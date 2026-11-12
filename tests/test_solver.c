/* Shared contract tests deliberately include public headers only. */
#include "solver.h"
#include "test_support.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static void test_preparation_errors(void)
{
    Solver solver;
    const char *error;
    SolverOption bad[] = {{"unknown", "1"}, {"theta", "nan"}, {"theta", "-1"}};
    assert(!solver_create(&solver, "missing", 1, NULL, 0, &error));
    assert(error && !solver.descriptor && !solver.context && !solver.acceleration);
    solver_destroy(&solver);
    size_t n2_count = 0;
    for (size_t i = 0; i < solver_count(); ++i) {
        n2_count += strcmp(solver_name_at(i), "n2") == 0;
    }
    assert(n2_count == 1);
    for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); ++i) {
        assert(!solver_create(&solver, "barnes-hut", 1, &bad[i], 1, &error));
        assert(error && !solver.context && !solver.descriptor);
        solver_destroy(&solver);
    }
    SolverOption duplicate[] = {{"theta", "0"}, {"theta", "1"}};
    assert(!solver_create(&solver, "barnes-hut", 1, duplicate, 2, &error));
    assert(!solver_create(&solver, "n2", 1, duplicate, 1, &error));
    assert(!solver_create(&solver, "barnes-hut", SIZE_MAX, NULL, 0, &error));
    assert(!solver.context && !solver.descriptor);
    solver_destroy(&solver);
    assert(solver_name_at(solver_count()) == NULL);
}

typedef struct { const char *name; MetadataValue value; bool found; } MetadataCheck;

static void find_metadata(const char *name, MetadataValue value, void *visitor_context)
{
    MetadataCheck *check = visitor_context;
    if (strcmp(name, check->name) == 0) {
        assert(!check->found); /* Keys must be unambiguous. */
        check->value = value;
        check->found = true;
    }
}

static void test_configuration_contract(void)
{
    SolverOption options[] = {{"theta", "0"}, {"bucket_capacity", "2"},
                              {"max_depth", "0"}, {"node_capacity", "1"}};
    const size_t expected[] = {0, 2, 0, 1};
    Solver solver;
    const char *error;
    assert(solver_create(&solver, "barnes-hut", 3, options, 4, &error));
    for (size_t i = 0; i < 4; ++i) {
        MetadataCheck check = {.name = options[i].name};
        solver_settings(&solver, find_metadata, &check);
        assert(check.found);
        if (i == 0) {
            assert(check.value.type == METADATA_DOUBLE && check.value.value.real == 0);
        } else {
            assert(check.value.type == METADATA_SIZE && check.value.value.size == expected[i]);
        }
    }
    solver_destroy(&solver);
    /* Omitted and explicitly automatic capacities report effective settings,
     * including the minimum allocation for an empty state. */
    const size_t counts[] = {0, 3};
    SolverOption automatic = {"node_capacity", "0"};
    for (size_t c = 0; c < 2; ++c) {
        for (size_t explicit_option = 0; explicit_option < 2; ++explicit_option) {
            assert(solver_create(&solver, "barnes-hut", counts[c], &automatic,
                                 explicit_option, &error));
            const size_t defaults[] = {0, 8, 32, counts[c] ? 4 * counts[c] : 1};
            for (size_t i = 0; i < 4; ++i) {
                MetadataCheck check = {.name = options[i].name};
                solver_settings(&solver, find_metadata, &check);
                assert(check.found);
                if (i == 0) {
                    assert(check.value.type == METADATA_DOUBLE && check.value.value.real == 0.5);
                } else {
                    assert(check.value.type == METADATA_SIZE && check.value.value.size == defaults[i]);
                }
            }
            solver_destroy(&solver);
        }
    }
    SolverOption invalid[] = {{"bucket_capacity", "0"}, {"max_depth", "-1"},
                             {"node_capacity", "-1"}};
    for (size_t i = 0; i < 3; ++i) {
        assert(!solver_create(&solver, "barnes-hut", 3, &invalid[i], 1, &error));
        solver_destroy(&solver);
    }
    SolverOption removed = {"reuse_tree", "1"};
    assert(!solver_create(&solver, "barnes-hut", 3, &removed, 1, &error));
    assert(error && !solver.context && !solver.descriptor);
    solver_destroy(&solver);
}

enum { CLOUD_SIZE = 129 };

static void fill_cloud(ReferenceBody *bodies)
{
    for (size_t i = 0; i < CLOUD_SIZE; ++i) {
        bodies[i] = (ReferenceBody){1 + (double)(i % 3),
            {(i % 2 ? -10 : 10) + 0.125 * (double)(i % 7),
             0.25 * (double)(i % 11), 0.25 * (double)(i % 13)}, {0}, {0}};
    }
    bodies[2].position = bodies[1].position; /* Unequal masses, coincident positions. */
}

/* Return the owned input so reuse tests can keep earlier states alive. */
static Simulation *check_force(Solver *solver, const ReferenceBody *bodies, size_t count, double tolerance)
{
    Simulation *s = from_bodies(bodies, count);
    SimulationView view = simulation_view(s);
    double ax[CLOUD_SIZE + 2], ay[CLOUD_SIZE + 2], az[CLOUD_SIZE + 2];
    Vec3 expected[CLOUD_SIZE];
    reference_accelerations(bodies, count, expected);
    for (size_t i = 0; i < count + 2; ++i) {
        ax[i] = ay[i] = az[i] = NAN;
    }
    ax[0] = ay[0] = az[0] = ax[count + 1] = ay[count + 1] = az[count + 1] = 12345;
    solver->acceleration(&view, ax + 1, ay + 1, az + 1, solver->context);
    double error = 0, norm = 0;
    for (size_t i = 0; i < count; ++i) {
        assert(isfinite(ax[i + 1]) && isfinite(ay[i + 1]) && isfinite(az[i + 1]));
        error = hypot(error, hypot(hypot(ax[i + 1] - expected[i].x,
                                        ay[i + 1] - expected[i].y), az[i + 1] - expected[i].z));
        norm = hypot(norm, hypot(hypot(expected[i].x, expected[i].y), expected[i].z));
        assert_body(body_at(s, i), bodies[i], 0); /* Input is immutable. */
    }
    assert_close(error, 0, tolerance * (1 + norm));
    assert(ax[0] == 12345 && ay[0] == 12345 && az[0] == 12345);
    assert(ax[count + 1] == 12345 && ay[count + 1] == 12345 && az[count + 1] == 12345);
    return s;
}

static void test_n2_first_row_initialization(void)
{
    const ReferenceBody bodies[] = {
        {2, {0, 0, 0}, {0}, {0}},
        {3, {0, 0, 0}, {0}, {0}},
        {5, {1, 2, 3}, {0}, {0}}
    };
    Solver solver;
    const char *error;
    assert(solver_create(&solver, "n2", 3, NULL, 0, &error));
    /* NaN-prefilled outputs must be initialized even when the first pair coincides. */
    simulation_destroy(check_force(&solver, bodies, 3, 1e-11));
    solver_destroy(&solver);
}

static void test_registered_solver_contract(void)
{
    ReferenceBody bodies[CLOUD_SIZE];
    for (size_t implementation = 0; implementation < solver_count(); ++implementation) {
        const char *name = solver_name_at(implementation);
        assert(name && *name);
        for (size_t previous = 0; previous < implementation; ++previous) {
            assert(strcmp(name, solver_name_at(previous)) != 0);
        }
        for (size_t count = 0; count <= 3; ++count) {
            fprintf(stderr, "solver contract: %s, count=%zu\n", name, count);
            fill_cloud(bodies);
            Solver solver;
            const char *error;
            assert(solver_create(&solver, name, count, NULL, 0, &error));
            Simulation *previous = check_force(&solver, bodies, count, 1e-11);
            /* Keep the first allocation alive: reuse must not read cached state. */
            ReferenceBody before[3];
            for (size_t i = 0; i < count; ++i) {
                before[i] = bodies[i];
                bodies[i].w *= 2;
                bodies[i].position = (Vec3){100, -10, 3};
            }
            simulation_destroy(check_force(&solver, bodies, count, 1e-11));
            for (size_t i = 0; i < count; ++i) {
                assert_body(body_at(previous, i), before[i], 0);
            }
            simulation_destroy(previous);
            solver_destroy(&solver);
            solver_destroy(&solver);
        }
    }
}

static void test_fmm_cutoff_settings(void)
{
    /* The default and explicit order four select the same optimization. */
    const char *orders[] = {NULL, "1", "2", "3", "4", "5", "6"};
    for (size_t i = 0; i < sizeof(orders) / sizeof(orders[0]); ++i) {
        Solver solver;
        const char *error;
        SolverOption order = {"order", orders[i]};
        assert(solver_create(&solver, "fmm", 1, &order, orders[i] ? 1 : 0, &error));
        MetadataCheck cutoff = {.name = "direct_pair_cutoff"};
        solver_settings(&solver, find_metadata, &cutoff);
        assert(cutoff.found && cutoff.value.type == METADATA_SIZE);
        assert(cutoff.value.value.size == (i == 0 || i == 4 ? 192u : 0u));
        solver_destroy(&solver);
    }
}

static void fill_distant_groups(ReferenceBody *bodies)
{
    /* Bodies at the same site exert zero acceleration on each other.
     * All of the reference force therefore comes from the other two sites. */
    const Vec3 sites[] = {{0, 0, 0}, {16, 9, -5}, {-4, 13, 21}};
    const size_t counts[] = {57, 43, 29};
    size_t count = 0;
    for (size_t site = 0; site < 3; ++site) {
        for (size_t j = 0; j < counts[site]; ++j) {
            bodies[count++] = (ReferenceBody){1 + (double)(j % 3), sites[site], {0}, {0}};
        }
    }
    assert(count == CLOUD_SIZE);
}

static void test_barnes_hut_topology_reuse(void)
{
    const char *thetas[] = {"0", "0.5"};
    for (size_t opening = 0; opening < 2; ++opening) {
        ReferenceBody bodies[CLOUD_SIZE], original[CLOUD_SIZE];
        fill_distant_groups(bodies);
        memcpy(original, bodies, sizeof(original));
        /* A depth limit leaves room for small movements of coincident
         * groups. Their far forces expose stale mass summaries. */
        SolverOption options[] = {{"theta", thetas[opening]}, {"max_depth", "4"}};
        Solver solver;
        const char *error;
        assert(solver_create(&solver, "barnes-hut", CLOUD_SIZE, options, 2, &error));
        double tolerance = opening ? 0.001 : 1e-11;
        Simulation *previous = check_force(&solver, bodies, CLOUD_SIZE, tolerance);
        for (size_t step = 0; step < 2; ++step) {
            for (size_t i = 0; i < CLOUD_SIZE; ++i) {
                bodies[i].w *= 2;
                bodies[i].position.x += 0.001 * (6 - bodies[i].position.x);
                bodies[i].position.y += 0.001 * (6.5 - bodies[i].position.y);
                bodies[i].position.z += 0.001 * (8 - bodies[i].position.z);
            }
            /* check_force supplies a different allocation of the same count. */
            simulation_destroy(check_force(&solver, bodies, CLOUD_SIZE, tolerance));
        }
        /* Move one late member, leaving its bucket's other bodies in place.
         * A representative-only containment check would miss this move. */
        bodies[CLOUD_SIZE - 1].position = bodies[0].position;
        bodies[CLOUD_SIZE - 1].w = 1000;
        simulation_destroy(check_force(&solver, bodies, CLOUD_SIZE, tolerance));
        for (size_t i = 0; i < CLOUD_SIZE; ++i) {
            assert_body(body_at(previous, i), original[i], 0);
        }
        simulation_destroy(previous);
        solver_destroy(&solver);
    }
}

static void test_barnes_hut_nonfinite_reuse(void)
{
    ReferenceBody bodies[] = {{1, {-1, -1, -1}, {0}, {0}},
                             {2, {0, 0, 0}, {0}, {0}}, {3, {1, 1, 1}, {0}, {0}}};
    Solver solver;
    const char *error;
    assert(solver_create(&solver, "barnes-hut", 3, NULL, 0, &error));
    Simulation *state = check_force(&solver, bodies, 3, 1e-11);
    MutableSimulationView mutable = simulation_mutable_view(state);
    double *coordinates[] = {mutable.x, mutable.y, mutable.z};
    const double invalid[] = {NAN, INFINITY, -INFINITY};
    for (size_t axis = 0; axis < 3; ++axis) {
        coordinates[axis][2] = invalid[axis];
        SimulationView view = simulation_view(state);
        double acceleration[3][3];
        solver.acceleration(&view, acceleration[0], acceleration[1], acceleration[2], solver.context);
        for (size_t component = 0; component < 3; ++component) {
            for (size_t i = 0; i < 3; ++i) {
                assert(isnan(acceleration[component][i]));
            }
        }
        coordinates[axis][2] = 1;
        simulation_destroy(check_force(&solver, bodies, 3, 1e-11));
    }
    simulation_destroy(state);
    solver_destroy(&solver);
}

static void test_fmm_distant_groups(void)
{
    ReferenceBody bodies[CLOUD_SIZE];
    fill_distant_groups(bodies);
    Solver solver;
    const char *error;
    assert(solver_create(&solver, "fmm", CLOUD_SIZE, NULL, 0, &error));
    simulation_destroy(check_force(&solver, bodies, CLOUD_SIZE, 0.01));
    solver_destroy(&solver);
}

static void test_fmm_sparse_cells(void)
{
    /* Two unevenly populated lines extend along different axes. Treating a
     * partially occupied cell as concentrated at its center loses far forces. */
    ReferenceBody bodies[CLOUD_SIZE];
    for (size_t i = 0; i < CLOUD_SIZE; ++i) {
        double k = (double)(i % 9);
        Vec3 point = i < 80 ? (Vec3){-1, k, 1} : (Vec3){4, 3, 2 * k};
        bodies[i] = (ReferenceBody){1 + (double)(i % 3), point, {0}, {0}};
    }
    Solver solver;
    const char *error;
    assert(solver_create(&solver, "fmm", CLOUD_SIZE, NULL, 0, &error));
    simulation_destroy(check_force(&solver, bodies, CLOUD_SIZE, 0.001));
    solver_destroy(&solver);
}

static void test_fmm_order_contract(void)
{
    const char *orders[] = {"1", "2", "3", "4", "5", "6"};
    ReferenceBody bodies[CLOUD_SIZE];
    fill_cloud(bodies);
    for (size_t order = 0; order < sizeof(orders) / sizeof(orders[0]); ++order) {
        for (size_t count = 0; count <= 3; ++count) {
            Solver solver;
            const char *error;
            SolverOption option = {"order", orders[order]};
            assert(solver_create(&solver, "fmm", count, &option, 1, &error));
            simulation_destroy(check_force(&solver, bodies, count, 1e-11));
            solver_destroy(&solver);
        }
    }
}

static void test_tree_summary_reuse(void)
{
    /* These fixed fixtures exercise approximations, including order four's
     * larger leaves. Budgets are fixture-specific, not general error bounds. */
    const char *orders[] = {NULL, "1", "2", "3", "4", "5", "6"};
    for (size_t variant = 0; variant < sizeof(orders) / sizeof(orders[0]); ++variant) {
        ReferenceBody bodies[CLOUD_SIZE];
        bool fmm = variant != 0;
        if (variant == 4) {
            fill_distant_groups(bodies);
        } else {
            fill_cloud(bodies);
        }
        double tolerance = variant == 4 ? 0.01 : 0.001;
        Solver solver;
        const char *error;
        SolverOption option = {"order", orders[variant]};
        assert(solver_create(&solver, fmm ? "fmm" : "barnes-hut", CLOUD_SIZE,
                             &option, fmm ? 1 : 0, &error));
        Simulation *previous = check_force(&solver, bodies, CLOUD_SIZE, tolerance);
        ReferenceBody before[CLOUD_SIZE];
        memcpy(before, bodies, sizeof(before));
        MetadataCheck approximations = {.name = fmm ? "translations" : "aggregate_evaluations"};
        solver_diagnostics(&solver, find_metadata, &approximations);
        assert(approximations.found && approximations.value.type == METADATA_SIZE);
        assert(approximations.value.value.size > 0);

        /* Rotate and translate the geometry, and change the source strengths.
         * Keep the old allocation alive to catch retained input pointers. */
        for (size_t i = 0; i < CLOUD_SIZE; ++i) {
            bodies[i].w *= 2;
            Vec3 point = bodies[i].position;
            bodies[i].position = (Vec3){point.y + 100, point.z - 20, point.x + 5};
        }
        simulation_destroy(check_force(&solver, bodies, CLOUD_SIZE, tolerance));
        for (size_t i = 0; i < CLOUD_SIZE; ++i) {
            assert_body(body_at(previous, i), before[i], 0);
        }
        simulation_destroy(previous);
        if (fmm) {
            /* Collapse and restore the cloud so rebuilt trees require very
             * different occupied extents in the same prepared storage. */
            for (size_t i = 0; i < CLOUD_SIZE; ++i) {
                bodies[i].position = (Vec3){0};
            }
            simulation_destroy(check_force(&solver, bodies, CLOUD_SIZE, 1e-11));
            simulation_destroy(check_force(&solver, before, CLOUD_SIZE, tolerance));
        }
        solver_destroy(&solver);
    }
}

static void test_gravity_accuracy(void)
{
    /* Known Newtonian answer anchors the independent reference's normalization. */
    ReferenceBody pair[] = {{2, {0,0,0}, {0}, {0}}, {3, {2,0,0}, {0}, {0}}};
    Vec3 exact[2];
    reference_accelerations(pair, 2, exact);
    assert_close(exact[0].x, 0.75, 0);
    assert_close(exact[1].x, -0.5, 0);

    const char *names[] = {"n2", "barnes-hut", "barnes-hut", "barnes-hut", "barnes-hut"};
    const SolverOption options[] = {{0}, {"theta", "0"}, {"node_capacity", "1"},
                                    {"max_depth", "0"}, {"bucket_capacity", "256"}};
    const size_t sizes[] = {2, 7, 33, CLOUD_SIZE};
    ReferenceBody bodies[CLOUD_SIZE];
    for (size_t variant = 0; variant < 5; ++variant) {
        for (size_t size = 0; size < sizeof(sizes) / sizeof(sizes[0]); ++size) {
            size_t count = sizes[size];
            fprintf(stderr, "exact gravity: %s, option=%s, count=%zu\n",
                    names[variant], variant ? options[variant].name : "none", count);
            fill_cloud(bodies);
            Solver solver;
            const char *error;
            assert(solver_create(&solver, names[variant], count, &options[variant], variant ? 1 : 0, &error));
            simulation_destroy(check_force(&solver, bodies, count, 1e-11));
            for (size_t i = 0; i < count; ++i) {
                bodies[i].position.x += 100;
            }
            simulation_destroy(check_force(&solver, bodies, count, 1e-11));
            /* Distinct coincident bodies still count as direct evaluations. */
            for (size_t i = 1; i < count; ++i) {
                bodies[i].position = bodies[0].position;
            }
            simulation_destroy(check_force(&solver, bodies, count, 1e-11));
            if (variant) {
                MetadataCheck direct = {.name = "direct_evaluations"};
                MetadataCheck aggregates = {.name = "aggregate_evaluations"};
                solver_diagnostics(&solver, find_metadata, &direct);
                solver_diagnostics(&solver, find_metadata, &aggregates);
                assert(direct.found && direct.value.type == METADATA_SIZE);
                assert(direct.value.value.size == 3 * count * (count - 1));
                assert(aggregates.found && aggregates.value.type == METADATA_SIZE);
                assert(aggregates.value.value.size == 0);
            }
            solver_destroy(&solver);
        }
    }
    /* A fixed regression budget for this cloud, not a universal theta error bound. */
    Solver barnes_hut;
    const char *error;
    fill_cloud(bodies);
    SolverOption opening = {"theta", "0.8"};
    assert(solver_create(&barnes_hut, "barnes-hut", CLOUD_SIZE, &opening, 1, &error));
    simulation_destroy(check_force(&barnes_hut, bodies, CLOUD_SIZE, 0.02));
    MetadataCheck calls = {.name = "calls"};
    solver_diagnostics(&barnes_hut, find_metadata, &calls);
    assert(calls.found && calls.value.type == METADATA_SIZE && calls.value.value.size == 1);
    MetadataCheck aggregates = {.name = "aggregate_evaluations"};
    solver_diagnostics(&barnes_hut, find_metadata, &aggregates);
    assert(aggregates.found && aggregates.value.type == METADATA_SIZE);
    assert(aggregates.value.value.size > 0); /* The approximation was exercised. */
    solver_destroy(&barnes_hut);
    /* Loose opening must still exclude the target's own mass. */
    const ReferenceBody separated[] = {{1, {0,0,0}, {0}, {0}},
                                       {2, {10,0,0}, {0}, {0}},
                                       {3, {11,0,0}, {0}, {0}}};
    SolverOption loose[] = {{"theta", "100"}, {"bucket_capacity", "1"}};
    assert(solver_create(&barnes_hut, "barnes-hut", 3, loose, 2, &error));
    simulation_destroy(check_force(&barnes_hut, separated, 3, 0.01));
    solver_destroy(&barnes_hut);

}

#ifdef _OPENMP
static void test_fmm_thread_consistency(void)
{
    enum { COUNT = 1536 };
    const struct { size_t count; int distribution; } cases[] = {
        {1023, 0}, {1024, 0}, {COUNT, 0}, {COUNT, 1}, {COUNT, 2}, {COUNT, 3}
    };
    const char *orders[] = {"1", "2", "3", "4", "5", "6"};
    int saved_threads = omp_get_max_threads(), saved_dynamic = omp_get_dynamic();
    omp_set_dynamic(0);
    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
        size_t count = cases[c].count;
        for (size_t order = 0; order < 6; ++order) {
            ReferenceBody bodies[COUNT];
            for (size_t i = 0; i < count; ++i) {
                Vec3 point = {(double)(i % 13) - 6, (double)(i / 13 % 13) - 6,
                              (double)(i / 169) - 4};
                if (cases[c].distribution == 1) {
                    double offset = i < count * 9 / 10 ? -20 : 20;
                    point = (Vec3){offset + point.x * 0.1, point.y * 0.1, point.z * 0.1};
                } else if (cases[c].distribution == 2) {
                    point = (Vec3){0.125 * (double)i, i % 2 ? -10 : 10, 0};
                } else if (cases[c].distribution == 3) {
                    point = (Vec3){0, 0, 0};
                }
                bodies[i] = (ReferenceBody){1 + (double)(i % 3), point, {0.01, 0, 0}, {0}};
            }
            Solver solver;
            SolverOption option = {"order", orders[order]};
            const char *error;
            /* FMM reads the runtime team setting on each evaluation. */
            omp_set_num_threads(1);
            assert(solver_create(&solver, "fmm", count, &option, 1, &error));
            for (size_t run = 0; run < 2; ++run) {
                Vec3 expected[COUNT];
                reference_accelerations(bodies, count, expected);
                Simulation *states[2];
                double actual[2][3][COUNT];
                for (size_t team = 0; team < 2; ++team) {
                    omp_set_num_threads(team ? 4 : 1);
                    states[team] = from_bodies(bodies, count);
                    SimulationView view = simulation_view(states[team]);
                    solver.acceleration(&view, actual[team][0], actual[team][1], actual[team][2],
                                        solver.context);
                    for (size_t i = 0; i < count; ++i) {
                        assert_body(body_at(states[team], i), bodies[i], 0);
                    }
                    simulate_with_solver(states[team], 0.00001, 2, solver.acceleration, solver.context);
                }
                double error_norm = 0, reference_norm = 0;
                for (size_t i = 0; i < count; ++i) {
                    for (size_t axis = 0; axis < 3; ++axis) {
                        assert_close(actual[1][axis][i], actual[0][axis][i], 0);
                    }
                    error_norm = hypot(error_norm, hypot(hypot(actual[0][0][i] - expected[i].x,
                                      actual[0][1][i] - expected[i].y), actual[0][2][i] - expected[i].z));
                    reference_norm = hypot(reference_norm, hypot(hypot(expected[i].x, expected[i].y),
                                                                 expected[i].z));
                    assert_body(body_at(states[0], i), body_at(states[1], i), 0);
                    /* Reuse the solver with a new allocation, masses, and positions. */
                    bodies[i].w *= 1.25;
                    bodies[i].position.x *= 1.1;
                    bodies[i].position.y += 0.25;
                }
                /* Fixture-specific accuracy budget, including order one. */
                assert_close(error_norm, 0, 0.1 * (1 + reference_norm));
                simulation_destroy(states[0]);
                simulation_destroy(states[1]);
            }
            solver_destroy(&solver);
        }
    }
    omp_set_num_threads(saved_threads);
    omp_set_dynamic(saved_dynamic);
}

static void test_fmm_parallel_fallback(void)
{
    enum { COUNT = 1024 };
    ReferenceBody bodies[COUNT];
    int saved_threads = omp_get_max_threads(), saved_dynamic = omp_get_dynamic();
    omp_set_dynamic(0);
    omp_set_num_threads(4);
    const char *orders[] = {"1", "2", "3", "4", "5", "6"};
    for (size_t order = 0; order < 6; ++order) {
        Solver solver;
        SolverOption option = {"order", orders[order]};
        const char *error;
        assert(solver_create(&solver, "fmm", COUNT, &option, 1, &error));
        for (size_t run = 0; run < 2; ++run) {
            /* Finite tree summaries, but overflowing distances between opposite
             * corners fail during interactions. The next call must recover. */
            double scale = run ? 1 : 8e307;
            for (size_t i = 0; i < COUNT; ++i) {
                bodies[i] = (ReferenceBody){1, {i & 1 ? scale : -scale,
                                                i & 2 ? scale : -scale,
                                                i & 4 ? scale : -scale}, {0}, {0}};
            }
            Simulation *state = from_bodies(bodies, COUNT);
            SimulationView view = simulation_view(state);
            double actual[3][COUNT], expected[3][COUNT];
            solver.acceleration(&view, actual[0], actual[1], actual[2], solver.context);
            omp_set_num_threads(1);
            /* Failed calls match direct summation; recovery matches a fresh
             * serial FMM, including the low-order approximation error. */
            Solver reference;
            assert(solver_create(&reference, run ? "fmm" : "n2", COUNT,
                                 run ? &option : NULL, run ? 1 : 0, &error));
            reference.acceleration(&view, expected[0], expected[1], expected[2], reference.context);
            solver_destroy(&reference);
            omp_set_num_threads(4);
            for (size_t axis = 0; axis < 3; ++axis) {
                for (size_t i = 0; i < COUNT; ++i) {
                    if (!run) {
                        assert((isnan(actual[axis][i]) && isnan(expected[axis][i]))
                               || actual[axis][i] == expected[axis][i]);
                    } else {
                        assert_close(actual[axis][i], expected[axis][i], 0);
                    }
                }
            }
            MetadataCheck fallback = {.name = "direct_fallbacks"};
            solver_diagnostics(&solver, find_metadata, &fallback);
            assert(fallback.found && fallback.value.value.size == 1);
            simulation_destroy(state);
        }
        solver_destroy(&solver);
    }
    omp_set_num_threads(saved_threads);
    omp_set_dynamic(saved_dynamic);
}

static void test_n2_thread_consistency(void)
{
    int saved_threads = omp_get_max_threads();
    int saved_dynamic = omp_get_dynamic();
    int saved_levels = omp_get_max_active_levels();
    omp_set_dynamic(0);
    omp_set_max_active_levels(1);
    const size_t counts[] = {0, 1, 2, 3, CLOUD_SIZE};
    for (size_t c = 0; c < sizeof(counts) / sizeof(counts[0]); ++c) {
        size_t count = counts[c];
        ReferenceBody bodies[CLOUD_SIZE];
        fill_cloud(bodies);
        Solver solvers[2];
        const char *error;
        for (size_t team = 0; team < 2; ++team) {
            omp_set_num_threads(team == 0 ? 1 : 4);
            assert(solver_create(&solvers[team], "n2", count, NULL, 0, &error));
        }
        Simulation *previous = NULL;
        for (size_t repeat = 0; repeat < 3; ++repeat) {
            Vec3 expected[CLOUD_SIZE];
            reference_accelerations(bodies, count, expected);
            for (size_t team = 0; team < 2; ++team) {
                /* Changing the runtime default cannot enlarge prepared storage
                 * or change the single-worker solver into a parallel solver. */
                omp_set_num_threads(team == 0 ? 4 : 1);
                Simulation *state = check_force(&solvers[team], bodies, count, 1e-11);
                SimulationView view = simulation_view(state);
                double actual[3][CLOUD_SIZE];
                solvers[team].acceleration(&view, actual[0], actual[1], actual[2], solvers[team].context);
                for (size_t i = 0; i < count; ++i) {
                    assert_close(actual[0][i], expected[i].x, 1e-11);
                    assert_close(actual[1][i], expected[i].y, 1e-11);
                    assert_close(actual[2][i], expected[i].z, 1e-11);
                }
                /* Keep the earlier allocation alive until after the next call. */
                simulation_destroy(previous);
                previous = state;
            }
            /* A nested region must serialize, even for a solver prepared for
             * four workers. Scratch from the preceding call must be ignored. */
#pragma omp parallel num_threads(2) default(none) shared(solvers, bodies, count)
            {
#pragma omp single
                {
                    simulation_destroy(check_force(&solvers[1], bodies, count, 1e-11));
                }
            }
            for (size_t i = 0; i < count; ++i) {
                bodies[i].w *= 1.5;
                bodies[i].position.x += 0.03125 * (double)(i % 5);
                bodies[i].position.z *= 1.125;
            }
        }
        simulation_destroy(previous);
        solver_destroy(&solvers[0]);
        solver_destroy(&solvers[1]);
    }
    omp_set_num_threads(4);
    Solver solver;
    const char *error;
    const size_t overflowing[] = {SIZE_MAX, SIZE_MAX / (3 * sizeof(double)) / 4 + 1};
    for (size_t i = 0; i < sizeof(overflowing) / sizeof(overflowing[0]); ++i) {
        assert(!solver_create(&solver, "n2", overflowing[i], NULL, 0, &error));
        assert(error && !solver.context && !solver.descriptor && !solver.acceleration);
        solver_destroy(&solver);
    }
    omp_set_max_active_levels(saved_levels);
    omp_set_num_threads(saved_threads);
    omp_set_dynamic(saved_dynamic);
}

static void test_trajectory_thread_consistency(void)
{
    int saved_threads = omp_get_max_threads();
    int saved_dynamic = omp_get_dynamic();
    omp_set_dynamic(0);
    ReferenceBody bodies[CLOUD_SIZE];
    fill_cloud(bodies);
    for (size_t i = 0; i < CLOUD_SIZE; ++i) {
        bodies[i].velocity = (Vec3){0.125 * (double)(i % 3), -0.25, 0.5};
        bodies[i].acceleration = (Vec3){1000, -2000, 3000};
    }
    const char *names[] = {"n2", "barnes-hut", "fmm"};
    const SolverOption options[] = {{"theta", "0.5"}, {"order", "4"}};
    for (size_t backend = 0; backend < 3; ++backend) {
        Simulation *states[2];
        for (size_t team = 0; team < 2; ++team) {
            omp_set_num_threads(team == 0 ? 1 : 4);
            states[team] = from_bodies(bodies, CLOUD_SIZE);
            Solver solver;
            const char *error;
            assert(solver_create(&solver, names[backend], CLOUD_SIZE,
                                backend ? &options[backend - 1] : NULL,
                                backend ? 1 : 0, &error));
            for (size_t run = 0; run < 2; ++run) {
                simulate_with_solver(states[team], 0.0001, 5,
                                     solver.acceleration, solver.context);
            }
            solver_destroy(&solver);
        }
        for (size_t i = 0; i < CLOUD_SIZE; ++i) {
            assert_body(body_at(states[0], i), body_at(states[1], i), backend == 0 ? 1e-11 : 0);
        }
        simulation_destroy(states[0]);
        simulation_destroy(states[1]);
    }
    omp_set_num_threads(saved_threads);
    omp_set_dynamic(saved_dynamic);
}

static void test_barnes_hut_thread_consistency(void)
{
    int saved_threads = omp_get_max_threads();
    int saved_dynamic = omp_get_dynamic();
    int saved_levels = omp_get_max_active_levels();
    omp_set_dynamic(0);
    omp_set_max_active_levels(1);
    ReferenceBody bodies[CLOUD_SIZE];
    fill_cloud(bodies);
    Simulation *state = from_bodies(bodies, CLOUD_SIZE);
    SimulationView view = simulation_view(state);
    const char *thetas[] = {"0", "0.8"};
    const double tolerances[] = {1e-11, 0.02};
    const char *counter_names[] = {
        "source_node_visits", "direct_evaluations", "aggregate_evaluations"
    };
    for (size_t variant = 0; variant < 2; ++variant) {
        Solver solvers[2];
        SolverOption opening = {"theta", thetas[variant]};
        const char *error;
        for (size_t team = 0; team < 2; ++team) {
            omp_set_num_threads(team == 0 ? 1 : 4);
            assert(solver_create(&solvers[team], "barnes-hut", CLOUD_SIZE, &opening, 1, &error));
            simulation_destroy(check_force(&solvers[team], bodies, CLOUD_SIZE, tolerances[variant]));
        }
        size_t per_call[3];
        for (size_t field = 0; field < 3; ++field) {
            for (size_t team = 0; team < 2; ++team) {
                MetadataCheck counter = {.name = counter_names[field]};
                solver_diagnostics(&solvers[team], find_metadata, &counter);
                assert(counter.found && counter.value.type == METADATA_SIZE);
                if (team == 0) {
                    per_call[field] = counter.value.value.size;
                } else {
                    assert(counter.value.value.size == per_call[field]);
                }
            }
        }
        double expected[3][CLOUD_SIZE], actual[3][CLOUD_SIZE];
        /* Runtime defaults cannot enlarge prepared storage. A nested call
         * serializes, so stale statistics from inactive workers must be ignored. */
        for (size_t repeat = 0; repeat < 3; ++repeat) {
            omp_set_num_threads(4);
            solvers[0].acceleration(&view, expected[0], expected[1], expected[2], solvers[0].context);
            omp_set_num_threads(1);
            if (repeat == 1) {
#pragma omp parallel num_threads(2) default(none) shared(solvers, view, actual)
                {
#pragma omp single
                    {
                        solvers[1].acceleration(&view, actual[0], actual[1], actual[2], solvers[1].context);
                    }
                }
            } else {
                solvers[1].acceleration(&view, actual[0], actual[1], actual[2], solvers[1].context);
            }
            for (size_t axis = 0; axis < 3; ++axis) {
                for (size_t i = 0; i < CLOUD_SIZE; ++i) {
                    assert(isfinite(actual[axis][i]));
                    assert(actual[axis][i] == expected[axis][i]);
                }
            }
            /* Identical inputs contribute the same counts on each call,
             * independently of team size and which worker evaluates a leaf. */
            for (size_t field = 0; field < 3; ++field) {
                for (size_t team = 0; team < 2; ++team) {
                    MetadataCheck counter = {.name = counter_names[field]};
                    solver_diagnostics(&solvers[team], find_metadata, &counter);
                    assert(counter.found && counter.value.type == METADATA_SIZE);
                    assert(counter.value.value.size == (repeat + 2) * per_call[field]);
                }
            }
        }
        for (size_t team = 0; team < 2; ++team) {
            MetadataCheck calls = {.name = "calls"};
            solver_diagnostics(&solvers[team], find_metadata, &calls);
            assert(calls.found && calls.value.type == METADATA_SIZE && calls.value.value.size == 4);
#ifdef NBODY_PROFILE_BH
            const char *profile_names[] = {"worker_slots", "min_threads", "max_threads"};
            const size_t expected_profile[] = {team ? 13 : 4, 1, team ? 4 : 1};
            for (size_t field = 0; field < 3; ++field) {
                MetadataCheck profile = {.name = profile_names[field]};
                solver_diagnostics(&solvers[team], find_metadata, &profile);
                assert(profile.found && profile.value.type == METADATA_SIZE);
                assert(profile.value.value.size == expected_profile[field]);
            }
#endif
            solver_destroy(&solvers[team]);
        }
    }
    simulation_destroy(state);
    omp_set_max_active_levels(saved_levels);
    omp_set_num_threads(saved_threads);
    omp_set_dynamic(saved_dynamic);
}
#endif

static void test_barnes_hut_monopole_acceleration(void)
{
    /* Two coincident pairs make each accepted source an exact point mass.
     * Small buckets allow aggregates without exposing private tree structure. */
    struct {
        const char *name;
        double mass;
        Vec3 displacement;
    } cases[] = {
        {"ordinary", 3, {2, -3, 4}},
        {"x axis", 1, {-2, 0, 0}},
        {"y axis", 1, {0, 2, 0}},
        {"z axis", 1, {0, 0, -2}},
    };
    SolverOption options[] = {{"theta", "0.5"}, {"bucket_capacity", "1"},
                              {"max_depth", "8"}, {"node_capacity", "128"}};
    Solver solver;
    const char *error;
    assert(solver_create(&solver, "barnes-hut", 4, options, 4, &error));
    size_t previous_aggregates = 0;
    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
        for (int repeat = 0; repeat < 3; ++repeat) {
            /* Refresh masses at fixed positions, then move the sources.
             * Reuse one solver across all cases. */
            double mass = cases[c].mass * (repeat == 1 ? 0.5 : 1.0);
            double direction = repeat == 2 ? -1.0 : 1.0;
            Vec3 position = {direction * cases[c].displacement.x,
                             direction * cases[c].displacement.y,
                             direction * cases[c].displacement.z};
            ReferenceBody bodies[4] = {
                {mass * 0.5, {0}, {0}, {0}}, {mass * 0.5, {0}, {0}, {0}},
                {mass * 0.5, position, {0}, {0}}, {mass * 0.5, position, {0}, {0}}
            };
            Simulation *state = from_bodies(bodies, 4);
            SimulationView view = simulation_view(state);
            Vec3 expected[4];
            double actual[3][4];
            reference_accelerations(bodies, 4, expected);
            solver.acceleration(&view, actual[0], actual[1], actual[2], solver.context);
            MetadataCheck aggregates = {.name = "aggregate_evaluations"};
            solver_diagnostics(&solver, find_metadata, &aggregates);
            assert(aggregates.found && aggregates.value.type == METADATA_SIZE);
            assert(aggregates.value.value.size > previous_aggregates);
            previous_aggregates = aggregates.value.value.size;
            for (size_t i = 0; i < 4; ++i) {
                const double components[] = {expected[i].x, expected[i].y, expected[i].z};
                for (size_t axis = 0; axis < 3; ++axis) {
                    double want = components[axis], got = actual[axis][i];
                    double tolerance = 64 * DBL_EPSILON * fabs(want);
                    if (!isfinite(got) || !isfinite(want)
                        || (want != 0.0 && got == 0.0)
                        || fabs(got - want) > tolerance) {
                        fprintf(stderr, "monopole %s repeat=%d axis=%zu: "
                                "actual %.17g expected %.17g\n",
                                cases[c].name, repeat, axis, got, want);
                        assert(false);
                    }
                }
                assert_body(body_at(state, i), bodies[i], 0);
            }
            simulation_destroy(state);
        }
    }
    solver_destroy(&solver);
}

static void test_barnes_hut_unrepresentable_tree(void)
{
    const ReferenceBody cases[][3] = {
        {{1e308, {0, 0, 0}, {0}, {0}},
         {1e308, {1e100, 0, 0}, {0}, {0}},
         {1e308, {2e100, 0, 0}, {0}, {0}}},
        {{1, {-1e308, 0, 0}, {0}, {0}},
         {1, {0, 0, 0}, {0}, {0}},
         {1, {1e308, 0, 0}, {0}, {0}}}
    };
    ReferenceBody ordinary[] = {{1, {0, 0, 0}, {0}, {0}},
                                {2, {1e100, 0, 0}, {0}, {0}},
                                {3, {2e100, 0, 0}, {0}, {0}}};
    Solver solver;
    const char *error;
    assert(solver_create(&solver, "barnes-hut", 3, NULL, 0, &error));
    simulation_destroy(check_force(&solver, ordinary, 3, 1e-11));
    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
        MetadataCheck contributions[] = {
            {.name = "direct_evaluations"}, {.name = "aggregate_evaluations"}
        };
        for (size_t field = 0; field < 2; ++field) {
            solver_diagnostics(&solver, find_metadata, &contributions[field]);
            assert(contributions[field].found && contributions[field].value.type == METADATA_SIZE);
        }
        Simulation *state = from_bodies(cases[c], 3);
        SimulationView view = simulation_view(state);
        double actual[3][3] = {{0}};
        solver.acceleration(&view, actual[0], actual[1], actual[2], solver.context);
        for (size_t axis = 0; axis < 3; ++axis) {
            for (size_t i = 0; i < 3; ++i) {
                assert(isnan(actual[axis][i]));
            }
        }
        for (size_t i = 0; i < 3; ++i) {
            assert_body(body_at(state, i), cases[c][i], 0);
        }
        /* A failed call must not recount contributions from the previous call. */
        for (size_t field = 0; field < 2; ++field) {
            MetadataCheck after = {.name = contributions[field].name};
            solver_diagnostics(&solver, find_metadata, &after);
            assert(after.found && after.value.type == METADATA_SIZE);
            assert(after.value.value.size == contributions[field].value.value.size);
        }
        simulation_destroy(state);
        /* Recover after both an invalid summary on reused topology and a
         * failed rebuild, without recreating the solver. */
        simulation_destroy(check_force(&solver, ordinary, 3, 1e-11));
    }
    MetadataCheck calls = {.name = "calls"};
    solver_diagnostics(&solver, find_metadata, &calls);
    assert(calls.found && calls.value.type == METADATA_SIZE && calls.value.value.size == 5);
    solver_destroy(&solver);
}

static void check_fmm_fallback_contract(const SolverOption *option)
{
    /* Overflow in aggregate mass must fall back to direct evaluation and recover
     * on a later ordinary call. All expected accelerations remain representable. */
    ReferenceBody bodies[] = {{1e308, {0,0,0}, {0}, {0}},
                             {1e308, {1e100,0,0}, {0}, {0}},
                             {1e308, {2e100,0,0}, {0}, {0}}};
    Solver solver;
    const char *error;
    assert(solver_create(&solver, "fmm", 3, option, 1, &error));
    for (size_t i = 0; i < 3; ++i) {
        bodies[i].w = 1;
    }
    simulation_destroy(check_force(&solver, bodies, 3, 1e-11));
    for (size_t i = 0; i < 3; ++i) {
        bodies[i].w = 1e308;
    }
    simulation_destroy(check_force(&solver, bodies, 3, 1e-11));
    MetadataCheck fallback = {.name = "direct_fallbacks"};
    solver_diagnostics(&solver, find_metadata, &fallback);
    assert(fallback.found && fallback.value.type == METADATA_SIZE && fallback.value.value.size == 1);
    /* Recover at the same positions. */
    for (size_t i = 0; i < 3; ++i) {
        bodies[i].w = 1;
    }
    simulation_destroy(check_force(&solver, bodies, 3, 1e-11));
    /* Unrepresentable bounds retain the n2 solver's overflow behavior. */
    bodies[0].position.x = -1e308;
    bodies[1].position.x = 0;
    bodies[2].position.x = 1e308;
    for (size_t i = 0; i < 3; ++i) {
        bodies[i].w = 1;
    }
    Simulation *state = from_bodies(bodies, 3);
    SimulationView view = simulation_view(state);
    Solver n2;
    assert(solver_create(&n2, "n2", 3, NULL, 0, &error));
    double actual[3][3], expected[3][3];
    solver.acceleration(&view, actual[0], actual[1], actual[2], solver.context);
    n2.acceleration(&view, expected[0], expected[1], expected[2], n2.context);
    for (size_t axis = 0; axis < 3; ++axis) {
        for (size_t i = 0; i < 3; ++i) {
            assert((isnan(actual[axis][i]) && isnan(expected[axis][i]))
                   || actual[axis][i] == expected[axis][i]);
        }
    }
    MetadataCheck bounds = {.name = "direct_fallbacks"};
    solver_diagnostics(&solver, find_metadata, &bounds);
    assert(bounds.found && bounds.value.type == METADATA_SIZE && bounds.value.value.size == 2);
    solver_destroy(&n2);
    simulation_destroy(state);
    for (size_t i = 0; i < 3; ++i) {
        bodies[i].w = i + 1;
        bodies[i].position.x = (double)i;
    }
    simulation_destroy(check_force(&solver, bodies, 3, 1e-11));
    MetadataCheck calls = {.name = "calls"};
    solver_diagnostics(&solver, find_metadata, &calls);
    assert(calls.found && calls.value.type == METADATA_SIZE && calls.value.value.size == 5);
    fallback.found = false;
    solver_diagnostics(&solver, find_metadata, &fallback);
    assert(fallback.found && fallback.value.value.size == 2);
    solver_destroy(&solver);
}

void test_solver(void)
{
    test_preparation_errors();
    test_configuration_contract();
    test_registered_solver_contract();
    test_n2_first_row_initialization();
    test_barnes_hut_topology_reuse();
    test_barnes_hut_nonfinite_reuse();
    test_fmm_cutoff_settings();
    test_fmm_distant_groups();
    test_fmm_sparse_cells();
    test_fmm_order_contract();
    test_tree_summary_reuse();
    test_gravity_accuracy();
    test_barnes_hut_monopole_acceleration();
#ifdef _OPENMP
    test_fmm_thread_consistency();
    test_fmm_parallel_fallback();
    test_n2_thread_consistency();
    test_trajectory_thread_consistency();
    test_barnes_hut_thread_consistency();
#endif
    test_barnes_hut_unrepresentable_tree();
    const char *orders[] = {"1", "2", "3", "4", "5", "6"};
    for (size_t i = 0; i < sizeof(orders) / sizeof(orders[0]); ++i) {
        SolverOption option = {"order", orders[i]};
        check_fmm_fallback_contract(&option);
    }
}
