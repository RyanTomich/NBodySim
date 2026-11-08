/* Shared contract tests deliberately include public headers only. */
#include "solver.h"
#include "test_support.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static void test_preparation_errors(void)
{
    Solver solver;
    const char *error;
    SolverOption bad[] = {{"unknown", "1"}, {"theta", "nan"}, {"theta", "-1"}};
    assert(!solver_create(&solver, "missing", 1, NULL, 0, &error));
    assert(error && !solver.descriptor && !solver.context && !solver.acceleration);
    solver_destroy(&solver);
    assert(!solver_create(&solver, "direct", 1, NULL, 0, &error));
    assert(error && !solver.descriptor && !solver.context && !solver.acceleration);
    solver_destroy(&solver);
    size_t n2_count = 0;
    for (size_t i = 0; i < solver_count(); ++i) {
        assert(strcmp(solver_name_at(i), "direct") != 0);
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
    SolverOption invalid[] = {{"bucket_capacity", "0"}, {"max_depth", "-1"},
                             {"node_capacity", "-1"}};
    for (size_t i = 0; i < 3; ++i) {
        assert(!solver_create(&solver, "barnes-hut", 3, &invalid[i], 1, &error));
        solver_destroy(&solver);
    }
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

static void check_fallback_contract(const char *name, const SolverOption *option)
{
    /* Overflow in aggregate mass must fall back to direct evaluation and recover
     * on a later ordinary call. All expected accelerations remain representable. */
    ReferenceBody bodies[] = {{1e308, {0,0,0}, {0}, {0}},
                             {1e308, {1e100,0,0}, {0}, {0}},
                             {1e308, {2e100,0,0}, {0}, {0}}};
    Solver solver;
    const char *error;
    bool fmm = strcmp(name, "fmm") == 0;
    assert(solver_create(&solver, name, 3, option, option ? 1 : 0, &error));
    simulation_destroy(check_force(&solver, bodies, 3, 1e-11));
    MetadataCheck fallback = {.name = fmm ? "direct_fallbacks" : "aggregate_fallbacks"};
    solver_diagnostics(&solver, find_metadata, &fallback);
    assert(fallback.found && fallback.value.type == METADATA_SIZE && fallback.value.value.size == 1);
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
    MetadataCheck bounds = {.name = fmm ? "direct_fallbacks" : "bounds_fallbacks"};
    solver_diagnostics(&solver, find_metadata, &bounds);
    assert(bounds.found && bounds.value.type == METADATA_SIZE && bounds.value.value.size == (fmm ? 2u : 1u));
    solver_destroy(&n2);
    simulation_destroy(state);
    for (size_t i = 0; i < 3; ++i) {
        bodies[i].w = i + 1;
        bodies[i].position.x = (double)i;
    }
    simulation_destroy(check_force(&solver, bodies, 3, 1e-11));
    MetadataCheck calls = {.name = "calls"};
    solver_diagnostics(&solver, find_metadata, &calls);
    assert(calls.found && calls.value.type == METADATA_SIZE && calls.value.value.size == 3);
    fallback.found = false;
    solver_diagnostics(&solver, find_metadata, &fallback);
    assert(fallback.found && fallback.value.value.size == (fmm ? 2u : 1u));
    solver_destroy(&solver);
}

void test_solver(void)
{
    test_preparation_errors();
    test_configuration_contract();
    test_registered_solver_contract();
    test_fmm_cutoff_settings();
    test_fmm_distant_groups();
    test_fmm_sparse_cells();
    test_fmm_order_contract();
    test_tree_summary_reuse();
    test_gravity_accuracy();
    check_fallback_contract("barnes-hut", NULL);
    const char *orders[] = {"1", "2", "3", "4", "5", "6"};
    for (size_t i = 0; i < sizeof(orders) / sizeof(orders[0]); ++i) {
        SolverOption option = {"order", orders[i]};
        check_fallback_contract("fmm", &option);
    }
}
