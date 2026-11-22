#include "test_support.h"
#include "io.h"
#include "parse.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static const ReferenceBody initial[] = {
    {2, {1, 2, 3}, {4, 5, 6}, {2, 4, 6}},
    {3, {-1, 3, 2}, {-2, 1, 4}, {-1, 3, 5}}
};

typedef struct {
    ReferenceBody expected[257];
    size_t count;
    size_t calls;
    bool initializing;
} IntegrationCheck;

/* Independent scalar half-kick; drift is a separate pass in the reference. */
static void reference_kick(ReferenceBody *bodies, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        bodies[i].velocity.x += 0.125 * bodies[i].acceleration.x;
        bodies[i].velocity.y += 0.125 * bodies[i].acceleration.y;
        bodies[i].velocity.z += 0.125 * bodies[i].acceleration.z;
    }
}

static void assert_snapshot(const SimulationView *v, const IntegrationCheck *check)
{
    assert(v->count == check->count);
    for (size_t i = 0; i < v->count; ++i) {
        ReferenceBody actual = {v->w[i], {v->x[i], v->y[i], v->z[i]},
            {v->vx[i], v->vy[i], v->vz[i]}, {v->ax[i], v->ay[i], v->az[i]}};
        assert_body(actual, check->expected[i], 0);
    }
}

/* A position-dependent coupled force. Check the entire callback snapshot before
 * and after writing outputs, including the initial evaluation and half velocities. */
static const char *check_step(const SimulationView *v, double *ax, double *ay, double *az, void *context)
{
    IntegrationCheck *check = context;
#ifdef _OPENMP
    assert(!omp_in_parallel());
#endif
    if (!check->initializing) {
        reference_kick(check->expected, check->count);
        for (size_t i = 0; i < check->count; ++i) {
            ReferenceBody *b = &check->expected[i];
            b->position.x += 0.25 * b->velocity.x;
            b->position.y += 0.25 * b->velocity.y;
            b->position.z += 0.25 * b->velocity.z;
        }
    }
    assert_snapshot(v, check);
    for (size_t i = 0; i < v->count; ++i) {
        size_t j = (i + 1) % v->count;
        ax[i] = v->w[j] * (v->x[j] - v->x[i]);
        assert_snapshot(v, check);
        ay[i] = v->w[j] * (v->y[j] - v->y[i]);
        assert_snapshot(v, check);
        az[i] = v->w[j] * (v->z[j] - v->z[i]);
        assert_snapshot(v, check);
    }
    for (size_t i = 0; i < v->count; ++i) {
        check->expected[i].acceleration = (Vec3){ax[i], ay[i], az[i]};
    }
    if (!check->initializing) {
        reference_kick(check->expected, check->count);
    }
    check->initializing = false;
    ++check->calls;
    return NULL;
}

static void check_integration_contract(size_t count)
{
    ReferenceBody bodies[257];
    for (size_t i = 0; i < count; ++i) {
        bodies[i] = (ReferenceBody){1 + (double)(i % 3),
            {0.125 * (double)i, -0.25 * (double)(i % 7), 0.5 * (double)(i % 11)},
            {0.25, -0.5, 0.125}, {1000, -2000, 3000}};
    }
    for (size_t steps = 0; steps <= 4; ++steps) {
        for (size_t split = 0; split <= 1; ++split) {
            Simulation *s = from_bodies(bodies, count);
            IntegrationCheck check = {.count = count, .initializing = true};
            for (size_t i = 0; i < count; ++i) {
                check.expected[i] = bodies[i];
            }
            if (split) {
                for (size_t i = 0; i < steps; ++i) {
                    check.initializing = true;
                    assert(!simulate_with_solver(s, 0.25, 1, check_step, &check));
                }
            } else {
                assert(!simulate_with_solver(s, 0.25, steps, check_step, &check));
            }
            assert(check.calls == (steps == 0 || count == 0 ? 0 : split ? 2 * steps : steps + 1));
            for (size_t i = 0; i < count; ++i) {
                assert_body(body_at(s, i), check.expected[i], 0);
            }
            simulation_destroy(s);
        }
    }
}

static void test_integration_contract(void)
{
    const size_t counts[] = {0, 1, 2, 257};
#ifdef _OPENMP
    int saved_threads = omp_get_max_threads();
    int saved_dynamic = omp_get_dynamic();
    omp_set_dynamic(0);
    const int threads[] = {1, 4};
    for (size_t team = 0; team < 2; ++team) {
        omp_set_num_threads(threads[team]);
#endif
        for (size_t i = 0; i < sizeof(counts) / sizeof(counts[0]); ++i) {
            check_integration_contract(counts[i]);
        }
#ifdef _OPENMP
    }
    omp_set_num_threads(saved_threads);
    omp_set_dynamic(saved_dynamic);
#endif
}

static void test_acceleration_reinitialization(void)
{
    Simulation *a = from_bodies(initial, 2), *b = from_bodies(initial, 2);
    for (size_t run = 0; run < 2; ++run) {
        MutableSimulationView av = simulation_mutable_view(a);
        MutableSimulationView bv = simulation_mutable_view(b);
        if (run != 0) {
            /* Reinitialization also follows external mass/position changes. */
            av.w[0] = bv.w[0] = 7;
            av.x[1] = bv.x[1] = 10;
        }
        for (size_t i = 0; i < 2; ++i) {
            bv.ax[i] = 1000;
            bv.ay[i] = -2000;
            bv.az[i] = 3000;
        }
        IntegrationCheck first = {.expected = {body_at(a, 0), body_at(a, 1)}, .count = 2, .initializing = true};
        IntegrationCheck second = {.expected = {body_at(b, 0), body_at(b, 1)}, .count = 2, .initializing = true};
        assert(!simulate_with_solver(a, 0.25, 2, check_step, &first));
        assert(!simulate_with_solver(b, 0.25, 2, check_step, &second));
        assert(first.calls == 3 && second.calls == 3);
        for (size_t i = 0; i < 2; ++i) {
            assert_body(body_at(a, i), first.expected[i], 0);
            assert_body(body_at(b, i), second.expected[i], 0);
            assert_body(body_at(a, i), body_at(b, i), 0);
        }
    }
    simulation_destroy(a);
    simulation_destroy(b);
}

static void test_state_lifecycle(void)
{
    Simulation *state = from_bodies(initial, 2);
    assert(simulation_view(state).count == 2);
    for (size_t i = 0; i < 2; ++i) {
        assert_body(body_at(state, i), initial[i], 0);
    }
    simulation_destroy(state);
    Simulation *empty = from_bodies(NULL, 0);
    assert(simulation_view(empty).count == 0);
    simulation_destroy(empty);
    Simulation *failed;
    assert(!simulation_create(&failed, SIZE_MAX) && failed == NULL);
    simulation_destroy(failed);
}

static FILE *text_file(const char *text)
{
    FILE *file = tmpfile();
    assert(file && fputs(text, file) >= 0);
    rewind(file);
    return file;
}

static void test_io_contract(void)
{
    for (size_t count = 0; count <= 2; count += 2) {
        for (size_t steps = 0; steps <= 2; ++steps) {
            Simulation *s = from_bodies(initial, count), *loaded;
            IntegrationCheck check = {.expected = {initial[0], initial[1]}, .count = count, .initializing = true};
            assert(!simulate_with_solver(s, 0.25, steps, check_step, &check));
            FILE *file = tmpfile();
            assert(file && write_bodies(file, s));
            rewind(file);
            assert(read_bodies(file, &loaded));
            assert(simulation_view(loaded).count == count);
            for (size_t i = 0; i < count; ++i) {
                assert_body(body_at(s, i), body_at(loaded, i), 0);
            }
            check.initializing = true;
            IntegrationCheck resumed = check;
            assert(!simulate_with_solver(s, 0.25, 1, check_step, &check));
            assert(!simulate_with_solver(loaded, 0.25, 1, check_step, &resumed));
            for (size_t i = 0; i < count; ++i) {
                assert_body(body_at(s, i), body_at(loaded, i), 0);
            }
            assert(fclose(file) == 0);
            simulation_destroy(s);
            simulation_destroy(loaded);
        }
    }
    FILE *file = text_file("1\n0x1p0 1 2 3 4 5 6 7 8 9");
    Simulation *s;
    assert(read_bodies(file, &s)); /* Hexadecimal strtod input; no trailing newline. */
    assert_body(body_at(s, 0), (ReferenceBody){1, {1,2,3}, {4,5,6}, {7,8,9}}, 0);
    simulation_destroy(s);
    assert(fclose(file) == 0);
}

static void test_invalid_input(void)
{
    const char *invalid[] = {
        "", " \t\n", "-1", "1.5", "no", "0 extra", "0 1",
        "9999999999999999999999999999999999999999999999999999999999",
        "1\n1 0 0 0 0 0 0 0 0", /* Missing field. */
        "1\n1 0 0 0 0 0 0 0 0 0 7", /* Extra field. */
        "2\n1 0 0 0 0 0 0 0 0 0", /* Missing body. */
        "1\n0 0 0 0 0 0 0 0 0 0", /* Zero mass. */
        "1\n-1 0 0 0 0 0 0 0 0 0",
        "1\ninf 0 0 0 0 0 0 0 0 0",
        "1\n1 nan 0 0 0 0 0 0 0 0",
        "1\n1 0 0 0 1e999 0 0 0 0 0",
        "1\n1 0 0 0 0 0 0 0 0 text",
        "1\n1 0-1 0 0 0 0 0 0 0"
    };
    for (size_t i = 0; i < sizeof(invalid) / sizeof(invalid[0]); ++i) {
        FILE *file = text_file(invalid[i]);
        Simulation *simulation;
        assert(!read_bodies(file, &simulation));
        assert(simulation == NULL);
        assert(fclose(file) == 0);
    }

    char maximum[32];
    size_t count;
    assert(snprintf(maximum, sizeof(maximum), "%zu", SIZE_MAX) > 0);
    assert(parse_size(maximum, &count) && count == SIZE_MAX);
    FILE *file = text_file(maximum); /* Allocation size overflow. */
    Simulation *simulation;
    assert(!read_bodies(file, &simulation));
    assert(simulation == NULL);
    assert(fclose(file) == 0);

    file = tmpfile();
    assert(file);
    for (int i = 0; i < 200; ++i) {
        assert(fputc('0', file) != EOF);
    }
    rewind(file);
    assert(!read_bodies(file, &simulation));
    assert(fclose(file) == 0);
}

static const char *integrator_force(const SimulationView *state, double *ax, double *ay,
                             double *az, void *context)
{
    size_t *calls = context;
    ++*calls;
    assert(state->count == 1);
    ax[0] = 4;
    ay[0] = az[0] = 0;
    return NULL;
}

static void test_integrator_selection(void)
{
    const ReferenceBody body = {1, {1, 0, 0}, {2, 0, 0}, {12, 0, 0}};
    const char *(*integrators[])(Simulation *, double, size_t, AccelerationSolver, void *) = {
        simulate_with_solver, simulate_leapfrog_with_solver, simulate_euler_with_solver
    };
    for (size_t i = 0; i < sizeof(integrators) / sizeof(integrators[0]); ++i) {
        size_t calls = 0;
        Simulation *state = from_bodies(&body, 1);
        Simulation *empty = from_bodies(NULL, 0);
        assert(!integrators[i](empty, 0.5, 2, integrator_force, &calls));
        assert(!integrators[i](state, 0.5, 0, integrator_force, &calls));
        assert(calls == 0);
        assert_body(body_at(state, 0), body, 0);
        assert(!integrators[i](state, 0.5, 1, integrator_force, &calls));
        const bool euler = i == 2;
        assert(calls == (euler ? 1u : 2u));
        assert_body(body_at(state, 0),
                    (ReferenceBody){1, {euler ? 2 : 2.5, 0, 0},
                                    {euler ? 8 : 4, 0, 0}, {4, 0, 0}}, 0);
        simulation_destroy(empty);
        simulation_destroy(state);
    }
}

static const char force_error[] = "force evaluation failed";

static const char *fail_second_force(const SimulationView *state, double *ax, double *ay,
                                     double *az, void *context)
{
    size_t *calls = context;
    if (++*calls == 2) {
        ax[0] = 999;
        return force_error;
    }
    assert(state->count == 1);
    ax[0] = 2;
    ay[0] = 4;
    az[0] = 6;
    return NULL;
}

static void test_partial_force_is_discarded(void)
{
    for (int euler = 0; euler <= 1; ++euler) {
        ReferenceBody body = {1, {0}, {0}, {0}};
        Simulation *state = from_bodies(&body, 1);
        size_t calls = 0;
        const char *error = (euler ? simulate_euler_with_solver : simulate_with_solver)
                            (state, 1, 4, fail_second_force, &calls);
        assert(error == force_error && calls == 2);
        /* Leapfrog stopped after the drift and before its second half-kick.
         * Euler reached its second position and velocity update. */
        ReferenceBody expected = {1, {1, 2, 3}, {1, 2, 3}, {2, 4, 6}};
        if (euler) {
            expected.position = (Vec3){0};
            expected.velocity = (Vec3){2, 4, 6};
        }
        assert_body(body_at(state, 0), expected, 0);
        simulation_destroy(state);
    }
}

void test_state(void)
{
    test_integration_contract();
    test_integrator_selection();
    test_partial_force_is_discarded();
    test_acceleration_reinitialization();
    test_state_lifecycle();
    test_io_contract();
    test_invalid_input();
}
