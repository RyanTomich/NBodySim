#include "test_support.h"
#include "io.h"
#include "parse.h"
#include "solver.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>

static const ReferenceBody initial[] = {
    {2, {1, 2, 3}, {4, 5, 6}, {2, 4, 6}},
    {3, {-1, 3, 2}, {-2, 1, 4}, {-1, 3, 5}}
};

typedef struct { ReferenceBody expected[2]; size_t calls; } IntegrationCheck;

/* Check the callback's entire snapshot before and after writing each output.
 * This detects premature publication without assuming how buffers are stored. */
static void check_step(const SimulationView *v, double *ax, double *ay, double *az, void *context)
{
    IntegrationCheck *check = context;
    assert(v->count == 2);
    for (size_t i = 0; i < 2; ++i) {
        ReferenceBody *b = &check->expected[i];
        b->position.x += 0.25 * b->velocity.x;
        b->position.y += 0.25 * b->velocity.y;
        b->position.z += 0.25 * b->velocity.z;
        b->velocity.x += 0.25 * b->acceleration.x;
        b->velocity.y += 0.25 * b->acceleration.y;
        b->velocity.z += 0.25 * b->acceleration.z;
        ReferenceBody actual = {v->w[i], {v->x[i], v->y[i], v->z[i]},
            {v->vx[i], v->vy[i], v->vz[i]}, {v->ax[i], v->ay[i], v->az[i]}};
        assert_body(actual, *b, 0);
    }
    for (size_t i = 0; i < 2; ++i) {
        ax[i] = v->ax[i] + v->ay[1 - i];
        ay[i] = v->ay[i] - v->az[1 - i];
        az[i] = v->az[i] + v->ax[1 - i];
        for (size_t j = 0; j < 2; ++j) {
            assert(v->ax[j] == check->expected[j].acceleration.x);
            assert(v->ay[j] == check->expected[j].acceleration.y);
            assert(v->az[j] == check->expected[j].acceleration.z);
        }
    }
    for (size_t i = 0; i < 2; ++i) {
        check->expected[i].acceleration = (Vec3){ax[i], ay[i], az[i]};
    }
    ++check->calls;
}

static void test_integration_contract(void)
{
    for (size_t steps = 0; steps <= 4; ++steps) {
        for (size_t split = 0; split <= 1; ++split) {
            Simulation *s = from_bodies(initial, 2);
            IntegrationCheck check = {{initial[0], initial[1]}, 0};
            if (split) {
                for (size_t i = 0; i < steps; ++i) {
                    simulate_with_solver(s, 0.25, 1, check_step, &check);
                }
            } else {
                simulate_with_solver(s, 0.25, steps, check_step, &check);
            }
            assert(check.calls == steps);
            for (size_t i = 0; i < 2; ++i) {
                assert_body(body_at(s, i), check.expected[i], 0);
            }
            simulation_destroy(s);
        }
    }
    Simulation *empty = from_bodies(NULL, 0);
    IntegrationCheck check = {0};
    simulate_with_solver(empty, 0.25, 3, check_step, &check);
    assert(check.calls == 0);
    simulation_destroy(empty);
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
            IntegrationCheck check = {{initial[0], initial[1]}, 0};
            simulate_with_solver(s, 0.25, steps, check_step, &check);
            FILE *file = tmpfile();
            assert(file && write_bodies(file, s));
            rewind(file);
            assert(read_bodies(file, &loaded));
            assert(simulation_view(loaded).count == count);
            IntegrationCheck resumed = check;
            simulate_with_solver(s, 0.25, 1, check_step, &check);
            simulate_with_solver(loaded, 0.25, 1, check_step, &resumed);
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

static void constant_force(const SimulationView *view, double *ax, double *ay, double *az, void *context)
{
    ++*(size_t *)context;
    for (size_t i = 0; i < view->count; ++i) { ax[i] = 2; ay[i] = -1; az[i] = 0; }
}

static void test_leapfrog(void)
{
    ReferenceBody body = {1, {1, 2, 3}, {4, 5, 6}, {999, 999, 999}};
    Simulation *s = from_bodies(&body, 1);
    size_t calls = 0;
    simulate_leapfrog_with_solver(s, 0.1, 0, constant_force, &calls);
    assert(calls == 0);
    assert_body(body_at(s, 0), body, 0);
    simulate_leapfrog_with_solver(s, 0.1, 10, constant_force, &calls);
    assert(calls == 11);
    ReferenceBody expected = {1, {6, 6.5, 9}, {6, 4, 6}, {2, -1, 0}};
    assert_body(body_at(s, 0), expected, 1e-13);
    simulation_destroy(s);
    ReferenceBody pair[] = {{1, {-1, 0, 0}, {0, .2, 0}, {0}},
                            {2, {1, 0, 0}, {0, -.1, 0}, {0}}};
    s = from_bodies(pair, 2);
    Solver solver;
    const char *error;
    assert(solver_create(&solver, "n2", 2, NULL, 0, &error));
    simulate_leapfrog_with_solver(s, .001, 100, solver.acceleration, solver.context);
    MutableSimulationView v = simulation_mutable_view(s);
    for (size_t i = 0; i < 2; ++i) { v.vx[i] = -v.vx[i]; v.vy[i] = -v.vy[i]; v.vz[i] = -v.vz[i]; }
    simulate_leapfrog_with_solver(s, .001, 100, solver.acceleration, solver.context);
    for (size_t i = 0; i < 2; ++i) {
        ReferenceBody actual = body_at(s, i);
        assert_close(actual.position.x, pair[i].position.x, 1e-12);
        assert_close(actual.position.y, pair[i].position.y, 1e-12);
        assert_close(actual.velocity.x, -pair[i].velocity.x, 1e-12);
        assert_close(actual.velocity.y, -pair[i].velocity.y, 1e-12);
    }
    solver_destroy(&solver); simulation_destroy(s);
}

void test_state(void)
{
    test_leapfrog();
    test_integration_contract();
    test_state_lifecycle();
    test_io_contract();
    test_invalid_input();
}
