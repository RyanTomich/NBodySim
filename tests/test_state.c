#include "test_support.h"
#include "io.h"
#include <assert.h>
#include <math.h>
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
    Simulation *source = from_bodies(initial, 2), *copy = from_bodies(initial, 2);
    for (size_t step = 0; step < 3; ++step) {
        IntegrationCheck check = {{initial[0], initial[1]}, 0};
        assert(simulation_copy(copy, source));
        simulate_with_solver(copy, 0.25, step, check_step, &check);
        MutableSimulationView v = simulation_mutable_view(copy);
        double *fields[] = {v.w, v.x, v.y, v.z, v.vx, v.vy, v.vz, v.ax, v.ay, v.az};
        for (size_t field = 0; field < 10; ++field) {
            fields[field][0] = NAN;
            assert(!simulation_is_finite(copy));
            assert(simulation_copy(copy, source));
            assert(simulation_is_finite(copy));
            for (size_t i = 0; i < 2; ++i) {
                assert_body(body_at(copy, i), initial[i], 0);
                assert_body(body_at(source, i), initial[i], 0);
            }
        }
        assert(simulation_copy(copy, copy));
    }
    Simulation *empty = from_bodies(NULL, 0);
    assert(simulation_view(empty).count == 0 && simulation_is_finite(empty));
    assert(simulation_copy(empty, empty));
    assert(!simulation_copy(copy, empty));
    assert_body(body_at(copy, 0), initial[0], 0);
    simulation_destroy(empty);
    simulation_destroy(copy);
    simulation_destroy(source);
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

void test_state(void)
{
    test_integration_contract();
    test_state_lifecycle();
    test_io_contract();
    test_invalid_input();
}
