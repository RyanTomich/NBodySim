#include "gravity.h"
#include "approximate.h"
#include "io.h"
#include "nbody.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* AoS exists only in test fixtures and the independent reference solver. */
typedef struct {
    double w;
    Vec3 position, velocity, acceleration;
} ReferenceBody;

static Simulation from_bodies(const ReferenceBody *bodies, size_t count)
{
    Simulation simulation;
    assert(simulation_create(&simulation, count));
    for (size_t i = 0; i < count; ++i) {
        simulation.w[i] = bodies[i].w;
        simulation.x[i] = bodies[i].position.x;
        simulation.y[i] = bodies[i].position.y;
        simulation.z[i] = bodies[i].position.z;
        simulation.vx[i] = bodies[i].velocity.x;
        simulation.vy[i] = bodies[i].velocity.y;
        simulation.vz[i] = bodies[i].velocity.z;
        simulation.ax[i] = bodies[i].acceleration.x;
        simulation.ay[i] = bodies[i].acceleration.y;
        simulation.az[i] = bodies[i].acceleration.z;
        /* Neither kernel may depend on previously initialized workspace. */
        simulation.next_ax[i] = NAN;
        simulation.next_ay[i] = NAN;
        simulation.next_az[i] = NAN;
    }
    return simulation;
}

static ReferenceBody view_body(const SimulationView *state, size_t i)
{
    return (ReferenceBody){
        state->w[i], {state->x[i], state->y[i], state->z[i]},
        {state->vx[i], state->vy[i], state->vz[i]},
        {state->ax[i], state->ay[i], state->az[i]}
    };
}

static ReferenceBody body_at(const Simulation *simulation, size_t i)
{
    SimulationView view = simulation_view(simulation);
    return view_body(&view, i);
}

static void assert_vec(Vec3 actual, Vec3 expected)
{
    assert(fabs(actual.x - expected.x) < 1e-12);
    assert(fabs(actual.y - expected.y) < 1e-12);
    assert(fabs(actual.z - expected.z) < 1e-12);
}

static void assert_body_equal(ReferenceBody actual, ReferenceBody expected)
{
    assert(actual.w == expected.w);
    assert(actual.position.x == expected.position.x);
    assert(actual.position.y == expected.position.y);
    assert(actual.position.z == expected.position.z);
    assert(actual.velocity.x == expected.velocity.x);
    assert(actual.velocity.y == expected.velocity.y);
    assert(actual.velocity.z == expected.velocity.z);
    assert(actual.acceleration.x == expected.acceleration.x);
    assert(actual.acceleration.y == expected.acceleration.y);
    assert(actual.acceleration.z == expected.acceleration.z);
}

static void test_gravity(void)
{
    const ReferenceBody initial[] = {
        {.w = 7.0, .position = {1.0, 1.0, 1.0}},
        {.w = 27.0, .position = {2.0, 3.0, 3.0}}
    };
    Simulation simulation = from_bodies(initial, 2);
    SimulationView state = simulation_view(&simulation);
    assert_vec(gravity(&state, 0, 1), (Vec3){1.0, 2.0, 2.0});

    simulation.w[0] = 99.0; /* Acceleration does not depend on target mass. */
    assert_vec(gravity(&state, 0, 1), (Vec3){1.0, 2.0, 2.0});
    simulation.w[1] = 54.0;
    assert_vec(gravity(&state, 0, 1), (Vec3){2.0, 4.0, 4.0});

    simulation.w[1] = 27.0;
    simulation.x[1] = 3.0; /* Twice the distance. */
    simulation.y[1] = 5.0;
    simulation.z[1] = 5.0;
    assert_vec(gravity(&state, 0, 1), (Vec3){0.25, 0.5, 0.5});
    simulation.x[1] = 0.0;
    simulation.y[1] = -1.0;
    simulation.z[1] = -1.0;
    assert_vec(gravity(&state, 0, 1), (Vec3){-1.0, -2.0, -2.0});

    simulation.x[1] = simulation.x[0];
    simulation.y[1] = simulation.y[0];
    simulation.z[1] = simulation.z[0];
    assert_vec(gravity(&state, 0, 1), (Vec3){0.0, 0.0, 0.0});
    free(simulation.storage);
}

static void test_euler(void)
{
    ReferenceBody body = {
        .w = 1.0,
        .position = {1.0, 2.0, 3.0},
        .velocity = {4.0, 5.0, 6.0},
        .acceleration = {2.0, 4.0, 6.0}
    };
    Simulation simulation = from_bodies(&body, 1);
    simulate(&simulation, 0.25, 1, gravity);
    body = body_at(&simulation, 0);
    assert_vec(body.position, (Vec3){2.0, 3.25, 4.5});
    assert_vec(body.velocity, (Vec3){4.5, 6.0, 7.5});
    assert_vec(body.acceleration, (Vec3){0.0, 0.0, 0.0});

    simulate(&simulation, 0.25, 2, gravity);
    body = body_at(&simulation, 0);
    free(simulation.storage);
    assert_vec(body.position, (Vec3){4.25, 6.25, 8.25});
    assert_vec(body.velocity, (Vec3){4.5, 6.0, 7.5});
}

static void test_input_acceleration(void)
{
    ReferenceBody bodies[] = {
        {.w = 1.0, .position = {0.0, 0.0, 0.0}, .acceleration = {1.0, 0.0, 0.0}},
        {.w = 1.0, .position = {2.0, 0.0, 0.0}, .acceleration = {-1.0, 0.0, 0.0}}
    };
    Simulation simulation = from_bodies(bodies, 2);
    Simulation together = from_bodies(bodies, 2);
    simulate(&simulation, 0.5, 1, gravity);
    bodies[0] = body_at(&simulation, 0);
    bodies[1] = body_at(&simulation, 1);
    assert_vec(bodies[0].position, (Vec3){0.0, 0.0, 0.0});
    assert_vec(bodies[1].position, (Vec3){2.0, 0.0, 0.0});
    assert_vec(bodies[0].velocity, (Vec3){0.5, 0.0, 0.0});
    assert_vec(bodies[1].velocity, (Vec3){-0.5, 0.0, 0.0});
    assert_vec(bodies[0].acceleration, (Vec3){0.25, 0.0, 0.0});
    assert_vec(bodies[1].acceleration, (Vec3){-0.25, 0.0, 0.0});

    simulate(&simulation, 0.5, 1, gravity);
    bodies[0] = body_at(&simulation, 0);
    bodies[1] = body_at(&simulation, 1);
    assert_vec(bodies[0].position, (Vec3){0.25, 0.0, 0.0});
    assert_vec(bodies[1].position, (Vec3){1.75, 0.0, 0.0});
    assert_vec(bodies[0].velocity, (Vec3){0.625, 0.0, 0.0});
    assert_vec(bodies[1].velocity, (Vec3){-0.625, 0.0, 0.0});
    assert_vec(bodies[0].acceleration, (Vec3){4.0 / 9.0, 0.0, 0.0});
    assert_vec(bodies[1].acceleration, (Vec3){-4.0 / 9.0, 0.0, 0.0});

    simulate(&together, 0.5, 2, gravity);
    assert_body_equal(body_at(&together, 0), bodies[0]);
    assert_body_equal(body_at(&together, 1), bodies[1]);
    free(simulation.storage);
    free(together.storage);
}

static const ReferenceBody snapshot[] = {
    {1.0, {0.5, 2.0, 3.5}, {6.0, 12.0, 18.0}, {10.0, 20.0, 30.0}},
    {2.0, {5.0, 6.5, 8.0}, {24.0, 30.0, 36.0}, {40.0, 50.0, 60.0}},
    {3.0, {9.5, 11.0, 12.5}, {42.0, 48.0, 54.0}, {70.0, 80.0, 90.0}}
};
static size_t calls;
static unsigned pair_visits[3][3];

static Vec3 custom_interaction(const SimulationView *state, size_t target, size_t source)
{
    assert(target != source);
    assert(state->count == 3);
    assert(calls == target * 2 + (source < target ? source : source - 1));
    ++calls;
    ++pair_visits[target][source];
    for (size_t i = 0; i < 3; ++i) {
        assert_body_equal(view_body(state, i), snapshot[i]);
    }
    return (Vec3){state->w[source], state->vy[source], state->az[source]};
}

static void test_callback_snapshot(void)
{
    ReferenceBody bodies[] = {
        {1.0, {0.0, 1.0, 2.0}, {1.0, 2.0, 3.0}, {10.0, 20.0, 30.0}},
        {2.0, {3.0, 4.0, 5.0}, {4.0, 5.0, 6.0}, {40.0, 50.0, 60.0}},
        {3.0, {6.0, 7.0, 8.0}, {7.0, 8.0, 9.0}, {70.0, 80.0, 90.0}}
    };
    Simulation simulation = from_bodies(bodies, 3);
    simulate(&simulation, 0.5, 1, custom_interaction);
    for (size_t i = 0; i < 3; ++i) {
        bodies[i] = body_at(&simulation, i);
    }
    free(simulation.storage);
    assert(calls == 3 * 2);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            assert(pair_visits[i][j] == (i == j ? 0u : 1u));
        }
    }
    assert_vec(bodies[0].acceleration, (Vec3){5.0, 78.0, 150.0});
    assert_vec(bodies[1].acceleration, (Vec3){4.0, 60.0, 120.0});
    assert_vec(bodies[2].acceleration, (Vec3){3.0, 42.0, 90.0});
}

static ReferenceBody acceleration_snapshots[4][3];
static size_t acceleration_calls;

static Vec3 acceleration_interaction(const SimulationView *state, size_t target, size_t source)
{
    size_t step = acceleration_calls / 6;
    assert(step < 4 && state->count == 3 && target != source);
    assert(acceleration_calls % 6 == target * 2 + (source < target ? source : source - 1));
    ++acceleration_calls;
    for (size_t i = 0; i < 3; ++i) {
        assert_body_equal(view_body(state, i), acceleration_snapshots[step][i]);
    }
    return (Vec3){state->ax[source], state->ay[source], state->az[source]};
}

static void test_callback_multiple_steps(void)
{
    const double dt = 0.5;
    for (size_t steps = 3; steps <= 4; ++steps) {
        ReferenceBody expected[3];
        for (size_t i = 0; i < 3; ++i) {
            expected[i] = snapshot[i];
        }
        /* Independently compute each full callback snapshot before any sums. */
        for (size_t step = 0; step < steps; ++step) {
            for (size_t i = 0; i < 3; ++i) {
                ReferenceBody *body = &expected[i];
                body->position.x += body->velocity.x * dt;
                body->position.y += body->velocity.y * dt;
                body->position.z += body->velocity.z * dt;
                body->velocity.x += body->acceleration.x * dt;
                body->velocity.y += body->acceleration.y * dt;
                body->velocity.z += body->acceleration.z * dt;
                acceleration_snapshots[step][i] = *body;
            }
            for (size_t i = 0; i < 3; ++i) {
                Vec3 sum = {0.0, 0.0, 0.0};
                for (size_t j = 0; j < 3; ++j) {
                    if (i != j) {
                        Vec3 contribution = acceleration_snapshots[step][j].acceleration;
                        sum.x += contribution.x;
                        sum.y += contribution.y;
                        sum.z += contribution.z;
                    }
                }
                expected[i].acceleration = sum;
            }
        }

        Simulation together = from_bodies(snapshot, 3);
        Simulation separate = from_bodies(snapshot, 3);
        acceleration_calls = 0;
        simulate(&together, dt, steps, acceleration_interaction);
        assert(acceleration_calls == 6 * steps);
        acceleration_calls = 0;
        for (size_t step = 0; step < steps; ++step) {
            simulate(&separate, dt, 1, acceleration_interaction);
        }
        assert(acceleration_calls == 6 * steps);
        for (size_t i = 0; i < 3; ++i) {
            assert_body_equal(body_at(&together, i), expected[i]);
            assert_body_equal(body_at(&separate, i), expected[i]);
        }
        free(together.storage);
        free(separate.storage);
    }
}

static void test_construction_and_noops(void)
{
    ReferenceBody original = {2.0, {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    Simulation simulation = from_bodies(&original, 1);
    size_t previous_calls = calls;
    simulate(&simulation, 0.1, 0, custom_interaction);
    assert_body_equal(body_at(&simulation, 0), original);
    free(simulation.storage);

    assert(simulation_create(&simulation, 0));
    assert(simulation.count == 0 && simulation.storage == NULL);
    assert(simulation.w == NULL && simulation.next_az == NULL);
    simulate(&simulation, 0.1, 10, custom_interaction);
    simulate(&simulation, 0.1, 10, gravity);
    assert(calls == previous_calls);
    free(simulation.storage);

    assert(!simulation_create(&simulation, SIZE_MAX));
    assert(simulation.count == 0 && simulation.storage == NULL);
    assert(!simulation_create(&simulation, SIZE_MAX / (13 * sizeof(double)) + 1));
    assert(simulation.count == 0 && simulation.storage == NULL);
}

static Vec3 reference_gravity(const ReferenceBody *target, const ReferenceBody *source)
{
    double dx = source->position.x - target->position.x;
    double dy = source->position.y - target->position.y;
    double dz = source->position.z - target->position.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 == 0.0) {
        return (Vec3){0.0, 0.0, 0.0};
    }
    double scale = source->w / (r2 * sqrt(r2));
    return (Vec3){scale * dx, scale * dy, scale * dz};
}

/* Force the generic callback path, independent of the unordered-pair kernel. */
static Vec3 gravity_callback(const SimulationView *state, size_t target, size_t source)
{
    return gravity(state, target, source);
}

enum { MAX_REFERENCE_BODIES = 69 };

static const size_t reference_counts[] = {
    0, 1, 2, 3, 31, 32, 33, 67, MAX_REFERENCE_BODIES
};

static void fill_reference_bodies(ReferenceBody *bodies, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        bodies[i] = (ReferenceBody){
            .w = 0.5 + 0.25 * (double)(i % 11),
            .position = {1.25 * (double)i, (double)(i % 7) - 3.0,
                         (double)(i % 5) - 2.0},
            .velocity = {0.01 * (double)(i % 3), -0.02, 0.03},
            .acceleration = {-0.04, 0.05 * (double)(i % 4), 0.06}
        };
    }
    if (count > 2) {
        /* Coincident unequal-mass bodies remain coincident through the update. */
        bodies[2].position = bodies[1].position;
        bodies[2].velocity = bodies[1].velocity;
        bodies[2].acceleration = bodies[1].acceleration;
    }
}

static void test_gravity_accelerations(void)
{
    ReferenceBody bodies[MAX_REFERENCE_BODIES];
    double output[3][MAX_REFERENCE_BODIES + 2];
    for (size_t case_index = 0; case_index < sizeof(reference_counts) / sizeof(reference_counts[0]); ++case_index) {
        const size_t count = reference_counts[case_index];
        fill_reference_bodies(bodies, count);
        Simulation simulation = from_bodies(bodies, count);
        SimulationView state = simulation_view(&simulation);
        for (size_t axis = 0; axis < 3; ++axis) {
            output[axis][0] = -12345.0;
            for (size_t i = 0; i < count; ++i) {
                output[axis][i + 1] = NAN;
            }
            output[axis][count + 1] = 67890.0;
        }
        gravity_accelerations(&state, output[0] + 1, output[1] + 1, output[2] + 1);
        for (size_t i = 0; i < count; ++i) {
            Vec3 expected = {0.0, 0.0, 0.0};
            for (size_t j = 0; j < count; ++j) {
                if (i == j) {
                    continue;
                }
                Vec3 contribution = reference_gravity(&bodies[i], &bodies[j]);
                expected.x += contribution.x;
                expected.y += contribution.y;
                expected.z += contribution.z;
            }
            assert_vec((Vec3){output[0][i + 1], output[1][i + 1], output[2][i + 1]}, expected);
            assert_body_equal(body_at(&simulation, i), bodies[i]);
            assert(isnan(simulation.next_ax[i]));
            assert(isnan(simulation.next_ay[i]));
            assert(isnan(simulation.next_az[i]));
        }
        for (size_t axis = 0; axis < 3; ++axis) {
            assert(output[axis][0] == -12345.0);
            assert(output[axis][count + 1] == 67890.0);
        }
        if (count == 0) {
            gravity_accelerations(&state, NULL, NULL, NULL);
        }
        free(simulation.storage);
    }
}

static void test_gravity_workspace_initialization(void)
{
    const ReferenceBody initial[] = {
        {.w = 1.0, .position = {0.0, -0.0, 0.0}},
        {.w = 2.0, .position = {0.0, 0.0, -0.0}},
        {.w = 4.0, .position = {2.0, 0.0, 0.0}}
    };
    for (size_t count = 1; count <= 3; ++count) {
        Simulation simulation = from_bodies(initial, count);
        SimulationView state = simulation_view(&simulation);
        gravity_accelerations(&state, simulation.next_ax, simulation.next_ay,
                              simulation.next_az);
        for (size_t i = 0; i < count; ++i) {
            double expected_x = count < 3 ? 0.0 : (i < 2 ? 1.0 : -0.75);
            assert(simulation.next_ax[i] == expected_x);
            assert(simulation.next_ay[i] == 0.0 && !signbit(simulation.next_ay[i]));
            assert(simulation.next_az[i] == 0.0 && !signbit(simulation.next_az[i]));
            assert_body_equal(body_at(&simulation, i), initial[i]);
            simulation.x[i] = 0.0;
        }
        /* All pairs now coincide. Overwrite the preceding finite results,
         * including the first target and the single-body case. */
        gravity_accelerations(&state, simulation.next_ax, simulation.next_ay,
                              simulation.next_az);
        for (size_t i = 0; i < count; ++i) {
            assert(simulation.next_ax[i] == 0.0 && !signbit(simulation.next_ax[i]));
            assert(simulation.next_ay[i] == 0.0 && !signbit(simulation.next_ay[i]));
            assert(simulation.next_az[i] == 0.0 && !signbit(simulation.next_az[i]));
        }
        free(simulation.storage);
    }
}

static void test_target_accumulation_order(void)
{
    const size_t first = 1;
    const size_t count = first + 3;
    ReferenceBody bodies[4];
    fill_reference_bodies(bodies, count);
    for (size_t i = 0; i < count; ++i) {
        bodies[i].position = (Vec3){0.0, 0.0, 0.0};
    }
    bodies[first].position.x = 1.0;
    bodies[first].w = 0x1p54;
    bodies[first + 1].position.x = -1.0;
    bodies[first + 1].w = 0x1p54;
    bodies[first + 2].position.x = 1.0;
    bodies[first + 2].w = 1.0;
    Simulation simulation = from_bodies(bodies, count);
    SimulationView state = simulation_view(&simulation);
    gravity_accelerations(&state, simulation.next_ax, simulation.next_ay, simulation.next_az);
    /* Target 0 receives (2^54 - 2^54) + 1. Combining the negative and unit
     * contributions first would lose the unit. */
    assert(simulation.next_ax[0] == 1.0);
    assert(simulation.next_ay[0] == 0.0 && simulation.next_az[0] == 0.0);
    free(simulation.storage);
}

static void test_source_accumulation_order(void)
{
    const size_t source = 4;
    const size_t count = source + 1;
    ReferenceBody bodies[5];
    const Vec3 directions[] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    for (size_t axis = 0; axis < 3; ++axis) {
        for (size_t i = 0; i < count; ++i) {
            bodies[i] = (ReferenceBody){.w = 1.0};
        }
        bodies[0].position = directions[axis];
        bodies[0].w = 0x1p54;
        bodies[2].position = (Vec3){-directions[axis].x, -directions[axis].y,
                                   -directions[axis].z};
        bodies[2].w = 0x1p54;
        bodies[3].position = directions[axis];
        Simulation simulation = from_bodies(bodies, count);
        SimulationView state = simulation_view(&simulation);
        gravity_accelerations(&state, simulation.next_ax, simulation.next_ay, simulation.next_az);
        /* The source must receive (2^54 - 2^54) + 1. Combining the two
         * contributions from targets 2 and 3 first would lose the unit. */
        assert(simulation.next_ax[source] == directions[axis].x);
        assert(simulation.next_ay[source] == directions[axis].y);
        assert(simulation.next_az[source] == directions[axis].z);
        free(simulation.storage);
    }
}

static void test_acceleration_buffer_swaps(void)
{
    const Interaction interactions[] = {gravity, gravity_callback};
    const size_t steps_per_call[] = {0, 1, 0, 1, 3, 3};
    for (size_t variant = 0; variant < 2; ++variant) {
        for (size_t count = 0; count <= 3; ++count) {
            Simulation simulation = from_bodies(snapshot, count);
            const Simulation original = simulation; /* Borrow pointers only. */
            size_t total_steps = 0;
            for (size_t call = 0; call < sizeof(steps_per_call) / sizeof(steps_per_call[0]); ++call) {
                simulate(&simulation, 0.001, steps_per_call[call], interactions[variant]);
                total_steps += steps_per_call[call];
                bool swapped = count != 0 && total_steps % 2 != 0;
                assert(simulation.ax == (swapped ? original.next_ax : original.ax));
                assert(simulation.ay == (swapped ? original.next_ay : original.ay));
                assert(simulation.az == (swapped ? original.next_az : original.az));
                assert(simulation.next_ax == (swapped ? original.ax : original.next_ax));
                assert(simulation.next_ay == (swapped ? original.ay : original.next_ay));
                assert(simulation.next_az == (swapped ? original.az : original.next_az));
                assert(simulation.storage == original.storage && simulation.count == count);
                assert(simulation.w == original.w);
                assert(simulation.x == original.x && simulation.y == original.y && simulation.z == original.z);
                assert(simulation.vx == original.vx && simulation.vy == original.vy && simulation.vz == original.vz);
                SimulationView state = simulation_view(&simulation);
                assert(state.ax == simulation.ax && state.ay == simulation.ay && state.az == simulation.az);
            }
            free(simulation.storage);
        }
    }
}

static void check_reference_simulation(const ReferenceBody *initial, size_t count, size_t steps)
{
    const double dt = 0.001;
    ReferenceBody *reference = malloc((count + 1) * sizeof(*reference));
    Vec3 *accelerations = malloc((count + 1) * sizeof(*accelerations));
    assert(reference && accelerations);
    for (size_t i = 0; i < count; ++i) {
        reference[i] = initial[i];
    }
    Simulation paired = from_bodies(initial, count);
    Simulation ordered = from_bodies(initial, count);
    Simulation batched = from_bodies(initial, count);
    const Simulation original = paired; /* Borrow pointers only. */
    for (size_t step = 0; step < steps; ++step) {
        for (size_t i = 0; i < count; ++i) {
            ReferenceBody *body = &reference[i];
            body->position.x += body->velocity.x * dt;
            body->position.y += body->velocity.y * dt;
            body->position.z += body->velocity.z * dt;
            body->velocity.x += body->acceleration.x * dt;
            body->velocity.y += body->acceleration.y * dt;
            body->velocity.z += body->acceleration.z * dt;
        }
        for (size_t i = 0; i < count; ++i) {
            accelerations[i] = (Vec3){0.0, 0.0, 0.0};
            for (size_t j = 0; j < count; ++j) {
                if (i == j) {
                    continue;
                }
                Vec3 contribution = reference_gravity(&reference[i], &reference[j]);
                accelerations[i].x += contribution.x;
                accelerations[i].y += contribution.y;
                accelerations[i].z += contribution.z;
            }
        }
        for (size_t i = 0; i < count; ++i) {
            reference[i].acceleration = accelerations[i];
        }
        simulate(&paired, dt, 1, gravity);
        simulate(&ordered, dt, 1, gravity_callback);
        bool swapped = count != 0 && (step + 1) % 2 != 0;
        assert(paired.ax == (swapped ? original.next_ax : original.ax));
        assert(paired.ay == (swapped ? original.next_ay : original.ay));
        assert(paired.az == (swapped ? original.next_az : original.az));
        for (size_t i = 0; i < count; ++i) {
            ReferenceBody actual[] = {body_at(&paired, i), body_at(&ordered, i)};
            for (size_t variant = 0; variant < 2; ++variant) {
                assert(actual[variant].w == reference[i].w);
                assert_vec(actual[variant].position, reference[i].position);
                assert_vec(actual[variant].velocity, reference[i].velocity);
                assert_vec(actual[variant].acceleration, reference[i].acceleration);
            }
        }
    }
    simulate(&batched, dt, steps, gravity);
    for (size_t i = 0; i < count; ++i) {
        assert_body_equal(body_at(&batched, i), body_at(&paired, i));
    }
    free(paired.storage);
    free(ordered.storage);
    free(batched.storage);
    free(reference);
    free(accelerations);
}

static void test_reference_simulation(void)
{
    const ReferenceBody initial[] = {
        {1.5, {0.0, 1.0, 2.0}, {0.1, -0.2, 0.3}, {0.2, -0.3, 0.4}},
        {7.0, {0.0, 1.0, 2.0}, {0.1, -0.2, 0.3}, {0.2, -0.3, 0.4}},
        {0.5, {-2.0, 3.0, 4.0}, {-0.2, 0.4, 0.1}, {-0.1, 0.2, 0.3}},
        {3.0, {4.0, -3.0, 1.0}, {0.3, 0.1, -0.2}, {0.3, 0.1, -0.4}},
        {2.5, {-1.0, -2.0, -3.0}, {-0.1, -0.3, 0.2}, {-0.2, 0.4, 0.1}}
    };
    for (size_t count = 1; count <= 5; ++count) {
        check_reference_simulation(initial, count, 12);
    }
}

static void test_reference_simulation_sizes(void)
{
    ReferenceBody bodies[MAX_REFERENCE_BODIES];
    for (size_t case_index = 0; case_index < sizeof(reference_counts) / sizeof(reference_counts[0]); ++case_index) {
        const size_t count = reference_counts[case_index];
        fill_reference_bodies(bodies, count);
        check_reference_simulation(bodies, count, 3);
        check_reference_simulation(bodies, count, 4);
    }
}

static void test_coincident_bodies(void)
{
    const size_t count = MAX_REFERENCE_BODIES;
    const size_t source = count - 1;
    ReferenceBody bodies[MAX_REFERENCE_BODIES];
    for (size_t coincident = 0; coincident < 3; ++coincident) {
        fill_reference_bodies(bodies, count);
        if (coincident == 2) {
            bodies[1].position = bodies[0].position;
            bodies[1].velocity = bodies[0].velocity;
            bodies[1].acceleration = bodies[0].acceleration;
        }
        /* The last body coincides with the first, the second, or both. */
        const size_t target = coincident == 1 ? 1 : 0;
        bodies[source].position = bodies[target].position;
        bodies[source].velocity = bodies[target].velocity;
        bodies[source].acceleration = bodies[target].acceleration;
        check_reference_simulation(bodies, count, 3);
    }
}

static FILE *text_file(const char *text)
{
    FILE *file = tmpfile();
    assert(file);
    assert(fputs(text, file) >= 0);
    rewind(file);
    return file;
}

static void test_file_roundtrip(void)
{
    ReferenceBody original[] = {
        {0.1, {1e-180, -1e180, 3.141592653589793}, {0.0, -0.5, 0.75}, {2.0, 3.0, -4.0}},
        {2.5, {nextafter(1.0, 2.0), 8.0, 9.0}, {10.0, 11.0, 12.0}, {13.0, 14.0, 15.0}}
    };
    FILE *file = tmpfile();
    assert(file);
    Simulation simulation = from_bodies(original, 2);
    assert(write_bodies(file, &simulation));
    free(simulation.storage);
    rewind(file);
    Simulation loaded;
    assert(read_bodies(file, &loaded));
    assert(loaded.count == 2);
    assert_body_equal(body_at(&loaded, 0), original[0]);
    assert_body_equal(body_at(&loaded, 1), original[1]);
    free(loaded.storage);
    assert(fclose(file) == 0);

    file = text_file(" \n0\t\n");
    assert(read_bodies(file, &loaded));
    assert(loaded.count == 0 && loaded.storage == NULL);
    assert(fclose(file) == 0);

    file = text_file("1\n1 2 3 4 5 6 7 8 9 10"); /* No trailing newline. */
    assert(read_bodies(file, &loaded));
    assert(loaded.count == 1);
    assert_body_equal(body_at(&loaded, 0), (ReferenceBody){1.0, {2.0, 3.0, 4.0}, {5.0, 6.0, 7.0}, {8.0, 9.0, 10.0}});
    free(loaded.storage);
    assert(fclose(file) == 0);

    file = fopen("tests/data/two_body.txt", "r");
    assert(file);
    assert(read_bodies(file, &loaded));
    assert(loaded.count == 2);
    SimulationView state = simulation_view(&loaded);
    assert_vec(body_at(&loaded, 0).acceleration, gravity(&state, 0, 1));
    assert_vec(body_at(&loaded, 1).acceleration, gravity(&state, 1, 0));
    free(loaded.storage);
    assert(fclose(file) == 0);
}

static void test_simulated_state_roundtrip(void)
{
    for (size_t steps = 1; steps <= 4; ++steps) {
        Simulation simulation = from_bodies(snapshot, 3);
        simulate(&simulation, 0.001, steps, gravity);
        for (size_t i = 0; i < simulation.count; ++i) {
            simulation.next_ax[i] = NAN;
            simulation.next_ay[i] = NAN;
            simulation.next_az[i] = NAN;
        }
        FILE *file = tmpfile();
        assert(file);
        assert(write_bodies(file, &simulation));
        rewind(file);
        Simulation loaded;
        assert(read_bodies(file, &loaded));
        assert(loaded.count == simulation.count);
        for (size_t i = 0; i < simulation.count; ++i) {
            assert_body_equal(body_at(&loaded, i), body_at(&simulation, i));
        }
        /* Restarting from saved state must also work with either buffer layout. */
        simulate(&simulation, 0.001, 1, gravity);
        simulate(&loaded, 0.001, 1, gravity);
        for (size_t i = 0; i < simulation.count; ++i) {
            assert_body_equal(body_at(&loaded, i), body_at(&simulation, i));
        }
        free(simulation.storage);
        free(loaded.storage);
        assert(fclose(file) == 0);
    }
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
        Simulation simulation;
        assert(!read_bodies(file, &simulation));
        assert(simulation.storage == NULL && simulation.count == 0);
        assert(fclose(file) == 0);
    }

    char maximum[32];
    size_t count;
    assert(snprintf(maximum, sizeof(maximum), "%zu", SIZE_MAX) > 0);
    assert(parse_size(maximum, &count) && count == SIZE_MAX);
    FILE *file = text_file(maximum); /* Allocation size overflow. */
    Simulation simulation;
    assert(!read_bodies(file, &simulation));
    assert(simulation.storage == NULL && simulation.count == 0);
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

static void check_tree(const SimulationView *state, const ApproximateWorkspace *workspace)
{
    size_t leaves = 0;
    assert(workspace->node_count <= workspace->node_capacity);
    for (size_t n = 0; n < workspace->node_count; ++n) {
        const OctreeNode *node = &workspace->nodes[n];
        size_t next = node->begin;
        bool internal = false;
        for (size_t c = 0; c < 8; ++c) {
            if (node->children[c] != SIZE_MAX) {
                internal = true;
                assert(node->children[c] > n && node->children[c] < workspace->node_count);
                const OctreeNode *child = &workspace->nodes[node->children[c]];
                assert(child->begin == next && child->end > child->begin);
                next = child->end;
            }
        }
        if (internal) {
            assert(next == node->end);
        } else {
            leaves += node->end - node->begin;
        }
        for (size_t p = node->begin; p < node->end; ++p) {
            size_t i = workspace->indices[p];
            assert(i < state->count && workspace->rank[i] == p);
            assert(fabs(state->x[i] - node->center.x) <= node->half_width);
            assert(fabs(state->y[i] - node->center.y) <= node->half_width);
            assert(fabs(state->z[i] - node->center.z) <= node->half_width);
        }
    }
    assert(leaves == state->count);
}

static void test_octree_build_and_rebuild(void)
{
    Simulation simulation;
    assert(simulation_create(&simulation, 27));
    for (size_t i = 0; i < simulation.count; ++i) {
        simulation.w[i] = 1.0;
        simulation.x[i] = (double)(i % 3) - 1.0;
        simulation.y[i] = (double)((i / 3) % 3) - 1.0;
        simulation.z[i] = (double)(i / 9) - 1.0;
    }
    SimulationView state = simulation_view(&simulation);
    ApproximateConfig config = approximate_default_config();
    config.bucket_capacity = 1;
    ApproximateWorkspace workspace;
    assert(approximate_workspace_create(&workspace, state.count, config));
    gravity_approximate(&state, simulation.next_ax, simulation.next_ay,
                        simulation.next_az, &workspace);
    assert(workspace.node_count > 1 && !workspace.bounds_fallback);
    check_tree(&state, &workspace);
    /* Rebuild a coincident cloud at a different location using the same pool. */
    OctreeNode *pool = workspace.nodes;
    for (size_t i = 0; i < state.count; ++i) {
        simulation.x[i] = simulation.y[i] = simulation.z[i] = 100.0;
    }
    gravity_approximate(&state, simulation.next_ax, simulation.next_ay,
                        simulation.next_az, &workspace);
    assert(workspace.nodes == pool && workspace.nodes[0].center.x == 100.0);
    check_tree(&state, &workspace);
    for (size_t i = 0; i < state.count; ++i) {
        assert(simulation.next_ax[i] == 0.0 && simulation.next_ay[i] == 0.0
               && simulation.next_az[i] == 0.0);
    }
    /* Unrepresentable full width selects direct fallback and clears old nodes. */
    simulation.x[0] = -DBL_MAX;
    simulation.x[1] = DBL_MAX;
    gravity_approximate(&state, simulation.next_ax, simulation.next_ay,
                        simulation.next_az, &workspace);
    assert(workspace.bounds_fallback && workspace.node_count == 0);
    approximate_workspace_destroy(&workspace);
    for (size_t i = 0; i < state.count; ++i) {
        simulation.x[i] = (double)i;
    }
    config.node_capacity = 1;
    assert(approximate_workspace_create(&workspace, state.count, config));
    gravity_approximate(&state, simulation.next_ax, simulation.next_ay,
                        simulation.next_az, &workspace);
    assert(workspace.node_count == 1 && workspace.capacity_limited_leaves == 1);
    check_tree(&state, &workspace);
    approximate_workspace_destroy(&workspace);
    config.max_depth = 0;
    assert(approximate_workspace_create(&workspace, state.count, config));
    gravity_approximate(&state, simulation.next_ax, simulation.next_ay,
                        simulation.next_az, &workspace);
    assert(workspace.node_count == 1 && workspace.capacity_limited_leaves == 0);
    check_tree(&state, &workspace);
    approximate_workspace_destroy(&workspace);
    free(simulation.storage);

    config = approximate_default_config();
    assert(approximate_workspace_create(&workspace, 0, config));
    SimulationView empty = {0};
    gravity_approximate(&empty, NULL, NULL, NULL, &workspace);
    assert(workspace.node_count == 0 && !workspace.bounds_fallback);
    approximate_workspace_destroy(&workspace);
    approximate_workspace_destroy(&workspace);
    assert(!approximate_workspace_create(&workspace, SIZE_MAX, config));
    config.max_depth = APPROXIMATE_MAX_DEPTH + 1;
    assert(!approximate_workspace_create(&workspace, 1, config));
    config = approximate_default_config();
    config.theta = NAN;
    assert(!approximate_workspace_create(&workspace, 1, config));
    assert(!workspace.nodes && !workspace.indices && !workspace.rank);
}

static void test_octree_mass_and_traversal(void)
{
    const ReferenceBody bodies[] = {
        {1, {0, 0, 0}, {0}, {0}},
        {2, {10, 0, 0}, {0}, {0}},
        {3, {11, 0, 0}, {0}, {0}}
    };
    double direct_ax[3], direct_ay[3], direct_az[3];
    Simulation simulation = from_bodies(bodies, 3);
    SimulationView state = simulation_view(&simulation);
    ApproximateConfig config = approximate_default_config();
    config.bucket_capacity = 1;
    config.theta = 100.0; /* Still must open every node containing the target. */
    ApproximateWorkspace workspace;
    assert(approximate_workspace_create(&workspace, 3, config));
    gravity_approximate(&state, simulation.next_ax, simulation.next_ay,
                        simulation.next_az, &workspace);
    assert(workspace.nodes[0].mass == 6.0);
    assert_vec(workspace.nodes[0].center_of_mass, (Vec3){53.0 / 6.0, 0, 0});
    assert(workspace.accepted_nodes > 0);
    assert(fabs(simulation.next_ax[0] - 5.0 / (10.6 * 10.6)) < 1e-14);
    /* Targets inside the cluster must see the other particle, not themselves. */
    assert(fabs(simulation.next_ax[1] - (3.0 - 0.01)) < 1e-14);
    assert(fabs(simulation.next_ax[2] - (-2.0 - 1.0 / 121.0)) < 1e-14);
    for (size_t i = 0; i < 3; ++i) {
        assert(simulation.next_ay[i] == 0.0 && simulation.next_az[i] == 0.0);
    }
    /* Full opening exercises the tree but accepts no internal aggregates. */
    workspace.config.theta = 0.0;
    gravity_direct(&state, direct_ax, direct_ay, direct_az, NULL);
    gravity_approximate(&state, simulation.next_ax, simulation.next_ay,
                        simulation.next_az, &workspace);
    assert(workspace.accepted_nodes == 0);
    for (size_t i = 0; i < 3; ++i) {
        assert_vec((Vec3){simulation.next_ax[i], simulation.next_ay[i], simulation.next_az[i]},
                   (Vec3){direct_ax[i], direct_ay[i], direct_az[i]});
        simulation.x[i] += 100.0;
        simulation.y[i] += 20.0;
        simulation.z[i] -= 30.0;
    }
    gravity_approximate(&state, simulation.next_ax, simulation.next_ay,
                        simulation.next_az, &workspace);
    assert_vec(workspace.nodes[0].center_of_mass, (Vec3){100.0 + 53.0 / 6.0, 20, -30});
    /* Finite masses may have an unrepresentable aggregate. Direct fallback
     * still produces representable forces for this separated cloud. */
    simulation.w[0] = simulation.w[1] = simulation.w[2] = DBL_MAX / 2.0;
    gravity_direct(&state, direct_ax, direct_ay, direct_az, NULL);
    gravity_approximate(&state, simulation.next_ax, simulation.next_ay,
                        simulation.next_az, &workspace);
    assert(workspace.aggregate_fallback && !workspace.bounds_fallback);
    assert(workspace.accepted_nodes == 0);
    for (size_t i = 0; i < 3; ++i) {
        assert(isfinite(simulation.next_ax[i]));
        assert(simulation.next_ax[i] == direct_ax[i]);
    }
    /* A fresh build clears failed aggregates and diagnostics. */
    simulation.w[0] = 1; simulation.w[1] = 2; simulation.w[2] = 3;
    gravity_approximate(&state, simulation.next_ax, simulation.next_ay,
                        simulation.next_az, &workspace);
    assert(!workspace.aggregate_fallback && workspace.nodes[0].mass == 6.0);
    /* Mass-times-position would overflow, although mass and center are finite. */
    for (size_t i = 0; i < 3; ++i) {
        simulation.w[i] = 1e200;
        simulation.x[i] = simulation.y[i] = simulation.z[i] = 1e150;
    }
    gravity_approximate(&state, simulation.next_ax, simulation.next_ay,
                        simulation.next_az, &workspace);
    assert(!workspace.aggregate_fallback);
    assert(isfinite(workspace.nodes[0].mass));
    assert(workspace.nodes[0].center_of_mass.x == 1e150);
    approximate_workspace_destroy(&workspace);
    free(simulation.storage);
}

static void test_octree_full_opening_cloud(void)
{
    double direct_ax[64], direct_ay[64], direct_az[64];
    Simulation simulation;
    assert(simulation_create(&simulation, 64));
    for (size_t i = 0; i < 64; ++i) {
        simulation.w[i] = 1.0 + (double)(i % 7);
        simulation.x[i] = (double)(i % 4);
        simulation.y[i] = (double)((i / 4) % 4);
        simulation.z[i] = (double)(i / 16);
    }
    SimulationView state = simulation_view(&simulation);
    ApproximateConfig config = approximate_default_config();
    config.theta = 0;
    config.bucket_capacity = 1;
    ApproximateWorkspace workspace;
    assert(approximate_workspace_create(&workspace, state.count, config));
    gravity_direct(&state, direct_ax, direct_ay, direct_az, NULL);
    gravity_approximate(&state, simulation.next_ax, simulation.next_ay,
                        simulation.next_az, &workspace);
    assert(workspace.node_count > 8 && workspace.accepted_nodes == 0);
    for (size_t i = 0; i < state.count; ++i) {
        assert_vec((Vec3){simulation.next_ax[i], simulation.next_ay[i], simulation.next_az[i]},
                   (Vec3){direct_ax[i], direct_ay[i], direct_az[i]});
    }
    approximate_workspace_destroy(&workspace);
    free(simulation.storage);
}

static void test_solver_equivalence(void)
{
    const ReferenceBody bodies[] = {
        {2, {1, 2, 3}, {0.1, 0.2, 0.3}, {1, -2, 3}},
        {3, {1, 2, 3}, {0.1, 0.2, 0.3}, {1, -2, 3}},
        {4, {-2, 1, 0}, {-0.2, 0.1, 0.4}, {-1, 3, 2}}
    };
    assert(gravity_solver_by_name("direct") == gravity_direct);
    assert(gravity_solver_by_name("approximate") == gravity_approximate);
    assert(gravity_solver_by_name("unknown") == NULL);
    for (size_t count = 0; count <= 3; ++count) {
        for (size_t steps = 0; steps <= 4; ++steps) {
            Simulation legacy = from_bodies(bodies, count);
            Simulation direct = from_bodies(bodies, count);
            Simulation approximate = from_bodies(bodies, count);
            simulate(&legacy, 0.01, steps, gravity);
            simulate_with_solver(&direct, 0.01, steps, gravity_direct, NULL);
            ApproximateWorkspace workspace;
            assert(approximate_workspace_create(&workspace, count, approximate_default_config()));
            simulate_with_solver(&approximate, 0.01, steps, gravity_approximate, &workspace);
            approximate_workspace_destroy(&workspace);
            for (size_t i = 0; i < count; ++i) {
                assert_body_equal(body_at(&legacy, i), body_at(&direct, i));
                ReferenceBody actual = body_at(&approximate, i);
                ReferenceBody expected = body_at(&direct, i);
                assert(actual.w == expected.w);
                assert_vec(actual.position, expected.position);
                assert_vec(actual.velocity, expected.velocity);
                assert_vec(actual.acceleration, expected.acceleration);
            }
            free(legacy.storage);
            free(direct.storage);
            free(approximate.storage);
        }
    }
}

static void check_solver_snapshot(const SimulationView *state, double *ax,
                                  double *ay, double *az, void *context)
{
    size_t *calls = context;
    assert(state->count == 1);
    assert(state->x[0] == (*calls == 0 ? 2.0 : 3.125));
    assert(state->vx[0] == (*calls == 0 ? 4.5 : 7.0));
    assert(state->ax[0] == (*calls == 0 ? 2.0 : 10.0));
    ax[0] = 10.0 + *calls;
    ay[0] = 20.0 + *calls;
    az[0] = 30.0 + *calls;
    ++*calls;
}

static void test_solver_context_and_order(void)
{
    ReferenceBody body = {1, {1, 2, 3}, {4, 5, 6}, {2, 4, 6}};
    Simulation simulation = from_bodies(&body, 1);
    size_t calls = 0;
    double *original_ax = simulation.ax;
    simulate_with_solver(&simulation, 0.25, 0, check_solver_snapshot, &calls);
    assert(calls == 0 && simulation.ax == original_ax);
    simulate_with_solver(&simulation, 0.25, 2, check_solver_snapshot, &calls);
    assert(calls == 2 && simulation.ax == original_ax);
    assert_vec(body_at(&simulation, 0).acceleration, (Vec3){11, 21, 31});
    free(simulation.storage);
    simulation = from_bodies(NULL, 0);
    simulate_with_solver(&simulation, 0.25, 2, check_solver_snapshot, &calls);
    assert(calls == 2);
    free(simulation.storage);
}

int main(void)
{
    test_octree_build_and_rebuild();
    test_octree_mass_and_traversal();
    test_octree_full_opening_cloud();
    test_solver_equivalence();
    test_solver_context_and_order();
    test_gravity();
    test_gravity_accelerations();
    test_gravity_workspace_initialization();
    test_target_accumulation_order();
    test_source_accumulation_order();
    test_euler();
    test_input_acceleration();
    test_callback_snapshot();
    test_callback_multiple_steps();
    test_construction_and_noops();
    test_acceleration_buffer_swaps();
    test_reference_simulation();
    test_reference_simulation_sizes();
    test_coincident_bodies();
    test_file_roundtrip();
    test_simulated_state_roundtrip();
    test_invalid_input();
    puts("All unit tests passed.");
    return 0;
}
