#include "nbody.h"
#include "gravity.h"

/* Acceleration buffers swap roles each step; the arrays remain disjoint. */
static void advance_state(size_t count, double dt,
                          double *restrict x, double *restrict y, double *restrict z,
                          double *restrict vx, double *restrict vy, double *restrict vz,
                          const double *restrict ax, const double *restrict ay,
                          const double *restrict az)
{
    for (size_t i = 0; i < count; ++i) {
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
        z[i] += vz[i] * dt;
        vx[i] += ax[i] * dt;
        vy[i] += ay[i] * dt;
        vz[i] += az[i] * dt;
    }
}

void simulate_with_solver(Simulation *simulation, double dt, size_t steps,
                          AccelerationSolver solver, void *context)
{
    const size_t count = simulation->count;
    if (steps == 0 || count == 0) {
        return;
    }
    for (size_t step = 0; step < steps; ++step) {
        advance_state(count, dt, simulation->x, simulation->y, simulation->z,
                      simulation->vx, simulation->vy, simulation->vz,
                      simulation->ax, simulation->ay, simulation->az);

        const SimulationView state = simulation_view(simulation);
        solver(&state, simulation->next_ax, simulation->next_ay,
               simulation->next_az, context);

        /* Publish new accelerations only after every interaction has read state. */
        double *previous_ax = simulation->ax;
        double *previous_ay = simulation->ay;
        double *previous_az = simulation->az;
        simulation->ax = simulation->next_ax;
        simulation->ay = simulation->next_ay;
        simulation->az = simulation->next_az;
        simulation->next_ax = previous_ax;
        simulation->next_ay = previous_ay;
        simulation->next_az = previous_az;
    }
}

static void interaction_accelerations(const SimulationView *state, double *ax,
                                     double *ay, double *az, void *context)
{
    Interaction interaction = *(const Interaction *)context;
    for (size_t i = 0; i < state->count; ++i) {
        double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
        for (size_t j = 0; j < state->count; ++j) {
            if (i == j) {
                continue;
            }
            Vec3 contribution = interaction(state, i, j);
            sum_x += contribution.x;
            sum_y += contribution.y;
            sum_z += contribution.z;
        }
        ax[i] = sum_x;
        ay[i] = sum_y;
        az[i] = sum_z;
    }
}

void simulate(Simulation *simulation, double dt, size_t steps, Interaction interaction)
{
    if (interaction == gravity) {
        simulate_with_solver(simulation, dt, steps, gravity_direct, NULL);
    } else {
        simulate_with_solver(simulation, dt, steps, interaction_accelerations, &interaction);
    }
}
