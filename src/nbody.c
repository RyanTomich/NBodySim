#include "nbody.h"
#include "gravity.h"

/* All arrays are disjoint slices of storage. Rebind them on each step because
 * acceleration arrays change roles after the interactions finish. */
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

void simulate(Simulation *simulation, double dt, size_t steps, Interaction f)
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
        if (f == gravity) {
            gravity_accelerations(&state, simulation->next_ax,
                                  simulation->next_ay, simulation->next_az);
        } else {
            for (size_t i = 0; i < count; ++i) {
                double ax = 0.0;
                double ay = 0.0;
                double az = 0.0;
                for (size_t j = 0; j < count; ++j) {
                    if (i == j) {
                        continue;
                    }
                    Vec3 contribution = f(&state, i, j);
                    ax += contribution.x;
                    ay += contribution.y;
                    az += contribution.z;
                }
                simulation->next_ax[i] = ax;
                simulation->next_ay[i] = ay;
                simulation->next_az[i] = az;
            }
        }

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
