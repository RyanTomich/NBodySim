#include "state_internal.h"

static void publish_acceleration(Simulation *s)
{
    double *old = s->ax; s->ax = s->next_ax; s->next_ax = old;
    old = s->ay; s->ay = s->next_ay; s->next_ay = old;
    old = s->az; s->az = s->next_az; s->next_az = old;
}

void simulate_leapfrog_with_solver(Simulation *s, double dt, size_t steps,
                                   AccelerationSolver solver, void *context)
{
    if (!steps || !s->count) return;
    /* Compute each backend's own starting acceleration, then reuse it between
     * kicks. Supplied accelerations are not assumed consistent with the solver. */
    SimulationView view = simulation_view(s);
    solver(&view, s->next_ax, s->next_ay, s->next_az, context);
    publish_acceleration(s);
    for (size_t step = 0; step < steps; ++step) {
        for (size_t i = 0; i < s->count; ++i) {
            s->vx[i] += 0.5 * dt * s->ax[i];
            s->vy[i] += 0.5 * dt * s->ay[i];
            s->vz[i] += 0.5 * dt * s->az[i];
            s->x[i] += dt * s->vx[i];
            s->y[i] += dt * s->vy[i];
            s->z[i] += dt * s->vz[i];
        }
        view = simulation_view(s);
        solver(&view, s->next_ax, s->next_ay, s->next_az, context);
        publish_acceleration(s);
        for (size_t i = 0; i < s->count; ++i) {
            s->vx[i] += 0.5 * dt * s->ax[i];
            s->vy[i] += 0.5 * dt * s->ay[i];
            s->vz[i] += 0.5 * dt * s->az[i];
        }
    }
}
