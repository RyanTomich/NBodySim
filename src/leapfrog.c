#include "state_internal.h"

static void kick_drift(size_t count, double dt, double half_dt,
                          double *restrict x, double *restrict y, double *restrict z,
                          double *restrict vx, double *restrict vy, double *restrict vz,
                          const double *restrict ax, const double *restrict ay,
                          const double *restrict az)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < count; ++i) {
        double vx_half = vx[i] + half_dt * ax[i];
        double vy_half = vy[i] + half_dt * ay[i];
        double vz_half = vz[i] + half_dt * az[i];
        vx[i] = vx_half;
        vy[i] = vy_half;
        vz[i] = vz_half;
        x[i] += dt * vx_half;
        y[i] += dt * vy_half;
        z[i] += dt * vz_half;
    }
}

static void finish_kick(size_t count, double half_dt,
                        double *restrict vx, double *restrict vy, double *restrict vz,
                        const double *restrict ax, const double *restrict ay,
                        const double *restrict az)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < count; ++i) {
        vx[i] += half_dt * ax[i];
        vy[i] += half_dt * ay[i];
        vz[i] += half_dt * az[i];
    }
}

static void refresh_acceleration(Simulation *simulation, AccelerationSolver solver, void *context)
{
    const SimulationView state = simulation_view(simulation);
    solver(&state, simulation->next_ax, simulation->next_ay,
           simulation->next_az, context);

    /* Publish only after every interaction has read the whole input state. */
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

void simulate_with_solver(Simulation *simulation, double dt, size_t steps,
                          AccelerationSolver solver, void *context)
{
    const size_t count = simulation->count;
    if (steps == 0 || count == 0) {
        return;
    }
    const double half_dt = 0.5 * dt;
    /* Mutable public views make accelerations from any earlier call untrusted. */
    refresh_acceleration(simulation, solver, context);
    for (size_t step = 0; step < steps; ++step) {
        kick_drift(count, dt, half_dt, simulation->x, simulation->y, simulation->z,
                   simulation->vx, simulation->vy, simulation->vz,
                   simulation->ax, simulation->ay, simulation->az);
        refresh_acceleration(simulation, solver, context);
        finish_kick(count, half_dt, simulation->vx, simulation->vy, simulation->vz,
                    simulation->ax, simulation->ay, simulation->az);
    }
}
