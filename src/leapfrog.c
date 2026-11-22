#include "state_internal.h"
#include "integrator_internal.h"
#include <string.h>
#ifdef NBODY_BENCH_STEPS
#include "step_timing.h"
#endif

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

const char *simulate_leapfrog_with_solver(Simulation *simulation, double dt, size_t steps,
                                          AccelerationSolver solver, void *context)
{
    const size_t count = simulation->count;
    if (steps == 0 || count == 0) {
        return NULL;
    }
    const double half_dt = 0.5 * dt;
    /* Mutable public views make accelerations from any earlier call untrusted. */
    const char *error = refresh_acceleration(simulation, solver, context);
    if (error) {
        return error;
    }
#ifdef NBODY_BENCH_STEPS
    benchmark_step_boundary(0);
#endif
    for (size_t step = 0; step < steps; ++step) {
        kick_drift(count, dt, half_dt, simulation->x, simulation->y, simulation->z,
                   simulation->vx, simulation->vy, simulation->vz,
                   simulation->ax, simulation->ay, simulation->az);
        error = refresh_acceleration(simulation, solver, context);
        if (error) {
            return error;
        }
        finish_kick(count, half_dt, simulation->vx, simulation->vy, simulation->vz,
                    simulation->ax, simulation->ay, simulation->az);
#ifdef NBODY_BENCH_STEPS
        benchmark_step_boundary(step + 1);
#endif
    }
    return NULL;
}

const char *simulate_with_solver(Simulation *simulation, double dt, size_t steps,
                                 AccelerationSolver solver, void *context)
{
    return simulate_leapfrog_with_solver(simulation, dt, steps, solver, context);
}

Integrate parse_integrator(const char *name)
{
    if (!name || strcmp(name, "leapfrog") == 0) {
        return simulate_with_solver;
    }
    if (strcmp(name, "euler") == 0) {
        return simulate_euler_with_solver;
    }
    return NULL;
}
