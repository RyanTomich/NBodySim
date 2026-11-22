#include "state_internal.h"
#ifdef NBODY_BENCH_STEPS
#include "step_timing.h"
#endif

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

const char *simulate_euler_with_solver(Simulation *simulation, double dt, size_t steps,
                                       AccelerationSolver solver, void *context)
{
    const size_t count = simulation->count;
    if (steps == 0 || count == 0) {
        return NULL;
    }
#ifdef NBODY_BENCH_STEPS
    benchmark_step_boundary(0);
#endif
    for (size_t step = 0; step < steps; ++step) {
        advance_state(count, dt, simulation->x, simulation->y, simulation->z,
                      simulation->vx, simulation->vy, simulation->vz,
                      simulation->ax, simulation->ay, simulation->az);

        const char *error = refresh_acceleration(simulation, solver, context);
        if (error) {
            return error;
        }
#ifdef NBODY_BENCH_STEPS
        benchmark_step_boundary(step + 1);
#endif
    }
    return NULL;
}
