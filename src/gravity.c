#include "gravity.h"

#include <math.h>
#include <string.h>

Vec3 gravity(const SimulationView *state, size_t target, size_t source)
{
    Vec3 displacement = {
        state->x[source] - state->x[target],
        state->y[source] - state->y[target],
        state->z[source] - state->z[target]
    };
    double r2 = displacement.x * displacement.x
              + displacement.y * displacement.y
              + displacement.z * displacement.z;
    if (r2 == 0.0) {
        return (Vec3){0.0, 0.0, 0.0};
    }
    double scale = state->w[source] / (r2 * sqrt(r2));
    return (Vec3){
        scale * displacement.x,
        scale * displacement.y,
        scale * displacement.z
    };
}

void gravity_accelerations(const SimulationView *state, double *restrict next_ax,
                           double *restrict next_ay, double *restrict next_az)
{
    const size_t count = state->count;
    if (count == 0) {
        return;
    }
    const double *restrict x = state->x;
    const double *restrict y = state->y;
    const double *restrict z = state->z;
    const double *restrict mass = state->w;
    /* The first target initializes every output with its first contribution,
     * without reading the old workspace or making a separate clearing pass. */
    const double x_0 = x[0];
    const double y_0 = y[0];
    const double z_0 = z[0];
    const double mass_0 = mass[0];
    double first_ax = 0.0;
    double first_ay = 0.0;
    double first_az = 0.0;
    for (size_t j = 1; j < count; ++j) {
        Vec3 displacement = {x[j] - x_0, y[j] - y_0, z[j] - z_0};
        double r2 = displacement.x * displacement.x
                  + displacement.y * displacement.y
                  + displacement.z * displacement.z;
        if (r2 == 0.0) {
            next_ax[j] = 0.0;
            next_ay[j] = 0.0;
            next_az[j] = 0.0;
            continue;
        }
        double inv_r3 = 1.0 / (r2 * sqrt(r2));
        double scale_i = mass_0 * inv_r3;
        double scale_j = mass[j] * inv_r3;
        first_ax += scale_j * displacement.x;
        first_ay += scale_j * displacement.y;
        first_az += scale_j * displacement.z;
        /* Subtract from positive zero, as in the remaining targets' updates. */
        next_ax[j] = 0.0 - scale_i * displacement.x;
        next_ay[j] = 0.0 - scale_i * displacement.y;
        next_az[j] = 0.0 - scale_i * displacement.z;
    }
    next_ax[0] = first_ax;
    next_ay[0] = first_ay;
    next_az[0] = first_az;

    for (size_t i = 1; i < count; ++i) {
        const double x_i = x[i];
        const double y_i = y[i];
        const double z_i = z[i];
        const double mass_i = mass[i];
        /* Continue the contributions already received from lower indices. */
        double ax = next_ax[i];
        double ay = next_ay[i];
        double az = next_az[i];
        for (size_t j = i + 1; j < count; ++j) {
            Vec3 displacement = {
                x[j] - x_i,
                y[j] - y_i,
                z[j] - z_i
            };
            double r2 = displacement.x * displacement.x
                      + displacement.y * displacement.y
                      + displacement.z * displacement.z;
            if (r2 == 0.0) {
                continue;
            }
            double inv_r3 = 1.0 / (r2 * sqrt(r2));
            double scale_i = mass_i * inv_r3;
            double scale_j = mass[j] * inv_r3;
            ax += scale_j * displacement.x;
            ay += scale_j * displacement.y;
            az += scale_j * displacement.z;
            next_ax[j] -= scale_i * displacement.x;
            next_ay[j] -= scale_i * displacement.y;
            next_az[j] -= scale_i * displacement.z;
        }
        next_ax[i] = ax;
        next_ay[i] = ay;
        next_az[i] = az;
    }
}

void gravity_direct(const SimulationView *state, double *ax, double *ay,
                    double *az, void *context)
{
    (void)context;
    gravity_accelerations(state, ax, ay, az);
}

AccelerationSolver gravity_solver_by_name(const char *name)
{
    if (strcmp(name, "direct") == 0) {
        return gravity_direct;
    }
    if (strcmp(name, "approximate") == 0) {
        return gravity_approximate;
    }
    return NULL;
}
