#include "gravity_internal.h"

#include <math.h>

/* Newtonian acceleration in normalized units (G = 1), with w as mass.
 * Coincident positions contribute zero.
 */
static Vec3 gravity(const SimulationView *state, size_t target, size_t source)
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

Vec3 gravity_accumulate_sources(const SimulationView *state, size_t target,
                                const size_t *indices, size_t begin, size_t end,
                                Vec3 sum)
{
    for (size_t p = begin; p < end; ++p) {
        size_t source = indices[p];
        if (source == target) {
            continue;
        }
        Vec3 contribution = gravity(state, target, source);
        sum.x += contribution.x;
        sum.y += contribution.y;
        sum.z += contribution.z;
    }
    return sum;
}
