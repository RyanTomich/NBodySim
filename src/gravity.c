#include "gravity.h"

#include <math.h>

Vec3 gravity(const Body *target, const Body *source)
{
    Vec3 displacement = {
        source->position.x - target->position.x,
        source->position.y - target->position.y,
        source->position.z - target->position.z
    };
    double r2 = displacement.x * displacement.x
              + displacement.y * displacement.y
              + displacement.z * displacement.z;
    if (r2 == 0.0) {
        return (Vec3){0.0, 0.0, 0.0};
    }
    double scale = source->w / (r2 * sqrt(r2));
    return (Vec3){
        scale * displacement.x,
        scale * displacement.y,
        scale * displacement.z
    };
}

void gravity_accelerations(const Body *bodies, size_t count, Vec3 *accelerations)
{
    for (size_t i = 0; i < count; ++i) {
        accelerations[i] = (Vec3){0.0, 0.0, 0.0};
    }

    for (size_t i = 0; i < count; ++i) {
        const Vec3 position_i = bodies[i].position;
        const double mass_i = bodies[i].w;
        /* Keep contributions from pairs whose first index was less than i. */
        double ax = accelerations[i].x;
        double ay = accelerations[i].y;
        double az = accelerations[i].z;
        for (size_t j = i + 1; j < count; ++j) {
            Vec3 displacement = {
                bodies[j].position.x - position_i.x,
                bodies[j].position.y - position_i.y,
                bodies[j].position.z - position_i.z
            };
            double r2 = displacement.x * displacement.x
                      + displacement.y * displacement.y
                      + displacement.z * displacement.z;
            if (r2 == 0.0) {
                continue;
            }
            double inv_r3 = 1.0 / (r2 * sqrt(r2));
            double scale_i = mass_i * inv_r3;
            double scale_j = bodies[j].w * inv_r3;
            ax += scale_j * displacement.x;
            ay += scale_j * displacement.y;
            az += scale_j * displacement.z;
            accelerations[j].x -= scale_i * displacement.x;
            accelerations[j].y -= scale_i * displacement.y;
            accelerations[j].z -= scale_i * displacement.z;
        }
        accelerations[i] = (Vec3){ax, ay, az};
    }
}
