#include "nbody.h"
#include "gravity.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

bool simulate(Body *bodies, size_t count, double dt, size_t steps, Interaction f)
{
    if ((!bodies && count != 0) || !f || !isfinite(dt) || dt <= 0.0) {
        return false;
    }
    if (steps == 0 || count == 0) {
        return true;
    }
    if (count > SIZE_MAX / sizeof(Vec3)) {
        return false;
    }
    Vec3 *accelerations = malloc(count * sizeof(*accelerations));
    if (!accelerations) {
        return false;
    }

    for (size_t step = 0; step < steps; ++step) {
        for (size_t i = 0; i < count; ++i) {
            Body *body = &bodies[i];
            body->position.x += body->velocity.x * dt;
            body->position.y += body->velocity.y * dt;
            body->position.z += body->velocity.z * dt;
            body->velocity.x += body->acceleration.x * dt;
            body->velocity.y += body->acceleration.y * dt;
            body->velocity.z += body->acceleration.z * dt;
        }

        if (f == gravity) {
            gravity_accelerations(bodies, count, accelerations);
        } else {
            for (size_t i = 0; i < count; ++i) {
                double ax = 0.0;
                double ay = 0.0;
                double az = 0.0;
                for (size_t j = 0; j < count; ++j) {
                    if (i == j) {
                        continue;
                    }
                    Vec3 contribution = f(&bodies[i], &bodies[j]);
                    ax += contribution.x;
                    ay += contribution.y;
                    az += contribution.z;
                }
                accelerations[i] = (Vec3){ax, ay, az};
            }
        }

        for (size_t i = 0; i < count; ++i) {
            bodies[i].acceleration = accelerations[i];
        }
    }

    free(accelerations);
    return true;
}
