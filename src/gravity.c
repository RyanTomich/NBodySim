#include "gravity.h"

#include <math.h>

Vec3 gravity(const Body *target, const Body *source)
{
    Vec3 displacement = {
        source->position.x - target->position.x,
        source->position.y - target->position.y,
        source->position.z - target->position.z
    };
    double distance = hypot(hypot(displacement.x, displacement.y), displacement.z);
    if (distance == 0.0) {
        return (Vec3){0.0, 0.0, 0.0};
    }
    double scale = source->w / (distance * distance * distance);
    return (Vec3){
        scale * displacement.x,
        scale * displacement.y,
        scale * displacement.z
    };
}
