#ifndef NBODY_H
#define NBODY_H

#include <stdbool.h>
#include <stddef.h>

typedef struct {
    double x, y, z;
} Vec3;

typedef struct {
    double w;
    Vec3 position;
    Vec3 velocity;
    Vec3 acceleration;
} Body;

/* Return the source's acceleration contribution to the target, without mutation. */
typedef Vec3 (*Interaction)(const Body *target, const Body *source);

/* Update in place using forward Euler, then recompute acceleration with f.
 * Input accelerations drive the first step. False means invalid arguments or
 * allocation failure; in either case the bodies are unchanged.
 */
bool simulate(Body *bodies, size_t count, double dt, size_t steps, Interaction f);

#endif
