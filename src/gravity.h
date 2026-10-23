#ifndef GRAVITY_H
#define GRAVITY_H

#include "nbody.h"

/* Newtonian acceleration in normalized units (G = 1), with w as mass.
 * Coincident positions contribute zero.
 */
Vec3 gravity(const Body *target, const Body *source);

#endif
