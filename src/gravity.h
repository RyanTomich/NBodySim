#ifndef GRAVITY_H
#define GRAVITY_H

#include "nbody.h"

/* Newtonian acceleration in normalized units (G = 1), with w as mass.
 * Coincident positions contribute zero.
 */
Vec3 gravity(const Body *target, const Body *source);

/* Fill a separate count-element workspace, evaluating each unordered pair once.
 * Bodies are read-only; coincident positions contribute zero.
 */
void gravity_accelerations(const Body *bodies, size_t count, Vec3 *accelerations);

#endif
