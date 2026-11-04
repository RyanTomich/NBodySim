#ifndef GRAVITY_INTERNAL_H
#define GRAVITY_INTERNAL_H

#include "nbody.h"

/* Newtonian acceleration in normalized units (G = 1), with w as mass.
 * Coincident positions contribute zero.
 */
Vec3 gravity(const SimulationView *state, size_t target, size_t source);

/* Evaluate each unordered pair once, accumulating in source-index order.
 * Outputs need not be initialized and must be disjoint from each other and state.
 * State is read-only; coincident positions contribute zero. No allocation.
 */
void gravity_accelerations(const SimulationView *state, double *restrict ax,
                           double *restrict ay, double *restrict az);

/* Whole-system direct callback; context is ignored. */
void gravity_direct(const SimulationView *state, double *ax, double *ay,
                    double *az, void *context);
#endif
