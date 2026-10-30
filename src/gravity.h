#ifndef GRAVITY_H
#define GRAVITY_H

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

/* Whole-system solvers. Direct ignores context. Approximate requires a prepared
 * ApproximateWorkspace (approximate.h) matching state->count for nonempty calls.
 * It rebuilds a bounded octree and shares source lists within each target leaf.
 * Unrepresentable bounds or aggregates fall back to direct gravity.
 * No allocation occurs in either callback. Empty approximate calls allow NULL.
 */
void gravity_direct(const SimulationView *state, double *ax, double *ay,
                    double *az, void *context);
void gravity_approximate(const SimulationView *state, double *ax, double *ay,
                         double *az, void *context);
/* Resolve "direct" or "approximate"; unknown names return NULL. */
AccelerationSolver gravity_solver_by_name(const char *name);

#endif
