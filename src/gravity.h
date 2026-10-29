#ifndef GRAVITY_H
#define GRAVITY_H

#include "nbody.h"

/* Newtonian acceleration in normalized units (G = 1), with w as mass.
 * Coincident positions contribute zero.
 */
Vec3 gravity(const SimulationView *state, size_t target, size_t source);

/* Fill three separate count-element workspace arrays using scalar i < j loops.
 * Each unordered pair is evaluated once, updating both bodies. Each body's
 * contributions are accumulated in increasing source-index order, with the
 * current target's running sums retained locally throughout its source loop.
 * The first target initializes the outputs; prior workspace values are ignored.
 * State is read-only; coincident positions contribute zero. No allocation.
 * Workspace arrays must not overlap each other or state.
 */
void gravity_accelerations(const SimulationView *state, double *restrict ax,
                           double *restrict ay, double *restrict az);

/* Whole-system solvers. Direct ignores context. Approximate requires a prepared
 * ApproximateWorkspace (approximate.h) matching state->count for nonempty calls.
 * It rebuilds a bounded octree, aggregates mass, and traverses for each target.
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
