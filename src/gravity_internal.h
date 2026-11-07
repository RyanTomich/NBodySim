#ifndef GRAVITY_INTERNAL_H
#define GRAVITY_INTERNAL_H

#include "nbody.h"

/* Newtonian acceleration in normalized units (G = 1), with w as mass.
 * Coincident positions contribute zero.
 */
Vec3 gravity(const SimulationView *state, size_t target, size_t source);

/* Add an indexed source range to an existing sum, in range order. Self is
 * excluded by index; distinct coincident particles contribute zero. */
Vec3 gravity_accumulate_sources(const SimulationView *state, size_t target,
                                const size_t *indices, size_t begin, size_t end,
                                Vec3 sum);

/* Positive-mass source summary, separate from spatial tree geometry.
 * Zero-initialize before adding particles or child monopoles. */
typedef struct {
    double mass;
    Vec3 center_of_mass;
} GravityMonopole;

/* Add a particle (mass, position) or merge a child (mass, center_of_mass).
 * False means the aggregate is unrepresentable; discard it and fall back. */
bool gravity_monopole_add(GravityMonopole *monopole, double mass, Vec3 position);
/* The caller must establish a nonzero, finite distance to this monopole. */
Vec3 gravity_monopole_acceleration(const GravityMonopole *monopole, Vec3 target);

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
