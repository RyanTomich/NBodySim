#ifndef GRAVITY_INTERNAL_H
#define GRAVITY_INTERNAL_H

#include "nbody.h"

/* Add an indexed source range to an existing sum, in range order. Self is
 * excluded by index; distinct coincident particles contribute zero. */
Vec3 gravity_accumulate_sources(const SimulationView *state, size_t target,
                                const size_t *indices, size_t begin, size_t end,
                                Vec3 sum);

/* Serial kernel: evaluate each unordered pair once in source-index order.
 * Outputs need not be initialized and must be disjoint from each other and state.
 * State is read-only; coincident positions contribute zero. No allocation.
 */
void n2_accelerations_serial(const SimulationView *state, double *restrict ax,
                             double *restrict ay, double *restrict az);

#endif
