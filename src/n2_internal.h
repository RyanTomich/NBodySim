#ifndef N2_INTERNAL_H
#define N2_INTERNAL_H

#include "nbody.h"

/* Serial kernel: evaluate each unordered pair once in source-index order.
 * Outputs need not be initialized and must be disjoint from each other and state.
 * State is read-only; coincident positions contribute zero. No allocation.
 */
void n2_accelerations(const SimulationView *state, double *restrict ax,
                           double *restrict ay, double *restrict az);

#endif
