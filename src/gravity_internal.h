#ifndef GRAVITY_INTERNAL_H
#define GRAVITY_INTERNAL_H

#include "nbody.h"

/* Add an indexed source range to an existing sum, in range order. Self is
 * excluded by index; distinct coincident particles contribute zero. */
Vec3 gravity_accumulate_sources(const SimulationView *state, size_t target,
                                const size_t *indices, size_t begin, size_t end,
                                Vec3 sum);

#endif
