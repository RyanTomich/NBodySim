#ifndef NBODY_IO_H
#define NBODY_IO_H

#include "nbody.h"

#include <stdio.h>

/* Read directly into array storage; destination must not own an allocation.
 * Caller owns the result and must call simulation_destroy(*simulation). Failure leaves an
 * NULL handle, including when input is malformed after allocation.
 * Format: count followed by count rows of w x y z vx vy vz ax ay az.
 */
bool read_bodies(FILE *input, Simulation **simulation);
bool write_bodies(FILE *output, const Simulation *simulation);

#endif
