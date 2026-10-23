#ifndef NBODY_IO_H
#define NBODY_IO_H

#include "nbody.h"

#include <stdio.h>

/* Parse a decimal nonnegative integer or a finite floating-point number. */
bool parse_size(const char *text, size_t *value);
bool parse_double(const char *text, double *value);

/* Read directly into array storage; destination must not own an allocation.
 * Caller owns the result and must free(simulation->storage). Failure leaves an
 * empty simulation, including when input is malformed after allocation.
 * Format: count followed by count rows of w x y z vx vy vz ax ay az.
 */
bool read_bodies(FILE *input, Simulation *simulation);
bool write_bodies(FILE *output, const Simulation *simulation);

#endif
