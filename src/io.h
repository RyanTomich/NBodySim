#ifndef NBODY_IO_H
#define NBODY_IO_H

#include "nbody.h"

#include <stdio.h>

/* Parse a decimal nonnegative integer or a finite floating-point number. */
bool parse_size(const char *text, size_t *value);
bool parse_double(const char *text, double *value);

/* Caller owns the returned array and must free it. Failure sets NULL and zero.
 * Format: count followed by count rows of w x y z vx vy vz ax ay az.
 */
bool read_bodies(FILE *input, Body **bodies, size_t *count);
bool write_bodies(FILE *output, const Body *bodies, size_t count);

#endif
