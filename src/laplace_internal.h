#ifndef LAPLACE_INTERNAL_H
#define LAPLACE_INTERNAL_H

#include "cartesian_internal.h"

/* Add a positive 1/r source expansion to a target's local polynomial.
 * Moments and locals use coordinates scaled by their cell half-widths.
 * unit is the source-to-target center displacement divided by distance;
 * distance must be positive and finite, with well-separated source/target cells.
 * Prepared tables truncate combined source/target degree at order (1..6).
 * kernel is reusable scratch of terms->terms doubles. All arrays are disjoint;
 * local is accumulated, kernel is overwritten, and no allocation occurs.
 */
void laplace_translate(const CartesianTerms *terms, const double *moments,
                       double source_half_width, double target_half_width,
                       const double unit[3], double distance, double *kernel,
                       double *local);

#endif
