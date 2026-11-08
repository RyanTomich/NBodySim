#ifndef LAPLACE_INTERNAL_H
#define LAPLACE_INTERNAL_H

#include <stddef.h>

enum { CARTESIAN_MAX_ORDER = 6, CARTESIAN_MAX_TERMS = 84, CARTESIAN_MAX_SPLITS = 924 };

/* Each split represents a=b+c with factor choose(a,b), componentwise.
 * 84 = choose(6+3,3) terms and 924 = choose(6+6,6) splits. */
typedef struct {
    int a, b, c, factor;
} CartesianSplit;

/* Polynomial bookkeeping only; no particle properties or spatial kernel. */
typedef struct {
    size_t order;
    int terms, splits;
    int degree[CARTESIAN_MAX_TERMS];
    int exponent[CARTESIAN_MAX_TERMS][3], lower[CARTESIAN_MAX_TERMS][3];
    CartesianSplit split[CARTESIAN_MAX_SPLITS];
} CartesianTerms;

/* Prepare a zero-initialized table whose order is already set to 1..6.
 * The caller validates the order. No allocation occurs. */
void laplace_prepare_terms(CartesianTerms *terms);

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
