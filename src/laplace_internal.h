#ifndef LAPLACE_INTERNAL_H
#define LAPLACE_INTERNAL_H

enum { LAPLACE_ORDER = 4, LAPLACE_TERMS = 35, LAPLACE_SPLITS = 210 };

typedef struct {
    int exponent[3];
    int lower[3]; /* Index after lowering the axis exponent, or -1 if zero. */
} CartesianTerm;

/* Multi-indices satisfy a = b + c, with factor choose(a,b) componentwise.
 * For an upward or downward shift, a is the parent term, b the child term,
 * and c the center-offset monomial. The divisor rescales b for the child.
 * In translation, a indexes the kernel, b the target local, and c the source moment. */
typedef struct {
    int a, b, c, factor, child_divisor;
} CartesianSplit;

extern const CartesianTerm laplace_terms[LAPLACE_TERMS];
extern const CartesianSplit laplace_splits[LAPLACE_SPLITS];

/* Add a positive 1/r source expansion to a target's local polynomial.
 * Moments and locals use coordinates scaled by their cell half-widths.
 * unit is the source-to-target center displacement divided by distance;
 * distance must be positive and finite, with well-separated source/target cells.
 * The translation truncates combined source/target degree at four.
 * All arrays are disjoint. The function accumulates into local without allocating.
 */
void laplace_translate(const double *moments,
                       double source_half_width, double target_half_width,
                       const double unit[3], double distance, double *local);

#endif
