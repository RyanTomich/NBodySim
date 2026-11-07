#ifndef CARTESIAN_INTERNAL_H
#define CARTESIAN_INTERNAL_H

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

#endif
