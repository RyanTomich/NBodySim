#ifndef APPROXIMATE_H
#define APPROXIMATE_H

#include "nbody.h"

/* Recursion is bounded even for caller-supplied configurations. */
#define APPROXIMATE_MAX_DEPTH 64

typedef struct {
    double theta; /* Zero opens every internal node. */
    size_t bucket_capacity;
    size_t max_depth;
    size_t node_capacity; /* Zero selects max(1, 4 * count). */
} ApproximateConfig;

typedef struct {
    Vec3 center;
    double half_width;
    size_t children[8]; /* SIZE_MAX means absent. */
    size_t begin, end;
    double mass;
    Vec3 center_of_mass;
} OctreeNode;

typedef struct {
    size_t calls, capacity_limited_leaves, bounds_fallbacks, aggregate_fallbacks;
    size_t accepted_nodes;
} ApproximateTotals;

typedef struct {
    ApproximateConfig config;
    size_t count, node_capacity, node_count;
    OctreeNode *nodes;
    size_t *indices, *scratch, *rank;
    /* Diagnostics describe the most recent solver call, including empty calls. */
    size_t capacity_limited_leaves;
    bool bounds_fallback;
    bool aggregate_fallback;
    size_t accepted_nodes; /* Accepted aggregates across all targets in last call. */
    ApproximateTotals totals; /* Cumulative since creation; caller may zero between runs. */
} ApproximateWorkspace;

ApproximateConfig approximate_default_config(void);
/* Destination must not own storage. Failure leaves it empty. No state pointers
 * are retained. count must match every subsequent solver invocation. */
bool approximate_workspace_create(ApproximateWorkspace *workspace,
                                  size_t count, ApproximateConfig config);
/* Accepts zero-initialized/partially initialized storage; clears ownership. */
void approximate_workspace_destroy(ApproximateWorkspace *workspace);

#endif
