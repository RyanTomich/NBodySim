#ifndef APPROXIMATE_INTERNAL_H
#define APPROXIMATE_INTERNAL_H

#include "nbody.h"

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
    size_t begin, end; /* Half-open range in workspace->indices. */
    double mass;
    Vec3 center_of_mass;
    bool is_leaf;
} OctreeNode;

typedef struct {
    size_t calls, capacity_limited_leaves, bounds_fallbacks, aggregate_fallbacks;
    size_t source_node_visits, direct_evaluations, aggregate_evaluations;
} ApproximateTotals;

typedef struct {
    ApproximateConfig config;
    size_t count, node_capacity, node_count;
    OctreeNode *nodes;
    size_t *indices;
    /* count entries, reused for tree partitioning and group interaction lists. */
    size_t *scratch;
    /* Diagnostics describe the most recent solver call, including empty calls. */
    size_t capacity_limited_leaves;
    bool bounds_fallback;
    bool aggregate_fallback;
    size_t source_node_visits; /* Source nodes visited during group discovery. */
    size_t direct_evaluations; /* Directed particle contributions, excluding self. */
    size_t aggregate_evaluations; /* Aggregate contributions across all targets. */
    ApproximateTotals totals; /* Cumulative; caller may zero between runs. */
} ApproximateWorkspace;

ApproximateConfig approximate_default_config(void);
/* Destination must not own storage. Failure leaves it empty. No state pointers
 * are retained. count must match every subsequent solver invocation. */
bool approximate_workspace_create(ApproximateWorkspace *workspace,
                                  size_t count, ApproximateConfig config);
/* Accepts zero-initialized/partially initialized storage; clears ownership. */
void approximate_workspace_destroy(ApproximateWorkspace *workspace);

/* Rebuild the bounded tree and share source lists within each target leaf.
 * Context is a prepared workspace matching state->count. No allocation.
 * Unrepresentable bounds/aggregates use direct fallback; empty calls allow NULL.
 */
void gravity_approximate(const SimulationView *state, double *ax, double *ay,
                         double *az, void *context);

#endif
