#ifndef OCTREE_INTERNAL_H
#define OCTREE_INTERNAL_H
#include "nbody.h"

#define OCTREE_MAX_DEPTH 64

/* Successful nonempty builds store the root at nodes[0] and give every node
 * a nonempty, contiguous range in indices, which is a permutation of body
 * indices. A node's range covers its whole subtree. An internal node's child
 * ranges partition its range in octant order; leaf ranges partition [0, count).
 * Thus unrelated subtrees have disjoint ranges, and a traversal that selects
 * a node or descends into its children selects at most count source nodes.
 * Children always have larger node indices than their parents, so reverse
 * node order supports bottom-up summaries and forward order top-down passes.
 * Subtree nodes need not be contiguous in the node pool.
 * Solver-specific numerical summaries are stored separately from geometry. */
typedef struct {
    Vec3 center;
    double half_width;
    size_t children[8]; /* SIZE_MAX means absent. */
    size_t begin, end; /* Half-open subtree range in tree->indices. */
    bool is_leaf;
} OctreeNode;

typedef struct {
    size_t bucket_capacity, max_depth, count, node_capacity, node_count;
    size_t capacity_limited_leaves;
    OctreeNode *nodes;
    size_t *indices, *scratch;
} Octree;

bool octree_create(Octree *tree, size_t count, size_t bucket_capacity,
                   size_t max_depth, size_t node_capacity);
void octree_destroy(Octree *tree);
/* Requires the prepared body count. Empty input succeeds with no nodes;
 * unrepresentable bounds fail with node_count zero. Capacity/depth limits
 * retain unsplit leaves and do not fail the build. */
bool octree_build(const SimulationView *state, Octree *tree);
/* Check every body's current membership in its previously assigned leaf.
 * Assumes an unmodified topology produced by octree_build; this is not a
 * structural validator. Conservative at rounded cell boundaries: false
 * requires rebuilding. True permits topology reuse, but numerical summaries
 * must still be refreshed from current masses and positions. */
bool octree_membership_is_valid(const SimulationView *state, const Octree *tree);
#endif
