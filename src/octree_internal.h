#ifndef OCTREE_INTERNAL_H
#define OCTREE_INTERNAL_H
#include "nbody.h"

#define OCTREE_MAX_DEPTH 64

typedef struct {
    Vec3 center;
    double half_width;
    size_t children[8]; /* SIZE_MAX means absent. */
    size_t begin, end; /* Half-open range in workspace->indices. */
    union {
        double mass; /* Barnes–Hut */
        double occupied_half_width; /* FMM, about the expansion center */
    };
    Vec3 center_of_mass;
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
bool octree_build(const SimulationView *state, Octree *tree);
#endif
