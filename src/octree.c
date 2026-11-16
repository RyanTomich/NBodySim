#include "octree_internal.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void octree_destroy(Octree *workspace)
{
    free(workspace->nodes);
    free(workspace->indices);
    free(workspace->scratch);
    *workspace = (Octree){0};
}

bool octree_create(Octree *workspace, size_t count, size_t bucket_capacity,
                   size_t max_depth, size_t node_capacity)
{
    *workspace = (Octree){0};
    if (bucket_capacity == 0 || max_depth > OCTREE_MAX_DEPTH
        || count > SIZE_MAX / sizeof(size_t)) {
        return false;
    }
    size_t capacity = node_capacity;
    if (capacity == 0) {
        if (count > SIZE_MAX / 4) {
            return false;
        }
        capacity = count ? 4 * count : 1;
    }
    if (capacity > SIZE_MAX / sizeof(OctreeNode)) {
        return false;
    }
    workspace->bucket_capacity = bucket_capacity;
    workspace->max_depth = max_depth;
    workspace->count = count;
    workspace->node_capacity = capacity;
    workspace->nodes = malloc(capacity * sizeof(*workspace->nodes));
    if (count != 0) {
        workspace->indices = malloc(count * sizeof(*workspace->indices));
        workspace->scratch = malloc(count * sizeof(*workspace->scratch));
    }
    if (!workspace->nodes || (count && (!workspace->indices || !workspace->scratch))) {
        octree_destroy(workspace);
        return false;
    }
    return true;
}

static OctreeNode make_node(Vec3 center, double half_width, size_t begin, size_t end)
{
    OctreeNode node = {
        .center = center,
        .half_width = half_width,
        .begin = begin,
        .end = end,
        .is_leaf = true
    };
    for (size_t i = 0; i < 8; ++i) {
        node.children[i] = SIZE_MAX;
    }
    return node;
}

static unsigned octant(const SimulationView *state, size_t i, Vec3 center)
{
    return (state->x[i] >= center.x)
         | ((unsigned)(state->y[i] >= center.y) << 1)
         | ((unsigned)(state->z[i] >= center.z) << 2);
}

static void subdivide(const SimulationView *state, Octree *workspace,
                      size_t index, size_t depth)
{
    OctreeNode *node = &workspace->nodes[index];
    if (node->end - node->begin <= workspace->bucket_capacity
        || depth >= workspace->max_depth) {
        return;
    }
    double half = node->half_width * 0.5;
    Vec3 centers[8];
    for (unsigned c = 0; c < 8; ++c) {
        centers[c] = (Vec3){node->center.x + (c & 1 ? half : -half),
                            node->center.y + (c & 2 ? half : -half),
                            node->center.z + (c & 4 ? half : -half)};
        if (!isfinite(centers[c].x) || !isfinite(centers[c].y)
            || !isfinite(centers[c].z) || centers[c].x == node->center.x
            || centers[c].y == node->center.y || centers[c].z == node->center.z) {
            return;
        }
    }
    size_t counts[8] = {0};
    for (size_t p = node->begin; p < node->end; ++p) {
        size_t i = workspace->indices[p];
        unsigned c = octant(state, i, node->center);
        /* Keep a leaf if rounding would place a particle outside its child. */
        if (fabs(state->x[i] - centers[c].x) > half
            || fabs(state->y[i] - centers[c].y) > half
            || fabs(state->z[i] - centers[c].z) > half) {
            return;
        }
        ++counts[c];
    }
    size_t occupied = 0;
    for (size_t c = 0; c < 8; ++c) {
        occupied += counts[c] != 0;
    }
    if (occupied > workspace->node_capacity - workspace->node_count) {
        ++workspace->capacity_limited_leaves;
        return;
    }
    size_t starts[8], cursors[8];
    size_t end = node->begin;
    for (size_t c = 0; c < 8; ++c) {
        starts[c] = cursors[c] = end;
        end += counts[c];
    }
    for (size_t p = node->begin; p < node->end; ++p) {
        size_t i = workspace->indices[p];
        workspace->scratch[cursors[octant(state, i, node->center)]++] = i;
    }
    memcpy(workspace->indices + node->begin, workspace->scratch + node->begin,
           (node->end - node->begin) * sizeof(*workspace->indices));
    /* Reserve siblings before a deeper subtree can exhaust the node pool. */
    for (size_t c = 0; c < 8; ++c) {
        if (counts[c]) {
            size_t child = workspace->node_count++;
            node->children[c] = child;
            node->is_leaf = false;
            workspace->nodes[child] = make_node(centers[c], half, starts[c], cursors[c]);
        }
    }
    for (size_t c = 0; c < 8; ++c) {
        if (node->children[c] != SIZE_MAX) {
            subdivide(state, workspace, node->children[c], depth + 1);
        }
    }
}

static bool compute_root_bounds(const SimulationView *state, Vec3 *center, double *half)
{
    const double *coordinates[] = {state->x, state->y, state->z};
    double mid[3];
    *half = 0.0;
    double scale = 1.0;
    for (size_t axis = 0; axis < 3; ++axis) {
        double low = coordinates[axis][0], high = low;
        for (size_t i = 0; i < state->count; ++i) {
            double value = coordinates[axis][i];
            if (!isfinite(value)) {
                return false;
            }
            low = fmin(low, value);
            high = fmax(high, value);
        }
        mid[axis] = low * 0.5 + high * 0.5;
        *half = fmax(*half, fmax(high - mid[axis], mid[axis] - low));
        scale = fmax(scale, fmax(fabs(low), fabs(high)));
    }
    if (*half == 0.0) {
        *half = scale * DBL_EPSILON;
    }
    *half = nextafter(*half, INFINITY);
    if (!isfinite(*half) || *half > DBL_MAX * 0.5) {
        return false; /* Full cube width must also be representable. */
    }
    for (size_t axis = 0; axis < 3; ++axis) {
        if (!isfinite(mid[axis] - *half) || !isfinite(mid[axis] + *half)) {
            return false;
        }
    }
    *center = (Vec3){mid[0], mid[1], mid[2]};
    return true;
}

bool octree_build(const SimulationView *state, Octree *workspace)
{
    workspace->capacity_limited_leaves = 0;
    workspace->node_count = 0;
    for (size_t i = 0; i < state->count; ++i) {
        workspace->indices[i] = i;
    }
    if (state->count == 0) {
        return true;
    }

    Vec3 center;
    double half_width;
    if (!compute_root_bounds(state, &center, &half_width)) {
        return false;
    }
    workspace->nodes[0] = make_node(center, half_width, 0, state->count);
    workspace->node_count = 1;
    subdivide(state, workspace, 0, 0);
    return true;
}

bool octree_membership_is_valid(const SimulationView *state, const Octree *tree)
{
    if (state->count != tree->count) {
        return false;
    }
    if (state->count == 0) {
        return tree->node_count == 0;
    }
    /* Build failures leave node_count zero, even after an earlier valid tree. */
    if (!tree->nodes || !tree->indices || tree->node_count == 0
        || tree->node_count > tree->node_capacity) {
        return false;
    }
    const double *coordinates[] = {state->x, state->y, state->z};
    for (size_t n = 0; n < tree->node_count; ++n) {
        const OctreeNode *node = &tree->nodes[n];
        if (!node->is_leaf) {
            continue;
        }
        const double center[] = {node->center.x, node->center.y, node->center.z};
        if (!isfinite(node->half_width) || node->half_width <= 0.0) {
            return false;
        }
        for (size_t p = node->begin; p < node->end; ++p) {
            size_t i = tree->indices[p];
            for (size_t axis = 0; axis < 3; ++axis) {
                double position = coordinates[axis][i];
                double displacement = position - center[axis];
                /* A rounded difference equal to half_width could hide an
                 * outside point. Rebuild on equality rather than widen a cell.
                 * Subtraction also avoids rounding expanded cube endpoints. */
                if (!isfinite(position) || !isfinite(center[axis])
                    || !isfinite(displacement) || fabs(displacement) >= node->half_width) {
                    return false;
                }
            }
        }
    }
    return true;
}
