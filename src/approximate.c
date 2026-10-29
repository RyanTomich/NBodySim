#include "approximate.h"
#include "gravity.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

ApproximateConfig approximate_default_config(void)
{
    return (ApproximateConfig){0.5, 8, 32, 0};
}

void approximate_workspace_destroy(ApproximateWorkspace *workspace)
{
    free(workspace->nodes);
    free(workspace->indices);
    free(workspace->scratch);
    free(workspace->rank);
    *workspace = (ApproximateWorkspace){0};
}

bool approximate_workspace_create(ApproximateWorkspace *workspace,
                                  size_t count, ApproximateConfig config)
{
    *workspace = (ApproximateWorkspace){0};
    if (!isfinite(config.theta) || config.theta < 0.0
        || config.bucket_capacity == 0 || config.max_depth > APPROXIMATE_MAX_DEPTH
        || count > SIZE_MAX / sizeof(size_t)) {
        return false;
    }
    size_t capacity = config.node_capacity;
    if (capacity == 0) {
        if (count > SIZE_MAX / 4) {
            return false;
        }
        capacity = count ? 4 * count : 1;
    }
    if (capacity > SIZE_MAX / sizeof(OctreeNode)) {
        return false;
    }
    workspace->config = config;
    workspace->count = count;
    workspace->node_capacity = capacity;
    workspace->nodes = malloc(capacity * sizeof(*workspace->nodes));
    if (count != 0) {
        workspace->indices = malloc(count * sizeof(*workspace->indices));
        workspace->scratch = malloc(count * sizeof(*workspace->scratch));
        workspace->rank = malloc(count * sizeof(*workspace->rank));
    }
    if (!workspace->nodes || (count && (!workspace->indices
        || !workspace->scratch || !workspace->rank))) {
        approximate_workspace_destroy(workspace);
        return false;
    }
    return true;
}

static OctreeNode make_node(Vec3 center, double half_width, size_t begin, size_t end)
{
    OctreeNode node = {.center = center, .half_width = half_width,
                       .begin = begin, .end = end};
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

static void subdivide(const SimulationView *state, ApproximateWorkspace *workspace,
                      size_t index, size_t depth)
{
    OctreeNode *node = &workspace->nodes[index];
    if (node->end - node->begin <= workspace->config.bucket_capacity
        || depth >= workspace->config.max_depth) {
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
        /* Rounded child centers must still enclose every assigned particle.
         * Otherwise retain the parent as a terminal bucket. */
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
    /* Reserve all siblings before descending so a deep child cannot consume
     * the slots required by its siblings. The pool never moves. */
    for (size_t c = 0; c < 8; ++c) {
        if (counts[c]) {
            size_t child = workspace->node_count++;
            node->children[c] = child;
            workspace->nodes[child] = make_node(centers[c], half, starts[c], cursors[c]);
        }
    }
    for (size_t c = 0; c < 8; ++c) {
        if (node->children[c] != SIZE_MAX) {
            subdivide(state, workspace, node->children[c], depth + 1);
        }
    }
}

static bool root_bounds(const SimulationView *state, Vec3 *center, double *half)
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

static bool is_leaf(const OctreeNode *node)
{
    for (size_t c = 0; c < 8; ++c) {
        if (node->children[c] != SIZE_MAX) {
            return false;
        }
    }
    return true;
}

/* Incremental weighted means avoid forming mass * position. */
static bool add_mass(OctreeNode *node, double mass, Vec3 position)
{
    double total = node->mass + mass;
    if (!isfinite(total) || mass <= 0.0) {
        return false;
    }
    if (node->mass == 0.0) {
        node->center_of_mass = position;
    } else {
        double fraction = mass / total;
        node->center_of_mass.x += fraction * (position.x - node->center_of_mass.x);
        node->center_of_mass.y += fraction * (position.y - node->center_of_mass.y);
        node->center_of_mass.z += fraction * (position.z - node->center_of_mass.z);
    }
    node->mass = total;
    return isfinite(node->center_of_mass.x) && isfinite(node->center_of_mass.y)
        && isfinite(node->center_of_mass.z);
}

static bool aggregate(const SimulationView *state, ApproximateWorkspace *workspace)
{
    /* Every child has a larger pool index than its parent, so reverse pool
     * order is a postorder even though siblings are reserved together. */
    for (size_t n = workspace->node_count; n-- > 0;) {
        OctreeNode *node = &workspace->nodes[n];
        if (is_leaf(node)) {
            for (size_t p = node->begin; p < node->end; ++p) {
                size_t i = workspace->indices[p];
                if (!add_mass(node, state->w[i],
                              (Vec3){state->x[i], state->y[i], state->z[i]})) {
                    return false;
                }
            }
        } else {
            for (size_t c = 0; c < 8; ++c) {
                if (node->children[c] != SIZE_MAX) {
                    const OctreeNode *child = &workspace->nodes[node->children[c]];
                    if (!add_mass(node, child->mass, child->center_of_mass)) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

static Vec3 point_mass(Vec3 displacement, double mass, double distance)
{
    /* Equivalent to M * displacement / r^3 without explicitly cubing r. */
    double scale = (mass / distance) / distance;
    return (Vec3){scale * (displacement.x / distance),
                  scale * (displacement.y / distance),
                  scale * (displacement.z / distance)};
}

static void traverse(const SimulationView *state, ApproximateWorkspace *workspace,
                     size_t index, size_t target, Vec3 *sum)
{
    const OctreeNode *node = &workspace->nodes[index];
    if (is_leaf(node)) {
        for (size_t p = node->begin; p < node->end; ++p) {
            size_t source = workspace->indices[p];
            if (source != target) {
                Vec3 contribution = gravity(state, target, source);
                sum->x += contribution.x;
                sum->y += contribution.y;
                sum->z += contribution.z;
            }
        }
        return;
    }
    size_t rank = workspace->rank[target];
    if (rank < node->begin || rank >= node->end) {
        Vec3 displacement = {node->center_of_mass.x - state->x[target],
                             node->center_of_mass.y - state->y[target],
                             node->center_of_mass.z - state->z[target]};
        double distance = hypot(hypot(displacement.x, displacement.y), displacement.z);
        if (distance > 0.0 && isfinite(distance)
            && (2.0 * node->half_width) / distance < workspace->config.theta) {
            Vec3 contribution = point_mass(displacement, node->mass, distance);
            sum->x += contribution.x;
            sum->y += contribution.y;
            sum->z += contribution.z;
            ++workspace->accepted_nodes;
            return;
        }
    }
    for (size_t c = 0; c < 8; ++c) {
        if (node->children[c] != SIZE_MAX) {
            traverse(state, workspace, node->children[c], target, sum);
        }
    }
}

static void evaluate(const SimulationView *state, double *ax, double *ay,
                     double *az, ApproximateWorkspace *workspace)
{
    workspace->node_count = 0;
    workspace->capacity_limited_leaves = 0;
    workspace->bounds_fallback = false;
    workspace->aggregate_fallback = false;
    workspace->accepted_nodes = 0;
    for (size_t i = 0; i < state->count; ++i) {
        workspace->indices[i] = workspace->rank[i] = i;
    }
    if (state->count) {
        Vec3 center;
        double half;
        if (root_bounds(state, &center, &half)) {
            workspace->nodes[0] = make_node(center, half, 0, state->count);
            workspace->node_count = 1;
            subdivide(state, workspace, 0, 0);
            for (size_t p = 0; p < state->count; ++p) {
                workspace->rank[workspace->indices[p]] = p;
            }
        } else {
            workspace->bounds_fallback = true;
        }
    }
    if (workspace->bounds_fallback) {
        gravity_direct(state, ax, ay, az, NULL);
        return;
    }
    if (!aggregate(state, workspace)) {
        workspace->aggregate_fallback = true;
        gravity_direct(state, ax, ay, az, NULL);
        return;
    }
    for (size_t i = 0; i < state->count; ++i) {
        Vec3 sum = {0};
        traverse(state, workspace, 0, i, &sum);
        ax[i] = sum.x;
        ay[i] = sum.y;
        az[i] = sum.z;
    }
}

void gravity_approximate(const SimulationView *state, double *ax, double *ay,
                         double *az, void *context)
{
    ApproximateWorkspace *workspace = context;
    if (state->count == 0 && !workspace) {
        return;
    }
    assert(workspace && workspace->count == state->count && workspace->nodes);
    evaluate(state, ax, ay, az, workspace);
    ++workspace->totals.calls;
    workspace->totals.capacity_limited_leaves += workspace->capacity_limited_leaves;
    workspace->totals.bounds_fallbacks += workspace->bounds_fallback;
    workspace->totals.aggregate_fallbacks += workspace->aggregate_fallback;
    workspace->totals.accepted_nodes += workspace->accepted_nodes;
}
