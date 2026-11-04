#include "parse.h"
#include <string.h>
#include "solver_internal.h"
#include "approximate_internal.h"
#include "gravity_internal.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

ApproximateConfig approximate_default_config(void)
{
    return (ApproximateConfig){
        .theta = 0.5,
        .bucket_capacity = 8,
        .max_depth = 32,
        .node_capacity = 0
    };
}

void approximate_workspace_destroy(ApproximateWorkspace *workspace)
{
    free(workspace->nodes);
    free(workspace->indices);
    free(workspace->scratch);
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
    }
    if (!workspace->nodes || (count && (!workspace->indices || !workspace->scratch))) {
        approximate_workspace_destroy(workspace);
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

static bool build_tree(const SimulationView *state, ApproximateWorkspace *workspace)
{
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

static bool aggregate_tree(const SimulationView *state, ApproximateWorkspace *workspace)
{
    /* Children follow their parents in the pool; visit them first. */
    for (size_t n = workspace->node_count; n-- > 0;) {
        OctreeNode *node = &workspace->nodes[n];
        if (node->is_leaf) {
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

static Vec3 point_mass_acceleration(Vec3 displacement, double mass, double distance)
{
    /* Equivalent to M * displacement / r^3 without explicitly cubing r. */
    double scale = (mass / distance) / distance;
    return (Vec3){scale * (displacement.x / distance),
                  scale * (displacement.y / distance),
                  scale * (displacement.z / distance)};
}

typedef struct {
    double low[3], high[3];
} TargetBounds;

static TargetBounds compute_target_bounds(const SimulationView *state,
                                          const ApproximateWorkspace *workspace,
                                          const OctreeNode *group)
{
    const double *coordinates[] = {state->x, state->y, state->z};
    TargetBounds bounds;
    for (size_t axis = 0; axis < 3; ++axis) {
        double low = coordinates[axis][workspace->indices[group->begin]];
        double high = low;
        for (size_t p = group->begin + 1; p < group->end; ++p) {
            double value = coordinates[axis][workspace->indices[p]];
            low = fmin(low, value);
            high = fmax(high, value);
        }
        bounds.low[axis] = low;
        bounds.high[axis] = high;
    }
    return bounds;
}

static bool accept_source(const OctreeNode *source, const TargetBounds *bounds,
                          double theta)
{
    if (theta == 0.0) {
        return false;
    }
    const double center[] = {source->center_of_mass.x, source->center_of_mass.y,
                             source->center_of_mass.z};
    double squared_distance = 0.0;
    double far_distance = 0.0;
    for (size_t axis = 0; axis < 3; ++axis) {
        /* An upper L1 bound ensures every target displacement and norm is
         * representable, not just the nearest target's distance. */
        double far = fmax(fabs(center[axis] - bounds->low[axis]),
                          fabs(center[axis] - bounds->high[axis]));
        far_distance += far;
        double gap = fmax(0.0, fmax(bounds->low[axis] - center[axis],
                                   center[axis] - bounds->high[axis]));
        double square = gap * gap;
        squared_distance += square;
    }
    /* Nonnegative sums retain any overflow. Zero distance also requires opening. */
    if (!isfinite(far_distance) || !isfinite(squared_distance) || squared_distance == 0.0) {
        return false;
    }
    double distance = sqrt(squared_distance);
    return (2.0 * source->half_width) / distance < theta;
}

static void collect_sources(ApproximateWorkspace *workspace, size_t index,
                            const OctreeNode *group, const TargetBounds *bounds,
                            size_t *length)
{
    const OctreeNode *node = &workspace->nodes[index];
    ++workspace->source_node_visits;
    bool disjoint = node->end <= group->begin || group->end <= node->begin;
    if (node->is_leaf || (disjoint && accept_source(node, bounds, workspace->config.theta))) {
        /* Disjoint, nonempty source ranges bound the list length by count. */
        assert(*length < workspace->count);
        workspace->scratch[(*length)++] = index;
        return;
    }
    for (size_t c = 0; c < 8; ++c) {
        if (node->children[c] != SIZE_MAX) {
            collect_sources(workspace, node->children[c], group, bounds, length);
        }
    }
}

static void evaluate_group(const SimulationView *state, double *ax, double *ay,
                           double *az, ApproximateWorkspace *workspace,
                           const OctreeNode *group, size_t length)
{
    const size_t *sources = workspace->scratch;
    for (size_t p = group->begin; p < group->end; ++p) {
        size_t target = workspace->indices[p];
        Vec3 sum = {0};
        for (size_t entry = 0; entry < length; ++entry) {
            const OctreeNode *node = &workspace->nodes[sources[entry]];
            if (node->is_leaf) {
                for (size_t q = node->begin; q < node->end; ++q) {
                    size_t source = workspace->indices[q];
                    if (source == target) {
                        continue;
                    }
                    Vec3 contribution = gravity(state, target, source);
                    sum.x += contribution.x;
                    sum.y += contribution.y;
                    sum.z += contribution.z;
                    ++workspace->direct_evaluations;
                }
            } else {
                Vec3 displacement = {node->center_of_mass.x - state->x[target],
                                     node->center_of_mass.y - state->y[target],
                                     node->center_of_mass.z - state->z[target]};
                double distance = hypot(hypot(displacement.x, displacement.y), displacement.z);
                Vec3 contribution = point_mass_acceleration(displacement, node->mass, distance);
                sum.x += contribution.x;
                sum.y += contribution.y;
                sum.z += contribution.z;
                ++workspace->aggregate_evaluations;
            }
        }
        ax[target] = sum.x;
        ay[target] = sum.y;
        az[target] = sum.z;
    }
}

static void evaluate(const SimulationView *state, double *ax, double *ay,
                     double *az, ApproximateWorkspace *workspace)
{
    workspace->capacity_limited_leaves = 0;
    workspace->source_node_visits = 0;
    workspace->direct_evaluations = 0;
    workspace->aggregate_evaluations = 0;

    workspace->bounds_fallback = !build_tree(state, workspace);
    workspace->aggregate_fallback = !workspace->bounds_fallback
                                 && !aggregate_tree(state, workspace);
    if (workspace->bounds_fallback || workspace->aggregate_fallback) {
        /* Count directed contributions, even though direct gravity shares pairs. */
        workspace->direct_evaluations = state->count * (state->count - 1);
        gravity_direct(state, ax, ay, az, NULL);
        return;
    }

    for (size_t n = 0; n < workspace->node_count; ++n) {
        const OctreeNode *group = &workspace->nodes[n];
        if (group->is_leaf) {
            TargetBounds bounds = compute_target_bounds(state, workspace, group);
            size_t length = 0;
            collect_sources(workspace, 0, group, &bounds, &length);
            evaluate_group(state, ax, ay, az, workspace, group, length);
        }
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
    workspace->totals.source_node_visits += workspace->source_node_visits;
    workspace->totals.direct_evaluations += workspace->direct_evaluations;
    workspace->totals.aggregate_evaluations += workspace->aggregate_evaluations;
}

static bool approximate_create(void **context, size_t count, const SolverOption *options,
                               size_t option_count, const char **error)
{
    ApproximateConfig config = approximate_default_config();
    for (size_t i = 0; i < option_count; ++i) {
        const char *name = options[i].name, *value = options[i].value;
        bool valid = false;
        if (strcmp(name, "theta") == 0) {
            valid = parse_double(value, &config.theta) && config.theta >= 0;
        } else if (strcmp(name, "bucket_capacity") == 0) {
            valid = parse_size(value, &config.bucket_capacity) && config.bucket_capacity > 0;
        } else if (strcmp(name, "max_depth") == 0) {
            valid = parse_size(value, &config.max_depth) && config.max_depth <= APPROXIMATE_MAX_DEPTH;
        } else if (strcmp(name, "node_capacity") == 0) {
            valid = parse_size(value, &config.node_capacity);
        }
        if (!valid) {
            *error = "Unknown solver option or invalid value.";
            return false;
        }
    }
    ApproximateWorkspace *workspace = malloc(sizeof(*workspace));
    if (!workspace || !approximate_workspace_create(workspace, count, config)) {
        free(workspace);
        *error = "Could not allocate approximate workspace.";
        return false;
    }
    *context = workspace;
    return true;
}

static void approximate_destroy(void *context)
{
    approximate_workspace_destroy(context);
    free(context);
}

static void approximate_settings(const void *context, MetadataVisitor visit, void *user)
{
    const ApproximateWorkspace *workspace = context;
    visit("theta", (MetadataValue){METADATA_DOUBLE, {.real = workspace->config.theta}}, user);
    visit("bucket_capacity", (MetadataValue){METADATA_SIZE, {.size = workspace->config.bucket_capacity}}, user);
    visit("max_depth", (MetadataValue){METADATA_SIZE, {.size = workspace->config.max_depth}}, user);
    visit("node_capacity", (MetadataValue){METADATA_SIZE, {.size = workspace->node_capacity}}, user);
}

static void approximate_diagnostics(const void *context, MetadataVisitor visit, void *user)
{
    const ApproximateTotals *totals = &((const ApproximateWorkspace *)context)->totals;
#define REPORT(field) visit(#field, (MetadataValue){METADATA_SIZE, {.size = totals->field}}, user)
    REPORT(calls);
    REPORT(capacity_limited_leaves);
    REPORT(bounds_fallbacks);
    REPORT(aggregate_fallbacks);
    /* Historical result field, derived rather than counted twice. */
    visit("accepted_nodes", (MetadataValue){METADATA_SIZE, {.size = totals->aggregate_evaluations}}, user);
    REPORT(source_node_visits);
    REPORT(direct_evaluations);
    REPORT(aggregate_evaluations);
#undef REPORT
}

const SolverDescriptor approximate_solver = {
    .name = "approximate",
    .acceleration = gravity_approximate,
    .create = approximate_create,
    .destroy = approximate_destroy,
    .settings = approximate_settings,
    .diagnostics = approximate_diagnostics
};
