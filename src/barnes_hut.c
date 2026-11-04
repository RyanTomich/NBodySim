#include "parse.h"
#include "solver_internal.h"
#include "octree_internal.h"
#include "gravity_internal.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double theta; /* Zero opens every internal node. */
    size_t bucket_capacity;
    size_t max_depth;
    size_t node_capacity; /* Zero selects max(1, 4 * count). */
} BarnesHutConfig;

typedef struct {
    size_t calls, capacity_limited_leaves, bounds_fallbacks, aggregate_fallbacks;
    size_t source_node_visits, direct_evaluations, aggregate_evaluations;
} BarnesHutTotals;

typedef struct {
    BarnesHutConfig config;
    Octree tree;
    bool bounds_fallback;
    bool aggregate_fallback;
    size_t source_node_visits; /* Source nodes visited during group discovery. */
    size_t direct_evaluations; /* Directed particle contributions, excluding self. */
    size_t aggregate_evaluations; /* Aggregate contributions across all targets. */
    BarnesHutTotals totals; /* Cumulative; caller may zero between runs. */
} BarnesHutWorkspace;

static BarnesHutConfig barnes_hut_default_config(void)
{
    return (BarnesHutConfig){
        .theta = 0.5,
        .bucket_capacity = 8,
        .max_depth = 32,
        .node_capacity = 0
    };
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

static bool aggregate_tree(const SimulationView *state, BarnesHutWorkspace *workspace)
{
    /* Children follow their parents in the pool; visit them first. */
    for (size_t n = workspace->tree.node_count; n-- > 0;) {
        OctreeNode *node = &workspace->tree.nodes[n];
        if (node->is_leaf) {
            for (size_t p = node->begin; p < node->end; ++p) {
                size_t i = workspace->tree.indices[p];
                if (!add_mass(node, state->w[i],
                              (Vec3){state->x[i], state->y[i], state->z[i]})) {
                    return false;
                }
            }
        } else {
            for (size_t c = 0; c < 8; ++c) {
                if (node->children[c] != SIZE_MAX) {
                    const OctreeNode *child = &workspace->tree.nodes[node->children[c]];
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
                                          const BarnesHutWorkspace *workspace,
                                          const OctreeNode *group)
{
    const double *coordinates[] = {state->x, state->y, state->z};
    TargetBounds bounds;
    for (size_t axis = 0; axis < 3; ++axis) {
        double low = coordinates[axis][workspace->tree.indices[group->begin]];
        double high = low;
        for (size_t p = group->begin + 1; p < group->end; ++p) {
            double value = coordinates[axis][workspace->tree.indices[p]];
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

static void collect_sources(BarnesHutWorkspace *workspace, size_t index,
                            const OctreeNode *group, const TargetBounds *bounds,
                            size_t *length)
{
    const OctreeNode *node = &workspace->tree.nodes[index];
    ++workspace->source_node_visits;
    bool disjoint = node->end <= group->begin || group->end <= node->begin;
    if (node->is_leaf || (disjoint && accept_source(node, bounds, workspace->config.theta))) {
        /* Disjoint, nonempty source ranges bound the list length by count. */
        assert(*length < workspace->tree.count);
        workspace->tree.scratch[(*length)++] = index;
        return;
    }
    for (size_t c = 0; c < 8; ++c) {
        if (node->children[c] != SIZE_MAX) {
            collect_sources(workspace, node->children[c], group, bounds, length);
        }
    }
}

static void evaluate_group(const SimulationView *state, double *ax, double *ay,
                           double *az, BarnesHutWorkspace *workspace,
                           const OctreeNode *group, size_t length)
{
    const size_t *sources = workspace->tree.scratch;
    for (size_t p = group->begin; p < group->end; ++p) {
        size_t target = workspace->tree.indices[p];
        Vec3 sum = {0};
        for (size_t entry = 0; entry < length; ++entry) {
            const OctreeNode *node = &workspace->tree.nodes[sources[entry]];
            if (node->is_leaf) {
                for (size_t q = node->begin; q < node->end; ++q) {
                    size_t source = workspace->tree.indices[q];
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
                     double *az, BarnesHutWorkspace *workspace)
{
    workspace->source_node_visits = 0;
    workspace->direct_evaluations = 0;
    workspace->aggregate_evaluations = 0;

    workspace->bounds_fallback = !octree_build(state, &workspace->tree);
    workspace->aggregate_fallback = !workspace->bounds_fallback
                                 && !aggregate_tree(state, workspace);
    if (workspace->bounds_fallback || workspace->aggregate_fallback) {
        /* Count directed contributions, even though direct gravity shares pairs. */
        workspace->direct_evaluations = state->count * (state->count - 1);
        gravity_direct(state, ax, ay, az, NULL);
        return;
    }

    for (size_t n = 0; n < workspace->tree.node_count; ++n) {
        const OctreeNode *group = &workspace->tree.nodes[n];
        if (group->is_leaf) {
            TargetBounds bounds = compute_target_bounds(state, workspace, group);
            size_t length = 0;
            collect_sources(workspace, 0, group, &bounds, &length);
            evaluate_group(state, ax, ay, az, workspace, group, length);
        }
    }
}

static void gravity_barnes_hut(const SimulationView *state, double *ax, double *ay,
                         double *az, void *context)
{
    BarnesHutWorkspace *workspace = context;
    if (state->count == 0 && !workspace) {
        return;
    }
    assert(workspace && workspace->tree.count == state->count && workspace->tree.nodes);
    evaluate(state, ax, ay, az, workspace);
    ++workspace->totals.calls;
    workspace->totals.capacity_limited_leaves += workspace->tree.capacity_limited_leaves;
    workspace->totals.bounds_fallbacks += workspace->bounds_fallback;
    workspace->totals.aggregate_fallbacks += workspace->aggregate_fallback;
    workspace->totals.source_node_visits += workspace->source_node_visits;
    workspace->totals.direct_evaluations += workspace->direct_evaluations;
    workspace->totals.aggregate_evaluations += workspace->aggregate_evaluations;
}

static bool barnes_hut_create(void **context, size_t count, const SolverOption *options,
                               size_t option_count, const char **error)
{
    BarnesHutConfig config = barnes_hut_default_config();
    for (size_t i = 0; i < option_count; ++i) {
        const char *name = options[i].name, *value = options[i].value;
        bool valid = false;
        if (strcmp(name, "theta") == 0) {
            valid = parse_double(value, &config.theta) && config.theta >= 0;
        } else if (strcmp(name, "bucket_capacity") == 0) {
            valid = parse_size(value, &config.bucket_capacity) && config.bucket_capacity > 0;
        } else if (strcmp(name, "max_depth") == 0) {
            valid = parse_size(value, &config.max_depth) && config.max_depth <= OCTREE_MAX_DEPTH;
        } else if (strcmp(name, "node_capacity") == 0) {
            valid = parse_size(value, &config.node_capacity);
        }
        if (!valid) {
            *error = "Unknown solver option or invalid value.";
            return false;
        }
    }
    BarnesHutWorkspace *workspace = calloc(1, sizeof(*workspace));
    if (!workspace || !octree_create(&workspace->tree, count, config.bucket_capacity,
                                     config.max_depth, config.node_capacity)) {
        free(workspace);
        *error = "Could not allocate Barnes-Hut workspace.";
        return false;
    }
    workspace->config = config;
    *context = workspace;
    return true;
}

static void barnes_hut_destroy(void *context)
{
    octree_destroy(&((BarnesHutWorkspace *)context)->tree);
    free(context);
}

static void barnes_hut_settings(const void *context, MetadataVisitor visit, void *user)
{
    const BarnesHutWorkspace *workspace = context;
    visit("theta", (MetadataValue){METADATA_DOUBLE, {.real = workspace->config.theta}}, user);
    visit("bucket_capacity", (MetadataValue){METADATA_SIZE, {.size = workspace->config.bucket_capacity}}, user);
    visit("max_depth", (MetadataValue){METADATA_SIZE, {.size = workspace->config.max_depth}}, user);
    visit("node_capacity", (MetadataValue){METADATA_SIZE, {.size = workspace->tree.node_capacity}}, user);
}

static void barnes_hut_diagnostics(const void *context, MetadataVisitor visit, void *user)
{
    const BarnesHutTotals *totals = &((const BarnesHutWorkspace *)context)->totals;
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

const SolverDescriptor barnes_hut_solver = {
    .name = "barnes-hut",
    .acceleration = gravity_barnes_hut,
    .create = barnes_hut_create,
    .destroy = barnes_hut_destroy,
    .settings = barnes_hut_settings,
    .diagnostics = barnes_hut_diagnostics
};
