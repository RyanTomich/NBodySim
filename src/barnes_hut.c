#include "parse.h"
#include "solver_internal.h"
#include "octree_internal.h"
#include "gravity_internal.h"
#include "n2_internal.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Positive-mass source summary, separate from spatial tree geometry.
 * Zero-initialize before adding particles or child monopoles. */
typedef struct {
    double mass;
    Vec3 center_of_mass;
} GravityMonopole;

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
    GravityMonopole *monopoles;
    bool bounds_fallback;
    bool aggregate_fallback;
    size_t source_node_visits; /* Source nodes visited during group discovery. */
    size_t direct_evaluations; /* Directed particle contributions, excluding self. */
    size_t aggregate_evaluations; /* Aggregate contributions across all targets. */
    BarnesHutTotals totals; /* Cumulative; caller may zero between runs. */
} BarnesHutWorkspace;

/* Add a particle (mass, position) or merge a child (mass, center_of_mass).
 * False means the aggregate is unrepresentable; discard it and fall back. */
/* Incremental weighted means avoid forming mass * position. */
static bool gravity_monopole_add(GravityMonopole *monopole, double mass, Vec3 position)
{
    double total = monopole->mass + mass;
    if (!isfinite(total) || mass <= 0.0) {
        return false;
    }
    if (monopole->mass == 0.0) {
        monopole->center_of_mass = position;
    } else {
        double fraction = mass / total;
        monopole->center_of_mass.x += fraction * (position.x - monopole->center_of_mass.x);
        monopole->center_of_mass.y += fraction * (position.y - monopole->center_of_mass.y);
        monopole->center_of_mass.z += fraction * (position.z - monopole->center_of_mass.z);
    }
    monopole->mass = total;
    return isfinite(monopole->center_of_mass.x) && isfinite(monopole->center_of_mass.y)
        && isfinite(monopole->center_of_mass.z);
}

/* The caller must establish a nonzero, finite distance to this monopole. */
static Vec3 gravity_monopole_acceleration(const GravityMonopole *monopole, Vec3 target)
{
    Vec3 displacement = {monopole->center_of_mass.x - target.x,
                         monopole->center_of_mass.y - target.y,
                         monopole->center_of_mass.z - target.z};
    double distance = hypot(hypot(displacement.x, displacement.y), displacement.z);
    /* Equivalent to M * displacement / r^3 without explicitly cubing r. */
    double scale = (monopole->mass / distance) / distance;
    return (Vec3){scale * (displacement.x / distance),
                  scale * (displacement.y / distance),
                  scale * (displacement.z / distance)};
}

static BarnesHutConfig barnes_hut_default_config(void)
{
    return (BarnesHutConfig){
        .theta = 0.5,
        .bucket_capacity = 8,
        .max_depth = 32,
        .node_capacity = 0
    };
}

static bool aggregate_tree(const SimulationView *state, BarnesHutWorkspace *workspace)
{
    memset(workspace->monopoles, 0,
           workspace->tree.node_count * sizeof(*workspace->monopoles));
    /* Children follow their parents in the pool; visit them first. */
    for (size_t n = workspace->tree.node_count; n-- > 0;) {
        const OctreeNode *node = &workspace->tree.nodes[n];
        GravityMonopole *monopole = &workspace->monopoles[n];
        if (node->is_leaf) {
            for (size_t p = node->begin; p < node->end; ++p) {
                size_t i = workspace->tree.indices[p];
                if (!gravity_monopole_add(monopole, state->w[i],
                                          (Vec3){state->x[i], state->y[i], state->z[i]})) {
                    return false;
                }
            }
        } else {
            for (size_t c = 0; c < 8; ++c) {
                if (node->children[c] != SIZE_MAX) {
                    const GravityMonopole *child = &workspace->monopoles[node->children[c]];
                    if (!gravity_monopole_add(monopole, child->mass, child->center_of_mass)) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
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

static bool accept_source(const OctreeNode *source, const GravityMonopole *monopole,
                          const TargetBounds *bounds, double theta)
{
    if (theta == 0.0) {
        return false;
    }
    const double center[] = {monopole->center_of_mass.x, monopole->center_of_mass.y,
                             monopole->center_of_mass.z};
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
    if (node->is_leaf || (disjoint && accept_source(node, &workspace->monopoles[index],
                                                  bounds, workspace->config.theta))) {
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
                sum = gravity_accumulate_sources(state, target, workspace->tree.indices,
                                                 node->begin, node->end, sum);
                workspace->direct_evaluations += (node->end - node->begin)
                                              - (node->begin <= p && p < node->end);
            } else {
                const GravityMonopole *monopole = &workspace->monopoles[sources[entry]];
                Vec3 contribution = gravity_monopole_acceleration(monopole,
                    (Vec3){state->x[target], state->y[target], state->z[target]});
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
        /* Count directed contributions, even though n2 shares pairs. */
        workspace->direct_evaluations = state->count * (state->count - 1);
        n2_accelerations(state, ax, ay, az);
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
                         double *az, void *solver_context)
{
    BarnesHutWorkspace *workspace = solver_context;
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

static bool barnes_hut_create(void **solver_context, size_t count, const SolverOption *options,
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
    if (workspace->tree.node_capacity > SIZE_MAX / sizeof(*workspace->monopoles) ||
        !(workspace->monopoles = malloc(workspace->tree.node_capacity * sizeof(*workspace->monopoles)))) {
        octree_destroy(&workspace->tree);
        free(workspace);
        *error = "Could not allocate Barnes-Hut workspace.";
        return false;
    }
    workspace->config = config;
    *solver_context = workspace;
    return true;
}

static void barnes_hut_destroy(void *solver_context)
{
    BarnesHutWorkspace *workspace = solver_context;
    octree_destroy(&workspace->tree);
    free(workspace->monopoles);
    free(workspace);
}

static void barnes_hut_settings(const void *solver_context, MetadataVisitor visit, void *visitor_context)
{
    const BarnesHutWorkspace *workspace = solver_context;
    visit("theta", (MetadataValue){METADATA_DOUBLE, {.real = workspace->config.theta}}, visitor_context);
    visit("bucket_capacity", (MetadataValue){METADATA_SIZE, {.size = workspace->config.bucket_capacity}}, visitor_context);
    visit("max_depth", (MetadataValue){METADATA_SIZE, {.size = workspace->config.max_depth}}, visitor_context);
    visit("node_capacity", (MetadataValue){METADATA_SIZE, {.size = workspace->tree.node_capacity}}, visitor_context);
}

static void barnes_hut_diagnostics(const void *solver_context, MetadataVisitor visit, void *visitor_context)
{
    const BarnesHutTotals *totals = &((const BarnesHutWorkspace *)solver_context)->totals;
#define REPORT(field) visit(#field, (MetadataValue){METADATA_SIZE, {.size = totals->field}}, visitor_context)
    REPORT(calls);
    REPORT(capacity_limited_leaves);
    REPORT(bounds_fallbacks);
    REPORT(aggregate_fallbacks);
    /* Historical result field, derived rather than counted twice. */
    visit("accepted_nodes", (MetadataValue){METADATA_SIZE, {.size = totals->aggregate_evaluations}}, visitor_context);
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
