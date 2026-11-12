#define _POSIX_C_SOURCE 200809L

#include "parse.h"
#include "solver_internal.h"
#include "octree_internal.h"
#include "gravity_internal.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(NBODY_REQUIRE_OPENMP) && !defined(_OPENMP)
#error "OPENMP=1 requires OpenMP compiler support; check OMPFLAGS and OMPLIBS."
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

/* Experiment-only instrumentation: no timers or extra metadata in normal builds.
 * Worker timers are private. After joining, divide their sum by the actual
 * team size to attribute elapsed time without double-counting overlap. */
#ifdef NBODY_PROFILE_BH
#include <time.h>
static uint64_t profile_now(void)
{
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) {
        abort();
    }
    return (uint64_t)now.tv_sec * UINT64_C(1000000000) + (uint64_t)now.tv_nsec;
}
#define PROFILE_START(name) uint64_t name = profile_now()
#define PROFILE_ADD(destination, start) ((destination) += (double)(profile_now() - (start)) * 1e-9)
#else
#define PROFILE_START(name) ((void)0)
#define PROFILE_ADD(destination, start) ((void)0)
#endif

/* Positive-mass source summary, separate from spatial tree geometry.
 * Zero-initialize before adding particles or child monopoles. */
typedef struct {
    double mass;
    Vec3 center_of_mass;
} GravityMonopole;

/* Requested settings, used only during preparation. */
typedef struct {
    double theta; /* Zero opens every internal node. */
    size_t bucket_capacity;
    size_t max_depth;
    size_t node_capacity; /* Zero selects max(1, 4 * count). */
} BarnesHutConfig;

typedef struct {
    size_t calls, capacity_limited_leaves;
    size_t source_node_visits, direct_evaluations, aggregate_evaluations;
#ifdef NBODY_PROFILE_BH
    size_t membership_checks, tree_builds, tree_reuses;
    size_t worker_slots, min_threads, max_threads;
    double membership_seconds, construction_seconds, summaries_seconds;
    double bounds_seconds, traversal_seconds, force_seconds, evaluation_seconds;
    double bounds_worker_seconds, traversal_worker_seconds, force_worker_seconds;
#endif
} BarnesHutTotals;

typedef struct {
    /* Node indices shared by all targets in one leaf: leaves are evaluated
     * directly, accepted internal nodes through their monopoles. */
    size_t *source_nodes; /* Borrowed slice of workspace source storage. */
    size_t source_count;
    size_t source_node_visits, direct_evaluations, aggregate_evaluations;
#ifdef NBODY_PROFILE_BH
    double bounds_seconds, traversal_seconds, force_seconds;
#endif
} BarnesHutWorker;

typedef struct {
    double theta;
    Octree tree;
    GravityMonopole *monopoles;
    BarnesHutWorker *workers;
    size_t *source_storage;
    int worker_count; /* Fixed at preparation, including for later runtime changes. */
    bool topology_usable;
    BarnesHutTotals totals; /* Cumulative over the prepared solver's lifetime. */
} BarnesHutWorkspace;

/* Add a particle (mass, position) or merge a child (mass, center_of_mass).
 * False means the aggregate is unrepresentable. */
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
    double r2 = displacement.x * displacement.x
              + displacement.y * displacement.y
              + displacement.z * displacement.z;
    double scale = monopole->mass / (r2 * sqrt(r2));
    return (Vec3){scale * displacement.x,
                  scale * displacement.y,
                  scale * displacement.z};
}

static bool ensure_topology(const SimulationView *state, BarnesHutWorkspace *workspace)
{
    bool reusable = workspace->topology_usable;
    if (reusable) {
        PROFILE_START(check_start);
        reusable = octree_membership_is_valid(state, &workspace->tree);
        PROFILE_ADD(workspace->totals.membership_seconds, check_start);
#ifdef NBODY_PROFILE_BH
        ++workspace->totals.membership_checks;
#endif
    }
    if (!reusable) {
        PROFILE_START(build_start);
        workspace->topology_usable = octree_build(state, &workspace->tree);
        PROFILE_ADD(workspace->totals.construction_seconds, build_start);
#ifdef NBODY_PROFILE_BH
        ++workspace->totals.tree_builds;
#endif
    }
#ifdef NBODY_PROFILE_BH
    else {
        ++workspace->totals.tree_reuses;
    }
#endif
    return workspace->topology_usable;
}

/* Only topology is cached: refresh summaries from current masses and positions.
 * Failure leaves the topology available for a later call with valid summaries. */
static bool refresh_monopoles(const SimulationView *state, BarnesHutWorkspace *workspace)
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
} TargetLeafBounds;

static TargetLeafBounds compute_target_leaf_bounds(const SimulationView *state,
                                                   const BarnesHutWorkspace *workspace,
                                                   const OctreeNode *target_leaf)
{
    const double *coordinates[] = {state->x, state->y, state->z};
    TargetLeafBounds bounds;
    for (size_t axis = 0; axis < 3; ++axis) {
        double low = coordinates[axis][workspace->tree.indices[target_leaf->begin]];
        double high = low;
        for (size_t p = target_leaf->begin + 1; p < target_leaf->end; ++p) {
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
                          const TargetLeafBounds *bounds, double theta)
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

static void collect_source_nodes(const BarnesHutWorkspace *workspace, BarnesHutWorker *worker,
                                 size_t index, const OctreeNode *target_leaf,
                                 const TargetLeafBounds *bounds)
{
    const OctreeNode *node = &workspace->tree.nodes[index];
    ++worker->source_node_visits;
    bool disjoint = node->end <= target_leaf->begin || target_leaf->end <= node->begin;
    if (node->is_leaf || (disjoint && accept_source(node, &workspace->monopoles[index],
                                                  bounds, workspace->theta))) {
        /* Disjoint, nonempty source ranges bound the list length by count. */
        assert(worker->source_count < workspace->tree.count);
        worker->source_nodes[worker->source_count++] = index;
        size_t targets = target_leaf->end - target_leaf->begin;
        if (node->is_leaf) {
            /* Source leaves are identical to or disjoint from the target leaf. */
            worker->direct_evaluations += targets * (node->end - node->begin)
                                       - (node == target_leaf ? targets : 0);
        } else {
            worker->aggregate_evaluations += targets;
        }
        return;
    }
    for (size_t c = 0; c < 8; ++c) {
        if (node->children[c] != SIZE_MAX) {
            collect_source_nodes(workspace, worker, node->children[c], target_leaf, bounds);
        }
    }
}

static void evaluate_target_leaf(const SimulationView *state, double *ax, double *ay,
                                 double *az, const BarnesHutWorkspace *workspace,
                                 const BarnesHutWorker *worker,
                                 const OctreeNode *target_leaf)
{
    const size_t *source_nodes = worker->source_nodes;
    for (size_t p = target_leaf->begin; p < target_leaf->end; ++p) {
        size_t target = workspace->tree.indices[p];
        Vec3 sum = {0};
        for (size_t entry = 0; entry < worker->source_count; ++entry) {
            const OctreeNode *node = &workspace->tree.nodes[source_nodes[entry]];
            if (node->is_leaf) {
                sum = gravity_accumulate_sources(state, target, workspace->tree.indices,
                                                 node->begin, node->end, sum);
            } else {
                const GravityMonopole *monopole = &workspace->monopoles[source_nodes[entry]];
                Vec3 contribution = gravity_monopole_acceleration(monopole,
                    (Vec3){state->x[target], state->y[target], state->z[target]});
                sum.x += contribution.x;
                sum.y += contribution.y;
                sum.z += contribution.z;
            }
        }
        ax[target] = sum.x;
        ay[target] = sum.y;
        az[target] = sum.z;
    }
}

static void accumulate_worker_totals(BarnesHutWorkspace *workspace, int actual_workers)
{
    BarnesHutTotals *totals = &workspace->totals;
    assert(actual_workers > 0 && actual_workers <= workspace->worker_count);
#ifdef NBODY_PROFILE_BH
    size_t team_size = (size_t)actual_workers;
    totals->worker_slots += team_size;
    if (totals->min_threads == 0 || team_size < totals->min_threads) {
        totals->min_threads = team_size;
    }
    if (team_size > totals->max_threads) {
        totals->max_threads = team_size;
    }
#endif
    for (int i = 0; i < actual_workers; ++i) {
        const BarnesHutWorker *worker = &workspace->workers[i];
        totals->source_node_visits += worker->source_node_visits;
        totals->direct_evaluations += worker->direct_evaluations;
        totals->aggregate_evaluations += worker->aggregate_evaluations;
#ifdef NBODY_PROFILE_BH
        /* Keep raw worker time and normalize each call by its actual team size.
         * The enclosing wall-time residual includes idle and scheduling time. */
        totals->bounds_worker_seconds += worker->bounds_seconds;
        totals->traversal_worker_seconds += worker->traversal_seconds;
        totals->force_worker_seconds += worker->force_seconds;
        totals->bounds_seconds += worker->bounds_seconds / actual_workers;
        totals->traversal_seconds += worker->traversal_seconds / actual_workers;
        totals->force_seconds += worker->force_seconds / actual_workers;
#endif
    }
}

static void evaluate_target_leaves(const SimulationView *state, double *ax, double *ay,
                                   double *az, BarnesHutWorkspace *workspace)
{
    int actual_workers = 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(workspace->worker_count) default(none) shared(state, ax, ay, az, workspace, actual_workers)
#endif
    {
        int worker_index = 0;
#ifdef _OPENMP
        worker_index = omp_get_thread_num();
        if (worker_index == 0) {
            actual_workers = omp_get_num_threads();
        }
#endif
        /* Keep hot counters local so workers do not repeatedly share cache lines. */
        BarnesHutWorker worker = {.source_nodes = workspace->workers[worker_index].source_nodes};
        /* Each leaf owns a disjoint target range; source order stays serial. */
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (size_t n = 0; n < workspace->tree.node_count; ++n) {
            const OctreeNode *target_leaf = &workspace->tree.nodes[n];
            if (target_leaf->is_leaf) {
                PROFILE_START(bounds_start);
                TargetLeafBounds bounds = compute_target_leaf_bounds(state, workspace, target_leaf);
                PROFILE_ADD(worker.bounds_seconds, bounds_start);
                worker.source_count = 0;
                PROFILE_START(traversal_start);
                collect_source_nodes(workspace, &worker, 0, target_leaf, &bounds);
                PROFILE_ADD(worker.traversal_seconds, traversal_start);
                PROFILE_START(force_start);
                evaluate_target_leaf(state, ax, ay, az, workspace, &worker, target_leaf);
                PROFILE_ADD(worker.force_seconds, force_start);
            }
        }
        workspace->workers[worker_index] = worker;
    }
    accumulate_worker_totals(workspace, actual_workers);
}

static void invalidate_accelerations(size_t count, double *ax, double *ay, double *az)
{
    /* The callback has no error return; publish an invalid result. */
    for (size_t i = 0; i < count; ++i) {
        ax[i] = ay[i] = az[i] = NAN;
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
    ++workspace->totals.calls;
    bool valid = ensure_topology(state, workspace);
    if (valid) {
        PROFILE_START(summaries_start);
        valid = refresh_monopoles(state, workspace);
        PROFILE_ADD(workspace->totals.summaries_seconds, summaries_start);
    }
    if (valid) {
        PROFILE_START(evaluation_start);
        evaluate_target_leaves(state, ax, ay, az, workspace);
        PROFILE_ADD(workspace->totals.evaluation_seconds, evaluation_start);
    } else {
        invalidate_accelerations(state->count, ax, ay, az);
    }
    /* Count the current topology's limited leaves on every evaluation,
     * including evaluations which reuse it. */
    workspace->totals.capacity_limited_leaves += workspace->tree.capacity_limited_leaves;
}

static void barnes_hut_destroy(void *solver_context)
{
    BarnesHutWorkspace *workspace = solver_context;
    if (!workspace) {
        return;
    }
    free(workspace->source_storage);
    free(workspace->workers);
    octree_destroy(&workspace->tree);
    free(workspace->monopoles);
    free(workspace);
}

static bool barnes_hut_parse_config(BarnesHutConfig *config, const SolverOption *options,
                                    size_t option_count, const char **error)
{
    *config = (BarnesHutConfig){
        .theta = 0.5,
        .bucket_capacity = 8,
        .max_depth = 32,
        .node_capacity = 0
    };
    for (size_t i = 0; i < option_count; ++i) {
        const char *name = options[i].name, *value = options[i].value;
        bool valid = false;
        if (strcmp(name, "theta") == 0) {
            valid = parse_double(value, &config->theta) && config->theta >= 0;
        } else if (strcmp(name, "bucket_capacity") == 0) {
            valid = parse_size(value, &config->bucket_capacity) && config->bucket_capacity > 0;
        } else if (strcmp(name, "max_depth") == 0) {
            valid = parse_size(value, &config->max_depth) && config->max_depth <= OCTREE_MAX_DEPTH;
        } else if (strcmp(name, "node_capacity") == 0) {
            valid = parse_size(value, &config->node_capacity);
        }
        if (!valid) {
            *error = "Unknown solver option or invalid value.";
            return false;
        }
    }
    return true;
}

static bool barnes_hut_create(void **solver_context, size_t count, const SolverOption *options,
                              size_t option_count, const char **error)
{
    BarnesHutConfig config;
    if (!barnes_hut_parse_config(&config, options, option_count, error)) {
        return false;
    }
    BarnesHutWorkspace *workspace = calloc(1, sizeof(*workspace));
    if (!workspace || !octree_create(&workspace->tree, count, config.bucket_capacity,
                                     config.max_depth, config.node_capacity)) {
        goto allocation_failure;
    }
    if (workspace->tree.node_capacity > SIZE_MAX / sizeof(*workspace->monopoles) ||
        !(workspace->monopoles = malloc(workspace->tree.node_capacity * sizeof(*workspace->monopoles)))) {
        goto allocation_failure;
    }
    workspace->worker_count = 1;
#ifdef _OPENMP
    workspace->worker_count = omp_get_max_threads();
#endif
    if (workspace->worker_count <= 0 ||
        (size_t)workspace->worker_count > SIZE_MAX / sizeof(*workspace->workers) ||
        count > SIZE_MAX / sizeof(*workspace->source_storage) / (size_t)workspace->worker_count) {
        goto allocation_failure;
    }
    workspace->workers = calloc((size_t)workspace->worker_count, sizeof(*workspace->workers));
    if (!workspace->workers) {
        goto allocation_failure;
    }
    if (count != 0) {
        workspace->source_storage = malloc((size_t)workspace->worker_count * count
                                            * sizeof(*workspace->source_storage));
        if (!workspace->source_storage) {
            goto allocation_failure;
        }
        for (int i = 0; i < workspace->worker_count; ++i) {
            workspace->workers[i].source_nodes = workspace->source_storage + (size_t)i * count;
        }
    }
    workspace->theta = config.theta;
    *solver_context = workspace;
    return true;

allocation_failure:
    barnes_hut_destroy(workspace);
    *error = "Could not allocate Barnes-Hut workspace.";
    return false;
}

static void barnes_hut_settings(const void *solver_context, MetadataVisitor visit, void *visitor_context)
{
    const BarnesHutWorkspace *workspace = solver_context;
    visit("theta", (MetadataValue){METADATA_DOUBLE, {.real = workspace->theta}}, visitor_context);
    visit("bucket_capacity", (MetadataValue){METADATA_SIZE, {.size = workspace->tree.bucket_capacity}}, visitor_context);
    visit("max_depth", (MetadataValue){METADATA_SIZE, {.size = workspace->tree.max_depth}}, visitor_context);
    visit("node_capacity", (MetadataValue){METADATA_SIZE, {.size = workspace->tree.node_capacity}}, visitor_context);
}

static void barnes_hut_diagnostics(const void *solver_context, MetadataVisitor visit, void *visitor_context)
{
    const BarnesHutTotals *totals = &((const BarnesHutWorkspace *)solver_context)->totals;
#define REPORT(field) visit(#field, (MetadataValue){METADATA_SIZE, {.size = totals->field}}, visitor_context)
    REPORT(calls);
    REPORT(capacity_limited_leaves);
    REPORT(source_node_visits);
    REPORT(direct_evaluations);
    REPORT(aggregate_evaluations);
#ifdef NBODY_PROFILE_BH
    REPORT(membership_checks);
    REPORT(tree_builds);
    REPORT(tree_reuses);
    REPORT(worker_slots);
    REPORT(min_threads);
    REPORT(max_threads);
#endif
#undef REPORT
#ifdef NBODY_PROFILE_BH
#define REPORT_TIME(field) visit(#field, (MetadataValue){METADATA_DOUBLE, {.real = totals->field}}, visitor_context)
    REPORT_TIME(membership_seconds);
    REPORT_TIME(construction_seconds);
    REPORT_TIME(summaries_seconds);
    REPORT_TIME(bounds_seconds);
    REPORT_TIME(traversal_seconds);
    REPORT_TIME(force_seconds);
    REPORT_TIME(evaluation_seconds);
    REPORT_TIME(bounds_worker_seconds);
    REPORT_TIME(traversal_worker_seconds);
    REPORT_TIME(force_worker_seconds);
#undef REPORT_TIME
#endif
}

const SolverDescriptor barnes_hut_solver = {
    .name = "barnes-hut",
    .acceleration = gravity_barnes_hut,
    .create = barnes_hut_create,
    .destroy = barnes_hut_destroy,
    .settings = barnes_hut_settings,
    .diagnostics = barnes_hut_diagnostics
};
