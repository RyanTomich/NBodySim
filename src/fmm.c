#include "gravity.h"
#include "octree.h"
#include "solver.h"
#include "nbody.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

enum { LAPLACE_ORDER = 4, LAPLACE_TERMS = 35, LAPLACE_SPLITS = 210 };

typedef struct {
    int exponent[3];
    int lower[3]; /* Index after lowering the axis exponent, or -1 if zero. */
} CartesianTerm;

/* Multi-indices satisfy a = b + c, with factor choose(a,b) componentwise.
 * For an upward or downward shift, a is the parent term, b the child term,
 * and c the center-offset monomial. The divisor rescales b for the child.
 * In translation, a indexes the kernel, b the target local, and c the source moment. */
typedef struct {
    int a, b, c, factor, child_divisor;
} CartesianSplit;

#include "generated/laplace_order4.h"

typedef struct {
    Octree tree;
    double *moments, *locals;
} FmmWorkspace;

static void fmm_workspace_destroy(void *solver_context)
{
    FmmWorkspace *workspace = solver_context;
    octree_destroy(&workspace->tree);
    free(workspace->moments);
    free(workspace->locals);
    free(workspace);
}

static bool fmm_workspace_create(void **solver_context, size_t count, const SolverOption *options,
                                 size_t option_count, const char **error)
{
    (void)options;
    if (option_count) {
        *error = "FMM uses order four and accepts no solver options.";
        return false;
    }
    *error = "Could not allocate FMM workspace.";
    FmmWorkspace *workspace = calloc(1, sizeof(*workspace));
    if (!workspace) {
        return false;
    }
    /* Allow spare nodes for sparse and clustered trees. Full leaves remain
     * valid even if subdivision exhausts this pool. */
    size_t capacity = count / 2 + 1;
    if (capacity > SIZE_MAX / sizeof(double) / LAPLACE_TERMS ||
        !octree_create(&workspace->tree, count, 32, 32, capacity)) {
        fmm_workspace_destroy(workspace);
        return false;
    }
    size_t bytes = capacity * LAPLACE_TERMS * sizeof(double);
    workspace->moments = malloc(bytes);
    workspace->locals = malloc(bytes);
    if (!workspace->moments || !workspace->locals) {
        fmm_workspace_destroy(workspace);
        return false;
    }
    *solver_context = workspace;
    *error = NULL;
    return true;
}

/* Coordinates in each expansion are scaled by its cell half-width. */
static void evaluate_monomials(Vec3 point, const OctreeNode *cell, double *polynomial)
{
    double offset[] = {(point.x - cell->center.x) / cell->half_width,
                       (point.y - cell->center.y) / cell->half_width,
                       (point.z - cell->center.z) / cell->half_width};
    polynomial[0] = 1;
    for (int a = 1; a < LAPLACE_TERMS; ++a) {
        const int *lower = laplace_terms[a].lower;
        int axis = lower[0] >= 0 ? 0 : (lower[1] >= 0 ? 1 : 2);
        polynomial[a] = offset[axis] * polynomial[lower[axis]];
    }
}

static bool finite_terms(const double *terms, int count)
{
    for (int a = 0; a < count; ++a) {
        if (!isfinite(terms[a])) {
            return false;
        }
    }
    return true;
}

/* Child coordinates have half the scale of their parent's coordinates. */
static double shift_weight(const double *polynomial, CartesianSplit split)
{
    return split.factor * polynomial[split.c] / split.child_divisor;
}

static void shift_moments_up(FmmWorkspace *workspace, const OctreeNode *parent,
                             size_t child, double *parent_moments)
{
    double polynomial[LAPLACE_TERMS];
    evaluate_monomials(workspace->tree.nodes[child].center, parent, polynomial);
    const double *child_moments = workspace->moments + child * LAPLACE_TERMS;
    for (int i = 0; i < LAPLACE_SPLITS; ++i) {
        CartesianSplit split = laplace_splits[i];
        parent_moments[split.a] += shift_weight(polynomial, split) * child_moments[split.b];
    }
}

static void shift_locals_down(FmmWorkspace *workspace, const OctreeNode *parent,
                              size_t child, const double *parent_local)
{
    double polynomial[LAPLACE_TERMS];
    evaluate_monomials(workspace->tree.nodes[child].center, parent, polynomial);
    double *child_local = workspace->locals + child * LAPLACE_TERMS;
    for (int i = 0; i < LAPLACE_SPLITS; ++i) {
        CartesianSplit split = laplace_splits[i];
        child_local[split.b] += shift_weight(polynomial, split) * parent_local[split.a];
    }
}

static bool upward_pass(FmmWorkspace *workspace, const SimulationView *state)
{
    for (size_t n = workspace->tree.node_count; n-- > 0;) {
        const OctreeNode *cell = &workspace->tree.nodes[n];
        double *moments = workspace->moments + n * LAPLACE_TERMS;
        if (cell->is_leaf) {
            double polynomial[LAPLACE_TERMS];
            for (size_t p = cell->begin; p < cell->end; ++p) {
                size_t i = workspace->tree.indices[p];
                evaluate_monomials(simulation_view_position(state, i), cell, polynomial);
                for (int a = 0; a < LAPLACE_TERMS; ++a) {
                    moments[a] += state->w[i] * polynomial[a];
                }
            }
        } else {
            for (int c = 0; c < 8; ++c) {
                if (cell->children[c] != SIZE_MAX) {
                    shift_moments_up(workspace, cell, cell->children[c], moments);
                }
            }
        }
        if (!finite_terms(moments, LAPLACE_TERMS)) {
            return false;
        }
    }
    return true;
}

typedef struct {
    const SimulationView *state;
    double *acceleration[3];
} Field;

static void evaluate_direct_range(const FmmWorkspace *workspace, const Field *field,
                                  const OctreeNode *source, size_t begin, size_t end)
{
    for (size_t p = begin; p < end; ++p) {
        size_t i = workspace->tree.indices[p];
        Vec3 sum = {field->acceleration[0][i], field->acceleration[1][i],
                    field->acceleration[2][i]};
        sum = gravity_accumulate_sources(field->state, i, workspace->tree.indices,
                                         source->begin, source->end, sum);
        field->acceleration[0][i] = sum.x;
        field->acceleration[1][i] = sum.y;
        field->acceleration[2][i] = sum.z;
    }
}

/* Add a fourth-order source expansion to the target local polynomial. */
static void laplace_translate(const double *moments,
                              double source_half_width, double target_half_width,
                              const double unit[3], double r, double *local)
{
    double kernel[LAPLACE_TERMS];
    double source_scale[LAPLACE_ORDER + 1] = {1}, target_scale[LAPLACE_ORDER + 1] = {1};
    laplace_kernel_order4(kernel, unit);
    for (int d = 1; d <= LAPLACE_ORDER; ++d) {
        source_scale[d] = source_scale[d - 1] * (-source_half_width / r);
        target_scale[d] = target_scale[d - 1] * (target_half_width / r);
    }
    laplace_contract_order4(local, moments, kernel, source_scale, target_scale, r);
}

static bool interact_cell_pair(const FmmWorkspace *workspace, const Field *field,
                               size_t source_index, size_t target_index, size_t grain)
{
    const OctreeNode *source_cell = &workspace->tree.nodes[source_index];
    const OctreeNode *target_cell = &workspace->tree.nodes[target_index];
    double dx = target_cell->center.x - source_cell->center.x;
    double dy = target_cell->center.y - source_cell->center.y;
    double dz = target_cell->center.z - source_cell->center.z;
    double r = hypot(hypot(dx, dy), dz);
    if (!isfinite(r)) {
        return false;
    }
    /* The sum of the cells' enclosing sphere radii must be at most r / 2. */
    if (source_cell->half_width + target_cell->half_width <= r * (0.5 / sqrt(3.0))) {
        double unit[] = {dx / r, dy / r, dz / r};
        laplace_translate(workspace->moments + source_index * LAPLACE_TERMS,
                          source_cell->half_width, target_cell->half_width, unit, r,
                          workspace->locals + target_index * LAPLACE_TERMS);
        return true;
    }
    if (source_cell->is_leaf && target_cell->is_leaf) {
#ifdef _OPENMP
        if (target_cell->end - target_cell->begin > grain) {
            for (size_t begin = target_cell->begin; begin < target_cell->end; begin += grain) {
                size_t end = begin + grain < target_cell->end ? begin + grain : target_cell->end;
#pragma omp task firstprivate(workspace, field, source_cell, begin, end)
                evaluate_direct_range(workspace, field, source_cell, begin, end);
            }
#pragma omp taskwait
            return true;
        }
#endif
        evaluate_direct_range(workspace, field, source_cell, target_cell->begin, target_cell->end);
        return true;
    }
    bool split_target = source_cell->is_leaf ||
                        (!target_cell->is_leaf && target_cell->half_width >= source_cell->half_width);
#ifdef _OPENMP
    if (split_target && target_cell->end - target_cell->begin > grain) {
        bool failed[8] = {false};
        for (int c = 0; c < 8; ++c) {
            size_t child = target_cell->children[c];
            if (child == SIZE_MAX) {
                continue;
            }
            const OctreeNode *node = &workspace->tree.nodes[child];
            if (node->begin == target_cell->begin && node->end == target_cell->end) {
                /* Only child. */
                return interact_cell_pair(workspace, field, source_index, child, grain);
            }
            /* Sibling targets own disjoint locals and acceleration entries. */
#pragma omp task firstprivate(workspace, field, source_index, child, grain, c) shared(failed)
            failed[c] = !interact_cell_pair(workspace, field, source_index, child, grain);
        }
#pragma omp taskwait
        for (int c = 0; c < 8; ++c) {
            if (failed[c]) {
                return false;
            }
        }
        return true;
    }
#endif
    /* Source siblings accumulate into the same target in octant order. */
    const OctreeNode *split = split_target ? target_cell : source_cell;
    for (int c = 0; c < 8; ++c) {
        size_t child = split->children[c];
        if (child == SIZE_MAX) {
            continue;
        }
        size_t next_source = split_target ? source_index : child;
        size_t next_target = split_target ? child : target_index;
        if (!interact_cell_pair(workspace, field, next_source, next_target, grain)) {
            return false;
        }
    }
    return true;
}

static bool interaction_pass(const FmmWorkspace *workspace, const Field *field)
{
#ifdef _OPENMP
    size_t count = field->state->count;
    if (count >= 1024 && omp_get_max_threads() > 1) {
        bool success;
#pragma omp parallel
        {
#pragma omp single
            {
                size_t workers = (size_t)omp_get_num_threads();
                /* Aim for four pieces per worker, with at least 64 targets per piece. */
                size_t grain = 1 + (count - 1) / workers / 4;
                if (grain < 64) {
                    grain = 64;
                }
                success = interact_cell_pair(workspace, field, 0, 0, workers > 1 ? grain : SIZE_MAX);
            }
        }
        return success;
    }
#endif
    return interact_cell_pair(workspace, field, 0, 0, SIZE_MAX);
}

static bool downward_pass(FmmWorkspace *workspace, Field *field)
{
    for (size_t n = 0; n < workspace->tree.node_count; ++n) {
        const OctreeNode *cell = &workspace->tree.nodes[n];
        const double *local = workspace->locals + n * LAPLACE_TERMS;
        if (!finite_terms(local, LAPLACE_TERMS)) {
            return false;
        }
        if (!cell->is_leaf) {
            for (int c = 0; c < 8; ++c) {
                if (cell->children[c] != SIZE_MAX) {
                    shift_locals_down(workspace, cell, cell->children[c], local);
                }
            }
            continue;
        }
        double polynomial[LAPLACE_TERMS];
        for (size_t p = cell->begin; p < cell->end; ++p) {
            size_t i = workspace->tree.indices[p];
            evaluate_monomials(simulation_view_position(field->state, i), cell, polynomial);
            for (int axis = 0; axis < 3; ++axis) {
                double gradient_component = 0;
                for (int a = 1; a < LAPLACE_TERMS; ++a) {
                    int b = laplace_terms[a].lower[axis];
                    if (b >= 0) {
                        gradient_component += local[a] * laplace_terms[a].exponent[axis]
                                              * polynomial[b];
                    }
                }
                /* Positive gradient of the positive 1/r expansion is gravity.
                 * The gravitational field is already target acceleration. */
                field->acceleration[axis][i] += gradient_component / cell->half_width;
                if (!isfinite(field->acceleration[axis][i])) {
                    return false;
                }
            }
        }
    }
    return true;
}

static const char *fmm_accelerations(const SimulationView *state, double *ax, double *ay, double *az,
                                     void *solver_context)
{
    FmmWorkspace *workspace = solver_context;
    if (!state->count) {
        return NULL;
    }
    if (!octree_build(state, &workspace->tree)) {
        return "fmm: cannot represent tree bounds.";
    }
    size_t bytes = workspace->tree.node_count * LAPLACE_TERMS * sizeof(double);
    memset(workspace->moments, 0, bytes);
    memset(workspace->locals, 0, bytes);
    Field field = {state, {ax, ay, az}};
    for (int axis = 0; axis < 3; ++axis) {
        memset(field.acceleration[axis], 0, state->count * sizeof(double));
    }
    if (!upward_pass(workspace, state)) {
        return "fmm: nonfinite multipole moments in upward pass.";
    }
    if (!interaction_pass(workspace, &field)) {
        return "fmm: nonfinite cell distance in interaction pass.";
    }
    if (!downward_pass(workspace, &field)) {
        return "fmm: nonfinite local expansion or acceleration in downward pass.";
    }
    return NULL;
}

static void fmm_settings(const void *solver_context, MetadataVisitor visit, void *visitor_context)
{
    const FmmWorkspace *workspace = solver_context;
    const Octree *tree = &workspace->tree;
    visit("order", (MetadataValue){METADATA_SIZE, {.size = LAPLACE_ORDER}}, visitor_context);
    visit("bucket_capacity", (MetadataValue){METADATA_SIZE, {.size = tree->bucket_capacity}},
          visitor_context);
    visit("max_depth", (MetadataValue){METADATA_SIZE, {.size = tree->max_depth}}, visitor_context);
    visit("node_capacity", (MetadataValue){METADATA_SIZE, {.size = tree->node_capacity}},
          visitor_context);
}

const SolverDescriptor fmm_solver = {
    .name = "fmm",
    .acceleration = fmm_accelerations,
    .create = fmm_workspace_create,
    .destroy = fmm_workspace_destroy,
    .settings = fmm_settings
};
