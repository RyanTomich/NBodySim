#include "laplace_internal.h"
#include "gravity_internal.h"
#include "n2_internal.h"
#include "octree_internal.h"
#include "solver_internal.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    Octree tree;
    double *moments, *locals;
    size_t direct_fallbacks;
} Fmm;

static void destroy(void *solver_context)
{
    Fmm *f = solver_context;
    octree_destroy(&f->tree);
    free(f->moments);
    free(f->locals);
    free(f);
}

static bool create(void **solver_context, size_t count, const SolverOption *options, size_t option_count,
                   const char **error)
{
    (void)options;
    if (option_count) {
        *error = "FMM uses order four and accepts no solver options.";
        return false;
    }
    *error = "Could not allocate FMM workspace.";
    Fmm *f = calloc(1, sizeof(*f));
    if (!f) {
        return false;
    }
    /* Allow spare nodes for sparse and clustered trees. Full leaves remain
     * valid even if subdivision exhausts this pool. */
    size_t capacity = count / 2 + 1;
    if (capacity > SIZE_MAX / sizeof(double) / LAPLACE_TERMS ||
        !octree_create(&f->tree, count, 32, 32, capacity)) {
        destroy(f);
        return false;
    }
    size_t bytes = capacity * LAPLACE_TERMS * sizeof(double);
    f->moments = malloc(bytes);
    f->locals = malloc(bytes);
    if (!f->moments || !f->locals) {
        destroy(f);
        return false;
    }
    *solver_context = f;
    *error = NULL;
    return true;
}

static Vec3 position(const SimulationView *state, size_t i)
{
    return (Vec3){state->x[i], state->y[i], state->z[i]};
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

static void shift_moments_up(Fmm *f, const OctreeNode *parent, size_t child, double *parent_moments)
{
    double polynomial[LAPLACE_TERMS];
    evaluate_monomials(f->tree.nodes[child].center, parent, polynomial);
    const double *child_moments = f->moments + child * LAPLACE_TERMS;
    for (int i = 0; i < LAPLACE_SPLITS; ++i) {
        CartesianSplit split = laplace_splits[i];
        parent_moments[split.a] += shift_weight(polynomial, split) * child_moments[split.b];
    }
}

static void shift_locals_down(Fmm *f, const OctreeNode *parent, size_t child, const double *parent_local)
{
    double polynomial[LAPLACE_TERMS];
    evaluate_monomials(f->tree.nodes[child].center, parent, polynomial);
    double *child_local = f->locals + child * LAPLACE_TERMS;
    for (int i = 0; i < LAPLACE_SPLITS; ++i) {
        CartesianSplit split = laplace_splits[i];
        child_local[split.b] += shift_weight(polynomial, split) * parent_local[split.a];
    }
}

static bool upward_pass(Fmm *f, const SimulationView *state)
{
    for (size_t n = f->tree.node_count; n-- > 0;) {
        const OctreeNode *cell = &f->tree.nodes[n];
        double *moments = f->moments + n * LAPLACE_TERMS;
        if (cell->is_leaf) {
            double polynomial[LAPLACE_TERMS];
            for (size_t p = cell->begin; p < cell->end; ++p) {
                size_t i = f->tree.indices[p];
                evaluate_monomials(position(state, i), cell, polynomial);
                for (int a = 0; a < LAPLACE_TERMS; ++a) {
                    moments[a] += state->w[i] * polynomial[a];
                }
            }
        } else {
            for (int c = 0; c < 8; ++c) {
                if (cell->children[c] != SIZE_MAX) {
                    shift_moments_up(f, cell, cell->children[c], moments);
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

static void evaluate_direct_pair(const Fmm *f, const Field *field,
                                 const OctreeNode *source, const OctreeNode *target)
{
    for (size_t p = target->begin; p < target->end; ++p) {
        size_t i = f->tree.indices[p];
        Vec3 sum = {field->acceleration[0][i], field->acceleration[1][i],
                    field->acceleration[2][i]};
        sum = gravity_accumulate_sources(field->state, i, f->tree.indices,
                                         source->begin, source->end, sum);
        field->acceleration[0][i] = sum.x;
        field->acceleration[1][i] = sum.y;
        field->acceleration[2][i] = sum.z;
    }
}

static bool interact(const Fmm *f, const Field *field, size_t source, size_t target)
{
    const OctreeNode *s = &f->tree.nodes[source], *t = &f->tree.nodes[target];
    double dx = t->center.x - s->center.x, dy = t->center.y - s->center.y;
    double dz = t->center.z - s->center.z;
    double r = hypot(hypot(dx, dy), dz);
    if (!isfinite(r)) {
        return false;
    }
    /* The sum of the cells' enclosing sphere radii must be at most r / 2. */
    if (s->half_width + t->half_width <= r * (0.5 / sqrt(3.0))) {
        double unit[] = {dx / r, dy / r, dz / r};
        laplace_translate(f->moments + source * LAPLACE_TERMS,
                          s->half_width, t->half_width, unit, r,
                          f->locals + target * LAPLACE_TERMS);
        return true;
    }
    if (s->is_leaf && t->is_leaf) {
        evaluate_direct_pair(f, field, s, t);
        return true;
    }
    bool split_target = s->is_leaf || (!t->is_leaf && t->half_width >= s->half_width);
    const OctreeNode *split = split_target ? t : s;
    for (int c = 0; c < 8; ++c) {
        size_t child = split->children[c];
        if (child != SIZE_MAX &&
            !interact(f, field, split_target ? source : child, split_target ? child : target)) {
            return false;
        }
    }
    return true;
}

static bool interaction_pass(const Fmm *f, const Field *field)
{
#ifdef _OPENMP
    int workers = omp_get_max_threads();
    if (workers < 2 || field->state->count < 1024) {
        return interact(f, field, 0, 0);
    }
    const OctreeNode *root = &f->tree.nodes[0];
    size_t targets[8];
    int jobs = 0;
    for (int c = 0; c < 8; ++c) {
        if (root->children[c] != SIZE_MAX) {
            targets[jobs++] = root->children[c];
        }
    }
    if (jobs < 2) {
        return interact(f, field, 0, 0);
    }
    if (workers > jobs) {
        workers = jobs;
    }
    /* Root-to-root traversal first splits the target. Its child subtrees
     * have disjoint outputs and each retains the original source order. */
    int failed = 0;
#pragma omp parallel for schedule(dynamic, 1) num_threads(workers) reduction(|:failed)
    for (int i = 0; i < jobs; ++i) {
        failed |= !interact(f, field, 0, targets[i]);
    }
    return !failed;
#else
    return interact(f, field, 0, 0);
#endif
}

static bool downward_pass(Fmm *f, Field *field)
{
    for (size_t n = 0; n < f->tree.node_count; ++n) {
        const OctreeNode *cell = &f->tree.nodes[n];
        const double *local = f->locals + n * LAPLACE_TERMS;
        if (!finite_terms(local, LAPLACE_TERMS)) {
            return false;
        }
        if (!cell->is_leaf) {
            for (int c = 0; c < 8; ++c) {
                if (cell->children[c] != SIZE_MAX) {
                    shift_locals_down(f, cell, cell->children[c], local);
                }
            }
            continue;
        }
        double polynomial[LAPLACE_TERMS];
        for (size_t p = cell->begin; p < cell->end; ++p) {
            size_t i = f->tree.indices[p];
            evaluate_monomials(position(field->state, i), cell, polynomial);
            for (int axis = 0; axis < 3; ++axis) {
                double g = 0;
                for (int a = 1; a < LAPLACE_TERMS; ++a) {
                    int b = laplace_terms[a].lower[axis];
                    if (b >= 0) {
                        g += local[a] * laplace_terms[a].exponent[axis] * polynomial[b];
                    }
                }
                /* Positive gradient of the positive 1/r expansion is gravity.
                 * The gravitational field is already target acceleration. */
                field->acceleration[axis][i] += g / cell->half_width;
                if (!isfinite(field->acceleration[axis][i])) {
                    return false;
                }
            }
        }
    }
    return true;
}

static void evaluate(const SimulationView *state, double *ax, double *ay, double *az, void *solver_context)
{
    Fmm *f = solver_context;
    if (!state->count) {
        return;
    }
    Field field = {state, {ax, ay, az}};
    if (octree_build(state, &f->tree)) {
        size_t bytes = f->tree.node_count * LAPLACE_TERMS * sizeof(double);
        memset(f->moments, 0, bytes);
        memset(f->locals, 0, bytes);
        for (int axis = 0; axis < 3; ++axis) {
            memset(field.acceleration[axis], 0, state->count * sizeof(double));
        }
        if (upward_pass(f, state) && interaction_pass(f, &field) && downward_pass(f, &field)) {
            return;
        }
    }
    ++f->direct_fallbacks;
    n2_accelerations(state, ax, ay, az);
}

static void settings(const void *solver_context, MetadataVisitor visit, void *visitor_context)
{
    const Fmm *f = solver_context;
    visit("order", (MetadataValue){METADATA_SIZE, {.size = LAPLACE_ORDER}}, visitor_context);
    visit("bucket_capacity", (MetadataValue){METADATA_SIZE, {.size = f->tree.bucket_capacity}},
          visitor_context);
    visit("max_depth", (MetadataValue){METADATA_SIZE, {.size = f->tree.max_depth}}, visitor_context);
    visit("node_capacity", (MetadataValue){METADATA_SIZE, {.size = f->tree.node_capacity}}, visitor_context);
}

static void diagnostics(const void *solver_context, MetadataVisitor visit, void *visitor_context)
{
    const Fmm *f = solver_context;
    visit("direct_fallbacks", (MetadataValue){METADATA_SIZE, {.size = f->direct_fallbacks}}, visitor_context);
}

const SolverDescriptor fmm_solver = {.name = "fmm",
                                     .acceleration = evaluate,
                                     .create = create,
                                     .destroy = destroy,
                                     .settings = settings,
                                     .diagnostics = diagnostics};
