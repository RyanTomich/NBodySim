#include "laplace_internal.h"
#include "gravity_internal.h"
#include "n2_internal.h"
#include "octree_internal.h"
#include "parse.h"
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
    size_t direct_pair_cutoff;
    CartesianTerms expansion;
    double *moments, *locals;
    /* Per-node occupied extent about the cell center, refreshed in the upward pass. */
    double *occupied_half_width;
    double polynomial[CARTESIAN_MAX_TERMS];
    size_t calls, translations, direct_evaluations, capacity_limited_leaves, direct_fallbacks;
} Fmm;

static void destroy(void *solver_context)
{
    Fmm *f = solver_context;
    octree_destroy(&f->tree);
    free(f->moments);
    free(f->locals);
    free(f->occupied_half_width);
    free(f);
}

static bool create(void **solver_context, size_t count, const SolverOption *options, size_t option_count,
                   const char **error)
{
    size_t order = 4;
    for (size_t i = 0; i < option_count; ++i) {
        if (strcmp(options[i].name, "order") || !parse_size(options[i].value, &order) ||
            order < 1 || order > CARTESIAN_MAX_ORDER) {
            *error = "FMM accepts only order, an integer from 1 through 6.";
            return false;
        }
    }
    *error = "Could not allocate FMM workspace.";
    Fmm *f = calloc(1, sizeof(*f));
    if (!f) {
        return false;
    }
    f->expansion.order = order;
    f->direct_pair_cutoff = order == 4 ? 192 : 0;
    laplace_prepare_terms(&f->expansion);
    /* Allow spare nodes for sparse and clustered trees. Full leaves remain
     * valid even if subdivision exhausts this pool. */
    size_t capacity = count / 2 + 1;
    size_t bucket_capacity = order == 4 ? 32 : 16;
    if (!octree_create(&f->tree, count, bucket_capacity, 32, capacity) ||
        capacity > SIZE_MAX / sizeof(double) / (size_t)f->expansion.terms) {
        destroy(f);
        return false;
    }
    size_t bytes = capacity * f->expansion.terms * sizeof(double);
    f->moments = malloc(bytes);
    f->locals = malloc(bytes);
    /* The expansion allocation bound above also covers one double per node. */
    f->occupied_half_width = malloc(capacity * sizeof(*f->occupied_half_width));
    if (!f->moments || !f->locals || !f->occupied_half_width) {
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
static void monomials(Fmm *f, Vec3 point, const OctreeNode *cell)
{
    double offset[] = {(point.x - cell->center.x) / cell->half_width,
                       (point.y - cell->center.y) / cell->half_width,
                       (point.z - cell->center.z) / cell->half_width};
    f->polynomial[0] = 1;
    for (int a = 1; a < f->expansion.terms; ++a) {
        int axis = f->expansion.exponent[a][0] ? 0 : (f->expansion.exponent[a][1] ? 1 : 2);
        f->polynomial[a] = offset[axis] * f->polynomial[f->expansion.lower[a][axis]];
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

/* Translate either child moments upward or parent local coefficients downward.
 * Child coordinates have half the scale of their parent's coordinates. */
static void shift(Fmm *f, size_t parent, size_t child, bool upward)
{
    monomials(f, f->tree.nodes[child].center, &f->tree.nodes[parent]);
    double *buffer = upward ? f->moments : f->locals;
    double *p = buffer + parent * f->expansion.terms;
    double *q = buffer + child * f->expansion.terms;
    for (int i = 0; i < f->expansion.splits; ++i) {
        CartesianSplit s = f->expansion.split[i];
        double weight = s.factor * f->polynomial[s.c] / (1 << f->expansion.degree[s.b]);
        if (upward) {
            p[s.a] += weight * q[s.b];
        } else {
            q[s.b] += weight * p[s.a];
        }
    }
}

/* Radius in the infinity norm: half-width of an axis-aligned bounding cube. */
static double max_axis_distance(Vec3 a, Vec3 b)
{
    return fmax(fabs(a.x - b.x), fmax(fabs(a.y - b.y), fabs(a.z - b.z)));
}

static bool upward_pass(Fmm *f, const SimulationView *state)
{
    for (size_t n = f->tree.node_count; n-- > 0;) {
        const OctreeNode *cell = &f->tree.nodes[n];
        double *moments = f->moments + n * f->expansion.terms;
        double extent = 0.0;
        if (cell->is_leaf) {
            for (size_t p = cell->begin; p < cell->end; ++p) {
                size_t i = f->tree.indices[p];
                Vec3 point = position(state, i);
                extent = fmax(extent, max_axis_distance(point, cell->center));
                monomials(f, point, cell);
                /* Gravity has one source component: particle mass. */
                for (int a = 0; a < f->expansion.terms; ++a) {
                    moments[a] += state->w[i] * f->polynomial[a];
                }
            }
        } else {
            for (int c = 0; c < 8; ++c) {
                if (cell->children[c] != SIZE_MAX) {
                    const OctreeNode *child = &f->tree.nodes[cell->children[c]];
                    double enclosing = max_axis_distance(child->center, cell->center)
                                     + f->occupied_half_width[cell->children[c]];
                    extent = fmax(extent, enclosing);
                    shift(f, n, cell->children[c], true);
                }
            }
        }
        f->occupied_half_width[n] = fmin(cell->half_width, nextafter(extent, INFINITY));
        if (!finite_terms(moments, f->expansion.terms)) {
            return false;
        }
    }
    return true;
}

typedef struct {
    const SimulationView *state;
    double *acceleration[3];
} Field;

/* One traversal owns its scratch and counts, and writes only its targets. */
typedef struct {
    double kernel[CARTESIAN_MAX_TERMS];
    size_t translations, direct_evaluations;
} FmmInteraction;

static void evaluate_direct_pair(const Fmm *f, const Field *field, FmmInteraction *work,
                                 const OctreeNode *s, const OctreeNode *t)
{
    for (size_t p = t->begin; p < t->end; ++p) {
        size_t i = f->tree.indices[p];
        Vec3 sum = {field->acceleration[0][i], field->acceleration[1][i],
                    field->acceleration[2][i]};
        sum = gravity_accumulate_sources(field->state, i, f->tree.indices,
                                         s->begin, s->end, sum);
        field->acceleration[0][i] = sum.x;
        field->acceleration[1][i] = sum.y;
        field->acceleration[2][i] = sum.z;
        work->direct_evaluations += (s->end - s->begin) - (s->begin <= p && p < s->end);
    }
}

static bool interact(const Fmm *f, const Field *field, FmmInteraction *work,
                     size_t source, size_t target)
{
    const OctreeNode *s = &f->tree.nodes[source], *t = &f->tree.nodes[target];
    /* Every occupied node is nonempty. Division keeps the pair count test
     * valid even for ranges whose product would overflow size_t. */
    size_t source_count = s->end - s->begin, target_count = t->end - t->begin;
    if (f->direct_pair_cutoff && source_count <= f->direct_pair_cutoff / target_count) {
        evaluate_direct_pair(f, field, work, s, t);
        return true;
    }
    double dx = t->center.x - s->center.x, dy = t->center.y - s->center.y;
    double dz = t->center.z - s->center.z;
    double r = hypot(hypot(dx, dy), dz);
    if (!isfinite(r)) {
        return false;
    }
    if (r > 0 && (f->occupied_half_width[source] + f->occupied_half_width[target])
                <= r * (0.5 / sqrt(3.0))) {
        double unit[] = {dx / r, dy / r, dz / r};
        laplace_translate(&f->expansion, f->moments + source * f->expansion.terms,
                          s->half_width, t->half_width, unit, r, work->kernel,
                          f->locals + target * f->expansion.terms);
        ++work->translations;
    } else if (s->is_leaf && t->is_leaf) {
        evaluate_direct_pair(f, field, work, s, t);
    } else {
        bool split_target = s->is_leaf || (!t->is_leaf && t->half_width >= s->half_width);
        const OctreeNode *split = split_target ? t : s;
        for (int c = 0; c < 8; ++c) {
            size_t child = split->children[c];
            if (child != SIZE_MAX &&
                !interact(f, field, work, split_target ? source : child, split_target ? child : target)) {
                return false;
            }
        }
    }
    return true;
}

static bool interaction_pass(Fmm *f, const Field *field)
{
    size_t targets[8] = {0};
    int jobs = 1;
#ifdef _OPENMP
    int workers = omp_get_max_threads();
    const OctreeNode *root = &f->tree.nodes[0];
    if (workers > 1 && field->state->count >= 1024 && !root->is_leaf) {
        jobs = 0;
        for (int c = 0; c < 8; ++c) {
            if (root->children[c] != SIZE_MAX) {
                targets[jobs++] = root->children[c];
            }
        }
        if (jobs < 2) {
            jobs = 1;
            targets[0] = 0;
        }
    }
    if (workers > jobs) {
        workers = jobs;
    }
#endif
    FmmInteraction work[8];
    bool valid[8];
    for (int i = 0; i < jobs; ++i) {
        work[i].translations = work[i].direct_evaluations = 0;
    }
#ifdef _OPENMP
    if (jobs > 1) {
        /* Root-to-root traversal first splits the target. Its child subtrees
         * have disjoint outputs; each retains the original source order. */
#pragma omp parallel for schedule(dynamic, 1) num_threads(workers)
        for (int i = 0; i < jobs; ++i) {
            valid[i] = interact(f, field, &work[i], 0, targets[i]);
        }
    } else
#endif
    {
        valid[0] = interact(f, field, &work[0], 0, targets[0]);
    }
    bool success = true;
    for (int i = 0; i < jobs; ++i) {
        f->translations += work[i].translations;
        f->direct_evaluations += work[i].direct_evaluations;
        success = valid[i] && success;
    }
    return success;
}

static bool downward_pass(Fmm *f, Field *field)
{
    for (size_t n = 0; n < f->tree.node_count; ++n) {
        const OctreeNode *cell = &f->tree.nodes[n];
        const double *local = f->locals + n * f->expansion.terms;
        if (!finite_terms(local, f->expansion.terms)) {
            return false;
        }
        if (!cell->is_leaf) {
            for (int c = 0; c < 8; ++c) {
                if (cell->children[c] != SIZE_MAX) {
                    shift(f, n, cell->children[c], false);
                }
            }
            continue;
        }
        for (size_t p = cell->begin; p < cell->end; ++p) {
            size_t i = f->tree.indices[p];
            monomials(f, position(field->state, i), cell);
            for (int axis = 0; axis < 3; ++axis) {
                double g = 0;
                for (int a = 1; a < f->expansion.terms; ++a) {
                    int b = f->expansion.lower[a][axis];
                    if (b >= 0) {
                        g += local[a] * f->expansion.exponent[a][axis] * f->polynomial[b];
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
    ++f->calls;
    if (!state->count) {
        return;
    }
    Field field = {state, {ax, ay, az}};
    if (octree_build(state, &f->tree)) {
        f->capacity_limited_leaves += f->tree.capacity_limited_leaves;
        size_t bytes = f->tree.node_count * f->expansion.terms * sizeof(double);
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
    f->direct_evaluations += state->count * (state->count - 1);
    n2_accelerations(state, ax, ay, az);
}

static void settings(const void *solver_context, MetadataVisitor visit, void *visitor_context)
{
    const Fmm *f = solver_context;
    visit("order", (MetadataValue){METADATA_SIZE, {.size = f->expansion.order}}, visitor_context);
    visit("direct_pair_cutoff", (MetadataValue){METADATA_SIZE, {.size = f->direct_pair_cutoff}}, visitor_context);
    visit("bucket_capacity", (MetadataValue){METADATA_SIZE, {.size = f->tree.bucket_capacity}},
          visitor_context);
    visit("max_depth", (MetadataValue){METADATA_SIZE, {.size = f->tree.max_depth}}, visitor_context);
    visit("node_capacity", (MetadataValue){METADATA_SIZE, {.size = f->tree.node_capacity}}, visitor_context);
}

static void diagnostics(const void *solver_context, MetadataVisitor visit, void *visitor_context)
{
    const Fmm *f = solver_context;
#define REPORT(field) visit(#field, (MetadataValue){METADATA_SIZE, {.size = f->field}}, visitor_context)
    REPORT(calls);
    REPORT(translations);
    REPORT(direct_evaluations);
    REPORT(capacity_limited_leaves);
    REPORT(direct_fallbacks);
#undef REPORT
}

const SolverDescriptor fmm_solver = {.name = "fmm",
                                     .acceleration = evaluate,
                                     .create = create,
                                     .destroy = destroy,
                                     .settings = settings,
                                     .diagnostics = diagnostics};
