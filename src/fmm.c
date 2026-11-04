#include "fmm_order4_internal.h"
#include "gravity_internal.h"
#include "octree_internal.h"
#include "parse.h"
#include "solver_internal.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { MAX_ORDER = 6, MAX_TERMS = 84, MAX_SPLITS = 924 };

/* Each split represents a=b+c with factor choose(a,b), componentwise.
 * 84 = choose(6+3,3) terms and 924 = choose(6+6,6) splits. */
typedef struct {
    int a, b, c, factor;
} Split;

typedef struct {
    Octree tree;
    size_t order, direct_pair_cutoff;
    int terms, splits;
    int degree[MAX_TERMS], exponent[MAX_TERMS][3], lower[MAX_TERMS][3];
    Split split[MAX_SPLITS];
    double *moments, *locals;
    double polynomial[MAX_TERMS], kernel[MAX_TERMS];
    size_t calls, translations, direct_evaluations, capacity_limited_leaves, direct_fallbacks;
} Fmm;

static int choose(int n, int k)
{
    int result = 1;
    for (int i = 1; i <= k; ++i) {
        result = result * (n - i + 1) / i;
    }
    return result;
}

static void prepare_terms(Fmm *f)
{
    int index[MAX_ORDER + 1][MAX_ORDER + 1][MAX_ORDER + 1] = {0};
    for (int d = 0; d <= (int)f->order; ++d) {
        for (int x = 0; x <= d; ++x) {
            for (int y = 0; y <= d - x; ++y) {
                int z = d - x - y, a = f->terms++;
                index[x][y][z] = a;
                f->degree[a] = d;
                f->exponent[a][0] = x;
                f->exponent[a][1] = y;
                f->exponent[a][2] = z;
                f->lower[a][0] = x ? index[x - 1][y][z] : -1;
                f->lower[a][1] = y ? index[x][y - 1][z] : -1;
                f->lower[a][2] = z ? index[x][y][z - 1] : -1;
                for (int i = 0; i <= x; ++i) {
                    for (int j = 0; j <= y; ++j) {
                        for (int k = 0; k <= z; ++k) {
                            f->split[f->splits++] =
                                (Split){a, index[i][j][k], index[x - i][y - j][z - k],
                                        choose(x, i) * choose(y, j) * choose(z, k)};
                        }
                    }
                }
            }
        }
    }
}

static void destroy(void *context)
{
    Fmm *f = context;
    octree_destroy(&f->tree);
    free(f->moments);
    free(f->locals);
    free(f);
}

static bool create(void **context, size_t count, const SolverOption *options, size_t option_count,
                   const char **error)
{
    size_t order = 4;
    for (size_t i = 0; i < option_count; ++i) {
        if (strcmp(options[i].name, "order") || !parse_size(options[i].value, &order) ||
            order < 1 || order > MAX_ORDER) {
            *error = "FMM accepts only order, an integer from 1 through 6.";
            return false;
        }
    }
    *error = "Could not allocate FMM workspace.";
    Fmm *f = calloc(1, sizeof(*f));
    if (!f) {
        return false;
    }
    f->order = order;
    f->direct_pair_cutoff = order == 4 ? 192 : 0;
    prepare_terms(f);
    /* Allow spare nodes for sparse and clustered trees. Full leaves remain
     * valid even if subdivision exhausts this pool. */
    size_t capacity = count / 2 + 1;
    size_t bucket_capacity = order == 4 ? 32 : 16;
    if (!octree_create(&f->tree, count, bucket_capacity, 32, capacity) ||
        capacity > SIZE_MAX / sizeof(double) / (size_t)f->terms) {
        destroy(f);
        return false;
    }
    size_t bytes = capacity * f->terms * sizeof(double);
    f->moments = malloc(bytes);
    f->locals = malloc(bytes);
    if (!f->moments || !f->locals) {
        destroy(f);
        return false;
    }
    *context = f;
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
    for (int a = 1; a < f->terms; ++a) {
        int axis = f->exponent[a][0] ? 0 : (f->exponent[a][1] ? 1 : 2);
        f->polynomial[a] = offset[axis] * f->polynomial[f->lower[a][axis]];
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
    double *p = buffer + parent * f->terms;
    double *q = buffer + child * f->terms;
    for (int i = 0; i < f->splits; ++i) {
        Split s = f->split[i];
        double weight = s.factor * f->polynomial[s.c] / (1 << f->degree[s.b]);
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
        OctreeNode *cell = &f->tree.nodes[n];
        double *moments = f->moments + n * f->terms;
        double extent = 0.0;
        if (cell->is_leaf) {
            for (size_t p = cell->begin; p < cell->end; ++p) {
                size_t i = f->tree.indices[p];
                Vec3 point = position(state, i);
                extent = fmax(extent, max_axis_distance(point, cell->center));
                monomials(f, point, cell);
                for (int a = 0; a < f->terms; ++a) {
                    moments[a] += state->w[i] * f->polynomial[a];
                }
            }
        } else {
            for (int c = 0; c < 8; ++c) {
                if (cell->children[c] != SIZE_MAX) {
                    const OctreeNode *child = &f->tree.nodes[cell->children[c]];
                    double enclosing = max_axis_distance(child->center, cell->center)
                                     + child->occupied_half_width;
                    extent = fmax(extent, enclosing);
                    shift(f, n, cell->children[c], true);
                }
            }
        }
        cell->occupied_half_width = fmin(cell->half_width, nextafter(extent, INFINITY));
        if (!finite_terms(moments, f->terms)) {
            return false;
        }
    }
    return true;
}

static void translation_kernel(Fmm *f, const double unit[3])
{
    if (f->order == 4) {
        fmm_kernel_order4(f->kernel, unit);
        return;
    }
    /* Taylor coefficients of 1/|unit+x| from the homogeneous Legendre
     * recurrence: d H_d = -(2d-1)(unit.x)H_(d-1) -(d-1)|x|²H_(d-2). */
    f->kernel[0] = 1;
    for (int a = 1; a < f->terms; ++a) {
        int d = f->degree[a];
        double linear = 0, quadratic = 0;
        for (int axis = 0; axis < 3; ++axis) {
            int b = f->lower[a][axis];
            if (b >= 0) {
                linear += unit[axis] * f->kernel[b];
                int c = f->lower[b][axis];
                if (c >= 0) {
                    quadratic += f->kernel[c];
                }
            }
        }
        f->kernel[a] = -((2 * d - 1) * linear + (d - 1) * quadratic) / d;
    }
}

static void translate(Fmm *f, size_t source, size_t target, const double unit[3], double r)
{
    translation_kernel(f, unit);
    double source_scale[MAX_ORDER + 1] = {1}, target_scale[MAX_ORDER + 1] = {1};
    for (size_t d = 1; d <= f->order; ++d) {
        source_scale[d] = source_scale[d - 1] * (-f->tree.nodes[source].half_width / r);
        target_scale[d] = target_scale[d - 1] * (f->tree.nodes[target].half_width / r);
    }
    const double *moments = f->moments + source * f->terms;
    double *local = f->locals + target * f->terms;
    if (f->order == 4) {
        fmm_contract_order4(local, moments, f->kernel, source_scale, target_scale, r);
    } else {
        for (int i = 0; i < f->splits; ++i) {
            Split s = f->split[i];
            local[s.b] += (moments[s.c] / r) * source_scale[f->degree[s.c]] *
                          target_scale[f->degree[s.b]] * s.factor * f->kernel[s.a];
        }
    }
    ++f->translations;
}

typedef struct {
    const SimulationView *state;
    double *acceleration[3];
} Field;

static void evaluate_direct_pair(Fmm *f, const Field *field,
                                 const OctreeNode *s, const OctreeNode *t)
{
    for (size_t p = t->begin; p < t->end; ++p) {
        size_t i = f->tree.indices[p];
        Vec3 sum = {field->acceleration[0][i], field->acceleration[1][i],
                    field->acceleration[2][i]};
        for (size_t q = s->begin; q < s->end; ++q) {
            size_t j = f->tree.indices[q];
            if (i == j) {
                continue;
            }
            Vec3 g = gravity(field->state, i, j);
            sum.x += g.x;
            sum.y += g.y;
            sum.z += g.z;
        }
        field->acceleration[0][i] = sum.x;
        field->acceleration[1][i] = sum.y;
        field->acceleration[2][i] = sum.z;
        f->direct_evaluations += (s->end - s->begin) - (s->begin <= p && p < s->end);
    }
}

static bool interact(Fmm *f, Field *field, size_t source, size_t target)
{
    const OctreeNode *s = &f->tree.nodes[source], *t = &f->tree.nodes[target];
    /* Every occupied node is nonempty. Division keeps the pair count test
     * valid even for ranges whose product would overflow size_t. */
    size_t source_count = s->end - s->begin, target_count = t->end - t->begin;
    if (f->direct_pair_cutoff && source_count <= f->direct_pair_cutoff / target_count) {
        evaluate_direct_pair(f, field, s, t);
        return true;
    }
    double dx = t->center.x - s->center.x, dy = t->center.y - s->center.y;
    double dz = t->center.z - s->center.z;
    double r = hypot(hypot(dx, dy), dz);
    if (!isfinite(r)) {
        return false;
    }
    if (r > 0 && (s->occupied_half_width + t->occupied_half_width) <= r * (0.5 / sqrt(3.0))) {
        double unit[] = {dx / r, dy / r, dz / r};
        translate(f, source, target, unit, r);
    } else if (s->is_leaf && t->is_leaf) {
        evaluate_direct_pair(f, field, s, t);
    } else {
        bool split_target = s->is_leaf || (!t->is_leaf && t->half_width >= s->half_width);
        const OctreeNode *split = split_target ? t : s;
        for (int c = 0; c < 8; ++c) {
            size_t child = split->children[c];
            if (child != SIZE_MAX &&
                !interact(f, field, split_target ? source : child, split_target ? child : target)) {
                return false;
            }
        }
    }
    return true;
}

static bool downward_pass(Fmm *f, Field *field)
{
    for (size_t n = 0; n < f->tree.node_count; ++n) {
        const OctreeNode *cell = &f->tree.nodes[n];
        const double *local = f->locals + n * f->terms;
        if (!finite_terms(local, f->terms)) {
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
                for (int a = 1; a < f->terms; ++a) {
                    int b = f->lower[a][axis];
                    if (b >= 0) {
                        g += local[a] * f->exponent[a][axis] * f->polynomial[b];
                    }
                }
                field->acceleration[axis][i] += g / cell->half_width;
                if (!isfinite(field->acceleration[axis][i])) {
                    return false;
                }
            }
        }
    }
    return true;
}

static void evaluate(const SimulationView *state, double *ax, double *ay, double *az, void *context)
{
    Fmm *f = context;
    ++f->calls;
    if (!state->count) {
        return;
    }
    Field field = {state, {ax, ay, az}};
    if (octree_build(state, &f->tree)) {
        f->capacity_limited_leaves += f->tree.capacity_limited_leaves;
        size_t bytes = f->tree.node_count * f->terms * sizeof(double);
        memset(f->moments, 0, bytes);
        memset(f->locals, 0, bytes);
        for (int axis = 0; axis < 3; ++axis) {
            memset(field.acceleration[axis], 0, state->count * sizeof(double));
        }
        if (upward_pass(f, state) && interact(f, &field, 0, 0) && downward_pass(f, &field)) {
            return;
        }
    }
    ++f->direct_fallbacks;
    f->direct_evaluations += state->count * (state->count - 1);
    gravity_direct(state, ax, ay, az, NULL);
}

static void settings(const void *context, MetadataVisitor visit, void *user)
{
    const Fmm *f = context;
    visit("order", (MetadataValue){METADATA_SIZE, {.size = f->order}}, user);
    visit("direct_pair_cutoff", (MetadataValue){METADATA_SIZE, {.size = f->direct_pair_cutoff}}, user);
    visit("bucket_capacity", (MetadataValue){METADATA_SIZE, {.size = f->tree.bucket_capacity}},
          user);
    visit("max_depth", (MetadataValue){METADATA_SIZE, {.size = f->tree.max_depth}}, user);
    visit("node_capacity", (MetadataValue){METADATA_SIZE, {.size = f->tree.node_capacity}}, user);
}

static void diagnostics(const void *context, MetadataVisitor visit, void *user)
{
    const Fmm *f = context;
#define REPORT(field) visit(#field, (MetadataValue){METADATA_SIZE, {.size = f->field}}, user)
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
