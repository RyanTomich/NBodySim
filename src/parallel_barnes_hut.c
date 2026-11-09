/* Rushil's Barnes-Hut algorithm as a prepared solver: reusable topology,
 * bottom-up task summaries, and independent parallel target walks.
 * All storage is private to a solver instance; evaluation allocates nothing. */
#include "solver_internal.h"
#include "n2_internal.h"
#include "parse.h"
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double lo[3], hi[3], mass, com[3];
    int child[8], body, count;
} ParallelNode;

typedef struct {
    ParallelNode *tree;
    size_t count, capacity, max_depth, task_cutoff;
    int nodes, actual_threads;
    double theta;
    bool failed;
    size_t calls, rebuilds, reuses, direct_fallbacks;
} ParallelBarnesHut;

static double coordinate(const SimulationView *s, int k, size_t i)
{
    const double *position[3] = {s->x, s->y, s->z};
    return position[k][i];
}

static int new_node(ParallelBarnesHut *w, const double lo[3], const double hi[3])
{
    if ((size_t)w->nodes == w->capacity) { w->failed = true; return -1; }
    int index = w->nodes++;
    ParallelNode *q = &w->tree[index];
    *q = (ParallelNode){.body = -1};
    for (int k = 0; k < 3; ++k) { q->lo[k] = lo[k]; q->hi[k] = hi[k]; }
    for (int c = 0; c < 8; ++c) q->child[c] = -1;
    return index;
}

static int octant(const ParallelNode *q, const SimulationView *s, int body)
{
    int o = 0;
    for (int k = 0; k < 3; ++k)
        o |= (coordinate(s, k, (size_t)body) >= (q->lo[k] + q->hi[k]) / 2) << k;
    return o;
}

static int child_node(ParallelBarnesHut *w, int parent, int o)
{
    ParallelNode *q = &w->tree[parent];
    if (q->child[o] >= 0) return q->child[o];
    double lo[3], hi[3];
    for (int k = 0; k < 3; ++k) {
        double mid = (q->lo[k] + q->hi[k]) / 2;
        lo[k] = (o & (1 << k)) ? mid : q->lo[k];
        hi[k] = (o & (1 << k)) ? q->hi[k] : mid;
    }
    q->child[o] = new_node(w, lo, hi);
    return q->child[o];
}

static void insert(ParallelBarnesHut *w, const SimulationView *s, int ni, int bi, size_t depth)
{
    if (w->failed) return;
    ParallelNode *q = &w->tree[ni];
    q->count++;
    bool has_child = false;
    for (int c = 0; c < 8; ++c) has_child |= q->child[c] >= 0;
    if (q->body < 0 && !has_child) { q->body = bi; return; }
    if (depth >= w->max_depth) { w->failed = true; return; }
    if (q->body >= 0) {
        int old = q->body;
        q->body = -1;
        int child = child_node(w, ni, octant(q, s, old));
        if (w->failed) return;
        insert(w, s, child, old, depth + 1);
    }
    if (w->failed) return;
    int child = child_node(w, ni, octant(q, s, bi));
    if (!w->failed) insert(w, s, child, bi, depth + 1);
}

static bool contains(const ParallelNode *q, const SimulationView *s, int bi)
{
    for (int k = 0; k < 3; ++k) {
        double x = coordinate(s, k, (size_t)bi);
        if (!(x >= q->lo[k] && x <= q->hi[k])) return false;
    }
    return true;
}

static bool tree_is_valid(const ParallelBarnesHut *w, const SimulationView *s)
{
    if (!w->nodes || w->failed) return false;
    for (int i = 0; i < w->nodes; ++i)
        if (w->tree[i].body >= 0 && !contains(&w->tree[i], s, w->tree[i].body)) return false;
    return true;
}

static void rebuild(ParallelBarnesHut *w, const SimulationView *s)
{
    double lo[3] = {INFINITY, INFINITY, INFINITY}, hi[3] = {-INFINITY, -INFINITY, -INFINITY};
    w->nodes = 0;
    w->failed = false;
    ++w->rebuilds;
    for (size_t i = 0; i < s->count; ++i) for (int k = 0; k < 3; ++k) {
        double x = coordinate(s, k, i);
        lo[k] = fmin(lo[k], x); hi[k] = fmax(hi[k], x);
    }
    double side = 0;
    for (int k = 0; k < 3; ++k) side = fmax(side, hi[k] - lo[k]);
    side = side ? side * 1.001 : 1;
    for (int k = 0; k < 3; ++k) {
        double mid = (lo[k] + hi[k]) / 2;
        lo[k] = mid - side / 2; hi[k] = mid + side / 2;
        if (!isfinite(lo[k]) || !isfinite(hi[k])) { w->failed = true; return; }
    }
    new_node(w, lo, hi);
    for (size_t i = 0; i < s->count && !w->failed; ++i) insert(w, s, 0, (int)i, 0);
}

static void summarize(ParallelBarnesHut *w, const SimulationView *s, int ni)
{
    ParallelNode *q = &w->tree[ni];
    q->mass = 0;
    for (int k = 0; k < 3; ++k) q->com[k] = 0;
    if (q->body >= 0) {
        q->mass = s->w[q->body];
        for (int k = 0; k < 3; ++k) q->com[k] = coordinate(s, k, (size_t)q->body);
        return;
    }
    for (int c = 0; c < 8; ++c) if (q->child[c] >= 0) {
        int child = q->child[c];
        if ((size_t)q->count > w->task_cutoff) {
            #pragma omp task firstprivate(w, s, child)
            summarize(w, s, child);
        } else {
            summarize(w, s, child);
        }
    }
    #pragma omp taskwait
    for (int c = 0; c < 8; ++c) if (q->child[c] >= 0) {
        ParallelNode *child = &w->tree[q->child[c]];
        q->mass += child->mass;
        for (int k = 0; k < 3; ++k) q->com[k] += child->mass * child->com[k];
    }
    for (int k = 0; k < 3; ++k) q->com[k] /= q->mass;
}

static void walk(const ParallelBarnesHut *w, const SimulationView *s, int bi, int ni, double a[3])
{
    const ParallelNode *q = &w->tree[ni];
    if (q->body == bi || !q->mass) return;
    double d[3], r2 = 0, width = q->hi[0] - q->lo[0];
    for (int k = 0; k < 3; ++k) {
        d[k] = q->com[k] - coordinate(s, k, (size_t)bi);
        r2 += d[k] * d[k];
    }
    if (q->body >= 0 || (!contains(q, s, bi) && width * width < w->theta * w->theta * r2)) {
        double scale = r2 == 0 ? 0 : q->mass / (r2 * sqrt(r2));
        for (int k = 0; k < 3; ++k) a[k] += scale * d[k];
        return;
    }
    for (int c = 0; c < 8; ++c) if (q->child[c] >= 0) walk(w, s, bi, q->child[c], a);
}

static void evaluate(const SimulationView *s, double *ax, double *ay, double *az, void *context)
{
    ParallelBarnesHut *w = context;
    ++w->calls;
    if (!s->count) return;
    if (!tree_is_valid(w, s)) rebuild(w, s); else ++w->reuses;
    if (!w->failed) {
        #pragma omp parallel
        {
            #pragma omp single
            {
                w->actual_threads = omp_get_num_threads();
                summarize(w, s, 0);
            }
        }
        for (int i = 0; i < w->nodes; ++i) {
            ParallelNode *q = &w->tree[i];
            if (!isfinite(q->mass)) w->failed = true;
            for (int k = 0; k < 3; ++k) if (!isfinite(q->com[k])) w->failed = true;
        }
    }
    if (w->failed) {
        ++w->direct_fallbacks;
        n2_accelerations(s, ax, ay, az);
        return;
    }
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < s->count; ++i) {
        double a[3] = {0};
        walk(w, s, (int)i, 0, a);
        ax[i] = a[0]; ay[i] = a[1]; az[i] = a[2];
    }
}

static void destroy(void *context)
{
    ParallelBarnesHut *w = context;
    free(w->tree); free(w);
}

static bool create(void **context, size_t count, const SolverOption *options, size_t option_count, const char **error)
{
    *context = NULL;
    double theta = 0.5;
    size_t cutoff = 256, depth = 64, capacity = 0;
    for (size_t i = 0; i < option_count; ++i) {
        bool valid = false;
        if (strcmp(options[i].name, "theta") == 0) valid = parse_double(options[i].value, &theta) && theta >= 0;
        if (strcmp(options[i].name, "task_cutoff") == 0) valid = parse_size(options[i].value, &cutoff);
        if (strcmp(options[i].name, "max_depth") == 0) valid = parse_size(options[i].value, &depth) && depth <= 64;
        if (strcmp(options[i].name, "node_capacity") == 0) valid = parse_size(options[i].value, &capacity);
        if (!valid) { *error = "Invalid parallel Barnes-Hut option."; return false; }
    }
    if (count > (INT_MAX - 1) / 16) { *error = "Parallel Barnes-Hut count too large."; return false; }
    if (!capacity) capacity = 16 * count + 1;
    if (capacity > INT_MAX || capacity > SIZE_MAX / sizeof(ParallelNode)) {
        *error = "Parallel Barnes-Hut capacity too large."; return false;
    }
    ParallelBarnesHut *w = calloc(1, sizeof *w);
    if (!w) { *error = "Parallel Barnes-Hut allocation failed."; return false; }
    *w = (ParallelBarnesHut){.count = count, .capacity = capacity, .theta = theta,
                            .max_depth = depth, .task_cutoff = cutoff};
    w->tree = calloc(capacity, sizeof *w->tree);
    if (!w->tree) { free(w); *error = "Parallel Barnes-Hut allocation failed."; return false; }
    *context = w;
    return true;
}

static void settings(const void *context, MetadataVisitor visit, void *visitor_context)
{
    const ParallelBarnesHut *w = context;
    visit("theta", (MetadataValue){METADATA_DOUBLE, {.real = w->theta}}, visitor_context);
    visit("task_cutoff", (MetadataValue){METADATA_SIZE, {.size = w->task_cutoff}}, visitor_context);
    visit("max_depth", (MetadataValue){METADATA_SIZE, {.size = w->max_depth}}, visitor_context);
    visit("node_capacity", (MetadataValue){METADATA_SIZE, {.size = w->capacity}}, visitor_context);
    visit("max_threads", (MetadataValue){METADATA_SIZE, {.size = (size_t)omp_get_max_threads()}}, visitor_context);
}

static void diagnostics(const void *context, MetadataVisitor visit, void *visitor_context)
{
    const ParallelBarnesHut *w = context;
    visit("actual_threads", (MetadataValue){METADATA_SIZE, {.size = (size_t)w->actual_threads}}, visitor_context);
    visit("calls", (MetadataValue){METADATA_SIZE, {.size = w->calls}}, visitor_context);
    visit("rebuilds", (MetadataValue){METADATA_SIZE, {.size = w->rebuilds}}, visitor_context);
    visit("reuses", (MetadataValue){METADATA_SIZE, {.size = w->reuses}}, visitor_context);
    visit("direct_fallbacks", (MetadataValue){METADATA_SIZE, {.size = w->direct_fallbacks}}, visitor_context);
}

const SolverDescriptor parallel_barnes_hut_solver = {
    .name = "parallel-barnes-hut", .acceleration = evaluate,
    .create = create, .destroy = destroy, .settings = settings, .diagnostics = diagnostics
};
