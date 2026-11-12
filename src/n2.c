#include "n2_internal.h"
#include "solver_internal.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    size_t count;
    int reserved_workers;
    double *scratch; /* Three count-element arrays per reserved worker. */
} N2Workspace;

static inline void evaluate_row(const SimulationView *state, size_t i,
                                double *restrict ax, double *restrict ay,
                                double *restrict az, bool initialize)
{
    const double *restrict x = state->x;
    const double *restrict y = state->y;
    const double *restrict z = state->z;
    const double *restrict mass = state->w;
    const double x_i = x[i];
    const double y_i = y[i];
    const double z_i = z[i];
    const double mass_i = mass[i];
    double sum_x = initialize ? 0.0 : ax[i];
    double sum_y = initialize ? 0.0 : ay[i];
    double sum_z = initialize ? 0.0 : az[i];
    for (size_t j = i + 1; j < state->count; ++j) {
        Vec3 displacement = {
            x[j] - x_i,
            y[j] - y_i,
            z[j] - z_i
        };
        double r2 = displacement.x * displacement.x
                  + displacement.y * displacement.y
                  + displacement.z * displacement.z;
        if (r2 == 0.0) {
            if (initialize) {
                ax[j] = 0.0;
                ay[j] = 0.0;
                az[j] = 0.0;
            }
            continue;
        }
        double inv_r3 = 1.0 / (r2 * sqrt(r2));
        double scale_i = mass_i * inv_r3;
        double scale_j = mass[j] * inv_r3;
        sum_x += scale_j * displacement.x;
        sum_y += scale_j * displacement.y;
        sum_z += scale_j * displacement.z;
        /* Initialize by subtracting from positive zero, without reading outputs. */
        ax[j] = (initialize ? 0.0 : ax[j]) - scale_i * displacement.x;
        ay[j] = (initialize ? 0.0 : ay[j]) - scale_i * displacement.y;
        az[j] = (initialize ? 0.0 : az[j]) - scale_i * displacement.z;
    }
    ax[i] = sum_x;
    ay[i] = sum_y;
    az[i] = sum_z;
}

void n2_accelerations(const SimulationView *state, double *restrict ax,
                     double *restrict ay, double *restrict az)
{
    const size_t count = state->count;
    if (count == 0) {
        return;
    }
    /* The first target initializes every output with its first contribution,
     * without reading the old workspace or making a separate clearing pass. */
    evaluate_row(state, 0, ax, ay, az, true);

    for (size_t i = 1; i < count; ++i) {
        evaluate_row(state, i, ax, ay, az, false);
    }
}

/* The serial kernel above is also used directly by FMM fallbacks. */
#ifdef _OPENMP
static void evaluate_parallel(const SimulationView *state, double *ax, double *ay,
                              double *az, N2Workspace *workspace)
{
    const size_t count = state->count;
    const size_t worker_stride = 3 * count;
#pragma omp parallel num_threads(workspace->reserved_workers) default(none) shared(state, ax, ay, az, workspace, count, worker_stride)
    {
        const int team_size = omp_get_num_threads();
        if (team_size == 1) {
            n2_accelerations(state, ax, ay, az);
        } else {
            double *local_ax = workspace->scratch + (size_t)omp_get_thread_num() * worker_stride;
            double *local_ay = local_ax + count;
            double *local_az = local_ay + count;
            memset(local_ax, 0, worker_stride * sizeof(*local_ax));
            /* Cyclic rows spread the long early rows across workers. Each pair
             * updates only this worker's buffers, including both endpoints. */
#pragma omp for schedule(static, 1)
            for (size_t i = 0; i < count - 1; ++i) {
                evaluate_row(state, i, local_ax, local_ay, local_az, false);
            }
            /* The pair-loop barrier completes every private array before any
             * output is combined. Ignore buffers of absent runtime workers. */
#pragma omp for schedule(static)
            for (size_t i = 0; i < count; ++i) {
                double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
                for (int worker = 0; worker < team_size; ++worker) {
                    const double *source = workspace->scratch + (size_t)worker * worker_stride;
                    sum_x += source[i];
                    sum_y += source[count + i];
                    sum_z += source[2 * count + i];
                }
                ax[i] = sum_x;
                ay[i] = sum_y;
                az[i] = sum_z;
            }
        }
    }
}
#endif

static void evaluate(const SimulationView *state, double *ax, double *ay,
                     double *az, void *solver_context)
{
    N2Workspace *workspace = solver_context;
    assert(workspace && workspace->count == state->count);
#ifdef _OPENMP
    if (workspace->scratch) {
        evaluate_parallel(state, ax, ay, az, workspace);
        return;
    }
#endif
    n2_accelerations(state, ax, ay, az);
}

static void destroy(void *solver_context)
{
    N2Workspace *workspace = solver_context;
    if (workspace) {
        free(workspace->scratch);
        free(workspace);
    }
}

static bool create(void **solver_context, size_t count, const SolverOption *options,
                   size_t option_count, const char **error)
{
    (void)options;
    if (option_count) {
        *error = "Solver does not accept options.";
        return false;
    }
    N2Workspace *workspace = calloc(1, sizeof(*workspace));
    if (!workspace) {
        goto allocation_failure;
    }
    workspace->count = count;
    workspace->reserved_workers = 1;
#ifdef _OPENMP
    workspace->reserved_workers = omp_get_max_threads();
#endif
    if (workspace->reserved_workers <= 0) {
        goto allocation_failure;
    }
    if (workspace->reserved_workers > 1 && count > 1) {
        if (count > SIZE_MAX / (3 * sizeof(double)) / (size_t)workspace->reserved_workers) {
            goto allocation_failure;
        }
        workspace->scratch = malloc((size_t)workspace->reserved_workers * 3 * count * sizeof(double));
        if (!workspace->scratch) {
            goto allocation_failure;
        }
    }
    *solver_context = workspace;
    return true;

allocation_failure:
    destroy(workspace);
    *error = "Could not allocate n2 workspace.";
    return false;
}

const SolverDescriptor n2_solver = {
    .name = "n2",
    .acceleration = evaluate,
    .create = create,
    .destroy = destroy
};
