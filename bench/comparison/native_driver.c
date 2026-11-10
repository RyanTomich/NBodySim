#define _POSIX_C_SOURCE 200809L

/* Include Rushil's unchanged implementation to call its actual static step.
 * The two src builds instead link their own production simulation entry point. */
#ifdef RUSHIL_SOURCE
#define main rushil_original_main
#include RUSHIL_SOURCE
#undef main
#else
#include "solver.h"
#endif
#include "io.h"
#include "parse.h"

#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void copy_state(Simulation *destination, const Simulation *source)
{
    MutableSimulationView out = simulation_mutable_view(destination);
    SimulationView in = simulation_view(source);
    double *outputs[] = {out.w, out.x, out.y, out.z, out.vx, out.vy, out.vz,
                         out.ax, out.ay, out.az};
    const double *inputs[] = {in.w, in.x, in.y, in.z, in.vx, in.vy, in.vz,
                             in.ax, in.ay, in.az};
    for (size_t field = 0; field < 10; ++field) {
        memcpy(outputs[field], inputs[field], in.count * sizeof(double));
    }
}

static bool valid_state(const Simulation *state, const Simulation *initial,
                        const Simulation *previous)
{
    SimulationView v = simulation_view(state), in = simulation_view(initial);
    SimulationView before = simulation_view(previous);
    const double *fields[] = {v.w, v.x, v.y, v.z, v.vx, v.vy, v.vz, v.ax, v.ay, v.az};
    const double *expected[] = {before.w, before.x, before.y, before.z,
                               before.vx, before.vy, before.vz, before.ax, before.ay, before.az};
    if (v.count != in.count || v.count != before.count) {
        return false;
    }
    for (size_t i = 0; i < v.count; ++i) {
        if (v.w[i] != in.w[i]) {
            return false;
        }
        for (size_t field = 0; field < 10; ++field) {
            if (!isfinite(fields[field][i]) || fields[field][i] != expected[field][i]) {
                return false;
            }
        }
    }
    return true;
}

/* Independent ordered-pair direct reference at final positions, outside timing.
 * Accumulate in long double and avoid cubing the distance. */
static bool report_force_error(const Simulation *state)
{
    SimulationView v = simulation_view(state);
    size_t targets = v.count < 128 ? v.count : 128;
    long double error = 0, norm = 0;
    for (size_t k = 0; k < targets; ++k) {
        size_t i = targets > 1 ? k * (v.count - 1) / (targets - 1) : 0;
        long double sum[3] = {0};
        for (size_t j = 0; j < v.count; ++j) {
            if (i == j) {
                continue;
            }
            long double d[] = {(long double)v.x[j] - v.x[i],
                               (long double)v.y[j] - v.y[i],
                               (long double)v.z[j] - v.z[i]};
            long double distance = hypotl(hypotl(d[0], d[1]), d[2]);
            if (distance == 0) {
                continue;
            }
            long double scale = v.w[j] / distance / distance;
            for (size_t axis = 0; axis < 3; ++axis) {
                sum[axis] += scale * (d[axis] / distance);
            }
        }
        error = hypotl(error, hypotl(hypotl(v.ax[i] - sum[0], v.ay[i] - sum[1]), v.az[i] - sum[2]));
        norm = hypotl(norm, hypotl(hypotl(sum[0], sum[1]), sum[2]));
    }
    long double relative = norm > 0 ? error / norm : (error == 0 ? 0 : INFINITY);
    long double rms = error / sqrtl((long double)targets);
    if (!isfinite(relative) || !isfinite(rms)) {
        return false;
    }
    fprintf(stderr, "force_error: {\"targets\": %zu, \"relative_l2\": %.17Lg, \"rms_absolute\": %.17Lg}\n",
            targets, relative, rms);
    return true;
}

int main(int argc, char **argv)
{
    double dt, opening;
    size_t steps, samples, warmups;
    if (argc != 8 || !parse_double(argv[2], &dt) || dt <= 0
        || !parse_size(argv[3], &steps) || steps == 0
        || !parse_size(argv[4], &samples) || samples == 0
        || !parse_size(argv[5], &warmups) || warmups > SIZE_MAX - samples
        || !parse_double(argv[6], &opening) || opening <= 0) {
        fprintf(stderr, "Usage: native_driver INPUT DT STEPS SAMPLES WARMUPS THETA FINAL\n");
        return 1;
    }
    int result = 1;
    Simulation *initial = NULL, *state = NULL, *previous = NULL;
    FILE *input = fopen(argv[1], "r");
    if (!input) {
        perror("input");
        return 1;
    }
    bool loaded = read_bodies(input, &initial);
    int closed = fclose(input);
    if (!loaded || closed != 0) {
        fprintf(stderr, "Invalid input.\n");
        goto cleanup;
    }
    size_t count = simulation_view(initial).count;
    if (count < 2 || count > (INT_MAX - 1) / 16
        || !simulation_create(&state, count) || !simulation_create(&previous, count)) {
        fprintf(stderr, "Invalid count or allocation failure.\n");
        goto cleanup;
    }
    omp_set_dynamic(0);
    int team = 1;
#ifndef ORIGINAL_SRC
#pragma omp parallel
    {
#pragma omp single
        team = omp_get_num_threads();
    }
#endif
    fprintf(stderr, "threads: {\"requested\": %d, \"effective\": %d}\n", omp_get_max_threads(), team);
#ifdef RUSHIL_SOURCE
    fprintf(stderr, "rushil_task_cutoff: %d\n", SUMMARIZE_BASE_CASE);
#endif
    puts("sample,seconds");
    for (size_t run = 0; run < warmups + samples; ++run) {
        copy_state(state, initial);
#ifdef RUSHIL_SOURCE
        simulation = state;
        bodies = simulation_mutable_view(state);
        n = (int)count;
        theta = opening;
        position[0] = bodies.x; position[1] = bodies.y; position[2] = bodies.z;
        velocity[0] = bodies.vx; velocity[1] = bodies.vy; velocity[2] = bodies.vz;
        acceleration[0] = bodies.ax; acceleration[1] = bodies.ay; acceleration[2] = bodies.az;
        capacity = 16 * n + 1;
        nodes = 0;
        tree = calloc((size_t)capacity, sizeof(*tree));
        if (!tree) {
            goto cleanup;
        }
#else
        Solver solver;
        const char *error;
#ifdef N2_BACKEND
        bool prepared = solver_create(&solver, "n2", count, NULL, 0, &error);
#else
        SolverOption option = {"theta", argv[6]};
        bool prepared = solver_create(&solver, "barnes-hut", count, &option, 1, &error);
#endif
        if (!prepared) {
            fprintf(stderr, "%s\n", error);
            goto cleanup;
        }
#endif
        struct timespec start, end;
        bool timed = clock_gettime(CLOCK_MONOTONIC, &start) == 0;
        if (timed) {
#ifdef RUSHIL_SOURCE
            for (size_t s = 0; s < steps; ++s) {
                (void)step(dt);
            }
#else
            simulate_with_solver(state, dt, steps, solver.acceleration, solver.context);
#endif
            timed = clock_gettime(CLOCK_MONOTONIC, &end) == 0;
        }
#ifdef RUSHIL_SOURCE
        free(tree);
        tree = NULL;
        nodes = 0;
#else
        solver_destroy(&solver);
#endif
        if (!timed) {
            perror("clock_gettime");
            goto cleanup;
        }
        if (run == 0) {
            copy_state(previous, state);
        }
        if (!valid_state(state, initial, previous)) {
            fprintf(stderr, "Invalid or nondeterministic final state at run %zu.\n", run);
            goto cleanup;
        }
        if (run >= warmups) {
            double seconds = (double)(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
            /* Tiny validation fixtures may finish within one clock tick. The
             * runner requires strictly positive times for benchmark samples. */
            if (!isfinite(seconds) || seconds < 0) {
                fprintf(stderr, "Invalid elapsed time.\n");
                goto cleanup;
            }
            printf("%zu,%.17g\n", run - warmups + 1, seconds);
            fflush(stdout);
        }
    }
    if (!report_force_error(state)) {
        fprintf(stderr, "Invalid sampled force error.\n");
        goto cleanup;
    }
    FILE *output = fopen(argv[7], "w");
    if (!output) {
        perror("final state");
        goto cleanup;
    }
    bool written = write_bodies(output, state);
    closed = fclose(output);
    if (!written || closed != 0 || fflush(stdout) != 0 || ferror(stdout)) {
        goto cleanup;
    }
    result = 0;
cleanup:
    simulation_destroy(initial);
    simulation_destroy(state);
    simulation_destroy(previous);
    return result;
}
