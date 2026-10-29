#define _POSIX_C_SOURCE 200809L

#include "gravity.h"
#include "approximate.h"
#include "io.h"
#include "perf_gate.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static bool finite_state(const Simulation *simulation)
{
    /* Acceleration arrays can move; inactive workspace is not saved state. */
    const double *fields[] = {
        simulation->w, simulation->x, simulation->y, simulation->z,
        simulation->vx, simulation->vy, simulation->vz,
        simulation->ax, simulation->ay, simulation->az
    };
    for (size_t field = 0; field < sizeof(fields) / sizeof(fields[0]); ++field) {
        for (size_t i = 0; i < simulation->count; ++i) {
            if (!isfinite(fields[field][i])) {
                return false;
            }
        }
    }
    return true;
}

static void reset_state(Simulation *simulation, const Simulation *initial)
{
    if (initial->count == 0) {
        return;
    }
    double *destination[] = {
        simulation->w, simulation->x, simulation->y, simulation->z,
        simulation->vx, simulation->vy, simulation->vz,
        simulation->ax, simulation->ay, simulation->az
    };
    const double *source[] = {
        initial->w, initial->x, initial->y, initial->z,
        initial->vx, initial->vy, initial->vz,
        initial->ax, initial->ay, initial->az
    };
    for (size_t field = 0; field < sizeof(destination) / sizeof(destination[0]); ++field) {
        memcpy(destination[field], source[field], initial->count * sizeof(double));
    }
}

/* Internal CLI used by run.py. stdout contains timing samples only. */
int main(int argc, char **argv)
{
    double dt;
    size_t steps, samples, warmups;
    if (argc < 8 || argc > 10) {
        fprintf(stderr, "Usage: %s INPUT DT STEPS SAMPLES WARMUPS timing|perf|reference FINAL_FILE [direct|approximate [THETA]]\n", argv[0]);
        return EXIT_FAILURE;
    }
    bool timing = strcmp(argv[6], "timing") == 0;
    bool perf = strcmp(argv[6], "perf") == 0;
    bool reference = strcmp(argv[6], "reference") == 0;
    AccelerationSolver solver = gravity_solver_by_name(argc >= 9 ? argv[8] : "direct");
    ApproximateConfig config = approximate_default_config();
    if ((argc == 10 && (solver != gravity_approximate
                       || !parse_double(argv[9], &config.theta) || config.theta < 0.0))
        || (!timing && !perf && !reference) || !solver
        || !parse_double(argv[2], &dt) || dt <= 0.0
        || !parse_size(argv[3], &steps) || steps == 0
        || !parse_size(argv[4], &samples) || samples == 0
        || samples > SIZE_MAX / sizeof(double)
        || !parse_size(argv[5], &warmups)
        || (reference && (solver != gravity_direct || warmups != 0))
        || (!timing && samples != 1)) {
        fprintf(stderr, "Invalid benchmark arguments.\n");
        return EXIT_FAILURE;
    }

    FILE *input = fopen(argv[1], "r");
    if (!input) {
        perror("input");
        return EXIT_FAILURE;
    }
    Simulation initial;
    bool loaded = read_bodies(input, &initial);
    int closed = fclose(input);
    if (!loaded || closed != 0) {
        fprintf(stderr, "Could not read instance.\n");
        free(initial.storage);
        return EXIT_FAILURE;
    }

    int result = EXIT_FAILURE;
    ApproximateWorkspace workspace = {0};
    void *context = NULL;
    Simulation simulation;
    bool allocated = simulation_create(&simulation, initial.count);
    double *durations = malloc(samples * sizeof(*durations));
    if (!allocated || !durations) {
        fprintf(stderr, "Could not allocate benchmark state.\n");
        goto cleanup;
    }

    if (solver == gravity_approximate) {
        if (!approximate_workspace_create(&workspace, initial.count,
                                          config)) {
            fprintf(stderr, "Could not allocate approximate workspace.\n");
            goto cleanup;
        }
        context = &workspace;
    }

    if (context) {
        fprintf(stderr, "approximate_config: {\"theta\": %.17g, \"bucket_capacity\": %zu, "
                "\"max_depth\": %zu, \"node_capacity\": %zu}\n",
                config.theta, config.bucket_capacity, config.max_depth, workspace.node_capacity);
    }
    for (size_t i = 0; i < warmups; ++i) {
        reset_state(&simulation, &initial);
        simulate_with_solver(&simulation, dt, steps, solver, context);
        if (!finite_state(&simulation)) {
            fprintf(stderr, "Warmup produced a nonfinite state.\n");
            goto cleanup;
        }
    }

    for (size_t i = 0; i < samples; ++i) {
        reset_state(&simulation, &initial);
        struct timespec start, end;
        if (timing) {
            if (clock_gettime(CLOCK_MONOTONIC, &start) != 0) {
                perror("clock_gettime");
                goto cleanup;
            }
        } else if (perf && !perf_gate(true)) {
            perror("enable perf events");
            goto cleanup;
        }

        simulate_with_solver(&simulation, dt, steps, solver, context);

        if (timing) {
            if (clock_gettime(CLOCK_MONOTONIC, &end) != 0) {
                perror("clock_gettime");
                goto cleanup;
            }
            durations[i] = (double)(end.tv_sec - start.tv_sec)
                + (double)(end.tv_nsec - start.tv_nsec) * 1e-9;
        } else if (perf && !perf_gate(false)) {
            perror("disable perf events");
            goto cleanup;
        }
        if (!finite_state(&simulation)) {
            fprintf(stderr, "Sample %zu produced a nonfinite state.\n", i + 1);
            goto cleanup;
        }
    }

    FILE *output = fopen(argv[7], "w");
    if (!output) {
        perror("final state");
        goto cleanup;
    }
    bool written = write_bodies(output, &simulation);
    closed = fclose(output);
    if (!written || closed != 0) {
        fprintf(stderr, "Could not write final state.\n");
        goto cleanup;
    }
    if (timing) {
        puts("sample,seconds");
        for (size_t i = 0; i < samples; ++i) {
            printf("%zu,%.17g\n", i + 1, durations[i]);
        }
    }
    if (fflush(stdout) != 0 || ferror(stdout)) {
        fprintf(stderr, "Could not write measurements.\n");
        goto cleanup;
    }
    result = EXIT_SUCCESS;

cleanup:
    if (context) {
        const ApproximateTotals *totals = &workspace.totals;
        fprintf(stderr, "approximate_totals: {\"scope\": \"warmups_and_samples\", "
                "\"calls\": %zu, \"capacity_limited_leaves\": %zu, "
                "\"bounds_fallbacks\": %zu, \"aggregate_fallbacks\": %zu, "
                "\"accepted_nodes\": %zu}\n", totals->calls,
                totals->capacity_limited_leaves, totals->bounds_fallbacks,
                totals->aggregate_fallbacks, totals->accepted_nodes);
    }
    approximate_workspace_destroy(&workspace);
    free(durations);
    free(simulation.storage);
    free(initial.storage);
    return result;
}
