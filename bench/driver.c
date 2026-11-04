#define _POSIX_C_SOURCE 200809L

#include "solver.h"
#include "metadata.h"
#include "io.h"
#include "perf_gate.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Internal CLI used by run.py. stdout contains timing samples only. */
int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--metadata") == 0) {
        return metadata_query(argc - 2, argv + 2);
    }
    double dt;
    size_t steps, samples, warmups;
    if (argc < 8 || argc > 10) {
        fprintf(stderr, "Usage: %s INPUT DT STEPS SAMPLES WARMUPS timing|perf|reference FINAL_FILE [direct|approximate [THETA]]\n", argv[0]);
        return EXIT_FAILURE;
    }
    bool timing = strcmp(argv[6], "timing") == 0;
    bool perf = strcmp(argv[6], "perf") == 0;
    bool reference = strcmp(argv[6], "reference") == 0;
    const char *solver_name = argc >= 9 ? argv[8] : "direct";
    SolverOption option = {"theta", argc == 10 ? argv[9] : NULL};
    if ((!timing && !perf && !reference)
        || !parse_double(argv[2], &dt) || dt <= 0.0
        || !parse_size(argv[3], &steps) || steps == 0
        || !parse_size(argv[4], &samples) || samples == 0
        || samples > SIZE_MAX / sizeof(double)
        || !parse_size(argv[5], &warmups)
        || (reference && (strcmp(solver_name, "direct") != 0 || warmups != 0))
        || (!timing && samples != 1)) {
        fprintf(stderr, "Invalid benchmark arguments.\n");
        return EXIT_FAILURE;
    }

    FILE *input = fopen(argv[1], "r");
    if (!input) {
        perror("input");
        return EXIT_FAILURE;
    }
    Simulation *initial;
    bool loaded = read_bodies(input, &initial);
    int closed = fclose(input);
    if (!loaded || closed != 0) {
        fprintf(stderr, "Could not read instance.\n");
        simulation_destroy(initial);
        return EXIT_FAILURE;
    }

    int result = EXIT_FAILURE;
    Solver solver = {0};
    Simulation *simulation;
    bool allocated = simulation_create(&simulation, simulation_view(initial).count);
    double *durations = malloc(samples * sizeof(*durations));
    if (!allocated || !durations) {
        fprintf(stderr, "Could not allocate benchmark state.\n");
        goto cleanup;
    }

    const char *error;
    if (!solver_create(&solver, solver_name, simulation_view(initial).count,
                       &option, argc == 10 ? 1 : 0, &error)) {
        fprintf(stderr, "%s\n", error);
        goto cleanup;
    }
    metadata_log(stderr, &solver, false);
    for (size_t i = 0; i < warmups; ++i) {
        (void)simulation_copy(simulation, initial);
        simulate_with_solver(simulation, dt, steps, solver.acceleration, solver.context);
        if (!simulation_is_finite(simulation)) {
            fprintf(stderr, "Warmup produced a nonfinite state.\n");
            goto cleanup;
        }
    }

    for (size_t i = 0; i < samples; ++i) {
        (void)simulation_copy(simulation, initial);
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

        simulate_with_solver(simulation, dt, steps, solver.acceleration, solver.context);

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
        if (!simulation_is_finite(simulation)) {
            fprintf(stderr, "Sample %zu produced a nonfinite state.\n", i + 1);
            goto cleanup;
        }
    }

    FILE *output = fopen(argv[7], "w");
    if (!output) {
        perror("final state");
        goto cleanup;
    }
    bool written = write_bodies(output, simulation);
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
    if (solver.descriptor) {
        metadata_log(stderr, &solver, true);
    }
    solver_destroy(&solver);
    free(durations);
    simulation_destroy(simulation);
    simulation_destroy(initial);
    return result;
}
