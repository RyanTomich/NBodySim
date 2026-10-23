#define _POSIX_C_SOURCE 200809L

#include "gravity.h"
#include "io.h"
#include "perf_gate.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static bool finite_state(const Body *bodies, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        const Body *b = &bodies[i];
        if (!isfinite(b->w)
            || !isfinite(b->position.x) || !isfinite(b->position.y) || !isfinite(b->position.z)
            || !isfinite(b->velocity.x) || !isfinite(b->velocity.y) || !isfinite(b->velocity.z)
            || !isfinite(b->acceleration.x) || !isfinite(b->acceleration.y) || !isfinite(b->acceleration.z)) {
            return false;
        }
    }
    return true;
}

static void reset_state(Body *bodies, const Body *initial, size_t count)
{
    if (count != 0) {
        memcpy(bodies, initial, count * sizeof(*bodies));
    }
}

/* Internal CLI used by run.py. stdout contains timing samples only. */
int main(int argc, char **argv)
{
    double dt;
    size_t steps, samples, warmups;
    if (argc != 8) {
        fprintf(stderr, "Usage: %s INPUT DT STEPS SAMPLES WARMUPS timing|perf FINAL_FILE\n", argv[0]);
        return EXIT_FAILURE;
    }
    bool timing = strcmp(argv[6], "timing") == 0;
    if ((!timing && strcmp(argv[6], "perf") != 0)
        || !parse_double(argv[2], &dt) || dt <= 0.0
        || !parse_size(argv[3], &steps) || steps == 0
        || !parse_size(argv[4], &samples) || samples == 0
        || samples > SIZE_MAX / sizeof(double)
        || !parse_size(argv[5], &warmups)
        || (!timing && samples != 1)) {
        fprintf(stderr, "Invalid benchmark arguments.\n");
        return EXIT_FAILURE;
    }

    FILE *input = fopen(argv[1], "r");
    if (!input) {
        perror("input");
        return EXIT_FAILURE;
    }
    Body *initial;
    size_t count;
    bool loaded = read_bodies(input, &initial, &count);
    int closed = fclose(input);
    if (!loaded || closed != 0) {
        fprintf(stderr, "Could not read instance.\n");
        free(initial);
        return EXIT_FAILURE;
    }

    int result = EXIT_FAILURE;
    Body *bodies = count == 0 ? NULL : malloc(count * sizeof(*bodies));
    double *durations = malloc(samples * sizeof(*durations));
    if ((count != 0 && !bodies) || !durations) {
        fprintf(stderr, "Could not allocate benchmark state.\n");
        goto cleanup;
    }

    for (size_t i = 0; i < warmups; ++i) {
        reset_state(bodies, initial, count);
        if (!simulate(bodies, count, dt, steps, gravity) || !finite_state(bodies, count)) {
            fprintf(stderr, "Warmup failed or produced a nonfinite state.\n");
            goto cleanup;
        }
    }

    for (size_t i = 0; i < samples; ++i) {
        reset_state(bodies, initial, count);
        struct timespec start, end;
        if (timing) {
            if (clock_gettime(CLOCK_MONOTONIC, &start) != 0) {
                perror("clock_gettime");
                goto cleanup;
            }
        } else if (!perf_gate(true)) {
            perror("enable perf events");
            goto cleanup;
        }

        bool completed = simulate(bodies, count, dt, steps, gravity);

        if (timing) {
            if (clock_gettime(CLOCK_MONOTONIC, &end) != 0) {
                perror("clock_gettime");
                goto cleanup;
            }
            durations[i] = (double)(end.tv_sec - start.tv_sec)
                + (double)(end.tv_nsec - start.tv_nsec) * 1e-9;
        } else if (!perf_gate(false)) {
            perror("disable perf events");
            goto cleanup;
        }
        if (!completed || !finite_state(bodies, count)) {
            fprintf(stderr, "Sample %zu failed or produced a nonfinite state.\n", i + 1);
            goto cleanup;
        }
    }

    FILE *output = fopen(argv[7], "w");
    if (!output) {
        perror("final state");
        goto cleanup;
    }
    bool written = write_bodies(output, bodies, count);
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
    free(durations);
    free(bodies);
    free(initial);
    return result;
}
