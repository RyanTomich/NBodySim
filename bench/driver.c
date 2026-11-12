#define _POSIX_C_SOURCE 200809L

#include "solver.h"
#include "io.h"
#include "parse.h"
#include "step_timing.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef NBODY_BENCH_OPENMP
#ifndef _OPENMP
#error "The parallel solver sweep requires OpenMP support"
#endif
#include <omp.h>
#endif

#ifdef NBODY_BENCH_STEPS
#if defined(_OPENMP) && !defined(NBODY_BENCH_OPENMP)
#error "The single-run solver sweep requires a serial build"
#endif
static struct timespec *step_boundaries;
static bool step_clock_failed;

void benchmark_step_boundary(size_t boundary)
{
    if (clock_gettime(CLOCK_MONOTONIC, &step_boundaries[boundary]) != 0) {
        step_clock_failed = true;
    }
}
#endif

/* Reset logical state through fresh views; inactive workspace is overwritten by
 * the solver. The driver owns distinct, equally sized initial and working states. */
static void reset_state(Simulation *destination, const Simulation *source)
{
    MutableSimulationView out = simulation_mutable_view(destination);
    SimulationView in = simulation_view(source);
    assert(destination != source && out.count == in.count);
    if (in.count == 0) {
        return;
    }
    double *outputs[] = {out.w, out.x, out.y, out.z, out.vx, out.vy, out.vz,
                         out.ax, out.ay, out.az};
    const double *inputs[] = {in.w, in.x, in.y, in.z, in.vx, in.vy, in.vz,
                             in.ax, in.ay, in.az};
    for (size_t field = 0; field < 10; ++field) {
        memcpy(outputs[field], inputs[field], in.count * sizeof(double));
    }
}

static bool state_is_finite(SimulationView state)
{
    const double *fields[] = {state.w, state.x, state.y, state.z,
                             state.vx, state.vy, state.vz,
                             state.ax, state.ay, state.az};
    for (size_t field = 0; field < 10; ++field) {
        for (size_t i = 0; i < state.count; ++i) {
            if (!isfinite(fields[field][i])) {
                return false;
            }
        }
    }
    return true;
}

/* Report settings from the solver backend, outside measurement. */
typedef struct {
    FILE *output;
    size_t fields;
    const char *prefix; /* Emit a log prefix only if a field is reported. */
} JsonWriter;

static void json_string(FILE *output, const char *text)
{
    fputc('"', output);
    for (const unsigned char *p = (const unsigned char *)text; *p; ++p) {
        if (*p == '"' || *p == '\\') {
            fputc('\\', output);
        }
        if (*p < 32) {
            fprintf(output, "\\u%04x", *p);
        } else {
            fputc(*p, output);
        }
    }
    fputc('"', output);
}

static void json_field(const char *name, MetadataValue value, void *visitor_context)
{
    JsonWriter *writer = visitor_context;
    FILE *output = writer->output;
    if (!writer->fields && writer->prefix) {
        fputs(writer->prefix, output);
    }
    if (writer->fields++) {
        fputs(", ", output);
    }
    json_string(output, name);
    fputs(": ", output);
    switch (value.type) {
    case METADATA_SIZE: fprintf(output, "%zu", value.value.size); break;
    case METADATA_DOUBLE: fprintf(output, "%.17g", value.value.real); break;
    case METADATA_BOOL: fputs(value.value.boolean ? "true" : "false", output); break;
    case METADATA_STRING: json_string(output, value.value.string); break;
    }
}

static void metadata_log(FILE *output, const Solver *solver, bool diagnostics)
{
    JsonWriter writer = {output, 0, diagnostics ? "solver_totals: {" : "solver_config: {"};
    if (diagnostics) {
        solver_diagnostics(solver, json_field, &writer);
    } else {
        solver_settings(solver, json_field, &writer);
    }
    if (writer.fields) {
        if (diagnostics) {
            json_field("scope", (MetadataValue){METADATA_STRING, {.string = "warmups_and_samples"}}, &writer);
        }
        fputs("}\n", output);
    }
}

#ifdef NBODY_PROFILE_BH
/* Snapshots stay outside the measured call; differences exclude warmups.
 * Use the public metadata interface, not backend workspace representation. */
static void profile_snapshot(const Solver *solver, size_t sample, const char *boundary)
{
    JsonWriter writer = {.output = stderr};
    fputs("phase_sample: {", stderr);
    json_field("sample", (MetadataValue){METADATA_SIZE, {.size = sample}}, &writer);
    json_field("boundary", (MetadataValue){METADATA_STRING, {.string = boundary}}, &writer);
    solver_diagnostics(solver, json_field, &writer);
    fputs("}\n", stderr);
}
#endif

static int metadata_query(int argc, char **argv)
{
    if (argc == 0) {
        fputc('[', stdout);
        for (size_t i = 0; i < solver_count(); ++i) {
            if (i) {
                fputc(',', stdout);
            }
            json_string(stdout, solver_name_at(i));
        }
        puts("]");
    } else {
        size_t count;
        if ((argc != 2 && argc != 4) || !parse_size(argv[1], &count)) {
            return EXIT_FAILURE;
        }
        Solver solver;
        const char *error;
        SolverOption option = {argc == 4 ? argv[2] : NULL, argc == 4 ? argv[3] : NULL};
        if (!solver_create(&solver, argv[0], count, &option, argc == 4 ? 1 : 0, &error)) {
            fprintf(stderr, "%s\n", error);
            return EXIT_FAILURE;
        }
        JsonWriter writer = {.output = stdout};
        fputc('{', stdout);
        solver_settings(&solver, json_field, &writer);
        puts("}");
        solver_destroy(&solver);
    }
    return fflush(stdout) == 0 && !ferror(stdout) ? EXIT_SUCCESS : EXIT_FAILURE;
}

/* Internal CLI used by benchmark runners. */
int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--build-info") == 0) {
        puts("{\"openmp\": "
#ifdef _OPENMP
             "true"
#else
             "false"
#endif
             ", \"instrumented\": "
#if defined(NBODY_BENCH_STEPS) || defined(NBODY_PROFILE_BH)
             "true"
#else
             "false"
#endif
             "}");
        return fflush(stdout) == 0 && !ferror(stdout) ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (argc >= 2 && strcmp(argv[1], "--metadata") == 0) {
        return metadata_query(argc - 2, argv + 2);
    }
    double dt;
    size_t steps, samples, warmups;
    if (argc < 8 || (argc > 9 && argc != 11)) {
        fprintf(stderr, "Usage: %s INPUT DT STEPS SAMPLES WARMUPS timing|reference FINAL_FILE [SOLVER [OPTION VALUE]]\n", argv[0]);
        return EXIT_FAILURE;
    }
    bool timing = strcmp(argv[6], "timing") == 0;
    bool reference = strcmp(argv[6], "reference") == 0;
    const char *solver_name = argc >= 9 ? argv[8] : "n2";
    SolverOption option = {argc == 11 ? argv[9] : NULL, argc == 11 ? argv[10] : NULL};
    if ((!timing && !reference)
        || !parse_double(argv[2], &dt) || dt <= 0.0
        || !parse_size(argv[3], &steps) || steps == 0
        || !parse_size(argv[4], &samples) || samples == 0
        || samples > SIZE_MAX / sizeof(double)
        || !parse_size(argv[5], &warmups)
        || (reference && (strcmp(solver_name, "n2") != 0 || warmups != 0))
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
#ifdef NBODY_BENCH_STEPS
    if (steps >= SIZE_MAX / sizeof(*step_boundaries)
        || !(step_boundaries = calloc(steps + 1, sizeof(*step_boundaries)))) {
        fprintf(stderr, "Could not allocate step timing storage.\n");
        goto cleanup;
    }
#endif
    if (!allocated || !durations) {
        fprintf(stderr, "Could not allocate benchmark state.\n");
        goto cleanup;
    }

    const char *error;
#ifdef NBODY_BENCH_OPENMP
    /* Check runtime team formation before preparation and outside timing. */
    const int requested_threads = omp_get_max_threads();
    int effective_threads = 1;
#pragma omp parallel shared(effective_threads)
    {
#pragma omp single
        effective_threads = omp_get_num_threads();
    }
    fprintf(stderr, "thread_check: {\"openmp\": true, \"requested_threads\": %d, \"effective_threads\": %d}\n",
            requested_threads, effective_threads);
    if (effective_threads != requested_threads || omp_get_dynamic()) {
        fprintf(stderr, "Expected a fixed OpenMP team of %d threads.\n", requested_threads);
        goto cleanup;
    }
#endif
    if (!solver_create(&solver, solver_name, simulation_view(initial).count,
                       &option, argc == 11 ? 1 : 0, &error)) {
        fprintf(stderr, "%s\n", error);
        goto cleanup;
    }
    metadata_log(stderr, &solver, false);
    for (size_t i = 0; i < warmups; ++i) {
        reset_state(simulation, initial);
        simulate_with_solver(simulation, dt, steps, solver.acceleration, solver.context);
        if (!state_is_finite(simulation_view(simulation))) {
            fprintf(stderr, "Warmup produced a nonfinite state.\n");
            goto cleanup;
        }
    }

    for (size_t i = 0; i < samples; ++i) {
        reset_state(simulation, initial);
        struct timespec start, end;
#ifdef NBODY_PROFILE_BH
        if (timing) {
            profile_snapshot(&solver, i + 1, "before");
        }
#endif
        if (timing) {
            if (clock_gettime(CLOCK_MONOTONIC, &start) != 0) {
                perror("clock_gettime");
                goto cleanup;
            }
        }

        simulate_with_solver(simulation, dt, steps, solver.acceleration, solver.context);

        if (timing) {
            if (clock_gettime(CLOCK_MONOTONIC, &end) != 0) {
                perror("clock_gettime");
                goto cleanup;
            }
            durations[i] = (double)(end.tv_sec - start.tv_sec)
                + (double)(end.tv_nsec - start.tv_nsec) * 1e-9;
        }
#ifdef NBODY_BENCH_STEPS
        if (step_clock_failed) {
            fprintf(stderr, "Step clock failed.\n");
            goto cleanup;
        }
        if (timing && simulation_view(simulation).count) {
            for (size_t step = 0; step < steps; ++step) {
                double seconds = (double)(step_boundaries[step + 1].tv_sec - step_boundaries[step].tv_sec)
                    + (double)(step_boundaries[step + 1].tv_nsec - step_boundaries[step].tv_nsec) * 1e-9;
                fprintf(stderr, "step_timing: {\"sample\": %zu, \"step\": %zu, \"seconds\": %.17g}\n",
                        i + 1, step + 1, seconds);
            }
        }
#endif
#ifdef NBODY_PROFILE_BH
        if (timing) {
            profile_snapshot(&solver, i + 1, "after");
        }
#endif
        if (!state_is_finite(simulation_view(simulation))) {
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
#ifdef NBODY_BENCH_STEPS
    free(step_boundaries);
#endif
    if (solver.descriptor) {
        metadata_log(stderr, &solver, true);
    }
    solver_destroy(&solver);
    free(durations);
    simulation_destroy(simulation);
    simulation_destroy(initial);
    return result;
}
