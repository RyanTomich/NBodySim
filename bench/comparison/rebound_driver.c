#define _POSIX_C_SOURCE 200809L

#include "io.h"
#include "parse.h"
#include <rebound.h>

#include <omp.h>
#include <stdlib.h>
#include <time.h>

#if !defined(_OPENMP) || !defined(OPENMP) || defined(QUADRUPOLE)
#error "The REBOUND comparison requires OpenMP and monopole tree gravity"
#endif

static void add_particles(struct reb_simulation *sim, SimulationView state)
{
    for (size_t i = 0; i < state.count; ++i) {
        char name[32];
        snprintf(name, sizeof(name), "%zu", i);
        struct reb_particle body = {.name = name, .m = state.w[i],
            .x = state.x[i], .y = state.y[i], .z = state.z[i],
            .vx = state.vx[i], .vy = state.vy[i], .vz = state.vz[i]};
        reb_simulation_add(sim, body);
    }
}

static bool write_final(const char *path, struct reb_simulation *sim, Simulation *state)
{
    MutableSimulationView out = simulation_mutable_view(state);
    for (size_t i = 0; i < out.count; ++i) {
        char name[32];
        snprintf(name, sizeof(name), "%zu", i);
        const struct reb_particle *p = reb_simulation_get_particle_by_name(sim, name);
        if (!p || p->m != out.w[i]) {
            fprintf(stderr, "Missing or changed particle %zu\n", i);
            return false;
        }
        out.x[i] = p->x; out.y[i] = p->y; out.z[i] = p->z;
        out.vx[i] = p->vx; out.vy[i] = p->vy; out.vz[i] = p->vz;
        out.ax[i] = p->ax; out.ay[i] = p->ay; out.az[i] = p->az;
    }
    FILE *output = fopen(path, "w");
    if (!output) {
        perror(path);
        return false;
    }
    bool written = write_bodies(output, state);
    return fclose(output) == 0 && written;
}

int main(int argc, char **argv)
{
    double dt;
    size_t steps;
    if (argc != 5 || !parse_double(argv[2], &dt) || dt <= 0
        || !parse_size(argv[3], &steps) || steps == 0) {
        fprintf(stderr, "Usage: %s INPUT DT STEPS FINAL_FILE\n", argv[0]);
        return EXIT_FAILURE;
    }
    FILE *input = fopen(argv[1], "r");
    if (!input) {
        perror(argv[1]);
        return EXIT_FAILURE;
    }
    Simulation *state;
    bool loaded = read_bodies(input, &state);
    if (fclose(input) != 0 || !loaded) {
        fprintf(stderr, "Could not read instance\n");
        simulation_destroy(state);
        return EXIT_FAILURE;
    }

    int result = EXIT_FAILURE;
    struct reb_simulation *sim = reb_simulation_create();
    sim->save_messages = 1;
    struct reb_integrator_leapfrog_state *leapfrog = reb_simulation_set_integrator(sim, "leapfrog");
    leapfrog->order = 2;
    sim->dt = dt;
    sim->G = 1;
    sim->softening = 0;
    sim->gravity = REB_GRAVITY_TREE;
    sim->opening_angle2 = 0.25;
    sim->boundary = REB_BOUNDARY_NONE;
    sim->collision = REB_COLLISION_NONE;
    sim->root_size = 32; /* One isolated cube enclosing the uniform [-10,10] input. */
    SimulationView initial = simulation_view(state);
    add_particles(sim, initial);

    int workers = 0;
#pragma omp parallel
    {
#pragma omp single
        workers = omp_get_num_threads();
    }
    fprintf(stderr, "thread_check: {\"requested_threads\": %d, \"effective_threads\": %d}\n",
            omp_get_max_threads(), workers);
    if (workers != omp_get_max_threads() || omp_get_dynamic()) {
        fprintf(stderr, "Expected a fixed OpenMP team\n");
        goto cleanup;
    }
    struct timespec start, end;
    if (clock_gettime(CLOCK_MONOTONIC, &start) != 0) {
        perror("clock_gettime");
        goto cleanup;
    }
    enum REB_STATUS status = reb_simulation_steps(sim, steps);
    if (clock_gettime(CLOCK_MONOTONIC, &end) != 0) {
        perror("clock_gettime");
        goto cleanup;
    }
    if (sim->messages && sim->messages[0]) {
        fprintf(stderr, "REBOUND: %s\n", sim->messages[0] + 1);
        goto cleanup;
    }
    if (status != REB_STATUS_SUCCESS || sim->steps_done != steps || sim->N != initial.count) {
        fprintf(stderr, "REBOUND failed to finish the requested trajectory\n");
        goto cleanup;
    }
    if (!write_final(argv[4], sim, state)) {
        goto cleanup;
    }
    /* Second-order DKD evaluates gravity once on each completed step. */
    fprintf(stderr, "solver_totals: {\"calls\": %zu}\n", (size_t)sim->steps_done);
    double seconds = (double)(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
    printf("sample,seconds\n1,%.17g\n", seconds);
    if (fflush(stdout) == 0 && !ferror(stdout)) {
        result = EXIT_SUCCESS;
    }
cleanup:
    reb_simulation_free(sim);
    simulation_destroy(state);
    return result;
}
