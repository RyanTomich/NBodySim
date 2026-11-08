#include "solver.h"
#include "io.h"
#include "parse.h"
#include "nbody.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    double dt;
    size_t steps;
    const char *positional[3];
    size_t positional_count = 0;
    const char *solver_name = "n2";
    SolverOption options[2];
    size_t option_count = 0;
    bool solver_set = false;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--solver") == 0) {
            if (solver_set || ++i == argc) {
                fprintf(stderr, "Expected --solver n2|barnes-hut|fmm once.\n");
                return EXIT_FAILURE;
            }
            solver_name = argv[i];
            solver_set = true;
        } else if (strcmp(argv[i], "--theta") == 0 || strcmp(argv[i], "--order") == 0) {
            const char *name = argv[i] + 2;
            if (++i == argc || option_count == 2) {
                fprintf(stderr, "Invalid solver options.\n");
                return EXIT_FAILURE;
            }
            options[option_count++] = (SolverOption){name, argv[i]};
        } else if (positional_count < 3) {
            positional[positional_count++] = argv[i];
        } else {
            fprintf(stderr, "Unexpected argument: %s\n", argv[i]);
            return EXIT_FAILURE;
        }
    }
    if (positional_count != 3) {
        fprintf(stderr, "Usage: %s [--solver n2|barnes-hut|fmm] [--theta VALUE | --order P] INPUT DT STEPS\n", argv[0]);
        return EXIT_FAILURE;
    }
    if (!parse_double(positional[1], &dt) || dt <= 0.0 || !parse_size(positional[2], &steps)) {
        fprintf(stderr, "DT must be positive and finite; STEPS must be a nonnegative integer.\n");
        return EXIT_FAILURE;
    }

    FILE *input = fopen(positional[0], "r");
    if (!input) {
        perror(positional[0]);
        return EXIT_FAILURE;
    }
    Simulation *simulation;
    bool loaded = read_bodies(input, &simulation);
    int closed = fclose(input);
    if (!loaded || closed != 0) {
        fprintf(stderr, "Could not read instance: invalid data, I/O error, or allocation failure.\n");
        simulation_destroy(simulation);
        return EXIT_FAILURE;
    }

    Solver solver;
    const char *error;
    if (!solver_create(&solver, solver_name, simulation_view(simulation).count,
                       options, option_count, &error)) {
        fprintf(stderr, "%s\n", error);
        simulation_destroy(simulation);
        return EXIT_FAILURE;
    }
    simulate_with_solver(simulation, dt, steps, solver.acceleration, solver.context);
    solver_destroy(&solver);
    bool written = write_bodies(stdout, simulation);
    simulation_destroy(simulation);
    if (!written || fflush(stdout) != 0) {
        fprintf(stderr, "Could not write final state.\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
