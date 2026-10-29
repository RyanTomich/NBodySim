#include "gravity.h"
#include "io.h"
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
    AccelerationSolver solver = gravity_direct;
    bool solver_set = false;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--solver") == 0) {
            if (solver_set || ++i == argc || !(solver = gravity_solver_by_name(argv[i]))) {
                fprintf(stderr, "Expected --solver direct|approximate once.\n");
                return EXIT_FAILURE;
            }
            solver_set = true;
        } else if (positional_count < 3) {
            positional[positional_count++] = argv[i];
        } else {
            fprintf(stderr, "Unexpected argument: %s\n", argv[i]);
            return EXIT_FAILURE;
        }
    }
    if (positional_count != 3) {
        fprintf(stderr, "Usage: %s [--solver direct|approximate] INPUT DT STEPS\n", argv[0]);
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
    Simulation simulation;
    bool loaded = read_bodies(input, &simulation);
    int closed = fclose(input);
    if (!loaded || closed != 0) {
        fprintf(stderr, "Could not read instance: invalid data, I/O error, or allocation failure.\n");
        free(simulation.storage);
        return EXIT_FAILURE;
    }

    simulate_with_solver(&simulation, dt, steps, solver, NULL);
    bool written = write_bodies(stdout, &simulation);
    free(simulation.storage);
    if (!written || fflush(stdout) != 0) {
        fprintf(stderr, "Could not write final state.\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
