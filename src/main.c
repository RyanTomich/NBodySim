#include "gravity.h"
#include "io.h"
#include "nbody.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    double dt;
    size_t steps;
    if (argc != 4) {
        fprintf(stderr, "Usage: %s INPUT DT STEPS\n", argv[0]);
        return EXIT_FAILURE;
    }
    if (!parse_double(argv[2], &dt) || dt <= 0.0 || !parse_size(argv[3], &steps)) {
        fprintf(stderr, "DT must be positive and finite; STEPS must be a nonnegative integer.\n");
        return EXIT_FAILURE;
    }

    FILE *input = fopen(argv[1], "r");
    if (!input) {
        perror(argv[1]);
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

    simulate(&simulation, dt, steps, gravity);
    bool written = write_bodies(stdout, &simulation);
    free(simulation.storage);
    if (!written || fflush(stdout) != 0) {
        fprintf(stderr, "Could not write final state.\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
