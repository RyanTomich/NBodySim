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
    SolverOption theta = {"theta", NULL};
    bool solver_set = false;
    const char *integrator_name = NULL;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--solver") == 0) {
            if (solver_set || ++i == argc) {
                fprintf(stderr, "Expected --solver NAME once.\n");
                return EXIT_FAILURE;
            }
            solver_name = argv[i];
            solver_set = true;
        } else if (strcmp(argv[i], "--integrator") == 0) {
            if (integrator_name || ++i == argc) {
                fprintf(stderr, "Expected --integrator leapfrog|euler once.\n");
                return EXIT_FAILURE;
            }
            integrator_name = argv[i];
        } else if (strcmp(argv[i], "--theta") == 0) {
            if (theta.value || ++i == argc) {
                fprintf(stderr, "Expected --theta VALUE once.\n");
                return EXIT_FAILURE;
            }
            theta.value = argv[i];
        } else if (positional_count < 3) {
            positional[positional_count++] = argv[i];
        } else {
            fprintf(stderr, "Unexpected argument: %s\n", argv[i]);
            return EXIT_FAILURE;
        }
    }
    if (positional_count != 3) {
        fprintf(stderr, "Usage: %s [--solver NAME] [--integrator leapfrog|euler] [--theta VALUE] INPUT DT STEPS\n",
                argv[0]);
        fputs("Solvers:", stderr);
        for (size_t i = 0; i < solver_count(); ++i) {
            fprintf(stderr, " %s", solver_name_at(i));
        }
        fputc('\n', stderr);
        return EXIT_FAILURE;
    }
    Integrate integrate = parse_integrator(integrator_name);
    if (!integrate) {
        fprintf(stderr, "Integrator must be leapfrog or euler.\n");
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
                       &theta, theta.value != NULL, &error)) {
        fprintf(stderr, "%s\n", error);
        simulation_destroy(simulation);
        return EXIT_FAILURE;
    }
    error = integrate(simulation, dt, steps, solver.acceleration, solver.context);
    solver_destroy(&solver);
    if (error) {
        fprintf(stderr, "%s\n", error);
        simulation_destroy(simulation);
        return EXIT_FAILURE;
    }
    bool written = write_bodies(stdout, simulation);
    simulation_destroy(simulation);
    if (!written || fflush(stdout) != 0) {
        fprintf(stderr, "Could not write final state.\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
