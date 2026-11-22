#include "state_internal.h"
#include <stdint.h>
#include <stdlib.h>

bool simulation_create(Simulation **result, size_t count)
{
    *result = NULL;
    if (count > (SIZE_MAX - sizeof(Simulation)) / (13 * sizeof(double))) {
        return false;
    }
    Simulation *simulation = malloc(sizeof(*simulation) + 13 * count * sizeof(double));
    if (!simulation) {
        return false;
    }
    *result = simulation;
    if (count == 0) {
        *simulation = (Simulation){0};
        return true;
    }
    double *storage = simulation->storage;
    *simulation = (Simulation){
        .count = count,
        .w = storage, .x = storage + count, .y = storage + 2 * count,
        .z = storage + 3 * count, .vx = storage + 4 * count,
        .vy = storage + 5 * count, .vz = storage + 6 * count,
        .ax = storage + 7 * count, .ay = storage + 8 * count,
        .az = storage + 9 * count, .next_ax = storage + 10 * count,
        .next_ay = storage + 11 * count, .next_az = storage + 12 * count
    };
    return true;
}

SimulationView simulation_view(const Simulation *simulation)
{
    return (SimulationView){
        .count = simulation->count,
        .w = simulation->w, .x = simulation->x,
        .y = simulation->y, .z = simulation->z,
        .vx = simulation->vx, .vy = simulation->vy, .vz = simulation->vz,
        .ax = simulation->ax, .ay = simulation->ay, .az = simulation->az
    };
}

void simulation_destroy(Simulation *simulation)
{
    free(simulation);
}

MutableSimulationView simulation_mutable_view(Simulation *simulation)
{
    return (MutableSimulationView){
        simulation->count, simulation->w, simulation->x, simulation->y, simulation->z,
        simulation->vx, simulation->vy, simulation->vz,
        simulation->ax, simulation->ay, simulation->az
    };
}

const char *refresh_acceleration(Simulation *simulation, AccelerationSolver solver, void *context)
{
    const SimulationView state = simulation_view(simulation);
    const char *error = solver(&state, simulation->next_ax, simulation->next_ay,
                               simulation->next_az, context);
    if (error) {
        return error;
    }

    /* Publish only after every interaction has read the whole input state. */
    double *previous_ax = simulation->ax;
    double *previous_ay = simulation->ay;
    double *previous_az = simulation->az;
    simulation->ax = simulation->next_ax;
    simulation->ay = simulation->next_ay;
    simulation->az = simulation->next_az;
    simulation->next_ax = previous_ax;
    simulation->next_ay = previous_ay;
    simulation->next_az = previous_az;
    return NULL;
}
