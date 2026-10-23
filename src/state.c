#include "nbody.h"

#include <stdint.h>
#include <stdlib.h>

bool simulation_create(Simulation *simulation, size_t count)
{
    *simulation = (Simulation){0};
    if (count > SIZE_MAX / (13 * sizeof(double))) {
        return false;
    }
    if (count == 0) {
        return true;
    }
    double *storage = malloc(13 * count * sizeof(*storage));
    if (!storage) {
        return false;
    }
    *simulation = (Simulation){
        .count = count, .storage = storage,
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
