#include "state_internal.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

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

bool simulation_copy(Simulation *destination, const Simulation *source)
{
    if (destination->count != source->count) {
        return false;
    }
    if (destination == source || source->count == 0) {
        return true;
    }
    double *out[] = {destination->w, destination->x, destination->y, destination->z,
                    destination->vx, destination->vy, destination->vz,
                    destination->ax, destination->ay, destination->az};
    const double *in[] = {source->w, source->x, source->y, source->z,
                         source->vx, source->vy, source->vz,
                         source->ax, source->ay, source->az};
    for (size_t i = 0; i < 10; ++i) {
        memcpy(out[i], in[i], source->count * sizeof(double));
    }
    return true;
}

bool simulation_is_finite(const Simulation *simulation)
{
    SimulationView view = simulation_view(simulation);
    const double *fields[] = {view.w, view.x, view.y, view.z, view.vx, view.vy, view.vz,
                             view.ax, view.ay, view.az};
    for (size_t field = 0; field < 10; ++field) {
        for (size_t i = 0; i < view.count; ++i) {
            if (!isfinite(fields[field][i])) {
                return false;
            }
        }
    }
    return true;
}
