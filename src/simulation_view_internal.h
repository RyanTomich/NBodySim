#ifndef SIMULATION_VIEW_INTERNAL_H
#define SIMULATION_VIEW_INTERNAL_H

#include "nbody.h"

static inline Vec3 simulation_view_position(const SimulationView *state, size_t body)
{
    return (Vec3){state->x[body], state->y[body], state->z[body]};
}

#endif
