#ifndef STATE_INTERNAL_H
#define STATE_INTERNAL_H
#include "nbody.h"

struct Simulation {
    size_t count;
    double *w, *x, *y, *z;
    double *vx, *vy, *vz;
    double *ax, *ay, *az;
    double *next_ax, *next_ay, *next_az;
    double storage[];
};

#endif
