#ifndef RUSHIL_PHYSICS_H
#define RUSHIL_PHYSICS_H

#include "nbody.h"

/* Initialize shared Simulation arrays: equal masses, random positions, zero v/a. */
void init_bodies(Simulation *simulation);

/* Overwrite position [0..2] and velocity [3..5]; requires nonzero total mass. */
void calc_cm(double cm_vals[6], const Simulation *simulation);

#endif
