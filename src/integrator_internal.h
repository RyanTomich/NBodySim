#ifndef INTEGRATOR_INTERNAL_H
#define INTEGRATOR_INTERNAL_H

#include "nbody.h"

typedef void (*Integrate)(Simulation *, double, size_t, AccelerationSolver, void *);
/* NULL selects the default. Unknown names give a null function pointer. */
Integrate parse_integrator(const char *name);

#endif
