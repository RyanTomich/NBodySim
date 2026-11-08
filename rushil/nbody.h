#ifndef NBODY_H
#define NBODY_H

#include <stdbool.h>
#include <stddef.h>

typedef struct {
    double x, y, z;
} Vec3;

typedef struct Simulation Simulation;

/* Borrowed read-only view of the separate coordinate arrays. */
typedef struct {
    size_t count;
    const double *w, *x, *y, *z;
    const double *vx, *vy, *vz;
    const double *ax, *ay, *az;
} SimulationView;

/* One allocation owns the private header and ten contiguous arrays. Initialize the ten state arrays before use.
 * Creation failure sets *simulation to NULL; zero count is valid.
 * Destroy accepts NULL. Handles must not be copied as independent owners.
 */
bool simulation_create(Simulation **simulation, size_t count);
void simulation_destroy(Simulation *simulation);
SimulationView simulation_view(const Simulation *simulation);

typedef struct {
    size_t count;
    double *w, *x, *y, *z;
    double *vx, *vy, *vz;
    double *ax, *ay, *az;
} MutableSimulationView;

/* Borrowed writable arrays, valid until simulation_destroy. */
MutableSimulationView simulation_mutable_view(Simulation *simulation);

#endif
