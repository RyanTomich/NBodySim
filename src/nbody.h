#ifndef NBODY_H
#define NBODY_H

#include <stdbool.h>
#include <stddef.h>

typedef struct {
    double x, y, z;
} Vec3;

typedef struct Simulation Simulation;

/* Read-only state, excluding ownership and acceleration workspace.
 * Reacquire after simulate_with_solver: acceleration pointers can change roles.
 */
typedef struct {
    size_t count;
    const double *w, *x, *y, *z;
    const double *vx, *vy, *vz;
    const double *ax, *ay, *az;
} SimulationView;

/* One allocation owns the private header and thirteen contiguous arrays. Initialize the ten state arrays before use.
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

/* Borrowed arrays; reacquire views after integration swaps accelerations. */
MutableSimulationView simulation_mutable_view(Simulation *simulation);

/* Fill all three count-element outputs without modifying state or allocating.
 * Outputs must be disjoint from each other and state. Context is caller-owned;
 * any solver workspace must be prepared before simulation and reused per step.
 */
typedef void (*AccelerationSolver)(const SimulationView *state, double *ax,
                                   double *ay, double *az, void *context);

/* Shared forward Euler loop with a whole-system acceleration solver.
 * Positive finite dt required. Input accelerations drive the first step.
 * Publish new acceleration only after evaluating all bodies; views then expire.
 * No allocation, conversion, or validation in this loop.
 * The solver must be nonnull; context may be NULL if the solver permits it.
 * Zero count or zero steps does not invoke the solver.
 */
void simulate_with_solver(Simulation *simulation, double dt, size_t steps,
                          AccelerationSolver solver, void *context);

/* Kick-drift-kick for position-dependent forces. Recomputes initial acceleration
 * and makes steps+1 solver calls. No allocation. Zero steps/count is a no-op.
 * Views expire after integration; supplied accelerations are ignored. */
void simulate_leapfrog_with_solver(Simulation *simulation, double dt, size_t steps,
                                   AccelerationSolver solver, void *context);

#endif
