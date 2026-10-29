#ifndef NBODY_H
#define NBODY_H

#include <stdbool.h>
#include <stddef.h>

typedef struct {
    double x, y, z;
} Vec3;

typedef struct {
    size_t count;
    double *storage;
    double *w, *x, *y, *z;
    double *vx, *vy, *vz;
    double *ax, *ay, *az;
    double *next_ax, *next_ay, *next_az;
} Simulation;

/* Read-only state, excluding ownership and acceleration workspace.
 * Reacquire after simulate: acceleration pointers can change roles.
 */
typedef struct {
    size_t count;
    const double *w, *x, *y, *z;
    const double *vx, *vy, *vz;
    const double *ax, *ay, *az;
} SimulationView;

/* Allocate thirteen contiguous count-element arrays initially in declaration order.
 * Initialize the ten state arrays before simulation; workspace is overwritten.
 * The destination must not own an allocation. Failure leaves it empty.
 * Zero count succeeds with null pointers. Release with free(simulation.storage).
 * Copies of Simulation share storage; they do not transfer or duplicate ownership.
 */
bool simulation_create(Simulation *simulation, size_t count);
/* Capture the current named state pointers, not a copy of the array contents. */
SimulationView simulation_view(const Simulation *simulation);

/* Return the source's acceleration contribution to the target, without mutation. */
typedef Vec3 (*Interaction)(const SimulationView *state, size_t target, size_t source);

/* Fill all three count-element outputs without modifying state or allocating.
 * Outputs must be disjoint from each other and state. Context is caller-owned;
 * any solver workspace must be prepared before simulation and reused per step.
 */
typedef void (*AccelerationSolver)(const SimulationView *state, double *ax,
                                   double *ay, double *az, void *context);

/* Shared forward Euler loop with a whole-system acceleration solver.
 * Same storage, timestep, update order, and pointer lifetime contracts as simulate.
 * The solver must be nonnull; context may be NULL if the solver permits it.
 * Zero count or zero steps does not invoke the solver.
 */
void simulate_with_solver(Simulation *simulation, double dt, size_t steps,
                          AccelerationSolver solver, void *context);

/* Update in place using forward Euler, then recompute acceleration with f.
 * Input accelerations drive the first step. No allocation or argument validation.
 * Each step swaps ax/ay/az with next_ax/next_ay/next_az after all interactions.
 * Read state through the named pointers; reacquire cached views and acceleration
 * pointers after this call. The owning storage pointer never changes.
 * Requires storage from simulation_create, initialized state, positive finite dt,
 * and a nonnull callback. Zero count or zero steps is a no-op.
 */
void simulate(Simulation *simulation, double dt, size_t steps, Interaction f);

#endif
