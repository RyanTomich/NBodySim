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

static inline Vec3 simulation_view_position(const SimulationView *state, size_t body)
{
    return (Vec3){state->x[body], state->y[body], state->z[body]};
}

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
 * Return NULL on success or a static error message on failure. Discard failed outputs.
 * Outputs must be disjoint from each other and state. Context is caller-owned;
 * any solver workspace must be prepared before simulation and reused per step.
 * Accelerations must depend only on current masses and positions, not velocities
 * or previous accelerations. After a drift, visible velocities are half-step values.
 */
typedef const char *(*AccelerationSolver)(const SimulationView *state, double *ax,
                                          double *ay, double *az, void *context);

typedef const char *(*Integrate)(Simulation *, double, size_t, AccelerationSolver, void *);
/* NULL selects the default. Unknown names give a null function pointer. */
Integrate parse_integrator(const char *name);

/* Shared kick-drift-kick leapfrog loop with a whole-system acceleration solver.
 * Positive finite dt required. Each nonempty, positive-step call recomputes
 * initial accelerations from masses and positions, ignoring supplied acceleration
 * values, then evaluates once per step (steps + 1 solver calls in total).
 * Publish new acceleration only after evaluating all bodies; views then expire.
 * All integrators return NULL on success, with full-step velocities and final forces.
 * On force failure, return the callback's error and retain the last published
 * accelerations. Positions and velocities may have advanced: reset before reuse.
 * No allocation, conversion, or validation in this loop.
 * The solver must be nonnull; context may be NULL if the solver permits it.
 * Zero count or zero steps leaves state untouched and does not invoke the solver.
 */
const char *simulate_with_solver(Simulation *simulation, double dt, size_t steps,
                                 AccelerationSolver solver, void *context);

/* Kick-drift-kick for position-dependent forces. Recomputes initial acceleration
 * and makes steps+1 solver calls. No allocation. Zero steps/count is a no-op.
 * Views expire after integration; supplied accelerations are ignored. */
const char *simulate_leapfrog_with_solver(Simulation *simulation, double dt, size_t steps,
                                          AccelerationSolver solver, void *context);

/* Forward Euler uses the supplied acceleration for the first step and invokes
 * the solver once after each step. The same dt, view and no-op rules apply. */
const char *simulate_euler_with_solver(Simulation *simulation, double dt, size_t steps,
                                       AccelerationSolver solver, void *context);

#endif
