#include "gravity.h"

/* Replace this body with the approximation algorithm when it is ready.
 * Keep the direct implementation independent for reference evaluation.
 */
void gravity_approximate(const SimulationView *state, double *ax, double *ay,
                         double *az, void *context)
{
    gravity_direct(state, ax, ay, az, context);
}
