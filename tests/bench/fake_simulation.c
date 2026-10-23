/* Test double: confirms warmups are ungated and every call sees initial state. */
#include "nbody.h"
#include "perf_gate.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static bool enabled;

bool perf_gate(bool enable)
{
    assert(enable != enabled);
    enabled = enable;
    fprintf(stderr, "%s\n", enable ? "enable" : "disable");
    return true;
}

void simulate(Simulation *simulation, double dt, size_t steps, Interaction f)
{
    (void)dt;
    (void)f;
    assert(simulation->count == 1);
    double *fields[] = {
        simulation->w, simulation->x, simulation->y, simulation->z,
        simulation->vx, simulation->vy, simulation->vz,
        simulation->ax, simulation->ay, simulation->az
    };
    for (size_t i = 0; i < 10; ++i) {
        assert(fields[i][0] == (i == 0 ? 1.0 : 0.0));
    }
    fprintf(stderr, "simulate:%s\n", enabled ? "measured" : "ungated");
    if (steps % 2 != 0) {
        simulation->ax = simulation->next_ax;
        simulation->ay = simulation->next_ay;
        simulation->az = simulation->next_az;
        simulation->next_ax = fields[7];
        simulation->next_ay = fields[8];
        simulation->next_az = fields[9];
        fields[7] = simulation->ax;
        fields[8] = simulation->ay;
        fields[9] = simulation->az;
    }
    /* Change every field so incomplete resets are observable on the next call. */
    for (size_t i = 0; i < 10; ++i) {
        fields[i][0] = (double)(i + 1);
    }
    const char *failure = getenv("BENCH_TEST_FAIL");
    /* A failure must come from active state, not an inactive NaN. */
    simulation->next_ax[0] = failure ? 0.0 : NAN;
    simulation->next_ay[0] = failure ? 0.0 : NAN;
    simulation->next_az[0] = failure ? 0.0 : NAN;
    if (failure) {
        assert(failure[0] >= '1' && failure[0] <= '3' && failure[1] == '\0');
        fields[7 + (failure[0] - '1')][0] = NAN;
    }
}
