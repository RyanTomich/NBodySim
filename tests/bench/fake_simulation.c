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

void simulate_with_solver(Simulation *simulation, double dt, size_t steps,
                          AccelerationSolver solver, void *context)
{
    (void)dt;
    (void)steps;
    (void)solver;
    (void)context;
    MutableSimulationView state = simulation_mutable_view(simulation);
    assert(state.count == 1);
    double *fields[] = {
        state.w, state.x, state.y, state.z,
        state.vx, state.vy, state.vz,
        state.ax, state.ay, state.az
    };
    for (size_t i = 0; i < 10; ++i) {
        assert(fields[i][0] == (i == 0 ? 1.0 : 0.0));
    }
    fprintf(stderr, "simulate:%s\n", enabled ? "measured" : "ungated");
    /* Change every field so incomplete resets are observable on the next call. */
    for (size_t i = 0; i < 10; ++i) {
        fields[i][0] = (double)(i + 1);
    }
    const char *failure = getenv("BENCH_TEST_FAIL");
    if (failure) {
        assert(failure[0] >= '1' && failure[0] <= '3' && failure[1] == '\0');
        fields[7 + (failure[0] - '1')][0] = NAN;
    }
}
