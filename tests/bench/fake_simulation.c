/* Test double: confirms warmups are ungated and every call sees initial state. */
#include "nbody.h"
#include "perf_gate.h"

#include <assert.h>
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

bool simulate(Body *bodies, size_t count, double dt, size_t steps, Interaction f)
{
    (void)dt;
    (void)steps;
    (void)f;
    assert(count == 1);
    assert(bodies[0].position.x == 0.0);
    assert(bodies[0].velocity.x == 0.0);
    assert(bodies[0].acceleration.x == 0.0);
    fprintf(stderr, "simulate:%s\n", enabled ? "measured" : "ungated");
    bodies[0].position.x = 1.0;
    bodies[0].velocity.x = 2.0;
    bodies[0].acceleration.x = 3.0;
    return getenv("BENCH_TEST_FAIL") == NULL;
}
