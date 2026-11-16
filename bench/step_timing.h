#ifndef BENCH_STEP_TIMING_H
#define BENCH_STEP_TIMING_H

#include <stddef.h>

/* Supplied only by the instrumented benchmark executable. Ordinary simulation
 * builds have no observer, clock calls, or extra branch in the stepping loop. */
#ifdef NBODY_BENCH_STEPS
void benchmark_step_boundary(size_t boundary);
#endif

#endif
