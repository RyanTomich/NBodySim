#ifndef PERF_GATE_H
#define PERF_GATE_H

#include <stdbool.h>

/* Toggle externally attached, initially disabled perf events for this task. */
bool perf_gate(bool enable);

#endif
