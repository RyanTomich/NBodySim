#ifndef SOLVER_INTERNAL_H
#define SOLVER_INTERNAL_H
#include "solver.h"

/* create cleans up partial allocations on failure and supplies a static error.
 * Metadata callbacks emit scalar values only; they perform no serialization.
 * Only name and acceleration are required. Without create, options are rejected
 * and solver_context is NULL. Other absent callbacks do nothing.
 * All lifecycle and metadata calls occur outside integration.
 */
struct SolverDescriptor {
    const char *name;
    AccelerationSolver acceleration;
    bool (*create)(void **solver_context, size_t count, const SolverOption *options,
                   size_t option_count, const char **error);
    void (*destroy)(void *solver_context);
    void (*settings)(const void *solver_context, MetadataVisitor visit, void *visitor_context);
    void (*diagnostics)(const void *solver_context, MetadataVisitor visit, void *visitor_context);
};
extern const SolverDescriptor n2_solver;
extern const SolverDescriptor barnes_hut_solver;
extern const SolverDescriptor fmm_solver;

#endif
