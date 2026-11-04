#ifndef SOLVER_INTERNAL_H
#define SOLVER_INTERNAL_H
#include "solver.h"

/* create cleans up partial allocations on failure and supplies a static error.
 * Metadata callbacks emit scalar values only; they perform no serialization.
 * Only name and acceleration are required. Without create, options are rejected
 * and context is NULL. Other absent callbacks do nothing.
 * All lifecycle and metadata calls occur outside integration.
 */
struct SolverDescriptor {
    const char *name;
    AccelerationSolver acceleration;
    bool (*create)(void **context, size_t count, const SolverOption *options,
                   size_t option_count, const char **error);
    void (*destroy)(void *context);
    void (*settings)(const void *workspace, MetadataVisitor visit, void *context);
    void (*diagnostics)(const void *workspace, MetadataVisitor visit, void *context);
};
extern const SolverDescriptor direct_solver;
extern const SolverDescriptor approximate_solver;

#endif
