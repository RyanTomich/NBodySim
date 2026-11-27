#ifndef SOLVER_H
#define SOLVER_H

#include "nbody.h"

typedef struct SolverDescriptor SolverDescriptor;
typedef struct {
    const char *name;
    const char *value;
} SolverOption;

typedef enum { METADATA_SIZE, METADATA_DOUBLE, METADATA_BOOL, METADATA_STRING } MetadataType;
typedef struct {
    MetadataType type;
    union { size_t size; double real; bool boolean; const char *string; } value;
} MetadataValue;
/* Names and string values are borrowed for the duration of the visit. */
typedef void (*MetadataVisitor)(const char *name, MetadataValue value, void *visitor_context);

typedef struct {
    AccelerationSolver acceleration;
    void *context;
    const SolverDescriptor *descriptor;
} Solver;

size_t solver_count(void);
const char *solver_name_at(size_t index); /* NULL when out of range. */
/* Destination must not own a solver. Failure leaves it empty and sets error to
 * a static message. error and name must be nonnull; option names/values must
 * be nonnull strings. Options may be NULL when option_count is zero.
 * Duplicate option names are rejected before backend preparation.
 * All allocation occurs here. Prepared solvers retain no state pointers and
 * must only evaluate views with the specified count. Do not copy ownership.
 * A prepared solver must not be called concurrently. OpenMP n2 and grouped Barnes-Hut
 * capture the maximum thread count here; set OMP_NUM_THREADS before preparation.
 */
bool solver_create(Solver *solver, const char *name, size_t count,
                   const SolverOption *options, size_t option_count, const char **error);
void solver_destroy(Solver *solver); /* Accepts an empty handle; clears it. */
/* Metadata requires a successfully prepared solver and a nonnull visitor. */
void solver_settings(const Solver *solver, MetadataVisitor visit, void *visitor_context);
void solver_diagnostics(const Solver *solver, MetadataVisitor visit, void *visitor_context);

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
extern const SolverDescriptor serial_barnes_hut_solver;
extern const SolverDescriptor grouped_barnes_hut_solver;
#ifdef _OPENMP
extern const SolverDescriptor parallel_barnes_hut_solver;
#endif
extern const SolverDescriptor fmm_solver;

#endif
