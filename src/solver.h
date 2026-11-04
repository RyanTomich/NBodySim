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
typedef void (*MetadataVisitor)(const char *name, MetadataValue value, void *context);

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
 */
bool solver_create(Solver *solver, const char *name, size_t count,
                   const SolverOption *options, size_t option_count, const char **error);
void solver_destroy(Solver *solver); /* Accepts an empty handle; clears it. */
/* Metadata requires a successfully prepared solver and a nonnull visitor. */
void solver_settings(const Solver *solver, MetadataVisitor visit, void *context);
void solver_diagnostics(const Solver *solver, MetadataVisitor visit, void *context);

#endif
