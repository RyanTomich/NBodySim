#include "solver_internal.h"
#include <string.h>

static const SolverDescriptor *const registry[] = {&direct_solver, &approximate_solver};

size_t solver_count(void)
{
    return sizeof(registry) / sizeof(registry[0]);
}

const char *solver_name_at(size_t index)
{
    return index < solver_count() ? registry[index]->name : NULL;
}

bool solver_create(Solver *solver, const char *name, size_t count,
                   const SolverOption *options, size_t option_count, const char **error)
{
    *solver = (Solver){0};
    *error = NULL;
    for (size_t i = 0; i < option_count; ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (strcmp(options[i].name, options[j].name) == 0) {
                *error = "Duplicate solver option.";
                return false;
            }
        }
    }
    for (size_t i = 0; i < solver_count(); ++i) {
        const SolverDescriptor *descriptor = registry[i];
        if (strcmp(name, descriptor->name) != 0) {
            continue;
        }
        void *context = NULL;
        if (descriptor->create) {
            if (!descriptor->create(&context, count, options, option_count, error)) {
                return false;
            }
        } else if (option_count) {
            *error = "Solver does not accept options.";
            return false;
        }
        *solver = (Solver){descriptor->acceleration, context, descriptor};
        return true;
    }
    *error = "Unknown solver.";
    return false;
}

void solver_destroy(Solver *solver)
{
    if (solver->descriptor && solver->descriptor->destroy) {
        solver->descriptor->destroy(solver->context);
    }
    *solver = (Solver){0};
}

void solver_settings(const Solver *solver, MetadataVisitor visit, void *context)
{
    if (solver->descriptor->settings) {
        solver->descriptor->settings(solver->context, visit, context);
    }
}

void solver_diagnostics(const Solver *solver, MetadataVisitor visit, void *context)
{
    if (solver->descriptor->diagnostics) {
        solver->descriptor->diagnostics(solver->context, visit, context);
    }
}
