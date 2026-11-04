#ifndef BENCH_METADATA_H
#define BENCH_METADATA_H
#include "solver.h"
#include <stdio.h>
void metadata_log(FILE *output, const Solver *solver, bool diagnostics);
int metadata_query(int argc, char **argv);
#endif
