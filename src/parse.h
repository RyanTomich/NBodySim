#ifndef PARSE_H
#define PARSE_H
#include <stdbool.h>
#include <stddef.h>
/* Decimal nonnegative integer or finite floating-point number. */
bool parse_size(const char *text, size_t *value);
bool parse_double(const char *text, double *value);
#endif
