#include "io.h"

#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

bool parse_size(const char *text, size_t *value)
{
    size_t result = 0;
    if (*text == '\0') {
        return false;
    }
    for (; *text; ++text) {
        if (*text < '0' || *text > '9') {
            return false;
        }
        unsigned digit = (unsigned)(*text - '0');
        if (result > (SIZE_MAX - digit) / 10) {
            return false;
        }
        result = result * 10 + digit;
    }
    *value = result;
    return true;
}

bool parse_double(const char *text, double *value)
{
    char *end;
    double result = strtod(text, &end);
    if (end == text || *end != '\0' || !isfinite(result)) {
        return false;
    }
    *value = result;
    return true;
}

/* Bound token length and require whitespace between numbers. */
static bool read_token(FILE *input, char token[128])
{
    if (fscanf(input, "%127s", token) != 1) {
        return false;
    }
    int next = fgetc(input);
    return next == EOF ? !ferror(input) : isspace((unsigned char)next) != 0;
}

bool read_bodies(FILE *input, Body **bodies, size_t *count)
{
    char token[128];
    size_t n;
    *bodies = NULL;
    *count = 0;
    if (!read_token(input, token) || !parse_size(token, &n)
        || n > SIZE_MAX / sizeof(Body)) {
        return false;
    }

    Body *result = n == 0 ? NULL : malloc(n * sizeof(*result));
    if (n != 0 && !result) {
        return false;
    }
    for (size_t i = 0; i < n; ++i) {
        Body *body = &result[i];
        double *fields[] = {
            &body->w,
            &body->position.x, &body->position.y, &body->position.z,
            &body->velocity.x, &body->velocity.y, &body->velocity.z,
            &body->acceleration.x, &body->acceleration.y, &body->acceleration.z
        };
        for (size_t field = 0; field < sizeof(fields) / sizeof(fields[0]); ++field) {
            if (!read_token(input, token) || !parse_double(token, fields[field])) {
                free(result);
                return false;
            }
        }
        if (body->w <= 0.0) {
            free(result);
            return false;
        }
    }

    int next;
    do {
        next = fgetc(input);
    } while (next != EOF && isspace((unsigned char)next));
    if (next != EOF || ferror(input)) {
        free(result);
        return false;
    }
    *bodies = result;
    *count = n;
    return true;
}

bool write_bodies(FILE *output, const Body *bodies, size_t count)
{
    if (fprintf(output, "%zu\n", count) < 0) {
        return false;
    }
    for (size_t i = 0; i < count; ++i) {
        const Body *body = &bodies[i];
        if (fprintf(output, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                    body->w,
                    body->position.x, body->position.y, body->position.z,
                    body->velocity.x, body->velocity.y, body->velocity.z,
                    body->acceleration.x, body->acceleration.y, body->acceleration.z) < 0) {
            return false;
        }
    }
    return true;
}
