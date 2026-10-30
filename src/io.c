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

bool read_bodies(FILE *input, Simulation *simulation)
{
    char token[128];
    size_t count;
    *simulation = (Simulation){0};
    if (!read_token(input, token) || !parse_size(token, &count)) {
        return false;
    }

    Simulation result;
    if (!simulation_create(&result, count)) {
        return false;
    }
    for (size_t i = 0; i < count; ++i) {
        double *fields[] = {
            &result.w[i], &result.x[i], &result.y[i], &result.z[i],
            &result.vx[i], &result.vy[i], &result.vz[i],
            &result.ax[i], &result.ay[i], &result.az[i]
        };
        for (size_t field = 0; field < sizeof(fields) / sizeof(fields[0]); ++field) {
            if (!read_token(input, token) || !parse_double(token, fields[field])) {
                free(result.storage);
                return false;
            }
        }
        if (result.w[i] <= 0.0) {
            free(result.storage);
            return false;
        }
    }

    int next;
    do {
        next = fgetc(input);
    } while (next != EOF && isspace((unsigned char)next));
    if (next != EOF || ferror(input)) {
        free(result.storage);
        return false;
    }
    *simulation = result;
    return true;
}

bool write_bodies(FILE *output, const Simulation *simulation)
{
    if (fprintf(output, "%zu\n", simulation->count) < 0) {
        return false;
    }
    for (size_t i = 0; i < simulation->count; ++i) {
        if (fprintf(output, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                    simulation->w[i], simulation->x[i], simulation->y[i], simulation->z[i],
                    simulation->vx[i], simulation->vy[i], simulation->vz[i],
                    simulation->ax[i], simulation->ay[i], simulation->az[i]) < 0) {
            return false;
        }
    }
    return true;
}
