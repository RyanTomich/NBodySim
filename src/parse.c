#include "parse.h"
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

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
