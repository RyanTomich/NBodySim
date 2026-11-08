#include "io.h"
#include "parse.h"

#include <ctype.h>

/* Bound token length and require whitespace between numbers. */
static bool read_token(FILE *input, char token[128])
{
    if (fscanf(input, "%127s", token) != 1) {
        return false;
    }
    int next = fgetc(input);
    return next == EOF ? !ferror(input) : isspace((unsigned char)next) != 0;
}

bool read_bodies(FILE *input, Simulation **simulation)
{
    char token[128];
    size_t count;
    *simulation = NULL;
    if (!read_token(input, token) || !parse_size(token, &count)) {
        return false;
    }

    Simulation *result;
    if (!simulation_create(&result, count)) {
        return false;
    }
    MutableSimulationView state = simulation_mutable_view(result);
    for (size_t i = 0; i < count; ++i) {
        double *fields[] = {
            &state.w[i], &state.x[i], &state.y[i], &state.z[i],
            &state.vx[i], &state.vy[i], &state.vz[i],
            &state.ax[i], &state.ay[i], &state.az[i]
        };
        for (size_t field = 0; field < sizeof(fields) / sizeof(fields[0]); ++field) {
            if (!read_token(input, token) || !parse_double(token, fields[field])) {
                simulation_destroy(result);
                return false;
            }
        }
        if (state.w[i] <= 0.0) {
            simulation_destroy(result);
            return false;
        }
    }

    int next;
    do {
        next = fgetc(input);
    } while (next != EOF && isspace((unsigned char)next));
    if (next != EOF || ferror(input)) {
        simulation_destroy(result);
        return false;
    }
    *simulation = result;
    return true;
}

bool write_bodies(FILE *output, const Simulation *simulation)
{
    SimulationView state = simulation_view(simulation);
    if (fprintf(output, "%zu\n", state.count) < 0) {
        return false;
    }
    for (size_t i = 0; i < state.count; ++i) {
        if (fprintf(output, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                    state.w[i], state.x[i], state.y[i], state.z[i],
                    state.vx[i], state.vy[i], state.vz[i],
                    state.ax[i], state.ay[i], state.az[i]) < 0) {
            return false;
        }
    }
    return true;
}
