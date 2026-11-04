#include "metadata.h"
#include "parse.h"
#include <stdlib.h>

typedef struct {
    FILE *output;
    size_t fields;
    const char *prefix; /* Emit a log prefix only if a field is reported. */
} JsonWriter;

static void json_string(FILE *output, const char *text)
{
    fputc('"', output);
    for (const unsigned char *p = (const unsigned char *)text; *p; ++p) {
        if (*p == '"' || *p == '\\') {
            fputc('\\', output);
        }
        if (*p < 32) {
            fprintf(output, "\\u%04x", *p);
        } else {
            fputc(*p, output);
        }
    }
    fputc('"', output);
}

static void json_field(const char *name, MetadataValue value, void *context)
{
    JsonWriter *writer = context;
    FILE *output = writer->output;
    if (!writer->fields && writer->prefix) {
        fputs(writer->prefix, output);
    }
    if (writer->fields++) {
        fputs(", ", output);
    }
    json_string(output, name);
    fputs(": ", output);
    switch (value.type) {
    case METADATA_SIZE: fprintf(output, "%zu", value.value.size); break;
    case METADATA_DOUBLE: fprintf(output, "%.17g", value.value.real); break;
    case METADATA_BOOL: fputs(value.value.boolean ? "true" : "false", output); break;
    case METADATA_STRING: json_string(output, value.value.string); break;
    }
}

void metadata_log(FILE *output, const Solver *solver, bool diagnostics)
{
    JsonWriter writer = {output, 0, diagnostics ? "solver_totals: {" : "solver_config: {"};
    if (diagnostics) {
        solver_diagnostics(solver, json_field, &writer);
    } else {
        solver_settings(solver, json_field, &writer);
    }
    if (writer.fields) {
        if (diagnostics) {
            json_field("scope", (MetadataValue){METADATA_STRING, {.string = "warmups_and_samples"}}, &writer);
        }
        fputs("}\n", output);
    }
}

int metadata_query(int argc, char **argv)
{
    if (argc == 0) {
        fputc('[', stdout);
        for (size_t i = 0; i < solver_count(); ++i) {
            if (i) {
                fputc(',', stdout);
            }
            json_string(stdout, solver_name_at(i));
        }
        puts("]");
    } else {
        size_t count;
        if (argc < 2 || argc > 3 || !parse_size(argv[1], &count)) {
            return EXIT_FAILURE;
        }
        Solver solver;
        const char *error;
        SolverOption option = {"theta", argc == 3 ? argv[2] : NULL};
        if (!solver_create(&solver, argv[0], count, &option, argc == 3 ? 1 : 0, &error)) {
            fprintf(stderr, "%s\n", error);
            return EXIT_FAILURE;
        }
        JsonWriter writer = {.output = stdout};
        fputc('{', stdout);
        solver_settings(&solver, json_field, &writer);
        puts("}");
        solver_destroy(&solver);
    }
    return fflush(stdout) == 0 && !ferror(stdout) ? EXIT_SUCCESS : EXIT_FAILURE;
}
