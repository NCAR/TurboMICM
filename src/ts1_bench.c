#include "chem.h"
#include <errno.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern const ChemMechanism generated_mechanism;
typedef struct { ChemConditions conditions; double *y, *rates; } Cell;

/* ISO C line reader; unlike a fixed buffer it accepts long mechanism rows. */
static char *read_line(FILE *file)
{
    size_t n = 0, capacity = 256;
    char *line = malloc(capacity);
    if (!line) return NULL;
    int ch;
    while ((ch = fgetc(file)) != EOF && ch != '\n') {
        if (n+1 == capacity) {
            capacity *= 2;
            char *bigger = realloc(line, capacity);
            if (!bigger) { free(line); return NULL; }
            line = bigger;
        }
        if (ch != '\r') line[n++] = (char)ch;
    }
    if (ch == EOF && n == 0) { free(line); return NULL; }
    line[n] = '\0';
    return line;
}

static char **split(char *line, size_t *count)
{
    size_t n = 0, capacity = 16;
    char **fields = malloc(capacity*sizeof(*fields));
    if (!fields) return NULL;
    char *read = line, *write = line;
    for (;;) {
        if (n == capacity) {
            capacity *= 2;
            char **bigger = realloc(fields, capacity*sizeof(*fields));
            if (!bigger) { free(fields); return NULL; }
            fields = bigger;
        }
        fields[n++] = write;
        int quoted = *read == '"';
        if (quoted) ++read;
        while (*read) {
            if (quoted && *read == '"') {
                if (read[1] == '"') { *write++ = '"'; read += 2; continue; }
                ++read; break;
            }
            if (!quoted && *read == ',') break;
            *write++ = *read++;
        }
        int more = *read == ',';
        if (more) ++read;
        else if (*read) { free(fields); return NULL; }
        *write++ = '\0';
        if (!more) break;
    }
    *count = n;
    return fields;
}

static int read_cells(const char *path, Cell *cells, size_t count, const ChemMechanism *m)
{
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); return 1; }
    char *line = read_line(f);
    size_t cols = 0;
    char **header = line ? split(line, &cols) : NULL;
    int *mapping = header ? malloc(cols*sizeof(*mapping)) : NULL;
    if (!mapping) { fprintf(stderr, "invalid CSV header\n"); free(header); free(line); fclose(f); return 1; }
    for (size_t j = 0; j < cols; ++j) {
        mapping[j] = -3;
        if (!strcmp(header[j], "temperature")) mapping[j] = -1;
        else if (!strcmp(header[j], "pressure")) mapping[j] = -2;
        else {
            for (size_t k = 0; k < m->species_count; ++k)
                if (!strcmp(header[j], m->species_names[k])) mapping[j] = (int)k;
            for (size_t k = 0; k < m->parameter_count; ++k)
                if (!strcmp(header[j], m->parameter_names[k])) mapping[j] = (int)(m->species_count+k);
        }
        if (mapping[j] == -3) {
            fprintf(stderr, "unknown column: %s\n", header[j]);
            free(mapping); free(header); free(line); fclose(f); return 1;
        }
    }
    free(header); free(line);
    double *params = malloc(m->parameter_count*sizeof(*params));
    if (!params) { free(mapping); fclose(f); return 1; }
    int failed = 0;
    for (size_t i = 0; i < count; ++i) {
        cells[i].conditions = (ChemConditions){287.45, 101319.9, 0, 0};
        memcpy(cells[i].y, m->initial_species, m->species_count*sizeof(double));
        memcpy(params, m->initial_parameters, m->parameter_count*sizeof(double));
        line = read_line(f);
        size_t n = 0;
        char **fields = line ? split(line, &n) : NULL;
        if (!fields || n != cols) { fprintf(stderr, "invalid CSV row %zu\n", i+2); failed = 1; }
        if (!failed) for (size_t j = 0; j < cols; ++j) {
            char *end;
            errno = 0;
            double v = strtod(fields[j], &end);
            if (errno || end == fields[j] || *end || !isfinite(v)) { fprintf(stderr, "invalid number in row %zu\n", i+2); failed = 1; break; }
            int k = mapping[j];
            if (k == -1) cells[i].conditions.temperature = v;
            else if (k == -2) cells[i].conditions.pressure = v;
            else if ((size_t)k < m->species_count) cells[i].y[k] = v;
            else params[(size_t)k-m->species_count] = v;
        }
        if (!failed && chem_rate_constants(m, cells[i].conditions, params, cells[i].rates)) {
            fprintf(stderr, "invalid rate constants in row %zu\n", i+2); failed = 1;
        }
        free(fields); free(line);
        if (failed) break;
    }
    free(params); free(mapping); fclose(f);
    return failed;
}

static int write_cells(const Cell *cells, size_t count, const ChemMechanism *m)
{
    FILE *f = fopen("ts1_final.csv", "w");
    if (!f) { perror("ts1_final.csv"); return 1; }
    fputs("temperature,pressure", f);
    for (size_t i = 0; i < m->species_count; ++i) fprintf(f, ",%s", m->species_names[i]);
    fputc('\n', f);
    for (size_t cell = 0; cell < count; ++cell) {
        fprintf(f, "%.17g,%.17g", cells[cell].conditions.temperature, cells[cell].conditions.pressure);
        for (size_t i = 0; i < m->species_count; ++i) fprintf(f, ",%.17g", cells[cell].y[i]);
        fputc('\n', f);
    }
    return fclose(f) != 0;
}

int main(int argc, char **argv)
{
    int arg = 1, threads = omp_get_max_threads();
    if (argc > 2 && !strcmp(argv[1], "--threads")) { threads = atoi(argv[2]); arg += 2; }
    if (argc != arg+4 || threads < 1) { fprintf(stderr, "usage: %s [--threads N] CELLS STEPS DT CELLS_CSV\n", argv[0]); return 1; }
    char *end;
    unsigned long count = strtoul(argv[arg], &end, 10); if (*end || !count) return 1;
    unsigned long steps = strtoul(argv[arg+1], &end, 10); if (*end) return 1;
    double dt = strtod(argv[arg+2], &end); if (*end || !(dt > 0) || !isfinite(dt)) return 1;
    const ChemMechanism *m = &generated_mechanism;
    Cell *cells = calloc(count, sizeof(*cells));
    double *ys = calloc(count*m->species_count, sizeof(*ys));
    double *rates = calloc(count*m->reaction_count, sizeof(*rates));
    if (!cells || !ys || !rates) { fprintf(stderr, "out of memory\n"); free(cells); free(ys); free(rates); return 1; }
    for (size_t i = 0; i < count; ++i) { cells[i].y = ys+i*m->species_count; cells[i].rates = rates+i*m->reaction_count; }
    if (read_cells(argv[arg+3], cells, count, m)) { free(cells); free(ys); free(rates); return 1; }
    ChemRosenbrockParameters p = chem_rosenbrock_three_stage();
    unsigned long internal=0, accepted=0, rejected=0, functions=0, jacobians=0, decompositions=0, solves=0;
    int failed = 0;
    omp_set_num_threads(threads);
    double start = 0, finish = 0;
#pragma omp parallel reduction(+:internal,accepted,rejected,functions,jacobians,decompositions,solves) reduction(|:failed)
    {
        ChemWorkspace *work = chem_workspace_create(m);
        if (!work) failed = 1;
#pragma omp for schedule(static)
        for (size_t i = 0; i < count; ++i) {
            ChemStatistics s;
            if (!work || chem_rosenbrock_solve(m, work, &p, cells[i].conditions, cells[i].rates, dt, cells[i].y, &s)) failed = 1;
        }
#pragma omp barrier
#pragma omp single
        start = omp_get_wtime();
#pragma omp for schedule(static)
        for (size_t i = 0; i < count; ++i) {
            for (size_t j = 0; j < steps && work; ++j) {
                ChemStatistics s;
                if (chem_rosenbrock_solve(m, work, &p, cells[i].conditions, cells[i].rates, dt, cells[i].y, &s)) { failed = 1; break; }
                internal += s.internal_steps; accepted += s.accepted; rejected += s.rejected;
                functions += s.function_calls; jacobians += s.jacobian_updates; decompositions += s.decompositions; solves += s.solves;
            }
        }
#pragma omp single
        finish = omp_get_wtime();
        chem_workspace_free(work);
    }
    int status = failed || write_cells(cells, count, m);
    if (failed) fprintf(stderr, "solver failed\n");
    if (!status) printf("elapsed_ms=%.3f\ninternal_steps=%lu accepted=%lu rejected=%lu function_calls=%lu jacobian_updates=%lu decompositions=%lu solves=%lu\n",
                        (finish-start)*1000, internal, accepted, rejected, functions, jacobians, decompositions, solves);
    free(cells); free(ys); free(rates);
    return status;
}
