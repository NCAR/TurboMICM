#include "chem.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern const ChemMechanism chem_generated_mechanism;

typedef struct {
    size_t first, count;
    ChemState *state;
    ChemStats stats;
    int status;
} Group;

typedef struct { int kind; size_t index; } Column;
enum { COLUMN_TEMPERATURE, COLUMN_PRESSURE, COLUMN_SPECIES, COLUMN_PARAMETER };

static int parse_size(const char *text, size_t *value, int allow_zero)
{
    char *end;
    if (!text || !*text || *text == '-') return 0;
    errno = 0;
    unsigned long long n = strtoull(text, &end, 10);
    if (errno || *end || n > SIZE_MAX || (!allow_zero && n == 0)) return 0;
    *value = (size_t)n;
    return 1;
}

/* ISO C has no getline; grow the input buffer without imposing a CSV row limit. */
static int read_line(FILE *file, char **line, size_t *capacity)
{
    size_t n = 0;
    int ch;
    while ((ch = fgetc(file)) != EOF && ch != '\n') {
        if (n+1 >= *capacity) {
            if (*capacity > SIZE_MAX/2) return -1;
            size_t next = *capacity ? *capacity*2 : 256;
            char *buffer = realloc(*line, next);
            if (!buffer) return -1;
            *line = buffer; *capacity = next;
        }
        (*line)[n++] = (char)ch;
    }
    if (ferror(file)) return -1;
    if (!n && ch == EOF) return 0;
    if (n && (*line)[n-1] == '\r') --n;
    if (!*line) {
        *line = malloc(1);
        if (!*line) return -1;
        *capacity = 1;
    }
    (*line)[n] = '\0';
    return 1;
}

static int read_cells(const char *path, Group *groups, size_t group_size, size_t cells,
                      const ChemMechanism *m)
{
    FILE *file = fopen(path, "r");
    if (!file) { fprintf(stderr, "cannot open %s\n", path); return 0; }
    char *line = NULL;
    size_t capacity = 0, count = 0;
    Column *columns = NULL;
    int ok = 0;
    if (read_line(file, &line, &capacity) != 1 || !*line) {
        fprintf(stderr, "%s: missing header\n", path); goto done;
    }
    size_t ncols = 1;
    for (const char *p = line; *p; ++p) ncols += *p == ',';
    columns = calloc(ncols, sizeof(*columns));
    if (!columns) goto done;
    char *field = line;
    for (size_t i = 0; i < ncols; ++i) {
        char *comma = strchr(field, ',');
        if (comma) *comma = '\0';
        int index;
        if (!strcmp(field, "temperature")) columns[i] = (Column){COLUMN_TEMPERATURE, 0};
        else if (!strcmp(field, "pressure")) columns[i] = (Column){COLUMN_PRESSURE, 0};
        else if ((index = chem_find_species(m, field)) >= 0) columns[i] = (Column){COLUMN_SPECIES, (size_t)index};
        else if ((index = chem_find_parameter(m, field)) >= 0) columns[i] = (Column){COLUMN_PARAMETER, (size_t)index};
        else { fprintf(stderr, "unknown CSV column '%s'\n", field); goto done; }
        for (size_t j = 0; j < i; ++j) {
            if (columns[i].kind == columns[j].kind && columns[i].index == columns[j].index) {
                fprintf(stderr, "duplicate CSV column '%s'\n", field); goto done;
            }
        }
        field = comma ? comma+1 : field+strlen(field);
    }
    while (count < cells) {
        int read = read_line(file, &line, &capacity);
        if (read <= 0) { fprintf(stderr, "%s: expected %zu cells, found %zu\n", path, cells, count); goto done; }
        if (!*line) continue;
        ChemState *state = groups[count/group_size].state;
        size_t cell = count % group_size;
        field = line;
        for (size_t i = 0; i < ncols; ++i) {
            char *comma = strchr(field, ',');
            if (comma) *comma = '\0';
            if ((i+1 < ncols) != (comma != NULL)) {
                fprintf(stderr, "%s: wrong field count on row %zu\n", path, count+2); goto done;
            }
            char *end;
            errno = 0;
            double value = strtod(field, &end);
            if (end == field || *end || errno == ERANGE || !isfinite(value)) {
                fprintf(stderr, "%s: invalid value on row %zu, column %zu\n", path, count+2, i+1); goto done;
            }
            Column col = columns[i];
            double *out = col.kind == COLUMN_TEMPERATURE ? chem_temperatures(state) :
                          col.kind == COLUMN_PRESSURE ? chem_pressures(state) :
                          col.kind == COLUMN_SPECIES ? chem_species_values(state, col.index) :
                          chem_parameter_values(state, col.index);
            out[cell] = value;
            field = comma ? comma+1 : end;
        }
        ++count;
    }
    ok = 1;
done:
    free(columns); free(line); fclose(file);
    return ok;
}

static int write_final(const Group *groups, size_t ngroups, const ChemMechanism *m)
{
    FILE *file = fopen("ts1_final.csv", "w");
    if (!file) { fprintf(stderr, "cannot write ts1_final.csv\n"); return 0; }
    fprintf(file, "temperature,pressure");
    for (size_t i = 0; i < m->species_count; ++i) fprintf(file, ",%s", m->species[i].name);
    fputc('\n', file);
    for (size_t g = 0; g < ngroups; ++g) {
        ChemState *s = groups[g].state;
        for (size_t c = 0; c < groups[g].count; ++c) {
            fprintf(file, "%.17g,%.17g", chem_temperatures(s)[c], chem_pressures(s)[c]);
            for (size_t i = 0; i < m->species_count; ++i) fprintf(file, ",%.17g", chem_species_values(s, i)[c]);
            fputc('\n', file);
        }
    }
    int success = !ferror(file);
    if (fclose(file)) success = 0;
    return success;
}

static void add_stats(ChemStats *sum, const ChemStats *x)
{
    sum->internal_steps += x->internal_steps; sum->accepted += x->accepted;
    sum->rejected += x->rejected; sum->function_calls += x->function_calls;
    sum->jacobian_updates += x->jacobian_updates; sum->decompositions += x->decompositions;
    sum->solves += x->solves;
}

int main(int argc, char **argv)
{
    const ChemMechanism *m = &chem_generated_mechanism;
    size_t threads = (size_t)omp_get_max_threads(), group_size = 64, cells, steps;
    const char *setting = getenv("CHEM_GROUP_SIZE");
    if (setting && !parse_size(setting, &group_size, 0)) {
        fprintf(stderr, "invalid CHEM_GROUP_SIZE\n"); return 2;
    }
    int arg = 1;
    while (arg < argc && !strncmp(argv[arg], "--", 2)) {
        if (arg+1 >= argc) goto usage;
        if (!strcmp(argv[arg], "--threads")) {
            if (!parse_size(argv[arg+1], &threads, 0) || threads > INT_MAX) goto usage;
        } else if (!strcmp(argv[arg], "--group-size")) {
            if (!parse_size(argv[arg+1], &group_size, 0)) goto usage;
        } else goto usage;
        arg += 2;
    }
    if (argc-arg != 4 || !parse_size(argv[arg], &cells, 0) || !parse_size(argv[arg+1], &steps, 1)) goto usage;
    char *end;
    errno = 0;
    double dt = strtod(argv[arg+2], &end);
    if (end == argv[arg+2] || *end || errno || !isfinite(dt) || dt <= 0) goto usage;
    size_t ngroups = (cells-1)/group_size+1;
    if (ngroups > SIZE_MAX/sizeof(Group)) return 1;
    Group *groups = calloc(ngroups, sizeof(*groups));
    if (!groups) return 1;
    int result = 1;
    for (size_t g = 0; g < ngroups; ++g) {
        groups[g].first = g*group_size;
        groups[g].count = cells-groups[g].first < group_size ? cells-groups[g].first : group_size;
        groups[g].state = g ? chem_state_create_like(groups[0].state, groups[g].count) : chem_state_create(m, groups[g].count);
        if (!groups[g].state) { fprintf(stderr, "could not allocate solver state\n"); goto done; }
    }
    if (!read_cells(argv[arg+3], groups, group_size, cells, m)) goto done;
    omp_set_dynamic(0);
    omp_set_num_threads((int)threads);
    if (!getenv("OMP_SCHEDULE")) omp_set_schedule(omp_sched_static,0);
    ChemRosenbrock p = chem_rosenbrock_three_stage();
    #pragma omp parallel for schedule(runtime)
    for (size_t g = 0; g < ngroups; ++g) {
        ChemStats warmup;
        groups[g].status = chem_update_rates(groups[g].state);
        if (groups[g].status == CHEM_OK) groups[g].status = chem_solve(groups[g].state, &p, dt, &warmup);
    }
    for (size_t g = 0; g < ngroups; ++g) if (groups[g].status) {
        fprintf(stderr, "warm-up group %zu: %s\n", g, chem_error_string(groups[g].status)); goto done;
    }
    double start = omp_get_wtime();
    #pragma omp parallel for schedule(runtime)
    for (size_t g = 0; g < ngroups; ++g) {
        for (size_t step = 0; step < steps && !groups[g].status; ++step) {
            ChemStats stats;
            groups[g].status = chem_solve(groups[g].state, &p, dt, &stats);
            add_stats(&groups[g].stats, &stats);
        }
    }
    double elapsed_ms = (omp_get_wtime()-start)*1000.0;
    ChemStats total = {0};
    for (size_t g = 0; g < ngroups; ++g) {
        if (groups[g].status) { fprintf(stderr, "group %zu: %s\n", g, chem_error_string(groups[g].status)); goto done; }
        add_stats(&total, &groups[g].stats);
    }
    if (!write_final(groups, ngroups, m)) goto done;
    printf("elapsed_ms=%.3f\ninternal_steps=%lu accepted=%lu rejected=%lu function_calls=%lu "
           "jacobian_updates=%lu decompositions=%lu solves=%lu\n", elapsed_ms,
           total.internal_steps, total.accepted, total.rejected, total.function_calls,
           total.jacobian_updates, total.decompositions, total.solves);
    result = 0;
done:
    for (size_t g = 0; g < ngroups; ++g) chem_state_free(groups[g].state);
    free(groups);
    return result;
usage:
    fprintf(stderr, "usage: %s [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV\n", argv[0]);
    return 2;
}
