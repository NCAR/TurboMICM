#include "chem.h"
#include "rates.h"
#include "rosenbrock.h"
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern const ChemMechanism generated_mechanism;
extern const char *const generated_species[], *const generated_custom_labels[];
extern const size_t generated_custom_count;
extern const ChemRate generated_rates[];
extern const ChemRateMultiplier generated_multipliers[];
extern const ChemParameterizedSpecies generated_multiplier_functions[];
extern const double generated_custom_defaults[], generated_initial_concentrations[];
extern const double generated_default_temperature, generated_default_pressure;

typedef struct {
    double temperature, pressure;
    double *y, *custom, *rates;
    const ChemSparsePlan *plan;
    ChemSolverStats stats;
    ChemSolveStatus status;
} Cell;

static char *read_line(FILE *file)
{
    size_t capacity = 256, length = 0;
    char *line = malloc(capacity);
    if (!line) return NULL;
    int c;
    while ((c = fgetc(file)) != EOF && c != '\n') {
        if (length+1 >= capacity) {
            if (capacity > SIZE_MAX/2) { free(line); return NULL; }
            capacity *= 2;
            char *tmp = realloc(line, capacity);
            if (!tmp) { free(line); return NULL; }
            line = tmp;
        }
        if (c != '\r') line[length++] = (char)c;
    }
    if (c == EOF && !length) { free(line); return NULL; }
    line[length] = '\0';
    return line;
}

static int find_name(const char *name, const char *const *names, size_t count)
{
    for (size_t i = 0; i < count; ++i) if (!strcmp(name, names[i])) return (int)i;
    return -1;
}

/* Header field codes: -1 temperature, -2 pressure, then species and custom. */
static int read_cells(const char *path, Cell *cells, size_t count)
{
    FILE *file = fopen(path, "r");
    if (!file) { fprintf(stderr, "cannot open %s\n", path); return 0; }
    char *line = read_line(file);
    if (!line) { fprintf(stderr, "missing CSV header\n"); fclose(file); return 0; }
    size_t columns = 1;
    for (char *p = line; *p; ++p) if (*p == ',') ++columns;
    int *map = malloc(columns*sizeof(*map));
    if (!map) { free(line); fclose(file); return 0; }
    char *field = line;
    int ok = 1;
    for (size_t i = 0; i < columns; ++i) {
        char *comma = strchr(field, ','); if (comma) *comma = 0;
        int index = find_name(field, generated_species, generated_mechanism.species_count);
        if (!strcmp(field,"temperature")) map[i] = -1;
        else if (!strcmp(field,"pressure")) map[i] = -2;
        else if (index >= 0) map[i] = index;
        else if ((index = find_name(field, generated_custom_labels, generated_custom_count)) >= 0)
            map[i] = (int)generated_mechanism.species_count + index;
        else { fprintf(stderr, "unknown CSV column: %s\n", field); ok = 0; break; }
        for (size_t j = 0; j < i; ++j) if (map[j] == map[i]) {
            fprintf(stderr, "duplicate CSV column: %s\n", field); ok = 0;
        }
        field = comma ? comma+1 : field+strlen(field);
    }
    free(line);
    for (size_t row = 0; ok && row < count; ++row) {
        do {
            line = read_line(file);
            if (line && !*line) { free(line); line = NULL; }
            else break;
        } while (!feof(file));
        if (!line) { fprintf(stderr, "not enough CSV rows\n"); ok = 0; break; }
        field = line;
        for (size_t col = 0; col < columns; ++col) {
            char *comma = strchr(field, ','); if (comma) *comma = 0;
            char *end; errno = 0; double value = strtod(field, &end);
            if (errno || end == field || *end || !isfinite(value) || (!!comma != (col+1 < columns))) {
                fprintf(stderr, "invalid CSV field at row %zu column %zu\n", row+1, col+1); ok = 0; break;
            }
            int key = map[col];
            if (key == -1) cells[row].temperature = value;
            else if (key == -2) cells[row].pressure = value;
            else if ((size_t)key < generated_mechanism.species_count) cells[row].y[key] = value;
            else cells[row].custom[(size_t)key-generated_mechanism.species_count] = value;
            if (comma) field = comma+1;
        }
        free(line);
    }
    free(map); fclose(file); return ok;
}

static void rhs(void *ctx, const double *y, double *out)
{ chem_rhs(&generated_mechanism, y, ((Cell *)ctx)->rates, out); }
static void jac(void *ctx, const double *y, double *out)
{
    Cell *cell=ctx;
    chem_sparse_clear(cell->plan,out);
    chem_accumulate_jacobian(&generated_mechanism,y,cell->rates,out);
}
static void add_stats(ChemSolverStats *a, const ChemSolverStats *b)
{
    a->steps += b->steps; a->accepted += b->accepted; a->rejected += b->rejected;
    a->function_calls += b->function_calls; a->jacobian_updates += b->jacobian_updates;
    a->decompositions += b->decompositions; a->solves += b->solves;
}
static int parse_count(const char *s, size_t *value)
{
    char *end; errno = 0;
    if (*s == '-' || !*s) return 0;
    unsigned long long n = strtoull(s, &end, 10);
    if (errno || *end || n > SIZE_MAX) return 0;
    *value = (size_t)n; return 1;
}
int main(int argc, char **argv)
{
    size_t threads = (size_t)omp_get_max_threads(), count, steps;
    int arg = 1;
    if (arg < argc && !strcmp(argv[arg], "--threads")) {
        if (arg+1 >= argc || !parse_count(argv[arg+1], &threads)) threads = 0;
        arg += 2;
    }
    char *end;
    if (argc-arg != 4 || !threads || threads > INT_MAX ||
        !parse_count(argv[arg], &count) || !count || !parse_count(argv[arg+1], &steps)) {
        fprintf(stderr, "usage: ts1_bench_cpu [--threads N] CELLS STEPS DT CELLS_CSV\n"); return 1;
    }
    errno = 0; double dt = strtod(argv[arg+2], &end);
    if (errno || *end || end == argv[arg+2] || !isfinite(dt) || dt <= 0) { fprintf(stderr, "invalid DT\n"); return 1; }
    const size_t n = generated_mechanism.species_count, nr = generated_mechanism.reaction_count;
    if (count > SIZE_MAX/sizeof(Cell) || n > SIZE_MAX-generated_custom_count || n+generated_custom_count > SIZE_MAX-nr ||
        count > SIZE_MAX/(n+generated_custom_count+nr)/sizeof(double)) return 1;
    Cell *cells = calloc(count, sizeof(*cells));
    double *storage = calloc(count*(n+generated_custom_count+nr), sizeof(*storage));
    double *atol = malloc(n*sizeof(*atol));
    if (!cells || !storage || !atol) { fprintf(stderr, "allocation failed\n"); free(cells); free(storage); free(atol); return 1; }
    for (size_t i = 0; i < n; ++i) atol[i] = 1e-3;
    for (size_t i = 0; i < count; ++i) {
        cells[i].y = storage+i*(n+generated_custom_count+nr);
        cells[i].custom = cells[i].y+n; cells[i].rates = cells[i].custom+generated_custom_count;
        memcpy(cells[i].y, generated_initial_concentrations, n*sizeof(double));
        memcpy(cells[i].custom, generated_custom_defaults, generated_custom_count*sizeof(double));
        cells[i].temperature = generated_default_temperature; cells[i].pressure = generated_default_pressure;
    }
    if (!read_cells(argv[arg+3], cells, count)) { free(cells); free(storage); free(atol); return 1; }
    ChemRosenbrockMethod method = chem_rosenbrock_three_stage();
    ChemRosenbrockOptions options = chem_rosenbrock_defaults(atol);
    ChemSparsePlan *plan = chem_sparse_plan_create(&generated_mechanism);
    if (!plan) { fprintf(stderr, "cannot create sparse plan\n"); free(cells); free(storage); free(atol); return 1; }
    for (size_t i=0;i<count;++i) cells[i].plan=plan;
    options.sparse_plan = plan;
    omp_set_num_threads((int)threads);
    /* Allow the standard OMP_SCHEDULE setting to choose policy and chunk
       size; absent that setting, preserve the benchmark's static default. */
    if (!getenv("OMP_SCHEDULE")) omp_set_schedule(omp_sched_static,0);
#pragma omp parallel for schedule(runtime)
    for (size_t i = 0; i < count; ++i) {
        ChemConditions c = {cells[i].temperature, cells[i].pressure,
            cells[i].pressure/(1.380649e-23*6.02214076e23*cells[i].temperature)};
        chem_rates_parameterized(nr, generated_rates, generated_multipliers, generated_multiplier_functions, c, cells[i].custom, cells[i].rates);
        cells[i].status = chem_rosenbrock_solve(n,cells[i].y,dt,rhs,jac,&cells[i],&method,&options,NULL);
    }
    double start = omp_get_wtime();
#pragma omp parallel for schedule(runtime)
    for (size_t i = 0; i < count; ++i) {
        for (size_t step = 0; step < steps && cells[i].status == CHEM_SOLVE_OK; ++step) {
            ChemSolverStats stats;
            cells[i].status = chem_rosenbrock_solve(n,cells[i].y,dt,rhs,jac,&cells[i],&method,&options,&stats);
            add_stats(&cells[i].stats, &stats);
        }
    }
    double elapsed = (omp_get_wtime()-start)*1000;
    ChemSolverStats total = {0}; int ok = 1;
    for (size_t i = 0; i < count; ++i) {
        if (cells[i].status != CHEM_SOLVE_OK) { fprintf(stderr, "cell %zu failed with status %d\n", i, cells[i].status); ok = 0; break; }
        add_stats(&total, &cells[i].stats);
    }
    if (ok) {
        FILE *out = fopen("ts1_final.csv", "w");
        if (!out) { perror("ts1_final.csv"); ok = 0; }
        else {
            fprintf(out, "temperature,pressure");
            for (size_t j = 0; j < n; ++j) fprintf(out, ",%s", generated_species[j]);
            fputc('\n', out);
            for (size_t i = 0; i < count; ++i) {
                fprintf(out, "%.17g,%.17g", cells[i].temperature,cells[i].pressure);
                for (size_t j = 0; j < n; ++j) fprintf(out, ",%.17g", cells[i].y[j]);
                fputc('\n',out);
            }
            if (fclose(out)) { perror("ts1_final.csv"); ok = 0; }
        }
    }
    if (ok) printf("elapsed_ms=%.6f\ninternal_steps=%zu accepted=%zu rejected=%zu function_calls=%zu jacobian_updates=%zu decompositions=%zu solves=%zu\n",
        elapsed,total.steps,total.accepted,total.rejected,total.function_calls,total.jacobian_updates,total.decompositions,total.solves);
    chem_sparse_plan_destroy(plan);
    free(cells); free(storage); free(atol); return ok ? 0 : 1;
}
