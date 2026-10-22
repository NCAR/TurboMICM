/* Portable TS1 scenario driver; all filesystem access uses ISO C streams. */
#include "chem_specialized.h"
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef CHEM_GROUP_SIZE
#define CHEM_GROUP_SIZE 32
#endif
#include "bench_io.h"

typedef struct {
    char *cells_path, *photo_path;
    size_t steps, snapshot_every;
    double dt;
} Scenario;
typedef struct { size_t count; double *times, *factors; } Photolysis;

static char *trim(char *s)
{
    while (isspace((unsigned char)*s)) ++s;
    char *end = s + strlen(s);
    while (end > s && isspace((unsigned char)end[-1])) --end;
    *end = 0;
    return s;
}
static char *relative_path(const char *scenario, const char *value)
{
    const char *slash = strrchr(scenario, '/');
    size_t prefix = value[0] == '/' || !slash ? 0 : (size_t)(slash - scenario + 1);
    size_t length = strlen(value);
    if (prefix > SIZE_MAX - length - 1) return NULL;
    char *path = malloc(prefix + length + 1);
    if (path) { memcpy(path, scenario, prefix); memcpy(path + prefix, value, length + 1); }
    return path;
}
static int read_scenario(const char *path, Scenario *s)
{
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    int ok = 1;
    unsigned seen = 0;
    char *line;
    while ((line = read_line(f))) {
        char *comment = strchr(line, '#');
        if (comment) *comment = 0;
        char *key = trim(line), *equal = strchr(key, '=');
        if (!*key) { free(line); continue; }
        if (!equal) { free(line); ok = 0; break; }
        *equal = 0;
        char *value = trim(equal + 1);
        key = trim(key);
        unsigned field = 0;
        if (!strcmp(key, "cells")) { field = 1; s->cells_path = relative_path(path, value); ok &= s->cells_path != NULL && *value; }
        else if (!strcmp(key, "photolysis")) { field = 2; s->photo_path = relative_path(path, value); ok &= s->photo_path != NULL && *value; }
        else if (!strcmp(key, "dt")) { field = 4; ok &= parse_double(value, &s->dt) && s->dt > 0; }
        else if (!strcmp(key, "steps")) { field = 8; ok &= parse_size(value, &s->steps) && s->steps > 0; }
        else if (!strcmp(key, "snapshot_every")) { field = 16; ok &= parse_size(value, &s->snapshot_every); }
        else ok = 0;
        if (seen & field) ok = 0;
        seen |= field;
        free(line);
        if (!ok) break;
    }
    ok &= !ferror(f) && (seen & 15) == 15;
    fclose(f);
    return ok;
}
static size_t count_cells(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    char *line = read_line(f);
    free(line);
    size_t count = 0;
    while ((line = read_line(f))) { if (*line) ++count; free(line); }
    if (ferror(f)) count = 0;
    fclose(f);
    return count;
}
static int read_photolysis(const char *path, Photolysis *p)
{
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    char *line = read_line(f);
    int ok = line && !strcmp(line, "time_s,factor");
    free(line);
    while (ok && (line = read_line(f))) {
        if (!*line) { free(line); continue; }
        char **fields = NULL;
        size_t count = split(line, &fields);
        double t = 0, factor = 0;
        ok = count == 2 && parse_double(fields[0], &t) && parse_double(fields[1], &factor)
            && factor >= 0 && (!p->count || t > p->times[p->count - 1]);
        if (ok) {
            double *times = realloc(p->times, (p->count + 1) * sizeof(double));
            if (!times) ok = 0;
            else p->times = times;
            double *factors = realloc(p->factors, (p->count + 1) * sizeof(double));
            if (!factors) ok = 0;
            else p->factors = factors;
            if (ok) { p->times[p->count] = t; p->factors[p->count++] = factor; }
        }
        free(fields); free(line);
    }
    ok &= !ferror(f) && p->count > 0;
    fclose(f);
    return ok;
}
static double photo_factor(const Photolysis *p, double t)
{
    if (t <= p->times[0]) return p->factors[0];
    if (t >= p->times[p->count - 1]) return p->factors[p->count - 1];
    size_t lo = 0, hi = p->count - 1;
    while (hi - lo > 1) {
        size_t mid = lo + (hi - lo) / 2;
        if (p->times[mid] <= t) lo = mid; else hi = mid;
    }
    double w = (t - p->times[lo]) / (p->times[hi] - p->times[lo]);
    return (1.0 - w) * p->factors[lo] + w * p->factors[hi];
}
static int write_snapshots(const Group *groups, size_t ng, size_t group_size, size_t cells,
    size_t count, const size_t *steps, const double *values)
{
    FILE *f = fopen("ts1_snapshots.csv", "w");
    if (!f) return 0;
    fprintf(f, "step,cell,temperature,pressure");
    size_t n = ts1_data.species_count;
    for (size_t i = 0; i < n; ++i) fprintf(f, ",%s", ts1_data.species_names[i]);
    fputc('\n', f);
    for (size_t k = 0; k < count; ++k)
        for (size_t g = 0; g < ng; ++g)
            for (size_t c = 0; c < groups[g].cells; ++c) {
                size_t cell = g * group_size + c;
                ChemState *s = groups[g].state;
                fprintf(f, "%zu,%zu,%.17g,%.17g", steps[k], cell, chem_temperatures(s)[c], chem_pressures(s)[c]);
                for (size_t i = 0; i < n; ++i) fprintf(f, ",%.17g", values[(k * cells + cell) * n + i]);
                fputc('\n', f);
            }
    int ok = !ferror(f);
    return fclose(f) == 0 && ok;
}
int main(int argc, char **argv)
{
    size_t threads = (size_t)omp_get_max_threads(), group_size = CHEM_GROUP_SIZE;
    double rtol = 1e-6, atol = -1;
    const char *path = NULL;
    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        if (!strcmp(arg, "--threads") || !strcmp(arg, "--group-size")) {
            size_t value;
            if (++i >= argc || !parse_size(argv[i], &value) || !value) goto usage;
            if (!strcmp(arg, "--threads")) threads = value; else group_size = value;
        } else if (!strcmp(arg, "--rtol") || !strcmp(arg, "--atol")) {
            double value;
            if (++i >= argc || !parse_double(argv[i], &value) || !(value > 0)) goto usage;
            if (!strcmp(arg, "--rtol")) rtol = value; else atol = value;
        } else if (path || arg[0] == '-') goto usage;
        else path = arg;
    }
    if (!path || !group_size || threads > INT_MAX) goto usage;
    int result = 1;
    Scenario scenario = {0};
    Photolysis photo = {0};
    ChemMechanism *m = NULL;
    Group *groups = NULL;
    double *base = NULL, *snapshots = NULL;
    size_t *snapshot_steps = NULL, ng = 0;
    if (!read_scenario(path, &scenario) || !read_photolysis(scenario.photo_path, &photo)) {
        fprintf(stderr, "invalid scenario or photolysis file\n"); goto cleanup;
    }
    size_t cells = count_cells(scenario.cells_path), n = ts1_data.species_count, np = ts1_data.parameter_count;
    if (!cells || cells > SIZE_MAX / sizeof(double) / n || cells > SIZE_MAX / sizeof(double) / (np + 1)) goto cleanup;
    if (group_size > cells) group_size = cells;
    ng = cells / group_size + (cells % group_size != 0);
    if (threads > ng) threads = ng;
    omp_set_num_threads((int)threads);
    if (chem_mechanism_create(&ts1_data, &m) != CHEM_OK) goto cleanup;
    groups = calloc(ng, sizeof(*groups));
    base = malloc(cells * (np + 1) * sizeof(double));
    size_t ns = scenario.snapshot_every ? scenario.steps / scenario.snapshot_every : 0;
    if (!scenario.snapshot_every || scenario.steps % scenario.snapshot_every) ++ns;
    if (ns > SIZE_MAX / cells / n / sizeof(double) || ns > SIZE_MAX / sizeof(size_t)) goto cleanup;
    snapshot_steps = malloc(ns * sizeof(size_t));
    snapshots = malloc(ns * cells * n * sizeof(double));
    if (!groups || !base || !snapshots || !snapshot_steps) goto cleanup;
    for (size_t k = 0; k < ns; ++k)
        snapshot_steps[k] = k + 1 == ns ? scenario.steps : (k + 1) * scenario.snapshot_every;
    for (size_t g = 0; g < ng; ++g) {
        groups[g].cells = cells - g * group_size < group_size ? cells - g * group_size : group_size;
        if (chem_state_create(m, groups[g].cells, &groups[g].state) != CHEM_OK) goto cleanup;
        chem_relative_tolerance(groups[g].state, rtol);
        if (atol > 0) for (size_t i = 0; i < n; ++i) chem_absolute_tolerances(groups[g].state)[i] = atol;
    }
    if (!read_cells(scenario.cells_path, groups, cells, group_size)) goto cleanup;
    for (size_t g = 0; g < ng; ++g) {
        if (chem_calculate_air_density(groups[g].state) != CHEM_OK) goto cleanup;
        memcpy(base + g * group_size * np, chem_parameters(groups[g].state), groups[g].cells * np * sizeof(double));
    }
    ChemRosenbrock parameters = chem_rosenbrock_parameters(3);
    double elapsed = 0;
    size_t begin = 0;
    for (size_t k = 0; k < ns; ++k) {
        double start = omp_get_wtime();
        #pragma omp parallel for schedule(static)
        for (size_t g = 0; g < ng; ++g) {
            Group *group = &groups[g];
            size_t nc = group->cells;
            double *custom = chem_parameters(group->state);
            for (size_t step = begin; step < snapshot_steps[k] && group->status == CHEM_OK; ++step) {
                double factor = photo_factor(&photo, (double)step * scenario.dt);
                for (size_t p = 0; p < np; ++p) if (ts1_data.parameter_names[p][0] == 'j')
                    for (size_t c = 0; c < nc; ++c)
                        custom[p * nc + c] = base[g * group_size * np + p * nc + c] * factor;
                group->status = ts1_update_rates(group->state);
                if (group->status == CHEM_OK) {
                    ChemStats stats;
                    group->status = ts1_solve(group->state, &parameters, scenario.dt, &stats);
                    add_stats(&group->stats, &stats);
                }
            }
        }
        elapsed += omp_get_wtime() - start;
        for (size_t g = 0; g < ng; ++g) {
            if (groups[g].status != CHEM_OK) {
                fprintf(stderr, "group %zu: %s\n", g, chem_status_string(groups[g].status)); goto cleanup;
            }
            size_t nc = groups[g].cells;
            const double *y = chem_concentrations(groups[g].state);
            for (size_t c = 0; c < nc; ++c)
                for (size_t i = 0; i < n; ++i)
                    snapshots[(k * cells + g * group_size + c) * n + i] = y[i * nc + c];
        }
        begin = snapshot_steps[k];
    }
    if (!write_final(groups, ng) || (scenario.snapshot_every &&
        !write_snapshots(groups, ng, group_size, cells, ns, snapshot_steps, snapshots))) goto cleanup;
    ChemStats sum = {0};
    for (size_t g = 0; g < ng; ++g) add_stats(&sum, &groups[g].stats);
    printf("elapsed_ms=%.3f\n", elapsed * 1000);
    printf("internal_steps=%zu accepted=%zu rejected=%zu function_calls=%zu jacobian_updates=%zu decompositions=%zu solves=%zu\n",
        sum.internal_steps, sum.accepted, sum.rejected, sum.function_calls, sum.jacobian_updates, sum.decompositions, sum.solves);
    result = 0;
cleanup:
    if (groups) for (size_t g = 0; g < ng; ++g) chem_state_free(groups[g].state);
    free(groups); chem_mechanism_free(m); free(base); free(snapshots); free(snapshot_steps);
    free(scenario.cells_path); free(scenario.photo_path); free(photo.times); free(photo.factors);
    return result;
usage:
    fprintf(stderr, "usage: ts1_scenario_specialized [--threads N] [--group-size N] [--rtol R] [--atol A] SCENARIO\n");
    return 2;
}
