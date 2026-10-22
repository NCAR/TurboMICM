#include "chem.h"
#include <errno.h>
#include <math.h>
#include <limits.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef CHEM_GROUP_SIZE
#define CHEM_GROUP_SIZE 32
#endif
extern const ChemMechanismData ts1_data;
#ifdef CHEM_SPECIALIZED
#include "chem_specialized.h"
#define solve_step ts1_solve
#define update_rates ts1_update_rates
#else
#define solve_step chem_solve
#define update_rates chem_update_rates
#endif
#include "bench_io.h"
int main(int argc, char **argv)
{
    size_t threads = (size_t)omp_get_max_threads(), group_size = CHEM_GROUP_SIZE;
    char *positional[4];
    int count = 0, explicit_group = 0;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--threads") || !strcmp(argv[i], "--group-size")) {
            size_t value;
            int is_threads = !strcmp(argv[i], "--threads");
            if (++i >= argc || !parse_size(argv[i], &value) || !value) goto usage;
            if (is_threads) threads = value; else { group_size = value; explicit_group = 1; }
        } else if (count < 4) positional[count++] = argv[i];
        else goto usage;
    }
    size_t cells, steps;
    double dt;
    if (count != 4 || !parse_size(positional[0], &cells) || !cells || !parse_size(positional[1], &steps)
        || !parse_double(positional[2], &dt) || !(dt > 0) || threads > INT_MAX || !group_size) goto usage;
    if (!explicit_group) {
        size_t share = cells / threads;
        if (share == 0) share = 1;
        if (group_size > share) group_size = share;
    }
    if (group_size > cells) group_size = cells;
    size_t ng = cells / group_size + (cells % group_size != 0);
    if (threads > ng) threads = ng;
    omp_set_num_threads((int)threads);
    ChemMechanism *m = NULL;
    ChemStatus status = chem_mechanism_create(&ts1_data, &m);
    if (status != CHEM_OK) { fprintf(stderr, "%s\n", chem_status_string(status)); return 1; }
    Group *groups = calloc(ng, sizeof(*groups));
    int result = 1;
    if (!groups) { chem_mechanism_free(m); return 1; }
    for (size_t g = 0; g < ng; ++g) {
        groups[g].cells = cells - g * group_size < group_size ? cells - g * group_size : group_size;
        status = chem_state_create(m, groups[g].cells, &groups[g].state);
        if (status != CHEM_OK) goto cleanup;
    }
    if (!read_cells(positional[3], groups, cells, group_size)) goto cleanup;
    ChemRosenbrock parameters = chem_rosenbrock_parameters(3);
    #pragma omp parallel for schedule(static)
    for (size_t g = 0; g < ng; ++g) {
        groups[g].status = chem_calculate_air_density(groups[g].state);
        if (groups[g].status == CHEM_OK) groups[g].status = update_rates(groups[g].state);
        ChemStats warmup;
        if (groups[g].status == CHEM_OK)
            groups[g].status = solve_step(groups[g].state, &parameters, dt, &warmup);
    }
    double start = omp_get_wtime();
    #pragma omp parallel for schedule(static)
    for (size_t g = 0; g < ng; ++g) {
        Group *group = &groups[g];
        for (size_t i = 0; i < steps && group->status == CHEM_OK; ++i) {
            ChemStats stats;
            group->status = solve_step(group->state, &parameters, dt, &stats);
            add_stats(&group->stats, &stats);
        }
    }
    double elapsed = (omp_get_wtime() - start) * 1000;
    ChemStats sum = {0};
    for (size_t g = 0; g < ng; ++g) {
        if (groups[g].status != CHEM_OK) {
            fprintf(stderr, "group %zu: %s\n", g, chem_status_string(groups[g].status)); goto cleanup;
        }
        add_stats(&sum, &groups[g].stats);
    }
    if (!write_final(groups, ng)) { fprintf(stderr, "failed writing ts1_final.csv\n"); goto cleanup; }
    printf("elapsed_ms=%.3f\n", elapsed);
    printf("internal_steps=%zu accepted=%zu rejected=%zu function_calls=%zu jacobian_updates=%zu decompositions=%zu solves=%zu\n",
        sum.internal_steps, sum.accepted, sum.rejected, sum.function_calls, sum.jacobian_updates, sum.decompositions, sum.solves);
    result = 0;
cleanup:
    for (size_t g = 0; g < ng; ++g) chem_state_free(groups[g].state);
    free(groups); chem_mechanism_free(m);
    return result;
usage:
    fprintf(stderr, "usage: ts1_bench_cpu [--threads N] [--group-size N] CELLS STEPS DT CELLS_CSV\n");
    return 2;
}
