#include "chem.h"
#include <errno.h>
#include <math.h>
#include <limits.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

ChemModel *chem_load_ts1(int specialized);
#ifdef CHEM_SPECIALIZED
ChemModel *chem_load_ts1_specialized(int specialized);
void chem_rates_ts1_specialized(ChemState *state);
ChemStatus chem_solve_ts1_specialized(double dt, ChemState *state, ChemStats *stats);
#define chem_load_ts1 chem_load_ts1_specialized
#define chem_calculate_rates chem_rates_ts1_specialized
#define BENCH_SOLVE(solver, dt, state, stats) chem_solve_ts1_specialized(dt, state, stats)
#define BENCH_SOLVER(name) ((void)0)
#else
#define BENCH_SOLVE(solver, dt, state, stats) chem_solve(solver, dt, state, stats)
#define BENCH_SOLVER(name) ChemSolver name = chem_solver(3)
#endif
#ifndef CHEM_GROUP_SIZE
#define CHEM_GROUP_SIZE 32
#endif

#include "bench_io.h"

int main(int argc, char **argv) {
  size_t threads = (size_t)omp_get_max_threads(), group_size = CHEM_GROUP_SIZE, cells, steps;
  int arg = 1;
  while (arg+1 < argc && !strncmp(argv[arg], "--", 2)) {
    size_t v; if (!positive_size(argv[arg+1], &v, 0)) goto usage;
    if (!strcmp(argv[arg], "--threads")) threads = v;
    else if (!strcmp(argv[arg], "--group") || !strcmp(argv[arg], "--group-size")) group_size = v;
    else goto usage;
    arg += 2;
  }
  if (argc-arg != 4 || !positive_size(argv[arg], &cells, 0) || !positive_size(argv[arg+1], &steps, 1) ||
      threads > INT_MAX || !group_size) goto usage;
  char *end; errno = 0; double dt = strtod(argv[arg+2], &end);
  if (errno || *end || !isfinite(dt) || dt <= 0.0) goto usage;
  omp_set_num_threads((int)threads);
  if (!getenv("OMP_SCHEDULE")) omp_set_schedule(omp_sched_static,0);
  Table table = {0};
  if (!table_read(argv[arg+3], &table, cells)) { table_free(&table); return 1; }
  ChemModel *m = chem_load_ts1(1);
  if (!m) { fprintf(stderr, "cannot build mechanism\n"); table_free(&table); return 1; }
  if (!state_columns(m,&table)) { table_free(&table); chem_model_free(m); return 1; }
  size_t groups = 1+(cells-1)/group_size;
  ChemState **states = calloc(groups, sizeof *states);
  ChemStats *statistics = calloc(groups, sizeof *statistics);
  ChemStatus *statuses = calloc(groups, sizeof *statuses);
  if (!states || !statistics || !statuses) { free(states); free(statistics); free(statuses); table_free(&table); chem_model_free(m); return 1; }
  int failed = 0;
  #pragma omp parallel for schedule(runtime)
  for (size_t g = 0; g < groups; ++g) {
    size_t first = g*group_size, count = cells-first < group_size ? cells-first : group_size;
    ChemState *s = states[g] = chem_state_create(m, count);
    if (!s) { statuses[g] = CHEM_MEMORY; continue; }
    state_initialize(s,&table,first);
    chem_calculate_rates(s);
    BENCH_SOLVER(solver); ChemStats warmup;
    statuses[g] = BENCH_SOLVE(&solver, dt, s, &warmup);
  }
  for (size_t g = 0; g < groups; ++g) if (statuses[g]) { fprintf(stderr, "warm-up group %zu: %s\n", g, chem_status_string(statuses[g])); failed = 1; }
  double start = omp_get_wtime();
  if (!failed) {
    #pragma omp parallel for schedule(runtime)
    for (size_t g = 0; g < groups; ++g) {
      BENCH_SOLVER(solver);
      for (size_t step = 0; step < steps; ++step) {
        ChemStats stats;
        statuses[g] = BENCH_SOLVE(&solver, dt, states[g], &stats);
        add_stats(&statistics[g], &stats);
        if (statuses[g]) break;
      }
    }
  }
  double elapsed = (omp_get_wtime()-start)*1000.0;
  ChemStats sum = {0};
  for (size_t g = 0; g < groups; ++g) {
    add_stats(&sum, &statistics[g]);
    if (statuses[g]) { fprintf(stderr, "group %zu: %s\n", g, chem_status_string(statuses[g])); failed = 1; }
  }
  if (!failed) {
    print_stats(elapsed,&sum);
    FILE *out = fopen("ts1_final.csv", "w");
    if (!out) { perror("ts1_final.csv"); failed = 1; }
    else {
      state_header(out,m,0);
      if (!state_write(out,states,groups,steps,0)) failed = 1;
      if (fclose(out)) failed = 1;
    }
  }
  for (size_t g = 0; g < groups; ++g) chem_state_free(states[g]);
  free(states); free(statistics); free(statuses); table_free(&table); chem_model_free(m);
  return failed ? 1 : 0;
usage:
  fprintf(stderr, "usage: %s [--threads N] [--group N] CELLS STEPS DT CELLS_CSV\n", argv[0]);
  return 2;
}
