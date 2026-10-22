/* SPDX-License-Identifier: Apache-2.0 */
#include "ts1_io.h"

int main(int argc, char **argv) {
  const char *pos[4]; int nargs = 0;
  size_t threads = (size_t)omp_get_max_threads(), group_size = CH_DEFAULT_GROUP_SIZE;
  const char *group_env = getenv("CH_GROUP_SIZE");
  if (group_env && (!size_parse(group_env,&group_size) || !group_size)) goto usage;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i],"--threads") || !strcmp(argv[i],"--group-size")) {
      int is_threads = !strcmp(argv[i],"--threads");
      size_t value;
      if (++i == argc || !size_parse(argv[i],&value) || !value) goto usage;
      if (is_threads) threads = value; else group_size = value;
    } else {
      if (nargs == 4) goto usage;
      pos[nargs++] = argv[i];
    }
  }
  size_t cells, steps; double dt;
  if (nargs != 4 || !size_parse(pos[0],&cells) || !cells || !size_parse(pos[1],&steps) ||
      !double_parse(pos[2],&dt) || dt <= 0 || !threads || threads > INT_MAX || !group_size) goto usage;
  if (group_size > cells) group_size = cells;
  size_t ngroups = cells/group_size + (cells%group_size != 0);
  if (threads > ngroups) threads = ngroups;
  omp_set_num_threads((int)threads);
  /* Groups have independent adaptive histories. Dynamic assignment balances
   * heterogeneous cores without changing any group's arithmetic. */
  if (!getenv("OMP_SCHEDULE")) omp_set_schedule(omp_sched_dynamic,1);
  ch_model *m = ch_build_ts1();
  group *groups = calloc(ngroups,sizeof(*groups));
  int status = CH_MEMORY;
  if (!m || !groups) goto cleanup;
  for (size_t g = 0; g < ngroups; ++g) {
    size_t count = cells-g*group_size;
    if (count > group_size) count = group_size;
    groups[g].first = g*group_size;
    groups[g].state = ch_state_create(m,count);
    if (!groups[g].state) { status = CH_MEMORY; goto cleanup; }
    status = ch_init_ts1(groups[g].state);
    if (status) goto cleanup;
  }
  status = initialize_csv(pos[3],m,groups,cells,group_size);
  if (status) goto cleanup;
  ch_rosenbrock_parameters p = ch_rosenbrock_default(3);
  #pragma omp parallel for schedule(runtime)
  for (size_t g = 0; g < ngroups; ++g) {
    ch_stats warmup;
    groups[g].status = ch_update_rates(groups[g].state);
    if (!groups[g].status) groups[g].status = ch_solve(groups[g].state,dt,&p,&warmup);
  }
  for (size_t g = 0; g < ngroups; ++g) if ((status = groups[g].status)) goto cleanup;
  double begin = omp_get_wtime();
  #pragma omp parallel for schedule(runtime)
  for (size_t g = 0; g < ngroups; ++g) {
    for (size_t step = 0; step < steps; ++step) {
      ch_stats stats;
      groups[g].status = ch_solve(groups[g].state,dt,&p,&stats);
      ch_stats_add(&groups[g].stats,&stats);
      if (groups[g].status) break;
    }
  }
  double elapsed_ms = (omp_get_wtime()-begin)*1000.0;
  ch_stats total = {0};
  for (size_t g = 0; g < ngroups; ++g) {
    if ((status = groups[g].status)) goto cleanup;
    ch_stats_add(&total,&groups[g].stats);
  }
  FILE *file = fopen("ts1_final.csv","w");
  if (!file) { status = CH_INVALID; goto cleanup; }
  fprintf(file,"temperature,pressure");
  size_t nv = ch_model_species_count(m);
  for (size_t i = 0; i < nv; ++i) fprintf(file,",%s",ch_model_species_name(m,i));
  fputc('\n',file);
  for (size_t g = 0; g < ngroups; ++g) {
    ch_state *s = groups[g].state;
    size_t count = ch_state_cells(s);
    double *y = ch_state_values(s);
    ch_conditions *c = ch_state_conditions(s);
    for (size_t cell = 0; cell < count; ++cell) {
      fprintf(file,"%.17g,%.17g",c[cell].temperature,c[cell].pressure);
      for (size_t i = 0; i < nv; ++i) fprintf(file,",%.17g",y[i*count+cell]);
      fputc('\n',file);
    }
  }
  int write_failed = ferror(file);
  if (fclose(file) || write_failed) { status = CH_INVALID; goto cleanup; }
  printf("elapsed_ms=%.3f\ninternal_steps=%zu accepted=%zu rejected=%zu function_calls=%zu jacobian_updates=%zu decompositions=%zu solves=%zu\n",
    elapsed_ms,total.internal_steps,total.accepted,total.rejected,total.function_calls,total.jacobian_updates,total.decompositions,total.solves);
  status = CH_OK;
cleanup:
  if (groups) for (size_t g = 0; g < ngroups; ++g) ch_state_free(groups[g].state);
  free(groups); ch_model_free(m);
  if (status) fprintf(stderr,"TS1 benchmark: %s\n",ch_status_string(status));
  return status ? 1 : 0;
usage:
  fprintf(stderr,"usage: ts1_bench_cpu [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV\n");
  return 2;
}
