/* SPDX-License-Identifier: Apache-2.0 */
#include "chem.h"
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef CH_DEFAULT_GROUP_SIZE
#define CH_DEFAULT_GROUP_SIZE 32
#endif

ch_model *ch_build_ts1(void);
int ch_init_ts1(ch_state *);

typedef struct { ch_state *state; size_t first; int status; ch_stats stats; } group;
typedef struct { char *name; int kind, index; } column;

static char *line_read(FILE *file) {
  size_t n = 0, cap = 1024;
  char *line = malloc(cap);
  if (!line) return NULL;
  int ch;
  while ((ch = fgetc(file)) != EOF && ch != '\n') {
    if (n+1 == cap) {
      char *next = realloc(line,cap*2);
      if (!next) { free(line); return NULL; }
      line = next; cap *= 2;
    }
    line[n++] = (char)ch;
  }
  if (ch == EOF && n == 0) { free(line); return NULL; }
  if (n && line[n-1] == '\r') --n;
  line[n] = 0;
  return line;
}

/* The reference format has unquoted names and decimal values, no embedded commas. */
static char *field_next(char **cursor) {
  if (!*cursor) return NULL;
  char *start = *cursor, *comma = strchr(start,',');
  if (comma) { *comma = 0; *cursor = comma+1; }
  else *cursor = NULL;
  return start;
}
static int size_parse(const char *text, size_t *value) {
  char *end; errno = 0;
  if (!text[0] || text[0] == '-') return 0;
  unsigned long long n = strtoull(text,&end,10);
  if (errno || *end || n != (size_t)n) return 0;
  *value = (size_t)n; return 1;
}
static int double_parse(const char *text, double *value) {
  char *end; errno = 0; *value = strtod(text,&end);
  return end != text && !*end && isfinite(*value) && errno != ERANGE;
}
static int initialize_csv(const char *path, const ch_model *m, group *groups, size_t cells, size_t group_size) {
  FILE *file = fopen(path,"r");
  if (!file) { fprintf(stderr,"Cannot open input: %s\n",path); return CH_INVALID; }
  char *header = line_read(file), *line = NULL;
  column *cols = NULL;
  size_t ncols = 0;
  int status = CH_INVALID;
  if (!header) goto done;
  char *cursor = header, *field;
  while ((field = field_next(&cursor))) {
    column col = {field,0,0};
    if (!strcmp(field,"temperature")) col.kind = 0;
    else if (!strcmp(field,"pressure")) col.kind = 1;
    else if ((col.index = ch_model_species_index(m,field)) >= 0) col.kind = 2;
    else if ((col.index = ch_model_parameter_index(m,field)) >= 0) col.kind = 3;
    else { fprintf(stderr,"Unknown input column: %s\n",field); goto done; }
    for (size_t i = 0; i < ncols; ++i) if (!strcmp(cols[i].name,field)) goto done;
    column *next = realloc(cols,(ncols+1)*sizeof(*cols));
    if (!next) { status = CH_MEMORY; goto done; }
    cols = next; cols[ncols++] = col;
  }
  for (size_t c = 0; c < cells; ++c) {
    line = line_read(file);
    if (!line) { fprintf(stderr,"Input has fewer than %zu cells\n",cells); goto done; }
    cursor = line;
    ch_state *s = groups[c/group_size].state;
    size_t local = c-groups[c/group_size].first, count = ch_state_cells(s);
    for (size_t i = 0; i < ncols; ++i) {
      field = field_next(&cursor);
      double value;
      if (!field || !double_parse(field,&value)) {
        fprintf(stderr,"Invalid value in cell %zu, column %s\n",c,cols[i].name); goto done;
      }
      switch (cols[i].kind) {
        case 0: ch_state_conditions(s)[local].temperature = value; break;
        case 1: ch_state_conditions(s)[local].pressure = value; break;
        case 2: ch_state_values(s)[(size_t)cols[i].index*count+local] = value; break;
        case 3: ch_state_parameters(s)[(size_t)cols[i].index*count+local] = value; break;
      }
    }
    if (cursor) { fprintf(stderr,"Extra value in cell %zu\n",c); goto done; }
    ch_conditions *condition = &ch_state_conditions(s)[local];
    if (condition->temperature <= 0 || condition->pressure <= 0) goto done;
    ch_ideal_air_density(condition);
    free(line); line = NULL;
  }
  status = CH_OK;
done:
  free(line); free(header); free(cols); fclose(file); return status;
}

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
  ch_model *m = ch_build_ts1();
  group *groups = calloc(ngroups,sizeof(*groups));
  int status = CH_MEMORY;
  if (!m || !groups) goto cleanup;
  for (size_t g = 0; g < ngroups; ++g) {
    size_t count = cells-g*group_size;
    if (count > group_size) count = group_size;
    groups[g].first = g*group_size;
    groups[g].state = ch_state_create(m,count);
    if (!groups[g].state) goto cleanup;
    status = ch_init_ts1(groups[g].state);
    if (status) goto cleanup;
  }
  status = initialize_csv(pos[3],m,groups,cells,group_size);
  if (status) goto cleanup;
  ch_rosenbrock_parameters p = ch_rosenbrock_default(3);
  #pragma omp parallel for schedule(static)
  for (size_t g = 0; g < ngroups; ++g) {
    ch_stats warmup;
    groups[g].status = ch_update_rates(groups[g].state);
    if (!groups[g].status) groups[g].status = ch_solve(groups[g].state,dt,&p,&warmup);
  }
  for (size_t g = 0; g < ngroups; ++g) if ((status = groups[g].status)) goto cleanup;
  double begin = omp_get_wtime();
  #pragma omp parallel for schedule(static)
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
