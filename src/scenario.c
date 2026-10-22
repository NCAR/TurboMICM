#include "chem.h"
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bench_io.h"

#ifndef CHEM_GROUP_SIZE
#define CHEM_GROUP_SIZE 32
#endif
ChemModel *chem_load_ts1_specialized(int specialized);
void chem_rates_ts1_specialized(ChemState *state);
ChemStatus chem_solve_ts1_specialized(double dt, ChemState *state, ChemStats *stats);

typedef struct { char *cells, *photolysis; double dt; size_t steps, every; } Scenario;
static char *trim(char *text) {
  while (isspace((unsigned char)*text)) ++text;
  size_t n = strlen(text);
  while (n && isspace((unsigned char)text[n-1])) text[--n] = '\0';
  return text;
}
static int positive_double(const char *text, double *value) {
  char *end; errno = 0; *value = strtod(text,&end);
  return end != text && !*end && !errno && isfinite(*value) && *value > 0.0;
}
static char *relative_path(const char *scenario, const char *value) {
  const char *slash = strrchr(scenario,'/');
  size_t prefix = slash && value[0] != '/' ? (size_t)(slash-scenario)+1 : 0;
  if (strlen(value) > SIZE_MAX-prefix-1) return NULL;
  char *path = malloc(prefix+strlen(value)+1);
  if (path) { memcpy(path,scenario,prefix); strcpy(path+prefix,value); }
  return path;
}
static int scenario_read(const char *path, Scenario *s) {
  FILE *file = fopen(path,"r");
  if (!file) { perror(path); return 0; }
  unsigned seen = 0; int valid = 1;
  char *line;
  while ((line = line_read(file))) {
    char *comment = strchr(line,'#'); if (comment) *comment = '\0';
    char *key = trim(line);
    if (!*key) { free(line); continue; }
    char *equals = strchr(key,'=');
    if (!equals) { free(line); valid = 0; break; }
    *equals = '\0'; key = trim(key); char *value = trim(equals+1);
    unsigned bit = 0;
    if (!strcmp(key,"cells")) { bit = 1; if (!(seen&bit)) s->cells = relative_path(path,value); valid = *value && s->cells; }
    else if (!strcmp(key,"photolysis")) { bit = 2; if (!(seen&bit)) s->photolysis = relative_path(path,value); valid = *value && s->photolysis; }
    else if (!strcmp(key,"dt")) { bit = 4; valid = positive_double(value,&s->dt); }
    else if (!strcmp(key,"steps")) { bit = 8; valid = positive_size(value,&s->steps,0); }
    else if (!strcmp(key,"snapshot_every")) { bit = 16; valid = positive_size(value,&s->every,1); }
    else valid = 0;
    if (seen&bit) valid = 0;
    seen |= bit; free(line);
    if (!valid) break;
  }
  if (ferror(file)) valid = 0;
  fclose(file);
  if (!valid || (seen&15) != 15 || !isfinite(s->dt*(double)s->steps)) {
    fprintf(stderr,"invalid or incomplete scenario: %s\n",path); return 0;
  }
  return 1;
}
static double photo_factor(const Table *photo, double time) {
  if (time <= photo->values[0]) return photo->values[1];
  size_t last = photo->rows-1;
  if (time >= photo->values[2*last]) return photo->values[2*last+1];
  size_t lo = 1, hi = last;
  while (lo < hi) {
    size_t mid = lo+(hi-lo)/2;
    if (photo->values[2*mid] <= time) lo = mid+1; else hi = mid;
  }
  double t0 = photo->values[2*(lo-1)], t1 = photo->values[2*lo];
  double w = (time-t0)/(t1-t0);
  return (1.0-w)*photo->values[2*(lo-1)+1]+w*photo->values[2*lo+1];
}

int main(int argc, char **argv) {
  size_t threads = (size_t)omp_get_max_threads(), group_size = CHEM_GROUP_SIZE;
  double rtol = 1e-6, atol = 0.0;
  int arg = 1;
  while (arg+1 < argc && !strncmp(argv[arg],"--",2)) {
    if (!strcmp(argv[arg],"--threads")) { if (!positive_size(argv[arg+1],&threads,0)) goto usage; }
    else if (!strcmp(argv[arg],"--group-size") || !strcmp(argv[arg],"--group")) { if (!positive_size(argv[arg+1],&group_size,0)) goto usage; }
    else if (!strcmp(argv[arg],"--rtol")) { if (!positive_double(argv[arg+1],&rtol)) goto usage; }
    else if (!strcmp(argv[arg],"--atol")) { if (!positive_double(argv[arg+1],&atol)) goto usage; }
    else goto usage;
    arg += 2;
  }
  if (argc-arg != 1 || threads > INT_MAX || !group_size) goto usage;
  omp_set_num_threads((int)threads);
  if (!getenv("OMP_SCHEDULE")) omp_set_schedule(omp_sched_static,0);
  Scenario scenario = {0}; Table table = {0}, photo = {0};
  ChemModel *m = NULL; ChemState **states = NULL;
  ChemStats *statistics = NULL; ChemStatus *statuses = NULL;
  double **base = NULL; size_t *photo_ids = NULL, nphoto = 0, groups = 0;
  FILE *snapshots = NULL, *final = NULL;
  int failed = 1;
  if (!scenario_read(argv[arg],&scenario) || !table_read(scenario.cells,&table,0) ||
      !table_read(scenario.photolysis,&photo,0)) goto cleanup;
  if (photo.columns != 2 || strcmp(photo.names[0],"time_s") || strcmp(photo.names[1],"factor")) {
    fprintf(stderr,"expected time_s,factor photolysis columns\n"); goto cleanup;
  }
  for (size_t i = 0; i < photo.rows; ++i) if ((i && photo.values[2*i] <= photo.values[2*(i-1)]) || photo.values[2*i+1] < 0.0) {
    fprintf(stderr,"photolysis times must increase and factors must be nonnegative\n"); goto cleanup;
  }
  m = chem_load_ts1_specialized(1);
  if (!m || !state_columns(m,&table)) goto cleanup;
  groups = 1+(table.rows-1)/group_size;
  states = calloc(groups,sizeof *states); statistics = calloc(groups,sizeof *statistics);
  statuses = calloc(groups,sizeof *statuses); base = calloc(groups,sizeof *base);
  photo_ids = malloc((m->ncustom+1)*sizeof *photo_ids);
  if (!states || !statistics || !statuses || !base || !photo_ids) goto cleanup;
  for (size_t i = 0; i < m->ncustom; ++i) if (m->custom_names[i][0] == 'j') photo_ids[nphoto++] = i;
  #pragma omp parallel for schedule(runtime)
  for (size_t g = 0; g < groups; ++g) {
    size_t first = g*group_size, count = table.rows-first < group_size ? table.rows-first : group_size;
    ChemState *s = states[g] = chem_state_create(m,count);
    if (!s) { statuses[g] = CHEM_MEMORY; continue; }
    state_initialize(s,&table,first);
    s->relative_tolerance = rtol;
    if (atol) for (size_t i = 0; i < m->nvariables; ++i) s->absolute_tolerance[i] = atol;
    if (nphoto && count > SIZE_MAX/nphoto/sizeof(double)) { statuses[g] = CHEM_MEMORY; continue; }
    base[g] = malloc((nphoto ? nphoto : 1)*count*sizeof(double));
    if (!base[g]) { statuses[g] = CHEM_MEMORY; continue; }
    for (size_t j = 0; j < nphoto; ++j) for (size_t c = 0; c < count; ++c)
      base[g][j*count+c] = s->custom[photo_ids[j]*s->stride+c];
  }
  for (size_t g = 0; g < groups; ++g) if (statuses[g]) goto cleanup;
  if (scenario.every) {
    snapshots = fopen("ts1_snapshots.csv","w");
    if (!snapshots) { perror("ts1_snapshots.csv"); goto cleanup; }
    state_header(snapshots,m,1);
  }
  double elapsed = 0.0;
  for (size_t first = 0; first < scenario.steps;) {
    size_t remaining = scenario.steps-first;
    size_t end = first + (scenario.every && scenario.every < remaining ? scenario.every : remaining);
    double start = omp_get_wtime();
    #pragma omp parallel for schedule(runtime)
    for (size_t g = 0; g < groups; ++g) {
      ChemState *s = states[g];
      for (size_t step = first; step < end; ++step) {
        double factor = photo_factor(&photo,(double)step*scenario.dt);
        for (size_t j = 0; j < nphoto; ++j) {
          double *restrict u = s->custom+photo_ids[j]*s->stride;
          const double *restrict b = base[g]+j*s->cells;
          #pragma omp simd
          for (size_t c = 0; c < s->cells; ++c) u[c] = b[c]*factor;
        }
        chem_rates_ts1_specialized(s);
        ChemStats stats;
        statuses[g] = chem_solve_ts1_specialized(scenario.dt,s,&stats);
        add_stats(statistics+g,&stats);
        if (statuses[g]) break;
      }
    }
    elapsed += (omp_get_wtime()-start)*1000.0;
    for (size_t g = 0; g < groups; ++g) if (statuses[g]) {
      fprintf(stderr,"group %zu, steps %zu-%zu: %s\n",g,first+1,end,chem_status_string(statuses[g])); goto cleanup;
    }
    /* Snapshot formatting and file I/O are outside all timed regions. */
    if (snapshots && !state_write(snapshots,states,groups,end,1)) goto cleanup;
    first = end;
  }
  if (snapshots) {
    int error = fclose(snapshots); snapshots = NULL; if (error) goto cleanup;
  }
  final = fopen("ts1_final.csv","w");
  if (!final) { perror("ts1_final.csv"); goto cleanup; }
  state_header(final,m,0);
  if (!state_write(final,states,groups,scenario.steps,0)) goto cleanup;
  { int error = fclose(final); final = NULL; if (error) goto cleanup; }
  ChemStats sum = {0};
  for (size_t g = 0; g < groups; ++g) add_stats(&sum,statistics+g);
  print_stats(elapsed,&sum); failed = 0;
cleanup:
  if (final) fclose(final);
  if (snapshots) fclose(snapshots);
  for (size_t g = 0; g < groups; ++g) { if (states) chem_state_free(states[g]); if (base) free(base[g]); }
  free(states); free(statistics); free(statuses); free(base); free(photo_ids); chem_model_free(m);
  table_free(&table); table_free(&photo); free(scenario.cells); free(scenario.photolysis);
  return failed;
usage:
  fprintf(stderr,"usage: %s [--threads N] [--group-size N] [--rtol R] [--atol A] SCENARIO\n",argv[0]);
  return 2;
}
