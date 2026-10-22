/* SPDX-License-Identifier: Apache-2.0 */
#include "ts1_io.h"
#include <ctype.h>

typedef struct { char *cells, *photolysis; size_t steps, every; double dt; } scenario;
typedef struct { double time, factor; } point;

static char *trim(char *s) {
  while (isspace((unsigned char)*s)) ++s;
  size_t n = strlen(s);
  while (n && isspace((unsigned char)s[n-1])) s[--n] = 0;
  return s;
}
static char *relative_path(const char *file, const char *value) {
  const char *slash = strrchr(file,'/');
  size_t prefix = value[0] == '/' || !slash ? 0 : (size_t)(slash-file)+1;
  if (strlen(value) >= SIZE_MAX-prefix) return NULL;
  char *result = malloc(prefix+strlen(value)+1);
  if (result) { memcpy(result,file,prefix); strcpy(result+prefix,value); }
  return result;
}
static int read_scenario(const char *path, scenario *s) {
  FILE *file = fopen(path,"r");
  if (!file) return CH_INVALID;
  unsigned seen = 0;
  int status = CH_INVALID;
  char *line;
  while ((line = line_read(file))) {
    char *comment = strchr(line,'#');
    if (comment) *comment = 0;
    char *key = trim(line), *eq = strchr(key,'=');
    if (!*key) { free(line); continue; }
    if (!eq) { free(line); goto done; }
    *eq = 0;
    key = trim(key);
    char *value = trim(eq+1);
    unsigned bit = 0;
    if (!strcmp(key,"cells")) {
      bit = 1;
      if (!(seen&bit)) s->cells = relative_path(path,value);
      if (!s->cells || !*value) { free(line); goto done; }
    } else if (!strcmp(key,"photolysis")) {
      bit = 2;
      if (!(seen&bit)) s->photolysis = relative_path(path,value);
      if (!s->photolysis || !*value) { free(line); goto done; }
    } else if (!strcmp(key,"dt")) {
      bit = 4;
      if (!double_parse(value,&s->dt) || s->dt <= 0) { free(line); goto done; }
    } else if (!strcmp(key,"steps")) {
      bit = 8;
      if (!size_parse(value,&s->steps) || !s->steps) { free(line); goto done; }
    } else if (!strcmp(key,"snapshot_every")) {
      bit = 16;
      if (!size_parse(value,&s->every)) { free(line); goto done; }
    }
    free(line);
    if (!bit || (seen&bit)) goto done;
    seen |= bit;
  }
  if ((seen&15) == 15 && isfinite(s->dt*s->steps) && !ferror(file)) status = CH_OK;
done:
  fclose(file); return status;
}
static size_t cell_count(const char *path) {
  FILE *file = fopen(path,"r");
  if (!file) return 0;
  char *line = line_read(file);
  free(line);
  size_t count = 0;
  while ((line = line_read(file))) { ++count; free(line); }
  if (ferror(file)) count = 0;
  fclose(file); return count;
}
static int read_photolysis(const char *path, point **points, size_t *count) {
  FILE *file = fopen(path,"r");
  if (!file) return CH_INVALID;
  char *line = line_read(file);
  int status = CH_INVALID;
  if (!line || strcmp(line,"time_s,factor")) { free(line); goto done; }
  free(line);
  while ((line = line_read(file))) {
    char *cursor = line, *time = field_next(&cursor), *factor = field_next(&cursor);
    point p;
    int valid = factor && !cursor && double_parse(time,&p.time) && double_parse(factor,&p.factor) &&
                p.factor >= 0 && (!*count || p.time > (*points)[*count-1].time);
    free(line);
    if (!valid) goto done;
    if (*count >= SIZE_MAX/sizeof(point)-1) { status = CH_MEMORY; goto done; }
    point *next = realloc(*points,(*count+1)*sizeof(point));
    if (!next) { status = CH_MEMORY; goto done; }
    *points = next; next[(*count)++] = p;
  }
  if (*count && !ferror(file)) status = CH_OK;
done:
  fclose(file); return status;
}
static double factor_at(const point *p, size_t n, double time) {
  if (time <= p[0].time) return p[0].factor;
  if (time >= p[n-1].time) return p[n-1].factor;
  size_t lo = 0, hi = n-1;
  while (hi-lo > 1) { size_t mid = lo+(hi-lo)/2; if (p[mid].time <= time) lo = mid; else hi = mid; }
  double weight = (time-p[lo].time)/(p[hi].time-p[lo].time);
  return (1.0-weight)*p[lo].factor + weight*p[hi].factor;
}
static int output(const char *path, const ch_model *m, group *groups, size_t ngroups,
                  size_t cells, const double *snapshots, const size_t *steps, size_t nsnap, int with_step) {
  FILE *file = fopen(path,"w");
  if (!file) return CH_INVALID;
  size_t nv = ch_model_species_count(m);
  if (with_step) fprintf(file,"step,cell,");
  fprintf(file,"temperature,pressure");
  for (size_t i = 0; i < nv; ++i) fprintf(file,",%s",ch_model_species_name(m,i));
  fputc('\n',file);
  for (size_t snap = with_step ? 0 : nsnap-1; snap < nsnap; ++snap) {
    for (size_t g = 0; g < ngroups; ++g) {
      ch_state *s = groups[g].state;
      size_t nc = ch_state_cells(s);
      const double *y = snapshots + snap*cells*nv + groups[g].first*nv;
      const ch_conditions *cond = ch_state_conditions(s);
      for (size_t c = 0; c < nc; ++c) {
        if (with_step) fprintf(file,"%zu,%zu,",steps[snap],groups[g].first+c);
        fprintf(file,"%.17g,%.17g",cond[c].temperature,cond[c].pressure);
        for (size_t i = 0; i < nv; ++i) fprintf(file,",%.17g",y[i*nc+c]);
        fputc('\n',file);
      }
    }
  }
  int failed = ferror(file);
  if (fclose(file)) failed = 1;
  return failed ? CH_INVALID : CH_OK;
}

int main(int argc, char **argv) {
  size_t threads = (size_t)omp_get_max_threads(), group_size = CH_DEFAULT_GROUP_SIZE;
  double rtol = 1e-6, atol = 1e-3;
  const char *path = NULL, *env = getenv("CH_GROUP_SIZE");
  if (env && (!size_parse(env,&group_size) || !group_size)) goto usage;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i],"--threads") || !strcmp(argv[i],"--group-size")) {
      int is_threads = !strcmp(argv[i],"--threads");
      size_t value;
      if (++i == argc || !size_parse(argv[i],&value) || !value) goto usage;
      if (is_threads) threads = value; else group_size = value;
    } else if (!strcmp(argv[i],"--rtol") || !strcmp(argv[i],"--atol")) {
      int is_rtol = !strcmp(argv[i],"--rtol");
      double value;
      if (++i == argc || !double_parse(argv[i],&value) || value <= 0) goto usage;
      if (is_rtol) rtol = value; else atol = value;
    } else if (path || argv[i][0] == '-') goto usage;
    else path = argv[i];
  }
  if (!path || !threads || threads > INT_MAX || !group_size) goto usage;
  scenario spec = {0};
  point *points = NULL;
  size_t npoints = 0, ngroups = 0;
  ch_model *m = NULL;
  group *groups = NULL;
  double *base = NULL, *snapshots = NULL;
  size_t *snapshot_steps = NULL, *photo = NULL;
  int status = read_scenario(path,&spec);
  if (status) goto cleanup;
  size_t cells = cell_count(spec.cells);
  if (!cells) { status = CH_INVALID; goto cleanup; }
  status = read_photolysis(spec.photolysis,&points,&npoints);
  if (status) goto cleanup;
  m = ch_build_ts1();
  status = CH_MEMORY;
  if (!m) goto cleanup;
  size_t nv = ch_model_species_count(m), np = ch_model_parameter_count(m);
  size_t nsnap = spec.every ? spec.steps/spec.every + (spec.steps%spec.every != 0) : 1;
  if (cells > SIZE_MAX/nv/sizeof(double) || nsnap > SIZE_MAX/(cells*nv)/sizeof(double) ||
      nsnap > SIZE_MAX/sizeof(size_t) || (np && cells > SIZE_MAX/np/sizeof(double))) goto cleanup;
  snapshot_steps = malloc(nsnap*sizeof(size_t));
  snapshots = malloc(nsnap*cells*nv*sizeof(double));
  base = malloc((np ? np : 1)*cells*sizeof(double));
  photo = malloc((np ? np : 1)*sizeof(size_t));
  if (!snapshot_steps || !snapshots || !base || !photo) goto cleanup;
  for (size_t i = 0; i < nsnap; ++i) snapshot_steps[i] = i+1 == nsnap ? spec.steps : (i+1)*spec.every;
  size_t nphoto = 0;
  for (size_t i = 0; i < np; ++i) if (ch_model_parameter_name(m,i)[0] == 'j') photo[nphoto++] = i;
  if (group_size > cells) group_size = cells;
  ngroups = cells/group_size+(cells%group_size != 0);
  groups = calloc(ngroups,sizeof(*groups));
  if (!groups) goto cleanup;
  for (size_t g = 0; g < ngroups; ++g) {
    size_t count = cells-g*group_size;
    if (count > group_size) count = group_size;
    groups[g].first = g*group_size;
    groups[g].state = ch_state_create(m,count);
    if (!groups[g].state) { status = CH_MEMORY; goto cleanup; }
    status = ch_init_ts1(groups[g].state);
    if (status) goto cleanup;
    status = ch_state_tolerances(groups[g].state,rtol,&atol,1);
    if (status) goto cleanup;
  }
  status = initialize_csv(spec.cells,m,groups,cells,group_size);
  if (status) goto cleanup;
  for (size_t g = 0; g < ngroups; ++g)
    memcpy(base+groups[g].first*np,ch_state_parameters(groups[g].state),np*ch_state_cells(groups[g].state)*sizeof(double));
  if (threads > ngroups) threads = ngroups;
  omp_set_num_threads((int)threads);
  if (!getenv("OMP_SCHEDULE")) omp_set_schedule(omp_sched_dynamic,1);
  ch_rosenbrock_parameters parameters = ch_rosenbrock_default(3);
  /* Stop the wall timer for snapshot copies; all parameter setting and solves
   * are timed. Segments end together at the next requested snapshot. */
  double elapsed_ms = 0.0;
  size_t first_step = 0;
  for (size_t snap = 0; snap < nsnap; ++snap) {
    double begin = omp_get_wtime();
    #pragma omp parallel for schedule(runtime)
    for (size_t g = 0; g < ngroups; ++g) {
      ch_state *s = groups[g].state;
      size_t nc = ch_state_cells(s);
      double *custom = ch_state_parameters(s);
      const double *initial = base+groups[g].first*np;
      for (size_t step = first_step; step < snapshot_steps[snap]; ++step) {
        double factor = factor_at(points,npoints,step*spec.dt);
        for (size_t j = 0; j < nphoto; ++j) {
          size_t offset = photo[j]*nc;
          #pragma omp simd
          for (size_t c = 0; c < nc; ++c) custom[offset+c] = initial[offset+c]*factor;
        }
        ch_stats stats = {0};
        groups[g].status = ch_update_rates(s);
        if (!groups[g].status) groups[g].status = ch_solve(s,spec.dt,&parameters,&stats);
        ch_stats_add(&groups[g].stats,&stats);
        if (groups[g].status) break;
      }
    }
    elapsed_ms += (omp_get_wtime()-begin)*1000.0;
    for (size_t g = 0; g < ngroups; ++g) {
      if ((status = groups[g].status)) goto cleanup;
      memcpy(snapshots+snap*cells*nv+groups[g].first*nv,ch_state_values(groups[g].state),
             nv*ch_state_cells(groups[g].state)*sizeof(double));
    }
    first_step = snapshot_steps[snap];
  }
  status = output("ts1_final.csv",m,groups,ngroups,cells,snapshots,snapshot_steps,nsnap,0);
  if (!status) status = output("ts1_snapshots.csv",m,groups,ngroups,cells,snapshots,snapshot_steps,nsnap,1);
  if (status) goto cleanup;
  ch_stats total = {0};
  for (size_t g = 0; g < ngroups; ++g) ch_stats_add(&total,&groups[g].stats);
  printf("elapsed_ms=%.3f\ninternal_steps=%zu accepted=%zu rejected=%zu function_calls=%zu jacobian_updates=%zu decompositions=%zu solves=%zu\n",
         elapsed_ms,total.internal_steps,total.accepted,total.rejected,total.function_calls,total.jacobian_updates,total.decompositions,total.solves);
cleanup:
  if (groups) for (size_t g = 0; g < ngroups; ++g) ch_state_free(groups[g].state);
  free(groups); ch_model_free(m); free(base); free(snapshots); free(snapshot_steps); free(photo);
  free(points); free(spec.cells); free(spec.photolysis);
  if (status) fprintf(stderr,"TS1 scenario: %s\n",ch_status_string(status));
  return status ? 1 : 0;
usage:
  fprintf(stderr,"usage: ts1_scenario_specialized [--threads N] [--group-size G] [--rtol R] [--atol A] SCENARIO\n");
  return 2;
}
