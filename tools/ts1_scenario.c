/* Scenario runner; no numerical method depends on TS1 or these input formats.
 * SPDX-License-Identifier: Apache-2.0 */
#include "chem.h"
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef CH_DEFAULT_GROUP_SIZE
#define CH_DEFAULT_GROUP_SIZE 32
#endif
ch_model *ch_build_ts1(void);
int ch_init_ts1(ch_state *);

typedef struct {
  ch_state *state;
  size_t first;
  int status;
  ch_stats stats;
  ch_integrator_stats extra;
  ch_integrator *integrator;
  double *base_parameters;
} group;
#include "ts1_io.h"

typedef struct {
  char *cells, *photolysis;
  double dt;
  size_t steps, every;
} scenario;
typedef struct { double time, factor; } sunlight;

static char *trim(char *s) {
  while (isspace((unsigned char)*s)) ++s;
  size_t n = strlen(s);
  while (n && isspace((unsigned char)s[n-1])) s[--n] = 0;
  return s;
}
static char *relative_path(const char *file, const char *value) {
  const char *slash = strrchr(file,'/');
  size_t prefix = value[0] == '/' || !slash ? 0 : (size_t)(slash-file)+1;
  char *result = malloc(prefix+strlen(value)+1);
  if (result) { memcpy(result,file,prefix); strcpy(result+prefix,value); }
  return result;
}
static int read_scenario(const char *path, scenario *s) {
  FILE *f = fopen(path,"r");
  if (!f) return 0;
  char *line;
  unsigned seen = 0;
  int ok = 1;
  while ((line = line_read(f))) {
    char *comment = strchr(line,'#'); if (comment) *comment = 0;
    char *key = trim(line);
    if (!*key) { free(line); continue; }
    char *value = strchr(key,'=');
    if (!value) { free(line); ok = 0; break; }
    *value++ = 0; key = trim(key); value = trim(value);
    unsigned bit = 0;
    if (!strcmp(key,"cells")) { bit = 1; if (!(seen&bit)) s->cells = relative_path(path,value); }
    else if (!strcmp(key,"photolysis")) { bit = 2; if (!(seen&bit)) s->photolysis = relative_path(path,value); }
    else if (!strcmp(key,"dt")) { bit = 4; ok = double_parse(value,&s->dt) && s->dt > 0; }
    else if (!strcmp(key,"steps")) { bit = 8; ok = size_parse(value,&s->steps) && s->steps > 0; }
    else if (!strcmp(key,"snapshot_every")) { bit = 16; ok = size_parse(value,&s->every); }
    else ok = 0;
    if ((seen&bit) || !*value) ok = 0;
    seen |= bit;
    free(line);
    if (!ok) break;
  }
  if (ferror(f)) ok = 0;
  fclose(f);
  return ok && (seen&15) == 15 && s->cells && s->photolysis && isfinite(s->dt*s->steps);
}
static size_t count_cells(const char *path) {
  FILE *f = fopen(path,"r");
  if (!f) return 0;
  char *line = line_read(f); free(line);
  size_t n = 0;
  while ((line = line_read(f))) { if (*line) ++n; free(line); }
  if (ferror(f)) n = 0;
  fclose(f); return n;
}
static int read_sunlight(const char *path, sunlight **points, size_t *count) {
  FILE *f = fopen(path,"r");
  if (!f) return 0;
  char *line = line_read(f);
  int ok = line && !strcmp(line,"time_s,factor"); free(line);
  while (ok && (line = line_read(f))) {
    char *cursor = line, *a = field_next(&cursor), *b = field_next(&cursor);
    sunlight p;
    ok = a && b && !cursor && double_parse(a,&p.time) && double_parse(b,&p.factor) && p.factor >= 0;
    if (ok && *count && p.time <= (*points)[*count-1].time) ok = 0;
    if (ok) {
      sunlight *next = realloc(*points,(*count+1)*sizeof(**points));
      if (!next) ok = 0;
      else { *points = next; (*points)[(*count)++] = p; }
    }
    free(line);
  }
  if (ferror(f)) ok = 0;
  fclose(f); return ok && *count;
}
static double sun_factor(const sunlight *p, size_t n, double time) {
  if (time <= p[0].time) return p[0].factor;
  if (time >= p[n-1].time) return p[n-1].factor;
  size_t low = 0, high = n-1;
  while (high-low > 1) { size_t mid = low+(high-low)/2; if (p[mid].time <= time) low = mid; else high = mid; }
  double w = (time-p[low].time)/(p[high].time-p[low].time);
  return (1-w)*p[low].factor+w*p[high].factor;
}

static int method_id(const char *name, int *legacy) {
  const char *old[] = {"rosenbrock2","rosenbrock3","rosenbrock4","rosenbrock4dae","rosenbrock6dae"};
  const char *new_names[] = {"ros2","ros3","ros4","rodas3","rodas4"};
  const int ids[] = {2,3,4,43,6};
  for (int i = 0; i < 5; ++i) {
    if (!strcmp(name,old[i])) { *legacy = 1; return ids[i]; }
    if (!strcmp(name,new_names[i])) { *legacy = 0; return ids[i]; }
  }
  *legacy = 0;
  if (!strcmp(name,"rodas31")) return 431;
  return !strcmp(name,"sdirk2") ? 102 : 0;
}
static int option(const char *arg, ch_integrator_options *o, size_t *group_size, double *rtol, double *atol) {
  char key[64];
  const char *eq = strchr(arg,'=');
  if (!eq || eq == arg || (size_t)(eq-arg) >= sizeof(key)) return 0;
  memcpy(key,arg,(size_t)(eq-arg)); key[eq-arg] = 0;
  const char *value = eq+1;
  if (!strcmp(key,"controller")) {
    if (!strcmp(value,"i")) o->controller = 0;
    else if (!strcmp(value,"pi")) o->controller = 1;
    else if (!strcmp(value,"h211b")) o->controller = 2;
    else return 0;
    return 1;
  }
  if (!strcmp(key,"norm")) {
    if (!strcmp(value,"rms")) o->norm = 0;
    else if (!strcmp(value,"cell")) o->norm = 1;
    else if (!strcmp(value,"max")) o->norm = 2;
    else return 0;
    return 1;
  }
  size_t integer;
  if (!strcmp(key,"group_size")) return size_parse(value,group_size) && *group_size;
  int *dest = NULL;
  if (!strcmp(key,"jacobian_steps")) dest = &o->jacobian_steps;
  else if (!strcmp(key,"warm_start")) dest = &o->warm_start;
  else if (!strcmp(key,"max_steps")) dest = &o->max_steps;
  else if (!strcmp(key,"newton_iterations")) dest = &o->newton_iterations;
  if (dest) {
    if (!size_parse(value,&integer) || integer > INT_MAX) return 0;
    *dest = (int)integer; return 1;
  }
  double *real = NULL;
  if (!strcmp(key,"rtol")) real = rtol;
  else if (!strcmp(key,"atol")) real = atol;
  else if (!strcmp(key,"h_start")) real = &o->h_start;
  else if (!strcmp(key,"h_max")) real = &o->h_max;
  else if (!strcmp(key,"safety")) real = &o->safety;
  else if (!strcmp(key,"factor_min")) real = &o->factor_min;
  else if (!strcmp(key,"factor_max")) real = &o->factor_max;
  else if (!strcmp(key,"pi_beta")) real = &o->pi_beta;
  else if (!strcmp(key,"lu_band")) real = &o->lu_band;
  else if (!strcmp(key,"positivity")) real = &o->positivity;
  else if (!strcmp(key,"newton_tolerance")) real = &o->newton_tolerance;
  else if (!strcmp(key,"filter_b")) real = &o->filter_b;
  else if (!strcmp(key,"filter_k")) real = &o->filter_k;
  return real && double_parse(value,real) && *real >= 0;
}

static int write_csv(const char *path, const ch_model *m, const group *groups, size_t ngroups,
                     size_t cells, const double *snapshots, const size_t *steps, size_t ns, int with_step) {
  FILE *f = fopen(path,"w");
  if (!f) return 0;
  if (with_step) fprintf(f,"step,cell,");
  fprintf(f,"temperature,pressure");
  size_t nv = ch_model_species_count(m);
  for (size_t i = 0; i < nv; ++i) fprintf(f,",%s",ch_model_species_name(m,i));
  fputc('\n',f);
  for (size_t k = with_step ? 0 : ns-1; k < ns; ++k) {
    for (size_t g = 0; g < ngroups; ++g) {
      size_t nc = ch_state_cells(groups[g].state);
      ch_conditions *conditions = ch_state_conditions(groups[g].state);
      for (size_t c = 0; c < nc; ++c) {
        size_t global = groups[g].first+c;
        if (with_step) fprintf(f,"%zu,%zu,",steps[k],global);
        fprintf(f,"%.17g,%.17g",conditions[c].temperature,conditions[c].pressure);
        for (size_t i = 0; i < nv; ++i) fprintf(f,",%.17g",snapshots[(k*cells+global)*nv+i]);
        fputc('\n',f);
      }
    }
  }
  int failed = ferror(f);
  return !fclose(f) && !failed;
}

int main(int argc, char **argv) {
  size_t threads = (size_t)omp_get_max_threads(), group_size = CH_DEFAULT_GROUP_SIZE;
  const char *path = NULL, *method = "rosenbrock3";
  double rtol = -1, atol = -1;
  int legacy = 1;
  ch_integrator_options o = ch_integrator_default(3);
  const char *env = getenv("CH_GROUP_SIZE");
  if (env && (!size_parse(env,&group_size) || !group_size)) goto usage;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i],"--threads")) {
      if (++i == argc || !size_parse(argv[i],&threads) || !threads || threads > INT_MAX) goto usage;
    } else if (!strcmp(argv[i],"--method")) {
      if (++i == argc) goto usage;
      method = argv[i];
    } else if (!strcmp(argv[i],"--option")) {
      if (++i == argc || !option(argv[i],&o,&group_size,&rtol,&atol)) goto usage;
    } else if (argv[i][0] == '-' || path) goto usage;
    else path = argv[i];
  }
  o.method = method_id(method,&legacy);
  if (!path || !o.method || !threads || !group_size) goto usage;
  scenario sc = {0}; sunlight *sun = NULL; size_t nsun = 0, ngroups = 0;
  group *groups = NULL; ch_model *m = NULL;
  double *snapshots = NULL; size_t *snapshot_steps = NULL;
  int status = CH_INVALID;
  if (!read_scenario(path,&sc) || !read_sunlight(sc.photolysis,&sun,&nsun)) goto cleanup;
  size_t cells = count_cells(sc.cells);
  if (!cells) goto cleanup;
  if (group_size > cells) group_size = cells;
  ngroups = cells/group_size+(cells%group_size != 0);
  if (threads > ngroups) threads = ngroups;
  omp_set_num_threads((int)threads);
  if (!getenv("OMP_SCHEDULE")) omp_set_schedule(omp_sched_dynamic,1);
  m = ch_build_ts1(); groups = calloc(ngroups,sizeof(*groups));
  if (!m || !groups) { status = CH_MEMORY; goto cleanup; }
  size_t nv = ch_model_species_count(m), np = ch_model_parameter_count(m);
  size_t ns = sc.every ? sc.steps/sc.every+(sc.steps%sc.every != 0) : 1;
  if (ns > SIZE_MAX/cells/nv/sizeof(double)) { status = CH_MEMORY; goto cleanup; }
  snapshots = malloc(ns*cells*nv*sizeof(double)); snapshot_steps = malloc(ns*sizeof(size_t));
  if (!snapshots || !snapshot_steps) { status = CH_MEMORY; goto cleanup; }
  for (size_t i = 0; i < ns; ++i) snapshot_steps[i] = i+1 == ns ? sc.steps : (i+1)*sc.every;
  for (size_t g = 0; g < ngroups; ++g) {
    size_t nc = cells-g*group_size;
    if (nc > group_size) nc = group_size;
    groups[g].first = g*group_size;
    ch_state *s = groups[g].state = ch_state_create(m,nc);
    if (!s) { status = CH_MEMORY; goto cleanup; }
    status = ch_init_ts1(s); if (status) goto cleanup;
    if (!legacy || rtol >= 0 || atol >= 0) {
      double r = rtol >= 0 ? rtol : legacy ? ch_state_relative_tolerance(s) : 1e-4;
      double a = atol >= 0 ? atol : 1e-18;
      status = ch_state_tolerances(s,r,&a,1); if (status) goto cleanup;
    }
    if (!legacy) {
      groups[g].integrator = ch_integrator_create(s,&o);
      if (!groups[g].integrator) { status = CH_INVALID; goto cleanup; }
    }
    groups[g].base_parameters = malloc(np*nc*sizeof(double));
    if (!groups[g].base_parameters) { status = CH_MEMORY; goto cleanup; }
  }
  status = initialize_csv(sc.cells,m,groups,cells,group_size);
  if (status) goto cleanup;
  for (size_t g = 0; g < ngroups; ++g)
    memcpy(groups[g].base_parameters,ch_state_parameters(groups[g].state),np*ch_state_cells(groups[g].state)*sizeof(double));
  ch_rosenbrock_parameters p = ch_rosenbrock_default(o.method);
  /* Legacy methods keep all defaults unless the option is explicitly given. */
  for (int i = 1; i < argc; ++i) if (!strcmp(argv[i],"--option") && i+1 < argc) {
    const char *v = argv[++i];
    if (!strncmp(v,"h_start=",8)) p.h_start = o.h_start;
    else if (!strncmp(v,"h_max=",6)) p.h_max = o.h_max;
    else if (!strncmp(v,"max_steps=",10)) p.max_steps = o.max_steps;
  }
  double elapsed = 0, minimum = INFINITY; size_t next_snapshot = 0;
  for (size_t step = 0; step < sc.steps; ++step) {
    double start = omp_get_wtime();
    double factor = sun_factor(sun,nsun,(double)step*sc.dt);
    #pragma omp parallel for schedule(runtime)
    for (size_t g = 0; g < ngroups; ++g) {
      ch_state *s = groups[g].state; size_t nc = ch_state_cells(s);
      double *params = ch_state_parameters(s);
      for (size_t j = 0; j < np; ++j) if (ch_model_parameter_name(m,j)[0] == 'j')
        for (size_t c = 0; c < nc; ++c) params[j*nc+c] = groups[g].base_parameters[j*nc+c]*factor;
      ch_stats st = {0}; ch_integrator_stats extra = {0};
      groups[g].status = ch_update_rates(s);
      if (!groups[g].status) groups[g].status = legacy ? ch_solve(s,sc.dt,&p,&st) :
        ch_integrator_step(groups[g].integrator,sc.dt,&st,&extra);
      ch_stats_add(&groups[g].stats,&st); ch_integrator_stats_add(&groups[g].extra,&extra);
    }
    elapsed += omp_get_wtime()-start;
    /* Diagnostics and snapshot copies are explicitly outside the timer. */
    for (size_t g = 0; g < ngroups; ++g) {
      if ((status = groups[g].status)) {
        fprintf(stderr,"group=%zu step=%zu status=%s\n",g,step+1,ch_status_string(status));
        goto report;
      }
      ch_state *s = groups[g].state; size_t nc = ch_state_cells(s); const double *y = ch_state_values(s);
      for (size_t i = 0; i < nv; ++i) for (size_t c = 0; c < nc; ++c) {
        double value = y[i*nc+c];
        if (!isfinite(value)) { status = CH_NAN; goto report; }
        minimum = fmin(minimum,value);
        if (step+1 == snapshot_steps[next_snapshot])
          snapshots[(next_snapshot*cells+groups[g].first+c)*nv+i] = value;
      }
    }
    if (step+1 == snapshot_steps[next_snapshot]) ++next_snapshot;
  }
  if (!write_csv("ts1_final.csv",m,groups,ngroups,cells,snapshots,snapshot_steps,ns,0) ||
      !write_csv("ts1_snapshots.csv",m,groups,ngroups,cells,snapshots,snapshot_steps,ns,1)) status = CH_INVALID;
report:;
  ch_stats total = {0}; ch_integrator_stats extra = {0};
  for (size_t g = 0; g < ngroups; ++g) {
    ch_stats_add(&total,&groups[g].stats); ch_integrator_stats_add(&extra,&groups[g].extra);
  }
  printf("elapsed_ms=%.3f\ninternal_steps=%zu accepted=%zu rejected=%zu function_calls=%zu jacobian_updates=%zu decompositions=%zu solves=%zu\n",
    elapsed*1000,total.internal_steps,total.accepted,total.rejected,total.function_calls,total.jacobian_updates,total.decompositions,total.solves);
  printf("jacobian_reuses=%zu factor_reuses=%zu positivity_rejections=%zu nonfinite_rejections=%zu clipped_values=%zu nonlinear_iterations=%zu nonlinear_failures=%zu minimum_concentration=%.17g\n",
    extra.jacobian_reuses,extra.factor_reuses,extra.positivity_rejections,extra.nonfinite_rejections,extra.clipped_values,extra.nonlinear_iterations,extra.nonlinear_failures,minimum);
cleanup:
  if (groups) for (size_t g = 0; g < ngroups; ++g) {
    ch_integrator_free(groups[g].integrator); ch_state_free(groups[g].state); free(groups[g].base_parameters);
  }
  free(groups); ch_model_free(m); free(sc.cells); free(sc.photolysis); free(sun); free(snapshots); free(snapshot_steps);
  if (status) fprintf(stderr,"TS1 scenario: %s\n",ch_status_string(status));
  return status ? 1 : 0;
usage:
  fprintf(stderr,"usage: ts1_scenario [--threads N] [--method NAME] [--option KEY=VALUE ...] SCENARIO\n");
  return 2;
}
