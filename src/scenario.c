#include "chem.h"
#include "chem_adaptive.h"
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
ChemModel *chem_load_ts1(int specialized);

typedef struct { char *cells, *light; double dt; size_t steps, snapshots; } Scenario;
typedef struct { double *time, *factor; size_t count; } Light;
static char *trim(char *s) {
  while (isspace((unsigned char)*s)) ++s;
  char *end = s+strlen(s);
  while (end > s && isspace((unsigned char)end[-1])) --end;
  *end = '\0'; return s;
}
static int number(const char *s, double *value) {
  char *end; errno = 0; *value = strtod(s,&end);
  return *s && end != s && !*end && !errno && isfinite(*value);
}
static char *relative_path(const char *path, const char *name) {
  size_t dir = 0;
  const char *slash = strrchr(path,'/');
  if (slash && *name != '/') dir = (size_t)(slash-path)+1;
  size_t n = strlen(name);
  if (dir > SIZE_MAX-n-1) return NULL;
  char *result = malloc(dir+n+1);
  if (result) { memcpy(result,path,dir); memcpy(result+dir,name,n+1); }
  return result;
}
static int scenario_read(const char *path, Scenario *s) {
  FILE *file = fopen(path,"r"); if (!file) { perror(path); return 0; }
  char *line; int ok = 1, dt_set = 0, steps_set = 0;
  while ((line = line_read(file))) {
    char *comment = strchr(line,'#'); if (comment) *comment = '\0';
    char *key = trim(line);
    if (!*key) { free(line); continue; }
    char *eq = strchr(key,'=');
    if (!eq) { free(line); ok = 0; break; }
    *eq = '\0'; char *value = trim(eq+1); key = trim(key);
    if (!strcmp(key,"cells") && !s->cells) s->cells = relative_path(path,value);
    else if (!strcmp(key,"photolysis") && !s->light) s->light = relative_path(path,value);
    else if (!strcmp(key,"dt") && !dt_set) { ok = number(value,&s->dt) && s->dt > 0.0; dt_set = 1; }
    else if (!strcmp(key,"steps") && !steps_set) { ok = positive_size(value,&s->steps,0); steps_set = 1; }
    else if (!strcmp(key,"snapshot_every") && !s->snapshots) ok = positive_size(value,&s->snapshots,0);
    else ok = 0;
    free(line); if (!ok) break;
  }
  if (ferror(file)) ok = 0;
  fclose(file);
  return ok && s->cells && s->light && dt_set && steps_set && isfinite(s->dt*s->steps);
}
static size_t row_count(const char *path) {
  FILE *f = fopen(path,"r"); if (!f) return 0;
  char *line = line_read(f); free(line); size_t count = 0;
  while ((line = line_read(f))) { if (*trim(line)) ++count; free(line); }
  if (ferror(f)) count = 0;
  fclose(f); return count;
}
static int light_read(const char *path, Light *l) {
  Table t = {0}; size_t count = row_count(path);
  if (!count || !table_read(path,&t,count)) { table_free(&t); return 0; }
  int ok = t.columns == 2 && !strcmp(t.names[0],"time_s") && !strcmp(t.names[1],"factor");
  l->time = malloc(count*sizeof(double)); l->factor = malloc(count*sizeof(double));
  if (!l->time || !l->factor) ok = 0;
  if (ok) for (size_t i = 0; i < count; ++i) {
    l->time[i] = t.values[2*i]; l->factor[i] = t.values[2*i+1];
    if (l->factor[i] < 0.0 || (i && l->time[i] <= l->time[i-1])) ok = 0;
  }
  l->count = count; table_free(&t); return ok;
}
static double light_at(const Light *l, double time) {
  if (time <= l->time[0]) return l->factor[0];
  for (size_t i = 1; i < l->count; ++i) if (time < l->time[i]) {
    double weight = (time-l->time[i-1])/(l->time[i]-l->time[i-1]);
    return l->factor[i-1]+weight*(l->factor[i]-l->factor[i-1]);
  }
  return l->factor[l->count-1];
}
static void stats_add(ChemAdaptiveStats *a, const ChemAdaptiveStats *b) {
  a->base.internal_steps += b->base.internal_steps; a->base.accepted += b->base.accepted;
  a->base.rejected += b->base.rejected; a->base.function_calls += b->base.function_calls;
  a->base.jacobian_updates += b->base.jacobian_updates; a->base.decompositions += b->base.decompositions;
  a->base.solves += b->base.solves; a->nonlinear_iterations += b->nonlinear_iterations;
  a->matrix_reuses += b->matrix_reuses; a->positivity_rejections += b->positivity_rejections;
  a->failed_iterations += b->failed_iterations;
}
static int header(FILE *f, const ChemModel *m, int snapshot) {
  fprintf(f,"%stemperature,pressure",snapshot ? "step,cell," : "");
  for (size_t v = 0; v < m->nvariables; ++v) fprintf(f,",%s",m->species[m->variable_species[v]].name);
  fputc('\n',f); return !ferror(f);
}
static int output(FILE *f, ChemState **states, size_t cells, size_t group, size_t step) {
  for (size_t c = 0; c < cells; ++c) {
    ChemState *s = states[c/group]; size_t k = c%group;
    if (step) fprintf(f,"%zu,%zu,",step,c);
    fprintf(f,"%.17g,%.17g",s->conditions[k].temperature,s->conditions[k].pressure);
    for (size_t v = 0; v < s->model->nvariables; ++v) fprintf(f,",%.17g",s->variables[v*s->stride+k]);
    fputc('\n',f);
  }
  return !ferror(f);
}
static int option(const char *text, ChemAdaptive *p, double *rtol, double *atol, size_t *group) {
  const char *eq = strchr(text,'='); if (!eq || (size_t)(eq-text) > 63) return 0;
  char key[64]; memcpy(key,text,(size_t)(eq-text)); key[eq-text] = '\0';
  double x; if (!number(eq+1,&x)) return 0;
  if (!strcmp(key,"rtol")) { if (!(x > 0.0)) return 0; *rtol = x; }
  else if (!strcmp(key,"atol")) { if (!(x > 0.0)) return 0; *atol = x; }
  else if (!strcmp(key,"h_start")) { if (!(x > 0.0)) return 0; p->h_start = x; }
  else if (!strcmp(key,"h_max")) { if (x < 0.0) return 0; p->h_max = x; }
  else if (!strcmp(key,"safety")) { if (!(x > 0.0 && x < 1.0)) return 0; p->safety = x; }
  else if (!strcmp(key,"growth")) { if (x < 1.0) return 0; p->growth = x; }
  else if (!strcmp(key,"positivity")) { if (x < 0.0) return 0; p->positivity = x; }
  else if (!strcmp(key,"reuse_ratio")) { if (!(x >= 0.0 && x < 1.0)) return 0; p->reuse_ratio = x; }
  else {
    size_t value; if (!positive_size(eq+1,&value,1)) return 0;
    if (!strcmp(key,"controller") && value <= 1) p->controller = (int)value;
    else if (!strcmp(key,"norm") && value <= 2) p->norm = (int)value;
    else if (!strcmp(key,"order") && value >= 1 && value <= 5) p->max_order = (int)value;
    else if (!strcmp(key,"reuse") && value <= 1) p->reuse = (int)value;
    else if (!strcmp(key,"jacobian_age") && value) p->jacobian_age = value;
    else if (!strcmp(key,"max_steps") && value) p->max_steps = value;
    else if (!strcmp(key,"group") && value) *group = value;
    else return 0;
  }
  return 1;
}
int main(int argc, char **argv) {
  size_t threads = (size_t)omp_get_max_threads(), group = CHEM_GROUP_SIZE;
  const char *method = "micm3", *path = NULL;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i],"--method") && i+1 < argc) method = argv[++i];
  }
  int legacy = !strncmp(method,"micm",4), id = 0;
  if (legacy && strlen(method) == 5 && method[4] >= '2' && method[4] <= '6') id = method[4]-'0';
  else if (!strcmp(method,"ros2")) id = 2;
  else if (!strcmp(method,"ros3")) id = 3;
  else if (!strcmp(method,"ros4")) id = 4;
  else if (!strcmp(method,"rodas3")) id = 5;
  else if (!strcmp(method,"rodas4")) id = 6;
  else if (!strcmp(method,"bdf")) id = 10;
  if (!id) goto usage;
  ChemAdaptive p = chem_adaptive_parameters(id);
  double rtol = legacy ? 0.0 : 1e-4, atol = legacy ? 0.0 : 1e-18;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i],"--threads") && i+1 < argc) { if (!positive_size(argv[++i],&threads,0) || threads > INT_MAX) goto usage; }
    else if (!strcmp(argv[i],"--method") && i+1 < argc) ++i;
    else if (!strcmp(argv[i],"--option") && i+1 < argc) { if (!option(argv[++i],&p,&rtol,&atol,&group)) goto usage; }
    else if (!strcmp(argv[i],"--group") && i+1 < argc) { if (!positive_size(argv[++i],&group,0)) goto usage; }
    else if (argv[i][0] == '-' || path) goto usage;
    else path = argv[i];
  }
  if (!path || !group || threads > INT_MAX) goto usage;
  omp_set_num_threads((int)threads);
  if (!getenv("OMP_SCHEDULE")) omp_set_schedule(omp_sched_static,0);
  Scenario scenario = {0}; Light light = {0}; Table table = {0};
  ChemModel *m = NULL; ChemState **states = NULL; ChemWorkspace **work = NULL;
  ChemStatus *statuses = NULL; ChemAdaptiveStats *totals = NULL;
  double **base_rates = NULL, **base_custom = NULL;
  FILE *snapshots = NULL; size_t groups = 0; int failed = 1;
  if (!scenario_read(path,&scenario) || !light_read(scenario.light,&light)) { fprintf(stderr,"invalid scenario or photolysis file\n"); goto cleanup; }
  size_t cells = row_count(scenario.cells);
  if (!cells || !table_read(scenario.cells,&table,cells)) goto cleanup;
  m = chem_load_ts1(1); if (!m) goto cleanup;
  for (size_t col = 0; col < table.columns; ++col) {
    const char *name = table.names[col];
    if (strcmp(name,"temperature") && strcmp(name,"pressure") && chem_variable_id(m,name) == SIZE_MAX && chem_custom_id(m,name) == SIZE_MAX) {
      fprintf(stderr,"unknown CSV column: %s\n",name); goto cleanup;
    }
  }
  groups = 1+(cells-1)/group;
  states = calloc(groups,sizeof *states); work = calloc(groups,sizeof *work);
  statuses = calloc(groups,sizeof *statuses); totals = calloc(groups,sizeof *totals);
  base_rates = calloc(groups,sizeof *base_rates); base_custom = calloc(groups,sizeof *base_custom);
  if (!states || !work || !statuses || !totals || !base_rates || !base_custom) goto cleanup;
  for (size_t g = 0; g < groups; ++g) {
    size_t first = g*group, count = cells-first < group ? cells-first : group;
    ChemState *s = states[g] = chem_state_create(m,count); if (!s) goto cleanup;
    for (size_t col = 0; col < table.columns; ++col) {
      const char *name = table.names[col]; size_t v = chem_variable_id(m,name), u = chem_custom_id(m,name);
      for (size_t c = 0; c < count; ++c) {
        double x = table.values[(first+c)*table.columns+col];
        if (!strcmp(name,"temperature")) s->conditions[c].temperature = x;
        else if (!strcmp(name,"pressure")) s->conditions[c].pressure = x;
        else if (v != SIZE_MAX) s->variables[v*count+c] = x;
        else s->custom[u*count+c] = x;
      }
    }
    for (size_t c = 0; c < count; ++c) chem_ideal_air_density(&s->conditions[c]);
    if (rtol) s->relative_tolerance = rtol;
    if (atol) for (size_t i = 0; i < m->nvariables; ++i) s->absolute_tolerance[i] = atol;
    chem_calculate_rates(s);
    base_rates[g] = malloc(m->nreactions*count*sizeof(double)); base_custom[g] = malloc(m->ncustom*count*sizeof(double));
    if (!base_rates[g] || !base_custom[g]) goto cleanup;
    memcpy(base_rates[g],s->rates,m->nreactions*count*sizeof(double)); memcpy(base_custom[g],s->custom,m->ncustom*count*sizeof(double));
    if (!legacy) { work[g] = chem_workspace_create(s); if (!work[g]) goto cleanup; }
  }
  snapshots = fopen("ts1_snapshots.csv","w");
  if (!snapshots || !header(snapshots,m,1)) goto cleanup;
  double elapsed = 0.0, previous_light = -1.0;
  failed = 0;
  for (size_t step = 0; step < scenario.steps; ++step) {
    double start = omp_get_wtime();
    double factor = light_at(&light,step*scenario.dt);
    #pragma omp parallel for schedule(runtime)
    for (size_t g = 0; g < groups; ++g) {
      ChemState *s = states[g]; size_t nc = s->cells;
      for (size_t u = 0; u < m->ncustom; ++u) if (m->custom_names[u][0] == 'j')
        for (size_t c = 0; c < nc; ++c) s->custom[u*nc+c] = base_custom[g][u*nc+c]*factor;
      for (size_t r = 0; r < m->nreactions; ++r) {
        const ChemRate *rate = &m->reactions[r].rate;
        if (rate->type == CHEM_USER && rate->label[0] == 'j')
          for (size_t c = 0; c < nc; ++c) s->rates[r*nc+c] = base_rates[g][r*nc+c]*factor;
      }
      ChemAdaptiveStats one = {0};
      if (legacy) {
        ChemSolver solver = chem_solver(id);
        statuses[g] = chem_solve(&solver,scenario.dt,s,&one.base);
      } else {
        if (factor != previous_light) chem_workspace_reset(work[g]);
        statuses[g] = chem_adaptive_solve(&p,scenario.dt,s,work[g],&one);
      }
      stats_add(totals+g,&one);
      if (!statuses[g]) for (size_t a = 0; a < m->nvariables*nc; ++a)
        if (!isfinite(s->variables[a]) || s->variables[a] < -1e-20) { statuses[g] = CHEM_NAN; break; }
    }
    elapsed += (omp_get_wtime()-start)*1000.0; previous_light = factor;
    for (size_t g = 0; g < groups; ++g) if (statuses[g]) {
      fprintf(stderr,"step=%zu group=%zu status=%s\n",step+1,g,chem_status_string(statuses[g])); failed = 1; break;
    }
    if (failed) break;
    if (scenario.snapshots && ((step+1)%scenario.snapshots == 0 || step+1 == scenario.steps))
      if (!output(snapshots,states,cells,group,step+1)) { failed = 1; break; }
  }
  ChemAdaptiveStats total = {0}; for (size_t g = 0; g < groups; ++g) stats_add(&total,totals+g);
  ChemStats *s = &total.base;
  printf("elapsed_ms=%.3f\ninternal_steps=%llu accepted=%llu rejected=%llu function_calls=%llu jacobian_updates=%llu decompositions=%llu solves=%llu\n",
    elapsed,(unsigned long long)s->internal_steps,(unsigned long long)s->accepted,(unsigned long long)s->rejected,
    (unsigned long long)s->function_calls,(unsigned long long)s->jacobian_updates,(unsigned long long)s->decompositions,(unsigned long long)s->solves);
  printf("nonlinear_iterations=%llu matrix_reuses=%llu positivity_rejections=%llu failed_iterations=%llu\n",
    (unsigned long long)total.nonlinear_iterations,(unsigned long long)total.matrix_reuses,
    (unsigned long long)total.positivity_rejections,(unsigned long long)total.failed_iterations);
  if (!failed) {
    FILE *final = fopen("ts1_final.csv","w");
    if (!final) failed = 1;
    else { if (!header(final,m,0) || !output(final,states,cells,group,0)) failed = 1; if (fclose(final)) failed = 1; }
  }
cleanup:
  if (snapshots && fclose(snapshots)) failed = 1;
  for (size_t g = 0; g < groups; ++g) {
    if (work) chem_workspace_free(work[g]); if (states) chem_state_free(states[g]);
    if (base_rates) free(base_rates[g]); if (base_custom) free(base_custom[g]);
  }
  free(states); free(work); free(statuses); free(totals); free(base_rates); free(base_custom);
  chem_model_free(m); table_free(&table); free(light.time); free(light.factor); free(scenario.cells); free(scenario.light);
  return failed ? 1 : 0;
usage:
  fprintf(stderr,"usage: %s [--threads N] [--method micm2..micm6|ros2|ros3|ros4|rodas3|rodas4|bdf] [--option KEY=VALUE ...] SCENARIO\n",argv[0]);
  return 2;
}
