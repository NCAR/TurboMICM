#include "chem.h"
#include "chem_integrator.h"
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef CK_DEFAULT_GROUP_SIZE
#define CK_DEFAULT_GROUP_SIZE 32
#endif
extern const ck_mechanism ts1;
#include "bench_io.h"

typedef struct { double time, factor; } photo_point;
typedef struct
{
  char *cells, *photo;
  double dt;
  size_t steps, snapshot;
  photo_point *points;
  size_t npoints;
} scenario;
static char *trim(char *s)
{
  while (isspace((unsigned char)*s)) ++s;
  char *end = s + strlen(s);
  while (end > s && isspace((unsigned char)end[-1])) *--end = 0;
  return s;
}
static int number(const char *s, double *v)
{
  char *end;
  errno = 0;
  *v = strtod(s, &end);
  return end != s && !*end && !errno && isfinite(*v);
}
static char *relative_path(const char *scenario_path, const char *value)
{
  const char *slash = strrchr(scenario_path, '/');
  size_t prefix = value[0] == '/' || !slash ? 0 : (size_t)(slash-scenario_path)+1;
  char *path = malloc(prefix+strlen(value)+1);
  if (path) { memcpy(path, scenario_path, prefix); strcpy(path+prefix, value); }
  return path;
}
static int read_scenario(const char *path, scenario *s)
{
  FILE *f = fopen(path, "r");
  if (!f) return 0;
  char *line;
  unsigned seen = 0;
  int ok = 1;
  while ((line = read_line(f)))
  {
    char *hash = strchr(line, '#');
    if (hash) *hash = 0;
    char *key = trim(line);
    if (*key)
    {
      char *eq = strchr(key, '=');
      if (!eq) { free(line); ok = 0; break; }
      *eq = 0;
      char *value = trim(eq+1);
      key = trim(key);
      unsigned bit = 0;
      if (!strcmp(key, "cells")) { bit=1; s->cells=relative_path(path,value); ok=s->cells && *value; }
      else if (!strcmp(key, "photolysis")) { bit=2; s->photo=relative_path(path,value); ok=s->photo && *value; }
      else if (!strcmp(key, "dt")) { bit=4; ok=number(value,&s->dt) && s->dt>0; }
      else if (!strcmp(key, "steps")) { bit=8; ok=size_argument(value,&s->steps) && s->steps>0; }
      else if (!strcmp(key, "snapshot_every")) { bit=16; ok=size_argument(value,&s->snapshot) && s->snapshot>0; }
      else ok=0;
      if (seen & bit) ok=0;
      seen |= bit;
    }
    free(line);
    if (!ok) break;
  }
  if (ferror(f)) ok=0;
  fclose(f);
  if (!ok || (seen & 15) != 15 || !isfinite(s->dt*s->steps)) return 0;
  f = fopen(s->photo,"r");
  if (!f) return 0;
  line = read_line(f);
  ok = line && !strcmp(trim(line),"time_s,factor");
  free(line);
  while (ok && (line=read_line(f)))
  {
    char *comma=strchr(line,',');
    photo_point p;
    if (!comma) ok=0;
    else
    {
      *comma=0;
      ok=number(trim(line),&p.time) && number(trim(comma+1),&p.factor) && p.factor>=0 &&
         (!s->npoints || p.time>s->points[s->npoints-1].time);
      if (ok)
      {
        photo_point *next=realloc(s->points,(s->npoints+1)*sizeof *next);
        if (!next) ok=0;
        else { s->points=next; s->points[s->npoints++]=p; }
      }
    }
    free(line);
  }
  if (ferror(f)) ok=0;
  fclose(f);
  return ok && s->npoints;
}
static double photo_factor(const scenario *s, double time)
{
  if (time <= s->points[0].time) return s->points[0].factor;
  size_t hi=1;
  while (hi<s->npoints && s->points[hi].time<time) ++hi;
  if (hi==s->npoints) return s->points[hi-1].factor;
  const photo_point *a=&s->points[hi-1], *b=&s->points[hi];
  return a->factor+(b->factor-a->factor)*((time-a->time)/(b->time-a->time));
}
static size_t cell_count(const char *path)
{
  FILE *f=fopen(path,"r");
  if (!f) return 0;
  char *line=read_line(f);
  free(line);
  size_t n=0;
  while ((line=read_line(f))) { if (*trim(line)) ++n; free(line); }
  if (ferror(f)) n=0;
  fclose(f);
  return n;
}
static int write_snapshots(const ck_model *m, const group *groups, size_t ngroups,
                           const double *snapshots, const size_t *steps, size_t nsnap, size_t cells)
{
  FILE *f=fopen("ts1_snapshots.csv","w");
  if (!f) return 0;
  fputs("step,cell,temperature,pressure",f);
  size_t n=ck_model_nspecies(m);
  for (size_t j=0;j<n;++j) fprintf(f,",%s",ck_species_name(m,j));
  fputc('\n',f);
  for (size_t k=0;k<nsnap;++k)
    for (size_t g=0;g<ngroups;++g)
    {
      const ck_state *s=groups[g].state;
      for (size_t c=0;c<s->cells;++c)
      {
        size_t cell=groups[g].first+c;
        fprintf(f,"%zu,%zu,%.17g,%.17g",steps[k],cell,s->conditions[c].temperature,s->conditions[c].pressure);
        for (size_t j=0;j<n;++j) fprintf(f,",%.17g",snapshots[(k*n+j)*cells+cell]);
        fputc('\n',f);
      }
    }
  int ok=!ferror(f);
  if (fclose(f)) ok=0;
  return ok;
}
static int method_id(const char *name, int *legacy)
{
  *legacy=0;
  if (!strcmp(name,"micm3")) { *legacy=1; return 3; }
  if (!strcmp(name,"micm2")) { *legacy=1; return 2; }
  if (!strcmp(name,"micm4")) { *legacy=1; return 4; }
  if (!strcmp(name,"micm43")) { *legacy=1; return 43; }
  if (!strcmp(name,"micm6")) { *legacy=1; return 6; }
  if (!strcmp(name,"ros2")) return 2;
  if (!strcmp(name,"ros3")) return 3;
  if (!strcmp(name,"ros4")) return 4;
  if (!strcmp(name,"rodas3")) return 43;
  if (!strcmp(name,"rodas4")) return 6;
  if (!strcmp(name,"extrap2")) return 12;
  if (!strcmp(name,"extrap3")) return 13;
  return 0;
}
int main(int argc,char **argv)
{
  scenario sc={0};
  const char *path=NULL, *method="micm3";
  ck_integrator_options options=ck_integrator_defaults();
  size_t threads=(size_t)omp_get_max_threads(), width=CK_DEFAULT_GROUP_SIZE;
  double rtol=-1, atol=-1;
  int legacy=1;
  const char *env=getenv("CK_GROUP_SIZE");
  int explicit_width=env!=NULL;
  if (env && (!size_argument(env,&width) || !width)) goto usage;
  for (int a=1;a<argc;++a)
  {
    if (!strcmp(argv[a],"--threads"))
    {
      if (++a==argc || !size_argument(argv[a],&threads) || !threads || threads>INT_MAX) goto usage;
    }
    else if (!strcmp(argv[a],"--method")) { if (++a==argc) goto usage; method=argv[a]; }
    else if (!strcmp(argv[a],"--option"))
    {
      if (++a==argc) goto usage;
      char key[80];
      const char *eq=strchr(argv[a],'=');
      double v;
      if (!eq || (size_t)(eq-argv[a])>=sizeof key || !number(eq+1,&v)) goto usage;
      memcpy(key,argv[a],(size_t)(eq-argv[a])); key[eq-argv[a]]=0;
      if (!strcmp(key,"rtol") && v>=0) rtol=v;
      else if (!strcmp(key,"atol") && v>0) atol=v;
      else if (!strcmp(key,"group_size") && size_argument(eq+1,&width) && width) explicit_width=1;
      else if (!strcmp(key,"jacobian_steps") && size_argument(eq+1,&options.jacobian_steps) && options.jacobian_steps) {}
      else if (!strcmp(key,"max_steps") && size_argument(eq+1,&options.max_steps) && options.max_steps) {}
      else if (!strcmp(key,"start") && v>0) options.start=v;
      else if (!strcmp(key,"safety") && v>0 && v<=1) options.safety=v;
      else if (!strcmp(key,"growth") && v>=1) options.growth=v;
      else if (!strcmp(key,"shrink") && v>0 && v<1) options.shrink=v;
      else if (!strcmp(key,"deadband") && v>=0 && v<1) options.deadband=v;
      else if (!strcmp(key,"warm_start") && (v==0 || v==1)) options.warm_start=(int)v;
      else if (!strcmp(key,"positivity") && (v==0 || v==1)) options.positivity=(int)v;
      else if (!strcmp(key,"reuse_rhs") && (v==0 || v==1)) options.reuse_rhs=(int)v;
      else if (!strcmp(key,"controller") && (v==0 || v==1)) options.controller=(int)v;
      else if (!strcmp(key,"norm") && (v==0 || v==1 || v==2)) options.norm=(int)v;
      else goto usage;
    }
    else if (!path && argv[a][0]!='-') path=argv[a];
    else goto usage;
  }
  options.method=method_id(method,&legacy);
  if (!path || !options.method || !width) goto usage;
  if (!read_scenario(path,&sc)) { fprintf(stderr,"Invalid scenario: %s\n",path); return 2; }
  size_t cells=cell_count(sc.cells);
  if (!cells) { fprintf(stderr,"Empty or unreadable cells\n"); return 2; }
  if (!explicit_width && width>(cells-1)/threads+1) width=(cells-1)/threads+1;
  if (width>cells) width=cells;
  size_t ngroups=(cells-1)/width+1;
  if (threads>ngroups) threads=ngroups;
  omp_set_num_threads((int)threads);
  if (!getenv("OMP_SCHEDULE")) omp_set_schedule(omp_sched_static,0);
  ck_model *model=NULL;
  group *groups=calloc(ngroups,sizeof *groups);
  ck_integrator **integrators=calloc(ngroups,sizeof *integrators);
  double **base=calloc(ngroups,sizeof *base);
  size_t *snapsteps=NULL;
  double *snapshots=NULL;
  int exit_code=1;
  ck_rosenbrock_parameters parameters;
  ck_rosenbrock_parameters_init(options.method,&parameters);
  if (!groups || !integrators || !base || ck_model_create(&ts1,&model)) goto done;
  size_t n=ck_model_nspecies(model), nc=ck_model_ncustom(model);
  size_t nsnap=sc.snapshot ? (sc.steps-1)/sc.snapshot+1 : 0;
  if (nsnap)
  {
    if (nsnap>SIZE_MAX/n/cells/sizeof(double)) goto done;
    snapshots=calloc(nsnap*n*cells,sizeof(double));
    snapsteps=calloc(nsnap,sizeof(size_t));
    if (!snapshots || !snapsteps) goto done;
  }
  for (size_t g=0;g<ngroups;++g)
  {
    groups[g].first=g*width;
    size_t count=cells-g*width;
    if (count>width) count=width;
    if (ck_state_create(model,count,&groups[g].state)) goto done;
  }
  if (!read_cells(sc.cells,model,groups,ngroups,cells,width)) goto done;
  for (size_t g=0;g<ngroups;++g)
  {
    ck_state *s=groups[g].state;
    if (rtol>=0) s->relative_tolerance=rtol;
    else if (!legacy) s->relative_tolerance=1e-4;
    if (atol>0 || !legacy)
      for (size_t j=0;j<n;++j) s->absolute_tolerance[j]=atol>0 ? atol : 1e-18;
    if (!legacy && ck_integrator_create(s,&options,&integrators[g])) goto done;
    base[g]=malloc(nc*s->cells*sizeof(double));
    if (!base[g]) goto done;
    memcpy(base[g],s->custom,nc*s->cells*sizeof(double));
  }
  double elapsed=0;
  size_t begin=0, snapshot_index=0;
  while (begin<sc.steps)
  {
    size_t end=sc.snapshot && sc.snapshot<sc.steps-begin ? begin+sc.snapshot : sc.steps;
    double t0=omp_get_wtime();
#pragma omp parallel for schedule(runtime)
    for (size_t g=0;g<ngroups;++g)
    {
      ck_state *s=groups[g].state;
      for (size_t step=begin;step<end && !groups[g].status;++step)
      {
        double factor=photo_factor(&sc,step*sc.dt);
        for (size_t j=0;j<nc;++j)
          if (ck_custom_name(model,j)[0]=='j')
            for (size_t c=0;c<s->cells;++c) s->custom[j*s->cells+c]=base[g][j*s->cells+c]*factor;
        groups[g].status=ck_update_rates(s);
        if (groups[g].status) break;
        ck_stats stats;
        groups[g].status=legacy ? ck_solve(s,sc.dt,&parameters,&stats) : ck_integrator_advance(integrators[g],sc.dt,&stats);
        ck_stats_add(&groups[g].stats,&stats);
        for (size_t i=0;i<n*s->cells;++i)
          if (!isfinite(s->y[i]) || s->y[i]<-1e-20) groups[g].status=CK_NAN;
      }
    }
    elapsed+=(omp_get_wtime()-t0)*1000;
    for (size_t g=0;g<ngroups;++g)
      if (groups[g].status)
      { fprintf(stderr,"group=%zu solve_failed=%s\n",g,ck_status_string(groups[g].status)); goto report; }
    if (nsnap)
    {
      snapsteps[snapshot_index]=end;
      for (size_t g=0;g<ngroups;++g)
      {
        ck_state *s=groups[g].state;
        for (size_t j=0;j<n;++j)
          memcpy(snapshots+(snapshot_index*n+j)*cells+groups[g].first,s->y+j*s->cells,s->cells*sizeof(double));
      }
      ++snapshot_index;
    }
    begin=end;
  }
  if (!write_final(model,groups,ngroups) ||
      (nsnap && !write_snapshots(model,groups,ngroups,snapshots,snapsteps,nsnap,cells))) goto done;
  exit_code=0;
report:
  {
    ck_stats st={0};
    for (size_t g=0;g<ngroups;++g) ck_stats_add(&st,&groups[g].stats);
    printf("elapsed_ms=%.3f\n",elapsed);
    printf("internal_steps=%llu accepted=%llu rejected=%llu function_calls=%llu jacobian_updates=%llu decompositions=%llu solves=%llu nonlinear_iterations=0\n",
      (unsigned long long)st.internal_steps,(unsigned long long)st.accepted,(unsigned long long)st.rejected,
      (unsigned long long)st.function_calls,(unsigned long long)st.jacobian_updates,
      (unsigned long long)st.decompositions,(unsigned long long)st.solves);
    printf("cells=%zu groups=%zu group_size=%zu threads=%zu method=%s status=%s\n",cells,ngroups,width,threads,method,exit_code ? "FAIL":"PASS");
  }
done:
  if (integrators) for (size_t g=0;g<ngroups;++g) ck_integrator_free(integrators[g]);
  if (groups) for (size_t g=0;g<ngroups;++g) ck_state_free(groups[g].state);
  if (base) for (size_t g=0;g<ngroups;++g) free(base[g]);
  free(groups); free(integrators); free(base); free(snapshots); free(snapsteps);
  free(sc.cells); free(sc.photo); free(sc.points); ck_model_free(model);
  return exit_code;
usage:
  fprintf(stderr,"usage: %s [--threads N] [--method NAME] [--option KEY=VALUE ...] SCENARIO\n",argv[0]);
  return 2;
}
