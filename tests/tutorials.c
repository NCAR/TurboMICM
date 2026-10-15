/* C ports of the eight tutorial programs and configuration equivalence checks.
 * Original fixtures: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
 * SPDX-License-Identifier: Apache-2.0 */
#include "chem.h"
#include "chem_array.h"
#include <inttypes.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

ChemModel *chem_load_ts1(int);
ChemModel *chem_load_example_no_user(int);
ChemModel *chem_load_example_user(int);
ChemModel *chem_load_example_robertson(int);
ChemModel *chem_load_validation(int);

static void check(int condition, const char *message) {
  if (!condition) { fprintf(stderr,"%s\n",message); exit(EXIT_FAILURE); }
}
static void status(ChemStatus code) { check(code == CHEM_OK,chem_status_string(code)); }
static void close_value(double a, double b) { check(isfinite(a) && isfinite(b) && fabs(a-b) <= 1e-30+1e-10*fabs(b),"configuration differs from code-built mechanism"); }
static ChemModel *new_model(const char *const *names, size_t n) {
  ChemModel *m = chem_model_create(); check(m != NULL,"allocate model");
  for (size_t i = 0; i < n; ++i) status(chem_model_species(m,names[i],0,0,NULL));
  return m;
}
static void reaction(ChemModel *m, const size_t *r, size_t nr, const size_t *p, const double *y, size_t np, ChemRate k) {
  status(chem_model_reaction(m,r,nr,p,y,np,k));
}
static void set(ChemState *s, const char *name, double value) {
  for (size_t c = 0; c < s->cells; ++c) s->variables[chem_variable_id(s->model,name)*s->stride+c] = value;
}
static double get(const ChemState *s, const char *name, size_t cell) { return s->variables[chem_variable_id(s->model,name)*s->stride+cell]; }
static void parameter(ChemState *s, const char *name, double value) {
  size_t i = chem_custom_id(s->model,name); check(i != SIZE_MAX,"unknown rate parameter");
  for (size_t c = 0; c < s->cells; ++c) s->custom[i*s->stride+c] = value;
}
static ChemModel *finish(ChemModel *m) { status(chem_model_reorder(m)); status(chem_model_finalize(m)); return m; }
static ChemModel *robertson(void) {
  const char *names[] = {"A","B","C"}; ChemModel *m = new_model(names,3); ChemRate k = chem_rate(CHEM_USER);
  k.label = "r1"; reaction(m,(size_t[]){0},1,(size_t[]){1},(double[]){1},1,k);
  k.label = "r2"; reaction(m,(size_t[]){1,1},2,(size_t[]){1,2},(double[]){1,1},2,k);
  k.label = "r3"; reaction(m,(size_t[]){1,2},2,(size_t[]){0,2},(double[]){1,1},2,k);
  return finish(m);
}
static ChemModel *rate_example(int user) {
  const char *names[] = {"A","B","C","D","E","F","G"}; ChemModel *m = new_model(names,7);
  double conversion = 1e-6*6.02214076e23;
  ChemRate k = chem_rate(CHEM_ARRHENIUS); k.A = 2.15e-4; k.C = 110;
  reaction(m,(size_t[]){0},1,(size_t[]){1},(double[]){1},1,k);
  k = chem_rate(CHEM_BRANCHED); k.X = 1.2; k.Y = 204.3; k.a0 = .001; k.n = 2;
  reaction(m,(size_t[]){1},1,(size_t[]){2},(double[]){1},1,k);
  k.branch = 1; reaction(m,(size_t[]){1},1,(size_t[]){3},(double[]){1},1,k);
  k = chem_rate(CHEM_SURFACE); k.label = "C"; k.diffusion = 2.3e2; k.molecular_weight = .025; k.probability = .9;
  reaction(m,(size_t[]){2},1,(size_t[]){4},(double[]){1},1,k);
  k = chem_rate(CHEM_TERNARY); k.k0_A = 1.2; k.k0_B = 2.3; k.k0_C = 302.3;
  k.kinf_A = 2.6/conversion; k.kinf_B = -3.1; k.kinf_C = 402.1; k.Fc = .9; k.N = 1.2;
  reaction(m,(size_t[]){3},1,(size_t[]){5},(double[]){2},1,k);
  k = chem_rate(CHEM_TROE); k.k0_A = 1.8e-30*conversion*conversion; k.k0_B = -3;
  k.kinf_A = 2.8e-11*conversion; k.Fc = .6; k.N = 1;
  reaction(m,(size_t[]){4,4},2,(size_t[]){6},(double[]){1},1,k);
  k = chem_rate(CHEM_TUNNELING); k.A = 1.2; k.B = 2.3; k.C = 302.3;
  reaction(m,(size_t[]){5},1,(size_t[]){6},(double[]){1},1,k);
  if (user) {
    k = chem_rate(CHEM_USER); k.label = "my photolysis rate";
    reaction(m,(size_t[]){2},1,(size_t[]){6},(double[]){1},1,k);
    k.label = "my emission rate"; reaction(m,NULL,0,(size_t[]){0},(double[]){1},1,k);
    k.label = "my loss rate"; reaction(m,(size_t[]){1},1,NULL,NULL,0,k);
  }
  return finish(m);
}
static ChemModel *validation_model(void) {
  const char *names[] = {"A","B","C"}; ChemModel *m = new_model(names,3); status(chem_model_species(m,"M",0,1,NULL));
  ChemRate k = chem_rate(CHEM_ARRHENIUS); k.A = .001; k.C = 5;
  reaction(m,(size_t[]){0,3},2,(size_t[]){1},(double[]){1},1,k);
  k = chem_rate(CHEM_TROE); k.k0_A = .0001; k.kinf_A = .002;
  reaction(m,(size_t[]){1,1},2,(size_t[]){2},(double[]){1},1,k);
  k = chem_rate(CHEM_TERNARY); k.k0_A = .002; k.kinf_A = .0001;
  reaction(m,(size_t[]){2},1,(size_t[]){0},(double[]){1},1,k);
  k = chem_rate(CHEM_JPL19); k.k0_A = 1e-5; k.kinf_A = .001; k.kint_A = .00002;
  reaction(m,(size_t[]){0,1},2,(size_t[]){2},(double[]){1},1,k);
  k = chem_rate(CHEM_BRANCHED); k.X = .0001; k.Y = 5; k.a0 = .2; k.n = 2;
  reaction(m,(size_t[]){1},1,(size_t[]){0},(double[]){.8},1,k);
  k.branch = 1; reaction(m,(size_t[]){1},1,(size_t[]){2},(double[]){1.2},1,k);
  k = chem_rate(CHEM_TUNNELING); k.A = .0004; k.B = 2; k.C = 30;
  reaction(m,(size_t[]){2},1,(size_t[]){1},(double[]){1},1,k);
  k = chem_rate(CHEM_TAYLOR); k.A = .0002; k.ncoefficients = 16; k.coefficients[1] = .0001; k.coefficients[15] = 1e-40;
  reaction(m,(size_t[]){0},1,(size_t[]){2},(double[]){1},1,k);
  k = chem_rate(CHEM_REVERSIBLE); k.A = .0001; k.C = 4; k.B = .3;
  reaction(m,(size_t[]){2},1,(size_t[]){0},(double[]){1},1,k);
  k = chem_rate(CHEM_USER); k.label = "photo";
  reaction(m,(size_t[]){1},1,(size_t[]){2},(double[]){1},1,k);
  k = chem_rate(CHEM_SURFACE); k.label = "aerosol"; k.diffusion = 1e-5; k.molecular_weight = .025; k.probability = .7;
  reaction(m,(size_t[]){0},1,(size_t[]){1},(double[]){1},1,k);
  m->species[0].initial = 1; m->species[1].initial = .2; m->species[2].initial = .3;
  m->default_temperature = 287; m->default_pressure = 100000; finish(m);
  m->custom_defaults[chem_custom_id(m,"photo")] = .0005;
  m->custom_defaults[chem_custom_id(m,"aerosol.effective radius [m]")] = 1e-7;
  m->custom_defaults[chem_custom_id(m,"aerosol.particle number concentration [# m-3]")] = 2.5e6;
  return m;
}
static void init_robertson(ChemState *s, int multiple) {
  for (size_t c = 0; c < s->cells; ++c) {
    s->variables[chem_variable_id(s->model,"A")*s->stride+c] = multiple ? (double[]){1,2,.5}[c] : 1;
    s->conditions[c] = (ChemConditions){272.5,101253.3,1e6};
  }
  parameter(s,"r1",.04); parameter(s,"r2",3e7); parameter(s,"r3",1e4);
}
static void init_rate_example(ChemState *s) {
  set(s,"A",1);
  for (size_t c = 0; c < s->cells; ++c) {
    s->conditions[c] = (ChemConditions){287.45,101319.9,0}; chem_ideal_air_density(s->conditions+c);
  }
  const char *label = chem_custom_id(s->model,"C.effective radius [m]") != SIZE_MAX ? "C" : "C surface";
  char name[128]; snprintf(name,sizeof name,"%s.effective radius [m]",label); parameter(s,name,1e-7);
  snprintf(name,sizeof name,"%s.particle number concentration [# m-3]",label); parameter(s,name,2.5e6);
  if (chem_custom_id(s->model,"my emission rate") != SIZE_MAX) {
    parameter(s,"my emission rate",1e-20); parameter(s,"my loss rate",1e-23);
  }
}
static void header(const char *const *names, size_t n, int grid) {
  printf("%6s","time"); if (grid) printf(",%6s","grid");
  for (size_t i = 0; i < n; ++i) printf(",%11s",names[i]);
  putchar('\n');
}
static void state_line(ChemState *s, double time, const char *const *names, size_t n) {
  for (size_t c = 0; c < s->cells; ++c) {
    printf("%6g",time); if (s->cells > 1) printf(",%6zu",c);
    for (size_t i = 0; i < n; ++i) printf(",%11.2e",get(s,names[i],c));
    putchar('\n');
  }
}
static void print_stats(const ChemStats *s) {
  printf("accepted: %"PRIu64"\nfunction_calls: %"PRIu64"\njacobian_updates: %"PRIu64"\nnumber_of_steps: %"PRIu64"\naccepted: %"PRIu64"\nrejected: %"PRIu64"\ndecompositions: %"PRIu64"\nsolves: %"PRIu64"\n",
         s->accepted,s->function_calls,s->jacobian_updates,s->internal_steps,s->accepted,s->rejected,s->decompositions,s->solves);
}
static void add_stats(ChemStats *total, const ChemStats *s) {
  total->accepted += s->accepted; total->function_calls += s->function_calls; total->jacobian_updates += s->jacobian_updates;
  total->internal_steps += s->internal_steps; total->rejected += s->rejected; total->decompositions += s->decompositions; total->solves += s->solves;
}
static void compare_states(ChemState *a, ChemState *b) {
  check(a && b && a->cells == b->cells && a->model->nvariables == b->model->nvariables,"state dimensions");
  for (size_t i = 0; i < a->model->nvariables; ++i) {
    const char *name = a->model->species[a->model->variable_species[i]].name;
    for (size_t c = 0; c < a->cells; ++c) close_value(get(a,name,c),get(b,name,c));
  }
}
static void compare_processes(ChemState *a, ChemState *b) {
  chem_calculate_rates(a); chem_calculate_rates(b);
  check(a->model->nreactions == b->model->nreactions,"reaction count");
  for (size_t r = 0; r < a->model->nreactions; ++r) {
    size_t br = r;
    if (a->model->reactions[r].rate.type == CHEM_USER && strcmp(a->model->reactions[r].rate.label,b->model->reactions[br].rate.label)) {
      for (br = 0; br < b->model->nreactions; ++br)
        if (b->model->reactions[br].rate.type == CHEM_USER && !strcmp(a->model->reactions[r].rate.label,b->model->reactions[br].rate.label)) break;
    }
    check(br < b->model->nreactions,"reaction label");
    for (size_t c = 0; c < a->cells; ++c) close_value(a->rates[r*a->stride+c],b->rates[br*b->stride+c]);
  }
  size_t n = a->model->nvariables;
  chem_forcing(a,a->variables,a->forcing); chem_forcing(b,b->variables,b->forcing);
  chem_jacobian(a,a->variables,a->jacobian); chem_jacobian(b,b->variables,b->jacobian);
  for (size_t i = 0; i < n; ++i) {
    const char *name = a->model->species[a->model->variable_species[i]].name; size_t bi = chem_variable_id(b->model,name);
    for (size_t c = 0; c < a->cells; ++c) close_value(a->forcing[i*a->stride+c],b->forcing[bi*b->stride+c]);
    for (size_t j = 0; j < n; ++j) {
      const char *other = a->model->species[a->model->variable_species[j]].name; size_t bj = chem_variable_id(b->model,other);
      int ai = a->model->index[i*n+j], bb = b->model->index[bi*n+bj];
      for (size_t c = 0; c < a->cells; ++c) close_value(ai < 0 ? 0 : a->jacobian[(size_t)ai*a->stride+c],bb < 0 ? 0 : b->jacobian[(size_t)bb*b->stride+c]);
    }
  }
}
static void rates_tutorial(int user) {
  ChemModel *m = rate_example(user);
  ChemModel *config[2] = {user ? chem_load_example_user(0) : chem_load_example_no_user(0),user ? chem_load_example_user(1) : chem_load_example_no_user(1)};
  ChemState *s = chem_state_create(m,1), *other[2]; check(s != NULL,"state allocation"); init_rate_example(s);
  for (int i = 0; i < 2; ++i) { check(config[i] != NULL,"load configuration"); other[i] = chem_state_create(config[i],1); check(other[i] != NULL,"state allocation"); init_rate_example(other[i]); }
  const char *names[] = {"A","B","C","D","E","F","G"}; header(names,7,0); state_line(s,0,names,7);
  double photo = 1e-10; ChemSolver solver = chem_solver(3);
  for (int t = 0; t < 10; ++t) {
    if (user) { parameter(s,"my photolysis rate",photo); for (int i = 0; i < 2; ++i) parameter(other[i],"my photolysis rate",photo); }
    for (int i = 0; i < 2; ++i) compare_processes(s,other[i]);
    ChemStats stats; status(chem_advance(chem_step_rosenbrock,&solver,500,1,s,&stats));
    for (int i = 0; i < 2; ++i) { status(chem_advance(chem_step_rosenbrock,&solver,500,1,other[i],&stats)); compare_states(s,other[i]); }
    state_line(s,500*(t+1),names,7); photo *= 1.5;
  }
  for (int i = 0; i < 2; ++i) { chem_state_free(other[i]); chem_model_free(config[i]); }
  chem_state_free(s); chem_model_free(m);
}
static void robertson_tutorial(int method, size_t cells, int stats_only, int timing) {
  ChemModel *m = robertson(), *config = chem_load_example_robertson(1); check(config != NULL,"load configuration");
  ChemState *s = chem_state_create(m,cells), *other = chem_state_create(config,cells); check(s && other,"state allocation");
  init_robertson(s,cells > 1); init_robertson(other,cells > 1);
  const char *names[] = {"A","B","C"}; ChemSolver solver = chem_solver(method); ChemStats total = {0};
  if (!stats_only) { header(names,3,cells > 1); state_line(s,0,names,3); }
  double start = omp_get_wtime();
  for (int t = 0; t < (stats_only ? 1 : 10); ++t) {
    ChemStats stats; compare_processes(s,other); status(chem_advance(chem_step_rosenbrock,&solver,200,1,s,&stats)); add_stats(&total,&stats);
    if (stats_only) { puts("Solver state: Converged"); print_stats(&stats); printf("final simulation time: %g\n",stats.final_time); }
    status(chem_advance(chem_step_rosenbrock,&solver,200,1,other,&stats)); compare_states(s,other);
    if (!stats_only) state_line(s,200*(t+1),names,3);
  }
  if (timing) { printf("Total solve time: %.0f nanoseconds\n",(omp_get_wtime()-start)*1e9); print_stats(&total); }
  chem_state_free(other); chem_model_free(config); chem_state_free(s); chem_model_free(m);
}
static void readme_tutorial(void) {
  const char *names[] = {"Foo","Bar","Baz"}; ChemModel *m = new_model(names,3); ChemRate k = chem_rate(CHEM_ARRHENIUS);
  k.A = .001; reaction(m,(size_t[]){0},1,(size_t[]){1,2},(double[]){.8,.2},2,k);
  k.A = 1e-5; k.C = 110; reaction(m,(size_t[]){0,1},2,(size_t[]){2},(double[]){1},1,k); finish(m);
  ChemState *s = chem_state_create(m,1); check(s != NULL,"state allocation"); set(s,"Foo",20);
  s->conditions[0] = (ChemConditions){287.45,101319.9,0}; chem_ideal_air_density(s->conditions);
  ChemSolver solver = chem_solver(3); header(names,3,0);
  for (int i = 0; i < 10; ++i) { ChemStats stats; chem_calculate_rates(s); status(chem_solve(&solver,500,s,&stats)); state_line(s,500*i,names,3); }
  chem_state_free(s); chem_model_free(m);
}
static void openmp_tutorial(void) {
  ChemModel *m = robertson(); double results[3][3];
  #pragma omp parallel for num_threads(3) schedule(static)
  for (int i = 0; i < 3; ++i) {
    ChemState *s = chem_state_create(m,1); check(s != NULL,"state allocation"); init_robertson(s,0); ChemSolver solver = chem_solver(3);
    #pragma omp critical
    printf("Running solver on thread %d\n",omp_get_thread_num());
    for (int t = 0; t < 10; ++t) { ChemStats stats; chem_calculate_rates(s); status(chem_advance(chem_step_rosenbrock,&solver,200,1,s,&stats)); }
    const char *names[] = {"A","B","C"}; for (int j = 0; j < 3; ++j) results[i][j] = get(s,names[j],0);
    chem_state_free(s);
  }
  for (int i = 0; i < 3; ++i) { printf("Thread %d\n%10s,%10s,%10s\n%10.2e,%10.2e,%10.2e\n",i+1,"A","B","C",results[i][0],results[i][1],results[i][2]); }
  chem_model_free(m);
}
static void matrix_tutorial(void) {
  ChemArray *a = chem_array_dense(3,3,1,0), *v = chem_array_dense(3,3,3,0); check(a && v,"array allocation");
  for (size_t c = 0; c < 3; ++c) for (size_t j = 0; j < 3; ++j) {
    double x = 1.1+(double)c+.1*(double)j; status(chem_array_set(a,c,j,x)); status(chem_array_set(v,c,j,x));
  }
  for (size_t i = 0; i < 9; ++i) printf("%.6g\n",a->values[i]);
  putchar('\n');
  for (size_t i = 0; i < 9; ++i) printf("%.6g\n",v->values[i]);
  putchar('\n');
  /* The original demonstration's while condition is false initially, so its
   * Solve helper only updates parameters and leaves these variables unchanged. */
  for (size_t c = 0; c < 3; ++c) {
    printf("Cell %zu\n%10s%20s%20s\n",c,"Species","Regular","Vectorized");
    for (size_t j = 0; j < 3; ++j) {
      double av,vv; status(chem_array_get(a,c,j,&av)); status(chem_array_get(v,c,j,&vv));
      printf("%10s%20.6g%20.6g\n",(const char *[]){"A","B","C"}[j],av,vv);
    }
    putchar('\n');
  }
  chem_array_free(a); chem_array_free(v);
  for (size_t group = 1; group <= 3; group += 2) {
    a = chem_array_sparse(3,4,group,0,(size_t[]){0,2,2,3},(size_t[]){1,1,3,2},4,0); check(a != NULL,"sparse allocation");
    for (size_t c = 0; c < 3; ++c) for (size_t r = 0; r < 4; ++r) for (size_t col = 0; col < 4; ++col)
      if (a->index[r*4+col] >= 0) a->values[chem_sparse_index(a,c,r,col)] = (double)(100*c+10*r+col);
    printf("\nSparse matrix %s ordering elements\n",group == 1 ? "standard" : "vector");
    for (size_t i = 0; i < a->allocated; ++i) { unsigned x = (unsigned)a->values[i]; printf("%u.%u.%u\n",x/100,x/10%10,x%10); }
    chem_array_free(a);
  }
}
static void configuration_test(void) {
  ChemModel *models[] = {validation_model(),chem_load_validation(0),chem_load_validation(1)};
  ChemState *states[3];
  for (int i = 0; i < 3; ++i) { check(models[i] != NULL,"load validation configuration"); states[i] = chem_state_create(models[i],5); check(states[i] != NULL,"state allocation"); }
  for (int i = 0; i < 3; ++i) for (size_t c = 0; c < 5; ++c) {
    states[i]->conditions[c].temperature += 2*(double)c; states[i]->conditions[c].pressure -= 100*(double)c;
    chem_ideal_air_density(states[i]->conditions+c);
    for (size_t v = 0; v < 3; ++v) states[i]->variables[v*states[i]->stride+c] *= 1+.05*(double)c;
  }
  for (int step = 0; step < 10; ++step) {
    for (int i = 1; i < 3; ++i) { compare_states(states[0],states[i]); compare_processes(states[0],states[i]); }
    for (int i = 0; i < 3; ++i) { ChemSolver solver = chem_solver(3); ChemStats stats; status(chem_solve(&solver,30,states[i],&stats)); }
  }
  for (int i = 0; i < 3; ++i) { chem_state_free(states[i]); chem_model_free(models[i]); }
  ChemModel *a = chem_load_ts1(0), *b = chem_load_ts1(1); check(a && b,"load TS1 configuration");
  check(a->nspecies == 210 && a->nreactions == 547,"TS1 dimensions");
  ChemState *x = chem_state_create(a,5), *y = chem_state_create(b,5); check(x && y,"state allocation");
  for (size_t c = 0; c < 5; ++c) {
    x->conditions[c].temperature += (double)c*3; x->conditions[c].pressure += (double)c*100;
    chem_ideal_air_density(x->conditions+c); y->conditions[c] = x->conditions[c];
  }
  for (size_t i = 0; i < a->nvariables; ++i) for (size_t c = 0; c < 5; ++c) {
    double value = x->variables[i*x->stride+c]*(1+.01*(double)c);
    x->variables[i*x->stride+c] = y->variables[i*y->stride+c] = value;
  }
  for (size_t i = 0; i < a->ncustom*5; ++i) x->custom[i] = y->custom[i] = x->custom[i]*.01;
  const char *dark[] = {"jo2_a","jh2o_b","jh2o_c","jch4_a","jch4_b","jhf"};
  for (size_t i = 0; i < 6; ++i) { parameter(x,dark[i],0); parameter(y,dark[i],0); }
  compare_processes(x,y); ChemSolver solver = chem_solver(3); ChemStats stats;
  status(chem_advance(chem_step_rosenbrock,&solver,30,1,x,&stats)); status(chem_advance(chem_step_rosenbrock,&solver,30,1,y,&stats)); compare_states(x,y);
  chem_state_free(x); chem_state_free(y); chem_model_free(a); chem_model_free(b);
}
int main(int argc, char **argv) {
  check(argc == 2,"usage: tutorials CASE"); const char *name = argv[1];
  if (!strcmp(name,"README_example")) readme_tutorial();
  else if (!strcmp(name,"solve_results")) robertson_tutorial(3,3,1,0);
  else if (!strcmp(name,"multiple_grid_cells")) robertson_tutorial(3,3,0,0);
  else if (!strcmp(name,"rate_constants_no_user_defined_example_by_hand")) rates_tutorial(0);
  else if (!strcmp(name,"rate_constants_user_defined_example_by_hand")) rates_tutorial(1);
  else if (!strcmp(name,"solver_configuration")) {
    const char *labels[] = {"Two stages: ","Three stages: ","Four stages: ","Four stages differential algebraic: ","Six stages differential algebraic: "};
    for (int i = 0; i < 5; ++i) { if (i) putchar('\n'); puts(labels[i]); robertson_tutorial(i+2,1,0,1); }
  } else if (!strcmp(name,"openmp_tutorial")) openmp_tutorial();
  else if (!strcmp(name,"vectorized_matrix_solver")) matrix_tutorial();
  else if (!strcmp(name,"configuration")) configuration_test();
  else check(0,"unknown tutorial");
  return EXIT_SUCCESS;
}
