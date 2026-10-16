/* Native ports of the supplied tutorial programs.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "chem.h"
#include "chem_array.h"
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>

ch_model *ch_build_example_no_user(void);
ch_model *ch_build_example_user(void);
ch_model *ch_build_example_robertson(void);

static void check(int result) {
  if (result < 0) { fprintf(stderr,"%s\n",ch_status_string(result)); exit(EXIT_FAILURE); }
}
static void *allocated(void *p) {
  if (!p) { fputs("Allocation failed\n",stderr); exit(EXIT_FAILURE); }
  return p;
}
static void set(ch_state *s, const char *name, double value) { check(ch_state_set(s,name,&value,1)); }
static void parameter(ch_state *s, const char *name, double value) { check(ch_state_set_parameter(s,name,&value,1)); }
static void reaction(ch_model *m, size_t nr, const int *r, size_t np, const int *p, const double *y, ch_rate rate) {
  check(ch_model_add_reaction(m,nr,r,np,p,y,&rate));
}
static ch_rate user_rate(ch_model *m, const char *label) {
  ch_rate rate = ch_rate_default(CH_USER);
  rate.custom_index = ch_model_add_parameter(m,label,0); check(rate.custom_index); return rate;
}
static ch_model *robertson(void) {
  ch_model *m = allocated(ch_model_create());
  const char *names[] = {"A","B","C"};
  for (int i = 0; i < 3; ++i) check(ch_model_add_species(m,names[i],1e-3,0));
  reaction(m,1,(int[]){0},1,(int[]){1},(double[]){1},user_rate(m,"r1"));
  reaction(m,2,(int[]){1,1},2,(int[]){1,2},(double[]){1,1},user_rate(m,"r2"));
  reaction(m,2,(int[]){1,2},2,(int[]){0,2},(double[]){1,1},user_rate(m,"r3"));
  check(ch_model_finalize(m,1,0)); return m;
}
static ch_model *readme_model(void) {
  ch_model *m = allocated(ch_model_create());
  const char *names[] = {"Foo","Bar","Baz"};
  for (int i = 0; i < 3; ++i) check(ch_model_add_species(m,names[i],1e-3,0));
  ch_rate rate = ch_rate_default(CH_ARRHENIUS); rate.p[0] = 1e-3;
  reaction(m,1,(int[]){0},2,(int[]){1,2},(double[]){.8,.2},rate);
  rate.p[0] = 1e-5; rate.p[2] = 110;
  reaction(m,2,(int[]){0,1},1,(int[]){2},(double[]){1},rate);
  check(ch_model_finalize(m,1,0)); return m;
}
static ch_model *rate_model(int user) {
  ch_model *m = allocated(ch_model_create());
  const char *names[] = {"A","B","C","D","E","F","G"};
  for (int i = 0; i < 7; ++i) check(ch_model_add_species(m,names[i],1e-3,0));
  check(ch_model_set_species_property(m,2,"molecular weight [kg mol-1]",.025));
  check(ch_model_set_species_property(m,2,"diffusion coefficient [m2 s-1]",230.));
  ch_rate r = ch_rate_default(CH_ARRHENIUS); r.p[0] = 2.15e-4; r.p[2] = 110;
  reaction(m,1,(int[]){0},1,(int[]){1},(double[]){1},r);
  r = ch_rate_default(CH_BRANCHED); r.p[0] = 1.2; r.p[1] = 204.3; r.p[2] = .001; r.p[3] = 2;
  reaction(m,1,(int[]){1},1,(int[]){2},(double[]){1},r);
  r.p[4] = 1; reaction(m,1,(int[]){1},1,(int[]){3},(double[]){1},r);
  r = ch_rate_default(CH_SURFACE); r.p[2] = .9;
  r.custom_index = ch_model_add_parameter(m,"C surface.effective radius [m]",0); check(r.custom_index);
  check(ch_model_add_parameter(m,"C surface.particle number concentration [# m-3]",0));
  reaction(m,1,(int[]){2},1,(int[]){4},(double[]){1},r);
  const double conversion = 1e-6*CH_AVOGADRO;
  r = ch_rate_default(CH_TERNARY);
  r.p[0] = 1.2; r.p[1] = 2.3; r.p[2] = 302.3; r.p[3] = 2.6/conversion;
  r.p[4] = -3.1; r.p[5] = 402.1; r.p[6] = .9; r.p[7] = 1.2;
  reaction(m,1,(int[]){3},1,(int[]){5},(double[]){2},r);
  r = ch_rate_default(CH_TROE); r.p[0] = 1.8e-30*conversion*conversion; r.p[1] = -3.; r.p[3] = 2.8e-11*conversion;
  reaction(m,2,(int[]){4,4},1,(int[]){6},(double[]){1},r);
  r = ch_rate_default(CH_TUNNELING); r.p[0] = 1.2; r.p[1] = 2.3; r.p[2] = 302.3;
  reaction(m,1,(int[]){5},1,(int[]){6},(double[]){1},r);
  if (user) {
    reaction(m,1,(int[]){2},1,(int[]){6},(double[]){1},user_rate(m,"my photolysis rate"));
    reaction(m,0,NULL,1,(int[]){0},(double[]){1},user_rate(m,"my emission rate"));
    reaction(m,1,(int[]){1},0,NULL,NULL,user_rate(m,"my loss rate"));
  }
  check(ch_model_finalize(m,1,0)); return m;
}
static ch_state *initial(ch_model *m, size_t cells, int kind) {
  ch_state *s = allocated(ch_state_create(m,cells));
  ch_conditions *cond = ch_state_conditions(s);
  for (size_t c = 0; c < cells; ++c) {
    cond[c] = kind < 2 ? (ch_conditions){287.45,101319.9,0} : (ch_conditions){272.5,101253.3,1e6};
    if (kind < 2) ch_ideal_air_density(&cond[c]);
  }
  set(s,"A",1.);
  if (kind < 2) {
    parameter(s,"C surface.effective radius [m]",1e-7);
    parameter(s,"C surface.particle number concentration [# m-3]",2.5e6);
    if (kind) { parameter(s,"my emission rate",1e-20); parameter(s,"my loss rate",1e-23); }
  } else {
    parameter(s,"r1",.04); parameter(s,"r2",3e7); parameter(s,"r3",1e4);
    if (cells == 3) check(ch_state_set(s,"A",(double[]){1,2,.5},3));
  }
  return s;
}
static void header(size_t cells, size_t n, const char *const *names) {
  printf("%6s","time"); if (cells > 1) printf(",%6s","grid");
  for (size_t i = 0; i < n; ++i) printf(",%11s",names[i]);
  putchar('\n');
}
static void print_state(ch_state *s, int time, size_t n, const char *const *names) {
  size_t nc = ch_state_cells(s); const ch_model *m = ch_state_model(s); double *y = ch_state_values(s);
  for (size_t c = 0; c < nc; ++c) {
    printf("%6d",time); if (nc > 1) printf(",%6zu",c);
    for (size_t i = 0; i < n; ++i) printf(",%11.2e",y[(size_t)ch_model_species_index(m,names[i])*nc+c]);
    putchar('\n');
  }
}
static void advance(ch_state *s, double dt, int method, ch_stats *total) {
  double elapsed = 0;
  ch_rosenbrock_parameters p = ch_rosenbrock_default(method);
  while (elapsed < dt) {
    ch_stats stats; check(ch_update_rates(s));
    int result = ch_solve(s,dt-elapsed,&p,&stats);
    if (result != CH_MAX_STEPS) check(result);
    if (!(stats.final_time > 0)) { fputs("Solver made no progress\n",stderr); exit(EXIT_FAILURE); }
    if (total) ch_stats_add(total,&stats);
    elapsed += stats.final_time;
  }
}
static void print_stats(const ch_stats *s) {
  printf("accepted: %zu\nfunction_calls: %zu\njacobian_updates: %zu\nnumber_of_steps: %zu\n",
         s->accepted,s->function_calls,s->jacobian_updates,s->internal_steps);
  printf("accepted: %zu\nrejected: %zu\ndecompositions: %zu\nsolves: %zu\n",s->accepted,s->rejected,s->decompositions,s->solves);
}

static int readme(void) {
  ch_model *m = readme_model(); ch_state *s = allocated(ch_state_create(m,1));
  const char *names[] = {"Foo","Bar","Baz"};
  ch_state_conditions(s)[0] = (ch_conditions){287.45,101319.9,0}; ch_ideal_air_density(ch_state_conditions(s));
  set(s,"Foo",20.); header(1,3,names);
  for (int i = 0; i < 10; ++i) { advance(s,500,3,NULL); print_state(s,i*500,3,names); }
  ch_state_free(s); ch_model_free(m); return 0;
}
static int rates_example(int user) {
  ch_model *m = rate_model(user); ch_state *s = initial(m,1,user);
  const char *names[] = {"A","B","C","D","E","F","G"};
  header(1,7,names); print_state(s,0,7,names); double photo = 1e-10;
  for (int i = 0; i < 10; ++i) {
    if (user) parameter(s,"my photolysis rate",photo);
    advance(s,500,3,NULL); print_state(s,500*(i+1),7,names); photo *= 1.5;
  }
  ch_state_free(s); ch_model_free(m); return 0;
}
static int multiple(int statistics_only) {
  ch_model *m = robertson(); ch_state *s = initial(m,3,2); const char *names[] = {"A","B","C"};
  if (statistics_only) {
    ch_stats stats; ch_rosenbrock_parameters p = ch_rosenbrock_default(3);
    check(ch_update_rates(s)); int status = ch_solve(s,200,&p,&stats); check(status);
    printf("Solver state: %s\n",ch_status_string(status)); print_stats(&stats);
    printf("final simulation time: %g\n",stats.final_time);
  } else {
    header(3,3,names); print_state(s,0,3,names);
    for (int i = 0; i < 10; ++i) { advance(s,200,3,NULL); print_state(s,200*(i+1),3,names); }
  }
  ch_state_free(s); ch_model_free(m); return 0;
}
static int configuration(void) {
  const int methods[] = {2,3,4,43,6};
  const char *labels[] = {"Two stages: ","Three stages: ","Four stages: ",
                         "Four stages differential algebraic: ","Six stages differential algebraic: "};
  const char *names[] = {"A","B","C"}; ch_model *m = robertson();
  for (int j = 0; j < 5; ++j) {
    if (j) putchar('\n'); puts(labels[j]);
    ch_state *s = initial(m,1,2); ch_stats total = {0}; double seconds = 0;
    header(1,3,names); print_state(s,0,3,names);
    for (int i = 0; i < 10; ++i) {
      double start = omp_get_wtime(); advance(s,200,methods[j],&total); seconds += omp_get_wtime()-start;
      print_state(s,200*(i+1),3,names);
    }
    printf("Total solve time: %.0f nanoseconds\n",seconds*1e9); print_stats(&total); ch_state_free(s);
  }
  ch_model_free(m); return 0;
}
static int openmp_tutorial(void) {
  ch_model *m = robertson(); ch_state *states[3]; ch_stats stats[3]; int statuses[3];
  ch_rosenbrock_parameters p = ch_rosenbrock_default(3);
  for (int i = 0; i < 3; ++i) { states[i] = initial(m,1,2); check(ch_update_rates(states[i])); }
  for (int i = 0; i < 3; ++i) printf("Running solver on thread %d\n",i);
  for (int i = 0; i < 10; ++i)
    if (ch_solve_many(states,3,200,&p,3,statuses,stats)) { fputs("Parallel solve failed\n",stderr); return 1; }
  for (int i = 0; i < 3; ++i) {
    printf("Thread %d\n%10s,%10s,%10s\n",i+1,"A","B","C");
    double *y = ch_state_values(states[i]);
    printf("%10.2e,%10.2e,%10.2e\n",y[ch_model_species_index(m,"A")],y[ch_model_species_index(m,"B")],y[ch_model_species_index(m,"C")]);
    ch_state_free(states[i]);
  }
  ch_model_free(m); return 0;
}
static int vectorized(void) {
  ch_array *a[2];
  for (int group = 0; group < 2; ++group) {
    a[group] = allocated(ch_array_dense(3,3,group ? 3 : 0,CH_ARRAY_REAL,0));
    for (int c = 0; c < 3; ++c) for (int s = 0; s < 3; ++s)
      check(ch_array_set(a[group],(size_t)c,(size_t)s,(double)(10*(c+1)+s+1)/10));
    double *v = ch_array_data(a[group]);
    for (size_t i = 0; i < ch_array_storage_size(a[group]); ++i) printf("%g\n",v[i]);
    putchar('\n');
  }
  /* The reference tutorial's loop condition skips Solve. Preserve its output;
   * actual vector integrations are covered by the named integration cases. */
  ch_model *m = robertson(); ch_state *state = initial(m,3,2);
  for (int i = 0; i < 10; ++i) check(ch_update_rates(state));
  ch_state_free(state); ch_model_free(m);
  for (int c = 0; c < 3; ++c) {
    printf("Cell %d\n%10s%20s%20s\n",c,"Species","Regular","Vectorized");
    for (int s = 0; s < 3; ++s) {
      double x,y; check(ch_array_get(a[0],(size_t)c,(size_t)s,&x)); check(ch_array_get(a[1],(size_t)c,(size_t)s,&y));
      printf("%10c%20g%20g\n",'A'+s,x,y);
    }
    putchar('\n');
  }
  for (int i = 0; i < 2; ++i) ch_array_free(a[i]);
  const int rows[] = {0,2,2,3},cols[] = {1,1,3,2}; unsigned char pattern[16] = {0};
  for (int i = 0; i < 4; ++i) pattern[rows[i]*4+cols[i]] = 1;
  char strings[12][16];
  for (int c = 0; c < 3; ++c) for (int e = 0; e < 4; ++e) snprintf(strings[c*4+e],16,"%d.%d.%d",c,rows[e],cols[e]);
  for (int group = 0; group < 2; ++group) {
    ch_array *s = allocated(ch_array_sparse(3,4,pattern,0,group ? 3 : 0,CH_ARRAY_INT,0));
    for (int c = 0; c < 3; ++c) for (int e = 0; e < 4; ++e) check(ch_array_set(s,(size_t)c,(size_t)e,c*4+e));
    printf("\nSparse matrix %s ordering elements\n",group ? "vector" : "standard");
    int *ids = ch_array_data(s);
    for (size_t i = 0; i < ch_array_storage_size(s); ++i) puts(strings[ids[i]]);
    ch_array_free(s);
  }
  return 0;
}

static int near(double a, double b) { return isfinite(a) && isfinite(b) && fabs(a-b) <= 1e-30+1e-12*fabs(b); }
static int parity(void) {
  for (int kind = 0; kind < 3; ++kind) {
    ch_model *hand = kind == 2 ? robertson() : rate_model(kind);
    ch_model *config = allocated(kind == 2 ? ch_build_example_robertson() : kind ? ch_build_example_user() : ch_build_example_no_user());
    if (ch_model_species_count(hand) != ch_model_species_count(config) || ch_model_reaction_count(hand) != ch_model_reaction_count(config) ||
        ch_model_parameter_count(hand) != ch_model_parameter_count(config)) return 1;
    ch_state *a = initial(hand,3,kind), *b = initial(config,3,kind);
    for (int c = 0; c < 3; ++c) {
      ch_state_conditions(a)[c].temperature += c*10;
      ch_state_conditions(b)[c] = ch_state_conditions(a)[c];
    }
    for (int step = 0; step < 10; ++step) {
      if (kind == 1) {
        parameter(a,"my photolysis rate",1e-10*pow(1.5,step)); parameter(b,"my photolysis rate",1e-10*pow(1.5,step));
      }
      advance(a,kind == 2 ? 200 : 500,3,NULL); advance(b,kind == 2 ? 200 : 500,3,NULL);
      for (size_t i = 0; i < ch_model_species_count(hand); ++i) {
        int j = ch_model_species_index(config,ch_model_species_name(hand,i)); if (j < 0) return 1;
        for (size_t c = 0; c < 3; ++c)
          if (!near(ch_state_values(a)[i*3+c],ch_state_values(b)[(size_t)j*3+c])) {
            fprintf(stderr,"Configuration mismatch: example %d, step %d, species %s\n",kind,step,ch_model_species_name(hand,i));
            return 1;
          }
      }
    }
    ch_state_free(a); ch_state_free(b); ch_model_free(hand); ch_model_free(config);
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc != 2) { fputs("usage: tutorial_cpu CASE\n",stderr); return 2; }
  if (!strcmp(argv[1],"README_example")) return readme();
  if (!strcmp(argv[1],"solve_results")) return multiple(1);
  if (!strcmp(argv[1],"multiple_grid_cells")) return multiple(0);
  if (!strcmp(argv[1],"rate_constants_no_user_defined_example_by_hand")) return rates_example(0);
  if (!strcmp(argv[1],"rate_constants_user_defined_example_by_hand")) return rates_example(1);
  if (!strcmp(argv[1],"solver_configuration")) return configuration();
  if (!strcmp(argv[1],"vectorized_matrix_solver")) return vectorized();
  if (!strcmp(argv[1],"openmp_tutorial")) return openmp_tutorial();
  if (!strcmp(argv[1],"config_parity")) return parity();
  return 2;
}
