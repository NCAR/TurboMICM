/* Native concurrency checks run without Python's numerical-library runtimes.
 * SPDX-License-Identifier: Apache-2.0 */
#include "chem.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct { size_t factors, solves, controls; } hook_counts;

static int alternative_factor(const ch_sparse *pattern, size_t cells, double *a, double *scratch, void *context) {
  (void)scratch;
  hook_counts *counts = context;
  ++counts->factors;
  size_t n = ch_sparse_dimension(pattern);
  const int *rows = ch_sparse_row_offsets(pattern), *columns = ch_sparse_column_indices(pattern);
  const int *diagonal = ch_sparse_diagonal_indices(pattern);
  if (n != 2 || (size_t)rows[n] != ch_sparse_size(pattern)) return CH_SIZE;
  for (size_t i = 0; i < n; ++i)
    if (columns[diagonal[i]] != (int)i || diagonal[i] < rows[i] || diagonal[i] >= rows[i+1]) return CH_SIZE;
  return ch_lu_factor(pattern,cells,a,1);
}
static void alternative_solve(const ch_sparse *pattern, size_t cells, const double *a, double *x, void *context) {
  ++((hook_counts *)context)->solves;
  ch_lu_solve(pattern,cells,a,x);
}
static double alternative_controller(double error, const ch_rosenbrock_parameters *p, void *context) {
  ++((hook_counts *)context)->controls;
  return ch_default_controller(error,p,NULL);
}

/* A caller-supplied parameter table, controller, and linear solver must work
 * together without any chemistry-engine changes or mechanism specialization. */
static int extension_case(void) {
  int result = 1;
  ch_model *m = ch_model_create(); ch_state *baseline = NULL, *modified = NULL;
  if (!m) return 1;
  if (ch_model_add_species(m,"A",1e-10,0) < 0 || ch_model_add_species(m,"B",1e-10,0) < 0) goto done;
  int reactant = 0, product = 1; double yield = 1;
  ch_rate rate = ch_rate_default(CH_ARRHENIUS); rate.p[0] = 2;
  if (ch_model_add_reaction(m,1,&reactant,1,&product,&yield,&rate) < 0 || ch_model_finalize(m,0,0)) goto done;
  baseline = ch_state_create(m,3);
  if (!baseline) goto done;
  for (size_t i = 0; i < 3; ++i) {
    ch_state_conditions(baseline)[i] = (ch_conditions){280+10*i,101325,40};
    ch_state_values(baseline)[i] = 1+i;
  }
  if (ch_update_rates(baseline)) goto done;
  modified = ch_state_clone(baseline);
  if (!modified) goto done;
  hook_counts counts = {0};
  if (ch_state_set_linear_solver(modified,alternative_factor,alternative_solve,&counts)) goto done;
  ch_rosenbrock_parameters p = ch_rosenbrock_default(3), custom = p;
  custom.controller = alternative_controller; custom.controller_context = &counts;
  ch_stats a, b;
  if (ch_solve(baseline,.1,&p,&a) || ch_solve(modified,.1,&custom,&b)) goto done;
  if (counts.factors != b.decompositions || counts.solves != b.solves ||
      counts.controls != b.internal_steps || !counts.controls) goto done;
  for (size_t i = 0; i < 6; ++i)
    if (fabs(ch_state_values(baseline)[i]-ch_state_values(modified)[i]) > 1e-12) goto done;
  result = 0;
done:
  ch_state_free(baseline); ch_state_free(modified); ch_model_free(m);
  if (result) fputs("Native extension hook check failed\n",stderr);
  return result;
}

static int openmp_case(void) {
  int status = 1;
  ch_model *m = ch_model_create();
  ch_state *states[8] = {0};
  if (!m) return 1;
  const char *names[] = {"A","B","C"};
  for (int i = 0; i < 3; ++i) if (ch_model_add_species(m,names[i],1e-3,0) < 0) goto done;
  const int nr[] = {1,2,2}, np[] = {1,2,2};
  const int reactants[3][2] = {{0,0},{1,1},{1,2}}, products[3][2] = {{1,0},{1,2},{0,2}};
  const double rates[] = {.04,3e7,1e4}, yields[] = {1,1};
  for (int i = 0; i < 3; ++i) {
    ch_rate r = ch_rate_default(CH_USER);
    r.custom_index = ch_model_add_parameter(m,names[i],rates[i]);
    if (r.custom_index < 0 || ch_model_add_reaction(m,(size_t)nr[i],reactants[i],(size_t)np[i],products[i],yields,&r)<0)
      goto done;
  }
  if (ch_model_finalize(m,0,0)) goto done;
  for (int i = 0; i < 8; ++i) {
    states[i] = ch_state_create(m,1);
    if (!states[i]) goto done;
    ch_state_values(states[i])[0] = 1;
    ch_state_conditions(states[i])[0] = (ch_conditions){272.5,101253.3,1e6};
    if (ch_update_rates(states[i])) goto done;
  }
  ch_rosenbrock_parameters p = ch_rosenbrock_default(3);
  ch_stats stats[8]; int results[8];
  for (int i = 0; i < 10; ++i) if (ch_solve_many(states,8,200,&p,8,results,stats)) goto done;
  for (int i = 1; i < 8; ++i)
    for (int j = 0; j < 3; ++j)
      if (ch_state_values(states[i])[j] != ch_state_values(states[0])[j]) goto done;
  if (!(ch_state_values(states[0])[0] < 1)) goto done;
  status = 0;
done:
  for (int i = 0; i < 8; ++i) ch_state_free(states[i]);
  ch_model_free(m);
  return status;
}

int main(int argc, char **argv) {
  if (argc == 2 && !strcmp(argv[1],"openmp")) return openmp_case();
  if (argc == 2 && !strcmp(argv[1],"extensions")) return extension_case();
  fputs("usage: test_native_cpu openmp|extensions\n",stderr);
  return 2;
}
