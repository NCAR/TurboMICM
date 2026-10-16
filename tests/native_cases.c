/* Native concurrency checks run without Python's numerical-library runtimes.
 * SPDX-License-Identifier: Apache-2.0 */
#include "chem.h"
#include <stdio.h>
#include <string.h>

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
  fputs("usage: test_native_cpu openmp\n",stderr);
  return 2;
}
