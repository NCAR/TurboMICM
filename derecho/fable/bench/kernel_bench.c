/* Per-kernel timing of the generated TS1 kernels on one state: forcing,
 * Jacobian, LU, and solve, in nanoseconds per cell per call.
 *
 *   kernel_bench [--group-size G] [--repeat N]
 */
#include "internal.h"
#include "ts1.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

int main(int argc, char **argv)
{
  int group = 8, repeat = 2000;
  for (int i = 1; i < argc; ++i)
  {
    if (strcmp(argv[i], "--group-size") == 0 && i + 1 < argc)
      group = atoi(argv[++i]);
    else if (strcmp(argv[i], "--repeat") == 0 && i + 1 < argc)
      repeat = atoi(argv[++i]);
  }
  tmc_system system = ts1_system();
  tmc_reaction *reactions;
  int n_reactions = ts1_reactions(&reactions);
  tmc_builder *builder = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(builder, &system);
  tmc_builder_set_reactions(builder, reactions, n_reactions);
  tmc_builder_set_ignore_unused_species(builder, true);
  tmc_builder_set_kernels(builder, ts1_kernels());
  tmc_error error;
  tmc_solver *solver = tmc_builder_build(builder, &error);
  if (!solver)
  {
    fprintf(stderr, "%s\n", error.message);
    return 1;
  }
  tmc_state *st = tmc_solver_create_state(solver, group);
  ts1_initial_state(st);
  tmc_conditions *cond = tmc_state_conditions(st);
  for (int c = 0; c < group; ++c)
  {
    cond[c].temperature = 260.0 + 50.0 * c / (group > 1 ? group - 1 : 1);
    cond[c].pressure = 101319.9 * (0.6 + 0.4 * c / (group > 1 ? group - 1 : 1));
    tmc_conditions_calculate_ideal_air_density(&cond[c]);
  }
  tmc_solver_update_state_parameters(solver, st);
  /* one real step so the state holds typical values */
  tmc_solver_solve(solver, 30.0, st);

  const tmc_kernels *k = solver->kernels;
  const int nb = st->n_blocks;
  double *y = st->variables->data, *rc = st->rate_constants->data, *f = st->initial_forcing->data;
  double *lu = st->jacobian->data, *recip = st->reciprocals, *x = st->K[0]->data;
  const double alpha = 1.0 / (30.0 / 9 * 0.4358665215);
  tmc_dense *x0 = tmc_dense_copy(st->K[0]);
  double cells = (double)group * repeat;

  for (int r = 0; r < repeat / 10; ++r)
  {
    k->forcing(y, rc, f, nb);
    k->jacobian(y, rc, alpha, lu, nb);
    k->lu_decompose(lu, recip, nb);
    k->lu_solve(x, lu, recip, nb);
  }
  double t0 = omp_get_wtime();
  for (int r = 0; r < repeat; ++r)
    k->forcing(y, rc, f, nb);
  double t1 = omp_get_wtime();
  for (int r = 0; r < repeat; ++r)
    k->jacobian(y, rc, alpha, lu, nb);
  double t2 = omp_get_wtime();
  for (int r = 0; r < repeat; ++r)
  {
    k->jacobian(y, rc, alpha, lu, nb);
    k->lu_decompose(lu, recip, nb);
  }
  double t3 = omp_get_wtime();
  for (int r = 0; r < repeat; ++r)
  {
    memcpy(x, x0->data, x0->size * sizeof(double));
    k->lu_solve(x, lu, recip, nb);
  }
  double t4 = omp_get_wtime();
  for (int r = 0; r < repeat; ++r)
    memcpy(x, x0->data, x0->size * sizeof(double));
  double t5 = omp_get_wtime();
  {
    long subnormal = 0, nonfinite = 0, tiny = 0;
    for (size_t i = 0; i < st->jacobian->size; ++i)
    {
      double v = fabs(lu[i]);
      if (v != 0.0 && v < DBL_MIN) ++subnormal;
      if (!isfinite(v)) ++nonfinite;
      if (v != 0.0 && v < 1e-290) ++tiny;
    }
    long sub_x = 0;
    for (size_t i = 0; i < x0->size; ++i) { double v = fabs(x[i]); if (v != 0.0 && v < DBL_MIN) ++sub_x; }
    printf("LU matrix: %zu values, %ld subnormal, %ld below 1e-290, %ld non-finite; solution: %ld subnormal\n", st->jacobian->size, subnormal, tiny, nonfinite, sub_x);
  }
  double jac = (t2 - t1) / cells * 1e9, luc = (t3 - t2) / cells * 1e9 - jac, cp = (t5 - t4) / cells * 1e9;
  double sv = (t4 - t3) / cells * 1e9 - cp, fo = (t1 - t0) / cells * 1e9;
  printf("group %d (%d blocks of %d): ns per cell: forcing %.0f, jacobian %.0f, lu %.0f, solve %.0f; per internal step (1 J + 1 LU + 2 F + 3 S): %.0f ns\n",
         group, nb, st->L, fo, jac, luc, sv, jac + luc + 2 * fo + 3 * sv);
  return 0;
}
