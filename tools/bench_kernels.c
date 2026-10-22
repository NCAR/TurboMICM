/* Development microbenchmark; each kernel uses the same initialized state. */
#include "internal.h"
#include "chem_specialized.h"
#include <omp.h>
#include <stdio.h>
extern const ck_specialized ts1_specialized;
int main(void)
{
  ck_model *m;
  ck_state *s;
  if (ck_model_create(ts1_specialized.mechanism, &m) || ck_state_create(m, 32, &s) || ck_update_rates(s))
    return 1;
  ck_workspace *w = s->workspace;
  double *matrix = calloc(m->nnz * s->cells, sizeof(double));
  if (!matrix)
    return 1;
  ck_jacobian(s, s->y, matrix);
  ck_alpha_minus_jacobian(m, s->cells, matrix, 100);
  for (int specialized = 0; specialized < 2; ++specialized)
  {
    const int repeats = 2000;
    double start = omp_get_wtime();
    for (int i = 0; i < repeats; ++i)
      if (specialized) ts1_specialized.forcing(s, s->y, w->initial);
      else ck_forcing(s, s->y, w->initial);
    printf("%s forcing %.3f us\n", specialized ? "generated" : "generic", (omp_get_wtime()-start)*1e6/repeats);
    start = omp_get_wtime();
    for (int i = 0; i < repeats; ++i)
      if (specialized) ts1_specialized.jacobian(s, s->y, w->jacobian);
      else ck_jacobian(s, s->y, w->jacobian);
    printf("%s jacobian %.3f us\n", specialized ? "generated" : "generic", (omp_get_wtime()-start)*1e6/repeats);
    start = omp_get_wtime();
    for (int i = 0; i < repeats; ++i)
    {
      memcpy(w->jacobian, matrix, m->nnz*s->cells*sizeof(double));
      if (specialized) ts1_specialized.factor(s->cells, w->jacobian);
      else ck_factor(m, s->cells, w->jacobian, CK_LU_MOZART);
    }
    printf("%s factor+copy %.3f us\n", specialized ? "generated" : "generic", (omp_get_wtime()-start)*1e6/repeats);
    start = omp_get_wtime();
    for (int i = 0; i < repeats; ++i)
    {
      memcpy(w->stage[0], s->y, m->n*s->cells*sizeof(double));
      if (specialized) ts1_specialized.linear_solve(s->cells, w->jacobian, w->stage[0]);
      else ck_linear_solve(m, s->cells, w->jacobian, w->stage[0]);
    }
    printf("%s solve+copy %.3f us\n", specialized ? "generated" : "generic", (omp_get_wtime()-start)*1e6/repeats);
  }
  free(matrix);
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
