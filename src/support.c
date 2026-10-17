/* Support functions for the ported MICM tests (see include/tmc_support.h). */
#include "internal.h"
#include "tmc_support.h"

/* ---- lambda rate constants by label ---- */

tmc_lambda_parameters *tmc_support_solver_lambda_rate_constant(tmc_solver *s, const char *label, tmc_error *err)
{
  for (int r = 0; r < s->n_reactions; ++r)
  {
    tmc_rate_constant *rc = &s->reactions[r].rate_constant;
    if (rc->type == TMC_RATE_LAMBDA && strcmp(rc->p.lambda.label, label) == 0)
      return &rc->p.lambda;
  }
  tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_RATE_CONSTANT_NOT_FOUND,
                "Lambda rate constant with name '%s' not found in any process", label);
  return NULL;
}

/* ---- unique names in a given order ---- */

int tmc_support_system_unique_names_ordered(const tmc_system *s, tmc_support_name_order order, void *context,
                                            char (*names)[TMC_NAME_MAX], int max_names)
{
  int n = tmc_system_unique_names(s, NULL, 0);
  char (*declared)[TMC_NAME_MAX] = tmc_malloc((size_t)(n > 0 ? n : 1) * sizeof *declared);
  tmc_system_unique_names(s, declared, n);
  for (int i = 0; i < n && i < max_names; ++i)
  {
    const char *name = order ? order((const char(*)[TMC_NAME_MAX])declared, n, i, context) : declared[i];
    tmc_strlcpy(names[i], name, TMC_NAME_MAX);
  }
  free(declared);
  return n;
}
