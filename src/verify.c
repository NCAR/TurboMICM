/* Central/one-sided differences and comparison semantics follow MICM's utility. */
#include "chem_verify.h"
#include <math.h>
#include <string.h>
#include <stdint.h>
#define VALUE(matrix, row, column) (*(double *)ck_matrix_at((matrix), (row), (column)))
#define CONST_VALUE(matrix, row, column) (*(const double *)ck_matrix_const_at((matrix), (row), (column)))
ck_status ck_finite_difference(ck_matrix_forcing function, void *context, const ck_matrix *variables,
                               double perturbation, ck_matrix *result)
{
  if (!function || !variables || !result || variables->element_size != sizeof(double) ||
      !(perturbation > 0) || (variables->columns && variables->columns > SIZE_MAX / variables->columns))
    return CK_INVALID;
  size_t n = variables->columns, b = variables->rows;
  double zero = 0;
  ck_matrix plus = {0}, minus = {0}, fp = {0}, fm = {0}, base = {0};
  ck_status status;
  status = ck_matrix_create(b, n * n, variables->width, sizeof(double), &zero, result);
  if (status)
    return status;
  if ((status = ck_matrix_clone(variables, &plus)) || (status = ck_matrix_clone(variables, &minus)) ||
      (status = ck_matrix_create(b, n, variables->width, sizeof(double), &zero, &fp)) ||
      (status = ck_matrix_create(b, n, variables->width, sizeof(double), &zero, &fm)) ||
      (status = ck_matrix_create(b, n, variables->width, sizeof(double), &zero, &base)))
    goto done;
  function(variables, &base, context);
  for (size_t col = 0; col < n; ++col)
  {
    memcpy(plus.data, variables->data, variables->size * sizeof(double));
    memcpy(minus.data, variables->data, variables->size * sizeof(double));
    for (size_t c = 0; c < b; ++c)
    {
      double x = CONST_VALUE(variables, c, col), h = perturbation * fmax(1, fabs(x));
      VALUE(&plus, c, col) = x + h;
      VALUE(&minus, c, col) = x - h;
    }
    ck_matrix_fill(&fp, &zero);
    ck_matrix_fill(&fm, &zero);
    function(&plus, &fp, context);
    function(&minus, &fm, context);
    for (size_t c = 0; c < b; ++c)
    {
      double x = CONST_VALUE(variables, c, col), h = perturbation * fmax(1, fabs(x));
      for (size_t row = 0; row < n; ++row)
        VALUE(result, c, row * n + col) = x - h >= 0 ? (VALUE(&fp, c, row) - VALUE(&fm, c, row)) / (2 * h)
                                                     : (VALUE(&fp, c, row) - VALUE(&base, c, row)) / h;
    }
  }
done:
  ck_matrix_free(&plus);
  ck_matrix_free(&minus);
  ck_matrix_free(&fp);
  ck_matrix_free(&fm);
  ck_matrix_free(&base);
  if (status)
    ck_matrix_free(result);
  return status;
}
static ck_status compare(const ck_sparse *a, const ck_matrix *fd, double atol, double rtol, int sparsity,
                         ck_jacobian_comparison *r)
{
  if (!a || !fd || !r || a->values.element_size != sizeof(double) || fd->element_size != sizeof(double) ||
      fd->rows != a->values.rows || fd->columns != a->order * a->order || !(atol >= 0) || !(rtol >= 0))
    return CK_INVALID;
  memset(r, 0, sizeof *r);
  r->passed = 1;
  for (size_t b = 0; b < fd->rows; ++b)
    for (size_t row = 0; row < a->order; ++row)
      for (size_t col = 0; col < a->order; ++col)
      {
        const double *entry = ck_sparse_const_at(a, b, row, col);
        if (sparsity ? entry != NULL : entry == NULL)
          continue;
        double f = CONST_VALUE(fd, b, row * a->order + col), analytic = sparsity ? 0 : -*entry;
        double error = fabs(analytic - f), scale = fmax(fabs(analytic), fabs(f)),
               rel = scale > 0 ? error / scale : 0;
        if (error > r->max_abs_error)
        {
          r->max_abs_error = error;
          r->worst_block = b;
          r->worst_row = row;
          r->worst_column = col;
          r->worst_analytical = analytic;
          r->worst_fd = sparsity ? fabs(f) : f;
        }
        if (rel > r->max_rel_error)
          r->max_rel_error = rel;
        if (error > atol + rtol * scale)
          r->passed = 0;
      }
  return CK_OK;
}
ck_status ck_compare_jacobian(const ck_sparse *a, const ck_matrix *fd, double atol, double rtol,
                              ck_jacobian_comparison *r)
{
  return compare(a, fd, atol, rtol, 0, r);
}
ck_status ck_check_sparsity(const ck_sparse *a, const ck_matrix *fd, double threshold,
                            ck_jacobian_comparison *r)
{
  return compare(a, fd, threshold, 0, 1, r);
}
