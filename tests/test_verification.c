#include "test.h"
#include "chem_verify.h"
static void simple(const ck_matrix *v, ck_matrix *f, void *context)
{
  (void)context;
  for (size_t b = 0; b < v->rows; ++b)
  {
    double x = *(const double *)ck_matrix_const_at(v, b, 0), y = *(const double *)ck_matrix_const_at(v, b, 1);
    *(double *)ck_matrix_at(f, b, 0) += x * y;
    *(double *)ck_matrix_at(f, b, 1) += x * x - y;
  }
}
int test_verification(const char *suite, const char *test)
{
  if (strcmp(suite, "JacobianVerification"))
    return 2;
  size_t b = IS("MultiBlockFiniteDifference") ? 3 : 1;
  double zero = 0;
  ck_matrix v = {0}, fd = {0};
  CHECK(ck_matrix_create(b, 2, 0, sizeof(double), &zero, &v) == CK_OK);
  double values[] = {3, 5, 0, 0, 0, 0};
  if (b == 3)
  {
    double a[] = {1, 2, .5, 4, 10, .1};
    memcpy(values, a, sizeof a);
  }
  if (IS("NearZeroVariableHandled"))
  {
    values[0] = 0;
    values[1] = 1;
  }
  for (size_t c = 0; c < b; ++c)
    CHECK(ck_matrix_set_row(&v, c, values + 2 * c, 2) == CK_OK);
  CHECK(ck_finite_difference(simple, NULL, &v, 1e-8, &fd) == CK_OK);
  for (size_t c = 0; c < b; ++c)
  {
    double expected[] = {values[c * 2 + 1], values[c * 2], 2 * values[c * 2], -1};
    for (size_t j = 0; j < 4; ++j)
      NEAR(*(double *)ck_matrix_at(&fd, c, j), expected[j], 1e-5);
  }
  if (strstr(test, "Compare") || strstr(test, "Sparsity"))
  {
    unsigned char pattern[] = {1, 1, 1, 1};
    if (IS("SparsityCompletenessFailsWhenMissingEntry"))
      pattern[2] = 0;
    ck_sparse a;
    CHECK(ck_sparse_create(2, b, 0, CK_SPARSE_ROW, pattern, sizeof(double), &zero, &a) == CK_OK);
    double expected[] = {IS("CompareFailsForWrongJacobian") ? -999 : -5, -3, -6, 1};
    for (size_t row = 0; row < 2; ++row)
      for (size_t col = 0; col < 2; ++col)
      {
        double *p = ck_sparse_at(&a, 0, row, col);
        if (p)
          *p = expected[row * 2 + col];
      }
    ck_jacobian_comparison result;
    if (strstr(test, "Sparsity"))
      CHECK(ck_check_sparsity(&a, &fd, 1e-6, &result) == CK_OK);
    else
      CHECK(ck_compare_jacobian(&a, &fd, 1e-7, 1e-7, &result) == CK_OK);
    if (strstr(test, "Fails"))
    {
      CHECK(!result.passed);
      CHECK(result.worst_row == (size_t)(strstr(test, "Sparsity") != NULL));
      CHECK(result.worst_column == 0);
    }
    else
    {
      CHECK(result.passed);
      if (strstr(test, "Compare"))
        CHECK(result.max_abs_error < 1e-5);
    }
    ck_sparse_free(&a);
  }
  ck_matrix_free(&fd);
  ck_matrix_free(&v);
  return 0;
}
