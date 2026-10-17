/* MICM test/unit/util/test_view_category_sfinae.cpp: four cases are
 * compile-time checks of C++ type traits and concepts (static_assert only),
 * with nothing to run in C. The fifth uses a plain vector where MICM accepts
 * vector-like arguments; in C that is a padded vector passed to an array
 * function. */
#include "harness.h"

#include "tmc_dense_ext.h"

TMC_TEST(vector_types)
{
  TMC_NA("compile-time only: static_asserts that std::vector and std::array have no view category and satisfy VectorLike; C has no type traits");
}

TMC_TEST(matrix_views)
{
  TMC_NA("compile-time only: static_asserts on the view-category tags of Matrix and VectorMatrix column views; C has no type traits");
}

TMC_TEST(sparse_matrix_views)
{
  TMC_NA("compile-time only: static_asserts on the view-category tags of SparseMatrix block views; C has no type traits");
}

TMC_TEST(matrix_types_excluded)
{
  TMC_NA("compile-time only: static_asserts that matrix types do not satisfy the VectorLike concept; C has no concepts");
}

static void sum_values(double *const *x, void *accumulator, void *ctx)
{
  (void)ctx;
  *(double *)accumulator += *x[0];
}

typedef struct vector_sums
{
  double strict, all;
} vector_sums;

/* sums vector argument 0, without and with its padding lanes */
static void sum_vector_kernel(tmc_array_group *g, void *ctx)
{
  vector_sums *sums = ctx;
  double *lanes[] = { tmc_array_group_vector(g, 0) };
  tmc_array_group_reduce_strict(g, sum_values, NULL, 1, lanes, &sums->strict);
  tmc_array_group_reduce(g, sum_values, NULL, 1, lanes, &sums->all);
}

TMC_TEST(runtime_vector_use)
{
  const double values[] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
  /* L = 1 is a plain vector (MICM's std::vector); L = 4 adds padding, which
   * the strict reduction leaves out and the full one sees as 0, the identity
   * of the sum */
  const int widths[] = { 1, 4 };
  for (int w = 0; w < 2; ++w)
  {
    tmc_padded_vector *v = tmc_padded_vector_from_values(values, 5, widths[w]);
    vector_sums sums = { 0.0, 0.0 };
    tmc_error err;
    tmc_error_clear(&err);
    tmc_array_arg args[] = { tmc_array_const_vector(v) };
    tmc_array_function *f = tmc_array_function_create(sum_vector_kernel, &sums, 1, args, &err);
    TMC_ASSERT_TRUE(f != NULL);
    TMC_EXPECT_EQ(tmc_array_function_apply(f, 1, args, &err), TMC_OK);
    TMC_EXPECT_DOUBLE_EQ(sums.strict, 15.0);
    TMC_EXPECT_DOUBLE_EQ(sums.all, 15.0);
    tmc_array_function_free(f);
    tmc_padded_vector_free(v);
  }
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("view_category_sfinae:ViewCategorySFINAE.VectorTypes", vector_types),
  TMC_CASE("view_category_sfinae:ViewCategorySFINAE.MatrixViews", matrix_views),
  TMC_CASE("view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews", sparse_matrix_views),
  TMC_CASE("view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded", matrix_types_excluded),
  TMC_CASE("view_category_sfinae:ViewCategorySFINAE.RuntimeVectorUse", runtime_vector_use),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_view_category_sfinae(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
