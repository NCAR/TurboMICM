/* MICM test/unit/util/test_jacobian_verification.cpp: the finite-difference
 * Jacobian check of micm/util/jacobian_verification.hpp, ported to C in
 * tmc_solver_ext.h. MICM runs it on its standard dense matrix (L = 1); the
 * multi-block case also runs at L = 4, where the last block is padded. */
#include "harness.h"
#include "tmc_solver_ext.h"

static const double fd_tolerance = 1e-5;

/* f(x, y) = [x y, x^2 - y], added to forcing for every block.
 * Jacobian: df0/dx = y, df0/dy = x, df1/dx = 2x, df1/dy = -1. */
static void simple_forcing(const tmc_dense *vars, tmc_dense *forcing, void *context)
{
  (void)context;
  for (int block = 0; block < vars->n_rows; ++block)
  {
    double x = tmc_dense_get(vars, block, 0);
    double y = tmc_dense_get(vars, block, 1);
    tmc_dense_set(forcing, block, 0, tmc_dense_get(forcing, block, 0) + x * y);
    tmc_dense_set(forcing, block, 1, tmc_dense_get(forcing, block, 1) + (x * x - y));
  }
}

static tmc_dense *fd_jacobian_at(double x, double y)
{
  tmc_dense *variables = tmc_dense_create(1, 2, 1, 0.0);
  tmc_dense_set(variables, 0, 0, x);
  tmc_dense_set(variables, 0, 1, y);
  tmc_dense *fd = tmc_finite_difference_jacobian(simple_forcing, NULL, variables, 2, TMC_FD_PERTURBATION);
  tmc_dense_free(variables);
  return fd;
}

/* A one-block 2 x 2 sparse matrix (standard ordering) with the given elements and values. */
static tmc_sparse *analytical_jacobian(const bool present[2][2], const double values[2][2])
{
  tmc_sparse_builder *builder = tmc_sparse_builder_create(2);
  tmc_sparse_builder_set_number_of_blocks(builder, 1);
  tmc_sparse_builder_set_initial_value(builder, 0.0);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      if (present[i][j])
        tmc_sparse_builder_with_element(builder, i, j, NULL);
  tmc_sparse *m = tmc_sparse_build(builder, TMC_SPARSE_CSR, 1, false);
  tmc_sparse_builder_free(builder);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      if (present[i][j])
        tmc_sparse_set(m, 0, i, j, values[i][j]);
  return m;
}

static const bool all_elements[2][2] = { { true, true }, { true, true } };

TMC_TEST(finite_difference_matches_analytical)
{
  const int n = 2;
  tmc_dense *fd = fd_jacobian_at(3.0, 5.0);
  TMC_EXPECT_EQ(fd->n_rows, 1);
  TMC_EXPECT_EQ(fd->n_cols, n * n);
  TMC_EXPECT_NEAR(tmc_dense_get(fd, 0, 0 * n + 0), 5.0, fd_tolerance);  /* df0/dx */
  TMC_EXPECT_NEAR(tmc_dense_get(fd, 0, 0 * n + 1), 3.0, fd_tolerance);  /* df0/dy */
  TMC_EXPECT_NEAR(tmc_dense_get(fd, 0, 1 * n + 0), 6.0, fd_tolerance);  /* df1/dx */
  TMC_EXPECT_NEAR(tmc_dense_get(fd, 0, 1 * n + 1), -1.0, fd_tolerance); /* df1/dy */
  tmc_dense_free(fd);
  return TMC_TEST_PASS;
}

TMC_TEST(multi_block_finite_difference)
{
  static const int widths[2] = { 1, 4 };
  for (int w = 0; w < 2; ++w)
  {
    tmc_dense *variables = tmc_dense_create(3, 2, widths[w], 0.0);
    tmc_dense_set(variables, 0, 0, 1.0);
    tmc_dense_set(variables, 0, 1, 2.0);
    tmc_dense_set(variables, 1, 0, 0.5);
    tmc_dense_set(variables, 1, 1, 4.0);
    tmc_dense_set(variables, 2, 0, 10.0);
    tmc_dense_set(variables, 2, 1, 0.1);
    tmc_dense *fd = tmc_finite_difference_jacobian(simple_forcing, NULL, variables, 2, TMC_FD_PERTURBATION);
    /* block 0: x = 1, y = 2 */
    TMC_EXPECT_NEAR(tmc_dense_get(fd, 0, 0), 2.0, fd_tolerance);
    TMC_EXPECT_NEAR(tmc_dense_get(fd, 0, 1), 1.0, fd_tolerance);
    TMC_EXPECT_NEAR(tmc_dense_get(fd, 0, 2), 2.0, fd_tolerance);
    TMC_EXPECT_NEAR(tmc_dense_get(fd, 0, 3), -1.0, fd_tolerance);
    /* block 1: x = 0.5, y = 4 */
    TMC_EXPECT_NEAR(tmc_dense_get(fd, 1, 0), 4.0, fd_tolerance);
    TMC_EXPECT_NEAR(tmc_dense_get(fd, 1, 1), 0.5, fd_tolerance);
    TMC_EXPECT_NEAR(tmc_dense_get(fd, 1, 2), 1.0, fd_tolerance);
    TMC_EXPECT_NEAR(tmc_dense_get(fd, 1, 3), -1.0, fd_tolerance);
    /* block 2: x = 10, y = 0.1 */
    TMC_EXPECT_NEAR(tmc_dense_get(fd, 2, 0), 0.1, fd_tolerance);
    TMC_EXPECT_NEAR(tmc_dense_get(fd, 2, 1), 10.0, fd_tolerance);
    TMC_EXPECT_NEAR(tmc_dense_get(fd, 2, 2), 20.0, fd_tolerance);
    TMC_EXPECT_NEAR(tmc_dense_get(fd, 2, 3), -1.0, fd_tolerance);
    tmc_dense_free(fd);
    tmc_dense_free(variables);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(compare_passes_for_correct_jacobian)
{
  tmc_dense *fd = fd_jacobian_at(3.0, 5.0);
  /* -df/dx, MICM's convention */
  const double values[2][2] = { { -5.0, -3.0 }, { -6.0, 1.0 } };
  tmc_sparse *analytical = analytical_jacobian(all_elements, values);
  tmc_jacobian_comparison result = tmc_compare_jacobian_to_finite_difference(analytical, fd, 2, TMC_FD_COMPARE_ATOL, TMC_FD_COMPARE_RTOL);
  TMC_EXPECT_TRUE(result.passed);
  TMC_EXPECT_TRUE(result.max_abs_error < fd_tolerance);
  tmc_sparse_free(analytical);
  tmc_dense_free(fd);
  return TMC_TEST_PASS;
}

TMC_TEST(compare_fails_for_wrong_jacobian)
{
  tmc_dense *fd = fd_jacobian_at(3.0, 5.0);
  const double values[2][2] = { { -999.0, -3.0 }, { -6.0, 1.0 } };
  tmc_sparse *analytical = analytical_jacobian(all_elements, values);
  tmc_jacobian_comparison result = tmc_compare_jacobian_to_finite_difference(analytical, fd, 2, TMC_FD_COMPARE_ATOL, TMC_FD_COMPARE_RTOL);
  TMC_EXPECT_FALSE(result.passed);
  TMC_EXPECT_EQ(result.worst_row, 0);
  TMC_EXPECT_EQ(result.worst_col, 0);
  tmc_sparse_free(analytical);
  tmc_dense_free(fd);
  return TMC_TEST_PASS;
}

TMC_TEST(sparsity_completeness_passes_when_complete)
{
  tmc_dense *fd = fd_jacobian_at(3.0, 5.0);
  const double values[2][2] = { { 0.0, 0.0 }, { 0.0, 0.0 } };
  tmc_sparse *full_sparsity = analytical_jacobian(all_elements, values);
  tmc_jacobian_comparison result = tmc_check_jacobian_sparsity_completeness(full_sparsity, fd, 2, TMC_FD_SPARSITY_THRESHOLD);
  TMC_EXPECT_TRUE(result.passed);
  tmc_sparse_free(full_sparsity);
  tmc_dense_free(fd);
  return TMC_TEST_PASS;
}

TMC_TEST(sparsity_completeness_fails_when_missing_entry)
{
  tmc_dense *fd = fd_jacobian_at(3.0, 5.0);
  /* (1, 0) is left out, but df1/dx = 2x = 6 */
  const bool present[2][2] = { { true, true }, { false, true } };
  const double values[2][2] = { { 0.0, 0.0 }, { 0.0, 0.0 } };
  tmc_sparse *incomplete_sparsity = analytical_jacobian(present, values);
  tmc_jacobian_comparison result = tmc_check_jacobian_sparsity_completeness(incomplete_sparsity, fd, 2, TMC_FD_SPARSITY_THRESHOLD);
  TMC_EXPECT_FALSE(result.passed);
  TMC_EXPECT_EQ(result.worst_row, 1);
  TMC_EXPECT_EQ(result.worst_col, 0);
  tmc_sparse_free(incomplete_sparsity);
  tmc_dense_free(fd);
  return TMC_TEST_PASS;
}

TMC_TEST(near_zero_variable_handled)
{
  const int n = 2;
  /* x = 0 takes a forward difference */
  tmc_dense *fd = fd_jacobian_at(0.0, 1.0);
  TMC_EXPECT_NEAR(tmc_dense_get(fd, 0, 0 * n + 0), 1.0, fd_tolerance);
  TMC_EXPECT_NEAR(tmc_dense_get(fd, 0, 0 * n + 1), 0.0, fd_tolerance);
  TMC_EXPECT_NEAR(tmc_dense_get(fd, 0, 1 * n + 0), 0.0, fd_tolerance);
  TMC_EXPECT_NEAR(tmc_dense_get(fd, 0, 1 * n + 1), -1.0, fd_tolerance);
  tmc_dense_free(fd);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("jacobian_verification:JacobianVerification.FiniteDifferenceMatchesAnalytical", finite_difference_matches_analytical),
  TMC_CASE("jacobian_verification:JacobianVerification.MultiBlockFiniteDifference", multi_block_finite_difference),
  TMC_CASE("jacobian_verification:JacobianVerification.ComparePassesForCorrectJacobian", compare_passes_for_correct_jacobian),
  TMC_CASE("jacobian_verification:JacobianVerification.CompareFailsForWrongJacobian", compare_fails_for_wrong_jacobian),
  TMC_CASE("jacobian_verification:JacobianVerification.SparsityCompletenessPassesWhenComplete", sparsity_completeness_passes_when_complete),
  TMC_CASE("jacobian_verification:JacobianVerification.SparsityCompletenessFailsWhenMissingEntry", sparsity_completeness_fails_when_missing_entry),
  TMC_CASE("jacobian_verification:JacobianVerification.NearZeroVariableHandled", near_zero_variable_handled),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_jacobian_verification(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
