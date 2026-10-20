#include "jacobian_verification.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

static int simple_forcing(size_t block, const double *state, size_t n, double *forcing, void *context) {
  (void)block;
  (void)context;
  if (n != 2) return -1;
  forcing[0] += state[0] * state[1];
  forcing[1] += state[0] * state[0] - state[1];
  return 0;
}

static int nonnegative_forcing(size_t block, const double *state, size_t n, double *forcing, void *context) {
  for (size_t i = 0; i < n; ++i)
    if (state[i] < 0.0) return -1;
  return simple_forcing(block, state, n, forcing, context);
}

static int close_to(double actual, double expected, double tolerance) {
  return isfinite(actual) && fabs(actual - expected) <= tolerance;
}

static int finite_difference_matches(void) {
  const double variables[] = { 3.0, 5.0 };
  double jacobian[4];
  if (chem_finite_difference_jacobian(simple_forcing, NULL, variables, 1, 2, 1.0e-8, jacobian) != 0)
    return 0;
  const double expected[] = { 5.0, 3.0, 6.0, -1.0 };
  for (size_t i = 0; i < 4; ++i)
    if (!close_to(jacobian[i], expected[i], 1.0e-5)) return 0;
  return 1;
}

static int multi_block_finite_difference(void) {
  const double variables[] = { 1.0, 2.0, 0.5, 4.0, 10.0, 0.1 };
  const double expected[] = { 2.0, 1.0, 2.0, -1.0, 4.0, 0.5, 1.0, -1.0,
                              0.1, 10.0, 20.0, -1.0 };
  double jacobian[12];
  if (chem_finite_difference_jacobian(simple_forcing, NULL, variables, 3, 2, 1.0e-8, jacobian) != 0)
    return 0;
  for (size_t i = 0; i < 12; ++i)
    if (!close_to(jacobian[i], expected[i], 1.0e-5)) return 0;
  return 1;
}

static int compare_passes_for_correct_jacobian(void) {
  const double variables[] = { 3.0, 5.0 };
  const double negative_jacobian[] = { -5.0, -3.0, -6.0, 1.0 };
  const unsigned char pattern[] = { 1, 1, 1, 1 };
  double fd[4];
  chem_jacobian_comparison result;
  if (chem_finite_difference_jacobian(simple_forcing, NULL, variables, 1, 2, 1.0e-8, fd) != 0 ||
      chem_compare_jacobian_to_finite_difference(pattern, negative_jacobian, fd, 1, 2,
                                                 1.0e-7, 1.0e-7, &result) != 0)
    return 0;
  return result.passed && result.max_absolute_error < 1.0e-7;
}

static int compare_fails_for_wrong_jacobian(void) {
  const double variables[] = { 3.0, 5.0 };
  const double negative_jacobian[] = { -999.0, -3.0, -6.0, 1.0 };
  const unsigned char pattern[] = { 1, 1, 1, 1 };
  double fd[4];
  chem_jacobian_comparison result;
  if (chem_finite_difference_jacobian(simple_forcing, NULL, variables, 1, 2, 1.0e-8, fd) != 0 ||
      chem_compare_jacobian_to_finite_difference(pattern, negative_jacobian, fd, 1, 2,
                                                 1.0e-7, 1.0e-7, &result) != 0)
    return 0;
  return !result.passed && result.worst_row == 0 && result.worst_column == 0;
}

static int sparsity_completeness(int expect_missing) {
  const double variables[] = { 3.0, 5.0 };
  const unsigned char full[] = { 1, 1, 1, 1 };
  const unsigned char missing[] = { 1, 1, 0, 1 };
  const unsigned char *pattern = expect_missing ? missing : full;
  double fd[4];
  chem_jacobian_comparison result;
  if (chem_finite_difference_jacobian(simple_forcing, NULL, variables, 1, 2, 1.0e-8, fd) != 0 ||
      chem_check_jacobian_sparsity(pattern, fd, 1, 2, 1.0e-6, &result) != 0)
    return 0;
  if (!expect_missing) return result.passed;
  return !result.passed && result.worst_row == 1 && result.worst_column == 0;
}

static int near_zero_variable(void) {
  const double variables[] = { 0.0, 1.0 };
  double jacobian[4];
  if (chem_finite_difference_jacobian(nonnegative_forcing, NULL, variables, 1, 2, 1.0e-8, jacobian) != 0)
    return 0;
  const double expected[] = { 1.0, 0.0, 0.0, -1.0 };
  for (size_t i = 0; i < 4; ++i)
    if (!close_to(jacobian[i], expected[i], 1.0e-5)) return 0;
  return 1;
}

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  if (strcmp(argv[1], "FiniteDifferenceMatchesAnalytical") == 0) return finite_difference_matches() ? 0 : 1;
  if (strcmp(argv[1], "MultiBlockFiniteDifference") == 0) return multi_block_finite_difference() ? 0 : 1;
  if (strcmp(argv[1], "ComparePassesForCorrectJacobian") == 0) return compare_passes_for_correct_jacobian() ? 0 : 1;
  if (strcmp(argv[1], "CompareFailsForWrongJacobian") == 0) return compare_fails_for_wrong_jacobian() ? 0 : 1;
  if (strcmp(argv[1], "SparsityCompletenessPassesWhenComplete") == 0) return sparsity_completeness(0) ? 0 : 1;
  if (strcmp(argv[1], "SparsityCompletenessFailsWhenMissingEntry") == 0) return sparsity_completeness(1) ? 0 : 1;
  if (strcmp(argv[1], "NearZeroVariableHandled") == 0) return near_zero_variable() ? 0 : 1;
  return 2;
}
