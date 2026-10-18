/* The checks of MICM's test/unit/util/test_matrix_policy.hpp, shared by the
 * Matrix suite (test_matrix.c, L = 1) and the VectorMatrix<L> suite
 * (test_vector_matrix.c). Each takes the width L and returns TMC_TEST_PASS or
 * TMC_TEST_FAIL; failures are also recorded through the harness. A check that
 * returns a matrix hands it to the caller, which checks its layout and frees
 * it. */
#ifndef TMC_TEST_DENSE_POLICY_H
#define TMC_TEST_DENSE_POLICY_H

#include "harness.h"

tmc_dense *dense_policy_small_matrix(int L);
tmc_dense *dense_policy_small_const_matrix(int L);
int dense_policy_initialize_matrix(int L);
int dense_policy_initialize_const_matrix(int L);
int dense_policy_loop_over_matrix(int L);
int dense_policy_loop_over_const_matrix(int L);
tmc_dense *dense_policy_strides(int L);
int dense_policy_conversion_to_vector(int L);
int dense_policy_const_conversion_to_vector(int L);
int dense_policy_conversion_from_vector(int L);
int dense_policy_assignment_from_vector(int L);
int dense_policy_axpy(int L);
int dense_policy_for_each(int L);
int dense_policy_set_scalar(int L);
int dense_policy_max(int L);
int dense_policy_min(int L);
int dense_policy_print(int L);
int dense_policy_array_function(int L);
int dense_policy_multi_matrix_array_function(int L);
int dense_policy_vector_in_matrix_function(int L);
int dense_policy_multi_matrix_different_rows_from_creation(int L);
int dense_policy_matrix_vector_different_rows_from_creation(int L);
int dense_policy_mismatched_rows_at_invocation(int L);
int dense_policy_multiple_matrices_mismatched_rows_at_invocation(int L);
int dense_policy_wrong_column_count_at_invocation(int L);
int dense_policy_mismatched_row_dimensions(int L);
int dense_policy_mismatched_column_dimensions(int L);
int dense_policy_wrong_matrix_dimensions(int L);
int dense_policy_multiple_temporaries(int L);
int dense_policy_column_view_reuse(int L);
int dense_policy_function_reusability(int L);
int dense_policy_const_matrix_function(int L);
int dense_policy_empty_matrix_function(int L);
int dense_policy_vector_too_small(int L);
int dense_policy_vector_too_large(int L);
int dense_policy_empty_vector_non_empty_matrix(int L);
int dense_policy_non_empty_vector_empty_matrix(int L);
int dense_policy_empty_vector_empty_matrix(int L);
int dense_policy_multiple_vectors_different_sizes(int L);
int dense_policy_multiple_vectors_same_size(int L);
int dense_policy_multiple_matrices_one_vector(int L);
int dense_policy_multiple_matrices_different_rows_vector(int L);
int dense_policy_vector_length_matches_one_matrix_only(int L);
int dense_policy_const_vector(int L);
int dense_policy_mutable_vector(int L);
int dense_policy_function_reusability_with_vectors(int L);
int dense_policy_function_invocation_with_wrong_sized_vector(int L);
int dense_policy_mixed_vector_column_view_row_variable(int L);
int dense_policy_integer_vector(int L);
int dense_policy_function_with_const_signature(int L);
int dense_policy_fill(int L);
int dense_policy_copy(int L);
int dense_policy_reduce_sum(int L);
int dense_policy_reduce_max(int L);
int dense_policy_reduce_lor(int L);
int dense_policy_reduce_land(int L);
int dense_policy_reduce_strict(int L);
int dense_policy_vector_capture(int L);

#endif
