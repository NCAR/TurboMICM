/* The checks of MICM's test/unit/util/test_sparse_matrix_policy.hpp, shared
 * by the four sparse matrix suites: standard compressed row and column
 * ordering (L = 1) and vector compressed row and column ordering (width L).
 * Each check takes the element numbering and the width and returns
 * TMC_TEST_PASS or TMC_TEST_FAIL; failures are also recorded through the
 * harness. A check that returns a matrix hands it to the caller, which checks
 * its data order and frees it. */
#ifndef TMC_TEST_SPARSE_POLICY_H
#define TMC_TEST_SPARSE_POLICY_H

#include "harness.h"

tmc_sparse *sparse_policy_zero_matrix(tmc_sparse_ordering ordering, int L);
tmc_sparse *sparse_policy_const_zero_matrix(tmc_sparse_ordering ordering, int L);
tmc_sparse *sparse_policy_single_block_matrix(tmc_sparse_ordering ordering, int L);
tmc_sparse *sparse_policy_const_single_block_matrix(tmc_sparse_ordering ordering, int L);
tmc_sparse *sparse_policy_multi_block_matrix(tmc_sparse_ordering ordering, int L);
int sparse_policy_set_scalar(tmc_sparse_ordering ordering, int L);
int sparse_policy_add_to_diagonal(tmc_sparse_ordering ordering, int L);
int sparse_policy_print(tmc_sparse_ordering ordering, int L);
int sparse_policy_print_non_zero(tmc_sparse_ordering ordering, int L);
int sparse_policy_builder_bad_configuration(void);

/* Function */
int sparse_policy_array_function(tmc_sparse_ordering ordering, int L);
int sparse_policy_multi_matrix_array_function(tmc_sparse_ordering ordering, int L);
int sparse_policy_mismatched_block_dimensions(tmc_sparse_ordering ordering, int L);
int sparse_policy_mismatched_element_dimensions(tmc_sparse_ordering ordering, int L);
int sparse_policy_multiple_temporaries(tmc_sparse_ordering ordering, int L);
int sparse_policy_block_view_reuse(tmc_sparse_ordering ordering, int L);
int sparse_policy_function_reusability(tmc_sparse_ordering ordering, int L);
int sparse_policy_two_sparse_matrices_different_structure(tmc_sparse_ordering ordering, int L);
/* a sparse matrix of width L with a dense matrix of width dense_L (MICM's
 * Matrix at 1, its VectorMatrix<L> otherwise) */
int sparse_policy_sparse_and_dense_matrix_function(tmc_sparse_ordering ordering, int L, int dense_L);
int sparse_policy_incompatible_ordering(tmc_sparse_ordering ordering, int L);
int sparse_policy_incompatible_vector_ordering(tmc_sparse_ordering ordering, int L, int dense_L);
int sparse_policy_incompatible_sparse_ordering(tmc_sparse_ordering ordering, int L, int other_L);
int sparse_policy_const_sparse_matrix_function(tmc_sparse_ordering ordering, int L);
int sparse_policy_empty_sparse_matrix_function(tmc_sparse_ordering ordering, int L);
int sparse_policy_multiple_sparse_matrices_different_blocks_from_creation(tmc_sparse_ordering ordering, int L);
int sparse_policy_sparse_matrix_vector_different_blocks_from_creation(tmc_sparse_ordering ordering, int L);
int sparse_policy_mismatched_blocks_at_invocation(tmc_sparse_ordering ordering, int L);
int sparse_policy_multiple_sparse_matrices_mismatched_blocks_at_invocation(tmc_sparse_ordering ordering, int L);

/* vectors */
int sparse_policy_vector_in_sparse_matrix_function(tmc_sparse_ordering ordering, int L);
int sparse_policy_vector_too_small(tmc_sparse_ordering ordering, int L);
int sparse_policy_vector_too_large(tmc_sparse_ordering ordering, int L);
int sparse_policy_empty_vector_non_empty_sparse_matrix(tmc_sparse_ordering ordering, int L);
int sparse_policy_non_empty_vector_empty_sparse_matrix(tmc_sparse_ordering ordering, int L);
int sparse_policy_empty_vector_empty_sparse_matrix(tmc_sparse_ordering ordering, int L);
int sparse_policy_multiple_vectors_different_sizes(tmc_sparse_ordering ordering, int L);
int sparse_policy_multiple_vectors_same_size(tmc_sparse_ordering ordering, int L);
int sparse_policy_multiple_sparse_matrices_one_vector(tmc_sparse_ordering ordering, int L);
int sparse_policy_multiple_sparse_matrices_different_blocks_vector(tmc_sparse_ordering ordering, int L);
int sparse_policy_vector_size_matches_one_sparse_matrix_only(tmc_sparse_ordering ordering, int L);
int sparse_policy_const_vector_sparse(tmc_sparse_ordering ordering, int L);
int sparse_policy_mutable_vector_sparse(tmc_sparse_ordering ordering, int L);
int sparse_policy_function_reusability_with_vectors_sparse(tmc_sparse_ordering ordering, int L);
int sparse_policy_function_invocation_with_wrong_sized_vector_sparse(tmc_sparse_ordering ordering, int L);
int sparse_policy_mixed_vector_block_view_block_variable(tmc_sparse_ordering ordering, int L);
int sparse_policy_integer_vector_sparse(tmc_sparse_ordering ordering, int L);
int sparse_policy_function_with_const_signature_sparse(tmc_sparse_ordering ordering, int L);
int sparse_policy_get_block_view_by_vector_index(tmc_sparse_ordering ordering, int L);
int sparse_policy_fill(tmc_sparse_ordering ordering, int L);
/* the dense matrix of the cross-type copies has the sparse matrix's width */
int sparse_policy_copy(tmc_sparse_ordering ordering, int L);

#endif
