#ifndef CHEM_VERIFY_H
#define CHEM_VERIFY_H
#include "chem_matrix.h"

typedef struct
{
  int passed;
  double max_abs_error, max_rel_error;
  size_t worst_block, worst_row, worst_column;
  double worst_analytical, worst_fd;
} ck_jacobian_comparison;
typedef void (*ck_matrix_forcing)(const ck_matrix *variables, ck_matrix *forcing, void *context);
/* result has rows=cells, columns=species*species and stores +df/dy. */
ck_status ck_finite_difference(ck_matrix_forcing forcing, void *context, const ck_matrix *variables,
                               double perturbation, ck_matrix *result);
ck_status ck_compare_jacobian(const ck_sparse *negative_jacobian, const ck_matrix *finite_difference,
                              double atol, double rtol, ck_jacobian_comparison *result);
ck_status ck_check_sparsity(const ck_sparse *jacobian, const ck_matrix *finite_difference, double threshold,
                            ck_jacobian_comparison *result);
#endif
