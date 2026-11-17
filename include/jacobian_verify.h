#ifndef CHEM_JACOBIAN_VERIFY_H
#define CHEM_JACOBIAN_VERIFY_H
#include "matrix.h"
#include "sparse.h"
typedef void (*ChemForcingFunction)(size_t block, const double *variables, double *forcing, void *context);
typedef struct {
    int passed;
    double max_abs_error, max_rel_error, worst_analytical, worst_fd;
    size_t worst_block, worst_row, worst_column;
} ChemJacobianComparison;
/* fd has one row per block and n*n row-major Jacobian columns. */
int chem_finite_difference_jacobian(const ChemMatrix *variables, double perturbation, ChemForcingFunction,
                                    void *context, ChemMatrix *fd);
/* Sparse analytical entries hold the negative forcing derivative. */
int chem_compare_jacobian(const ChemSparseMatrix *, const ChemMatrix *fd, double atol, double rtol, ChemJacobianComparison *);
int chem_check_jacobian_sparsity(const ChemSparseMatrix *, const ChemMatrix *fd, double threshold, ChemJacobianComparison *);
#endif
