#ifndef CHEM_JACOBIAN_H
#define CHEM_JACOBIAN_H

#include "chem_matrix.h"
#include "chem_sparse.h"

typedef int (*ChemForcingFunction)(const ChemMatrix *variables, ChemMatrix *forcing, void *context);
typedef struct {
    int passed;
    double max_absolute_error, max_relative_error;
    size_t worst_block, worst_row, worst_column;
    double worst_analytical, worst_finite_difference;
} ChemJacobianComparison;

/* A positive perturbation scales max(1,abs(x)). Pass a zero-initialized result;
 * the returned matrix has one row per cell and columns row*species_count+column.
 * The caller destroys it. */
int chem_finite_difference_jacobian(ChemForcingFunction forcing, void *context,
                                    const ChemMatrix *variables, double perturbation, ChemMatrix *result);
/* The sparse analytical matrix stores the negative Jacobian. */
int chem_compare_jacobian(const ChemSparseMatrix *analytical, const ChemMatrix *finite_difference,
                          double absolute_tolerance, double relative_tolerance, ChemJacobianComparison *result);
int chem_check_jacobian_sparsity(const ChemSparseMatrix *analytical, const ChemMatrix *finite_difference,
                                 double threshold, ChemJacobianComparison *result);

#endif
