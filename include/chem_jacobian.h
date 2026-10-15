/* SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_JACOBIAN_H
#define CHEM_JACOBIAN_H
#include "chem_solver.h"
typedef void (*ChemForcingFunction)(const double *values,double *forcing,void *context);
typedef struct {
    int passed;
    double max_absolute_error,max_relative_error;
    size_t worst_cell,worst_row,worst_column;
    double worst_analytical,worst_fd;
} ChemJacobianComparison;
/* Dense buffers use [(row*n+column)*cells+cell]. Default perturbation is 1e-8. */
ChemStatus chem_finite_difference(size_t n,size_t cells,ChemForcingFunction,void *context,
                                 const double *values,double perturbation,double *jacobian);
/* Analytic values use the sparse layout and store -df/dy. Default atol=rtol=1e-7. */
ChemJacobianComparison chem_compare_jacobian(const ChemSparse *,size_t cells,const double *analytic,
                                             const double *fd,double atol,double rtol);
/* Default threshold is 1e-6. */
ChemJacobianComparison chem_check_sparsity(const ChemSparse *,size_t cells,const double *fd,double threshold);
#endif
