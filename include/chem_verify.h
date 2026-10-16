#ifndef CHEM_VERIFY_H
#define CHEM_VERIFY_H
#include "chem.h"
typedef void (*chem_forcing_function)(void *context, const double *variables, double *forcing);
typedef struct {
    int passed;
    double max_absolute_error, max_relative_error;
    size_t worst_cell, worst_row, worst_column;
    double worst_analytical, worst_finite_difference;
} chem_jacobian_comparison;
chem_status chem_finite_difference(size_t species, size_t cells, const double *variables,
                                   chem_forcing_function function, void *context,
                                   double perturbation, double *result);
chem_jacobian_comparison chem_compare_jacobian(size_t species, size_t cells, const size_t *sparse_map,
    const double *negative_jacobian, const double *finite_difference, double atol, double rtol);
chem_jacobian_comparison chem_check_sparsity(size_t species, size_t cells, const size_t *sparse_map,
    const double *finite_difference, double threshold);
#endif
