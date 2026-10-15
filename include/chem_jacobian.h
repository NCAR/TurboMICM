#ifndef CHEM_JACOBIAN_H
#define CHEM_JACOBIAN_H
#include "chem_sparse.h"

/* Concentrations/forcing are species-major, with contiguous cells. The dense
   Jacobian is (row * species + column) * cells + cell. The forcing callback
   adds to its zero-initialized output, as does chem_add_forcing. */
typedef void (*ChemForcingFunction)(size_t cells, size_t species,
    const double *variables, double *forcing, void *context);
ChemStatus chem_finite_difference(size_t cells, size_t species, ChemForcingFunction,
    void *context, const double *variables, double perturbation, double *jacobian);
typedef struct {
    int passed;
    double max_absolute_error, max_relative_error;
    size_t worst_block, worst_row, worst_column;
    double worst_analytical, worst_finite_difference;
} ChemJacobianComparison;
/* Analytical sparse entries use the solver's negative-Jacobian convention;
   finite-difference entries are positive derivatives. */
ChemStatus chem_compare_jacobian(const ChemSparse *analytical, const double *finite_difference,
    double absolute_tolerance, double relative_tolerance, ChemJacobianComparison *);
ChemStatus chem_check_jacobian_sparsity(const ChemSparse *analytical, const double *finite_difference,
    double threshold, ChemJacobianComparison *);
#endif
