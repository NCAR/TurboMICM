#ifndef CHEM_VERIFY_H
#define CHEM_VERIFY_H
#include "chem_array.h"
typedef void (*ChemArrayForcing)(const ChemArray *variables, ChemArray *forcing, void *context);
typedef struct {
  int passed;
  double max_abs_error, max_rel_error;
  size_t worst_block, worst_row, worst_col;
  double worst_analytical, worst_fd;
} ChemJacobianComparison;
ChemStatus chem_finite_difference(const ChemArray *variables, size_t species, double perturbation,
                                  ChemArrayForcing function, void *context, ChemArray **result);
ChemJacobianComparison chem_compare_jacobian(const ChemArray *negative_analytical, const ChemArray *fd,
                                             size_t species, double absolute, double relative);
ChemJacobianComparison chem_check_sparsity(const ChemArray *pattern, const ChemArray *fd, size_t species, double threshold);
#endif
