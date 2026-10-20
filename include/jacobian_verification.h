#ifndef JACOBIAN_VERIFICATION_H
#define JACOBIAN_VERIFICATION_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*chem_forcing_callback)(
    size_t block, const double *state, size_t species_count, double *forcing, void *context);

typedef struct {
  int passed;
  double max_absolute_error;
  double max_relative_error;
  size_t worst_block;
  size_t worst_row;
  size_t worst_column;
  double worst_analytical;
  double worst_finite_difference;
} chem_jacobian_comparison;

/* Approximate df/dy with central differences, using a forward difference when
   y-h would be negative. The forcing output is zeroed before each callback.
   Output is [block][row * species_count + column]. */
int chem_finite_difference_jacobian(
    chem_forcing_callback forcing,
    void *context,
    const double *state,
    size_t block_count,
    size_t species_count,
    double perturbation,
    double *jacobian);

/* Compare a sparse-pattern row-major matrix storing -df/dy against the dense
   finite-difference result. A nonzero pattern byte selects a matrix entry. */
int chem_compare_jacobian_to_finite_difference(
    const unsigned char *pattern,
    const double *negative_jacobian,
    const double *finite_difference_jacobian,
    size_t block_count,
    size_t species_count,
    double absolute_tolerance,
    double relative_tolerance,
    chem_jacobian_comparison *result);

/* Detect finite-difference derivatives outside the declared sparsity pattern. */
int chem_check_jacobian_sparsity(
    const unsigned char *pattern,
    const double *finite_difference_jacobian,
    size_t block_count,
    size_t species_count,
    double threshold,
    chem_jacobian_comparison *result);

#ifdef __cplusplus
}
#endif

#endif
