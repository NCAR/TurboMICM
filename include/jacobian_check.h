#ifndef JACOBIAN_CHECK_H
#define JACOBIAN_CHECK_H
#include <stddef.h>

/* The callback adds its forcing to an initially zeroed output vector. */
typedef void (*ChemForcingCallback)(void *context, size_t cell, const double *state, double *forcing);
typedef struct {
    int passed;
    double max_abs_error, max_rel_error;
    size_t worst_cell, worst_row, worst_col;
    double worst_analytical, worst_fd;
} ChemJacobianComparison;
/* states: cells*n; output: cells*n*n, row-major. Returns 0 on success. */
int chem_finite_difference_jacobian(size_t cells, size_t n, const double *states,
        ChemForcingCallback forcing, void *context, double perturbation, double *output);
/* pattern has n*n entries shared by all cells; analytical contains -df/dy. */
ChemJacobianComparison chem_compare_jacobian(size_t cells, size_t n, const unsigned char *pattern,
        const double *analytical, const double *fd, double atol, double rtol);
ChemJacobianComparison chem_check_jacobian_sparsity(size_t cells, size_t n, const unsigned char *pattern,
        const double *fd, double threshold);
#endif
