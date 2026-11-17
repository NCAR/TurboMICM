#ifndef LINEAR_H
#define LINEAR_H
#include <stddef.h>
/* Dense row-major, no-pivot Mozart LU. Factorization is in place. */
int chem_dense_factor(double *matrix, size_t n);
void chem_dense_solve(const double *matrix, size_t n, double *rhs);
#endif
