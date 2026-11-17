#ifndef LINEAR_H
#define LINEAR_H
#include <stddef.h>
typedef enum { CHEM_LU_MOZART, CHEM_LU_DOOLITTLE } ChemLuAlgorithm;
/* Dense row-major, no-pivot LU. Factorization is in place. */
int chem_dense_factor(double *matrix, size_t n);
int chem_dense_factor_copy(const double *matrix, double *factors, size_t n);
int chem_dense_factor_doolittle(double *matrix, size_t n);
/* Separate output must have n*n entries and not overlap the input. */
int chem_dense_factor_doolittle_copy(const double *matrix, double *factors, size_t n);
void chem_dense_solve(const double *matrix, size_t n, double *rhs);
/* Row-major matrix-vector product; output must not overlap either input. */
void chem_dense_matvec(size_t rows, size_t columns, const double *matrix,
                       const double *vector, double *output);
#endif
