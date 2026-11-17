#include "linear.h"
#include <math.h>
#include <stdint.h>
#include <string.h>

/* Dense form of the no-pivot Mozart elimination: scale a column, then
   update each trailing column. The factor is stored in place. */
int chem_dense_factor(double *a, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        if (a[i*n+i] == 0 || !isfinite(a[i*n+i])) return 0;
        double inverse = 1.0 / a[i*n+i];
        for (size_t j = i+1; j < n; ++j) a[j*n+i] *= inverse;
        for (size_t k = i+1; k < n; ++k)
            for (size_t j = i+1; j < n; ++j) a[j*n+k] -= a[j*n+i]*a[i*n+k];
    }
    return 1;
}

void chem_dense_solve(const double *a, size_t n, double *b)
{
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < i; ++j) b[i] -= a[i*n+j]*b[j];
    for (size_t i = n; i-- > 0;) {
        for (size_t j = i+1; j < n; ++j) b[i] -= a[i*n+j]*b[j];
        b[i] /= a[i*n+i];
    }
}


int chem_dense_factor_doolittle(double *a, size_t n)
{
    for (size_t i=0; i<n; ++i) {
        for (size_t k=i; k<n; ++k)
            for (size_t j=0; j<i; ++j) a[i*n+k] -= a[i*n+j]*a[j*n+k];
        if (a[i*n+i] == 0 || !isfinite(a[i*n+i])) return 0;
        for (size_t k=i+1; k<n; ++k) {
            for (size_t j=0; j<i; ++j) a[k*n+i] -= a[k*n+j]*a[j*n+i];
            a[k*n+i] /= a[i*n+i];
        }
    }
    return 1;
}

int chem_dense_factor_doolittle_copy(const double *matrix, double *factors, size_t n)
{
    if (!matrix || !factors || !n || n > SIZE_MAX/n || n*n > SIZE_MAX/sizeof(double)) return 0;
    memcpy(factors,matrix,n*n*sizeof(double));
    return chem_dense_factor_doolittle(factors,n);
}

int chem_dense_factor_copy(const double *matrix, double *factors, size_t n)
{
    if (!matrix || !factors || !n || n > SIZE_MAX/n || n*n > SIZE_MAX/sizeof(double)) return 0;
    memcpy(factors,matrix,n*n*sizeof(double));
    return chem_dense_factor(factors,n);
}
