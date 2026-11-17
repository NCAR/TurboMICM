#include "linear.h"
#include <math.h>

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

