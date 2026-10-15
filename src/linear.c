/* Sparse in-place Mozart elimination, with the reference's pivot and update
   ordering. There is no numerical pivoting. Symbolic fill is computed once. */
#include "internal.h"
#include <math.h>

ChemStatus chem_factor(const ChemMechanism *m, size_t nc, double *matrix)
{
    for (size_t i = 0; i < m->data->species_count; ++i) {
        const double *restrict pivot = matrix + m->diagonal[i] * nc;
        for (size_t c = 0; c < nc; ++c) if (pivot[c] == 0 || !isfinite(pivot[c])) return CHEM_SINGULAR;
        for (size_t j = m->lower_start[i]; j < m->lower_start[i + 1]; ++j) {
            double *restrict dst = matrix + m->lower[j] * nc;
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) dst[c] *= 1.0 / pivot[c];
        }
        for (size_t j = m->operation_start[i]; j < m->operation_start[i + 1]; ++j) {
            ChemElimination op = m->operations[j];
            double *restrict dst = matrix + op.destination * nc;
            const double *restrict left = matrix + op.left * nc;
            const double *restrict right = matrix + op.right * nc;
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) dst[c] -= left[c] * right[c];
        }
    }
    return CHEM_OK;
}
void chem_linear_solve(const ChemMechanism *m, size_t nc, const double *matrix, double *x)
{
    size_t n = m->data->species_count;
    for (size_t i = 0; i < n; ++i) {
        double *restrict xi = x + i * nc;
        for (size_t j = m->row_start[i]; j < m->diagonal[i]; ++j) {
            const double *restrict lij = matrix + j * nc;
            const double *restrict xj = x + m->entries[j].column * nc;
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) xi[c] -= lij[c] * xj[c];
        }
    }
    for (size_t ii = n; ii > 0; --ii) {
        size_t i = ii - 1;
        double *restrict xi = x + i * nc;
        for (size_t j = m->diagonal[i] + 1; j < m->row_start[i + 1]; ++j) {
            const double *restrict uij = matrix + j * nc;
            const double *restrict xj = x + m->entries[j].column * nc;
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) xi[c] -= uij[c] * xj[c];
        }
        const double *restrict uii = matrix + m->diagonal[i] * nc;
        #pragma omp simd
        for (size_t c = 0; c < nc; ++c) xi[c] /= uii[c];
    }
}
