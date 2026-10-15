/* Sparse in-place Mozart elimination, with the reference's pivot and update
   ordering. There is no numerical pivoting. Symbolic fill is computed once. */
#include "internal.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

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

static ChemStatus factor_doolittle(const ChemMechanism *m, size_t nc, double *matrix)
{
    size_t n = m->data->species_count;
    for (size_t i = 0; i < n; ++i) {
        for (size_t k = i; k < n; ++k) {
            size_t ik = m->index[i * n + k];
            if (ik == SIZE_MAX) continue;
            double *restrict dst = matrix + ik * nc;
            for (size_t j = 0; j < i; ++j) {
                size_t ij = m->index[i * n + j], jk = m->index[j * n + k];
                if (ij == SIZE_MAX || jk == SIZE_MAX) continue;
                const double *restrict a = matrix + ij * nc, *restrict b = matrix + jk * nc;
                #pragma omp simd
                for (size_t c = 0; c < nc; ++c) dst[c] -= a[c] * b[c];
            }
        }
        const double *restrict pivot = matrix + m->diagonal[i] * nc;
        for (size_t c = 0; c < nc; ++c) if (pivot[c] == 0 || !isfinite(pivot[c])) return CHEM_SINGULAR;
        for (size_t k = i + 1; k < n; ++k) {
            size_t ki = m->index[k * n + i];
            if (ki == SIZE_MAX) continue;
            double *restrict dst = matrix + ki * nc;
            for (size_t j = 0; j < i; ++j) {
                size_t kj = m->index[k * n + j], ji = m->index[j * n + i];
                if (kj == SIZE_MAX || ji == SIZE_MAX) continue;
                const double *restrict a = matrix + kj * nc, *restrict b = matrix + ji * nc;
                #pragma omp simd
                for (size_t c = 0; c < nc; ++c) dst[c] -= a[c] * b[c];
            }
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) dst[c] /= pivot[c];
        }
    }
    return CHEM_OK;
}
ChemLinearSolver chem_linear_mozart(void) { return (ChemLinearSolver){chem_factor, chem_linear_solve}; }
ChemLinearSolver chem_linear_doolittle(void) { return (ChemLinearSolver){factor_doolittle, chem_linear_solve}; }

ChemStatus chem_factor_separate(const ChemMechanism *m, size_t nc,
    const double *matrix, double *lower, double *upper, ChemLinearSolver solver)
{
    if (!m || !nc || !matrix || !lower || !upper || lower == upper
        || matrix == lower || matrix == upper || !solver.factor
        || m->nnz > SIZE_MAX / nc / sizeof(double)) return CHEM_INVALID;
    size_t length = m->nnz * nc;
    double *work = malloc(length * sizeof(*work));
    if (!work) return CHEM_MEMORY;
    memcpy(work, matrix, length * sizeof(*work));
    ChemStatus status = solver.factor(m, nc, work);
    if (status == CHEM_OK) {
        memset(lower, 0, length * sizeof(*lower));
        memset(upper, 0, length * sizeof(*upper));
        for (size_t row = 0; row < m->data->species_count; ++row)
            for (size_t j = m->row_start[row]; j < m->row_start[row + 1]; ++j) {
                size_t column = m->entries[j].column;
                for (size_t c = 0; c < nc; ++c) {
                    if (column < row) lower[j * nc + c] = work[j * nc + c];
                    else upper[j * nc + c] = work[j * nc + c];
                    if (column == row) lower[j * nc + c] = 1.0;
                }
            }
    }
    free(work);
    return status;
}

void chem_linear_solve_separate(const ChemMechanism *m, size_t nc,
    const double *lower, const double *upper, double *x)
{
    size_t n = m->data->species_count;
    for (size_t i = 0; i < n; ++i) {
        double *restrict xi = x + i * nc;
        for (size_t j = m->row_start[i]; j < m->diagonal[i]; ++j) {
            const double *restrict lij = lower + j * nc;
            const double *restrict xj = x + m->entries[j].column * nc;
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) xi[c] -= lij[c] * xj[c];
        }
    }
    for (size_t ii = n; ii > 0; --ii) {
        size_t i = ii - 1;
        double *restrict xi = x + i * nc;
        for (size_t j = m->diagonal[i] + 1; j < m->row_start[i + 1]; ++j) {
            const double *restrict uij = upper + j * nc;
            const double *restrict xj = x + m->entries[j].column * nc;
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) xi[c] -= uij[c] * xj[c];
        }
        const double *restrict uii = upper + m->diagonal[i] * nc;
        #pragma omp simd
        for (size_t c = 0; c < nc; ++c) xi[c] /= uii[c];
    }
}

ChemStatus chem_diagonal_markowitz(size_t n, size_t count,
    const ChemMatrixIndex *entries, size_t *permutation)
{
    if (!n || !permutation || (count && !entries) || n > SIZE_MAX / n) return CHEM_INVALID;
    unsigned char *pattern = calloc(n * n, 1), *alive = malloc(n);
    if (!pattern || !alive) { free(pattern); free(alive); return CHEM_MEMORY; }
    memset(alive, 1, n);
    for (size_t k = 0; k < count; ++k) {
        if (entries[k].row >= n || entries[k].column >= n) {
            free(pattern); free(alive); return CHEM_INVALID;
        }
        pattern[entries[k].row * n + entries[k].column] = 1;
    }
    for (size_t step = 0; step < n; ++step) {
        size_t pivot = SIZE_MAX, best = SIZE_MAX;
        for (size_t i = 0; i < n; ++i) if (alive[i]) {
            size_t rows = 0, columns = 0;
            for (size_t j = 0; j < n; ++j) if (alive[j]) {
                rows += pattern[i * n + j]; columns += pattern[j * n + i];
            }
            size_t cost = (rows - 1) * (columns - 1);
            if (pivot == SIZE_MAX || cost < best) { pivot = i; best = cost; }
        }
        permutation[step] = pivot;
        alive[pivot] = 0;
        for (size_t i = 0; i < n; ++i) if (alive[i] && pattern[i * n + pivot])
            for (size_t j = 0; j < n; ++j) if (alive[j] && pattern[pivot * n + j]) pattern[i * n + j] = 1;
    }
    free(pattern); free(alive);
    return CHEM_OK;
}
ChemStatus chem_set_linear_solver(ChemState *s, ChemLinearSolver linear)
{
    if (!s || !linear.factor || !linear.solve) return CHEM_INVALID;
    s->linear = linear;
    return CHEM_OK;
}

void chem_alpha_minus_jacobian(const ChemMechanism *m, size_t nc, double *matrix, double alpha)
{
    for (size_t i = 0; i < m->data->species_count; ++i) if (!m->algebraic[i]) {
        double *restrict diagonal = matrix + m->diagonal[i] * nc;
        #pragma omp simd
        for (size_t c = 0; c < nc; ++c) diagonal[c] += alpha;
    }
}
