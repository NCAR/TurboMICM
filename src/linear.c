#include "internal.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void *allocate(size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size) return NULL;
    return calloc(n ? n : 1, size);
}

/* Symbolic elimination is independent of the concentrations. It is kept
 * separate from numeric factorization so another linear solver can replace it. */
int chem_make_symbolic(ChemState *s)
{
    const ChemMechanism *m = s->mechanism;
    size_t n = m->species_count;
    if (n > INT_MAX || n > SIZE_MAX / n || n*n > INT_MAX) return CHEM_INVALID;
    unsigned char *pattern = allocate(n*n, sizeof(*pattern));
    if (!pattern) return CHEM_MEMORY;
    for (size_t i = 0; i < n; ++i) pattern[i*n+i] = 1;
    for (size_t r = 0; r < m->reaction_count; ++r) {
        const ChemReaction *reaction = &m->reactions[r];
        for (size_t k = 0; k < reaction->reactant_count; ++k) {
            int col = m->reactants[reaction->reactant_start+k];
            if (col < 0 || (size_t)col >= n) { free(pattern); return CHEM_INVALID; }
            for (size_t j = 0; j < reaction->reactant_count; ++j) {
                int row = m->reactants[reaction->reactant_start+j];
                if (row < 0 || (size_t)row >= n) { free(pattern); return CHEM_INVALID; }
                pattern[(size_t)row*n+(size_t)col] = 1;
            }
            for (size_t j = 0; j < reaction->product_count; ++j) {
                int row = m->products[reaction->product_start+j].species;
                if (row < 0 || (size_t)row >= n) { free(pattern); return CHEM_INVALID; }
                pattern[(size_t)row*n+(size_t)col] = 1;
            }
        }
    }
    for (size_t i = 0; i < n; ++i)
        for (size_t k = i+1; k < n; ++k)
            if (pattern[i*n+k])
                for (size_t j = i+1; j < n; ++j)
                    if (pattern[j*n+i]) pattern[j*n+k] = 1;

    s->index = allocate(n*n, sizeof(*s->index));
    s->row_start = allocate(n+1, sizeof(*s->row_start));
    s->pivots = allocate(n, sizeof(*s->pivots));
    if (!s->index || !s->row_start || !s->pivots) { free(pattern); return CHEM_MEMORY; }
    for (size_t i = 0; i < n*n; ++i) s->index[i] = -1;
    for (size_t i = 0; i < n; ++i) {
        s->row_start[i] = s->nonzeros;
        for (size_t j = 0; j < n; ++j)
            if (pattern[i*n+j]) s->index[i*n+j] = (int)s->nonzeros++;
    }
    s->row_start[n] = s->nonzeros;
    s->columns = allocate(s->nonzeros, sizeof(*s->columns));
    if (!s->columns) { free(pattern); return CHEM_MEMORY; }
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (pattern[i*n+j]) s->columns[s->index[i*n+j]] = (int)j;

    size_t nl = 0, nu = 0, nz = 0;
    for (size_t i = 0; i < n; ++i) {
        size_t count_l = 0, count_u = 0;
        for (size_t j = i+1; j < n; ++j) count_l += pattern[j*n+i];
        for (size_t k = i+1; k < n; ++k) count_u += pattern[i*n+k];
        if (count_l && count_u > (SIZE_MAX - nz) / count_l) { free(pattern); return CHEM_MEMORY; }
        nz += count_l * count_u; nl += count_l; nu += count_u;
    }
    s->lower = allocate(nl, sizeof(*s->lower));
    s->upper = allocate(nu, sizeof(*s->upper));
    s->updates = allocate(nz, sizeof(*s->updates));
    if (!s->lower || !s->upper || !s->updates) { free(pattern); return CHEM_MEMORY; }
    nl = nu = nz = 0;
    for (size_t i = 0; i < n; ++i) {
        ChemPivot *pivot = &s->pivots[i];
        pivot->diagonal = s->index[i*n+i];
        pivot->first_lower = nl;
        for (size_t j = i+1; j < n; ++j)
            if (pattern[j*n+i]) { s->lower[nl++] = s->index[j*n+i]; ++pivot->lowers; }
        pivot->first_upper = nu;
        for (size_t k = i+1; k < n; ++k) if (pattern[i*n+k]) {
            ChemUpper *u = &s->upper[nu++];
            u->position = s->index[i*n+k];
            u->first_update = nz;
            for (size_t j = i+1; j < n; ++j) if (pattern[j*n+i]) {
                s->updates[nz++] = (ChemUpdate){s->index[j*n+k], s->index[j*n+i]};
                ++u->updates;
            }
            ++pivot->uppers;
        }
    }
    free(pattern);
    return CHEM_OK;
}

int chem_factor(ChemState *s, double alpha)
{
    size_t n = s->mechanism->species_count, cells = s->cells;
    double *a = s->matrix, *inverse = s->scratch;
    for (size_t i = 0; i < n; ++i) {
        double *diag = a + (size_t)s->pivots[i].diagonal*cells;
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) diag[c] += alpha;
    }
    for (size_t i = 0; i < n; ++i) {
        ChemPivot *p = &s->pivots[i];
        double *diag = a + (size_t)p->diagonal*cells;
        for (size_t c = 0; c < cells; ++c)
            if (diag[c] == 0.0 || !isfinite(diag[c])) return CHEM_SINGULAR;
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) inverse[c] = 1.0 / diag[c];
        for (size_t j = 0; j < p->lowers; ++j) {
            double *lower = a + (size_t)s->lower[p->first_lower+j]*cells;
            #pragma omp simd
            for (size_t c = 0; c < cells; ++c) lower[c] *= inverse[c];
        }
        for (size_t k = 0; k < p->uppers; ++k) {
            ChemUpper *u = &s->upper[p->first_upper+k];
            const double *upper = a + (size_t)u->position*cells;
            for (size_t j = 0; j < u->updates; ++j) {
                ChemUpdate *op = &s->updates[u->first_update+j];
                double *target = a + (size_t)op->target*cells;
                const double *lower = a + (size_t)op->lower*cells;
                #pragma omp simd
                for (size_t c = 0; c < cells; ++c) target[c] -= lower[c] * upper[c];
            }
        }
    }
    return CHEM_OK;
}

void chem_substitute(const ChemState *s, double *x)
{
    size_t n = s->mechanism->species_count, cells = s->cells;
    const double *a = s->matrix;
    for (size_t i = 0; i < n; ++i) {
        double *xi = x + i*cells;
        for (size_t k = s->row_start[i]; k < s->row_start[i+1] && (size_t)s->columns[k] < i; ++k) {
            const double *lij = a + k*cells, *xj = x + (size_t)s->columns[k]*cells;
            #pragma omp simd
            for (size_t c = 0; c < cells; ++c) xi[c] -= lij[c]*xj[c];
        }
    }
    for (size_t i = n; i-- > 0;) {
        double *xi = x + i*cells;
        size_t diagonal = (size_t)s->pivots[i].diagonal;
        for (size_t k = diagonal+1; k < s->row_start[i+1]; ++k) {
            const double *uij = a + k*cells, *xj = x + (size_t)s->columns[k]*cells;
            #pragma omp simd
            for (size_t c = 0; c < cells; ++c) xi[c] -= uij[c]*xj[c];
        }
        const double *diag = a + diagonal*cells;
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) xi[c] /= diag[c];
    }
}
