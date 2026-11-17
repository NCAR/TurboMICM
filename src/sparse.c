#include "sparse.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

struct ChemSparsePlan {
    size_t n;
    size_t *row_offsets, *row_columns;
    size_t *col_offsets, *col_rows;
};

void chem_sparse_plan_destroy(ChemSparsePlan *p)
{
    if (p) {
        free(p->row_offsets); free(p->row_columns);
        free(p->col_offsets); free(p->col_rows); free(p);
    }
}
size_t chem_sparse_plan_size(const ChemSparsePlan *p) { return p ? p->n : 0; }

ChemSparsePlan *chem_sparse_plan_create(const ChemMechanism *m)
{
    if (!m || !m->species_count || m->species_count > SIZE_MAX/m->species_count ||
        m->species_count >= SIZE_MAX/sizeof(size_t)) return NULL;
    size_t n = m->species_count;
    unsigned char *pattern = calloc(n*n, 1);
    ChemSparsePlan *p = calloc(1, sizeof(*p));
    if (!pattern || !p) { free(pattern); free(p); return NULL; }
    p->n = n;
    for (size_t i = 0; i < n; ++i) pattern[i*n+i] = 1;
    for (size_t r = 0; r < m->reaction_count; ++r) {
        const ChemReaction *reaction = m->reactions+r;
        for (size_t c = 0; c < reaction->reactant_count; ++c) {
            size_t col = m->reactants[reaction->reactant_begin+c];
            if (col >= n) goto fail;
            for (size_t j = 0; j < reaction->reactant_count; ++j) {
                size_t row = m->reactants[reaction->reactant_begin+j];
                if (row >= n) goto fail;
                pattern[row*n+col] = 1;
            }
            for (size_t j = 0; j < reaction->product_count; ++j) {
                size_t row = m->products[reaction->product_begin+j].species;
                if (row >= n) goto fail;
                pattern[row*n+col] = 1;
            }
        }
    }
    /* Complete the structural fill for no-pivot elimination. */
    for (size_t pivot = 0; pivot < n; ++pivot)
        for (size_t row = pivot+1; row < n; ++row)
            if (pattern[row*n+pivot])
                for (size_t col = pivot+1; col < n; ++col)
                    if (pattern[pivot*n+col]) pattern[row*n+col] = 1;
    size_t entries = 0;
    for (size_t i = 0; i < n*n; ++i) entries += pattern[i] != 0;
    if (entries > SIZE_MAX/sizeof(size_t)) goto fail;
    p->row_offsets = malloc((n+1)*sizeof(size_t));
    p->col_offsets = malloc((n+1)*sizeof(size_t));
    p->row_columns = malloc(entries*sizeof(size_t));
    p->col_rows = malloc(entries*sizeof(size_t));
    if (!p->row_offsets || !p->col_offsets || !p->row_columns || !p->col_rows) goto fail;
    size_t row_pos = 0, col_pos = 0;
    for (size_t i = 0; i < n; ++i) {
        p->row_offsets[i] = row_pos; p->col_offsets[i] = col_pos;
        for (size_t j = 0; j < n; ++j) {
            if (pattern[i*n+j]) p->row_columns[row_pos++] = j;
            if (pattern[j*n+i]) p->col_rows[col_pos++] = j;
        }
    }
    p->row_offsets[n] = row_pos; p->col_offsets[n] = col_pos;
    free(pattern); return p;
fail:
    free(pattern); chem_sparse_plan_destroy(p); return NULL;
}

int chem_sparse_factor(const ChemSparsePlan *p, double *a)
{
    size_t n = p->n;
    for (size_t i = 0; i < n; ++i) {
        if (a[i*n+i] == 0 || !isfinite(a[i*n+i])) return 0;
        double inverse = 1.0/a[i*n+i];
        size_t cb = p->col_offsets[i], ce = p->col_offsets[i+1];
        for (size_t q = cb; q < ce; ++q) {
            size_t row = p->col_rows[q];
            if (row > i) a[row*n+i] *= inverse;
        }
        for (size_t q = p->row_offsets[i]; q < p->row_offsets[i+1]; ++q) {
            size_t col = p->row_columns[q];
            if (col <= i) continue;
            for (size_t z = cb; z < ce; ++z) {
                size_t row = p->col_rows[z];
                if (row > i) a[row*n+col] -= a[row*n+i]*a[i*n+col];
            }
        }
    }
    return 1;
}
void chem_sparse_solve(const ChemSparsePlan *p, const double *a, double *b)
{
    size_t n = p->n;
    for (size_t i = 0; i < n; ++i)
        for (size_t q = p->row_offsets[i]; q < p->row_offsets[i+1]; ++q) {
            size_t col = p->row_columns[q];
            if (col >= i) break;
            b[i] -= a[i*n+col]*b[col];
        }
    for (size_t i = n; i-- > 0;) {
        for (size_t q = p->row_offsets[i]; q < p->row_offsets[i+1]; ++q) {
            size_t col = p->row_columns[q];
            if (col > i) b[i] -= a[i*n+col]*b[col];
        }
        b[i] /= a[i*n+i];
    }
}
