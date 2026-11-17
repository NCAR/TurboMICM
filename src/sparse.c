#include "sparse.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

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

void chem_sparse_clear(const ChemSparsePlan *p, double *a)
{
    for (size_t row=0;row<p->n;++row)
        for (size_t q=p->row_offsets[row];q<p->row_offsets[row+1];++q)
            a[row*p->n+p->row_columns[q]]=0;
}
void chem_sparse_negate(const ChemSparsePlan *p, double *a)
{
    for (size_t row=0;row<p->n;++row)
        for (size_t q=p->row_offsets[row];q<p->row_offsets[row+1];++q) {
            size_t index=row*p->n+p->row_columns[q];
            a[index]=-a[index];
        }
}

ChemSparsePlan *chem_sparse_plan_from_pattern(size_t n, const unsigned char *input)
{
    if (!n || !input || n > SIZE_MAX/n || n >= SIZE_MAX/sizeof(size_t)) return NULL;
    unsigned char *pattern = malloc(n*n);
    ChemSparsePlan *p = calloc(1, sizeof(*p));
    if (!pattern || !p) { free(pattern); free(p); return NULL; }
    memcpy(pattern,input,n*n);
    p->n = n;
    for (size_t i=0; i<n; ++i) pattern[i*n+i] = 1;
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

ChemSparsePlan *chem_sparse_plan_create(const ChemMechanism *m)
{
    if (!m || !m->species_count || m->species_count > SIZE_MAX/m->species_count) return NULL;
    size_t n=m->species_count;
    unsigned char *pattern=calloc(n*n,1);
    if (!pattern) return NULL;
    if (!chem_jacobian_pattern(m,pattern)) { free(pattern); return NULL; }
    ChemSparsePlan *p=chem_sparse_plan_from_pattern(n,pattern);
    free(pattern); return p;
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
int chem_sparse_factor_copy(const ChemSparsePlan *p, const double *a, double *factors)
{
    if (!p || !a || !factors || p->n > SIZE_MAX/p->n || p->n*p->n > SIZE_MAX/sizeof(double)) return 0;
    memcpy(factors,a,p->n*p->n*sizeof(double));
    return chem_sparse_factor(p,factors);
}
/* Intersect a factored row and column up to the current pivot. */
static double subtract_products(const ChemSparsePlan *p, const double *a,
                                size_t row, size_t col, size_t pivot, double value)
{
    size_t r=p->row_offsets[row], re=p->row_offsets[row+1];
    size_t c=p->col_offsets[col], ce=p->col_offsets[col+1], n=p->n;
    while (r<re && c<ce) {
        size_t rj=p->row_columns[r], cj=p->col_rows[c];
        if (rj >= pivot || cj >= pivot) break;
        if (rj < cj) ++r;
        else if (cj < rj) ++c;
        else { value -= a[row*n+rj]*a[rj*n+col]; ++r; ++c; }
    }
    return value;
}

int chem_sparse_factor_doolittle(const ChemSparsePlan *p, double *a)
{
    size_t n=p->n;
    for (size_t i=0; i<n; ++i) {
        for (size_t q=p->row_offsets[i]; q<p->row_offsets[i+1]; ++q) {
            size_t col=p->row_columns[q];
            if (col >= i) a[i*n+col]=subtract_products(p,a,i,col,i,a[i*n+col]);
        }
        if (a[i*n+i] == 0 || !isfinite(a[i*n+i])) return 0;
        for (size_t q=p->col_offsets[i]; q<p->col_offsets[i+1]; ++q) {
            size_t row=p->col_rows[q];
            if (row > i) a[row*n+i]=subtract_products(p,a,row,i,i,a[row*n+i])/a[i*n+i];
        }
    }
    return 1;
}

int chem_sparse_factor_doolittle_copy(const ChemSparsePlan *p, const double *a, double *factors)
{
    if (!p || !a || !factors || p->n > SIZE_MAX/p->n || p->n*p->n > SIZE_MAX/sizeof(double)) return 0;
    memcpy(factors,a,p->n*p->n*sizeof(double));
    return chem_sparse_factor_doolittle(p,factors);
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
