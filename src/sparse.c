/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Symbolic sparse plans and the MICM Mozart and Doolittle LU algorithms. */
#include "chem_solver.h"
#include "internal.h"
#include <math.h>
#include <string.h>

typedef struct { size_t dst, left, right; } Update;
struct ChemSparse {
    size_t n, nnz, update_count;
    size_t *map, *row_start, *columns, *diagonal;
    size_t *pivot_start, *pivot_end;
    Update *updates;
};
size_t chem_sparse_order(const ChemSparse *a) { return a->n; }
size_t chem_sparse_elements(const ChemSparse *a) { return a->nnz; }
size_t chem_sparse_index(const ChemSparse *a, size_t r, size_t c)
{ return r < a->n && c < a->n ? a->map[r*a->n+c] : SIZE_MAX; }
void chem_sparse_free(ChemSparse *a)
{
    if (!a) return;
    free(a->map); free(a->row_start); free(a->columns); free(a->diagonal);
    free(a->pivot_start); free(a->pivot_end); free(a->updates); free(a);
}
ChemStatus chem_sparse_create(size_t n, const unsigned char *pattern, ChemSparse **out)
{
    if (!out || (n && !pattern)) return CHEM_INVALID;
    *out = NULL;
    if (n && n > SIZE_MAX/n) return CHEM_MEMORY;
    unsigned char *fill = chem_alloc(n*n, 1);
    ChemSparse *a = calloc(1, sizeof(*a));
    if (!a || !fill) { free(a); free(fill); return CHEM_MEMORY; }
    a->n=n;
    memcpy(fill,pattern,n*n);
    for (size_t i=0;i<n;++i) fill[i*n+i]=1;
    for (size_t i=0;i<n;++i)
        for (size_t k=i+1;k<n;++k) if (fill[i*n+k])
            for (size_t j=i+1;j<n;++j) if (fill[j*n+i]) fill[j*n+k]=1;
    a->map=chem_alloc(n*n,sizeof(size_t));
    a->row_start=chem_alloc(n+1,sizeof(size_t));
    a->diagonal=chem_alloc(n,sizeof(size_t));
    a->pivot_start=chem_alloc(n,sizeof(size_t));
    a->pivot_end=chem_alloc(n,sizeof(size_t));
    if (!a->map || !a->row_start || !a->diagonal || !a->pivot_start || !a->pivot_end) goto fail;
    for (size_t i=0;i<n;++i) {
        a->row_start[i]=a->nnz;
        for (size_t j=0;j<n;++j) a->map[i*n+j]=fill[i*n+j]?a->nnz++:SIZE_MAX;
        a->diagonal[i]=a->map[i*n+i];
    }
    a->row_start[n]=a->nnz;
    a->columns=chem_alloc(a->nnz,sizeof(size_t));
    if (!a->columns) goto fail;
    for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j)
        if (fill[i*n+j]) a->columns[a->map[i*n+j]]=j;
    for (size_t i=0;i<n;++i)
        for (size_t k=i+1;k<n;++k) if (fill[i*n+k])
            for (size_t j=i+1;j<n;++j) if (fill[j*n+i]) ++a->update_count;
    a->updates=chem_alloc(a->update_count,sizeof(Update));
    if (!a->updates) goto fail;
    size_t u=0;
    for (size_t i=0;i<n;++i) {
        a->pivot_start[i]=u;
        for (size_t k=i+1;k<n;++k) if (fill[i*n+k])
            for (size_t j=i+1;j<n;++j) if (fill[j*n+i])
                a->updates[u++]=(Update){a->map[j*n+k],a->map[j*n+i],a->map[i*n+k]};
        a->pivot_end[i]=u;
    }
    free(fill); *out=a; return CHEM_OK;
fail:
    free(fill); chem_sparse_free(a); return CHEM_MEMORY;
}
ChemStatus chem_lu_factor(const ChemSparse *a, size_t cells, double *v, ChemLuAlgorithm algorithm)
{
    size_t n=a->n;
    if (algorithm == CHEM_LU_MOZART) {
        for (size_t i=0;i<n;++i) {
            const double *diag=v+a->diagonal[i]*cells;
            for (size_t c=0;c<cells;++c) if (diag[c]==0 || !isfinite(diag[c])) return CHEM_SINGULAR;
            for (size_t j=i+1;j<n;++j) {
                size_t ij=a->map[j*n+i];
                if (ij==SIZE_MAX) continue;
                double *restrict dst=v+ij*cells;
                #pragma omp simd
                for (size_t c=0;c<cells;++c) dst[c] *= 1.0/diag[c];
            }
            for (size_t u=a->pivot_start[i];u<a->pivot_end[i];++u) {
                Update op=a->updates[u];
                double *restrict dst=v+op.dst*cells;
                const double *left=v+op.left*cells, *right=v+op.right*cells;
                #pragma omp simd
                for (size_t c=0;c<cells;++c) dst[c] -= left[c]*right[c];
            }
        }
    } else if (algorithm == CHEM_LU_DOOLITTLE) {
        for (size_t i=0;i<n;++i) {
            for (size_t k=i;k<n;++k) {
                size_t ik=a->map[i*n+k];
                if (ik==SIZE_MAX) continue;
                double *dst=v+ik*cells;
                for (size_t j=0;j<i;++j) {
                    size_t ij=a->map[i*n+j], jk=a->map[j*n+k];
                    if (ij==SIZE_MAX || jk==SIZE_MAX) continue;
                    #pragma omp simd
                    for (size_t c=0;c<cells;++c) dst[c]-=v[ij*cells+c]*v[jk*cells+c];
                }
            }
            const double *diag=v+a->diagonal[i]*cells;
            for (size_t c=0;c<cells;++c) if (diag[c]==0 || !isfinite(diag[c])) return CHEM_SINGULAR;
            for (size_t k=i+1;k<n;++k) {
                size_t ki=a->map[k*n+i];
                if (ki==SIZE_MAX) continue;
                double *dst=v+ki*cells;
                for (size_t j=0;j<i;++j) {
                    size_t kj=a->map[k*n+j], ji=a->map[j*n+i];
                    if (kj==SIZE_MAX || ji==SIZE_MAX) continue;
                    #pragma omp simd
                    for (size_t c=0;c<cells;++c) dst[c]-=v[kj*cells+c]*v[ji*cells+c];
                }
                #pragma omp simd
                for (size_t c=0;c<cells;++c) dst[c]/=diag[c];
            }
        }
    } else return CHEM_INVALID;
    return CHEM_OK;
}
void chem_lu_solve(const ChemSparse *a, size_t cells, const double *v, double *x)
{
    for (size_t i=0;i<a->n;++i) {
        double *restrict dst=x+i*cells;
        for (size_t k=a->row_start[i];k<a->diagonal[i];++k) {
            const double *lij=v+k*cells, *yj=x+a->columns[k]*cells;
            #pragma omp simd
            for (size_t c=0;c<cells;++c) dst[c]-=lij[c]*yj[c];
        }
    }
    for (size_t i=a->n;i-->0;) {
        double *restrict dst=x+i*cells;
        for (size_t k=a->diagonal[i]+1;k<a->row_start[i+1];++k) {
            const double *uij=v+k*cells, *xj=x+a->columns[k]*cells;
            #pragma omp simd
            for (size_t c=0;c<cells;++c) dst[c]-=uij[c]*xj[c];
        }
        const double *diag=v+a->diagonal[i]*cells;
        #pragma omp simd
        for (size_t c=0;c<cells;++c) dst[c]/=diag[c];
    }
}
ChemStatus chem_markowitz(size_t n, const unsigned char *pattern, size_t *perm)
{
    if (n && (!pattern || !perm)) return CHEM_INVALID;
    if (n && n>SIZE_MAX/n) return CHEM_MEMORY;
    unsigned char *work=chem_alloc(n*n,1), *alive=chem_alloc(n,1);
    if (!work || !alive) { free(work); free(alive); return CHEM_MEMORY; }
    memcpy(work,pattern,n*n); memset(alive,1,n);
    for (size_t step=0;step<n;++step) {
        size_t pivot=n,best=SIZE_MAX;
        for (size_t v=0;v<n;++v) if (alive[v]) {
            size_t row=0,col=0;
            for (size_t j=0;j<n;++j) if (alive[j]) { row+=!!work[v*n+j]; col+=!!work[j*n+v]; }
            /* Unsigned arithmetic matches MICM, including structurally empty rows. */
            size_t cost=(row-1)*(col-1);
            if (pivot==n || cost<best) { best=cost; pivot=v; }
        }
        perm[step]=pivot; alive[pivot]=0;
        for (size_t i=0;i<n;++i) if (alive[i] && work[i*n+pivot])
            for (size_t j=0;j<n;++j) if (alive[j] && work[pivot*n+j]) work[i*n+j]=1;
    }
    free(work); free(alive); return CHEM_OK;
}
