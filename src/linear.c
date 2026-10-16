/* Sparse Mozart and Doolittle algorithms: MICM/UCAR, Apache-2.0. */
#include "chem.h"
#include <stdlib.h>
#include <string.h>

void chem_lu_destroy(chem_lu_plan *p)
{
    free(p->map); free(p->row_start); free(p->column); free(p->diagonal); free(p->ops);
    memset(p, 0, sizeof(*p));
}

static int append(chem_lu_plan *p, size_t *capacity, unsigned kind, size_t target, size_t left, size_t right)
{
    if (p->n_ops == *capacity) {
        size_t next = *capacity ? *capacity * 2 : 256;
        if (next < *capacity || next > SIZE_MAX / sizeof(*p->ops)) return 0;
        chem_lu_op *ops = realloc(p->ops, next * sizeof(*ops));
        if (!ops) return 0;
        p->ops = ops; *capacity = next;
    }
    p->ops[p->n_ops++] = (chem_lu_op){target, left, right, kind};
    return 1;
}

chem_status chem_lu_create(chem_lu_plan *p, size_t n, const unsigned char *pattern, chem_lu_algorithm algorithm)
{
    memset(p, 0, sizeof(*p));
    if (!n || !pattern || n > SIZE_MAX / n / sizeof(size_t) || algorithm > CHEM_LU_DOOLITTLE) return CHEM_INVALID;
    p->n=n; p->algorithm=algorithm;
    unsigned char *fill = malloc(n*n);
    if (!fill) return CHEM_MEMORY;
    memcpy(fill, pattern, n*n);
    for (size_t i=0; i<n; ++i) {
        if (!fill[i*n+i]) { free(fill); return CHEM_SINGULAR; }
        for (size_t k=i+1; k<n; ++k) if (fill[i*n+k])
            for (size_t j=i+1; j<n; ++j) if (fill[j*n+i]) fill[j*n+k]=1;
    }
    p->map=malloc(n*n*sizeof(size_t));
    p->row_start=calloc(n+1,sizeof(size_t));
    p->diagonal=malloc(n*sizeof(size_t));
    if (!p->map || !p->row_start || !p->diagonal) goto memory_error;
    for (size_t i=0; i<n; ++i) {
        p->row_start[i]=p->nnz;
        for (size_t j=0; j<n; ++j) p->map[i*n+j] = fill[i*n+j] ? p->nnz++ : CHEM_NONE;
        p->diagonal[i]=p->map[i*n+i];
    }
    p->row_start[n]=p->nnz;
    p->column=malloc(p->nnz*sizeof(size_t));
    if (!p->column) goto memory_error;
    for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j)
        if (fill[i*n+j]) p->column[p->map[i*n+j]]=j;
    size_t capacity=0;
#define EMIT(kind,t,l,r) do { if (!append(p,&capacity,kind,t,l,r)) goto memory_error; } while (0)
#define AT(i,j) p->map[(i)*n+(j)]
    for (size_t i=0; i<n; ++i) {
        if (algorithm==CHEM_LU_MOZART) {
            EMIT(0,AT(i,i),0,0); /* scratch = reciprocal of pivot */
            for (size_t j=i+1; j<n; ++j) if (fill[j*n+i]) EMIT(1,AT(j,i),0,0);
            for (size_t k=i+1; k<n; ++k) if (fill[i*n+k])
                for (size_t j=i+1; j<n; ++j) if (fill[j*n+i]) EMIT(2,AT(j,k),AT(j,i),AT(i,k));
        } else {
            for (size_t k=i; k<n; ++k) if (fill[i*n+k])
                for (size_t j=0; j<i; ++j) if (fill[i*n+j] && fill[j*n+k]) EMIT(2,AT(i,k),AT(i,j),AT(j,k));
            for (size_t k=i+1; k<n; ++k) if (fill[k*n+i]) {
                for (size_t j=0; j<i; ++j) if (fill[k*n+j] && fill[j*n+i]) EMIT(2,AT(k,i),AT(k,j),AT(j,i));
                EMIT(3,AT(k,i),AT(i,i),0);
            }
        }
    }
#undef AT
#undef EMIT
    free(fill);
    return CHEM_OK;
memory_error:
    free(fill); chem_lu_destroy(p); return CHEM_MEMORY;
}

void chem_lu_factor(const chem_lu_plan *p, size_t cells, double *matrix, double *scratch)
{
    for (size_t op=0; op<p->n_ops; ++op) {
        const chem_lu_op *o=&p->ops[op];
        double *a=matrix+o->target*cells;
        const double *b=matrix+o->left*cells, *c=matrix+o->right*cells;
        switch (o->kind) {
        case 0:
            #pragma omp simd
            for (size_t v=0; v<cells; ++v) scratch[v]=1.0/a[v];
            break;
        case 1:
            #pragma omp simd
            for (size_t v=0; v<cells; ++v) a[v]*=scratch[v];
            break;
        case 2:
            #pragma omp simd
            for (size_t v=0; v<cells; ++v) a[v]-=b[v]*c[v];
            break;
        case 3:
            #pragma omp simd
            for (size_t v=0; v<cells; ++v) a[v]/=b[v];
            break;
        }
    }
}

void chem_lu_solve(const chem_lu_plan *p, size_t cells, const double *restrict a, double *restrict x)
{
    for (size_t i=0; i<p->n; ++i) {
        double *xi=x+i*cells;
        for (size_t e=p->row_start[i]; e<p->diagonal[i]; ++e) {
            const double *aij=a+e*cells, *xj=x+p->column[e]*cells;
            #pragma omp simd
            for (size_t c=0; c<cells; ++c) xi[c]-=aij[c]*xj[c];
        }
    }
    for (size_t i=p->n; i-- > 0;) {
        double *xi=x+i*cells;
        for (size_t e=p->diagonal[i]+1; e<p->row_start[i+1]; ++e) {
            const double *aij=a+e*cells, *xj=x+p->column[e]*cells;
            #pragma omp simd
            for (size_t c=0; c<cells; ++c) xi[c]-=aij[c]*xj[c];
        }
        const double *diag=a+p->diagonal[i]*cells;
        #pragma omp simd
        for (size_t c=0; c<cells; ++c) xi[c]/=diag[c];
    }
}

chem_status chem_markowitz(size_t n, const unsigned char *pattern, size_t *permutation)
{
    if (!n || !pattern || !permutation || n>SIZE_MAX/n) return CHEM_INVALID;
    unsigned char *graph=malloc(n*n), *alive=malloc(n);
    size_t *row=calloc(n,sizeof(size_t)), *col=calloc(n,sizeof(size_t));
    if (!graph || !alive || !row || !col) { free(graph); free(alive); free(row); free(col); return CHEM_MEMORY; }
    memcpy(graph,pattern,n*n); memset(alive,1,n);
    for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j)
        if (graph[i*n+j]) { ++row[i]; ++col[j]; }
    for (size_t step=0; step<n; ++step) {
        size_t pivot=CHEM_NONE, cost=SIZE_MAX;
        for (size_t i=0; i<n; ++i) if (alive[i]) {
            /* Unsigned arithmetic also matches the reference's isolated-node handling. */
            size_t candidate=(row[i]-1)*(col[i]-1);
            if (pivot==CHEM_NONE || candidate<cost) { pivot=i; cost=candidate; }
        }
        permutation[step]=pivot; alive[pivot]=0;
        for (size_t i=0; i<n; ++i) if (alive[i] && graph[i*n+pivot])
            for (size_t j=0; j<n; ++j) if (alive[j] && graph[pivot*n+j] && !graph[i*n+j]) {
                graph[i*n+j]=1; ++row[i]; ++col[j];
            }
        for (size_t i=0; i<n; ++i) if (alive[i]) {
            if (graph[pivot*n+i]) --col[i];
            if (graph[i*n+pivot]) --row[i];
        }
    }
    free(graph); free(alive); free(row); free(col); return CHEM_OK;
}
