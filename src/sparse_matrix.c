#include "sparse_matrix.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void chem_sparse_matrix_destroy(ChemSparseMatrix *m)
{
    if (m) { free(m->offsets); free(m->indices); chem_matrix_destroy(&m->values); *m=(ChemSparseMatrix){0}; }
}
static int compare_size(const void *a,const void *b)
{
    size_t x=*(const size_t *)a,y=*(const size_t *)b;
    return (x>y)-(x<y);
}
ChemMatrixStatus chem_sparse_matrix_init(ChemSparseMatrix *m,size_t n,size_t blocks,size_t width,
        ChemSparseOrdering ordering,const ChemSparseEntry *entries,size_t count,double initial)
{
    if (!m || (count && !entries) || (ordering!=CHEM_SPARSE_CSR && ordering!=CHEM_SPARSE_CSC)) return CHEM_MATRIX_INVALID;
    if (n==SIZE_MAX || (n && n>SIZE_MAX/n) || n+1>SIZE_MAX/sizeof(size_t) || count>SIZE_MAX/sizeof(size_t)) return CHEM_MATRIX_SIZE;
    size_t *keys=count ? malloc(count*sizeof(*keys)) : NULL;
    if (count && !keys) return CHEM_MATRIX_ALLOCATION;
    for (size_t i=0;i<count;++i) {
        if (entries[i].row>=n || entries[i].column>=n) { free(keys); return CHEM_MATRIX_INVALID; }
        size_t major=ordering==CHEM_SPARSE_CSR ? entries[i].row : entries[i].column;
        size_t minor=ordering==CHEM_SPARSE_CSR ? entries[i].column : entries[i].row;
        keys[i]=major*n+minor;
    }
    if (count>1) qsort(keys,count,sizeof(*keys),compare_size);
    size_t unique=0;
    for (size_t i=0;i<count;++i) if (!i || keys[i]!=keys[i-1]) keys[unique++]=keys[i];
    ChemSparseMatrix result={.dimension=n,.nonzeros=unique,.ordering=ordering};
    result.offsets=calloc(n+1,sizeof(*result.offsets));
    result.indices=unique ? malloc(unique*sizeof(*result.indices)) : NULL;
    if (!result.offsets || (unique && !result.indices)) { free(keys); chem_sparse_matrix_destroy(&result); return CHEM_MATRIX_ALLOCATION; }
    for (size_t i=0;i<unique;++i) { ++result.offsets[keys[i]/n+1]; result.indices[i]=keys[i]%n; }
    for (size_t i=0;i<n;++i) result.offsets[i+1]+=result.offsets[i];
    free(keys);
    ChemMatrixStatus status=chem_matrix_init(&result.values,blocks,unique,width,initial);
    if (status!=CHEM_MATRIX_OK) { chem_sparse_matrix_destroy(&result); return status; }
    chem_sparse_matrix_destroy(m); *m=result; return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_sparse_matrix_clone(ChemSparseMatrix *dst,const ChemSparseMatrix *src)
{
    if (!dst || !src || !src->offsets) return CHEM_MATRIX_INVALID;
    if (dst==src) return CHEM_MATRIX_OK;
    ChemSparseMatrix result={.dimension=src->dimension,.nonzeros=src->nonzeros,.ordering=src->ordering};
    result.offsets=malloc((src->dimension+1)*sizeof(*result.offsets));
    result.indices=src->nonzeros ? malloc(src->nonzeros*sizeof(*result.indices)) : NULL;
    if (!result.offsets || (src->nonzeros && !result.indices)) { chem_sparse_matrix_destroy(&result); return CHEM_MATRIX_ALLOCATION; }
    memcpy(result.offsets,src->offsets,(src->dimension+1)*sizeof(*result.offsets));
    if (src->nonzeros) memcpy(result.indices,src->indices,src->nonzeros*sizeof(*result.indices));
    ChemMatrixStatus status=chem_matrix_clone(&result.values,&src->values);
    if (status!=CHEM_MATRIX_OK) { chem_sparse_matrix_destroy(&result); return status; }
    chem_sparse_matrix_destroy(dst); *dst=result; return CHEM_MATRIX_OK;
}
size_t chem_sparse_matrix_entry(const ChemSparseMatrix *m,size_t row,size_t col)
{
    if (!m || !m->offsets || row>=m->dimension || col>=m->dimension) return SIZE_MAX;
    size_t major=m->ordering==CHEM_SPARSE_CSR ? row : col, minor=m->ordering==CHEM_SPARSE_CSR ? col : row;
    size_t lo=m->offsets[major],hi=m->offsets[major+1];
    while (lo<hi) { size_t mid=lo+(hi-lo)/2; if (m->indices[mid]<minor) lo=mid+1; else hi=mid; }
    return lo<m->offsets[major+1] && m->indices[lo]==minor ? lo : SIZE_MAX;
}
size_t chem_sparse_matrix_index(const ChemSparseMatrix *m,size_t block,size_t row,size_t col)
{
    if (!m || block>=m->values.rows) return SIZE_MAX;
    size_t entry=chem_sparse_matrix_entry(m,row,col);
    return entry==SIZE_MAX ? SIZE_MAX : chem_matrix_index(&m->values,block,entry);
}
double *chem_sparse_matrix_at(ChemSparseMatrix *m,size_t block,size_t row,size_t col)
{
    size_t index=chem_sparse_matrix_index(m,block,row,col);
    return index==SIZE_MAX ? NULL : m->values.data+index;
}
const double *chem_sparse_matrix_const_at(const ChemSparseMatrix *m,size_t block,size_t row,size_t col)
{
    size_t index=chem_sparse_matrix_index(m,block,row,col);
    return index==SIZE_MAX ? NULL : m->values.data+index;
}
ChemMatrixStatus chem_sparse_matrix_diagonal_indices(const ChemSparseMatrix *m,size_t block,
        size_t *indices,size_t capacity,size_t *count)
{
    if (!m || !count || block>=m->values.rows || (capacity && !indices)) return CHEM_MATRIX_INVALID;
    size_t needed=0;
    for (size_t i=0;i<m->dimension;++i) if (chem_sparse_matrix_entry(m,i,i)!=SIZE_MAX) ++needed;
    *count=needed;
    if (capacity<needed) return CHEM_MATRIX_SIZE;
    size_t pos=0;
    for (size_t i=0;i<m->dimension;++i) {
        size_t index=chem_sparse_matrix_index(m,block,i,i);
        if (index!=SIZE_MAX) indices[pos++]=index;
    }
    return CHEM_MATRIX_OK;
}
void chem_sparse_matrix_fill(ChemSparseMatrix *m,double value) { chem_matrix_fill(&m->values,value); }
void chem_sparse_matrix_add_diagonal(ChemSparseMatrix *m,double value)
{
    for (size_t i=0;i<m->dimension;++i) {
        size_t entry=chem_sparse_matrix_entry(m,i,i);
        if (entry!=SIZE_MAX) for (size_t block=0;block<m->values.rows;++block)
            *chem_matrix_at(&m->values,block,entry)+=value;
    }
}
ChemMatrixStatus chem_sparse_matrix_print(const ChemSparseMatrix *m,FILE *stream,int nonzeros_only)
{
    if (!m || !stream || !m->offsets) return CHEM_MATRIX_INVALID;
    for (size_t block=0;block<m->values.rows;++block) {
        if (fprintf(stream,"Block %zu\n",block)<0) return CHEM_MATRIX_INVALID;
        for (size_t row=0;row<m->dimension;++row) for (size_t col=0;col<m->dimension;++col) {
            const double *p=chem_sparse_matrix_const_at(m,block,row,col);
            if (nonzeros_only) {
                if (p && fprintf(stream,"%zu, %zu, %g\n",row,col,*p)<0) return CHEM_MATRIX_INVALID;
            } else if (fprintf(stream,col ? ",%g" : "%g",p ? *p : 0.0)<0) return CHEM_MATRIX_INVALID;
            if (!nonzeros_only && col+1==m->dimension && fputc('\n',stream)==EOF) return CHEM_MATRIX_INVALID;
        }
    }
    return CHEM_MATRIX_OK;
}
