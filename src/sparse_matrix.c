/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "chem_sparse_matrix.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
struct ChemSparseBuilder {
    size_t order,blocks,count,capacity;
    ChemSparseElement *elements;
    double initial;
};
ChemStatus chem_sparse_builder_create(size_t order,ChemSparseBuilder **out)
{
    if(!out)return CHEM_INVALID;
    ChemSparseBuilder *b=calloc(1,sizeof(*b));if(!b)return CHEM_MEMORY;
    b->order=order;b->blocks=1;*out=b;return CHEM_OK;
}
void chem_sparse_builder_free(ChemSparseBuilder *b){if(b){free(b->elements);free(b);}}
ChemStatus chem_sparse_builder_add(ChemSparseBuilder *b,size_t r,size_t c)
{
    if(!b)return CHEM_INVALID;
    if(r>=b->order||c>=b->order)return CHEM_DIMENSION;
    for(size_t i=0;i<b->count;++i)if(b->elements[i].row==r&&b->elements[i].column==c)return CHEM_OK;
    if(b->count==SIZE_MAX||(b->blocks && b->count+1>SIZE_MAX/b->blocks))return CHEM_MEMORY;
    if(b->count==b->capacity){
        size_t cap=b->capacity?b->capacity*2:8;
        if(cap<b->capacity||cap>SIZE_MAX/sizeof(ChemSparseElement))return CHEM_MEMORY;
        ChemSparseElement *p=realloc(b->elements,cap*sizeof(*p));if(!p)return CHEM_MEMORY;
        b->elements=p;b->capacity=cap;
    }
    b->elements[b->count++]=(ChemSparseElement){r,c};return CHEM_OK;
}
ChemStatus chem_sparse_builder_blocks(ChemSparseBuilder *b,size_t n)
{
    if(!b)return CHEM_INVALID;
    if(b->count&&n>SIZE_MAX/b->count)return CHEM_MEMORY;
    b->blocks=n;return CHEM_OK;
}
void chem_sparse_builder_initial(ChemSparseBuilder *b,double v){if(b)b->initial=v;}
size_t chem_sparse_builder_elements(const ChemSparseBuilder *b){return b?b->count*b->blocks:0;}
static int compare_row(const void *lhs,const void *rhs)
{
    const ChemSparseElement *a=lhs,*b=rhs;
    if(a->row!=b->row)return a->row<b->row?-1:1;
    return a->column<b->column?-1:a->column>b->column;
}
static int compare_column(const void *lhs,const void *rhs)
{
    const ChemSparseElement *a=lhs,*b=rhs;
    if(a->column!=b->column)return a->column<b->column?-1:1;
    return a->row<b->row?-1:a->row>b->row;
}
ChemStatus chem_sparse_matrix_build(const ChemSparseBuilder *b,size_t group,ChemSparseOrdering ordering,ChemMatrixType type,ChemSparseMatrix *out)
{
    if(!b||!out||(ordering!=CHEM_SPARSE_ROW&&ordering!=CHEM_SPARSE_COLUMN))return CHEM_INVALID;
    ChemSparseMatrix m={.order=b->order,.ordering=ordering};
    if(b->count){m.elements=malloc(b->count*sizeof(*m.elements));if(!m.elements)return CHEM_MEMORY;
        memcpy(m.elements,b->elements,b->count*sizeof(*m.elements));
        qsort(m.elements,b->count,sizeof(*m.elements),ordering==CHEM_SPARSE_ROW?compare_row:compare_column);}
    ChemStatus s=chem_matrix_create(b->blocks,b->count,group,type,b->initial,&m.values);
    if(s!=CHEM_OK){free(m.elements);return s;}
    *out=m;return CHEM_OK;
}
void chem_sparse_matrix_free(ChemSparseMatrix *m)
{if(m){free(m->elements);chem_matrix_free(&m->values);*m=(ChemSparseMatrix){0};}}
ChemStatus chem_sparse_matrix_copy(const ChemSparseMatrix *m,ChemSparseMatrix *out)
{
    if(!m||!out||m==out)return CHEM_INVALID;
    ChemSparseMatrix copy={.order=m->order,.ordering=m->ordering};
    if(m->values.columns){copy.elements=malloc(m->values.columns*sizeof(*copy.elements));if(!copy.elements)return CHEM_MEMORY;
        memcpy(copy.elements,m->elements,m->values.columns*sizeof(*copy.elements));}
    ChemStatus s=chem_matrix_copy(&m->values,&copy.values);
    if(s!=CHEM_OK){free(copy.elements);return s;}
    *out=copy;return CHEM_OK;
}
static ChemStatus entry(const ChemSparseMatrix *m,size_t row,size_t column,size_t *out)
{
    if(!m||!out)return CHEM_INVALID;
    if(row>=m->order||column>=m->order)return CHEM_DIMENSION;
    ChemSparseElement key={row,column};size_t lo=0,hi=m->values.columns;
    int (*compare)(const void *,const void *)=m->ordering==CHEM_SPARSE_ROW?compare_row:compare_column;
    while(lo<hi){size_t mid=lo+(hi-lo)/2;int c=compare(&key,m->elements+mid);
        if(c==0){*out=mid;return CHEM_OK;}if(c<0)hi=mid;else lo=mid+1;}
    return CHEM_NOT_FOUND;
}
ChemStatus chem_sparse_matrix_index(const ChemSparseMatrix *m,size_t b,size_t r,size_t c,size_t *out)
{
    if(!m||!out)return CHEM_INVALID;
    if(b>=m->values.rows)return CHEM_DIMENSION;
    size_t i;ChemStatus s=entry(m,r,c,&i);if(s!=CHEM_OK)return s;
    *out=chem_matrix_index(&m->values,b,i);return CHEM_OK;
}
ChemStatus chem_sparse_matrix_single_index(const ChemSparseMatrix *m,size_t r,size_t c,size_t *out)
{
    if(!m)return CHEM_INVALID;
    if(m->values.rows!=1)return CHEM_MISSING_BLOCK;
    return chem_sparse_matrix_index(m,0,r,c,out);
}
ChemStatus chem_sparse_matrix_is_zero(const ChemSparseMatrix *m,size_t r,size_t c,int *out)
{
    if(!out)return CHEM_INVALID;
    size_t i;ChemStatus s=entry(m,r,c,&i);
    if(s!=CHEM_OK&&s!=CHEM_NOT_FOUND)return s;
    *out=s==CHEM_NOT_FOUND;return CHEM_OK;
}
double chem_sparse_matrix_get(const ChemSparseMatrix *m,size_t b,size_t r,size_t c)
{
    size_t i;if(!m||b>=m->values.rows||entry(m,r,c,&i)!=CHEM_OK)return NAN;
    return chem_matrix_get(&m->values,b,i);
}
ChemStatus chem_sparse_matrix_set(ChemSparseMatrix *m,size_t b,size_t r,size_t c,double value)
{
    if(!m)return CHEM_INVALID;
    if(b>=m->values.rows)return CHEM_DIMENSION;
    size_t i;ChemStatus s=entry(m,r,c,&i);return s==CHEM_OK?chem_matrix_set(&m->values,b,i,value):s;
}
ChemStatus chem_sparse_matrix_diagonal_indices(const ChemSparseMatrix *m,size_t b,size_t *out,size_t cap,size_t *count)
{
    if(!m||!count)return CHEM_INVALID;
    if(b>=m->values.rows)return CHEM_DIMENSION;
    size_t n=0;for(size_t i=0;i<m->values.columns;++i)n+=m->elements[i].row==m->elements[i].column;
    *count=n;if(cap<n)return CHEM_DIMENSION;if(n&&!out)return CHEM_INVALID;
    size_t j=0;for(size_t i=0;i<m->values.columns;++i)if(m->elements[i].row==m->elements[i].column)out[j++]=chem_matrix_index(&m->values,b,i);
    return CHEM_OK;
}
ChemStatus chem_sparse_matrix_add_diagonal(ChemSparseMatrix *m,double v)
{
    if(!m)return CHEM_INVALID;
    for(size_t i=0;i<m->values.columns;++i)if(m->elements[i].row==m->elements[i].column)
        for(size_t b=0;b<m->values.rows;++b){ChemStatus s=chem_matrix_set(&m->values,b,i,chem_matrix_get(&m->values,b,i)+v);if(s!=CHEM_OK)return s;}
    return CHEM_OK;
}
ChemStatus chem_sparse_matrix_print(const ChemSparseMatrix *m,FILE *out,int nonzero)
{
    if(!m||!out)return CHEM_INVALID;
    for(size_t b=0;b<m->values.rows;++b){
        if(fprintf(out,"Block %zu\n",b)<0)return CHEM_INVALID;
        for(size_t r=0;r<m->order;++r){
            for(size_t c=0;c<m->order;++c){
                size_t i;int found=entry(m,r,c,&i)==CHEM_OK;double v=found?chem_matrix_get(&m->values,b,i):0;
                if(nonzero){if(found&&fprintf(out,"%zu, %zu, %g\n",r,c,v)<0)return CHEM_INVALID;}
                else if(fprintf(out,c?",%g":"%g",v)<0)return CHEM_INVALID;
            }
            if(!nonzero&&fputc('\n',out)==EOF)return CHEM_INVALID;
        }
    }
    return CHEM_OK;
}
ChemArray chem_array_sparse_matrix(ChemSparseMatrix *m)
{ChemArray a=chem_array_matrix(&m->values);a.kind=CHEM_ARRAY_SPARSE;return a;}
ChemArray chem_array_const_sparse_matrix(const ChemSparseMatrix *m)
{ChemArray a=chem_array_const_matrix(&m->values);a.kind=CHEM_ARRAY_SPARSE;return a;}
ChemStatus chem_sparse_block_view(ChemArray a,size_t index,ChemView *out)
{
    if(a.kind!=CHEM_ARRAY_SPARSE)return CHEM_INVALID;
    size_t group=a.group?a.group:1;
    if(index%group)return CHEM_DIMENSION;
    return chem_array_column(a,index/group,out);
}
