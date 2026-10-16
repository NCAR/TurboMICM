#include "chem_sparse.h"
#include <stdlib.h>
#include <string.h>
chem_status chem_sparse_builder_create(chem_sparse_builder *b,size_t n)
{
    *b=(chem_sparse_builder){0};
    if (n && n>SIZE_MAX/n) return CHEM_SIZE;
    b->pattern=n ? calloc(n*n,1) : NULL;
    if (n && !b->pattern) return CHEM_MEMORY;
    b->n=n; b->blocks=1; return CHEM_OK;
}
chem_status chem_sparse_builder_add(chem_sparse_builder *b,size_t r,size_t c)
{
    if (r>=b->n || c>=b->n) return CHEM_SIZE;
    if (!b->pattern[r*b->n+c]) { b->pattern[r*b->n+c]=1; ++b->nnz; }
    return CHEM_OK;
}
void chem_sparse_builder_destroy(chem_sparse_builder *b)
{ free(b->pattern); *b=(chem_sparse_builder){0}; }
chem_status chem_sparse_create(chem_sparse *s,const chem_sparse_builder *b,chem_sparse_order order,size_t v,chem_element_type type)
{
    *s=(chem_sparse){0};
    if ((order!=CHEM_CSR && order!=CHEM_CSC) || (b->n && b->n>SIZE_MAX/b->n/sizeof(size_t))) return CHEM_INVALID;
    s->n=b->n; s->nnz=b->nnz; s->order=order;
    s->map=b->n ? malloc(b->n*b->n*sizeof(size_t)) : NULL;
    s->starts=calloc(b->n+1,sizeof(size_t)); s->indices=b->nnz ? malloc(b->nnz*sizeof(size_t)) : NULL;
    if ((b->n && !s->map) || !s->starts || (b->nnz && !s->indices)) { chem_sparse_destroy(s);return CHEM_MEMORY; }
    for(size_t i=0;i<b->n*b->n;++i) s->map[i]=CHEM_NONE;
    size_t k=0;
    for(size_t major=0;major<b->n;++major) {
        s->starts[major]=k;
        for(size_t minor=0;minor<b->n;++minor) {
            size_t r=order==CHEM_CSR ? major : minor,c=order==CHEM_CSR ? minor : major;
            if(b->pattern[r*b->n+c]) { s->map[r*b->n+c]=k; s->indices[k++]=minor; }
        }
    }
    s->starts[b->n]=k;
    chem_status status=chem_array_create(&s->values,b->blocks,b->nnz,v,type,b->initial);
    if(status) { chem_sparse_destroy(s);return status; }
    s->values.kind=CHEM_SPARSE_ARRAY;return CHEM_OK;
}
void chem_sparse_destroy(chem_sparse *s)
{ free(s->map);free(s->starts);free(s->indices);chem_array_destroy(&s->values);*s=(chem_sparse){0}; }
chem_status chem_sparse_clone(chem_sparse *d,const chem_sparse *s)
{
    chem_sparse_builder b; chem_status status=chem_sparse_builder_create(&b,s->n);
    if(status) return status;
    b.blocks=s->values.rows;
    for(size_t i=0;i<s->n;++i) for(size_t j=0;j<s->n;++j) if(s->map[i*s->n+j]!=CHEM_NONE) chem_sparse_builder_add(&b,i,j);
    status=chem_sparse_create(d,&b,s->order,s->values.vector_length,s->values.type);
    chem_sparse_builder_destroy(&b);
    if(!status && s->values.size) memcpy(d->values.data,s->values.data,s->values.size*chem_element_size(s->values.type));
    return status;
}
chem_status chem_sparse_is_zero(const chem_sparse *s,size_t r,size_t c,int *zero)
{
    if(r>=s->n || c>=s->n) return CHEM_SIZE;
    *zero=s->map[r*s->n+c]==CHEM_NONE;return CHEM_OK;
}
chem_status chem_sparse_index(const chem_sparse *s,size_t b,size_t r,size_t c,size_t *index)
{
    *index=CHEM_NONE;
    if(b==CHEM_NONE) { if(s->values.rows!=1) return CHEM_INVALID; b=0; }
    if(b>=s->values.rows || r>=s->n || c>=s->n) return CHEM_SIZE;
    size_t k=s->map[r*s->n+c];if(k==CHEM_NONE) return CHEM_NOT_FOUND;
    *index=chem_array_index(&s->values,b,k);return CHEM_OK;
}
chem_status chem_sparse_get(const chem_sparse *s,size_t b,size_t r,size_t c,double *value)
{
    size_t index; chem_status status=chem_sparse_index(s,b,r,c,&index);if(status) return status;
    (void)index;*value=chem_array_get(&s->values,b==CHEM_NONE ? 0 : b,s->map[r*s->n+c]);return CHEM_OK;
}
chem_status chem_sparse_set(chem_sparse *s,size_t b,size_t r,size_t c,double value)
{
    size_t index; chem_status status=chem_sparse_index(s,b,r,c,&index);if(status) return status;
    (void)index;return chem_array_set(&s->values,b==CHEM_NONE ? 0 : b,s->map[r*s->n+c],value);
}
chem_status chem_sparse_diagonal_indices(const chem_sparse *s,size_t b,size_t *ids,size_t capacity,size_t *count)
{
    *count=0;if(b>=s->values.rows) return CHEM_SIZE;
    for(size_t i=0;i<s->n;++i) if(s->map[i*s->n+i]!=CHEM_NONE) ++*count;
    if(capacity<*count) return CHEM_SIZE;
    size_t k=0;for(size_t i=0;i<s->n;++i) if(s->map[i*s->n+i]!=CHEM_NONE) chem_sparse_index(s,b,i,i,&ids[k++]);
    return CHEM_OK;
}
void chem_sparse_add_diagonal(chem_sparse *s,double value)
{
    for(size_t i=0;i<s->n;++i) if(s->map[i*s->n+i]!=CHEM_NONE)
        for(size_t b=0;b<s->values.rows;++b) {
            size_t k=s->map[i*s->n+i]; chem_array_set(&s->values,b,k,chem_array_get(&s->values,b,k)+value);
        }
}
chem_status chem_sparse_matvec(const chem_sparse *s,const chem_array *x,chem_array *y)
{
    if (!s || !x || !y || x->data==y->data || s->values.type!=CHEM_DOUBLE || x->type!=CHEM_DOUBLE || y->type!=CHEM_DOUBLE)
        return CHEM_INVALID;
    if (x->rows!=s->values.rows || y->rows!=x->rows || x->columns!=s->n || y->columns!=s->n) return CHEM_SIZE;
    for (size_t b=0;b<x->rows;++b) for (size_t r=0;r<s->n;++r) {
        double sum=0;
        for (size_t c=0;c<s->n;++c) {
            size_t k=s->map[r*s->n+c];
            if (k!=CHEM_NONE) sum+=chem_array_get(&s->values,b,k)*chem_array_get(x,b,c);
        }
        chem_array_set(y,b,r,sum);
    }
    return CHEM_OK;
}
chem_status chem_sparse_print(FILE *f,const chem_sparse *s,int only)
{
    for(size_t b=0;b<s->values.rows;++b) {
        if(fprintf(f,"Block %zu\n",b)<0) return CHEM_INVALID;
        for(size_t r=0;r<s->n;++r) {
            for(size_t c=0;c<s->n;++c) {
                size_t k=s->map[r*s->n+c];double value=k==CHEM_NONE ? 0 : chem_array_get(&s->values,b,k);
                if(only) { if(k!=CHEM_NONE && fprintf(f,"%zu, %zu, %g\n",r,c,value)<0) return CHEM_INVALID; }
                else if(fprintf(f,"%s%g",c ? "," : "",value)<0) return CHEM_INVALID;
            }
            if(!only && fputc('\n',f)==EOF) return CHEM_INVALID;
        }
    }
    return CHEM_OK;
}
chem_span chem_sparse_block_span(chem_array_group g,size_t index)
{
    size_t width=g.array->vector_length ? g.array->vector_length : 1;
    if(g.array->kind!=CHEM_SPARSE_ARRAY || index%width) return (chem_span){.status=CHEM_INVALID};
    chem_span span=chem_group_column(g,index/width);span.kind=CHEM_ELEMENT_VIEW;return span;
}
chem_span chem_sparse_const_block_span(chem_array_group g,size_t index)
{ g.writable=0;return chem_sparse_block_span(g,index); }
chem_span chem_sparse_coordinate_span(const chem_sparse *s,chem_array_group g,size_t r,size_t c)
{
    if(g.array!=&s->values) return (chem_span){.status=CHEM_INVALID};
    size_t index;chem_status status=chem_sparse_index(s,0,r,c,&index);
    return status ? (chem_span){.status=status} : chem_sparse_block_span(g,index);
}
chem_raw_view chem_sparse_block(const chem_sparse *s,size_t index,int writable)
{ return (chem_raw_view){&s->values,index,writable,CHEM_RAW_ELEMENT_VIEW}; }
