#include "sparse.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int chem_sparse_builder_init(ChemSparseBuilder *b,size_t n)
{
    if (!b || (n && n>SIZE_MAX/n)) return 1;
    unsigned char *pattern=n ? calloc(n*n,1) : NULL;
    if (n && !pattern) return 1;
    *b=(ChemSparseBuilder){n,0,pattern}; return 0;
}
void chem_sparse_builder_free(ChemSparseBuilder *b)
{
    if (!b) return;
    free(b->pattern); *b=(ChemSparseBuilder){0};
}
int chem_sparse_builder_add(ChemSparseBuilder *b,size_t row,size_t column)
{
    if (!b || row>=b->order || column>=b->order) return 1;
    size_t index=row*b->order+column;
    if (!b->pattern[index]) { b->pattern[index]=1; ++b->count; }
    return 0;
}
int chem_sparse_init(ChemSparseMatrix *m,const ChemSparseBuilder *b,size_t blocks,size_t width,int column_major,double initial)
{
    if (!m || !b || !width || blocks>SIZE_MAX-(width-1)) return 1;
    size_t n=b->order, groups=(blocks+width-1)/width;
    if ((n && (n>SIZE_MAX/n || n*n>SIZE_MAX/sizeof(size_t))) ||
        (groups && (b->count>SIZE_MAX/groups || groups*b->count>SIZE_MAX/width || groups*b->count*width>SIZE_MAX/sizeof(double)))) return 1;
    size_t *indices=n ? malloc(n*n*sizeof(*indices)) : NULL;
    size_t count=groups*b->count*width;
    double *data=count ? malloc(count*sizeof(*data)) : NULL;
    if ((n && !indices) || (count && !data)) { free(indices); free(data); return 1; }
    for(size_t i=0;i<n*n;++i) indices[i]=SIZE_MAX;
    size_t index=0;
    for(size_t outer=0;outer<n;++outer) for(size_t inner=0;inner<n;++inner) {
        size_t row=column_major ? inner : outer, col=column_major ? outer : inner;
        if (b->pattern[row*n+col]) indices[row*n+col]=index++;
    }
    for(size_t i=0;i<count;++i) data[i]=initial;
    *m=(ChemSparseMatrix){n,blocks,width,groups,b->count,count,column_major,indices,data}; return 0;
}
void chem_sparse_free(ChemSparseMatrix *m)
{
    if (!m) return;
    free(m->indices); free(m->data); *m=(ChemSparseMatrix){0};
}
int chem_sparse_index(const ChemSparseMatrix *m,size_t block,size_t row,size_t column,size_t *index)
{
    if (!m || !index || block>=m->blocks || row>=m->order || column>=m->order) return 1;
    size_t flat=m->indices[row*m->order+column];
    if (flat==SIZE_MAX) return 1;
    *index=(block/m->width)*m->nonzeros*m->width+flat*m->width+block%m->width;
    return 0;
}
int chem_sparse_index_single(const ChemSparseMatrix *m,size_t row,size_t column,size_t *index)
{
    return !m || m->blocks!=1 ? 1 : chem_sparse_index(m,0,row,column,index);
}
int chem_sparse_is_zero(const ChemSparseMatrix *m,size_t row,size_t column,int *zero)
{
    if (!m || !zero || row>=m->order || column>=m->order) return 1;
    *zero=m->indices[row*m->order+column]==SIZE_MAX; return 0;
}
int chem_sparse_diagonal_indices(const ChemSparseMatrix *m,size_t block,size_t *indices,size_t capacity,size_t *count)
{
    if (!m || block>=m->blocks || !count) return 1;
    size_t n=0;
    for (size_t i=0;i<m->order;++i) {
        size_t index;
        if (!chem_sparse_index(m,block,i,i,&index)) {
            if (n>=capacity || !indices) return 1;
            indices[n++]=index;
        }
    }
    *count=n; return 0;
}

double *chem_sparse_at(ChemSparseMatrix *m,size_t block,size_t row,size_t column)
{
    size_t i; return chem_sparse_index(m,block,row,column,&i) ? NULL : &m->data[i];
}
const double *chem_sparse_const_at(const ChemSparseMatrix *m,size_t block,size_t row,size_t column)
{
    size_t i; return chem_sparse_index(m,block,row,column,&i) ? NULL : &m->data[i];
}
int chem_sparse_group_block(ChemSparseMatrix *m,size_t group,size_t row,size_t column,ChemSparseGroup *view)
{
    if (!m || !view || group>=m->groups) return 1;
    size_t first=group*m->width, count=m->blocks-first;
    if(count>m->width) count=m->width;
    double *base=chem_sparse_at(m,first,row,column);
    if(!base) return 1;
    *view=(ChemSparseGroup){base,count}; return 0;
}
int chem_sparse_const_group_block(const ChemSparseMatrix *m,size_t group,size_t row,size_t column,ChemConstSparseGroup *view)
{
    if (!m || !view || group>=m->groups) return 1;
    size_t first=group*m->width, count=m->blocks-first;
    if(count>m->width) count=m->width;
    const double *base=chem_sparse_const_at(m,first,row,column);
    if(!base) return 1;
    *view=(ChemConstSparseGroup){base,count}; return 0;
}
void chem_sparse_group_fill(ChemSparseGroup view,double value)
{
    for(size_t i=0;i<view.count;++i) view.base[i]=value;
}
int chem_sparse_group_copy(ChemSparseGroup out,ChemConstSparseGroup in)
{
    if(out.count!=in.count || (out.count && (!out.base || !in.base))) return 1;
    for(size_t i=0;i<out.count;++i) out.base[i]=in.base[i];
    return 0;
}
void chem_sparse_fill(ChemSparseMatrix *m,double value)
{
    for(size_t i=0;i<m->storage_size;++i) m->data[i]=value;
}
void chem_sparse_add_diagonal(ChemSparseMatrix *m,double value)
{
    for(size_t b=0;b<m->blocks;++b) for(size_t i=0;i<m->order;++i) {
        double *p=chem_sparse_at(m,b,i,i); if (p) *p+=value;
    }
}
int chem_sparse_print(const ChemSparseMatrix *m,FILE *file,int nonzero_only)
{
    if (!m || !file) return 1;
    for(size_t b=0;b<m->blocks;++b) {
        if (fprintf(file,"Block %zu\n",b)<0) return 1;
        for(size_t r=0;r<m->order;++r) {
            for(size_t c=0;c<m->order;++c) {
                const double *p=chem_sparse_const_at(m,b,r,c);
                if (nonzero_only) { if (p && fprintf(file,"%zu, %zu, %g\n",r,c,*p)<0) return 1; }
                else if (fprintf(file,"%s%g",c ? "," : "",p ? *p : 0)<0) return 1;
            }
            if (!nonzero_only && fputc('\n',file)==EOF) return 1;
        }
    }
    return 0;
}

double *chem_sparse_at_vector_index(ChemSparseMatrix *m,size_t block,size_t index)
{
    if (!m || block>=m->blocks || index % m->width || index/m->width>=m->nonzeros) return NULL;
    return &m->data[(block/m->width)*m->nonzeros*m->width+index+block%m->width];
}
int chem_sparse_fill_block_view(ChemSparseMatrix *m,size_t row,size_t column,double value)
{
    if (!m || row>=m->order || column>=m->order || m->indices[row*m->order+column]==SIZE_MAX) return 1;
    for(size_t b=0;b<m->blocks;++b) *chem_sparse_at(m,b,row,column)=value;
    return 0;
}
int chem_sparse_block_to_vector(const ChemSparseMatrix *m,size_t row,size_t column,double *values,size_t size)
{
    if (!m || size!=m->blocks || (size && !values) || row>=m->order || column>=m->order || m->indices[row*m->order+column]==SIZE_MAX) return 1;
    for(size_t b=0;b<size;++b) values[b]=*chem_sparse_const_at(m,b,row,column);
    return 0;
}
int chem_sparse_vector_to_block(ChemSparseMatrix *m,size_t row,size_t column,const double *values,size_t size)
{
    if (!m || size!=m->blocks || (size && !values) || row>=m->order || column>=m->order || m->indices[row*m->order+column]==SIZE_MAX) return 1;
    for(size_t b=0;b<size;++b) *chem_sparse_at(m,b,row,column)=values[b];
    return 0;
}
int chem_sparse_copy_block_view(ChemSparseMatrix *out,size_t orow,size_t ocol,const ChemSparseMatrix *in,size_t irow,size_t icol)
{
    if (!out || !in || out->blocks!=in->blocks || orow>=out->order || ocol>=out->order || irow>=in->order || icol>=in->order ||
        out->indices[orow*out->order+ocol]==SIZE_MAX || in->indices[irow*in->order+icol]==SIZE_MAX) return 1;
    for(size_t b=0;b<out->blocks;++b) *chem_sparse_at(out,b,orow,ocol)=*chem_sparse_const_at(in,b,irow,icol);
    return 0;
}
int chem_sparse_function_init(ChemSparseFunction *f,ChemSparseMatrix *const *sparse,size_t ns,
                              ChemMatrix *const *dense,size_t nd,size_t nv,ChemBlockFunction fn,void *context)
{
    if (!f || !fn || !ns || !sparse || (nd && !dense) || !sparse[0]) return 1;
    size_t width=sparse[0]->width;
    for(size_t i=0;i<ns;++i) if (!sparse[i] || sparse[i]->width!=width) return 1;
    for(size_t i=0;i<nd;++i) if (!dense[i] || dense[i]->width!=width) return 1;
    *f=(ChemSparseFunction){width,ns,nd,nv,fn,context}; return 0;
}
int chem_sparse_function_call(const ChemSparseFunction *f,ChemSparseMatrix *const *sparse,size_t ns,
                              ChemMatrix *const *dense,size_t nd,ChemVector *const *vectors,size_t nv)
{
    if (!f || !f->function || ns!=f->sparse_count || nd!=f->dense_count || nv!=f->vector_count ||
        !ns || !sparse || (nd && !dense) || (nv && !vectors) || !sparse[0]) return 1;
    size_t blocks=sparse[0]->blocks;
    for(size_t i=0;i<ns;++i) if (!sparse[i] || sparse[i]->width!=f->width || sparse[i]->blocks!=blocks) return 1;
    for(size_t i=0;i<nd;++i) if (!dense[i] || dense[i]->width!=f->width || dense[i]->rows!=blocks) return 1;
    for(size_t i=0;i<nv;++i) if (!vectors[i] || vectors[i]->size!=blocks || (blocks && !vectors[i]->data)) return 1;
    for(size_t b=0;b<blocks;++b) if (f->function(b,sparse,dense,vectors,f->context)) return 1;
    return 0;
}
