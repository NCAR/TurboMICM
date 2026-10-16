#include "chem_factor.h"
#include <stdlib.h>
#include <string.h>
chem_status chem_sparse_lu_create(chem_lu_plan *p,const chem_sparse *s,chem_lu_algorithm algorithm)
{
    if(!s->n || s->values.type!=CHEM_DOUBLE) return CHEM_INVALID;
    unsigned char *pattern=malloc(s->n*s->n);if(!pattern) return CHEM_MEMORY;
    for(size_t i=0;i<s->n*s->n;++i) pattern[i]=s->map[i]!=CHEM_NONE;
    chem_status status=chem_lu_create(p,s->n,pattern,algorithm);free(pattern);return status;
}
static chem_status make_factor(const chem_lu_plan *p,const chem_sparse *s,chem_sparse *out,double initial,int triangle)
{
    chem_sparse_builder b;chem_status status=chem_sparse_builder_create(&b,p->n);if(status) return status;
    b.blocks=s->values.rows;b.initial=initial;
    for(size_t i=0;i<p->n;++i)for(size_t j=0;j<p->n;++j)
        if(p->map[i*p->n+j]!=CHEM_NONE && (!triangle || (triangle==1 ? j<=i : j>=i))) chem_sparse_builder_add(&b,i,j);
    status=chem_sparse_create(out,&b,s->order,s->values.vector_length,CHEM_DOUBLE);
    chem_sparse_builder_destroy(&b);return status;
}
chem_status chem_sparse_lu_matrix(const chem_lu_plan *p,const chem_sparse *s,chem_sparse *combined,double initial)
{return make_factor(p,s,combined,initial,0);}
chem_status chem_sparse_lu_matrices(const chem_lu_plan *p,const chem_sparse *s,chem_sparse *lower,chem_sparse *upper,double initial)
{
    *lower=(chem_sparse){0};*upper=(chem_sparse){0};
    chem_status status=make_factor(p,s,lower,initial,1);if(status) return status;
    status=make_factor(p,s,upper,initial,2);if(status) chem_sparse_destroy(lower);return status;
}
static chem_status check_factor(const chem_lu_plan *p,const chem_sparse *s,int triangle)
{
    if(s->n!=p->n || s->values.type!=CHEM_DOUBLE) return CHEM_SIZE;
    for(size_t i=0;i<p->n;++i)for(size_t j=0;j<p->n;++j) {
        int want=p->map[i*p->n+j]!=CHEM_NONE && (!triangle || (triangle==1 ? j<=i : j>=i));
        if((s->map[i*p->n+j]!=CHEM_NONE)!=want) return CHEM_SIZE;
    }
    return CHEM_OK;
}
/* Copy only real blocks. Tail lanes are never factored, avoiding artificial
 * zero pivots in padding and keeping every group width valid. */
static void pack(const chem_lu_plan *p,const chem_sparse *s,size_t start,size_t cells,double *data)
{
    memset(data,0,p->nnz*cells*sizeof(double));
    for(size_t i=0;i<p->n;++i)for(size_t j=0;j<p->n;++j) {
        size_t k=p->map[i*p->n+j],source=s->map[i*p->n+j];
        if(k!=CHEM_NONE && source!=CHEM_NONE)
            for(size_t c=0;c<cells;++c) data[k*cells+c]=chem_array_get(&s->values,start+c,source);
    }
}
static void unpack(const chem_lu_plan *p,chem_sparse *s,size_t start,size_t cells,const double *data,int lower)
{
    for(size_t i=0;i<p->n;++i)for(size_t j=0;j<p->n;++j) {
        size_t target=s->map[i*p->n+j];if(target==CHEM_NONE) continue;
        size_t k=p->map[i*p->n+j];
        for(size_t c=0;c<cells;++c) chem_array_set(&s->values,start+c,target,lower && i==j ? 1 : data[k*cells+c]);
    }
}
static chem_status factor(const chem_lu_plan *p,const chem_sparse *a,chem_sparse *l,chem_sparse *u)
{
    if(a->n!=p->n || a->values.type!=CHEM_DOUBLE) return CHEM_SIZE;
    if(check_factor(p,l,u ? 1 : 0) || (u && check_factor(p,u,2))) return CHEM_SIZE;
    if(l->values.rows!=a->values.rows || (u && u->values.rows!=a->values.rows)) return CHEM_SIZE;
    size_t width=a->values.vector_length ? a->values.vector_length : 1;
    if(p->nnz==SIZE_MAX || width>SIZE_MAX/(p->nnz+1)/sizeof(double)) return CHEM_SIZE;
    double *data=malloc((p->nnz+1)*width*sizeof(double));if(!data) return CHEM_MEMORY;
    for(size_t start=0;start<a->values.rows;start+=width) {
        size_t cells=a->values.rows-start<width ? a->values.rows-start : width;
        pack(p,a,start,cells,data);chem_lu_factor(p,cells,data,data+p->nnz*width);
        unpack(p,l,start,cells,data,u!=NULL);if(u) unpack(p,u,start,cells,data,0);
    }
    free(data);return CHEM_OK;
}
chem_status chem_sparse_lu_factor_inplace(const chem_lu_plan *p,chem_sparse *a)
{return factor(p,a,a,NULL);}
chem_status chem_sparse_lu_factor(const chem_lu_plan *p,const chem_sparse *a,chem_sparse *l,chem_sparse *u)
{return factor(p,a,l,u);}
static chem_status solve(const chem_lu_plan *p,const chem_sparse *l,const chem_sparse *u,chem_array *rhs)
{
    if(check_factor(p,l,u ? 1 : 0) || (u && check_factor(p,u,2))) return CHEM_SIZE;
    if(rhs->rows!=l->values.rows || rhs->columns!=p->n || rhs->type!=CHEM_DOUBLE || (u && u->values.rows!=rhs->rows)) return CHEM_SIZE;
    size_t width=rhs->vector_length ? rhs->vector_length : 1;
    if(p->nnz>SIZE_MAX-p->n || width>SIZE_MAX/(p->nnz+p->n)/sizeof(double)) return CHEM_SIZE;
    double *data=malloc((p->nnz+p->n)*width*sizeof(double));if(!data) return CHEM_MEMORY;
    double *x=data+p->nnz*width;
    for(size_t start=0;start<rhs->rows;start+=width) {
        size_t cells=rhs->rows-start<width ? rhs->rows-start : width;
        pack(p,l,start,cells,data);
        if(u) for(size_t i=0;i<p->n;++i)for(size_t j=i;j<p->n;++j) {
            size_t k=p->map[i*p->n+j];if(k!=CHEM_NONE)
                for(size_t c=0;c<cells;++c)data[k*cells+c]=chem_array_get(&u->values,start+c,u->map[i*p->n+j]);
        }
        for(size_t i=0;i<p->n;++i)for(size_t c=0;c<cells;++c)x[i*cells+c]=chem_array_get(rhs,start+c,i);
        chem_lu_solve(p,cells,data,x);
        for(size_t i=0;i<p->n;++i)for(size_t c=0;c<cells;++c)chem_array_set(rhs,start+c,i,x[i*cells+c]);
    }
    free(data);return CHEM_OK;
}
chem_status chem_sparse_lu_solve_inplace(const chem_lu_plan *p,const chem_sparse *a,chem_array *rhs)
{return solve(p,a,NULL,rhs);}
chem_status chem_sparse_lu_solve(const chem_lu_plan *p,const chem_sparse *l,const chem_sparse *u,chem_array *rhs)
{return solve(p,l,u,rhs);}
