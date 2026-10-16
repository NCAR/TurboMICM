#ifndef TEST_BACKEND_H
#define TEST_BACKEND_H
#include "chem_factor.h"
typedef struct { chem_sparse a,l,u; chem_array rhs; int inplace; } integration_backend;
static inline void integration_backend_destroy(integration_backend *b)
{ chem_sparse_destroy(&b->a); chem_sparse_destroy(&b->l); chem_sparse_destroy(&b->u); chem_array_destroy(&b->rhs); }
static inline chem_status integration_factor(chem_state *s,void *context)
{
    integration_backend *b=context; size_t n=s->model->n_species;
    for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) {
        size_t k=s->model->lu.map[i*n+j];
        if (k!=CHEM_NONE) for (size_t c=0;c<s->cells;++c) chem_sparse_set(&b->a,c,i,j,s->jacobian[k*s->cells+c]);
    }
    return b->inplace ? chem_sparse_lu_factor_inplace(&s->model->lu,&b->a) : chem_sparse_lu_factor(&s->model->lu,&b->a,&b->l,&b->u);
}
static inline chem_status integration_linear_solve(const chem_state *s,double *rhs,void *context)
{
    integration_backend *b=context; size_t n=s->model->n_species;
    for (size_t i=0;i<n;++i) for (size_t c=0;c<s->cells;++c) chem_array_set(&b->rhs,c,i,rhs[i*s->cells+c]);
    chem_status status=b->inplace ? chem_sparse_lu_solve_inplace(&s->model->lu,&b->a,&b->rhs) : chem_sparse_lu_solve(&s->model->lu,&b->l,&b->u,&b->rhs);
    if (status!=CHEM_OK) return status;
    for (size_t i=0;i<n;++i) for (size_t c=0;c<s->cells;++c) rhs[i*s->cells+c]=chem_array_get(&b->rhs,c,i);
    return CHEM_OK;
}
static inline chem_status integration_backend_create(integration_backend *b,chem_state *s,size_t width,chem_sparse_order order,int inplace)
{
    *b=(integration_backend){.inplace=inplace}; chem_sparse_builder builder={0};
    chem_status status=chem_sparse_builder_create(&builder,s->model->n_species); if (status!=CHEM_OK) return status;
    builder.blocks=s->cells;
    for (size_t i=0;i<builder.n;++i) for (size_t j=0;j<builder.n;++j)
        if (s->model->lu.map[i*builder.n+j]!=CHEM_NONE) chem_sparse_builder_add(&builder,i,j);
    status=chem_sparse_create(&b->a,&builder,order,width,CHEM_DOUBLE); chem_sparse_builder_destroy(&builder);
    if (status==CHEM_OK && !inplace) status=chem_sparse_lu_matrices(&s->model->lu,&b->a,&b->l,&b->u,0);
    if (status==CHEM_OK) status=chem_array_create(&b->rhs,s->cells,s->model->n_species,width,CHEM_DOUBLE,0);
    if (status==CHEM_OK) status=chem_state_set_linear_backend(s,integration_factor,integration_linear_solve,b);
    if (status!=CHEM_OK) integration_backend_destroy(b);
    return status;
}
typedef struct { int method; size_t width; chem_lu_algorithm lu; chem_sparse_order order; int inplace; double tolerance; } integration_variant;
#endif
