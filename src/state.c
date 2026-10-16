#include "chem.h"
#include <stdlib.h>
#include <string.h>

void chem_state_destroy(chem_state *s)
{
    free(s->y); free(s->rates); free(s->custom); free(s->atol); free(s->conditions);
    free(s->work); free(s->jacobian); free(s->scratch);
    memset(s,0,sizeof(*s));
}

chem_status chem_state_create(chem_state *s, const chem_model *m, size_t cells)
{
    memset(s,0,sizeof(*s));
    if (!m || !cells || cells>SIZE_MAX/(m->lu.nnz+10*m->n_species+m->n_custom+m->n_reactions+1)/sizeof(double)) return CHEM_INVALID;
    size_t n=m->n_species;
    s->model=m; s->cells=cells; s->rtol=1e-6;
    s->y=calloc(n*cells,sizeof(double)); s->rates=calloc((m->n_reactions+1)*cells,sizeof(double));
    s->custom=calloc((m->n_custom+1)*cells,sizeof(double)); s->atol=malloc(n*sizeof(double));
    s->conditions=calloc(cells,sizeof(chem_conditions)); s->work=calloc(9*n*cells,sizeof(double));
    s->jacobian=calloc(m->lu.nnz*cells,sizeof(double)); s->scratch=calloc(cells,sizeof(double));
    if (!s->y || !s->rates || !s->custom || !s->atol || !s->conditions || !s->work || !s->jacobian || !s->scratch) {
        chem_state_destroy(s); return CHEM_MEMORY;
    }
    for (size_t i=0; i<n; ++i) {
        const chem_species *sp=&m->mechanism->species[m->var_to_spec[i]];
        s->atol[i]=sp->absolute_tolerance > 0 ? sp->absolute_tolerance : 1e-3;
        for (size_t c=0; c<cells; ++c) s->y[i*cells+c]=sp->initial;
    }
    for (size_t c=0; c<cells; ++c) s->conditions[c]=m->mechanism->conditions;
    for (size_t i=0; i<m->mechanism->n_defaults; ++i) {
        const chem_named_value *v=&m->mechanism->defaults[i];
        chem_status status=chem_state_set_custom(s,v->name,&v->value,1);
        if (status!=CHEM_OK) { chem_state_destroy(s); return status; }
    }
    return CHEM_OK;
}

chem_status chem_state_copy(chem_state *dst, const chem_state *src)
{
    if (dst==src) return CHEM_OK;
    chem_state copy;
    chem_status rc=chem_state_create(&copy,src->model,src->cells);
    if (rc!=CHEM_OK) return rc;
    size_t c=src->cells, n=src->model->n_species;
    memcpy(copy.y,src->y,n*c*sizeof(double)); memcpy(copy.rates,src->rates,src->model->n_reactions*c*sizeof(double));
    memcpy(copy.custom,src->custom,src->model->n_custom*c*sizeof(double)); memcpy(copy.atol,src->atol,n*sizeof(double));
    memcpy(copy.conditions,src->conditions,c*sizeof(chem_conditions)); copy.rtol=src->rtol;
    copy.factor=src->factor; copy.solve_linear=src->solve_linear; copy.linear_context=src->linear_context;
    chem_state_destroy(dst); *dst=copy;
    return CHEM_OK;
}

chem_status chem_state_set_linear_backend(chem_state *s,
    chem_status (*factor)(chem_state *,void *),
    chem_status (*solve)(const chem_state *,double *,void *),void *context)
{
    if (!s || (!factor != !solve)) return CHEM_INVALID;
    s->factor=factor; s->solve_linear=solve; s->linear_context=context;
    return CHEM_OK;
}
chem_status chem_state_factor(chem_state *s)
{
    if (s->factor) return s->factor(s,s->linear_context);
    chem_lu_factor(&s->model->lu,s->cells,s->jacobian,s->scratch);
    return CHEM_OK;
}
chem_status chem_state_solve_linear(const chem_state *s,double *rhs)
{
    if (s->solve_linear) return s->solve_linear(s,rhs,s->linear_context);
    chem_lu_solve(&s->model->lu,s->cells,s->jacobian,rhs);
    return CHEM_OK;
}

static chem_status set_column(double *data, size_t cells, size_t index, const double *values, size_t count)
{
    if (index==CHEM_NONE) return CHEM_NOT_FOUND;
    if (!values || (count!=1 && count!=cells)) return CHEM_SIZE;
    for (size_t c=0; c<cells; ++c) data[index*cells+c]=values[count==1 ? 0 : c];
    return CHEM_OK;
}
chem_status chem_state_set_species(chem_state *s, const char *name, const double *values, size_t count)
{ return set_column(s->y,s->cells,chem_species_index(s->model,name),values,count); }
chem_status chem_state_set_custom(chem_state *s, const char *name, const double *values, size_t count)
{
    size_t first=chem_custom_index(s->model,name);
    chem_status status=set_column(s->custom,s->cells,first,values,count);
    if (status!=CHEM_OK) return status;
    /* A parameter label can be shared by several reactions. */
    for (size_t i=first+1; i<s->model->n_custom; ++i)
        if (!strcmp(s->model->custom_names[i],name)) set_column(s->custom,s->cells,i,values,count);
    return CHEM_OK;
}

void chem_calculate_rates(chem_state *s)
{
    const chem_model *m=s->model;
    size_t cells=s->cells;
    for (size_t r=0; r<m->n_reactions; ++r) {
        const chem_rate *rate=&m->rates[r];
        const chem_reaction *rx=&m->mechanism->reactions[r];
        for (size_t c=0; c<cells; ++c) {
            double custom[2]={0,0};
            if (rate->type==CHEM_USER || rate->type==CHEM_SURFACE) custom[0]=s->custom[m->custom_offset[r]*cells+c];
            if (rate->type==CHEM_SURFACE) custom[1]=s->custom[(m->custom_offset[r]+1)*cells+c];
            double k=chem_rate_calculate(rate,&s->conditions[c],custom);
            for (size_t j=0; j<rx->n_reactants; ++j) {
                const chem_species *sp=&m->mechanism->species[rx->reactants[j]];
                if (sp->type!=CHEM_SPECIES_VARIABLE) k*=chem_species_evaluate(sp,&s->conditions[c]);
            }
            s->rates[r*cells+c]=k;
        }
    }
    chem_update_constraint_parameters(s);
}

void chem_stats_add(chem_stats *a, const chem_stats *b)
{
    a->internal_steps+=b->internal_steps; a->accepted+=b->accepted; a->rejected+=b->rejected;
    a->function_calls+=b->function_calls; a->jacobian_updates+=b->jacobian_updates;
    a->constraint_init_iterations+=b->constraint_init_iterations;
    a->decompositions+=b->decompositions; a->solves+=b->solves; a->final_time+=b->final_time;
}
