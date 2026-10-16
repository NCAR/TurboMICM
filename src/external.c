#include "external_internal.h"
#include <stdlib.h>
#include <string.h>

typedef struct {
    chem_external_kernel forcing,jacobian;
    chem_status (*finalize)(const chem_model *,void *);
    void (*update)(chem_state *,void *);
    void *context;
    void (*destroy)(void *);
} external_binding;
typedef struct { size_t row,column; } index_pair;
typedef struct chem_composition {
    chem_mechanism mechanism;
    size_t base_species,base_constraints,base_parameters,n_bindings,n_elements;
    chem_species *species;
    chem_constraint *constraints;
    chem_named_value *parameters;
    external_binding *bindings;
    index_pair *elements;
} chem_composition;
static char *copy_string(const char *s)
{ if (!s) return NULL; size_t n=strlen(s)+1; char *p=malloc(n); if (p) memcpy(p,s,n); return p; }
static size_t species_index(const chem_composition *c,const char *name)
{
    if (name) for (size_t i=0;i<c->mechanism.n_species;++i) if (!strcmp(c->species[i].name,name)) return i;
    return CHEM_NONE;
}
static void empty_constraint(const chem_state *s,const double *y,double *out,void *context)
{ (void)s; (void)y; (void)out; (void)context; }
void chem_external_destroy(chem_model *m)
{
    chem_composition *c=m->composition; if (!c) return;
    if (c->species) for (size_t i=c->base_species;i<c->mechanism.n_species;++i) free((void *)c->species[i].name);
    if (c->constraints) for (size_t i=c->base_constraints;i<c->mechanism.n_constraints;++i) free((void *)c->constraints[i].dependencies);
    if (c->parameters) for (size_t i=c->base_parameters;i<c->mechanism.n_parameters;++i) free((void *)c->parameters[i].name);
    if (c->bindings) for (size_t i=0;i<c->n_bindings;++i) if (c->bindings[i].destroy && c->bindings[i].context) c->bindings[i].destroy(c->bindings[i].context);
    free(c->species); free(c->constraints); free(c->parameters); free(c->bindings); free(c->elements); free(c); m->composition=NULL;
}
chem_status chem_external_prepare(chem_model *m,const chem_mechanism **description)
{
    const chem_mechanism *base=*description;
    if (!base->n_external) return CHEM_OK;
    if (!base->external || (base->n_constraints && !base->constraints) || (base->n_parameters && !base->parameters)) return CHEM_INVALID;
    chem_composition *c=calloc(1,sizeof(*c)); if (!c) return CHEM_MEMORY; m->composition=c;
    c->mechanism=*base; c->base_species=base->n_species; c->base_constraints=base->n_constraints; c->base_parameters=base->n_parameters;
    c->n_bindings=base->n_external;
    for (size_t i=0;i<base->n_external;++i) {
        const chem_external_model *e=base->external+i;
        if ((e->n_species && !e->species) || (e->n_parameters && !e->parameters) ||
            (e->n_process_elements && (!e->process_elements || !e->forcing || !e->jacobian)) ||
            (e->n_algebraic && (!e->algebraic || !e->residual || !e->constraint_jacobian)) ||
            (e->n_constraint_elements && !e->constraint_elements) || (e->clone && !e->destroy) ||
            (e->context_size && !e->context)) return CHEM_INVALID;
        if (SIZE_MAX-c->mechanism.n_species<e->n_species || SIZE_MAX-c->mechanism.n_constraints<e->n_algebraic ||
            SIZE_MAX-c->mechanism.n_parameters<e->n_parameters || SIZE_MAX-c->n_elements<e->n_process_elements) return CHEM_SIZE;
        c->mechanism.n_species+=e->n_species; c->mechanism.n_constraints+=e->n_algebraic;
        c->mechanism.n_parameters+=e->n_parameters; c->n_elements+=e->n_process_elements;
    }
    c->species=calloc(c->mechanism.n_species,sizeof(*c->species));
    c->constraints=calloc(c->mechanism.n_constraints+1,sizeof(*c->constraints));
    c->parameters=calloc(c->mechanism.n_parameters+1,sizeof(*c->parameters));
    c->bindings=calloc(c->n_bindings,sizeof(*c->bindings)); c->elements=calloc(c->n_elements+1,sizeof(*c->elements));
    if (!c->species || !c->constraints || !c->parameters || !c->bindings || !c->elements) return CHEM_MEMORY;
    memcpy(c->species,base->species,base->n_species*sizeof(*c->species));
    if (base->n_constraints) memcpy(c->constraints,base->constraints,base->n_constraints*sizeof(*c->constraints));
    if (base->n_parameters) memcpy(c->parameters,base->parameters,base->n_parameters*sizeof(*c->parameters));
    size_t ns=c->base_species,nc=c->base_constraints,np=c->base_parameters,ne=0;
    for (size_t i=0;i<base->n_external;++i) {
        const chem_external_model *e=base->external+i;
        for (size_t j=0;j<e->n_species;++j) {
            if (!e->species[j].name || e->species[j].type!=CHEM_SPECIES_VARIABLE) return CHEM_INVALID;
            c->species[ns]=e->species[j]; c->species[ns++].name=copy_string(e->species[j].name);
            if (!c->species[ns-1].name) return CHEM_MEMORY;
        }
    }
    for (size_t i=0;i<c->mechanism.n_species;++i) {
        if (!c->species[i].name) return CHEM_INVALID;
        for (size_t j=0;j<i;++j) if (!strcmp(c->species[i].name,c->species[j].name)) return CHEM_INVALID;
    }
    for (size_t i=0;i<base->n_external;++i) {
        const chem_external_model *e=base->external+i; external_binding *binding=c->bindings+i;
        *binding=(external_binding){e->forcing,e->jacobian,e->finalize,e->update,(void *)e->context,NULL};
        if (e->clone) { binding->context=e->clone(e->context); binding->destroy=e->destroy; if (!binding->context) return CHEM_MEMORY; }
        else if (e->context_size) {
            binding->context=malloc(e->context_size); binding->destroy=free; if (!binding->context) return CHEM_MEMORY;
            memcpy(binding->context,e->context,e->context_size);
        }
        for (size_t j=0;j<e->n_parameters;++j) {
            if (!e->parameters[j].name) return CHEM_INVALID;
            c->parameters[np]=e->parameters[j]; c->parameters[np++].name=copy_string(e->parameters[j].name);
            if (!c->parameters[np-1].name) return CHEM_MEMORY;
        }
        for (size_t j=0;j<e->n_process_elements;++j) {
            size_t row=species_index(c,e->process_elements[j].row),col=species_index(c,e->process_elements[j].column);
            if (row==CHEM_NONE || col==CHEM_NONE) return CHEM_NOT_FOUND;
            if (c->species[row].type!=CHEM_SPECIES_VARIABLE || c->species[col].type!=CHEM_SPECIES_VARIABLE) return CHEM_INVALID;
            c->elements[ne++]=(index_pair){row,col};
        }
        for (size_t j=0;j<e->n_algebraic;++j) {
            size_t row=species_index(c,e->algebraic[j]); if (row==CHEM_NONE) return CHEM_NOT_FOUND;
            chem_constraint *con=c->constraints+nc++;
            *con=(chem_constraint){.name=c->species[row].name,.type=CHEM_CUSTOM_CONSTRAINT,.algebraic_species=row,
                .residual=j ? empty_constraint : e->residual,.jacobian=j ? empty_constraint : e->constraint_jacobian,.context=binding->context};
            size_t *deps=calloc(e->n_constraint_elements+1,sizeof(size_t)); if (!deps) return CHEM_MEMORY; con->dependencies=deps;
            for (size_t k=0;k<e->n_constraint_elements;++k) {
                size_t rr=species_index(c,e->constraint_elements[k].row),col=species_index(c,e->constraint_elements[k].column);
                if (rr==CHEM_NONE || col==CHEM_NONE) return CHEM_NOT_FOUND;
                if (rr==row) deps[con->n_dependencies++]=col;
            }
        }
    }
    c->mechanism.species=c->species; c->mechanism.constraints=c->constraints; c->mechanism.parameters=c->parameters;
    /* A generated base kernel's species and sparse offsets change when composed. */
    c->mechanism.forcing=NULL; c->mechanism.jacobian=NULL; c->mechanism.n_external=0; c->mechanism.external=NULL;
    *description=&c->mechanism; return CHEM_OK;
}
void chem_external_pattern(const chem_model *m,unsigned char *pattern)
{
    const chem_composition *c=m->composition; if (!c) return;
    for (size_t i=0;i<c->n_elements;++i) {
        size_t row=m->spec_to_var[c->elements[i].row],col=m->spec_to_var[c->elements[i].column];
        pattern[row*m->n_species+col]=1;
    }
}
int chem_external_species_used(const chem_model *m,size_t species)
{
    const chem_composition *c=m->composition; if (!c) return 0;
    if (species>=c->base_species) return 1;
    for (size_t i=0;i<c->n_elements;++i) if (c->elements[i].row==species || c->elements[i].column==species) return 1;
    return 0;
}
chem_status chem_external_finalize(chem_model *m)
{
    const chem_composition *c=m->composition; if (!c) return CHEM_OK;
    for (size_t i=0;i<c->n_bindings;++i) if (c->bindings[i].finalize) {
        chem_status status=c->bindings[i].finalize(m,c->bindings[i].context); if (status!=CHEM_OK) return status;
    }
    return CHEM_OK;
}
void chem_external_update(chem_state *s)
{
    const chem_composition *c=s->model->composition; if (!c) return;
    for (size_t i=0;i<c->n_bindings;++i) if (c->bindings[i].update) c->bindings[i].update(s,c->bindings[i].context);
}
const void *chem_model_external_context(const chem_model *m,size_t i)
{ return m && m->composition && i<m->composition->n_bindings ? m->composition->bindings[i].context : NULL; }
static void external_evaluate(const chem_state *s,const double *y,double *out,int jacobian)
{
    const chem_model *m=s->model; const chem_composition *c=m->composition; if (!c) return;
    /* Preserve caller-provided values in masked rows while allowing callbacks
     * to use every declared sparse element, including masked process entries. */
    if (s->external_work) for (size_t row=0;row<m->n_species;++row) if (m->algebraic[row]) {
        size_t begin=jacobian ? m->lu.row_start[row] : row,end=jacobian ? m->lu.row_start[row+1] : row+1;
        memcpy(s->external_work+begin*s->cells,out+begin*s->cells,(end-begin)*s->cells*sizeof(double));
    }
    for (size_t i=0;i<c->n_bindings;++i) {
        chem_external_kernel kernel=jacobian ? c->bindings[i].jacobian : c->bindings[i].forcing;
        if (kernel) kernel(s,y,out,c->bindings[i].context);
    }
    if (s->external_work) for (size_t row=0;row<m->n_species;++row) if (m->algebraic[row]) {
        size_t begin=jacobian ? m->lu.row_start[row] : row,end=jacobian ? m->lu.row_start[row+1] : row+1;
        memcpy(out+begin*s->cells,s->external_work+begin*s->cells,(end-begin)*s->cells*sizeof(double));
    }
}
void chem_state_forcing(const chem_state *s,const double *y,double *out)
{ s->model->forcing(s->model,s->cells,s->rates,y,out); external_evaluate(s,y,out,0); }
void chem_state_jacobian(const chem_state *s,const double *y,double *out)
{ s->model->jacobian(s->model,s->cells,s->rates,y,out); external_evaluate(s,y,out,1); }
