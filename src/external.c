/* SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <string.h>

static ChemStatus resolve_names(const ChemMechanism *m,size_t count,const char *const *names,
                                const size_t **out,int parameters)
{
    size_t *indices=chem_alloc(count,sizeof(*indices));
    if (!indices) return CHEM_MEMORY;
    *out=indices;
    for (size_t i=0;i<count;++i) {
        ChemStatus status=parameters?chem_custom_index(m,names[i],indices+i):chem_species_index(m,names[i],indices+i);
        if (status) return status;
    }
    return CHEM_OK;
}
static ChemStatus resolve_elements(const ChemMechanism *m,size_t count,const ChemNamedJacobianElement *names,
                                   const ChemJacobianElement **out)
{
    ChemJacobianElement *elements=chem_alloc(count,sizeof(*elements));
    if (!elements) return CHEM_MEMORY;
    *out=elements;
    for (size_t i=0;i<count;++i) {
        ChemStatus status=chem_species_index(m,names[i].row,&elements[i].row);
        if (!status) status=chem_species_index(m,names[i].column,&elements[i].column);
        if (status) return status;
    }
    return CHEM_OK;
}
void chem_models_free(ChemMechanism *m)
{
    if (!m->models) return;
    for (size_t i=0;i<m->model_count;++i) {
        ChemBoundModel *b=m->models+i;
        if (b->context) {
            if (b->callbacks.destroy) b->callbacks.destroy(b->context);
            else free(b->context);
        }
        free((void *)b->binding.variables);free((void *)b->binding.parameters);free((void *)b->binding.algebraic);
        free((void *)b->binding.process_elements);free((void *)b->binding.constraint_elements);
    }
    free(m->models);
}
static int valid(const ChemExternalModel *d)
{
    return !(d->additional_species_count && !d->additional_species) &&
        !(d->variable_count && !d->variables) && !(d->parameter_count && !d->parameters) &&
        !(d->algebraic_count && (!d->algebraic || !d->constraint_residual || !d->constraint_jacobian)) &&
        !(d->process_element_count && !d->process_elements) &&
        !(d->constraint_element_count && !d->constraint_elements) &&
        !(d->context && !d->clone && !d->context_size) &&
        !(!d->context && (d->clone || d->context_size)) && (!d->destroy || d->clone);
}
ChemStatus chem_mechanism_create_with_models(size_t ns,const ChemSpecies *species,size_t nr,const ChemReaction *reactions,
                                             size_t count,const ChemExternalModel *definitions,ChemMechanism **out)
{
    if (!out || (ns && !species) || (count && !definitions)) return CHEM_INVALID;
    *out=NULL;
    size_t total=ns,parameters=0;
    for (size_t i=0;i<count;++i) {
        const ChemExternalModel *d=definitions+i;
        if (!valid(d)) return CHEM_INVALID;
        if (d->additional_species_count>SIZE_MAX-total || d->parameter_count>SIZE_MAX-parameters) return CHEM_MEMORY;
        total+=d->additional_species_count;parameters+=d->parameter_count;
        for (size_t j=0;j<d->additional_species_count;++j) if (d->additional_species[j].parameter) return CHEM_INVALID;
    }
    ChemSpecies *all=chem_alloc(total,sizeof(*all));
    if (!all) return CHEM_MEMORY;
    if (ns) memcpy(all,species,ns*sizeof(*all));
    size_t at=ns;
    for (size_t i=0;i<count;++i) {
        const ChemExternalModel *d=definitions+i;
        if (d->additional_species_count) memcpy(all+at,d->additional_species,d->additional_species_count*sizeof(*all));
        at+=d->additional_species_count;
    }
    ChemMechanism *m=NULL;
    ChemStatus status=chem_mechanism_create(total,all,nr,reactions,&m);
    free(all);
    if (status) return status;
    status=CHEM_MEMORY;
    if (parameters>SIZE_MAX-m->custom_count) goto fail;
    char **names=chem_alloc(m->custom_count+parameters,sizeof(*names));
    if (!names) goto fail;
    memcpy(names,m->custom_names,m->custom_count*sizeof(*names));
    free(m->custom_names);m->custom_names=names;
    for (size_t i=0;i<count;++i) for (size_t j=0;j<definitions[i].parameter_count;++j) {
        const char *name=definitions[i].parameters[j];
        status=CHEM_INVALID;
        if (!name || !*name) goto fail;
        size_t existing;
        if (!chem_custom_index(m,name,&existing)) goto fail;
        size_t length=strlen(name)+1;
        char *copy=malloc(length);status=CHEM_MEMORY;
        if (!copy) goto fail;
        memcpy(copy,name,length);m->custom_names[m->custom_count++]=copy;
    }
    m->base_custom_count=m->custom_count;
    m->models=chem_alloc(count,sizeof(*m->models));
    if (!m->models) goto fail;
    m->model_count=count;
    for (size_t i=0;i<count;++i) {
        const ChemExternalModel *d=definitions+i;
        ChemBoundModel *b=m->models+i;
        ChemModelBinding *v=&b->binding;
        b->callbacks=(ChemExternalModel){.destroy=d->destroy,.update_process_parameters=d->update_process_parameters,
            .update_constraint_parameters=d->update_constraint_parameters,.add_forcing=d->add_forcing,
            .constraint_residual=d->constraint_residual,.process_jacobian=d->process_jacobian,
            .constraint_jacobian=d->constraint_jacobian,.finalize=d->finalize,.finalize_sparse=d->finalize_sparse};
        if (d->context) {
            b->context=d->clone?d->clone(d->context):malloc(d->context_size);
            status=CHEM_MEMORY;if (!b->context) goto fail;
            if (!d->clone) memcpy(b->context,d->context,d->context_size);
        }
        v->variable_count=d->variable_count;v->parameter_count=d->parameter_count;v->algebraic_count=d->algebraic_count;
        v->process_element_count=d->process_element_count;v->constraint_element_count=d->constraint_element_count;
        status=resolve_names(m,d->variable_count,d->variables,&v->variables,0);if (status) goto fail;
        status=resolve_names(m,d->parameter_count,d->parameters,&v->parameters,1);if (status) goto fail;
        status=resolve_names(m,d->algebraic_count,d->algebraic,&v->algebraic,0);if (status) goto fail;
        status=resolve_elements(m,d->process_element_count,d->process_elements,&v->process_elements);if (status) goto fail;
        status=resolve_elements(m,d->constraint_element_count,d->constraint_elements,&v->constraint_elements);if (status) goto fail;
        for (size_t j=0;j<v->algebraic_count;++j) {
            status=CHEM_INVALID;if (m->algebraic[v->algebraic[j]]) goto fail;
            m->algebraic[v->algebraic[j]]=1;++m->external_constraint_count;
        }
        for (size_t j=0;j<v->constraint_element_count;++j) {
            int declared=0;
            for (size_t k=0;k<v->algebraic_count;++k) declared|=v->algebraic[k]==v->constraint_elements[j].row;
            status=CHEM_INVALID;if (!declared) goto fail;
        }
    }
    *out=m;return CHEM_OK;
fail:
    chem_mechanism_free(m);return status;
}
size_t chem_external_model_count(const ChemMechanism *m) { return m->model_count; }
void chem_models_reorder(ChemMechanism *m,const size_t *inverse)
{
    for (size_t i=0;i<m->model_count;++i) {
        ChemModelBinding *b=&m->models[i].binding;
        size_t *variables=(size_t *)b->variables,*rows=(size_t *)b->algebraic;
        for (size_t j=0;j<b->variable_count;++j) variables[j]=inverse[variables[j]];
        for (size_t j=0;j<b->algebraic_count;++j) { rows[j]=inverse[rows[j]];m->algebraic[rows[j]]=1; }
        for (size_t kind=0;kind<2;++kind) {
            size_t count=kind?b->constraint_element_count:b->process_element_count;
            ChemJacobianElement *e=(ChemJacobianElement *)(kind?b->constraint_elements:b->process_elements);
            for (size_t j=0;j<count;++j) { e[j].row=inverse[e[j].row];e[j].column=inverse[e[j].column]; }
        }
    }
}
int chem_external_has_element(const ChemMechanism *m,size_t row,size_t col,int constraints)
{
    for (size_t i=0;i<m->model_count;++i) {
        const ChemModelBinding *b=&m->models[i].binding;
        const ChemJacobianElement *e=constraints?b->constraint_elements:b->process_elements;
        size_t count=constraints?b->constraint_element_count:b->process_element_count;
        for (size_t j=0;j<count;++j) if (e[j].row==row && e[j].column==col) return 1;
    }
    return 0;
}
void chem_models_update(ChemState *s,int constraints)
{
    const ChemMechanism *m=s->mechanism;
    for (size_t i=0;i<m->model_count;++i) {
        const ChemBoundModel *b=m->models+i;
        ChemModelUpdate fn=constraints?b->callbacks.update_constraint_parameters:b->callbacks.update_process_parameters;
        if (fn) fn(&b->binding,b->context,s);
    }
}
void chem_external_add_forcing(const ChemState *s,const double *y,double *f)
{
    const ChemMechanism *m=s->mechanism;
    for (size_t i=0;i<m->model_count;++i) {
        const ChemBoundModel *b=m->models+i;
        if (b->callbacks.add_forcing) b->callbacks.add_forcing(&b->binding,b->context,s,y,f);
    }
}
void chem_models_residual(const ChemState *s,const double *y,double *f)
{
    const ChemMechanism *m=s->mechanism;
    for (size_t i=0;i<m->model_count;++i) {
        const ChemBoundModel *b=m->models+i;
        if (b->callbacks.constraint_residual) b->callbacks.constraint_residual(&b->binding,b->context,s,y,f);
    }
}
ChemStatus chem_model_jacobian_add(ChemModelJacobian *j,size_t row,size_t col,size_t cell,double derivative)
{
    if (!j || row>=j->order || col>=j->order || cell>=j->cells) return CHEM_DIMENSION;
    size_t index=j->sparse?chem_sparse_index(j->sparse,row,col):row*j->order+col;
    if (index==SIZE_MAX) return CHEM_NOT_FOUND;
    j->values[index*j->cells+cell]+=j->sign*derivative;
    return CHEM_OK;
}
void chem_models_jacobian(const ChemState *s,const double *y,double *jac,const ChemSparse *sparse,double sign,int constraints)
{
    const ChemMechanism *m=s->mechanism;
    ChemModelJacobian view={m->variable_count,s->cells,sparse,jac,sign};
    for (size_t i=0;i<m->model_count;++i) {
        const ChemBoundModel *b=m->models+i;
        ChemModelDerivative fn=constraints?b->callbacks.constraint_jacobian:b->callbacks.process_jacobian;
        if (fn) fn(&b->binding,b->context,s,y,&view);
    }
}
void chem_external_subtract_jacobian(const ChemState *s,const double *y,double *jac)
{ chem_models_jacobian(s,y,jac,NULL,-1,0); }
void chem_clear_algebraic_jacobian(const ChemState *s,double *jac,const ChemSparse *sparse)
{
    const ChemMechanism *m=s->mechanism;size_t n=m->variable_count;
    for (size_t r=0;r<n;++r) if (m->algebraic[r]) for (size_t c=0;c<n;++c) {
        size_t index=sparse?chem_sparse_index(sparse,r,c):r*n+c;
        if (index!=SIZE_MAX) memset(jac+index*s->cells,0,s->cells*sizeof(double));
    }
}
