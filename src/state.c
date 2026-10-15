/* SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <string.h>

void chem_state_free(ChemState *s)
{
    if (!s) return;
    free(s->conditions); free(s->values); free(s->custom); free(s->rates);
    free(s->absolute_tolerance); free(s->scratch); free(s);
}
ChemStatus chem_state_create(const ChemMechanism *m, size_t cells, ChemState **out)
{
    if (!m || !out) return CHEM_INVALID;
    *out = NULL;
    ChemState *s = calloc(1, sizeof(*s));
    if (!s) return CHEM_MEMORY;
    s->mechanism = m; s->cells = cells; s->relative_tolerance = 1e-6;
    if ((m->variable_count && cells > SIZE_MAX / m->variable_count) ||
        (m->reaction_count && cells > SIZE_MAX / m->reaction_count) ||
        (m->custom_count && cells > SIZE_MAX / m->custom_count)) {
        free(s); return CHEM_MEMORY;
    }
    s->scratch = chem_alloc(cells, sizeof(double));
    s->conditions = chem_alloc(cells, sizeof(*s->conditions));
    s->values = chem_alloc(cells * m->variable_count, sizeof(double));
    s->rates = chem_alloc(cells * m->reaction_count, sizeof(double));
    s->custom = chem_alloc(cells * m->custom_count, sizeof(double));
    s->absolute_tolerance = chem_alloc(m->variable_count, sizeof(double));
    if (!s->conditions || !s->values || !s->rates || !s->custom || !s->absolute_tolerance || !s->scratch) {
        chem_state_free(s); return CHEM_MEMORY;
    }
    for (size_t i = 0; i < m->variable_count; ++i) {
        double a = m->species[m->variable_species[i]].absolute_tolerance;
        s->absolute_tolerance[i] = a > 0 ? a : 1e-3;
    }
    *out = s; return CHEM_OK;
}
ChemStatus chem_state_copy(const ChemState *s, ChemState **out)
{
    if (!s) return CHEM_INVALID;
    ChemStatus status = chem_state_create(s->mechanism, s->cells, out);
    if (status) return status;
    ChemState *d = *out;
    const ChemMechanism *m = s->mechanism;
    memcpy(d->conditions, s->conditions, s->cells * sizeof(*s->conditions));
    memcpy(d->values, s->values, s->cells * m->variable_count * sizeof(double));
    memcpy(d->custom, s->custom, s->cells * m->custom_count * sizeof(double));
    memcpy(d->rates, s->rates, s->cells * m->reaction_count * sizeof(double));
    memcpy(d->absolute_tolerance, s->absolute_tolerance, m->variable_count * sizeof(double));
    d->relative_tolerance = s->relative_tolerance;
    return CHEM_OK;
}
static ChemStatus set_field(double *field, size_t cells, size_t index, const double *v, size_t count)
{
    if (count != cells) return CHEM_DIMENSION;
    if (count && !v) return CHEM_INVALID;
    for (size_t c = 0; c < cells; ++c) field[index * cells + c] = v[c];
    return CHEM_OK;
}
ChemStatus chem_state_set_concentration(ChemState *s, const char *name, const double *v, size_t count)
{
    size_t i;
    if (!s) return CHEM_INVALID;
    ChemStatus status = chem_species_index(s->mechanism, name, &i);
    return status ? status : set_field(s->values, s->cells, i, v, count);
}
ChemStatus chem_state_set_custom(ChemState *s, const char *name, const double *v, size_t count)
{
    if (!s) return CHEM_INVALID;
    size_t i;
    ChemStatus status = chem_custom_index(s->mechanism, name, &i);
    if (status) return status;
    /* A label can intentionally drive several reactions. */
    for (; i < s->mechanism->custom_count; ++i) {
        if (!strcmp(s->mechanism->custom_names[i], name)) {
            status = set_field(s->custom, s->cells, i, v, count);
            if (status) return status;
        }
    }
    return CHEM_OK;
}
void chem_update_rates(ChemState *s)
{
    const ChemMechanism *m = s->mechanism;
    size_t n = s->cells;
    for (size_t i = 0; i < m->reaction_count; ++i) {
        const ChemProcess *p = &m->processes[i];
        const ChemRate *r = &p->rate;
        for (size_t c = 0; c < n; ++c) {
            double custom[2] = {0, 0};
            if (r->type == CHEM_USER || r->type == CHEM_SURFACE) custom[0] = s->custom[r->custom_index * n + c];
            if (r->type == CHEM_SURFACE) custom[1] = s->custom[(r->custom_index + 1) * n + c];
            s->rates[i * n + c] = chem_rate_evaluate(r, &s->conditions[c], custom);
        }
        if (p->parameter_count) for (size_t c = 0; c < n; ++c) {
            double mult = 1;
            for (size_t j = 0; j < p->parameter_count; ++j) {
                const ChemSpecies *sp = &m->species[p->parameters[j]];
                mult *= sp->parameter(&s->conditions[c], sp->parameter_context);
            }
            s->rates[i * n + c] *= mult;
        }
    }
    chem_models_update(s,0);
    chem_constraints_update_parameters(s);
}
void chem_add_kinetic_forcing(const ChemState *s, const double *restrict y, double *restrict f)
{
    const ChemMechanism *m = s->mechanism;
    size_t n = s->cells;
    double *restrict work = s->scratch;
    for (size_t i = 0; i < m->reaction_count; ++i) {
        const ChemProcess *p = &m->processes[i];
        memcpy(work, s->rates + i*n, n*sizeof(double));
        for (size_t j=0; j<p->reactant_count; ++j) {
            const double *v=y+p->reactants[j]*n;
            #pragma omp simd
            for (size_t c=0; c<n; ++c) work[c]*=v[c];
        }
        for (size_t j=0; j<p->reactant_count; ++j) {
            if (m->algebraic[p->reactants[j]]) continue;
            double *restrict dst=f+p->reactants[j]*n;
            #pragma omp simd
            for (size_t c=0; c<n; ++c) dst[c]-=work[c];
        }
        for (size_t j=0; j<p->product_count; ++j) {
            if (m->algebraic[p->products[j].species]) continue;
            double *restrict dst=f+p->products[j].species*n;
            double yield=p->products[j].yield;
            #pragma omp simd
            for (size_t c=0; c<n; ++c) dst[c]+=yield*work[c];
        }
    }
}
void chem_forcing(const ChemState *s,const double *y,double *f)
{
    memset(f,0,s->mechanism->variable_count*s->cells*sizeof(double));
    chem_add_kinetic_forcing(s,y,f);
    if (s->mechanism->model_count) {
        chem_external_add_forcing(s,y,f);
        for (size_t r=0;r<s->mechanism->variable_count;++r) if (s->mechanism->algebraic[r])
            memset(f+r*s->cells,0,s->cells*sizeof(double));
    }
    chem_constraints_forcing(s,y,f);
}
static void kinetic_jacobian_add(const ChemState *s, const double *restrict y, double *restrict jac,double sign)
{
    const ChemMechanism *m = s->mechanism;
    size_t n = s->cells, nv = m->variable_count;
    double *restrict work=s->scratch;
    for (size_t i = 0; i < m->reaction_count; ++i) {
        const ChemProcess *p = &m->processes[i];
        for (size_t j = 0; j < p->reactant_count; ++j) {
            size_t independent = p->reactants[j];
            memcpy(work,s->rates+i*n,n*sizeof(double));
            int skipped=0;
            for (size_t k=0;k<p->reactant_count;++k) {
                if (!skipped && p->reactants[k]==independent) { skipped=1; continue; }
                const double *v=y+p->reactants[k]*n;
                #pragma omp simd
                for (size_t c=0;c<n;++c) work[c]*=v[c];
            }
            for (size_t k=0;k<p->reactant_count;++k) {
                if (m->algebraic[p->reactants[k]]) continue;
                double *restrict dst=jac+(p->reactants[k]*nv+independent)*n;
                #pragma omp simd
                for (size_t c=0;c<n;++c) dst[c]-=sign*work[c];
            }
            for (size_t k=0;k<p->product_count;++k) {
                if (m->algebraic[p->products[k].species]) continue;
                double *restrict dst=jac+(p->products[k].species*nv+independent)*n;
                double yield=p->products[k].yield;
                #pragma omp simd
                for (size_t c=0;c<n;++c) dst[c]+=sign*yield*work[c];
            }
        }
    }
}
void chem_subtract_kinetic_jacobian(const ChemState *s,const double *y,double *jac)
{ kinetic_jacobian_add(s,y,jac,-1); }
void chem_jacobian(const ChemState *s,const double *y,double *jac)
{
    size_t nv=s->mechanism->variable_count;
    memset(jac,0,nv*nv*s->cells*sizeof(double));
    kinetic_jacobian_add(s,y,jac,1);
    if (s->mechanism->model_count) {
        chem_models_jacobian(s,y,jac,NULL,1,0);
        chem_clear_algebraic_jacobian(s,jac,NULL);
    }
    chem_constraints_jacobian_add(s,y,jac,NULL,1);
}

ChemStatus chem_state_assign(ChemState **destination,const ChemState *source)
{
    if (!destination || !source) return CHEM_INVALID;
    if (*destination==source) return CHEM_OK;
    ChemState *copy=NULL;ChemStatus status=chem_state_copy(source,&copy);
    if (status) return status;
    chem_state_free(*destination);*destination=copy;return CHEM_OK;
}
void chem_state_move(ChemState **destination,ChemState **source)
{
    if (!destination || !source || destination==source) return;
    chem_state_free(*destination);*destination=*source;*source=NULL;
}
ChemStatus chem_state_fill_concentration(ChemState *s,const char *name,double value)
{
    if (!s) return CHEM_INVALID;
    size_t index;ChemStatus status=chem_species_index(s->mechanism,name,&index);if (status) return status;
    for (size_t c=0;c<s->cells;++c) s->values[index*s->cells+c]=value;
    return CHEM_OK;
}
ChemStatus chem_state_fill_custom(ChemState *s,const char *name,double value)
{
    if (!s) return CHEM_INVALID;
    size_t index;ChemStatus status=chem_custom_index(s->mechanism,name,&index);if (status) return status;
    for (;index<s->mechanism->custom_count;++index) if (!strcmp(s->mechanism->custom_names[index],name))
        for (size_t c=0;c<s->cells;++c) s->custom[index*s->cells+c]=value;
    return CHEM_OK;
}
ChemStatus chem_state_set_concentrations(ChemState *s,size_t count,const ChemNamedValues *values)
{
    if (!s || (count && !values)) return CHEM_INVALID;
    for (size_t i=0;i<count;++i) {
        ChemStatus status=chem_state_set_concentration(s,values[i].name,values[i].values,values[i].count);
        if (status) return status;
    }
    return CHEM_OK;
}
ChemStatus chem_state_set_custom_parameters(ChemState *s,size_t count,const ChemNamedValues *values)
{
    if (!s || (count && !values)) return CHEM_INVALID;
    for (size_t i=0;i<count;++i) {
        ChemStatus status=chem_state_set_custom(s,values[i].name,values[i].values,values[i].count);
        if (status) return status;
    }
    return CHEM_OK;
}
ChemStatus chem_state_set_custom_matrix(ChemState *s,size_t cells,size_t parameters,const double *values)
{
    if (!s) return CHEM_INVALID;
    if (cells!=s->cells || parameters!=s->mechanism->custom_count) return CHEM_DIMENSION;
    if (cells && parameters && !values) return CHEM_INVALID;
    for (size_t c=0;c<cells;++c) for (size_t p=0;p<parameters;++p) s->custom[p*cells+c]=values[c*parameters+p];
    return CHEM_OK;
}
const double *chem_state_concentration_const_data(const ChemState *s,const char *name)
{
    if (!s) return NULL;
    size_t index;if (chem_species_index(s->mechanism,name,&index)) return NULL;
    return s->values+index*s->cells;
}
double *chem_state_concentration_data(ChemState *s,const char *name)
{
    if (!s) return NULL;
    size_t index;if (chem_species_index(s->mechanism,name,&index)) return NULL;
    return s->values+index*s->cells;
}
