/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Linear and Van't Hoff equilibrium constraints translated to C. */
#include "internal.h"
#include <math.h>
#include <string.h>

ChemVantHoff chem_vant_hoff(double k,double h)
{ return (ChemVantHoff){k,h,CHEM_GAS_CONSTANT,298.15}; }
ChemStatus chem_constraint_validate(const ChemConstraint *c)
{
    if (!c || !c->name || !c->algebraic_species ||
        (c->type!=CHEM_LINEAR_CONSTRAINT && c->type!=CHEM_EQUILIBRIUM_CONSTRAINT)) return CHEM_INVALID;
    if (!c->reactant_count) return CHEM_EMPTY_REACTANTS;
    if (!c->reactants || (c->product_count && !c->products)) return CHEM_INVALID;
    if (c->type==CHEM_EQUILIBRIUM_CONSTRAINT && !c->product_count) return CHEM_EMPTY_PRODUCTS;
    for (size_t j=0;j<c->reactant_count+c->product_count;++j) {
        const ChemConstraintTerm *t=j<c->reactant_count?c->reactants+j:c->products+j-c->reactant_count;
        if (!t->species) return CHEM_INVALID;
        if (c->type==CHEM_EQUILIBRIUM_CONSTRAINT && !(t->coefficient>0)) return CHEM_BAD_STOICHIOMETRY;
    }
    if (c->type==CHEM_EQUILIBRIUM_CONSTRAINT && !(c->equilibrium.reference_constant>0)) return CHEM_BAD_EQUILIBRIUM;
    return CHEM_OK;
}
static char *duplicate(const char *s)
{
    size_t n=strlen(s)+1;char *p=malloc(n);
    if (p) memcpy(p,s,n);
    return p;
}
void chem_constraints_free(ChemBoundConstraint *constraints,size_t n)
{
    if (!constraints) return;
    for (size_t i=0;i<n;++i) {
        ChemConstraint *c=&constraints[i].definition;
        free((void *)c->name);free((void *)c->algebraic_species);
        if (c->reactants) for (size_t j=0;j<c->reactant_count;++j) free((void *)c->reactants[j].species);
        if (c->products) for (size_t j=0;j<c->product_count;++j) free((void *)c->products[j].species);
        free((void *)c->reactants);free((void *)c->products);
        free(constraints[i].reactants);free(constraints[i].products);
    }
    free(constraints);
}
ChemStatus chem_mechanism_set_constraints(ChemMechanism *m,size_t count,const ChemConstraint *definitions)
{
    if (!m || (count && !definitions)) return CHEM_INVALID;
    ChemBoundConstraint *constraints=chem_alloc(count,sizeof(*constraints));
    unsigned char *algebraic=chem_alloc(m->variable_count,1);
    size_t custom_count=m->base_custom_count;
    for (size_t i=0;i<count;++i) custom_count+=definitions[i].type==CHEM_EQUILIBRIUM_CONSTRAINT;
    char **names=chem_alloc(custom_count,sizeof(*names));
    ChemStatus status=CHEM_MEMORY;
    if (!constraints || !algebraic || !names) goto fail;
    for (size_t i=0;i<m->base_custom_count;++i) if (!(names[i]=duplicate(m->custom_names[i]))) goto fail;
    for (size_t i=0;i<m->model_count;++i)
        for (size_t j=0;j<m->models[i].binding.algebraic_count;++j)
            algebraic[m->models[i].binding.algebraic[j]]=1;
    size_t parameter=m->base_custom_count;
    for (size_t i=0;i<count;++i) {
        const ChemConstraint *source=&definitions[i];
        status=chem_constraint_validate(source);if (status) goto fail;
        ChemBoundConstraint *bound=&constraints[i];ChemConstraint *c=&bound->definition;
        status=chem_species_index(m,source->algebraic_species,&bound->row);if (status) goto fail;
        if (algebraic[bound->row]) { status=CHEM_INVALID;goto fail; }
        algebraic[bound->row]=1;
        *c=*source;c->name=NULL;c->algebraic_species=NULL;c->reactants=NULL;c->products=NULL;
        c->name=duplicate(source->name);c->algebraic_species=duplicate(source->algebraic_species);
        ChemConstraintTerm *lhs=chem_alloc(c->reactant_count,sizeof(*lhs)),*rhs=chem_alloc(c->product_count,sizeof(*rhs));
        c->reactants=lhs;c->products=rhs;
        bound->reactants=chem_alloc(c->reactant_count,sizeof(size_t));bound->products=chem_alloc(c->product_count,sizeof(size_t));
        status=CHEM_MEMORY;
        if (!c->name || !c->algebraic_species || !lhs || !rhs || !bound->reactants || !bound->products) goto fail;
        for (size_t j=0;j<c->reactant_count+c->product_count;++j) {
            int reactant=j<c->reactant_count;size_t k=reactant?j:j-c->reactant_count;
            const ChemConstraintTerm *from=reactant?source->reactants+k:source->products+k;
            ChemConstraintTerm *to=reactant?lhs+k:rhs+k;
            to->coefficient=from->coefficient;to->species=duplicate(from->species);
            status=CHEM_MEMORY;if (!to->species) goto fail;
            status=chem_species_index(m,from->species,(reactant?bound->reactants:bound->products)+k);
            if (status) goto fail;
        }
        bound->parameter=SIZE_MAX;
        if (c->type==CHEM_EQUILIBRIUM_CONSTRAINT) {
            bound->parameter=parameter;
            names[parameter++]=duplicate(c->name);
            status=CHEM_MEMORY;if (!names[parameter-1]) goto fail;
        }
    }
    chem_constraints_free(m->constraints,m->constraint_count);free(m->algebraic);
    for (size_t i=0;i<m->custom_count;++i) free(m->custom_names[i]);
    free(m->custom_names);
    m->constraints=constraints;m->constraint_count=count;m->algebraic=algebraic;
    m->custom_names=names;m->custom_count=custom_count;
    return CHEM_OK;
fail:
    chem_constraints_free(constraints,count);free(algebraic);
    if (names) for (size_t i=0;i<custom_count;++i) free(names[i]);
    free(names);
    return status;
}
size_t chem_constraint_count(const ChemMechanism *m) { return m->constraint_count+m->external_constraint_count; }
const ChemConstraint *chem_constraint_definition(const ChemMechanism *m,size_t i)
{ return i<m->constraint_count?&m->constraints[i].definition:NULL; }
size_t chem_constraint_row(const ChemMechanism *m,size_t i)
{
    if (i<m->constraint_count) return m->constraints[i].row;
    i-=m->constraint_count;
    for (size_t model=0;model<m->model_count;++model) {
        const ChemModelBinding *b=&m->models[model].binding;
        if (i<b->algebraic_count) return b->algebraic[i];
        i-=b->algebraic_count;
    }
    return SIZE_MAX;
}
int chem_is_algebraic(const ChemMechanism *m,size_t v)
{ return v<m->variable_count && m->algebraic[v]; }
int chem_constraint_has_element(const ChemMechanism *m,size_t row,size_t col)
{
    for (size_t i=0;i<m->constraint_count;++i) {
        const ChemBoundConstraint *b=&m->constraints[i];if (b->row!=row) continue;
        if (col==row) return 1;
        for (size_t j=0;j<b->definition.reactant_count;++j) if (b->reactants[j]==col) return 1;
        for (size_t j=0;j<b->definition.product_count;++j) if (b->products[j]==col) return 1;
    }
    return chem_external_has_element(m,row,col,1);
}
void chem_constraints_update_parameters(ChemState *s)
{
    const ChemMechanism *m=s->mechanism;
    chem_models_update(s,1);
    for (size_t i=0;i<m->constraint_count;++i) {
        const ChemBoundConstraint *b=&m->constraints[i];
        if (b->definition.type!=CHEM_EQUILIBRIUM_CONSTRAINT) continue;
        ChemVantHoff p=b->definition.equilibrium;
        for (size_t c=0;c<s->cells;++c)
            s->custom[b->parameter*s->cells+c]=p.reference_constant*exp((p.enthalpy/p.gas_constant)*
                (1/s->conditions[c].temperature-1/p.reference_temperature));
    }
}
static double product(const ChemConstraintTerm *terms,const size_t *indices,size_t count,
                      const double *y,size_t cells,size_t cell,size_t skip,double start)
{
    for (size_t j=0;j<count;++j) if (j!=skip) {
        double value=y[indices[j]*cells+cell];
        start*=pow(value>0?value:0,terms[j].coefficient);
    }
    return start;
}
void chem_constraints_forcing(const ChemState *s,const double *y,double *f)
{
    const ChemMechanism *m=s->mechanism;size_t cells=s->cells;
    for (size_t i=0;i<m->constraint_count;++i) {
        const ChemBoundConstraint *b=&m->constraints[i];const ChemConstraint *d=&b->definition;
        for (size_t c=0;c<cells;++c) {
            double value=0;
            if (d->type==CHEM_LINEAR_CONSTRAINT) {
                for (size_t j=0;j<d->reactant_count;++j) value+=d->reactants[j].coefficient*y[b->reactants[j]*cells+c];
                value-=d->constant;
            } else {
                double rp=product(d->reactants,b->reactants,d->reactant_count,y,cells,c,SIZE_MAX,s->custom[b->parameter*cells+c]);
                double pp=product(d->products,b->products,d->product_count,y,cells,c,SIZE_MAX,1);
                value=rp-pp;
            }
            f[b->row*cells+c]=value;
        }
    }
    chem_models_residual(s,y,f);
}
void chem_constraints_jacobian_add(const ChemState *s,const double *y,double *jac,const ChemSparse *sparse,double sign)
{
    const ChemMechanism *m=s->mechanism;size_t cells=s->cells,n=m->variable_count;
    for (size_t i=0;i<m->constraint_count;++i) {
        const ChemBoundConstraint *b=&m->constraints[i];const ChemConstraint *d=&b->definition;
        for (size_t j=0;j<d->reactant_count+d->product_count;++j) {
            int reactant=j<d->reactant_count;size_t k=reactant?j:j-d->reactant_count;
            const ChemConstraintTerm *terms=reactant?d->reactants:d->products;
            const size_t *indices=reactant?b->reactants:b->products;
            size_t count=reactant?d->reactant_count:d->product_count;
            size_t index=sparse?chem_sparse_index(sparse,b->row,indices[k]):b->row*n+indices[k];
            for (size_t c=0;c<cells;++c) {
                double partial=terms[k].coefficient;
                if (d->type==CHEM_EQUILIBRIUM_CONSTRAINT) {
                    double p=product(terms,indices,count,y,cells,c,k,reactant?s->custom[b->parameter*cells+c]:1);
                    double conc=y[indices[k]*cells+c],stoich=terms[k].coefficient;
                    partial=stoich==1?p:conc>0?stoich*p*pow(conc,stoich-1):0;
                    if (!reactant) partial=-partial;
                }
                jac[index*cells+c]+=sign*partial;
            }
        }
    }
    chem_models_jacobian(s,y,jac,sparse,sign,1);
}
void chem_constraints_subtract_jacobian(const ChemState *s,const double *y,double *jac)
{ chem_constraints_jacobian_add(s,y,jac,NULL,-1); }
