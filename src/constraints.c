#include "internal.h"

#include <math.h>
#include <string.h>

static int validate_one(const ChemConstraint *c, size_t species_count)
{
    if (!c || c->algebraic_species < 0 || (size_t)c->algebraic_species >= species_count ||
        (c->kind != CHEM_LINEAR_CONSTRAINT && c->kind != CHEM_EQUILIBRIUM_CONSTRAINT) ||
        !c->reactant_count || !c->reactants) return CHEM_INVALID;
    if (c->kind == CHEM_EQUILIBRIUM_CONSTRAINT &&
        (!c->product_count || !c->products || !isfinite(c->equilibrium_constant) || c->equilibrium_constant <= 0))
        return CHEM_INVALID;
    if (c->kind == CHEM_LINEAR_CONSTRAINT && c->product_count) return CHEM_INVALID;
    for (int side=0; side<2; ++side) {
        const ChemTerm *terms=side ? c->products : c->reactants;
        size_t count=side ? c->product_count : c->reactant_count;
        for (size_t i=0; i<count; ++i)
            if (terms[i].species < 0 || (size_t)terms[i].species >= species_count ||
                !isfinite(terms[i].coefficient) ||
                (c->kind == CHEM_EQUILIBRIUM_CONSTRAINT && terms[i].coefficient <= 0)) return CHEM_INVALID;
    }
    return CHEM_OK;
}

int chem_validate_constraints(const ChemMechanism *m, const ChemConstraint *constraints, size_t count)
{
    if (!m || (count && !constraints)) return CHEM_INVALID;
    for (size_t i=0; i<count; ++i) {
        if (validate_one(&constraints[i],m->species_count)) return CHEM_INVALID;
        for (size_t j=0; j<i; ++j)
            if (constraints[i].algebraic_species == constraints[j].algebraic_species) return CHEM_INVALID;
    }
    return CHEM_OK;
}

static double equilibrium_constant(const ChemConstraint *c, double temperature)
{
    double gas=c->gas_constant == 0 ? 1.380649e-23*6.02214076e23 : c->gas_constant;
    double ref=c->reference_temperature == 0 ? 298.15 : c->reference_temperature;
    return c->equilibrium_constant*exp((c->enthalpy/gas)*(1.0/temperature-1.0/ref));
}

double chem_constraint_value(const ChemConstraint *c, const double *y, size_t stride,
                             size_t cell, double k, double *gradient, size_t n)
{
    if (gradient) memset(gradient,0,n*sizeof(*gradient));
    if (c->kind == CHEM_LINEAR_CONSTRAINT) {
        double sum=0;
        for (size_t i=0; i<c->reactant_count; ++i) {
            ChemTerm t=c->reactants[i];
            sum+=t.coefficient*y[(size_t)t.species*stride+cell];
            if (gradient) gradient[t.species]+=t.coefficient;
        }
        return sum-c->constant;
    }
    double left=k, right=1;
    for (size_t i=0; i<c->reactant_count; ++i) {
        ChemTerm t=c->reactants[i]; double v=y[(size_t)t.species*stride+cell];
        left*=pow(v > 0 ? v : 0,t.coefficient);
    }
    for (size_t i=0; i<c->product_count; ++i) {
        ChemTerm t=c->products[i]; double v=y[(size_t)t.species*stride+cell];
        right*=pow(v > 0 ? v : 0,t.coefficient);
    }
    if (gradient) for (int side=0; side<2; ++side) {
        const ChemTerm *terms=side ? c->products : c->reactants;
        size_t count=side ? c->product_count : c->reactant_count;
        for (size_t i=0; i<count; ++i) {
            double partial=side ? 1 : k;
            for (size_t j=0; j<count; ++j) if (i!=j) {
                double v=y[(size_t)terms[j].species*stride+cell];
                partial*=pow(v > 0 ? v : 0,terms[j].coefficient);
            }
            double v=y[(size_t)terms[i].species*stride+cell];
            if (terms[i].coefficient != 1)
                partial = v > 0 ? terms[i].coefficient*partial*pow(v,terms[i].coefficient-1) : 0;
            gradient[terms[i].species] += side ? -partial : partial;
        }
    }
    return left-right;
}

int chem_constraint_evaluate(const ChemConstraint *c, const double *y, size_t n,
                             double temperature, double *residual, double *gradient)
{
    if (!y || !residual || validate_one(c,n)) return CHEM_INVALID;
    double k=c->kind == CHEM_EQUILIBRIUM_CONSTRAINT ? equilibrium_constant(c,temperature) : 0;
    *residual=chem_constraint_value(c,y,1,0,k,gradient,n);
    return isfinite(*residual) ? CHEM_OK : CHEM_NONFINITE;
}

int chem_update_constraint_constants(ChemState *s)
{
    for (size_t i=0; i<s->constraint_count; ++i)
        for (size_t cell=0; cell<s->cells; ++cell) {
            const ChemConstraint *c=&s->constraints[i];
            double k=c->kind == CHEM_EQUILIBRIUM_CONSTRAINT ? equilibrium_constant(c,s->temperatures[cell]) : 0;
            if (!isfinite(k)) return CHEM_NONFINITE;
            s->constraint_constants[i*s->cells+cell]=k;
        }
    return CHEM_OK;
}

void chem_constraint_residuals(const ChemState *s, const double *y, double *out)
{
    size_t n=s->mechanism->species_count, cells=s->cells;
    for (size_t i=0; i<s->constraint_count; ++i) {
        const ChemConstraint *c=&s->constraints[i];
        for (size_t cell=0; cell<cells; ++cell)
            out[(size_t)c->algebraic_species*cells+cell]=
                chem_constraint_value(c,y,cells,cell,s->constraint_constants[i*cells+cell],NULL,n);
    }
}

void chem_constraint_jacobian(ChemState *s, const double *y)
{
    size_t n=s->mechanism->species_count, cells=s->cells;
    for (size_t i=0; i<s->constraint_count; ++i) {
        const ChemConstraint *c=&s->constraints[i];
        size_t row=(size_t)c->algebraic_species;
        for (size_t cell=0; cell<cells; ++cell) {
            chem_constraint_value(c,y,cells,cell,s->constraint_constants[i*cells+cell],s->constraint_gradient,n);
            for (size_t col=0; col<n; ++col) if (s->index[row*n+col]>=0)
                s->matrix[(size_t)s->index[row*n+col]*cells+cell] -= s->constraint_gradient[col];
        }
    }
}
