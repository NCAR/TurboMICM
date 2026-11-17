#include "constraint.h"
#include <math.h>

ChemConstraintStatus chem_constraint_validate(const ChemConstraintSet *set)
{
    if (!set || !set->variables || (set->count && !set->constraints)) return CHEM_CONSTRAINT_INVALID;
    for (size_t i=0;i<set->count;++i) {
        const ChemConstraint *c=set->constraints+i;
        if (c->row >= set->variables) return CHEM_CONSTRAINT_UNKNOWN_SPECIES;
        for (size_t j=0;j<i;++j) if (set->constraints[j].row == c->row) return CHEM_CONSTRAINT_DUPLICATE_ROW;
        if (c->type != CHEM_CONSTRAINT_LINEAR && c->type != CHEM_CONSTRAINT_EQUILIBRIUM) return CHEM_CONSTRAINT_INVALID;
        if (!c->reactant_count) return CHEM_CONSTRAINT_EMPTY_REACTANTS;
        if (!c->reactants || (c->product_count && !c->products)) return CHEM_CONSTRAINT_INVALID;
        if (c->type == CHEM_CONSTRAINT_EQUILIBRIUM && !c->product_count) return CHEM_CONSTRAINT_EMPTY_PRODUCTS;
        for (size_t side=0;side<2;++side) {
            const ChemConstraintTerm *terms=side ? c->products : c->reactants;
            size_t count=side ? c->product_count : c->reactant_count;
            for (size_t j=0;j<count;++j) {
                if (terms[j].species >= set->variables) return CHEM_CONSTRAINT_UNKNOWN_SPECIES;
                if (!isfinite(terms[j].coefficient) ||
                    (c->type == CHEM_CONSTRAINT_EQUILIBRIUM && terms[j].coefficient <= 0))
                    return CHEM_CONSTRAINT_BAD_STOICHIOMETRY;
            }
        }
        if (!isfinite(c->constant) || (c->type == CHEM_CONSTRAINT_EQUILIBRIUM && c->constant <= 0))
            return CHEM_CONSTRAINT_BAD_EQUILIBRIUM;
        if (c->type == CHEM_CONSTRAINT_EQUILIBRIUM &&
            (!isfinite(c->enthalpy) || !isfinite(c->gas_constant) || c->gas_constant <= 0 ||
             !isfinite(c->reference_temperature) || c->reference_temperature <= 0)) return CHEM_CONSTRAINT_INVALID;
    }
    return CHEM_CONSTRAINT_OK;
}
static double keq(const ChemConstraintSet *set, const ChemConstraint *c, size_t cell, size_t index,
                  const double *temperature, const double *parameters)
{
    if (parameters) return parameters[cell*set->count+index];
    return c->constant*exp(c->enthalpy/c->gas_constant*(1/temperature[cell]-1/c->reference_temperature));
}
static double product(const ChemConstraintTerm *terms, size_t count, const double *state, size_t excluded)
{
    double value=1;
    for (size_t i=0;i<count;++i) if (i != excluded)
        value *= pow(state[terms[i].species] > 0 ? state[terms[i].species] : 0,terms[i].coefficient);
    return value;
}
static double derivative(const ChemConstraintTerm *terms, size_t count, const double *state, size_t i)
{
    double value=product(terms,count,state,i), coefficient=terms[i].coefficient;
    double concentration=state[terms[i].species];
    if (coefficient == 1) return value;
    return concentration > 0 ? coefficient*value*pow(concentration,coefficient-1) : 0;
}
void chem_constraint_residuals(const ChemConstraintSet *set, size_t cells, const double *state,
                               const double *temperatures, const double *parameters, double *forcing)
{
    size_t n=set->variables;
    for (size_t cell=0;cell<cells;++cell) {
        const double *y=state+cell*n;
        for (size_t i=0;i<set->count;++i) {
            const ChemConstraint *c=set->constraints+i;
            double value;
            if (c->type == CHEM_CONSTRAINT_LINEAR) {
                value=0;
                for (size_t j=0;j<c->reactant_count;++j) value += c->reactants[j].coefficient*y[c->reactants[j].species];
                value -= c->constant;
            } else {
                value=keq(set,c,cell,i,temperatures,parameters)*product(c->reactants,c->reactant_count,y,(size_t)-1)
                     -product(c->products,c->product_count,y,(size_t)-1);
            }
            forcing[cell*n+c->row]=value;
        }
    }
}
void chem_constraint_subtract_jacobian(const ChemConstraintSet *set, size_t cells, const double *state,
                                       const double *temperatures, const double *parameters, double *jacobian)
{
    size_t n=set->variables;
    for (size_t cell=0;cell<cells;++cell) {
        const double *y=state+cell*n;
        for (size_t i=0;i<set->count;++i) {
            const ChemConstraint *c=set->constraints+i;
            double *row=jacobian+cell*n*n+c->row*n;
            if (c->type == CHEM_CONSTRAINT_LINEAR) {
                for (size_t j=0;j<c->reactant_count;++j) row[c->reactants[j].species]-=c->reactants[j].coefficient;
            } else {
                double k=keq(set,c,cell,i,temperatures,parameters);
                for (size_t j=0;j<c->reactant_count;++j) row[c->reactants[j].species]-=k*derivative(c->reactants,c->reactant_count,y,j);
                for (size_t j=0;j<c->product_count;++j) row[c->products[j].species]+=derivative(c->products,c->product_count,y,j);
            }
        }
    }
}
void chem_constraint_sparsity(const ChemConstraintSet *set, unsigned char *pattern)
{
    size_t n=set->variables;
    for (size_t i=0;i<set->count;++i) {
        const ChemConstraint *c=set->constraints+i;
        pattern[c->row*n+c->row]=1;
        for (size_t j=0;j<c->reactant_count;++j) pattern[c->row*n+c->reactants[j].species]=1;
        for (size_t j=0;j<c->product_count;++j) pattern[c->row*n+c->products[j].species]=1;
    }
}
void chem_constraint_step_errors(const ChemConstraintSet *set, size_t cells,
                                 const double *old_state, const double *new_state, double *errors)
{
    for (size_t cell=0;cell<cells;++cell)
        for (size_t i=0;i<set->count;++i) {
            size_t k=cell*set->variables+set->constraints[i].row;
            errors[k]=new_state[k]-old_state[k];
        }
}
