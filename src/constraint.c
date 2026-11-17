#include "chem.h"
#include <math.h>
#include <stdint.h>
#include <string.h>

ChemConstraintStatus chem_constraint_validate(const ChemConstraint *c, size_t n)
{
    if (!c || c->row >= n) return CHEM_CONSTRAINT_BAD_ROW;
    if (!c->reactant_count || !c->reactants) return CHEM_CONSTRAINT_EMPTY_REACTANTS;
    if (c->type == CHEM_EQUILIBRIUM_CONSTRAINT && (!c->product_count || !c->products)) return CHEM_CONSTRAINT_EMPTY_PRODUCTS;
    if (c->type != CHEM_LINEAR_CONSTRAINT && c->type != CHEM_EQUILIBRIUM_CONSTRAINT) return CHEM_CONSTRAINT_BAD_ROW;
    for (size_t side = 0; side < 2; ++side) {
        const ChemTerm *terms = side ? c->products : c->reactants;
        size_t count = side ? c->product_count : c->reactant_count;
        if (count && !terms) return CHEM_CONSTRAINT_BAD_SPECIES;
        for (size_t i = 0; i < count; ++i) {
            if (terms[i].species < 0 || (size_t)terms[i].species >= n) return CHEM_CONSTRAINT_BAD_SPECIES;
            if (c->type == CHEM_EQUILIBRIUM_CONSTRAINT && !(terms[i].yield > 0)) return CHEM_CONSTRAINT_BAD_STOICHIOMETRY;
        }
    }
    if (c->type == CHEM_EQUILIBRIUM_CONSTRAINT && !(c->constant > 0)) return CHEM_CONSTRAINT_BAD_EQUILIBRIUM;
    return CHEM_CONSTRAINT_OK;
}

static double product(const ChemTerm *terms, size_t count, const double *y, size_t omit)
{
    double value = 1;
    for (size_t i = 0; i < count; ++i) if (i != omit) {
        double concentration = y[terms[i].species];
        value *= pow(concentration > 0 ? concentration : 0, terms[i].yield);
    }
    return value;
}

ChemConstraintStatus chem_constraint_equilibrium_constant(const ChemConstraint *c,double temperature,double *value)
{
    if(!c || !value || c->type!=CHEM_EQUILIBRIUM_CONSTRAINT || !(c->constant>0) || !(temperature>0) ||
       !(c->reference_temperature>0) || !(c->gas_constant>0)) return CHEM_CONSTRAINT_BAD_EQUILIBRIUM;
    *value=c->constant*exp((c->enthalpy/c->gas_constant)*(1/temperature-1/c->reference_temperature));
    return isfinite(*value)?CHEM_CONSTRAINT_OK:CHEM_CONSTRAINT_BAD_EQUILIBRIUM;
}

ChemConstraintStatus chem_constraint_apply(const ChemConstraint *c, size_t n, double temperature,
                                          const double *y, double *forcing, double *negative_jacobian)
{
    ChemConstraintStatus status = chem_constraint_validate(c, n);
    if (status != CHEM_CONSTRAINT_OK) return status;
    if (!y || !forcing || !negative_jacobian) return CHEM_CONSTRAINT_BAD_ROW;
    double *row = negative_jacobian + c->row*n;
    memset(row, 0, n*sizeof(*row));
    if (c->type == CHEM_LINEAR_CONSTRAINT) {
        double residual = 0;
        for (size_t i = 0; i < c->reactant_count; ++i) {
            ChemTerm term = c->reactants[i];
            residual += term.yield*y[term.species];
            row[term.species] -= term.yield;
        }
        forcing[c->row] = residual-c->constant;
    } else {
        double k;
        if (chem_constraint_equilibrium_constant(c,temperature,&k)) return CHEM_CONSTRAINT_BAD_EQUILIBRIUM;
        forcing[c->row] = k*product(c->reactants,c->reactant_count,y,(size_t)-1)
                        - product(c->products,c->product_count,y,(size_t)-1);
        for (size_t side = 0; side < 2; ++side) {
            const ChemTerm *terms = side ? c->products : c->reactants;
            size_t count = side ? c->product_count : c->reactant_count;
            for (size_t i = 0; i < count; ++i) {
                double exponent = terms[i].yield, concentration = y[terms[i].species];
                double partial = product(terms,count,y,i)*(side ? 1 : k);
                if (exponent != 1) partial = concentration > 0 ? exponent*partial*pow(concentration, exponent-1) : 0;
                row[terms[i].species] += (side ? partial : -partial);
            }
        }
    }
    return CHEM_CONSTRAINT_OK;
}

ChemConstraintStatus chem_constraint_set_validate(ChemConstraintSet set, size_t n)
{
    if (set.count && !set.items) return CHEM_CONSTRAINT_BAD_ROW;
    for (size_t i=0; i<set.count; ++i) {
        ChemConstraintStatus status = chem_constraint_validate(&set.items[i],n);
        if (status) return status;
        for (size_t j=0; j<i; ++j)
            if (set.items[i].row == set.items[j].row) return CHEM_CONSTRAINT_DUPLICATE_ROW;
    }
    return CHEM_CONSTRAINT_OK;
}

ChemConstraintStatus chem_constraint_set_pattern(ChemConstraintSet set, size_t n, unsigned char *pattern)
{
    ChemConstraintStatus status = chem_constraint_set_validate(set,n);
    if (status) return status;
    if ((!pattern && n) || (n && n > SIZE_MAX/n)) return CHEM_CONSTRAINT_BAD_ROW;
    if (n) memset(pattern,0,n*n);
    for (size_t i=0; i<set.count; ++i) {
        const ChemConstraint *c=&set.items[i];
        pattern[c->row*n+c->row]=1;
        for (size_t j=0;j<c->reactant_count;++j) pattern[c->row*n+(size_t)c->reactants[j].species]=1;
        for (size_t j=0;j<c->product_count;++j) pattern[c->row*n+(size_t)c->products[j].species]=1;
    }
    return CHEM_CONSTRAINT_OK;
}

ChemConstraintStatus chem_constraint_mass_diagonal(ChemConstraintSet set,size_t n,unsigned char *diagonal)
{
    ChemConstraintStatus status=chem_constraint_set_validate(set,n);
    if(status) return status;
    if(n && !diagonal) return CHEM_CONSTRAINT_BAD_ROW;
    if(n) memset(diagonal,1,n);
    for(size_t i=0;i<set.count;++i) diagonal[set.items[i].row]=0;
    return CHEM_CONSTRAINT_OK;
}

ChemConstraintStatus chem_constraint_set_apply(ChemConstraintSet set, size_t n, size_t cells,
                                               const double *temperature, const double *y,
                                               double *forcing, double *jacobians)
{
    ChemConstraintStatus status=chem_constraint_set_validate(set,n);
    if (status) return status;
    if (set.count && cells && (!temperature || !y || !forcing || !jacobians)) return CHEM_CONSTRAINT_BAD_ROW;
    if (n && (n > SIZE_MAX/n || cells > SIZE_MAX/(n*n))) return CHEM_CONSTRAINT_BAD_ROW;
    for (size_t cell=0; cell<cells; ++cell) for (size_t i=0; i<set.count; ++i) {
        status=chem_constraint_apply(&set.items[i],n,temperature[cell],y+cell*n,forcing+cell*n,jacobians+cell*n*n);
        if (status) return status;
    }
    return CHEM_CONSTRAINT_OK;
}
