#include "chem.h"
#include <math.h>
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
        if (!(temperature > 0) || !(c->reference_temperature > 0) || !(c->gas_constant > 0)) return CHEM_CONSTRAINT_BAD_EQUILIBRIUM;
        double k = c->constant * exp((c->enthalpy/c->gas_constant)*(1/temperature-1/c->reference_temperature));
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
