/* SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_CONSTRAINTS_H
#define CHEM_CONSTRAINTS_H
#include "chem.h"

typedef enum { CHEM_LINEAR_CONSTRAINT, CHEM_EQUILIBRIUM_CONSTRAINT } ChemConstraintType;
typedef struct { const char *species; double coefficient; } ChemConstraintTerm;
typedef struct {
    double reference_constant, enthalpy, gas_constant, reference_temperature;
} ChemVantHoff;
typedef struct {
    ChemConstraintType type;
    const char *name, *algebraic_species;
    /* Linear: reactants are the weighted terms; products is empty. */
    size_t reactant_count, product_count;
    const ChemConstraintTerm *reactants, *products;
    double constant;
    ChemVantHoff equilibrium;
} ChemConstraint;
ChemVantHoff chem_vant_hoff(double reference_constant, double enthalpy);
ChemStatus chem_constraint_validate(const ChemConstraint *);
/* Copies all definitions. Call before creating states or solvers. Replacement
 * is transactional: an invalid definition leaves the old constraints intact. */
ChemStatus chem_mechanism_set_constraints(ChemMechanism *, size_t count, const ChemConstraint *);
size_t chem_constraint_count(const ChemMechanism *);
const ChemConstraint *chem_constraint_definition(const ChemMechanism *, size_t index);
size_t chem_constraint_row(const ChemMechanism *, size_t index);
int chem_is_algebraic(const ChemMechanism *, size_t variable);
int chem_constraint_has_element(const ChemMechanism *, size_t row, size_t column);
void chem_constraints_update_parameters(ChemState *);
/* These routines affect only constraint rows. Residuals overwrite their rows;
 * Jacobian contributions are subtracted from existing dense matrix entries. */
void chem_constraints_forcing(const ChemState *, const double *values, double *forcing);
void chem_constraints_subtract_jacobian(const ChemState *, const double *values, double *jacobian);
#endif
