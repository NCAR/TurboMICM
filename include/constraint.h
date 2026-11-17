#ifndef CONSTRAINT_H
#define CONSTRAINT_H
#include <stddef.h>

typedef struct { size_t species; double coefficient; } ChemConstraintTerm;
typedef enum { CHEM_CONSTRAINT_LINEAR, CHEM_CONSTRAINT_EQUILIBRIUM } ChemConstraintType;
typedef struct {
    const char *name;
    ChemConstraintType type;
    size_t row;
    const ChemConstraintTerm *reactants, *products;
    size_t reactant_count, product_count;
    double constant; /* linear total or equilibrium reference constant */
    double enthalpy, gas_constant, reference_temperature;
} ChemConstraint;
typedef struct { size_t variables, count; const ChemConstraint *constraints; } ChemConstraintSet;
typedef enum {
    CHEM_CONSTRAINT_OK, CHEM_CONSTRAINT_INVALID, CHEM_CONSTRAINT_EMPTY_REACTANTS,
    CHEM_CONSTRAINT_EMPTY_PRODUCTS, CHEM_CONSTRAINT_BAD_STOICHIOMETRY,
    CHEM_CONSTRAINT_BAD_EQUILIBRIUM, CHEM_CONSTRAINT_UNKNOWN_SPECIES,
    CHEM_CONSTRAINT_DUPLICATE_ROW
} ChemConstraintStatus;

/* Structures and their arrays remain owned by the caller. */
ChemConstraintStatus chem_constraint_validate(const ChemConstraintSet *set);
/* State and forcing are cells-by-variables, Jacobian is cells-by-n-by-n.
   parameters, when non-NULL, is cells-by-constraint and overrides K(T).
   Otherwise equilibrium constraints require temperatures, one per cell. */
void chem_constraint_residuals(const ChemConstraintSet *set, size_t cells, const double *state,
                               const double *temperatures, const double *parameters, double *forcing);
void chem_constraint_subtract_jacobian(const ChemConstraintSet *set, size_t cells, const double *state,
                                       const double *temperatures, const double *parameters, double *jacobian);
/* OR constraint positions, including each algebraic diagonal, into a row-major pattern. */
void chem_constraint_sparsity(const ChemConstraintSet *set, unsigned char *pattern);
void chem_constraint_step_errors(const ChemConstraintSet *set, size_t cells,
                                 const double *old_state, const double *new_state, double *errors);
#endif
