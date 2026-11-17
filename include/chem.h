#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>

/* Terms are ordered; repeated reactants have separate entries.  A reactant
 * index of -1 denotes the third body, parameterized by air density. */
typedef struct { int species; double yield; } ChemTerm;
typedef enum {
    CHEM_ARRHENIUS, CHEM_TROE, CHEM_TERNARY, CHEM_TUNNELING,
    CHEM_USER, CHEM_SURFACE, CHEM_BRANCH_ALKOXY, CHEM_BRANCH_NITRATE
} ChemRateType;
typedef struct {
    ChemRateType type;
    size_t reactant_begin, reactant_count, product_begin, product_count;
    int parameter;                 /* first custom parameter, or -1 */
    double p[10];                  /* interpreted by the rate type */
} ChemReaction;
typedef struct {
    size_t species_count, reaction_count, parameter_count;
    const char *const *species_names, *const *parameter_names;
    const double *initial_species, *initial_parameters;
    const ChemReaction *reactions;
    const ChemTerm *reactants, *products;
} ChemMechanism;
typedef struct { double temperature, pressure; } ChemConditions;

/* Return nonzero for invalid conditions, arguments, or non-finite results. */
int chem_rate_constants(const ChemMechanism *, ChemConditions, const double *parameters, double *rates);
void chem_forcing(const ChemMechanism *, ChemConditions, const double *rates, const double *y, double *f);
void chem_jacobian(const ChemMechanism *, ChemConditions, const double *rates, const double *y, double *j);

#endif
