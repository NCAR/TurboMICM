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

/* Rosenbrock parameters are data, so the step controller and formula can be
 * changed independently. Arrays use row-wise strictly lower triangular order. */
typedef struct {
    int stages, new_function[6], max_steps;
    double a[15], c[15], m[6], e[6], gamma, order;
    double factor_min, factor_max, rejection_factor, safety;
    double h_min, h_max, h_start, relative_tolerance, absolute_tolerance;
} ChemRosenbrockParameters;
typedef struct {
    unsigned long internal_steps, accepted, rejected, function_calls;
    unsigned long jacobian_updates, decompositions, solves;
} ChemStatistics;
typedef struct ChemWorkspace ChemWorkspace;

ChemRosenbrockParameters chem_rosenbrock_three_stage(void);
ChemWorkspace *chem_workspace_create(const ChemMechanism *);
void chem_workspace_free(ChemWorkspace *);
/* On failure returns nonzero; y contains the most recently accepted state. */
int chem_rosenbrock_solve(const ChemMechanism *, ChemWorkspace *,
                          const ChemRosenbrockParameters *, ChemConditions,
                          const double *rates, double dt, double *y, ChemStatistics *);

#endif
