#ifndef CHEM_INTERNAL_H
#define CHEM_INTERNAL_H

#include "chem.h"
#include <stdatomic.h>

typedef struct { int target, lower; } ChemUpdate;
typedef struct { int position; size_t first_update, updates; } ChemUpper;
typedef struct {
    int diagonal;
    size_t first_lower, lowers, first_upper, uppers;
} ChemPivot;

struct ChemState {
    const ChemMechanism *mechanism;
    const ChemConstraint *constraints;
    size_t constraint_count;
    size_t cells, nonzeros;
    double *values, *parameters, *temperatures, *pressures, *air_densities, *absolute_tolerances, *rates;
    unsigned char *explicit_air_density;
    ChemRateFunction *rate_functions;
    void **rate_contexts;
    double *matrix, *scratch, *forcing, *candidate, *error, *stages;
    double *constraint_constants, *constraint_gradient;
    unsigned char *algebraic;
    int *index, *columns;
    size_t *row_start;
    ChemPivot *pivots;
    int *lower;
    ChemUpper *upper;
    ChemUpdate *updates;
    atomic_size_t *symbolic_references;
};

int chem_make_symbolic(ChemState *state);
void chem_forcing(const ChemState *state, const double *values, double *out, double *scratch);
void chem_negative_jacobian(ChemState *state, const double *values);
void chem_constraint_residuals(const ChemState *state, const double *values, double *out);
void chem_constraint_jacobian(ChemState *state, const double *values);
int chem_update_constraint_constants(ChemState *state);
int chem_project_constraints(ChemState *state, const ChemRosenbrock *parameters, ChemStats *stats);
double chem_constraint_value(const ChemConstraint *constraint, const double *values, size_t stride,
                             size_t cell, double equilibrium_constant, double *gradient, size_t species_count);
int chem_factor(ChemState *state, double alpha);
void chem_substitute(const ChemState *state, double *values);

#endif
