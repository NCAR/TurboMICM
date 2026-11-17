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
    size_t cells, nonzeros;
    double *values, *parameters, *temperatures, *pressures, *rates;
    double *matrix, *scratch, *forcing, *candidate, *error, *stages;
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
int chem_factor(ChemState *state, double alpha);
void chem_substitute(const ChemState *state, double *values);

#endif
