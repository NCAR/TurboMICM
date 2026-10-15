#ifndef CHEM_INTERNAL_H
#define CHEM_INTERNAL_H
#include "chem.h"
#include <stdint.h>

typedef struct { size_t destination, left, right; } ChemElimination;
typedef struct { size_t index, column; } ChemEntry;
typedef struct { size_t reaction, independent, skipped; } ChemDerivative;
struct ChemMechanism {
    const ChemMechanismData *data;
    ChemRate *rate_specs;
    size_t nnz, *index, *row_start;
    ChemEntry *entries;
    size_t *diagonal;
    size_t *lower_start, *lower;
    size_t *operation_start;
    ChemElimination *operations;
    size_t derivative_count;
    ChemDerivative *derivatives;
    unsigned char *algebraic, *constraint_pattern;
};
struct ChemState {
    const ChemMechanism *mechanism;
    size_t cells;
    double *y, *custom, *temperature, *pressure, *air_density, *rates, *atol;
    double rtol;
    ChemLinearSolver linear;
    double *jacobian, *forcing, *candidate, *error, *stages, *scratch;
};
void chem_post_solve_clamp(ChemState *);
ChemStatus chem_factor(const ChemMechanism *, size_t cells, double *matrix);
void chem_linear_solve(const ChemMechanism *, size_t cells, const double *matrix, double *x);
#endif
