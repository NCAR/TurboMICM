#ifndef CHEM_INTERNAL_H
#define CHEM_INTERNAL_H
#include "chem.h"
#include <stdint.h>

typedef struct { size_t destination, left, right; } ChemElimination;
typedef struct { size_t index, column; } ChemEntry;
typedef struct { size_t reaction, independent, skipped; } ChemDerivative;
struct ChemMechanism {
    const ChemMechanismData *data;
    size_t nnz, *index, *row_start;
    ChemEntry *entries;
    size_t *diagonal;
    size_t *lower_start, *lower;
    size_t *operation_start;
    ChemElimination *operations;
    size_t derivative_count;
    ChemDerivative *derivatives;
};
struct ChemState {
    const ChemMechanism *mechanism;
    size_t cells;
    double *y, *custom, *temperature, *pressure, *rates, *atol;
    double rtol;
    double *jacobian, *forcing, *candidate, *error, *stages, *scratch;
};
ChemStatus chem_factor(const ChemMechanism *, size_t cells, double *matrix);
void chem_linear_solve(const ChemMechanism *, size_t cells, const double *matrix, double *x);
#endif
