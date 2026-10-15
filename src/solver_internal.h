#ifndef CHEM_SOLVER_INTERNAL_H
#define CHEM_SOLVER_INTERNAL_H
#include "chem_solver.h"
#include "internal.h"
typedef struct {
    size_t reaction, independent, dependent_count, product_count;
    size_t *dependent, *updates;
    double *yields;
} Derivative;
struct ChemSolver {
    const ChemMechanism *mechanism;
    size_t cells, count, derivative_count;
    ChemSparse *sparse;
    Derivative *derivatives;
    double *matrix, *next, *initial, *error, *stages, *work;
    ChemStepController controller;
    void *controller_context;
    ChemLinearSolver linear;
};
void chem_solver_negative_jacobian(ChemSolver *, const ChemState *, const double *);
#endif
