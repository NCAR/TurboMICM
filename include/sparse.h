#ifndef SPARSE_H
#define SPARSE_H
#include "chem.h"

/* Immutable symbolic Mozart LU plan; safe to share across solver threads.
   The matrix passed to factor and solve uses dense row-major storage, and
   every possible nonzero must be represented by the mechanism Jacobian. */
typedef struct ChemSparsePlan ChemSparsePlan;
ChemSparsePlan *chem_sparse_plan_create(const ChemMechanism *mechanism);
void chem_sparse_plan_destroy(ChemSparsePlan *plan);
size_t chem_sparse_plan_size(const ChemSparsePlan *plan);
int chem_sparse_factor(const ChemSparsePlan *plan, double *matrix);
void chem_sparse_solve(const ChemSparsePlan *plan, const double *matrix, double *rhs);
#endif
