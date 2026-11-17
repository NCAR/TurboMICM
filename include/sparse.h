#ifndef SPARSE_H
#define SPARSE_H
#include "chem.h"

/* Immutable symbolic Mozart LU plan; safe to share across solver threads.
   The matrix passed to factor and solve uses dense row-major storage, and
   every possible nonzero must be represented by the mechanism Jacobian. */
typedef struct ChemSparsePlan ChemSparsePlan;
/* The pattern is an n-by-n row-major array of zero/nonzero bytes. */
ChemSparsePlan *chem_sparse_plan_from_pattern(size_t n, const unsigned char *pattern);
ChemSparsePlan *chem_sparse_plan_create(const ChemMechanism *mechanism);
void chem_sparse_plan_destroy(ChemSparsePlan *plan);
size_t chem_sparse_plan_size(const ChemSparsePlan *plan);
/* Touch only the plan's original and fill positions in a dense matrix. */
void chem_sparse_clear(const ChemSparsePlan *plan, double *matrix);
void chem_sparse_negate(const ChemSparsePlan *plan, double *matrix);
int chem_sparse_factor(const ChemSparsePlan *plan, double *matrix);
int chem_sparse_factor_copy(const ChemSparsePlan *plan, const double *matrix, double *factors);
int chem_sparse_factor_doolittle(const ChemSparsePlan *plan, double *matrix);
int chem_sparse_factor_doolittle_copy(const ChemSparsePlan *plan, const double *matrix, double *factors);
void chem_sparse_solve(const ChemSparsePlan *plan, const double *matrix, double *rhs);
#endif
