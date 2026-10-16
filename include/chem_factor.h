#ifndef CHEM_FACTOR_H
#define CHEM_FACTOR_H
#include "chem_sparse.h"
/* Layout adapters for the same symbolic plan and arithmetic used by ODE solves.
 * A plan can be reused with any number of blocks or cell-group width. */
chem_status chem_sparse_lu_create(chem_lu_plan *plan, const chem_sparse *matrix, chem_lu_algorithm algorithm);
chem_status chem_sparse_lu_matrix(const chem_lu_plan *plan, const chem_sparse *prototype,
                                  chem_sparse *combined, double initial);
chem_status chem_sparse_lu_matrices(const chem_lu_plan *plan, const chem_sparse *prototype,
                                    chem_sparse *lower, chem_sparse *upper, double initial);
chem_status chem_sparse_lu_factor_inplace(const chem_lu_plan *plan, chem_sparse *combined);
chem_status chem_sparse_lu_factor(const chem_lu_plan *plan, const chem_sparse *matrix,
                                  chem_sparse *lower, chem_sparse *upper);
chem_status chem_sparse_lu_solve_inplace(const chem_lu_plan *plan, const chem_sparse *combined, chem_array *rhs);
chem_status chem_sparse_lu_solve(const chem_lu_plan *plan, const chem_sparse *lower,
                                 const chem_sparse *upper, chem_array *rhs);
#endif
