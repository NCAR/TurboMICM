#ifndef CHEM_CONSTRAINT_MATRIX_H
#define CHEM_CONSTRAINT_MATRIX_H

#include "chem.h"
#include "matrix.h"
#include "sparse.h"

/* Apply constraints to dense states and forcing, with a sparse negative
 * Jacobian. Standard and grouped matrix layouts may be mixed; all matrices
 * must have one row/block per temperature and the sparse pattern must contain
 * every constraint dependency. Only algebraic rows are overwritten. */
ChemConstraintStatus chem_constraint_set_apply_matrices(ChemConstraintSet set,
                                                        const double *temperatures, size_t cells,
                                                        const ChemMatrix *state, ChemMatrix *forcing,
                                                        ChemSparseMatrix *negative_jacobian);

#endif
