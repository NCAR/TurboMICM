#ifndef CHEM_INTERNAL_H
#define CHEM_INTERNAL_H
#include "chem.h"
/* Add a Jacobian to a matrix whose structural entries have been cleared. */
void chem_jacobian_add(const ChemMechanism *, ChemConditions, const double *, const double *, double *);
#endif
