#ifndef CHEM_INTERNAL_H
#define CHEM_INTERNAL_H
#include "chem.h"
static inline double chem_shift_diagonal(double negative_jacobian, double alpha, int differential)
{
    return negative_jacobian + (differential ? alpha : 0);
}
/* Add a Jacobian to a matrix whose structural entries have been cleared. */
void chem_jacobian_add(const ChemMechanism *, ChemConditions, const double *, const double *, double *);
#endif
