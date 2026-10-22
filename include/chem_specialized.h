#ifndef CHEM_SPECIALIZED_H
#define CHEM_SPECIALIZED_H
#include "chem.h"

/* Generated code uses the default Markowitz/Mozart model for this descriptor.
 * Conditions, custom parameters, concentrations, and tolerances remain mutable.
 * Rates must be refreshed after changing conditions or custom parameters.
 * solve uses MICM's default three-stage Rosenbrock and does not clamp its result;
 * like the generic low-level solver, callers clamp at external step boundaries.
 */
typedef struct {
    const chem_mechanism *mechanism;
    void (*rates)(chem_state *);
    chem_kernel forcing, jacobian;
    void (*factor)(size_t cells, double *matrix, double *scratch);
    void (*linear)(size_t cells, const double *matrix, double *rhs);
    chem_result (*solve)(chem_state *, double dt);
} chem_specialized;
#endif
