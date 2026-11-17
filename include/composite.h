#ifndef CHEM_COMPOSITE_H
#define CHEM_COMPOSITE_H
#include "rates.h"
/* Each component adds to the forcing and positive dense row-major Jacobian.
   Component contexts and arrays remain owned by the caller. */
typedef struct {
    void *context;
    void (*update)(void *context,ChemConditions conditions,double *parameters);
    void (*add_forcing)(void *context,const double *state,double *forcing);
    void (*add_jacobian)(void *context,const double *state,double *jacobian);
} ChemOdeComponent;
typedef struct { size_t variables,count; const ChemOdeComponent *components; } ChemCompositeOde;
void chem_composite_update(const ChemCompositeOde *system,ChemConditions conditions,double *parameters);
/* Solver-compatible callbacks; the composite clears the output once before
   calling its components. Both outputs must not alias the input state. */
void chem_composite_rhs(void *system,const double *state,double *forcing);
void chem_composite_jacobian(void *system,const double *state,double *jacobian);
#endif
