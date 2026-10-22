#ifndef CHEM_SPECIALIZED_H
#define CHEM_SPECIALIZED_H
#include "chem.h"
extern const ChemMechanismData ts1_data;
/* Generated functions require a state created from the matching descriptor.
   Temperature, pressure, air density, custom parameters and tolerances are data.
   Fixed rate descriptors must be changed in the configuration and regenerated. */
ChemStatus ts1_update_rates(ChemState *);
ChemStatus ts1_solve(ChemState *, const ChemRosenbrock *, double, ChemStats *);
#endif
