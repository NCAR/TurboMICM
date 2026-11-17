#include "chapman_names.h"
static const char *const reactions[]={"O2_1","O3_1","O3_2","N2_O1D_1","O1D_O2_1","O_O3_1","M_O_O2_1"};
static const char *const species[]={"M","Ar","CO2","H2O","N2","O1D","O","O2","O3"};
const char *const *chem_chapman_reaction_names(size_t *count)
{ if (count) *count=sizeof(reactions)/sizeof(reactions[0]); return reactions; }
const char *const *chem_chapman_photolysis_names(size_t *count)
{ if (count) *count=3; return reactions; }
const char *const *chem_chapman_species_names(size_t *count)
{ if (count) *count=sizeof(species)/sizeof(species[0]); return species; }
