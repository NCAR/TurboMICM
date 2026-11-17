#include "composite.h"
#include <string.h>
void chem_composite_update(const ChemCompositeOde *system,ChemConditions conditions,double *parameters)
{
    for (size_t i=0;i<system->count;++i) {
        const ChemOdeComponent *component=system->components+i;
        if (component->update) component->update(component->context,conditions,parameters);
    }
}
void chem_composite_rhs(void *v,const double *state,double *forcing)
{
    const ChemCompositeOde *system=v;
    memset(forcing,0,system->variables*sizeof(*forcing));
    for (size_t i=0;i<system->count;++i) {
        const ChemOdeComponent *component=system->components+i;
        if (component->add_forcing) component->add_forcing(component->context,state,forcing);
    }
}
void chem_composite_jacobian(void *v,const double *state,double *jacobian)
{
    const ChemCompositeOde *system=v;
    memset(jacobian,0,system->variables*system->variables*sizeof(*jacobian));
    for (size_t i=0;i<system->count;++i) {
        const ChemOdeComponent *component=system->components+i;
        if (component->add_jacobian) component->add_jacobian(component->context,state,jacobian);
    }
}
