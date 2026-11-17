#include "chem.h"
#include "rates.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
extern const ChemMechanism generated_mechanism;
extern const ChemRate generated_rates[];
extern const char *const generated_species[];
int main(void)
{
    assert(generated_mechanism.reaction_count == 3);
    ChemConditions c = {250,101253.3,30};
    double actual[3]; chem_rates(3,generated_rates,c,NULL,actual);
#ifdef REFERENCE_ORDER
    const size_t reversible=2, jpl=0;
#else
    const size_t reversible=0, jpl=2;
#endif
    assert(fabs(actual[reversible]-exp(.4)) < 1e-13);
    assert(fabs(actual[1]-3*(1+.01*250+.001*250*250)) < 1e-12);
    double k0=.003*pow(250.0/298,-2), ki=1.5*pow(250.0/298,-1), ratio=k0*30/ki;
    double expected=k0*30/(1+ratio)*pow(.6,1/(1+pow(log10(ratio),2)));
    assert(fabs(actual[jpl]-expected) < 1e-14);
    size_t a=0,b=0;
    for (size_t i=0;i<generated_mechanism.species_count;++i) {
        if (!strcmp(generated_species[i],"A")) a=i;
        if (!strcmp(generated_species[i],"B")) b=i;
    }
    assert(a!=b);
    double y[2]={0},forcing[2]; y[a]=1;
    chem_rhs(&generated_mechanism,y,actual,forcing);
    assert(fabs(forcing[a]+actual[reversible]+actual[1]+actual[jpl])<1e-11);
    assert(fabs(forcing[b]-actual[reversible]-2*actual[1]-3*actual[jpl])<1e-11);
    puts("generated rate extensions PASS");
    return 0;
}
