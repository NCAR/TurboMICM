#include "chem.h"
#include "rates.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
extern const ChemMechanism generated_mechanism;
extern const ChemRate generated_rates[];
int main(void)
{
    assert(generated_mechanism.reaction_count == 3);
    ChemConditions c = {250,101253.3,30};
    double actual[3]; chem_rates(3,generated_rates,c,NULL,actual);
    assert(fabs(actual[0]-exp(.4)) < 1e-13);
    assert(fabs(actual[1]-3*(1+.01*250+.001*250*250)) < 1e-12);
    double k0=.003*pow(250.0/298,-2), ki=1.5*pow(250.0/298,-1), ratio=k0*30/ki;
    double expected=k0*30/(1+ratio)*pow(.6,1/(1+pow(log10(ratio),2)));
    assert(fabs(actual[2]-expected) < 1e-14);
    puts("generated rate extensions PASS");
    return 0;
}
