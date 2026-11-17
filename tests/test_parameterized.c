#include "chem.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

extern const ChemMechanism generated_mechanism;
int main(void)
{
    const ChemMechanism *m = &generated_mechanism;
    if (m->species_count != 2 || m->reaction_count != 2 || m->parameterized_species_count != 2) return 1;
    size_t a = !strcmp(m->species_names[0], "A") ? 0 : 1, b = 1-a;
    if (strcmp(m->species_names[a], "A") || strcmp(m->species_names[b], "B")) return 1;
    ChemConditions c = {300,100000,5,1};
    double rates[2], y[2], f[2], j[4];
    y[a]=2; y[b]=3;
    if (chem_rate_constants(m,c,NULL,rates) || fabs(rates[0]-10)>1e-14 || fabs(rates[1]-5.5)>1e-14) return 2;
    chem_forcing(m,c,rates,y,f); chem_jacobian(m,c,rates,y,j);
    if (fabs(f[a]+3.5)>1e-14 || fabs(f[b]-3.5)>1e-14 ||
        fabs(j[a*2+a]+10)>1e-14 || fabs(j[a*2+b]-5.5)>1e-14 ||
        fabs(j[b*2+a]-10)>1e-14 || fabs(j[b*2+b]+5.5)>1e-14) return 3;
    puts("generated parameterized reactants: PASS");
    return 0;
}
