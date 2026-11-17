#include "chem.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>
extern const ChemMechanism generated_mechanism;
extern const char *const generated_species[];
int main(void)
{
    const double y[] = {1, 2, 3}, k[] = {0.04, 3e7, 1e4};
    double f[3];
    assert(generated_mechanism.species_count == 3);
    assert(generated_species[0][0] == 'A');
    chem_rhs(&generated_mechanism, y, k, f);
    assert(f[0] == -0.04 + 60000);
    assert(fabs(f[1] - (0.04 - 120000000 - 60000)) < 1e-7);
    assert(f[2] == 120000000);
    puts("generated mechanism PASS");
    return 0;
}
