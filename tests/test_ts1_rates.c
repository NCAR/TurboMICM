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
extern const unsigned generated_air_orders[];
extern const double generated_custom_defaults[];
extern const double generated_initial_concentrations[];
int main(void)
{
    const ChemConditions conditions = {287.45,101319.9,101319.9/(1.380649e-23*6.02214076e23*287.45)};
    double k[547], f[209], j[209*209];
    assert(generated_mechanism.species_count == 209 && generated_mechanism.reaction_count == 547);
    chem_rates_with_air(547, generated_rates, generated_air_orders, conditions, generated_custom_defaults, k);
    unsigned third_bodies = 0;
    for (size_t i=0; i<547; ++i) {
        assert(isfinite(k[i]) && k[i] >= 0);
        third_bodies += generated_air_orders[i];
    }
    assert(third_bodies == 4);
    chem_rhs(&generated_mechanism, generated_initial_concentrations, k, f);
    chem_jacobian(&generated_mechanism, generated_initial_concentrations, k, j);
    for (size_t i=0; i<209; ++i) assert(isfinite(f[i]));
    for (size_t i=0; i<209*209; ++i) assert(isfinite(j[i]));
    puts("TS1 generated rates and forcing PASS");
    return 0;
}
