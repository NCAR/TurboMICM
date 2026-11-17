#include "chem.h"
#include "rates.h"
#include "sparse.h"
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
extern const double generated_default_temperature, generated_default_pressure;
int main(void)
{
    assert(generated_default_temperature==287.45 && generated_default_pressure==101319.9);
    const ChemConditions conditions = {generated_default_temperature,generated_default_pressure,generated_default_pressure/(1.380649e-23*6.02214076e23*generated_default_temperature)};
    double k[547], f[209], j[209*209], incremental[209*209]={0};
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
    ChemSparsePlan *plan=chem_sparse_plan_create(&generated_mechanism); assert(plan);
    for (size_t pass=0;pass<2;++pass) {
        chem_sparse_clear(plan,incremental);
        chem_accumulate_jacobian(&generated_mechanism,generated_initial_concentrations,k,incremental);
        for (size_t i=0;i<209*209;++i) assert(incremental[i]==j[i]);
        chem_sparse_negate(plan,incremental);
        for (size_t i=0;i<209*209;++i) assert(incremental[i]==-j[i]);
    }
    chem_sparse_plan_destroy(plan);
    for (size_t i=0; i<209; ++i) assert(isfinite(f[i]));
    for (size_t i=0; i<209*209; ++i) assert(isfinite(j[i]));
    puts("TS1 generated rates and forcing PASS");
    return 0;
}
