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
extern const ChemSpeciesCatalog generated_species_catalog;
extern const ChemRate generated_rates[];
extern const ChemRateMultiplier generated_multipliers[];
extern const ChemParameterizedSpecies generated_multiplier_functions[];
extern const char *const generated_species[], *const generated_custom_labels[];
extern const size_t generated_custom_count;
int main(void)
{
    assert(generated_mechanism.species_count==1 && generated_mechanism.reaction_count==5);
    assert(chem_species_catalog_used_count(&generated_species_catalog)==4);
    assert(chem_species_catalog_is_used(&generated_species_catalog,"foo"));
    assert(chem_species_catalog_is_used(&generated_species_catalog,"bar"));
    assert(chem_species_catalog_is_used(&generated_species_catalog,"probe"));
    assert(chem_species_catalog_is_used(&generated_species_catalog,"Q"));
    assert(!strcmp(generated_species[0],"foo") && generated_custom_count==3);
    assert(!strcmp(generated_custom_labels[0],"bar_user"));
    assert(!strcmp(generated_custom_labels[1],"foo_surf.effective radius [m]"));
    assert(!strcmp(generated_custom_labels[2],"foo_surf.particle number concentration [# m-3]"));
    assert(generated_multipliers[0].count==1 && generated_multipliers[1].count==1 &&
           generated_multipliers[2].count==0 && generated_multipliers[3].count==1 && generated_multipliers[4].count==2);
    for (size_t cell=0;cell<5;++cell) {
        ChemConditions c={285+cell,101100+20*cell,10+.3*cell};
        double custom[]={.01+.001*cell,1e-8,1e5}, rates[5];
        chem_rates_parameterized(5,generated_rates,generated_multipliers,generated_multiplier_functions,c,custom,rates);
        double expected=12.2*exp(300/c.temperature)*(.82*c.air_density);
        assert(fabs(rates[0]-expected)<=1e-12*expected);
        expected=custom[0]*(.82*c.air_density);
        assert(fabs(rates[1]-expected)<=1e-12*expected);
        const double pi=3.14159265358979323846, gas=1.380649e-23*6.02214076e23;
        double speed=sqrt(8*gas/(pi*.025)*c.temperature);
        expected=4*custom[2]*pi*custom[1]*custom[1]/(custom[1]/230+4/speed);
        assert(fabs(rates[2]-expected)<=1e-10*expected);
        assert(rates[3]==2*c.air_density);
        double q=1+.01*c.temperature+.0001*c.pressure+.5*c.air_density;
        assert(fabs(rates[4]-q*q)<1e-12);
        double y[]={.5}, f[1], j[1];
        chem_rhs(&generated_mechanism,y,rates,f); chem_jacobian(&generated_mechanism,y,rates,j);
        expected=-rates[0]*.5-rates[2]*.5+rates[3]-rates[4]*.5;
        assert(fabs(f[0]-expected)<1e-12*fmax(1,fabs(expected)));
        expected=-rates[0]-rates[2]-rates[4];
        assert(fabs(j[0]-expected)<1e-12*fabs(expected));
    }
    ChemConditions probe={300,101325,5}; double custom[]={.01,1e-8,1e5}, k[5];
    chem_rates_parameterized(5,generated_rates,generated_multipliers,generated_multiplier_functions,probe,custom,k);
    assert(k[3]==10);
    const ChemRate simple_rates[]={{CHEM_RATE_ARRHENIUS,{1,0,0,300,0},0},{CHEM_RATE_ARRHENIUS,{1,0,0,300,0},0}};
    const ChemRateMultiplier simple_multipliers[]={{0,0},{0,1}};
    const ChemParameterizedSpecies density_double[]={{0,0,0,2}};
    double simple_out[2];
    chem_rates_parameterized(2,simple_rates,simple_multipliers,density_double,probe,NULL,simple_out);
    assert(simple_multipliers[0].count==0 && simple_multipliers[1].count==1 && simple_out[1]==10);
    puts("parameterized species and process rates PASS");
    return 0;
}
