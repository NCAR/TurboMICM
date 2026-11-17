#undef NDEBUG
#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "state.h"

static double equilibrium(double k,double enthalpy,double temperature)
{ return k*exp(enthalpy/(1.380649e-23*6.02214076e23)*(1/temperature-1/298.15)); }
static void run(size_t count)
{
    const char *const names[]={"A","B","C","D","E","F"};
    const double atol[]={.001,.001,.001,.001,.001,.001};
    const ChemConstraintTerm r1[]={{1,1}},p1[]={{2,1}},r2[]={{4,1}},p2[]={{5,1}};
    const double gas=1.380649e-23*6.02214076e23;
    ChemConstraint rows[]={
        {.name="B_C_eq",.type=CHEM_CONSTRAINT_EQUILIBRIUM,.row=2,.reactants=r1,.products=p1,
         .reactant_count=1,.product_count=1,.constant=.034,.enthalpy=-2400,.gas_constant=gas,.reference_temperature=298.15},
        {.name="E_F_eq",.type=CHEM_CONSTRAINT_EQUILIBRIUM,.row=5,.reactants=r2,.products=p2,
         .reactant_count=1,.product_count=1,.constant=.012,.enthalpy=-2000,.gas_constant=gas,.reference_temperature=298.15}};
    size_t n=count==1 ? 3 : 6; ChemConstraintSet set={n,count,rows};
    ChemStateConfig config={1,n,0,0,names,NULL,1e-6,atol};
    ChemState state={0};
    assert(chem_state_init_constrained(&state,&config,0,NULL,&set)==CHEM_STATE_OK);
    assert(state.species_count==n && state.constraint_count==count && state.parameter_count==count);
    for (size_t i=0;i<n;++i) {
        assert(chem_state_species_index(&state,names[i])==i);
        assert(state.algebraic[i]==(i==2 || (count==2 && i==5)));
    }
    assert(chem_state_parameter_index(&state,"B_C_eq")==0);
    if (count==2) assert(chem_state_parameter_index(&state,"E_F_eq")==1);
    state.concentrations[0]=1; state.concentrations[1]=.1; state.concentrations[2]=.034*.1;
    if (count==2) { state.concentrations[3]=.5; state.concentrations[4]=.05; state.concentrations[5]=.012*.05; }
    double temperature=count==1 ? 300 : 310;
    state.conditions[0]=(ChemConditions){temperature,101325,0};
    assert(chem_state_update_constraint_parameters(&state,&set)==CHEM_STATE_OK);
    double k1=equilibrium(.034,-2400,temperature);
    assert(fabs(state.parameters[0]-k1)<=fabs(k1)*100*DBL_EPSILON);
    if (count==2) {
        double k2=equilibrium(.012,-2000,temperature);
        assert(fabs(state.parameters[0]-k1)<1e-10);
        assert(fabs(state.parameters[1]-k2)<=fabs(k2)*100*DBL_EPSILON);
    }
    ChemState copy={0}; assert(chem_state_clone(&copy,&state)==CHEM_STATE_OK);
    assert(copy.constraint_count==count && copy.algebraic!=state.algebraic && !memcmp(copy.algebraic,state.algebraic,n));
    assert(copy.parameters[0]==state.parameters[0]);
    double old=state.parameters[0]; state.conditions[0].temperature=0;
    assert(chem_state_update_constraint_parameters(&state,&set)==CHEM_STATE_INVALID && state.parameters[0]==old);
    chem_state_destroy(&copy); chem_state_destroy(&state);
    config.cells=3;
    assert(chem_state_init_constrained(&state,&config,0,NULL,&set)==CHEM_STATE_OK);
    for (size_t c=0;c<3;++c) state.conditions[c].temperature=290+10*c;
    assert(chem_state_update_constraint_parameters(&state,&set)==CHEM_STATE_OK);
    for (size_t c=0;c<3;++c) {
        double expected=equilibrium(.034,-2400,290+10*c);
        assert(fabs(state.parameters[c*count]-expected)<=fabs(expected)*100*DBL_EPSILON);
    }
    old=state.parameters[0]; state.conditions[2].temperature=-1;
    assert(chem_state_update_constraint_parameters(&state,&set)==CHEM_STATE_INVALID && state.parameters[0]==old);
    chem_state_destroy(&state);
}
int main(void) { run(1); run(2); return 0; }
