#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <string.h>
#include "process.h"
#include "rosenbrock.h"
#include "state.h"

static void init_species(ChemSpecies *s,const char *name) { assert(chem_species_init(s,name)==CHEM_SPECIES_OK); }
typedef struct { const ChemBuiltMechanism *built; const double *rates; } Ode;
static void rhs(void *v,const double *y,double *out)
{ Ode *s=v; chem_rhs(&s->built->mechanism,y,s->rates,out); }
static void jac(void *v,const double *y,double *out)
{ Ode *s=v; chem_jacobian(&s->built->mechanism,y,s->rates,out); }
static void init_phase(ChemPhase *phase,const char *name,ChemSpecies *species,size_t count)
{
    ChemPhaseSpecies entries[5]={{0}}; assert(count<=5);
    for (size_t i=0;i<count;++i) assert(chem_phase_species_init(entries+i,species+i)==CHEM_SPECIES_OK);
    assert(chem_phase_init(phase,name,count,entries)==CHEM_SPECIES_OK);
    for (size_t i=0;i<count;++i) chem_phase_species_destroy(entries+i);
}
int main(void)
{
    ChemSpecies oxides[4]={{0}}; const char *names[]={"O3","NO","NO2","O2"};
    for (size_t i=0;i<4;++i) init_species(oxides+i,names[i]);
    ChemPhase gas={0}; init_phase(&gas,"gas",oxides,4);
    const ChemSpecies *reactants[]={oxides,oxides+1},*products[]={oxides+2,oxides+3};
    const double yields[]={1,1};
    ChemProcessConfig config={2,2,reactants,products,yields,&gas,{CHEM_RATE_ARRHENIUS,{.5,0,0,300,0},0}};
    ChemProcess reaction={0},copy={0};
    assert(chem_process_init(&reaction,&config)==CHEM_SPECIES_OK);
    assert(reaction.reactant_count==2 && !strcmp(reaction.products[0].species.name,"NO2") && !strcmp(reaction.phase.name,"gas"));
    assert(chem_process_clone(&copy,&reaction)==CHEM_SPECIES_OK);
    assert(copy.reactants!=reaction.reactants && !strcmp(copy.reactants[0].name,reaction.reactants[0].name));
    ChemSystem system={0}; assert(chem_system_init(&system,&gas)==CHEM_SPECIES_OK);
    ChemBuiltMechanism built={0};
    assert(chem_build_mechanism(&built,NULL,1,&reaction,0)==CHEM_BUILD_MISSING_SYSTEM);
    assert(chem_build_mechanism(&built,&system,1,&reaction,0)==CHEM_BUILD_OK);
    double k[1],y[]={2,3,0,0},f[4];
    chem_rates_parameterized(1,built.rates,built.multipliers,built.functions,(ChemConditions){300,0,0},NULL,k);
    chem_rhs(&built.mechanism,y,k,f);
    assert(k[0]==.5 && f[0]==-3 && f[1]==-3 && f[2]==3 && f[3]==3);
    {
        const double atol[]={1e-10,1e-10,1e-10,1e-10};
        ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol); options.relative_tolerance=1e-9;
        ChemRosenbrockMethod method=chem_rosenbrock_three_stage(); Ode system_ode={&built,k};
        assert(chem_rosenbrock_solve(4,y,1,rhs,jac,&system_ode,&method,&options,NULL)==CHEM_SOLVE_OK);
        double a=2/(3*exp(.5)-2);
        assert(fabs(y[0]-a)<1e-7 && fabs(y[1]-(a+1))<1e-7);
        assert(fabs(y[2]-(2-a))<1e-7 && fabs(y[3]-(2-a))<1e-7);
    }
    chem_built_mechanism_destroy(&built); chem_process_destroy(&copy); chem_process_destroy(&reaction);
    chem_system_destroy(&system); chem_phase_destroy(&gas);
    for (size_t i=0;i<4;++i) chem_species_destroy(oxides+i);

    ChemSpecies foo={0},bar={0}; init_species(&foo,"foo"); init_species(&bar,"bar");
    assert(chem_species_set_property(&foo,"molecular weight [kg mol-1]",CHEM_PROPERTY_DOUBLE,(ChemPropertyValue){.real=.025})==CHEM_SPECIES_OK);
    assert(chem_species_set_property(&foo,"diffusion coefficient [m2 s-1]",CHEM_PROPERTY_DOUBLE,(ChemPropertyValue){.real=230})==CHEM_SPECIES_OK);
    chem_species_set_parameterization(&bar,(ChemParameterizedSpecies){0,0,0,.82});
    ChemPhaseSpecies phase_species[2]={{0}};
    assert(chem_phase_species_init(phase_species,&foo)==CHEM_SPECIES_OK && chem_phase_species_init(phase_species+1,&bar)==CHEM_SPECIES_OK);
    assert(chem_phase_init(&gas,"gas",2,phase_species)==CHEM_SPECIES_OK);
    const ChemSpecies *r[]={&foo,&bar};
    config=(ChemProcessConfig){2,0,r,NULL,NULL,&gas,{CHEM_RATE_ARRHENIUS,{12.2,0,300,300,0},0}};
    assert(chem_process_init(&reaction,&config)==CHEM_SPECIES_OK && chem_process_clone(&copy,&reaction)==CHEM_SPECIES_OK);
    assert(!strcmp(copy.reactants[0].name,reaction.reactants[0].name) && copy.product_count==0);
    assert(!strcmp(copy.phase.name,reaction.phase.name) && copy.rate.type==CHEM_RATE_ARRHENIUS && copy.rate.p[0]==12.2);
    assert(chem_system_init(&system,&gas)==CHEM_SPECIES_OK);
    assert(chem_build_mechanism(&built,&system,1,&reaction,0)==CHEM_BUILD_OK && built.mechanism.species_count==1);
    chem_rates_parameterized(1,built.rates,built.multipliers,built.functions,(ChemConditions){300,0,5},NULL,k);
    assert(fabs(k[0]-(12.2*exp(1)*.82*5))<1e-12);
    chem_built_mechanism_destroy(&built); chem_system_destroy(&system); chem_phase_destroy(&gas);
    chem_process_destroy(&copy); chem_process_destroy(&reaction);
    for (size_t i=0;i<2;++i) chem_phase_species_destroy(phase_species+i);
    chem_species_destroy(&foo); chem_species_destroy(&bar);

    ChemSpecies abc[4]={{0}}; const char *abc_names[]={"A","B","C","D"};
    for (size_t i=0;i<4;++i) init_species(abc+i,abc_names[i]);
    init_phase(&gas,"gas",abc,4); assert(chem_system_init(&system,&gas)==CHEM_SPECIES_OK);
    ChemProcess chain[2]={{0}}; const ChemSpecies *from[]={abc,abc+1},*to[]={abc+1,abc+2};
    double unit[]={1};
    for (size_t i=0;i<2;++i) {
        ChemProcessConfig c={1,1,from+i,to+i,unit,&gas,{CHEM_RATE_ARRHENIUS,{1,0,0,300,0},0}};
        assert(chem_process_init(chain+i,&c)==CHEM_SPECIES_OK);
    }
    assert(chem_build_mechanism(&built,&system,2,chain,0)==CHEM_BUILD_UNUSED_SPECIES);
    assert(chem_build_mechanism(&built,&system,2,chain,1)==CHEM_BUILD_OK);
    assert(built.mechanism.species_count==4 && !strcmp(built.system.gas.name,"gas"));
    chem_built_mechanism_destroy(&built); chem_system_destroy(&system); chem_phase_destroy(&gas);
    for (size_t i=0;i<2;++i) chem_process_destroy(chain+i);
    for (size_t i=0;i<4;++i) chem_species_destroy(abc+i);

    ChemSpecies tolerances[2]={{0}}; init_species(tolerances,"foo"); init_species(tolerances+1,"bar");
    assert(chem_species_set_property(tolerances,"absolute tolerance",CHEM_PROPERTY_DOUBLE,(ChemPropertyValue){.real=1e-7})==CHEM_SPECIES_OK);
    assert(chem_species_set_property(tolerances+1,"absolute tolerance",CHEM_PROPERTY_DOUBLE,(ChemPropertyValue){.real=1e-8})==CHEM_SPECIES_OK);
    init_phase(&gas,"gas",tolerances,2); assert(chem_system_init(&system,&gas)==CHEM_SPECIES_OK);
    const ChemSpecies *one_r[]={tolerances},*one_p[]={tolerances+1};
    config=(ChemProcessConfig){1,1,one_r,one_p,unit,&gas,{CHEM_RATE_ARRHENIUS,{2e-11,0,110,300,0},0}};
    assert(chem_process_init(&reaction,&config)==CHEM_SPECIES_OK);
    assert(chem_build_mechanism(&built,&system,1,&reaction,0)==CHEM_BUILD_OK);
    for (size_t cells=1;cells<=10;++cells) {
        ChemState state={0};
        ChemStateConfig state_config={cells,2,0,1,built.state_names,NULL,1e-6,built.absolute_tolerances};
        assert(chem_state_init(&state,&state_config)==CHEM_STATE_OK && state.species_count==2);
        assert(state.absolute_tolerances[0]==1e-7 && state.absolute_tolerances[1]==1e-8);
        chem_state_destroy(&state);
    }
    chem_built_mechanism_destroy(&built); chem_system_destroy(&system); chem_phase_destroy(&gas); chem_process_destroy(&reaction);
    chem_species_destroy(tolerances); chem_species_destroy(tolerances+1);
    return 0;
}
