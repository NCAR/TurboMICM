#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "backward_euler.h"
#include "process.h"

typedef struct { const ChemMechanism *mechanism; const double *rates; } Context;
static void rhs(void *v,const double *y,double *out)
{ Context *c=v; chem_rhs(c->mechanism,y,c->rates,out); }
static void jac(void *v,const double *y,double *out)
{ Context *c=v; chem_jacobian(c->mechanism,y,c->rates,out); }
static void run(const ChemBuiltMechanism *built,int backward)
{
    const char *path=backward ? "../micm/test/regression/expected_results/backward_euler_wall_loss.csv" :
                               "../micm/test/regression/expected_results/rosenbrock_wall_loss.csv";
    FILE *file=fopen(path,"r"); assert(file);
    char line[512]; assert(fgets(line,sizeof(line),file));
    assert(!strcmp(line,"time,SOA1,SOA2,a-pinene,O3\n"));
    double custom[]={.01,.05},rates[3],y[]={0,0,8e-8,2e-5};
    ChemConditions conditions={298.15,101325,101325/(1.380649e-23*6.02214076e23*298.15)};
    chem_rates_parameterized(3,built->rates,built->multipliers,built->functions,conditions,custom,rates);
    Context context={&built->mechanism,rates};
    ChemRosenbrockMethod method=chem_rosenbrock_two_stage();
    ChemRosenbrockOptions ro=chem_rosenbrock_defaults(built->absolute_tolerances);
    ChemBackwardEulerOptions be=chem_backward_euler_defaults(built->absolute_tolerances);
    for (size_t step=0;step<=3600;++step) {
        double time,expected[4];
        assert(fgets(line,sizeof(line),file));
        assert(sscanf(line,"%lf,%lf,%lf,%lf,%lf",&time,expected,expected+1,expected+2,expected+3)==5);
        assert(time==(double)step);
        if (step<10) for (size_t j=0;j<4;++j) assert(fabs(y[j]-expected[j])<=1e-10);
        if (step==3600) break;
        double elapsed=0;
        while (elapsed<1) {
            ChemSolverStats stats;
            ChemSolveStatus status=backward ? chem_backward_euler_solve(4,y,1-elapsed,rhs,jac,&context,&be,&stats) :
                chem_rosenbrock_solve(4,y,1-elapsed,rhs,jac,&context,&method,&ro,&stats);
            assert((status==CHEM_SOLVE_OK || status==CHEM_SOLVE_MAX_STEPS) && stats.final_time>0);
            elapsed+=stats.final_time;
        }
        for (size_t j=0;j<4;++j) assert(isfinite(y[j]) && y[j]>=0);
        assert(fabs((y[3]-y[2])-(2e-5-8e-8))<1e-12);
    }
    assert(!fgets(line,sizeof(line),file)); fclose(file);
}
int main(void)
{
    const char *names[]={"SOA1","SOA2","a-pinene","O3"};
    ChemSpecies species[4]={{0}}; ChemPhaseSpecies phase_species[4]={{0}};
    for (size_t i=0;i<4;++i) {
        assert(chem_species_init(species+i,names[i])==CHEM_SPECIES_OK);
        assert(chem_phase_species_init(phase_species+i,species+i)==CHEM_SPECIES_OK);
    }
    ChemPhase phase={0}; ChemSystem system={0};
    assert(chem_phase_init(&phase,"gas",4,phase_species)==CHEM_SPECIES_OK);
    assert(chem_system_init(&system,&phase)==CHEM_SPECIES_OK);
    const ChemSpecies *inputs0[]={species},*inputs1[]={species+1},*inputs2[]={species+2,species+3};
    const ChemSpecies *outputs[]={species,species+1}; const double yields[]={.18,.09};
    ChemProcessConfig configs[]={
        {1,0,inputs0,NULL,NULL,&phase,{CHEM_RATE_USER,{1},0}},
        {1,0,inputs1,NULL,NULL,&phase,{CHEM_RATE_USER,{1},1}},
        {2,2,inputs2,outputs,yields,&phase,{CHEM_RATE_ARRHENIUS,{8.8e-17*6.02214076e17,0,0,300,0},0}}};
    ChemProcess processes[3]={{0}};
    for (size_t i=0;i<3;++i) assert(chem_process_init(processes+i,configs+i)==CHEM_SPECIES_OK);
    ChemBuiltMechanism built={0}; assert(chem_build_mechanism(&built,&system,3,processes,0)==CHEM_BUILD_OK);
    assert(built.mechanism.species_count==4 && built.mechanism.reaction_count==3);
    for (size_t i=0;i<4;++i) assert(!strcmp(built.state_names[i],names[i]));
    run(&built,0); run(&built,1);
    chem_built_mechanism_destroy(&built);
    for (size_t i=0;i<3;++i) chem_process_destroy(processes+i);
    chem_system_destroy(&system); chem_phase_destroy(&phase);
    for (size_t i=0;i<4;++i) { chem_phase_species_destroy(phase_species+i); chem_species_destroy(species+i); }
    return 0;
}
