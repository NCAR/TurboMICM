#include "chem.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>

static const double pi=3.14159265358979323846;
static unsigned long long random_state=1;
static double uniform_random(void) { random_state=random_state*16807u%2147483647u; return (double)random_state/2147483647.0; }
static double lognormal_random(void)
{
    double a=uniform_random(),b=uniform_random();
    return exp(-2+2*sqrt(-2*log(a))*cos(2*pi*b));
}
static int terminator(size_t cells)
{
    const ChemSpecies species[]={{"Cl2",0,1e-20,0,0},{"Cl",0,1e-20,0,0}};
    const int reactants[]={0,1,1}; const ChemProduct products[]={{1,2},{0,1}};
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_USER,.parameter=0,.a=1}},
        {1,2,1,1,0,{.kind=CHEM_ARRHENIUS,.a=1,.d=300}}};
    const char *names[]={"toy_k1"}; const double defaults[]={0};
    const ChemMechanism m={2,2,1,species,reactions,reactants,products,names,defaults,298,101300};
    ChemState *s=chem_state_create(&m,cells); if (!s) return 0;
    double initial_cl2[4],initial_cl[4]; random_state=1;
    for (size_t c=0;c<cells;++c) { initial_cl2[c]=lognormal_random()*1e-6; initial_cl[c]=lognormal_random()*1e-10; }
    ChemRosenbrock p=chem_rosenbrock_three_stage(); p.relative_tolerance=1e-8; p.max_steps=100000;
    int ok=1;
    for (int step=0;step<(int)floor(2*pi/.3) && ok;++step) {
        double lon=step*.3;
        for (size_t c=0;c<cells;++c) {
            double lat=pi/180*(c*(90.0/cells));
            double rate=fmax(0,sin(lat)*sin(pi*20/180)+cos(lat)*cos(pi*20/180)*cos(lon-pi*300/180));
            chem_species_values(s,0)[c]=initial_cl2[c]; chem_species_values(s,1)[c]=initial_cl[c];
            chem_parameter_values(s,0)[c]=rate; chem_set_air_density(s,c,42);
        }
        ChemStats stats; ok=chem_update_rates(s)==CHEM_OK && chem_solve(s,&p,30,&stats)==CHEM_OK;
        for (size_t c=0;c<cells;++c) {
            double r=chem_parameter_values(s,0)[c]/4,cly=initial_cl[c]+2*initial_cl2[c];
            double det=sqrt(r*r+2*r*cly),e=exp(-4*det*30),l=det*30>1e-16 ? (1-e)/det/30 : 4;
            double df=-l*(initial_cl[c]-det+r)*(initial_cl[c]+det+r)/(1+e+30*l*(initial_cl[c]+r));
            double cl=initial_cl[c]+30*df,cl2=initial_cl2[c]-15*df;
            ok &= fabs(chem_species_values(s,1)[c]-cl)<=fabs(cl)*1e-6+1e-14 &&
                  fabs(chem_species_values(s,0)[c]-cl2)<=fabs(cl2)*1e-6+1e-14;
        }
    }
    chem_state_free(s); return ok;
}
static int terminator_all(void)
{ return terminator(1) && terminator(2) && terminator(3) && terminator(4); }
static int integrated_rates(void)
{
    const ChemSpecies species[]={{"a",1,1e-3,0,0},{"b",0,1e-3,0,0},{"c",0,1e-3,0,0},
                                 {"irr_1",0,1e-3,0,0},{"irr_2",0,1e-3,0,0}};
    const int reactants[]={0,1}; const ChemProduct products[]={{1,1},{3,1},{2,1},{4,1}};
    const ChemReaction reactions[]={
        {0,1,0,2,0,{.kind=CHEM_ARRHENIUS,.a=2.15e-11,.c=110,.d=300}},
        {1,1,2,2,0,{.kind=CHEM_USER,.parameter=0,.a=1}}};
    const char *names[]={"r2"}; const double defaults[]={1};
    const ChemMechanism m={5,2,1,species,reactions,reactants,products,names,defaults,273,1000};
    ChemState *s=chem_state_create(&m,1); if (!s) return 0;
    ChemRosenbrock p=chem_rosenbrock_three_stage(); int ok=1;
    for (int i=0;i<100 && ok;++i) {
        if (i>50) chem_parameter_values(s,0)[0]=0;
        double before1=chem_species_values(s,3)[0],before2=chem_species_values(s,4)[0]; ChemStats stats;
        ok=chem_update_rates(s)==CHEM_OK && chem_solve(s,&p,30,&stats)==CHEM_OK &&
           chem_species_values(s,3)[0]>=before1 && chem_species_values(s,4)[0]>=before2;
    }
    chem_state_free(s); return ok;
}
static int flow_tube(int backward)
{
    const ChemSpecies species[]={{"SOA1",0,1e-3,0,0},{"SOA2",0,1e-3,0,0},
                                 {"a-pinene",8e-8,1e-3,0,0},{"O3",2e-5,1e-3,0,0}};
    const int reactants[]={0,1,2,3}; const ChemProduct products[]={{0,.18},{1,.09}};
    const ChemReaction reactions[]={
        {0,1,0,0,0,{.kind=CHEM_USER,.parameter=0,.a=1}},
        {1,1,0,0,0,{.kind=CHEM_USER,.parameter=1,.a=1}},
        {2,2,0,2,0,{.kind=CHEM_ARRHENIUS,.a=8.8e-17*1e-6*6.02214076e23,.d=300}}};
    const char *names[]={"r1","r2"}; const double defaults[]={.01,.05};
    const ChemMechanism m={4,3,2,species,reactions,reactants,products,names,defaults,298.15,101325};
    FILE *file=fopen(backward ? "tests/fixtures/backward_euler_wall_loss.csv" : "tests/fixtures/rosenbrock_wall_loss.csv","r");
    if (!file) return 0;
    char line[256]; int ok=fgets(line,sizeof(line),file) && !strcmp(line,"time,SOA1,SOA2,a-pinene,O3\n");
    ChemState *s=chem_state_create(&m,1); if (!s) { fclose(file); return 0; }
    ok &= chem_update_rates(s)==CHEM_OK;
    ChemRosenbrock rp=chem_rosenbrock_two_stage(); ChemBackwardEuler bp=chem_backward_euler_defaults();
    for (int step=0;step<=3600 && ok;++step) {
        double time,expected[4];
        ok=fgets(line,sizeof(line),file) && sscanf(line,"%lf,%lf,%lf,%lf,%lf",&time,&expected[0],&expected[1],&expected[2],&expected[3])==5 && time==step;
        if (!ok) break;
        if (step<10) for (size_t i=0;i<4;++i) ok &= fabs(chem_species_values(s,i)[0]-expected[i])<=1e-10;
        if (step==3600) break;
        double elapsed=0;
        while (elapsed<1 && ok) {
            ChemStats stats; int status=backward ? chem_solve_backward_euler(s,&bp,1-elapsed,&stats) : chem_solve(s,&rp,1-elapsed,&stats);
            ok=(status==CHEM_OK || status==CHEM_TOO_MANY_STEPS) && stats.final_time>0;
            elapsed+=stats.final_time;
        }
    }
    ok &= !fgets(line,sizeof(line),file);
    chem_state_free(s); fclose(file); return ok;
}
static int parallel_states(void)
{
    const ChemSpecies species[]={{"A",0,1e-3,0,0},{"B",0,1e-3,0,0},{"C",0,1e-3,0,0}};
    const int reactants[]={0,1,1,1,2}; const ChemProduct products[]={{1,1},{1,1},{2,1},{0,1},{2,1}};
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_USER,.parameter=0,.a=1}},
        {1,2,1,2,0,{.kind=CHEM_USER,.parameter=1,.a=1}},
        {3,2,3,2,0,{.kind=CHEM_USER,.parameter=2,.a=1}}};
    const char *names[]={"r1","r2","r3"}; const double defaults[]={.04,3e7,1e4};
    const ChemMechanism m={3,3,3,species,reactions,reactants,products,names,defaults,272.5,101253.3};
    ChemState *base=chem_state_create(&m,1); if (!base) return 0;
    double results[8][3]={{0}}; int good[8]={0};
    #pragma omp parallel for num_threads(8) schedule(static)
    for (int thread=0;thread<8;++thread) {
        ChemState *s=chem_state_create_like(base,1); int ok=s!=NULL;
        if (s) {
            chem_species_values(s,0)[0]=1; chem_set_air_density(s,0,1e6);
            ChemRosenbrock p=chem_rosenbrock_three_stage(); ok &= chem_update_rates(s)==CHEM_OK;
            for (int i=0;i<10 && ok;++i) {
                double elapsed=0;
                while (elapsed<200 && ok) {
                    ChemStats stats; int status=chem_solve(s,&p,200-elapsed,&stats);
                    ok=(status==CHEM_OK || status==CHEM_TOO_MANY_STEPS) && stats.final_time>0;
                    elapsed+=stats.final_time;
                }
            }
            for (size_t i=0;i<3;++i) results[thread][i]=chem_species_values(s,i)[0];
            chem_state_free(s);
        }
        good[thread]=ok;
    }
    int ok=1;
    for (int thread=0;thread<8;++thread) {
        ok &= good[thread];
        for (size_t i=0;i<3;++i) ok &= results[thread][i]==results[0][i];
    }
    chem_state_free(base); return ok;
}
static int failures;
static void report(const char *id,int ok)
{ printf("%s %s %s\n",id,ok ? "PASS" : "FAIL",ok ? "reference integration matched" : "integration differed"); failures+=!ok; }
int main(void)
{
    report("terminator:RosenbrockSolver.Terminator",terminator_all());
    report("terminator:RosenbrockSolver.VectorTerminator",terminator_all());
    report("integrated_reaction_rates:ChapmanIntegration.CanBuildChapmanSystem",integrated_rates());
    report("regression:RosenbrockRegression.FlowTube",flow_tube(0));
    report("regression:BackwardEulerRegression.FlowTube",flow_tube(1));
    report("openmp:OpenMP.OneSolverManyStates",parallel_states());
    return failures ? 1 : 0;
}
