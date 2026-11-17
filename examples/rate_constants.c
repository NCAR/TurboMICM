#include "chem.h"
#include "rates.h"
#include "rosenbrock.h"
#include <stdio.h>

/* Standalone, hand-built counterpart of the rate-constant tutorials. */
typedef struct { const ChemMechanism *mechanism; double rates[10]; } System;
static void rhs(void *v,const double *y,double *f)
{ System *s=v; chem_rhs(s->mechanism,y,s->rates,f); }
static void jacobian(void *v,const double *y,double *j)
{ System *s=v; chem_jacobian(s->mechanism,y,s->rates,j); }
static void print_state(double time,const double *y)
{
    printf("%6g,",time);
    for (size_t i=0;i<7;++i) printf(i ? ",%11.2e" : "%11.2e",y[i]);
    putchar('\n');
}
int main(void)
{
    const double conversion=6.02214076e17;
    const size_t reactants[]={0,1,1,2,3,4,4,5,2,1};
    const ChemProduct products[]={{1,1},{2,1},{3,1},{4,1},{5,2},{6,1},{6,1},{6,1},{0,1}};
    const ChemReaction reactions[]={
        {0,1,0,1},{1,1,1,1},{2,1,2,1},{3,1,3,1},{4,1,4,1},
        {5,2,5,1},{7,1,6,1},{8,1,7,1},{9,1,8,0},{10,0,8,1}};
    const ChemRate rates[]={
        {CHEM_RATE_ARRHENIUS,{2.15e-4,0,110,300,0},0},
        {CHEM_RATE_ALKOXY,{1.2,204.3,1e-3,2},0},
        {CHEM_RATE_NITRATE,{1.2,204.3,1e-3,2},0},
        {CHEM_RATE_SURFACE,{2.3e2,.025,.9},0},
        {CHEM_RATE_TERNARY,{1.2,2.3,302.3,2.6/conversion,-3.1,402.1,.9,1.2},0},
        {CHEM_RATE_TROE,{1.8e-30*conversion*conversion,-3,0,2.8e-11*conversion,0,0,.6,1},0},
        {CHEM_RATE_TUNNELING,{1.2,2.3,302.3},0},
        {CHEM_RATE_USER,{1},2},{CHEM_RATE_USER,{1},3},{CHEM_RATE_USER,{1},4}};
#ifdef WITH_CUSTOM
    const size_t count=10;
#else
    const size_t count=7;
#endif
    const ChemMechanism mechanism={7,count,reactions,reactants,products};
    const ChemConditions conditions={287.45,101319.9,101319.9/(1.380649e-23*6.02214076e23*287.45)};
    const double atol[]={1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3};
    double custom[]={1e-7,2.5e6,1e-10,1e-23,1e-20}, y[]={1,0,0,0,0,0,0};
    System system={&mechanism,{0}};
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
    ChemRosenbrockMethod method=chem_rosenbrock_three_stage();
    printf("%6s","time");
    for (int i=0;i<7;++i) printf(",%11c",'A'+i);
    putchar('\n');
    print_state(0,y);
    for (size_t step=0;step<10;++step) {
        chem_rates(count,rates,conditions,custom,system.rates);
        double elapsed=0;
        while (elapsed<500) {
            ChemSolverStats stats;
            ChemSolveStatus status=chem_rosenbrock_solve(7,y,500-elapsed,rhs,jacobian,&system,&method,&options,&stats);
            if ((status!=CHEM_SOLVE_OK && status!=CHEM_SOLVE_MAX_STEPS) || stats.final_time<=0) {
                fprintf(stderr,"Solver made no progress: %s\n",chem_solve_status_string(status));
                return 1;
            }
            elapsed+=stats.final_time;
        }
        print_state(500*(step+1),y);
        custom[2]*=1.5;
    }
    return 0;
}
