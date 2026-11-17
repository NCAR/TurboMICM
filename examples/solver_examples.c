#include "chem.h"
#include "rates.h"
#include "rosenbrock.h"
#include <stdio.h>

#ifndef EXAMPLE
#define EXAMPLE 0
#endif
typedef struct { const ChemMechanism *mechanism; double rates[3]; } System;
static void rhs(void *context,const double *y,double *out)
{ System *s=context; chem_rhs(s->mechanism,y,s->rates,out); }
static void jacobian(void *context,const double *y,double *out)
{ System *s=context; chem_jacobian(s->mechanism,y,s->rates,out); }
static void header(size_t cells,const char *const *names)
{
    printf("%6s","time");
    if (cells>1) printf(",%6s","grid");
    for (size_t i=0;i<3;++i) printf(",%11s",names[i]);
    putchar('\n');
}
static void print_state(double time,size_t cells,const double *y)
{
    for (size_t c=0;c<cells;++c) {
        printf("%6g,",time);
        if (cells>1) printf("%6zu,",c);
        for (size_t i=0;i<3;++i) printf(i ? ",%11.2e" : "%11.2e",y[c*3+i]);
        putchar('\n');
    }
}
int main(void)
{
    const size_t readme_reactants[]={0,0,1};
    const ChemProduct readme_products[]={{1,.8},{2,.2},{2,1}};
    const ChemReaction readme_reactions[]={{0,1,0,2},{1,2,2,1}};
    const ChemMechanism readme={3,2,readme_reactions,readme_reactants,readme_products};
    const size_t robertson_reactants[]={0,1,1,1,2};
    const ChemProduct robertson_products[]={{1,1},{1,1},{2,1},{0,1},{2,1}};
    const ChemReaction robertson_reactions[]={{0,1,0,1},{1,2,1,2},{3,2,3,2}};
    const ChemMechanism robertson={3,3,robertson_reactions,robertson_reactants,robertson_products};
    const char *const readme_names[]={"Foo","Bar","Baz"};
    const char *const robertson_names[]={"A","B","C"};
    const size_t cells=EXAMPLE==0 ? 1 : 3;
    System systems[3]; void *contexts[3];
    for (size_t i=0;i<cells;++i) {
        systems[i]=(System){EXAMPLE==0 ? &readme : &robertson,{0}};
        contexts[i]=&systems[i];
        if (EXAMPLE==0) {
            ChemConditions conditions={287.45,101319.9,0};
            ChemRate rates[]={{CHEM_RATE_ARRHENIUS,{.001,0,0,300,0},0},
                              {CHEM_RATE_ARRHENIUS,{1e-5,0,110,300,0},0}};
            chem_rates(2,rates,conditions,NULL,systems[i].rates);
        } else {
            systems[i].rates[0]=.04; systems[i].rates[1]=3e7; systems[i].rates[2]=1e4;
        }
    }
    double y[9]={0};
    if (EXAMPLE==0) y[0]=20; else { y[0]=1; y[3]=2; y[6]=.5; }
    const double atol[]={.001,.001,.001};
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
    ChemRosenbrockMethod method=chem_rosenbrock_three_stage();
    if (EXAMPLE!=2) {
        header(cells,EXAMPLE==0 ? readme_names : robertson_names);
        if (EXAMPLE==1) print_state(0,cells,y);
    }
    for (size_t step=0;step<(EXAMPLE==2 ? 1u : 10u);++step) {
        double duration=EXAMPLE==0 ? 500 : 200,elapsed=0;
        ChemSolverStats stats;
        do {
            ChemSolveStatus status=chem_rosenbrock_solve_grouped(cells,3,y,duration-elapsed,rhs,jacobian,
                                     contexts,NULL,&method,&options,&stats);
            if ((status!=CHEM_SOLVE_OK && status!=CHEM_SOLVE_MAX_STEPS) || stats.final_time<=0) {
                fprintf(stderr,"Solver made no progress: %s\n",chem_solve_status_string(status));
                return 1;
            }
            elapsed+=stats.final_time;
        } while (EXAMPLE!=2 && elapsed<duration);
        if (EXAMPLE==0) print_state(step*500,cells,y);
        else if (EXAMPLE==1) print_state((step+1)*200,cells,y);
        else {
            printf("Solver state: Converged\naccepted: %zu\nfunction_calls: %zu\njacobian_updates: %zu\n"
                   "number_of_steps: %zu\naccepted: %zu\nrejected: %zu\ndecompositions: %zu\nsolves: %zu\n"
                   "final simulation time: %g\n",stats.accepted,stats.function_calls,stats.jacobian_updates,
                   stats.steps,stats.accepted,stats.rejected,stats.decompositions,stats.solves,stats.final_time);
        }
    }
    return 0;
}
