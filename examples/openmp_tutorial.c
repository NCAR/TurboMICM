#include "chem.h"
#include "rates.h"
#include "rosenbrock.h"
#include "state.h"
#include <omp.h>
#include <stdio.h>

#ifdef TEST_MODE
enum { THREADS=8 };
#else
enum { THREADS=3 };
#endif
typedef struct { const ChemMechanism *mechanism; double rates[3]; } Context;
static void rhs(void *v,const double *y,double *out)
{ Context *c=v; chem_rhs(c->mechanism,y,c->rates,out); }
static void jacobian(void *v,const double *y,double *out)
{ Context *c=v; chem_jacobian(c->mechanism,y,c->rates,out); }
int main(void)
{
    const size_t reactants[]={0,1,1,1,2};
    const ChemProduct products[]={{1,1},{1,1},{2,1},{0,1},{2,1}};
    const ChemReaction reactions[]={{0,1,0,1},{1,2,1,2},{3,2,3,2}};
    const ChemMechanism mechanism={3,3,reactions,reactants,products};
    const ChemRate rate_descriptors[]={{CHEM_RATE_USER,{1},0},{CHEM_RATE_USER,{1},1},{CHEM_RATE_USER,{1},2}};
    const char *const names[]={"A","B","C"}, *const labels[]={"r1","r2","r3"};
    const double atol[]={.001,.001,.001};
    const ChemStateConfig config={1,3,3,3,names,labels,1e-6,atol};
    ChemRosenbrockMethod method=chem_rosenbrock_three_stage();
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
    ChemSparsePlan *plan=chem_sparse_plan_create(&mechanism);
    if (!plan) return 1;
    options.sparse_plan=plan;
    double results[THREADS][3]={{0}};
    int success[THREADS]={0};
    omp_set_dynamic(0);
#pragma omp parallel num_threads(THREADS)
    {
        int thread=omp_get_thread_num();
        ChemState state={0};
        int ok=chem_state_init(&state,&config)==CHEM_STATE_OK;
        if (ok) {
#ifndef TEST_MODE
#pragma omp critical
            printf("Running solver on thread %d\n",thread);
#endif
            const double one=1, k1=.04,k2=3e7,k3=1e4;
            ok=chem_state_set_concentration(&state,"A",&one,1)==CHEM_STATE_OK &&
               chem_state_set_parameter(&state,"r1",&k1,1)==CHEM_STATE_OK &&
               chem_state_set_parameter(&state,"r2",&k2,1)==CHEM_STATE_OK &&
               chem_state_set_parameter(&state,"r3",&k3,1)==CHEM_STATE_OK;
            state.conditions[0]=(ChemConditions){272.5,101253.3,1e6};
            Context context={&mechanism,{0}};
            for (size_t step=0;step<10 && ok;++step) {
                chem_rates(3,rate_descriptors,state.conditions[0],state.parameters,context.rates);
                double elapsed=0;
                while (elapsed<200 && ok) {
                    ChemSolverStats stats;
                    ChemSolveStatus status=chem_rosenbrock_solve(3,state.concentrations,200-elapsed,
                                             rhs,jacobian,&context,&method,&options,&stats);
                    if ((status!=CHEM_SOLVE_OK && status!=CHEM_SOLVE_MAX_STEPS) || stats.final_time<=0) ok=0;
                    else elapsed+=stats.final_time;
                }
            }
            for (size_t i=0;i<3;++i) results[thread][i]=state.concentrations[i];
        }
        success[thread]=ok;
        chem_state_destroy(&state);
    }
    chem_sparse_plan_destroy(plan);
    for (size_t thread=0;thread<THREADS;++thread) {
        if (!success[thread]) return 1;
        for (size_t i=0;i<3;++i) if (results[thread][i]!=results[0][i]) return 1;
    }
#ifndef TEST_MODE
    for (size_t thread=0;thread<THREADS;++thread) {
        printf("Thread %zu\n%10s,%10s,%10s\n",thread+1,"A","B","C");
        printf("%10.2e,%10.2e,%10.2e\n",results[thread][0],results[thread][1],results[thread][2]);
    }
#endif
    return 0;
}
