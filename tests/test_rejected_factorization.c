#undef NDEBUG
#include <assert.h>
#include <math.h>
#include "chem.h"
#include "rates.h"
#include "rosenbrock.h"

typedef struct { const ChemMechanism *mechanism; double rates[3]; } System;
static void rhs(void *v,const double *y,double *out)
{ System *s=v; chem_rhs(s->mechanism,y,s->rates,out); }
static void jacobian(void *v,const double *y,double *out)
{ System *s=v; chem_jacobian(s->mechanism,y,s->rates,out); }
int main(void)
{
    const size_t reactants[]={0,2,1,3};
    const ChemProduct products[]={{1,1},{4,2.4},{0,1},{3,1.4}};
    const ChemReaction reactions[]={{0,2,0,2},{2,1,2,2},{3,1,4,0}};
    const ChemMechanism mechanism={5,3,reactions,reactants,products};
    const ChemRate descriptors[]={
        {CHEM_RATE_ARRHENIUS,{2e-11,0,110,300,0},0},
        {CHEM_RATE_ARRHENIUS,{1e-6,0,0,300,0},0},
        {CHEM_RATE_ARRHENIUS,{3.5e-6,0,0,300,0},0}};
    System system={.mechanism=&mechanism};
    chem_rates(3,descriptors,(ChemConditions){298.15,0,0},NULL,system.rates);
    const double atol[]={.001,.001,.001,.001,.001};
    ChemRosenbrockMethod method=chem_rosenbrock_three_stage();
    ChemSparsePlan *plan=chem_sparse_plan_create(&mechanism); assert(plan);
    for (int sparse=0;sparse<2;++sparse) for (int algorithm=0;algorithm<2;++algorithm) {
        ChemSolverStats reference={0}; double first[5]={0};
        for (int separate=0;separate<2;++separate) {
            double y[]={1e12,1e12,1e12,1e12,1e12};
            ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
            options.h_start=1000; options.factor_out_of_place=separate;
            options.sparse_plan=sparse ? plan : NULL; options.lu_algorithm=(ChemLuAlgorithm)algorithm;
            ChemSolverStats stats;
            assert(chem_rosenbrock_solve(5,y,1000,rhs,jacobian,&system,&method,&options,&stats)==CHEM_SOLVE_OK);
            assert(stats.steps-stats.accepted>=2);
            if (!separate) { reference=stats; for (size_t i=0;i<5;++i) first[i]=y[i]; }
            else {
                assert(reference.steps==stats.steps && reference.accepted==stats.accepted && reference.rejected==stats.rejected);
                for (size_t i=0;i<5;++i) assert(fabs(first[i]-y[i])<=fabs(first[i])*1e-14);
            }
        }
    }
    chem_sparse_plan_destroy(plan);
    return 0;
}
