#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "rosenbrock.h"

typedef struct { const ChemMechanism *mechanism; const double *rates; } Context;
static void rhs(void *v,const double *y,double *out)
{ Context *c=v; chem_rhs(c->mechanism,y,c->rates,out); }
static void jac(void *v,const double *y,double *out)
{ Context *c=v; chem_jacobian(c->mechanism,y,c->rates,out); }
int main(void)
{
    const size_t reactants[]={0,1,1,2,3,3,3,4};
    const ChemProduct products[]={{1,2},{0,1},{3,1},{3,1},{4,1},{2,1},{4,1}};
    const ChemReaction reactions[]={{0,1,0,1},{1,2,1,1},{3,1,2,1},{4,2,3,2},{6,2,5,2}};
    const ChemMechanism mechanism={5,5,reactions,reactants,products};
    const double rates[]={0,1,.04,3e7,1e4},atol[]={1e-20,1e-20,.001,.001,.001};
    const ChemConstraintTerm terms[]={{2,1},{3,1},{4,1}};
    const ChemConstraint balance={.name="mass_conservation",.type=CHEM_CONSTRAINT_LINEAR,.row=4,
        .reactants=terms,.reactant_count=3,.constant=1};
    const ChemConstraintSet constraints={5,1,&balance};
    Context context={&mechanism,rates};
    ChemRosenbrockMethod method=chem_rosenbrock_four_stage_da();
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
    options.relative_tolerance=1e-8; options.constraints=&constraints; options.constraint_temperature=298;
    unsigned char pattern[25]={0}; assert(chem_jacobian_pattern(&mechanism,pattern));
    chem_constraint_sparsity(&constraints,pattern);
    ChemSparsePlan *plan=chem_sparse_plan_from_pattern(5,pattern); assert(plan);
    options.sparse_plan=plan;
    const double chlorine_initial=2*1.8e-10+1.2e-6;
    double y[]={1.8e-10,1.2e-6,1,0,0},duration=1;
    for (size_t step=0;step<12;++step) {
        double advanced=0;
        while (advanced<duration) {
            ChemSolverStats stats;
            ChemSolveStatus status=chem_rosenbrock_solve(5,y,duration-advanced,rhs,jac,&context,&method,&options,&stats);
            if (!((status==CHEM_SOLVE_OK || status==CHEM_SOLVE_MAX_STEPS) && stats.final_time>0))
                fprintf(stderr,"step %zu status %s advanced=%g of %g\n",step,chem_solve_status_string(status),advanced,duration);
            assert((status==CHEM_SOLVE_OK || status==CHEM_SOLVE_MAX_STEPS) && stats.final_time>0);
            double next=advanced+stats.final_time; assert(next>advanced); advanced=next;
        }
        assert(fabs(y[2]+y[3]+y[4]-1)<=1e-10);
        assert(fabs(2*y[0]+y[1]-chlorine_initial)<=1e-6*chlorine_initial);
        for (size_t i=0;i<5;++i) assert(isfinite(y[i]) && y[i]>=0);
        if (!step) assert(y[2]<.99 && y[3]>0 && y[4]>1e-3);
        duration*=10;
    }
    assert(y[4]>.95 && y[2]<.05 && y[0]>.4*chlorine_initial && y[1]<.1*chlorine_initial);
    chem_sparse_plan_destroy(plan);
    return 0;
}
