#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "rosenbrock.h"
typedef struct { size_t n,from,to; double k; } Reaction;
static void rhs(void *v,const double *y,double *f)
{
    Reaction *r=v; memset(f,0,r->n*sizeof(double)); f[r->from]-=r->k*y[r->from]; f[r->to]+=r->k*y[r->from];
}
static void jac(void *v,const double *y,double *j)
{
    Reaction *r=v; (void)y; memset(j,0,r->n*r->n*sizeof(double));
    j[r->from*r->n+r->from]-=r->k; j[r->to*r->n+r->from]+=r->k;
}
static void single_balance(void)
{
    const ChemConstraintTerm terms[]={{0,1},{1,1},{2,1}};
    const ChemConstraint c={.name="mass",.type=CHEM_CONSTRAINT_LINEAR,.row=2,.reactants=terms,.reactant_count=3,.constant=1e-6};
    const ChemConstraintSet set={3,1,&c}; const double atol[]={1e-12,1e-12,1e-12};
    ChemRosenbrockMethod methods[]={chem_rosenbrock_three_stage(),chem_rosenbrock_four_stage(),
        chem_rosenbrock_four_stage_da(),chem_rosenbrock_six_stage_da()};
    Reaction reaction={3,0,1,1e4};
    for (size_t m=0;m<4;++m) {
        double y[]={.9e-6,0,.1e-6}; ChemSolverStats stats;
        ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol); options.constraints=&set;
        ChemSolveStatus status=chem_rosenbrock_solve(3,y,30,rhs,jac,&reaction,&methods[m],&options,&stats);
        if (status!=CHEM_SOLVE_OK) fprintf(stderr,"balance method %zu status %d t %g\n",m,status,stats.final_time);
        assert(status==CHEM_SOLVE_OK);
        assert(fabs(y[0]+y[1]+y[2]-1e-6)<1e-12 && y[2]>=-1e-18);
        assert(y[0]<1e-12 && fabs(y[1]-.9e-6)<1e-10 && fabs(y[2]-.1e-6)<1e-10);
    }
}
static void equilibrium_balance(void)
{
    const ChemConstraintTerm a[]={{0,1}},b[]={{1,1}},all[]={{1,1},{2,1},{0,1}};
    const ChemConstraint c[]={
        {.name="eq",.type=CHEM_CONSTRAINT_EQUILIBRIUM,.row=1,.reactants=a,.products=b,.reactant_count=1,.product_count=1,.constant=10,.gas_constant=8.31446261815324,.reference_temperature=298.15},
        {.name="mass",.type=CHEM_CONSTRAINT_LINEAR,.row=0,.reactants=all,.reactant_count=3,.constant=1e-6}};
    const ChemConstraintSet set={3,2,c}; const double atol[]={1e-8,1e-8,1e-12};
    double y[]={1e-6/11,10e-6/11,0}; Reaction reaction={3,1,2,1e3};
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol); options.constraints=&set;
    ChemRosenbrockMethod method=chem_rosenbrock_four_stage_da(); ChemSolverStats stats;
    ChemSolveStatus status=chem_rosenbrock_solve(3,y,30,rhs,jac,&reaction,&method,&options,&stats);
    if (status!=CHEM_SOLVE_OK) fprintf(stderr,"eq balance status %d t %g\n",status,stats.final_time);
    assert(status==CHEM_SOLVE_OK && fabs(y[0]+y[1]+y[2]-1e-6)<1e-12);
    assert(y[0]>=-1e-18 && y[1]>=-1e-18 && fabs(y[2]-1e-6)<1e-8);
}
typedef struct { double gas,product; size_t accepted; } Result;
static Result cascade(double balance_atol,double k,double k1,double k2)
{
    const ChemConstraintTerm gas[]={{0,1}},aq[]={{1,1}},bq[]={{2,1}},all[]={{1,1},{2,1},{3,1},{0,1}};
    const ChemConstraint c[]={
        {.name="eq1",.type=CHEM_CONSTRAINT_EQUILIBRIUM,.row=1,.reactants=gas,.products=aq,.reactant_count=1,.product_count=1,.constant=k1,.gas_constant=8.31446261815324,.reference_temperature=298.15},
        {.name="eq2",.type=CHEM_CONSTRAINT_EQUILIBRIUM,.row=2,.reactants=aq,.products=bq,.reactant_count=1,.product_count=1,.constant=k2,.gas_constant=8.31446261815324,.reference_temperature=298.15},
        {.name="mass",.type=CHEM_CONSTRAINT_LINEAR,.row=0,.reactants=all,.reactant_count=4,.constant=1e-6}};
    const ChemConstraintSet set={4,3,c}; const double atol[]={balance_atol,1e-8,1e-8,1e-12};
    double ag=1e-6/(1+k1+k1*k2),y[]={ag,k1*ag,k1*k2*ag,0}; Reaction reaction={4,2,3,k};
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol); options.constraints=&set;
    ChemRosenbrockMethod method=chem_rosenbrock_four_stage_da(); ChemSolverStats stats;
    ChemSolveStatus status=chem_rosenbrock_solve(4,y,30,rhs,jac,&reaction,&method,&options,&stats);
    if (status!=CHEM_SOLVE_OK) fprintf(stderr,"cascade status %d t %g\n",status,stats.final_time);
    assert(status==CHEM_SOLVE_OK); return (Result){y[0],y[3],stats.accepted};
}
int main(void)
{
    single_balance(); equilibrium_balance();
    Result loose=cascade(1e-3,100,2,2),tight=cascade(1e-8,100,2,2);
    assert(tight.accepted>loose.accepted);
    Result stiff=cascade(1e-8,1e4,100,50);
    assert(stiff.gas>=-1e-8 && fabs(stiff.product-1e-6)<1e-8 && fabs(stiff.gas)<1e-10);
    return 0;
}
