#include "chem.h"

#include <math.h>
#include <stdio.h>

static int failures;
static void report(const char *id, int ok)
{
    printf("%s %s %s\n", id, ok ? "PASS" : "FAIL", ok ? "reference quantities matched" : "reference quantity differed");
    failures += !ok;
}

static int test_convergence(void)
{
    double residual[12]={0}, y[12]={0}, atol[3]={1e-6,1e-6,1e-6};
    double small=1e-6, rtol=1e-3;
    if (!chem_backward_euler_converged(3,4,residual,y,atol,rtol,small)) return 0;
    residual[1*4+0]=1e-5;
    if (chem_backward_euler_converged(3,4,residual,y,atol,rtol,small)) return 0;
    small=1e-4;
    if (!chem_backward_euler_converged(3,4,residual,y,atol,rtol,small)) return 0;
    residual[2*4+3]=1e-3;
    if (chem_backward_euler_converged(3,4,residual,y,atol,rtol,small)) return 0;
    y[2*4+3]=10;
    if (!chem_backward_euler_converged(3,4,residual,y,atol,rtol,small)) return 0;
    residual[2*4+3]=1e-1;
    if (chem_backward_euler_converged(3,4,residual,y,atol,rtol,small)) return 0;
    atol[2]=1;
    if (!chem_backward_euler_converged(3,4,residual,y,atol,rtol,small)) return 0;
    residual[2*4+3]=INFINITY;
    return !chem_backward_euler_converged(3,4,residual,y,atol,rtol,small);
}

static int test_backward_euler(void)
{
    const ChemSpecies species[]={{"A",1,1e-6,0,0},{"B",0,1e-6,0,0},{"C",0,1e-6,0,0}};
    const int reactants[]={0,1}; const ChemProduct products[]={{1,1},{2,1}};
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=2.15e-11,.c=110,.d=300}},
        {1,1,1,1,0,{.kind=CHEM_ARRHENIUS,.a=3.3e-11,.c=55,.d=300}}
    };
    const ChemMechanism m={3,2,0,species,reactions,reactants,products,NULL,NULL,272.5,101253.3};
    ChemState *s=chem_state_create(&m,1); if (!s) return 0;
    ChemBackwardEuler p=chem_backward_euler_defaults(); ChemStats stats;
    int ok=chem_set_air_density(s,0,1e6)==CHEM_OK && chem_update_rates(s)==CHEM_OK &&
           chem_solve_backward_euler(s,&p,1,&stats)==CHEM_OK && stats.final_time==1;
    chem_state_free(s);
    return ok;
}

static double callback(double t, double p, double d, void *context)
{ (void)p; (void)d; return t * *(const double *)context; }

static int test_override(void)
{
    const double initial=.001, override=.002;
    const ChemSpecies species[]={{"A",1,1e-6,0,0},{"B",0,1e-6,0,0}};
    const int reactants[]={0}; const ChemProduct products[]={{1,1}};
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_CALLBACK,.label="lambda_rc",.function=callback,.context=(void *)&initial}}
    };
    const ChemMechanism m={2,1,0,species,reactions,reactants,products,NULL,NULL,300,101253.3};
    ChemState *s=chem_state_create(&m,1); if (!s) return 0;
    ChemBackwardEuler p=chem_backward_euler_defaults(); ChemStats stats; double f[2], jac[4];
    int ok=chem_override_rate_function(s,"lambda_rc",callback,(void *)&override)==CHEM_OK &&
           chem_set_air_density(s,0,1e6)==CHEM_OK && chem_update_rates(s)==CHEM_OK &&
           chem_evaluate(s,0,f,jac)==CHEM_OK && fabs(f[0]+.6)<=1e-12 &&
           chem_solve_backward_euler(s,&p,1,&stats)==CHEM_OK && stats.final_time==1;
    chem_state_free(s);
    return ok;
}

int main(void)
{
    report("backward_euler:BackwardEuler.CanCallSolve",test_backward_euler());
    report("backward_euler:BackwardEuler.IsConverged",test_convergence());
    report("solver_lambda_rate_constant:Solver.GetLambdaRateConstantByNameCanOverrideLambda",test_override());
    const double vector[]={1,2,3,4,5};
    double sum=0;
    for (size_t i=0;i<sizeof vector/sizeof vector[0];++i) sum+=vector[i];
    report("view_category_sfinae:ViewCategorySFINAE.RuntimeVectorUse",sum==15);
    return failures ? 1 : 0;
}
