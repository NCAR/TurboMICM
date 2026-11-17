#include "backward_euler.h"
#include "rates.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>

typedef struct { const ChemMechanism *mechanism; double k[2]; } Chain;
static void chain_rhs(void *context, const double *y, double *out)
{
    Chain *chain = context; chem_rhs(chain->mechanism,y,chain->k,out);
}
static void chain_jac(void *context, const double *y, double *out)
{
    Chain *chain = context; chem_jacobian(chain->mechanism,y,chain->k,out);
}
static void rhs(void *context, const double *y, double *out)
{
    double k = *(double *)context;
    out[0] = -k*y[0]; out[1] = k*y[0];
}
static void jac(void *context, const double *y, double *out)
{
    double k = *(double *)context; (void)y;
    out[0] = -k; out[1] = 0; out[2] = k; out[3] = 0;
}
static void nonfinite_rhs(void *context,const double *y,double *out)
{ (void)y; out[0]=*(const double *)context; }
static void zero_jac(void *context,const double *y,double *out)
{ (void)context; (void)y; out[0]=0; }
int main(void)
{
    double residual[12] = {0}, state[12] = {0}, abs_tol[3] = {1e-6,1e-6,1e-6};
    double small = 1e-6, rtol = 1e-3;
#define CONVERGED() chem_backward_euler_converged(4,3,residual,state,small,abs_tol,rtol)
    assert(CONVERGED()); residual[1] = 1e-5; assert(!CONVERGED());
    small = 1e-4; assert(CONVERGED()); residual[11] = 1e-3; assert(!CONVERGED());
    state[11] = 10; assert(CONVERGED()); residual[11] = 1e-1; assert(!CONVERGED());
    abs_tol[2] = 1; assert(CONVERGED()); residual[11] = INFINITY; assert(!CONVERGED());
#undef CONVERGED
    const size_t reactants[] = {0};
    const ChemProduct products[] = {{1,1}};
    const ChemReaction reactions[] = {{0,1,0,1}};
    const ChemMechanism mechanism = {2,1,reactions,reactants,products};
    const double atol[] = {1e-12,1e-12};
    double y[] = {1,0}, k = 2;
    ChemBackwardEulerOptions options = chem_backward_euler_defaults(atol);
    ChemSolverStats stats;
    assert(chem_backward_euler_solve(2,y,1,rhs,jac,&k,&options,&stats) == CHEM_SOLVE_OK);
    assert(fabs(y[0]-1.0/3) < 1e-15 && fabs(y[1]-2.0/3) < 1e-15);
    assert(stats.steps == 2 && stats.accepted == 1 && stats.rejected == 0 && stats.final_time == 1);
    ChemSparsePlan *plan = chem_sparse_plan_create(&mechanism);
    assert(plan); options.sparse_plan = plan; y[0]=1; y[1]=0;
    assert(chem_backward_euler_solve(2,y,1,rhs,jac,&k,&options,&stats) == CHEM_SOLVE_OK);
    assert(fabs(y[0]-1.0/3) < 1e-15 && fabs(y[1]-2.0/3) < 1e-15);
    options.lu_algorithm = CHEM_LU_DOOLITTLE; y[0]=1; y[1]=0;
    assert(chem_backward_euler_solve(2,y,1,rhs,jac,&k,&options,&stats) == CHEM_SOLVE_OK);
    assert(fabs(y[0]-1.0/3) < 1e-15 && fabs(y[1]-2.0/3) < 1e-15);
    options.sparse_plan = NULL; y[0]=1; y[1]=0;
    assert(chem_backward_euler_solve(2,y,1,rhs,jac,&k,&options,&stats) == CHEM_SOLVE_OK);
    assert(fabs(y[0]-1.0/3) < 1e-15 && fabs(y[1]-2.0/3) < 1e-15);
    options.sparse_plan = plan;
    options.max_iterations = 1; y[0]=1; y[1]=0;
    assert(chem_backward_euler_solve(2,y,1,rhs,jac,&k,&options,&stats) == CHEM_SOLVE_UNCONVERGED);
    assert(stats.rejected == 6 && stats.accepted == 0 && fabs(stats.final_time-.00625) < 1e-15);
    {
        double rates[]={2,4,8}; void *contexts[]={rates,rates+1,rates+2};
        const unsigned char pattern[]={1,0,1,0};
        for (size_t sparse=0;sparse<2;++sparse) for (size_t lu=0;lu<2;++lu) {
            double group_y[]={1,0,2,0,3,0};
            ChemBackwardEulerOptions group=chem_backward_euler_defaults(atol); group.lu_algorithm=(ChemLuAlgorithm)lu;
            assert(chem_backward_euler_solve_grouped(3,2,group_y,1,rhs,jac,contexts,sparse ? pattern : NULL,&group,&stats)==CHEM_SOLVE_OK);
            assert(stats.steps==2 && stats.accepted==1 && stats.final_time==1);
            for (size_t c=0;c<3;++c) {
                double expected=(c+1)/(1+rates[c]);
                assert(fabs(group_y[2*c]-expected)<1e-14 && fabs(group_y[2*c+1]-((c+1)-expected))<1e-14);
            }
            group.max_iterations=1;
            assert(chem_backward_euler_solve_grouped(3,2,group_y,1,rhs,jac,contexts,sparse ? pattern : NULL,&group,&stats)==CHEM_SOLVE_UNCONVERGED);
            assert(stats.rejected==6 && stats.accepted==0 && fabs(stats.final_time-.00625)<1e-15);
        }
    }
    chem_sparse_plan_destroy(plan);
    /* The reference CanCallSolve case uses this chain and conditions. */
    const size_t chain_reactants[] = {0,1};
    const ChemProduct chain_products[] = {{1,1},{2,1}};
    const ChemReaction chain_reactions[] = {{0,1,0,1},{1,1,1,1}};
    const ChemMechanism chain_mechanism = {3,2,chain_reactions,chain_reactants,chain_products};
    ChemRate chain_rates[] = {{CHEM_RATE_ARRHENIUS,{2.15e-11,0,110,300,0},0},
                              {CHEM_RATE_ARRHENIUS,{3.3e-11,0,55,300,0},0}};
    ChemConditions conditions = {272.5,101253.3,1e6};
    Chain chain = {&chain_mechanism,{0,0}};
    chem_rates(2,chain_rates,conditions,NULL,chain.k);
    const double chain_atol[] = {1e-6,1e-6,1e-6};
    options = chem_backward_euler_defaults(chain_atol);
    double chain_y[] = {1,0,0};
    assert(chem_backward_euler_solve(3,chain_y,1,chain_rhs,chain_jac,&chain,&options,&stats) == CHEM_SOLVE_OK);
    assert(fabs(chain_y[0]+chain_y[1]+chain_y[2]-1) < 1e-14);
    const double one_atol[]={1};
    ChemBackwardEulerOptions extreme=chem_backward_euler_defaults(one_atol); extreme.max_iterations=1;
    double value=1,forcing=INFINITY;
    assert(chem_backward_euler_solve(1,&value,1,nonfinite_rhs,zero_jac,&forcing,&extreme,NULL)==CHEM_SOLVE_INFINITY);
    value=1; forcing=NAN;
    assert(chem_backward_euler_solve(1,&value,1,nonfinite_rhs,zero_jac,&forcing,&extreme,NULL)==CHEM_SOLVE_NAN);
    puts("backward Euler convergence and solve PASS");
    return 0;
}
