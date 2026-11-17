#include "chem.h"
#include "rosenbrock.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>

typedef struct { const ChemMechanism *mechanism; const double *rates; } Context;
static void rhs(void *context, const double *y, double *out)
{ Context *c = context; chem_rhs(c->mechanism,y,c->rates,out); }
static void jac(void *context, const double *y, double *out)
{ Context *c = context; chem_jacobian(c->mechanism,y,c->rates,out); }
int main(void)
{
    const size_t reactants[] = {0,1};
    const ChemProduct products[] = {{1,1},{2,1}};
    const ChemReaction reactions[] = {{0,1,0,1},{1,1,1,1}};
    const ChemMechanism m = {3,2,reactions,reactants,products};
    const double k[] = {.002,.003}, atol[] = {1e-3,1e-3,1e-3};
    Context context = {&m,k};
    ChemRosenbrockMethod methods[] = {chem_rosenbrock_two_stage(),chem_rosenbrock_three_stage(),
        chem_rosenbrock_four_stage(),chem_rosenbrock_four_stage_da(),chem_rosenbrock_six_stage_da()};
    ChemSparsePlan *plan = chem_sparse_plan_create(&m);
    assert(plan);
    for (size_t method=0; method<5; ++method) {
        for (size_t sparse=0; sparse<2; ++sparse) {
            ChemRosenbrockOptions options = chem_rosenbrock_defaults(atol);
            options.sparse_plan = sparse ? plan : NULL;
            for (size_t cell=0; cell<3; ++cell) {
                double a0=1.0-(double)cell/3.0, y[] = {a0,0,0};
                for (size_t step=1; step<1000; ++step) {
                    ChemSolverStats stats;
                    assert(chem_rosenbrock_solve(3,y,1,rhs,jac,&context,&methods[method],&options,&stats)==CHEM_SOLVE_OK);
                    double t=(double)step;
                    double e1=exp(-k[0]*t), e2=exp(-k[1]*t);
                    double expected[] = {a0*e1,a0*k[0]/(k[1]-k[0])*(e1-e2),
                        a0*(1+(k[0]*e2-k[1]*e1)/(k[1]-k[0]))};
                    for (size_t i=0; i<3; ++i) assert(fabs(y[i]-expected[i]) <= 2e-6);
                }
            }
        }
    }
    chem_sparse_plan_destroy(plan);

    const size_t stiff_reactants[] = {0,1,2,0,1};
    const ChemProduct stiff_products[] = {{2,1},{2,1},{3,1},{1,1},{0,1}};
    const ChemReaction stiff_reactions[] = {{0,1,0,1},{1,1,1,1},{2,1,2,1},{3,1,3,1},{4,1,4,1}};
    const ChemMechanism stiff = {4,5,stiff_reactions,stiff_reactants,stiff_products};
    const double stiff_rates[] = {.002,.002,.003,4e10,3.6e10};
    const double stiff_atol[] = {1e-3,1e-3,1e-3,1e-3};
    Context stiff_context = {&stiff,stiff_rates};
    plan = chem_sparse_plan_create(&stiff);
    assert(plan);
    for (size_t method=0; method<5; ++method)
        for (size_t sparse=0; sparse<2; ++sparse) {
            ChemRosenbrockOptions options=chem_rosenbrock_defaults(stiff_atol);
            options.sparse_plan=sparse ? plan : NULL;
            for (size_t cell=0; cell<3; ++cell) {
                double a0=1.0-(double)cell/3, y[] = {.5*a0,.5*a0,0,0};
                for (size_t step=1; step<1000; ++step) {
                    assert(chem_rosenbrock_solve(4,y,1,rhs,jac,&stiff_context,&methods[method],&options,NULL)==CHEM_SOLVE_OK);
                    double t=(double)step, e1=exp(-.002*t), e2=exp(-.003*t);
                    double expected[] = {a0*e1,a0*(.002/.001)*(e1-e2),
                        a0*(1+(.002*e2-.003*e1)/.001)};
                    assert(fabs(y[0]+y[1]-expected[0]) <= 2e-5);
                    assert(fabs(y[2]-expected[1]) <= 2e-5);
                    assert(fabs(y[3]-expected[2]) <= 2e-5);
                }
            }
        }
    chem_sparse_plan_destroy(plan);
    puts("five Rosenbrock methods: analytical and stiff photolysis PASS");
    return 0;
}
