#include "chem.h"
#include "rates.h"
#include "rosenbrock.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

extern const ChemMechanism generated_mechanism;
extern const ChemRate generated_rates[];
extern const char *const generated_custom_labels[];
extern const size_t generated_custom_count;
typedef struct { const ChemMechanism *mechanism; const double *rates; } System;
static void rhs(void *context, const double *y, double *out)
{
    const System *s = context;
    chem_rhs(s->mechanism, y, s->rates, out);
}
static void jac(void *context, const double *y, double *out)
{
    const System *s = context;
    chem_jacobian(s->mechanism, y, s->rates, out);
}
int main(void)
{
    /* Robertson kinetics: first step has independently tabulated values. */
    const double rates[] = {0.04, 3e7, 1e4};
    const size_t reactants[]={0,1,1,1,2};
    const ChemProduct products[]={{1,1},{1,1},{2,1},{0,1},{2,1}};
    const ChemReaction reactions[]={{0,1,0,1},{1,2,1,2},{3,2,3,2}};
    const ChemMechanism hand={3,3,reactions,reactants,products};
    const ChemRate hand_rates[]={{CHEM_RATE_USER,{1},0},{CHEM_RATE_USER,{1},1},{CHEM_RATE_USER,{1},2}};
    const char *labels[]={"r1","r2","r3"};
    double generated_k[3],hand_k[3],generated_f[3],hand_f[3],generated_j[9],hand_j[9];
    const ChemConditions conditions={300,101325,1};
    assert(generated_custom_count==3);
    for (size_t i=0;i<3;++i) assert(!strcmp(generated_custom_labels[i],labels[i]));
    chem_rates(3,generated_rates,conditions,rates,generated_k);
    chem_rates(3,hand_rates,conditions,rates,hand_k);
    const double probe[]={.9,.003,.097};
    for (size_t i=0;i<3;++i) assert(generated_k[i]==hand_k[i]);
    chem_rhs(&generated_mechanism,probe,generated_k,generated_f);
    chem_rhs(&hand,probe,hand_k,hand_f);
    chem_jacobian(&generated_mechanism,probe,generated_k,generated_j);
    chem_jacobian(&hand,probe,hand_k,hand_j);
    for (size_t i=0;i<3;++i) assert(generated_f[i]==hand_f[i]);
    for (size_t i=0;i<9;++i) assert(generated_j[i]==hand_j[i]);
    const double atol[] = {1e-12, 1e-16, 1e-12};
    double y[] = {1, 0, 0};
    System system = {&generated_mechanism, rates};
    ChemRosenbrockMethod method = chem_rosenbrock_three_stage();
    ChemRosenbrockOptions options = chem_rosenbrock_defaults(atol);
    options.relative_tolerance = 1e-7;
    options.max_steps = 100000;
    ChemSolverStats stats;
    ChemSolveStatus result = chem_rosenbrock_solve(3, y, 0.4, rhs, jac, &system, &method, &options, &stats);
    assert(result == CHEM_SOLVE_OK);
    assert(fabs(y[0] - 0.985172114) < 2e-8);
    assert(fabs(y[1] - 0.00003386395) < 2e-10);
    assert(fabs(y[0] + y[1] + y[2] - 1) < 1e-12);
    double sparse_y[] = {1,0,0};
    ChemSparsePlan *plan = chem_sparse_plan_create(&generated_mechanism);
    assert(plan);
    options.sparse_plan = plan;
    assert(chem_rosenbrock_solve(3, sparse_y, 0.4, rhs, jac, &system, &method, &options, NULL) == CHEM_SOLVE_OK);
    for (size_t i = 0; i < 3; ++i) assert(fabs(sparse_y[i]-y[i]) < 1e-14);
    chem_sparse_plan_destroy(plan);
    options.sparse_plan=NULL;
    System generated_system={&generated_mechanism,generated_k},hand_system={&hand,hand_k};
    double generated_y[]={1,0,0},hand_y[]={1,0,0};
    for (size_t step=0;step<10;++step) {
        assert(chem_rosenbrock_solve(3,generated_y,.4,rhs,jac,&generated_system,&method,&options,NULL)==CHEM_SOLVE_OK);
        assert(chem_rosenbrock_solve(3,hand_y,.4,rhs,jac,&hand_system,&method,&options,NULL)==CHEM_SOLVE_OK);
        for (size_t i=0;i<3;++i) assert(generated_y[i]==hand_y[i]);
    }
    puts("generated Robertson integration PASS");
    return 0;
}
