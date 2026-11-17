#include "chem.h"
#include "rosenbrock.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>

extern const ChemMechanism generated_mechanism;
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
    puts("generated Robertson integration PASS");
    return 0;
}
