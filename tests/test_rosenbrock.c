#include "rosenbrock.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>

typedef struct { double rate; } Decay;
static void rhs(void *ctx, const double *y, double *out)
{
    const Decay *d = ctx;
    out[0] = -d->rate*y[0];
    out[1] = d->rate*y[0];
}
static void jac(void *ctx, const double *y, double *out)
{
    const Decay *d = ctx;
    (void)y;
    out[0] = -d->rate; out[1] = 0;
    out[2] = d->rate; out[3] = 0;
}
int main(void)
{
    Decay d = {2};
    double y[2] = {1, 0}, atol[2] = {1e-12, 1e-12};
    ChemRosenbrockMethod method = chem_rosenbrock_three_stage();
    ChemRosenbrockOptions options = chem_rosenbrock_defaults(atol);
    ChemSolverStats stats;
    options.relative_tolerance = 1e-9;
    assert(chem_rosenbrock_solve(2, y, 1, rhs, jac, &d, &method, &options, &stats) == CHEM_SOLVE_OK);
    assert(fabs(y[0]-exp(-2)) < 2e-9);
    assert(fabs(y[0]+y[1]-1) < 1e-12);
    assert(stats.steps >= stats.accepted && stats.accepted > 0);
    assert(stats.solves == 3*stats.decompositions);
    assert(stats.jacobian_updates == stats.decompositions);
    assert(chem_rosenbrock_solve(2, y, -1, rhs, jac, &d, &method, &options, &stats) == CHEM_SOLVE_INVALID);
    ChemRosenbrockMethod methods[] = {
        chem_rosenbrock_two_stage(), chem_rosenbrock_three_stage(),
        chem_rosenbrock_four_stage(), chem_rosenbrock_four_stage_da(),
        chem_rosenbrock_six_stage_da()
    };
    options.relative_tolerance = 1e-6;
    options.max_steps = 100000;
    for (size_t i=0; i<sizeof(methods)/sizeof(methods[0]); ++i) {
        y[0]=1; y[1]=0;
        assert(chem_rosenbrock_solve(2,y,1,rhs,jac,&d,&methods[i],&options,&stats) == CHEM_SOLVE_OK);
        assert(fabs(y[0]-exp(-2)) < 2e-6);
        assert(fabs(y[0]+y[1]-1) < 1e-12);
        assert(stats.solves == methods[i].stages*stats.decompositions);
    }
    puts("Rosenbrock decay PASS");
    return 0;
}
