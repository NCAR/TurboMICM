#undef NDEBUG
#include <assert.h>
#include <math.h>
#include "jacobian_check.h"

static void forcing(void *ctx, size_t cell, const double *x, double *f)
{
    (void)ctx; (void)cell;
    f[0] += x[0] * x[1];
    f[1] += x[0] * x[0] - x[1];
}
static void near(double a, double b) { assert(fabs(a-b) < 1e-5); }
int main(void)
{
    const double states[] = { 3,5, 1,2, .5,4, 10,.1, 0,1 };
    const double expected[] = { 5,3,6,-1, 2,1,2,-1, 4,.5,1,-1, .1,10,20,-1, 1,0,0,-1 };
    double fd[20];
    assert(!chem_finite_difference_jacobian(5, 2, states, forcing, NULL, 1e-8, fd));
    for (size_t i = 0; i < 20; ++i) near(fd[i], expected[i]);
    unsigned char full[] = {1,1,1,1}, missing[] = {1,1,0,1};
    double analytical[] = {-5,-3,-6,1};
    ChemJacobianComparison c = chem_compare_jacobian(1,2,full,analytical,fd,1e-7,1e-7);
    assert(c.passed && c.max_abs_error < 1e-5);
    analytical[0] = -999;
    c = chem_compare_jacobian(1,2,full,analytical,fd,1e-7,1e-7);
    assert(!c.passed && c.worst_row == 0 && c.worst_col == 0);
    c = chem_check_jacobian_sparsity(1,2,full,fd,1e-6);
    assert(c.passed);
    c = chem_check_jacobian_sparsity(1,2,missing,fd,1e-6);
    assert(!c.passed && c.worst_row == 1 && c.worst_col == 0);
    return 0;
}
