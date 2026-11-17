#undef NDEBUG
#include <assert.h>
#include "chem.h"
#include "jacobian_check.h"

typedef struct { const ChemMechanism *m; const double *rates; } Context;
static void forcing(void *v, size_t cell, const double *x, double *f)
{
    Context *ctx = v;
    chem_add_forcing(ctx->m, x, ctx->rates + 2 * cell, NULL, f);
}
int main(void)
{
    const ChemReaction rx[] = {{0,2,0,2}};
    const size_t react[] = {0,1};
    const ChemProduct products[] = {{2,1},{3,1}};
    const ChemMechanism m = {4,1,rx,react,products};
    const unsigned char mask[] = {0,1,0,1};
    const double x[] = {1,2,3,4,5,6,7,8}, k[] = {10,20};
    for (size_t cell=0; cell<2; ++cell) {
        double f[] = {1000,1000,1000,1000}, j[16];
        for (size_t i=0; i<16; ++i) j[i]=500;
        chem_add_forcing(&m,x+4*cell,k+cell,mask,f);
        chem_subtract_jacobian(&m,x+4*cell,k+cell,mask,j);
        double rate = cell ? 600 : 20, da = cell ? 120 : 20, db = cell ? 100 : 10;
        assert(f[0]==1000-rate && f[1]==1000 && f[2]==1000+rate && f[3]==1000);
        assert(j[0]==500+da && j[1]==500+db && j[8]==500-da && j[9]==500-db);
        assert(j[4]==500 && j[5]==500 && j[12]==500 && j[13]==500);
    }
    const ChemReaction rx2[] = {{0,2,0,1},{2,1,1,1}};
    const size_t react2[] = {0,1,2};
    const ChemProduct products2[] = {{2,2},{0,1}};
    const ChemMechanism m2 = {3,2,rx2,react2,products2};
    const double state[] = {1,2,3,.5,4,1.5}, rates[] = {5,10,8,12};
    unsigned char pattern[] = {1,1,1,1,1,0,1,1,1};
    double fd[18], analytic[18] = {0};
    Context ctx = {&m2,rates};
    for (size_t cell=0; cell<2; ++cell)
        chem_subtract_jacobian(&m2,state+3*cell,rates+2*cell,NULL,analytic+9*cell);
    assert(!chem_finite_difference_jacobian(2,3,state,forcing,&ctx,1e-8,fd));
    assert(chem_compare_jacobian(2,3,pattern,analytic,fd,1e-7,1e-7).passed);
    assert(chem_check_jacobian_sparsity(2,3,pattern,fd,1e-6).passed);
    return 0;
}
