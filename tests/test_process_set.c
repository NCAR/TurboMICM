#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <float.h>
#include "chem.h"
#include "jacobian_check.h"
#include "rates.h"
extern const ChemSpeciesCatalog generated_species_catalog;
extern const ChemMechanism generated_mechanism;
extern const ChemRate generated_rates[];
extern const ChemRateMultiplier generated_multipliers[];
extern const ChemParameterizedSpecies generated_multiplier_functions[];

typedef struct { const ChemMechanism *m; const double *rates; } Context;
static void forcing(void *v, size_t cell, const double *x, double *f)
{
    Context *ctx = v;
    chem_add_forcing(ctx->m, x, ctx->rates + 2 * cell, NULL, f);
}
static void test_reference_fixture(void)
{
    /* Generated metadata includes qux, which has a density parameterization. */
    const ChemMechanism m=generated_mechanism;
    const char *names[]={"foo","bar","baz","quz","quuz","qux"};
    assert(chem_species_catalog_used_count(&generated_species_catalog)==6);
    for (size_t i=0;i<6;++i) assert(chem_species_catalog_is_used(&generated_species_catalog,names[i]));
    assert(!chem_species_catalog_is_used(&generated_species_catalog,"corge"));
    assert(!chem_species_catalog_is_used(&generated_species_catalog,"missing"));
    const double states[]={.1,.2,.3,.4,.5,0,1.1,1.2,1.3,1.4,1.5,0};
    const double custom[]={10,20,30,40,110,120,130,140};
    const unsigned char expected_pattern[]={
        1,1,1,0,0,0, 1,1,1,0,0,0, 1,0,1,0,0,0,
        0,1,1,1,0,0, 1,0,1,0,0,0, 0,0,0,0,0,0};
    unsigned char pattern[36]={0}, used[6]={0};
    assert(chem_jacobian_pattern(&m,pattern) && chem_used_species(&m,used));
    for (size_t i=0;i<36;++i) assert(pattern[i]==expected_pattern[i]);
    for (size_t i=0;i<6;++i) assert(used[i]==(i<5));
    for (size_t cell=0;cell<2;++cell) {
        const double *x=states+cell*6;
        double k[4];
        ChemConditions conditions={300,101325,cell ? 80 : 70};
        chem_rates_parameterized(4,generated_rates,generated_multipliers,generated_multiplier_functions,
                                 conditions,custom+4*cell,k);
        double f[]={1000,1000,1000,1000,1000,1000}, jac[36];
        for (size_t i=0;i<36;++i) jac[i]=100;
        chem_add_forcing(&m,x,k,NULL,f);
        chem_subtract_jacobian(&m,x,k,NULL,jac);
        double r1=k[0]*x[0]*x[2],r2=k[1]*x[1],r3=k[2]*x[3],r4=k[3]*x[2];
        double expected[]={1000-r1+r2,1000+r1-r2+r4,1000-r1-r4,
                           1000+1.4*r2-r3+2.5*r4,1000+2.4*r1,1000};
        for (size_t i=0;i<6;++i) assert(fabs(f[i]-expected[i])<=(cell==1 && i==1 ? 1e-12 : 4*DBL_EPSILON)*fmax(1,fabs(expected[i])));
        double reference[36]; for (size_t i=0;i<36;++i) reference[i]=100;
        reference[0]+=k[0]*x[2]; reference[1]-=k[1]; reference[2]+=k[0]*x[0];
        reference[6]-=k[0]*x[2]; reference[7]+=k[1]; reference[8]-=k[0]*x[0]+k[3];
        reference[12]+=k[0]*x[2]; reference[14]+=k[0]*x[0]+k[3];
        reference[19]-=1.4*k[1]; reference[20]-=2.5*k[3]; reference[21]+=k[2];
        reference[24]-=2.4*k[0]*x[2]; reference[26]-=2.4*k[0]*x[0];
        for (size_t i=0;i<36;++i) assert(fabs(jac[i]-reference[i])<=4*DBL_EPSILON*fmax(1,fabs(reference[i])));
    }
}
int main(void)
{
    test_reference_fixture();
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
