#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "chem.h"
#include "rates.h"
#include "rosenbrock.h"

/* Reference species order: M, Ar, CO2, H2O, N2, O1D, O, O2, O3. */
static const size_t reactants[]={4,5, 5,7, 6,8, 0,6,7, 7, 8, 8};
static const ChemProduct products[]={
    {4,1},{6,1}, {6,1},{7,1}, {7,2}, {0,1},{8,1},
    {6,2}, {5,1},{7,1}, {6,1},{7,1}};
static const ChemReaction reactions[]={
    {0,2,0,2},{2,2,2,2},{4,2,4,1},{6,3,5,2},
    {9,1,7,1},{10,1,8,2},{11,1,10,2}};
static const ChemMechanism mechanism={9,7,reactions,reactants,products};
typedef struct { double rates[7]; } System;
static void rhs(void *context,const double *y,double *out)
{ chem_rhs(&mechanism,y,((System *)context)->rates,out); }
static void jac(void *context,const double *y,double *out)
{ chem_jacobian(&mechanism,y,((System *)context)->rates,out); }

/* The reference's 23 stored entries are in this column-major order. */
static const size_t rows[]={0,6,7,8,1,2,3,4,5,6,5,6,6,7,8,5,6,7,8,5,6,7,8};
static const size_t cols[]={0,0,0,0,1,2,3,4,4,4,5,5,6,6,6,7,7,7,7,8,8,8,8};
static void matrix_checks(void)
{
    double lu[81]={0},b[9];
    for (size_t i=0;i<23;++i) lu[rows[i]*9+cols[i]]=1;
    for (size_t i=0;i<9;++i) b[i]=.5;
    chem_dense_solve(lu,9,b);
    const double solution[]={.5,.5,.5,.5,.5,-.5,-1,.5,0};
    for (size_t i=0;i<9;++i) assert(b[i]==solution[i]);

    double a[81]={0};
    for (size_t i=0;i<23;++i) a[rows[i]*9+cols[i]]=-1;
    for (size_t i=0;i<9;++i) a[i*9+i]+=2;
    assert(chem_dense_factor(a,9));
    const double expected[]={1,-1,-1,-1,1,1,1,1,-1,-1,1,-1,1,-1,-1,-1,-2,-1,3,-1,-2,-3,.125};
    for (size_t i=0;i<23;++i) {
        double value=a[rows[i]*9+cols[i]];
        if (rows[i]==cols[i]) value=1/value;
        assert(fabs(value-expected[i])<=.01);
    }
    /* Five inert rows do not participate in the Jacobian-vector product. */
    double matrix[81]={0}, vector[9], product[23]={0};
    for (size_t i=0;i<23;++i) if (rows[i]>=5) matrix[rows[i]*9+cols[i]]=1;
    for (size_t i=0;i<9;++i) vector[i]=.5;
    chem_dense_matvec(9,9,matrix,vector,product);
    const double vector_expected[]={0,0,0,0,0,2,3,2,2};
    for (size_t i=0;i<23;++i) assert(fabs(product[i]-(i<9 ? vector_expected[i] : 0))<=.01);
}
static void forcing_checks(void)
{
    double y[9],k[7],f[9];
    for (size_t i=0;i<9;++i) y[i]=1;
    for (size_t i=0;i<7;++i) k[i]=1;
    chem_rhs(&mechanism,y,k,f);
    const double expected[]={0,0,0,0,0,-1,3,2,-2};
    for (size_t i=0;i<9;++i) assert(f[i]==expected[i]);
    for (size_t i=0;i<9;++i) y[i]=5e-6;
    for (size_t i=0;i<7;++i) k[i]=3e-8;
    chem_rhs(&mechanism,y,k,f);
    const double small[]={0,0,0,0,0,1.49e-13,4.55e-13,1.5e-13,-3e-13};
    for (size_t i=0;i<9;++i)
        if (i<5) assert(f[i]==0); else assert(fabs(f[i]-small[i])<=.01);
}
static void solve_checks(void)
{
    const ChemRate r[]={
        {CHEM_RATE_ARRHENIUS,{2.15e-11,0,110,300,0},0},
        {CHEM_RATE_ARRHENIUS,{3.3e-11,0,55,300,0},0},
        {CHEM_RATE_ARRHENIUS,{8e-12,0,-2060,300,0},0},
        {CHEM_RATE_ARRHENIUS,{6e-34,2.4,0,300,0},0}};
    System s={{0}}; ChemConditions conditions={273.15,100000,2.7e19};
    for (size_t i=0;i<4;++i) s.rates[i]=chem_rate(r+i,conditions,NULL);
    s.rates[4]=1e-4; s.rates[5]=1e-5; s.rates[6]=1e-6;
    const double initial[]={1,.392,.0169,0,32.9,0,0,8.84,0};
    const double expected[]={1,.392,.0169,0,32.9,1.8039e-41,.00176789,8.83912,4.5031e-36};
    double atol[9]; for (size_t i=0;i<9;++i) atol[i]=.001;
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
    options.relative_tolerance=1e-4; options.max_steps=100;
    ChemRosenbrockMethod method=chem_rosenbrock_three_stage(); assert(method.stages==3);
    unsigned random_state=1;
    for (size_t trial=0;trial<4;++trial) {
        double y[9],scale=trial==0 ? 1 : (trial==1 ? 10 : .1);
        for (size_t i=0;i<9;++i) {
            random_state=random_state*1664525u+1013904223u;
            y[i]=trial==3 ? (random_state>>8)/16777216.0 : initial[i]*scale;
        }
        ChemSolveStatus status=chem_rosenbrock_solve(9,y,1,rhs,jac,&s,&method,&options,NULL);
        assert(status==CHEM_SOLVE_OK);
        if (trial<3) for (size_t i=0;i<9;++i) {
            if (fabs(y[i]-scale*expected[i])>1e-4)
                fprintf(stderr,"trial %zu species %zu: %.17g expected %.17g\n",trial,i,y[i],scale*expected[i]);
            assert(fabs(y[i]-scale*expected[i])<=1e-4);
        }
    }
}
int main(void) { matrix_checks(); forcing_checks(); solve_checks(); return 0; }
