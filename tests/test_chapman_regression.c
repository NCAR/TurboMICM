#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "chem.h"
#include "chapman_names.h"
#include "matrix.h"
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
/* Independent rate equations for the four reactive species. */
static void reference_rhs(const System *s,const double *y,double *f)
{
    for (size_t i=0;i<9;++i) f[i]=0;
    double q0=s->rates[0]*y[4]*y[5], q1=s->rates[1]*y[5]*y[7];
    double q2=s->rates[2]*y[6]*y[8], q3=s->rates[3]*y[0]*y[6]*y[7];
    double q4=s->rates[4]*y[7],q5=s->rates[5]*y[8],q6=s->rates[6]*y[8];
    f[5]=q5-q0-q1;
    f[6]=2*q4+q6+q0+q1-q2-q3;
    f[7]=-q4+q5+q6+2*q2-q3;
    f[8]=-q5-q6-q2+q3;
}
static void fixed_rhs(void *v,const double *y,double *f) { reference_rhs(v,y,f); }
static void fixed_jacobian(void *v,const double *y,double *j)
{
    System *s=v; const double *k=s->rates; memset(j,0,81*sizeof(*j));
    j[5*9+4]=-k[0]*y[5]; j[5*9+5]=-k[0]*y[4]-k[1]*y[7];
    j[5*9+7]=-k[1]*y[5]; j[5*9+8]=k[5];
    j[6*9+0]=-k[3]*y[6]*y[7]; j[6*9+4]=k[0]*y[5];
    j[6*9+5]=k[0]*y[4]+k[1]*y[7]; j[6*9+6]=-k[2]*y[8]-k[3]*y[0]*y[7];
    j[6*9+7]=2*k[4]+k[1]*y[5]-k[3]*y[0]*y[6]; j[6*9+8]=k[6]-k[2]*y[6];
    j[7*9+0]=-k[3]*y[6]*y[7]; j[7*9+6]=2*k[2]*y[8]-k[3]*y[0]*y[7];
    j[7*9+7]=-k[4]-k[3]*y[0]*y[6]; j[7*9+8]=k[5]+k[6]+2*k[2]*y[6];
    j[8*9+0]=k[3]*y[6]*y[7]; j[8*9+6]=-k[2]*y[8]+k[3]*y[0]*y[7];
    j[8*9+7]=k[3]*y[0]*y[6]; j[8*9+8]=-k[5]-k[6]-k[2]*y[6];
}
static void reference_integrate(const System *s,double *y)
{
    /* RK4 at 0.5 s resolves these slowly varying photolysis trajectories. */
    const double h=.5;
    for (size_t step=0;step<1000;++step) {
        double k1[9],k2[9],k3[9],k4[9],v[9];
        reference_rhs(s,y,k1);
        for (size_t i=0;i<9;++i) v[i]=y[i]+h*k1[i]/2;
        reference_rhs(s,v,k2);
        for (size_t i=0;i<9;++i) v[i]=y[i]+h*k2[i]/2;
        reference_rhs(s,v,k3);
        for (size_t i=0;i<9;++i) v[i]=y[i]+h*k3[i];
        reference_rhs(s,v,k4);
        for (size_t i=0;i<9;++i) y[i]+=h*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
    }
}
/* Deterministic lognormal samples with the reference mean and deviation. */
static double uniform(uint64_t *seed)
{
    const uint64_t modulus=2147483647, multiplier=48271;
    *seed=(*seed*multiplier)%modulus; double low=(double)(*seed-1);
    *seed=(*seed*multiplier)%modulus; double high=(double)(*seed-1),base=2147483646.0;
    return (low+high*base)/(base*base);
}
static void reference_order(size_t *order)
{
    unsigned char pattern[81]={0},alive[9];
    assert(chem_jacobian_pattern(&mechanism,pattern)); memset(alive,1,sizeof(alive));
    for (size_t step=0;step<9;++step) {
        size_t pivot=9,best=(size_t)-1;
        for (size_t i=0;i<9;++i) if (alive[i]) {
            size_t row=0,col=0;
            for (size_t j=0;j<9;++j) if (alive[j]) { row+=pattern[i*9+j]; col+=pattern[j*9+i]; }
            size_t cost=(row-1)*(col-1);
            if (pivot==9 || cost<best) { best=cost; pivot=i; }
        }
        order[step]=pivot; alive[pivot]=0;
        for (size_t i=0;i<9;++i) if (alive[i] && pattern[i*9+pivot])
            for (size_t j=0;j<9;++j) if (alive[j] && pattern[pivot*9+j]) pattern[i*9+j]=1;
    }
}
static void reference_trajectories(void)
{
    const ChemRate r[]={
        {CHEM_RATE_ARRHENIUS,{2.15e-11,0,110,300,0},0},
        {CHEM_RATE_ARRHENIUS,{3.3e-11,0,55,300,0},0},
        {CHEM_RATE_ARRHENIUS,{8e-12,0,-2060,300,0},0},
        {CHEM_RATE_ARRHENIUS,{6e-34,2.4,0,300,0},0}};
    const double temperatures[]={284.19,215.02,299.31},pressures[]={101245,100789.2,101398};
    const double photos[3][3]={{1e-4,1e-5,1e-6},{3.2e-4,7.3e-5,3.2e-6},{5.2e-4,8.2e-5,4.6e-6}};
    System systems[3]={0}; void *contexts[]={systems,systems+1,systems+2};
    double initial[27],expected[27],fixed[27],max_value=0,spare=0; int hot=0; uint64_t seed=1;
    size_t order[9]; reference_order(order);
    for (size_t i=0;i<27;++i) {
        double z;
        if (hot) { z=spare; hot=0; }
        else {
            double u,v,s;
            do { u=2*uniform(&seed)-1; v=2*uniform(&seed)-1; s=u*u+v*v; } while (s>1 || s==0);
            double scale=sqrt(-2*log(s)/s); z=u*scale; spare=v*scale; hot=1;
        }
        double sample=exp(-2+2*z); initial[(i/9)*9+order[i%9]]=sample;
        if (sample>max_value) max_value=sample;
    }
    memcpy(expected,initial,sizeof(expected)); memcpy(fixed,initial,sizeof(fixed));
    for (size_t c=0;c<3;++c) {
        ChemConditions conditions={temperatures[c],pressures[c],0};
        for (size_t i=0;i<4;++i) systems[c].rates[i]=chem_rate(r+i,conditions,NULL);
        for (size_t i=0;i<3;++i) systems[c].rates[4+i]=photos[c][i];
        reference_integrate(systems+c,expected+9*c);
    }
    double atol[9]; for (size_t i=0;i<9;++i) atol[i]=.001;
    ChemRosenbrockMethod methods[]={chem_rosenbrock_two_stage(),chem_rosenbrock_three_stage(),
        chem_rosenbrock_four_stage(),chem_rosenbrock_four_stage_da(),chem_rosenbrock_six_stage_da()};
    unsigned char pattern[81]={0}; assert(chem_jacobian_pattern(&mechanism,pattern));
    ChemRosenbrockOptions fixed_options=chem_rosenbrock_defaults(atol);
    fixed_options.relative_tolerance=1e-4; fixed_options.h_start=1e-6; fixed_options.max_steps=100;
    for (size_t c=0;c<3;++c)
        assert(chem_rosenbrock_solve(9,fixed+9*c,500,fixed_rhs,fixed_jacobian,systems+c,
                                     methods+1,&fixed_options,NULL)==CHEM_SOLVE_OK);
    for (size_t m=0;m<5;++m) for (size_t sparse=0;sparse<2;++sparse) for (size_t lu=0;lu<2;++lu)
      for (size_t width=0;width<=4;++width) {
        double y[27]; ChemMatrix state={0};
        assert(chem_matrix_init(&state,3,9,width,0)==CHEM_MATRIX_OK);
        for (size_t c=0;c<3;++c) {
            assert(chem_matrix_set_row(&state,c,initial+9*c,9)==CHEM_MATRIX_OK);
            assert(chem_matrix_get_row(&state,c,y+9*c,9)==CHEM_MATRIX_OK);
        }
        ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol); options.lu_algorithm=(ChemLuAlgorithm)lu;
        assert(chem_rosenbrock_solve_grouped(3,9,y,500,rhs,jac,contexts,sparse ? pattern : NULL,
                                            methods+m,&options,NULL)==CHEM_SOLVE_OK);
        for (size_t c=0;c<3;++c) {
            assert(chem_matrix_set_row(&state,c,y+9*c,9)==CHEM_MATRIX_OK);
            assert(chem_matrix_get_row(&state,c,y+9*c,9)==CHEM_MATRIX_OK);
        }
        double relative=m==0 ? 1e-2 : 1e-4;
        for (size_t i=0;i<27;++i) {
            double tolerance=(fabs(y[i])+fabs(expected[i]))*relative+max_value*1e-12;
            if (fabs(y[i]-expected[i])>tolerance)
                fprintf(stderr,"Chapman method=%zu sparse=%zu lu=%zu i=%zu error=%g tolerance=%g\n",
                        m,sparse,lu,i,fabs(y[i]-expected[i]),tolerance);
            assert(fabs(y[i]-expected[i])<=tolerance);
            tolerance=(fabs(y[i])+fabs(fixed[i]))*relative+max_value*1e-12;
            assert(fabs(y[i]-fixed[i])<=tolerance);
        }
        chem_matrix_destroy(&state);
    }
}
static void name_checks(void)
{
    const char *const expected_reactions[]={"O2_1","O3_1","O3_2","N2_O1D_1","O1D_O2_1","O_O3_1","M_O_O2_1"};
    const char *const expected_species[]={"M","Ar","CO2","H2O","N2","O1D","O","O2","O3"};
    size_t count=0; const char *const *names=chem_chapman_reaction_names(&count);
    assert(count==7);
    for (size_t i=0;i<count;++i) assert(!strcmp(names[i],expected_reactions[i]));
    names=chem_chapman_photolysis_names(&count); assert(count==3);
    for (size_t i=0;i<count;++i) assert(!strcmp(names[i],expected_reactions[i]));
    names=chem_chapman_species_names(&count); assert(count==9);
    for (size_t i=0;i<count;++i) assert(!strcmp(names[i],expected_species[i]));
    assert(chem_chapman_species_names(NULL)==names);
}
int main(void) { name_checks(); matrix_checks(); forcing_checks(); solve_checks(); reference_trajectories(); return 0; }
