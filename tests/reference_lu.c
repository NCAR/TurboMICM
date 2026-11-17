#include "chem.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

enum { MAX_N=10, MAX_CELLS=5, MAX_VALUES=MAX_N*MAX_N*MAX_CELLS };
typedef struct {
    size_t n, cells, nonzeros, rows[MAX_N+1];
    int columns[MAX_N*MAX_N], lookup[MAX_N*MAX_N];
    unsigned char algebraic[MAX_N];
    double input[MAX_VALUES], factor[MAX_VALUES];
} Fixture;

static uint32_t random_state;
static double uniform(void)
{
    random_state=random_state*1664525u+1013904223u;
    return ((double)random_state+1.0)/4294967297.0;
}
static double lognormal(double sigma)
{
    double z=sqrt(-2.0*log(uniform()))*cos(6.2831853071795864769*uniform());
    return exp(-2.0+sigma*z);
}

/* family: 0 is the reference's 3x3 dense example, 1 is a random sparse
 * 10x10 matrix, 2 is a six-element diagonal, 3 exercises initial values. */
static void setup(Fixture *fixture, int family, size_t cells, int separate, double initial)
{
    memset(fixture,0,sizeof(*fixture));
    fixture->n=family==0 ? 3 : family==2 ? 6 : family==3 && !separate ? 5 : 10;
    fixture->cells=cells;
    size_t n=fixture->n;
    unsigned char original[MAX_N*MAX_N]={0}, pattern[MAX_N*MAX_N]={0};
    random_state=12973u+(unsigned)family;
    for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) {
        int present=family==0 || i==j || (family!=2 && uniform()>.5);
        original[i*n+j]=pattern[i*n+j]=(unsigned char)present;
    }
    for (size_t i=0;i<n;++i) for (size_t k=i+1;k<n;++k) if (pattern[i*n+k])
        for (size_t j=i+1;j<n;++j) if (pattern[j*n+i]) pattern[j*n+k]=1;
    for (size_t i=0;i<n*n;++i) fixture->lookup[i]=-1;
    for (size_t i=0;i<n;++i) {
        fixture->rows[i]=fixture->nonzeros;
        for (size_t j=0;j<n;++j) if (pattern[i*n+j]) {
            size_t index=fixture->nonzeros++;
            fixture->lookup[i*n+j]=(int)index; fixture->columns[index]=(int)j;
        }
    }
    fixture->rows[n]=fixture->nonzeros;
    const double dense[9]={2,-1,-2,-4,6,3,-4,-2,8};
    for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) if (pattern[i*n+j]) {
        size_t pos=(size_t)fixture->lookup[i*n+j];
        double shared=family==0 ? dense[i*n+j] : original[i*n+j] ? lognormal(family==2 ? 4 : 2) : 0;
        for (size_t cell=0;cell<cells;++cell) {
            double value=shared;
            if (separate && family==1 && original[i*n+j]) value=lognormal(2);
            if (family==2 && cell && i==j) value=lognormal(4);
            fixture->input[pos*cells+cell]=value;
            fixture->factor[pos*cells+cell]=initial;
        }
    }
}

static ChemLinearSystem view(Fixture *fixture, double *values)
{
    return (ChemLinearSystem){fixture->n,fixture->cells,fixture->nonzeros,fixture->rows,fixture->columns,
                              fixture->lookup,fixture->algebraic,values};
}

static double element(const Fixture *f, const double *values, size_t i, size_t j, size_t cell)
{
    int k=f->lookup[i*f->n+j];
    return k<0 ? 0 : values[(size_t)k*f->cells+cell];
}

static int factor_and_check(int family, size_t cells, int separate, int doolittle, double initial, int solve)
{
    Fixture fixture;
    setup(&fixture,family,cells,separate,initial);
    ChemLinearSystem input=view(&fixture,fixture.input), result=view(&fixture,fixture.factor);
    ChemLinearSolver solver=doolittle ? chem_doolittle_linear_solver() : chem_mozart_linear_solver();
    if (separate) {
        if (chem_linear_factor_copy(&solver,&input,&result,0,NULL)!=CHEM_OK) return 0;
    } else {
        memcpy(fixture.factor,fixture.input,fixture.nonzeros*cells*sizeof(double));
        if (solver.factor(&result,0,NULL)!=CHEM_OK) return 0;
    }
    double tolerance=family==0 || family==2 ? 1e-10 : separate ? 1e-9 : 1e-12;
    for (size_t cell=0;cell<cells;++cell) for (size_t i=0;i<fixture.n;++i) for (size_t j=0;j<fixture.n;++j) {
        double sum=0;
        for (size_t k=0;k<fixture.n;++k) {
            double l=i==k ? 1 : i>k ? element(&fixture,fixture.factor,i,k,cell) : 0;
            double u=k<=j ? element(&fixture,fixture.factor,k,j,cell) : 0;
            sum+=l*u;
        }
        if (!isfinite(sum) || fabs(sum-element(&fixture,fixture.input,i,j,cell))>tolerance) return 0;
    }
    if (solve) {
        double rhs[MAX_N*MAX_CELLS], expected[MAX_N*MAX_CELLS];
        random_state=56317;
        for (size_t i=0;i<fixture.n;++i) for (size_t cell=0;cell<cells;++cell) {
            double b=family==0 ? (double[]){23,42,9}[i] : lognormal(family==2 ? 4 : 2);
            rhs[i*cells+cell]=expected[i*cells+cell]=b;
        }
        if (solver.solve(&result,rhs,NULL)!=CHEM_OK) return 0;
        for (size_t i=0;i<fixture.n;++i) for (size_t cell=0;cell<cells;++cell) {
            double sum=0;
            for (size_t j=0;j<fixture.n;++j) sum+=element(&fixture,fixture.input,i,j,cell)*rhs[j*cells+cell];
            double tol=family==2 ? 1e-12*fmax(1,fabs(expected[i*cells+cell])) : family==0 ? 1e-5 : 1e-6;
            if (!isfinite(sum) || fabs(sum-expected[i*cells+cell])>tol) return 0;
        }
    }
    return 1;
}

static int check(int family, int vector, int separate, int doolittle, int solve)
{
    const double initial[]={-INFINITY,-1,0,1,INFINITY};
    size_t cases=family==3 ? 5 : 1;
    for (size_t init=0;init<cases;++init) {
        size_t variants=vector ? 4 : 1;
        for (size_t width=1;width<=variants;++width) {
            size_t cells=family==0 ? 1 : family==1 && !vector && !solve ? 1 : 5;
            /* Exercise the same cell count as the reference as well as each
             * width-sized final group used by the interleaved layout. */
            if (!factor_and_check(family,cells,separate,doolittle,initial[init],solve)) return 0;
            if (vector && cells>width && !factor_and_check(family,width,separate,doolittle,initial[init],solve)) return 0;
        }
        if (family==1 && !vector && !solve && !factor_and_check(family,5,separate,doolittle,initial[init],solve)) return 0;
    }
    return 1;
}

int main(void)
{
    int failures=0;
    const char *prefixes[]={"lu_decomposition_doolittle","lu_decomposition_mozart",
                            "lu_decomposition_doolittle_in_place","lu_decomposition_mozart_in_place"};
    const char *suites[]={"LuDecompositionDoolittle","LuDecompositionMozart",
                          "LuDecompositionDoolittleInPlace","LuDecompositionMozartInPlace"};
    const char *names[]={"DenseMatrixStandardOrdering","RandomMatrixStandardOrdering","DiagonalMatrixStandardOrdering",
                         "AgnosticToInitialValueStandardOrdering","DenseMatrixVectorOrdering","RandomMatrixVectorOrdering",
                         "DiagonalMatrixVectorOrdering","AgnosticToInitialValueVectorOrdering"};
    for (int a=0;a<4;++a) for (int i=0;i<8;++i) {
        int ok=check(i%4,i>=4,a<2,a%2==0,0);
        printf("%s:%s.%s %s %s\n",prefixes[a],suites[a],names[i],ok ? "PASS" : "FAIL",
               ok ? "factor reconstruction matched" : "factor reconstruction differed");
        failures+=!ok;
    }
    const char *lp[]={"linear_solver","linear_solver_in_place"}, *ls[]={"LinearSolver","LinearSolverInPlace"};
    const char *ln[]={"DenseMatrixStandardOrdering","RandomMatrixStandardOrdering","DiagonalMatrixStandardOrdering",
                      "StandardOrderingAgnosticToInitialValue","DenseMatrixVectorOrdering","RandomMatrixVectorOrdering",
                      "VectorOrderingAgnosticToInitialValue","DiagonalMatrixVectorOrdering"};
    for (int a=0;a<2;++a) for (int i=0;i<8;++i) {
        int family=i==6 ? 3 : i==7 ? 2 : i%4;
        int ok=check(family,i>=4,a==0,a==0,1);
        printf("%s:%s.%s %s %s\n",lp[a],ls[a],ln[i],ok ? "PASS" : "FAIL",
               ok ? "linear solve reconstruction matched" : "linear solve reconstruction differed");
        failures+=!ok;
    }
    return failures ? 1 : 0;
}
