#include "linear.h"
#include "sparse.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static uint32_t next(uint32_t *seed) { *seed=*seed*1664525u+1013904223u; return *seed; }
static void check(size_t n, const double *a, const double *b, double tolerance)
{
    unsigned char pattern[900];
    for (size_t i=0;i<n*n;++i) pattern[i]=a[i]!=0;
    ChemSparsePlan *plan=chem_sparse_plan_from_pattern(n,pattern); assert(plan);
    const double initial[] = {-INFINITY,-1,0,1,INFINITY};
    for (size_t algorithm=0;algorithm<2;++algorithm)
      for (size_t sparse=0;sparse<2;++sparse)
        for (size_t separate=0;separate<2;++separate)
          for (size_t t=0;t<5;++t) {
            double lu[900], x[30];
            for (size_t i=0;i<n*n;++i) lu[i]=initial[t];
            int ok;
            if (separate) {
                if (algorithm) ok=sparse ? chem_sparse_factor_doolittle_copy(plan,a,lu) : chem_dense_factor_doolittle_copy(a,lu,n);
                else ok=sparse ? chem_sparse_factor_copy(plan,a,lu) : chem_dense_factor_copy(a,lu,n);
            } else {
                memcpy(lu,a,n*n*sizeof(double));
                if (algorithm) ok=sparse ? chem_sparse_factor_doolittle(plan,lu) : chem_dense_factor_doolittle(lu,n);
                else ok=sparse ? chem_sparse_factor(plan,lu) : chem_dense_factor(lu,n);
            }
            assert(ok); memcpy(x,b,n*sizeof(double));
            if (sparse) chem_sparse_solve(plan,lu,x); else chem_dense_solve(lu,n,x);
            for (size_t i=0;i<n;++i) {
                double result=0; for(size_t j=0;j<n;++j) result+=a[i*n+j]*x[j];
                assert(fabs(result-b[i]) <= tolerance);
            }
          }
    chem_sparse_plan_destroy(plan);
}
int main(void)
{
    const double a[]={2,-1,-2,-4,6,3,-4,-2,8}, b[]={23,42,9};
    check(3,a,b,1e-5);
    uint32_t seed=4723;
    for (size_t block=0;block<5;++block) {
        for (size_t kind=0;kind<3;++kind) {
            size_t n=kind==0 ? 10 : (kind==1 ? 30 : 6);
            double matrix[900]={0}, rhs[30];
            for(size_t i=0;i<n;++i) {
                rhs[i]=.1+(next(&seed)%1000)/100.0;
                for(size_t j=0;j<n;++j) {
                    uint32_t x=next(&seed);
                    if (i==j) matrix[i*n+j]=kind==2 ? exp(((double)(x%1600)-800)/100.0) : 30+(x%100)/100.0;
                    else if (kind!=2 && (x&0x1000u)) matrix[i*n+j]=(x%300)/200.0;
                }
            }
            check(n,matrix,rhs,kind==0 ? 1e-6 : (kind==1 ? 2e-6 : 1e-5));
        }
    }
    puts("linear factor and solve checks PASS"); return 0;
}
