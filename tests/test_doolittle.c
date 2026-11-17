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

/* Reconstruction checks the same L*U entries as the reference tests. */
static void check_matrix(size_t n, const double *original, double tolerance)
{
    unsigned char pattern[100] = {0};
    for (size_t i=0;i<n*n;++i) pattern[i]=original[i]!=0;
    ChemSparsePlan *plan=chem_sparse_plan_from_pattern(n,pattern);
    assert(plan);
    const double initial_values[] = {-INFINITY,-1,0,1,INFINITY};
    for (size_t sparse=0;sparse<2;++sparse)
        for (size_t separate=0;separate<2;++separate)
            for (size_t init=0;init<5;++init) {
                double a[100], original_copy[100];
                memcpy(original_copy,original,n*n*sizeof(double));
                for (size_t i=0;i<n*n;++i) a[i]=initial_values[init];
                int ok;
                if (separate)
                    ok=sparse ? chem_sparse_factor_doolittle_copy(plan,original,a) : chem_dense_factor_doolittle_copy(original,a,n);
                else {
                    memcpy(a,original,n*n*sizeof(double));
                    ok=sparse ? chem_sparse_factor_doolittle(plan,a) : chem_dense_factor_doolittle(a,n);
                }
                assert(ok);
                assert(!memcmp(original,original_copy,n*n*sizeof(double)));
                for (size_t i=0;i<n;++i)
                    for (size_t j=0;j<n;++j) {
                        double value=0;
                        for (size_t k=0;k<n;++k) {
                            double l = i==k ? 1 : (i>k ? a[i*n+k] : 0);
                            double u = k<=j ? a[k*n+j] : 0;
                            value += l*u;
                        }
                        assert(fabs(value-original[i*n+j]) <= tolerance);
                    }
            }
    chem_sparse_plan_destroy(plan);
}

static uint32_t next(uint32_t *seed) { *seed=*seed*1664525u+1013904223u; return *seed; }
int main(void)
{
    const double dense[] = {2,-1,-2,-4,6,3,-4,-2,8};
    check_matrix(3,dense,1e-10);
    uint32_t seed=12345;
    for (size_t block=0;block<5;++block) {
        double a[100];
        for (size_t i=0;i<10;++i)
            for (size_t j=0;j<10;++j) {
                uint32_t x=next(&seed);
                a[i*10+j] = i==j ? 8+(x%100)/100.0 : (x&0x1000u ? (x%300)/200.0 : 0);
            }
        check_matrix(10,a,1e-12);
        double diagonal[36]={0};
        for (size_t i=0;i<6;++i) diagonal[i*6+i]=exp(((double)(next(&seed)%1600)-800)/100.0);
        check_matrix(6,diagonal,1e-10);
        double small[25];
        for (size_t i=0;i<5;++i)
            for (size_t j=0;j<5;++j) {
                uint32_t x=next(&seed);
                small[i*5+j] = i==j ? 8+(x%100)/100.0 : (x&0x1000u ? (x%300)/200.0 : 0);
            }
        check_matrix(5,small,1e-12);
    }
    puts("Doolittle decomposition checks PASS");
    return 0;
}
