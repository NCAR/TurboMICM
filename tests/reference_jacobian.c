#include "chem.h"
#include "chem_jacobian.h"

#include <math.h>
#include <stdio.h>

#define D(m,r,c) (*(double *)chem_matrix_at((m),(r),(c)))
#define S(m,b,r,c) (*(double *)chem_sparse_at((m),(b),(r),(c)))
static int failures;
static void report(const char *name, int ok)
{
    printf("jacobian_verification:JacobianVerification.%s %s %s\n",name,ok ? "PASS" : "FAIL",
           ok ? "reference quantities matched" : "reference quantity differed");
    failures+=!ok;
}
static int forcing(const ChemMatrix *vars, ChemMatrix *out, void *context)
{
    (void)context;
    for (size_t i=0;i<vars->rows;++i) {
        double x=*(const double *)chem_matrix_at_const(vars,i,0), y=*(const double *)chem_matrix_at_const(vars,i,1);
        D(out,i,0)+=x*y; D(out,i,1)+=x*x-y;
    }
    return CHEM_OK;
}
static int finite_difference(size_t cells, int zero)
{
    ChemMatrix vars={0},fd={0}; double z=0;
    int ok=chem_matrix_init(&vars,cells,2,0,sizeof z,&z)==CHEM_OK;
    if (!ok) return 0;
    const double x[]={1,.5,10}, y[]={2,4,.1};
    for (size_t i=0;i<cells;++i) { D(&vars,i,0)=cells==1?(zero?0:3):x[i]; D(&vars,i,1)=cells==1?(zero?1:5):y[i]; }
    ok=chem_finite_difference_jacobian(forcing,NULL,&vars,1e-8,&fd)==CHEM_OK;
    for (size_t i=0;i<cells && ok;++i) {
        double a=D(&vars,i,0),b=D(&vars,i,1),expected[]={b,a,2*a,-1};
        for (size_t j=0;j<4;++j) ok &= fabs(D(&fd,i,j)-expected[j])<=1e-5;
    }
    chem_matrix_destroy(&vars); chem_matrix_destroy(&fd); return ok;
}
static int comparison(int completeness, int wrong)
{
    ChemMatrix vars={0},fd={0}; ChemSparseMatrix a={0}; double z=0;
    const ChemSparseCoordinate full[]={{0,0},{0,1},{1,0},{1,1}},missing[]={{0,0},{0,1},{1,1}};
    int ok=chem_matrix_init(&vars,1,2,0,sizeof z,&z)==CHEM_OK;
    if (!ok) goto done;
    D(&vars,0,0)=3; D(&vars,0,1)=5;
    ok=chem_finite_difference_jacobian(forcing,NULL,&vars,1e-8,&fd)==CHEM_OK &&
       chem_sparse_init(&a,2,1,completeness && wrong ? missing : full,completeness && wrong ? 3 : 4,
                        CHEM_SPARSE_ROW,0,sizeof z,&z)==CHEM_OK;
    if (!ok) goto done;
    if (!completeness) {
        S(&a,0,0,0)=wrong?-999:-5; S(&a,0,0,1)=-3; S(&a,0,1,0)=-6; S(&a,0,1,1)=1;
    }
    ChemJacobianComparison result;
    if (completeness) ok=chem_check_jacobian_sparsity(&a,&fd,1e-6,&result)==CHEM_OK;
    else ok=chem_compare_jacobian(&a,&fd,1e-7,1e-7,&result)==CHEM_OK;
    ok &= result.passed==!wrong;
    if (!wrong && !completeness) ok &= result.max_absolute_error<1e-5;
    if (wrong) ok &= result.worst_row==(completeness?1u:0u) && result.worst_column==0;
done:
    chem_matrix_destroy(&vars); chem_matrix_destroy(&fd); chem_sparse_destroy(&a); return ok;
}
int main(void)
{
    report("FiniteDifferenceMatchesAnalytical",finite_difference(1,0));
    report("MultiBlockFiniteDifference",finite_difference(3,0));
    report("ComparePassesForCorrectJacobian",comparison(0,0));
    report("CompareFailsForWrongJacobian",comparison(0,1));
    report("SparsityCompletenessPassesWhenComplete",comparison(1,0));
    report("SparsityCompletenessFailsWhenMissingEntry",comparison(1,1));
    report("NearZeroVariableHandled",finite_difference(1,1));
    return failures ? 1 : 0;
}
