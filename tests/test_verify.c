#include "test.h"
#include "chem_verify.h"
static void simple_forcing(void *context,const double *y,double *f)
{
    size_t cells=*(size_t *)context;
    for (size_t c=0;c<cells;++c) { double x=y[c],z=y[cells+c]; f[c]+=x*z; f[cells+c]+=x*x-z; }
}
static void fd_case(size_t cells,int zero)
{
    double y[6]={3,5},fd[12];
    if (cells==3) memcpy(y,(double[]){1,.5,10,2,4,.1},sizeof(y));
    if (zero) { y[0]=0; y[1]=1; }
    REQUIRE(chem_finite_difference(2,cells,y,simple_forcing,&cells,1e-8,fd)==CHEM_OK);
    for (size_t c=0;c<cells;++c) {
        NEAR(fd[c],y[cells+c],1e-5); NEAR(fd[cells+c],y[c],1e-5);
        NEAR(fd[2*cells+c],2*y[c],1e-5); NEAR(fd[3*cells+c],-1,1e-5);
    }
}
TEST(jacobian_verification,JacobianVerification,FiniteDifferenceMatchesAnalytical) { fd_case(1,0); }
TEST(jacobian_verification,JacobianVerification,MultiBlockFiniteDifference) { fd_case(3,0); }
TEST(jacobian_verification,JacobianVerification,NearZeroVariableHandled) { fd_case(1,1); }
static void compare_case(int wrong,int sparsity)
{
    double y[]={3,5},fd[4],a[]={wrong ? -999 : -5,-3,-6,1}; size_t cells=1,map[]={0,1,2,3};
    REQUIRE(chem_finite_difference(2,1,y,simple_forcing,&cells,1e-8,fd)==CHEM_OK);
    if (sparsity && wrong) map[2]=CHEM_NONE;
    chem_jacobian_comparison result=sparsity ? chem_check_sparsity(2,1,map,fd,1e-6) : chem_compare_jacobian(2,1,map,a,fd,1e-7,1e-7);
    CHECK(result.passed==!wrong);
    if (wrong) { CHECK(result.worst_row==(sparsity ? 1u : 0u)); CHECK(result.worst_column==0); }
    else if (!sparsity) CHECK(result.max_absolute_error<1e-5);
}
TEST(jacobian_verification,JacobianVerification,ComparePassesForCorrectJacobian) { compare_case(0,0); }
TEST(jacobian_verification,JacobianVerification,CompareFailsForWrongJacobian) { compare_case(1,0); }
TEST(jacobian_verification,JacobianVerification,SparsityCompletenessPassesWhenComplete) { compare_case(0,1); }
TEST(jacobian_verification,JacobianVerification,SparsityCompletenessFailsWhenMissingEntry) { compare_case(1,1); }
