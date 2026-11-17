#include "jacobian_verify.h"
#include <math.h>
#include <stdio.h>

static void forcing(size_t block,const double *x,double *f,void *context)
{
    (void)block; (void)context; f[0]+=x[0]*x[1]; f[1]+=x[0]*x[0]-x[1];
}
static int run(int kind)
{
    size_t blocks=kind==1?3:1; ChemMatrix vars={0},fd={0}; ChemSparseBuilder b={0}; ChemSparseMatrix a={0};
    if(chem_matrix_init(&vars,blocks,2,1,0,0) || chem_matrix_init(&fd,blocks,4,1,0,0)) return 0;
    double x[]={3,.5,10},y[]={5,4,.1}; if(kind==1) {x[0]=1;y[0]=2;} if(kind==6) {x[0]=0;y[0]=1;}
    for(size_t i=0;i<blocks;++i) {*chem_matrix_at(&vars,i,0)=x[i];*chem_matrix_at(&vars,i,1)=y[i];}
    int ok=!chem_finite_difference_jacobian(&vars,1e-8,forcing,NULL,&fd);
    if(kind==0 || kind==1 || kind==6) {
        for(size_t i=0;i<blocks;++i) { double expected[]={y[i],x[i],2*x[i],-1};
            for(size_t j=0;j<4;++j) if(fabs(*chem_matrix_at(&fd,i,j)-expected[j])>1e-5) ok=0;
        }
    } else {
        if(chem_sparse_builder_init(&b,2)) ok=0;
        for(size_t i=0;i<2;++i) for(size_t j=0;j<2;++j) if(!(kind==5 && i==1 && j==0)) if(chem_sparse_builder_add(&b,i,j)) ok=0;
        if(chem_sparse_init(&a,&b,1,1,0,0)) ok=0;
        chem_sparse_builder_free(&b);
        if(kind==2 || kind==3) {
            *chem_sparse_at(&a,0,0,0)=kind==3?-999:-5; *chem_sparse_at(&a,0,0,1)=-3;
            *chem_sparse_at(&a,0,1,0)=-6; *chem_sparse_at(&a,0,1,1)=1;
            ChemJacobianComparison result;
            if(chem_compare_jacobian(&a,&fd,1e-7,1e-7,&result)) ok=0;
            if(kind==2 && (!result.passed || result.max_abs_error>=1e-5)) ok=0;
            if(kind==3 && (result.passed || result.worst_row!=0 || result.worst_column!=0)) ok=0;
        } else {
            ChemJacobianComparison result;
            if(chem_check_jacobian_sparsity(&a,&fd,1e-6,&result)) ok=0;
            if(kind==4 && !result.passed) ok=0;
            if(kind==5 && (result.passed || result.worst_row!=1 || result.worst_column!=0)) ok=0;
        }
    }
    chem_sparse_free(&a); chem_matrix_free(&fd); chem_matrix_free(&vars); return ok;
}
int main(void)
{
    const char *names[]={"FiniteDifferenceMatchesAnalytical","MultiBlockFiniteDifference","ComparePassesForCorrectJacobian","CompareFailsForWrongJacobian","SparsityCompletenessPassesWhenComplete","SparsityCompletenessFailsWhenMissingEntry","NearZeroVariableHandled"};
    for(int i=0;i<7;++i) printf("jacobian_verification:JacobianVerification.%s %d\n",names[i],run(i));
    return 0;
}
