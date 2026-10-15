/* Copyright (C) 2024-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "test.h"
#include "chem_jacobian.h"
#include <stdlib.h>
static void forcing(const double *v,double *f,void *context)
{
    size_t cells=*(size_t *)context;
    for (size_t c=0;c<cells;++c) { double x=v[c],y=v[cells+c];f[c]+=x*y;f[cells+c]+=x*x-y; }
}
static void verify(int which)
{
    size_t cells=which==1?3:1;double values[]={3,5,0,0,0,0},fd[12];
    if (which==1) { double input[]={1,.5,10,2,4,.1};memcpy(values,input,sizeof(input)); }
    if (which==6) { values[0]=0;values[1]=1; }
    CHECK(chem_finite_difference(2,cells,forcing,&cells,values,1e-8,fd)==CHEM_OK);
    if (which<=1 || which==6) for (size_t c=0;c<cells;++c) {
        NEAR(fd[c],values[cells+c],1e-5);NEAR(fd[cells+c],values[c],1e-5);
        NEAR(fd[2*cells+c],2*values[c],1e-5);NEAR(fd[3*cells+c],-1,1e-5);
    } else {
        unsigned char pattern[]={1,1,which==5?0:1,1};ChemSparse *sparse=NULL;
        CHECK(chem_sparse_create(2,pattern,&sparse)==CHEM_OK);
        double *analytic=calloc(chem_sparse_elements(sparse),sizeof(double));CHECK(analytic!=NULL);
        double expected[]={which==3?-999:-5,-3,-6,1};
        for (size_t i=0;i<2;++i) for (size_t j=0;j<2;++j)
            if (pattern[2*i+j]) analytic[chem_sparse_index(sparse,i,j)]=expected[2*i+j];
        if (which<=3) {
            ChemJacobianComparison result=chem_compare_jacobian(sparse,1,analytic,fd,1e-7,1e-7);
            CHECK(result.passed==(which==2));
            if (which==2) CHECK(result.max_absolute_error<1e-5);
            else { CHECK(result.worst_row==0);CHECK(result.worst_column==0); }
        } else {
            ChemJacobianComparison result=chem_check_sparsity(sparse,1,fd,1e-6);CHECK(result.passed==(which==4));
            if (which==5) { CHECK(result.worst_row==1);CHECK(result.worst_column==0); }
        }
        free(analytic);chem_sparse_free(sparse);
    }
}
void test_jacobian_verification(void)
{
#define J(name,which) RUN("jacobian_verification:JacobianVerification." name,verify(which))
    J("FiniteDifferenceMatchesAnalytical",0);J("MultiBlockFiniteDifference",1);J("ComparePassesForCorrectJacobian",2);
    J("CompareFailsForWrongJacobian",3);J("SparsityCompletenessPassesWhenComplete",4);
    J("SparsityCompletenessFailsWhenMissingEntry",5);J("NearZeroVariableHandled",6);
#undef J
}
