/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Port of the in-place LU and linear solver policy tests. */
#include "test.h"
#include "chem_solver.h"
#include <stdint.h>
#include <stdlib.h>

/* Fixed portable uniform and lognormal draws replace implementation-defined
 * C++ random engines. Every matrix dimension, distribution and bound is kept. */
static uint64_t rng;
static double uniform(void)
{
    rng=(rng*48271)%2147483647;
    return (double)rng/2147483647;
}
static double lognormal(double sigma)
{
    double z=sqrt(-2*log(uniform()))*cos(2*CHEM_PI*uniform());
    return exp(-2+sigma*z);
}
static void matrix_test(int kind,int linear,ChemLuAlgorithm algorithm,size_t blocks,size_t group,double initial)
{
    size_t n=kind==0?3:kind==2?6:kind==3?(linear?30:5):10;
    unsigned char *pattern=calloc(n*n,1);
    double *original=calloc(n*n*blocks,sizeof(double)),*rhs=calloc(n*blocks,sizeof(double));
    CHECK(pattern && original && rhs);
    if (!pattern || !original || !rhs) goto done;
    rng=kind==3 && linear?12345:1;
    for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j)
        pattern[i*n+j]=kind==0 || i==j || (kind!=2 && uniform()<.5);
    const double dense[]={2,-1,-2,-4,6,3,-4,-2,8}, b[]={23,42,9};
    rng=kind==3 && linear?12345:1;
    for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) if (pattern[i*n+j]) {
        double shared=kind==0?dense[i*n+j]:lognormal(kind==2?4:2);
        for (size_t c=0;c<blocks;++c)
            original[(i*n+j)*blocks+c]=linear && kind!=0?lognormal(kind==2?4:2):shared;
    }
    for (size_t i=0;i<n;++i) for (size_t c=0;c<blocks;++c) rhs[i*blocks+c]=kind==0?b[i]:lognormal(kind==2?4:2);
    ChemSparse *a=NULL;
    CHECK(chem_sparse_create(n,pattern,&a)==CHEM_OK);
    if (!a) goto done;
    CHECK(chem_sparse_order(a)==n);
    for (size_t start=0;start<blocks;start+=group) {
        size_t cells=blocks-start<group?blocks-start:group;
        double *lu=malloc(chem_sparse_elements(a)*cells*sizeof(double)),*x=malloc(n*cells*sizeof(double));
        CHECK(lu && x);
        if (!lu || !x) { free(lu);free(x);break; }
        for (size_t i=0;i<chem_sparse_elements(a)*cells;++i) lu[i]=initial;
        for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) {
            size_t k=chem_sparse_index(a,i,j);
            if (k==SIZE_MAX) continue;
            for (size_t c=0;c<cells;++c) lu[k*cells+c]=original[(i*n+j)*blocks+start+c];
        }
        CHECK(chem_lu_factor(a,cells,lu,algorithm)==CHEM_OK);
        if (linear) {
            for (size_t i=0;i<n;++i) for (size_t c=0;c<cells;++c) x[i*cells+c]=rhs[i*blocks+start+c];
            chem_lu_solve(a,cells,lu,x);
            for (size_t c=0;c<cells;++c) for (size_t i=0;i<n;++i) {
                double result=0;
                for (size_t j=0;j<n;++j) result+=original[(i*n+j)*blocks+start+c]*x[j*cells+c];
                double expected=rhs[i*blocks+start+c];
                /* GoogleTest double equality allows four ULPs; two is tighter. */
                double tol=kind==0?1e-5:kind==1?1e-6:kind==3?2e-6:2*fabs(nextafter(expected,INFINITY)-expected);
                NEAR(result,expected,tol);
            }
        } else {
            for (size_t c=0;c<cells;++c) for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) {
                double sum=0;
                for (size_t k=0;k<=j;++k) {
                    size_t ukj=chem_sparse_index(a,k,j),lik=chem_sparse_index(a,i,k);
                    if (k==i && ukj!=SIZE_MAX) sum+=lu[ukj*cells+c];
                    else if (k<i && lik!=SIZE_MAX && ukj!=SIZE_MAX) sum+=lu[lik*cells+c]*lu[ukj*cells+c];
                }
                NEAR(sum,original[(i*n+j)*blocks+start+c],kind==0 || kind==2?1e-10:1e-12);
            }
        }
        free(lu);free(x);
    }
    chem_sparse_free(a);
done:
    free(pattern);free(original);free(rhs);
}
static void matrix_cases(int kind,int linear,ChemLuAlgorithm algorithm,int vector)
{
    double initials[]={-INFINITY,-1,0,1,INFINITY};
    size_t max_group=vector?4:5;
    for (size_t group=vector?1:5;group<=max_group;++group) {
        if (kind==3) for (size_t v=0;v<5;++v) {
            size_t blocks=linear && vector && group<=2?group:5;
            matrix_test(kind,linear,algorithm,blocks,group,initials[v]);
        } else {
            matrix_test(kind,linear,algorithm,kind==0?1:5,group,0);
            if (!linear && kind==1 && !vector) matrix_test(kind,linear,algorithm,1,1,0);
        }
    }
}
static void markowitz_test(void)
{
    unsigned char pattern[2500],reordered[2500];size_t perm[50];
    rng=1;
    for (size_t i=0;i<50;++i) for (size_t j=0;j<50;++j) pattern[i*50+j]=i==j || uniform()<.5;
    CHECK(chem_markowitz(50,pattern,perm)==CHEM_OK);
    size_t sum=0;unsigned char seen[50]={0};
    for (size_t i=0;i<50;++i) { CHECK(perm[i]<50);CHECK(!seen[perm[i]]);seen[perm[i]]=1;sum+=perm[i]; }
    CHECK(sum==50*49/2);
    for (size_t i=0;i<50;++i) for (size_t j=0;j<50;++j) reordered[i*50+j]=pattern[perm[i]*50+perm[j]];
    ChemSparse *a=NULL,*b=NULL;
    CHECK(chem_sparse_create(50,pattern,&a)==CHEM_OK);CHECK(chem_sparse_create(50,reordered,&b)==CHEM_OK);
    if (a && b) CHECK(chem_sparse_elements(a)>chem_sparse_elements(b));
    chem_sparse_free(a);chem_sparse_free(b);
}
void test_sparse(void)
{
    const char *suffix[]={"DenseMatrixStandardOrdering","RandomMatrixStandardOrdering","DiagonalMatrixStandardOrdering","AgnosticToInitialValueStandardOrdering",
        "DenseMatrixVectorOrdering","RandomMatrixVectorOrdering","DiagonalMatrixVectorOrdering","AgnosticToInitialValueVectorOrdering"};
    char id[200];
    for (int algorithm=0;algorithm<2;++algorithm) for (int test=0;test<8;++test) {
        snprintf(id,sizeof(id),"lu_decomposition_%s_in_place:LuDecomposition%sInPlace.%s",algorithm?"doolittle":"mozart",algorithm?"Doolittle":"Mozart",suffix[test]);
        RUN(id,matrix_cases(test%4,0,(ChemLuAlgorithm)algorithm,test>=4));
    }
    const char *linear_suffix[]={"DenseMatrixStandardOrdering","RandomMatrixStandardOrdering","DiagonalMatrixStandardOrdering","StandardOrderingAgnosticToInitialValue",
        "DenseMatrixVectorOrdering","RandomMatrixVectorOrdering","DiagonalMatrixVectorOrdering","VectorOrderingAgnosticToInitialValue"};
    for (int test=0;test<8;++test) {
        snprintf(id,sizeof(id),"linear_solver_in_place:LinearSolverInPlace.%s",linear_suffix[test]);
        RUN(id,matrix_cases(test%4,1,CHEM_LU_MOZART,test>=4));
    }
    RUN("linear_solver_in_place:LinearSolverInPlace.DiagonalMarkowitzReorder",markowitz_test());
    RUN("linear_solver_in_place:LinearSolverInPlace.VectorDiagonalMarkowitzReordering",markowitz_test());
}
