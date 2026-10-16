#include "test.h"
#include "chem_factor.h"

#include "test_random.h"
static double entry(const chem_sparse *s,size_t b,size_t i,size_t j)
{double x=0;if(s->map[i*s->n+j]!=CHEM_NONE) CHECK(chem_sparse_get(s,b,i,j,&x)==CHEM_OK);return x;}
static void ulp_equal(double a,double b)
{double low=b,high=b;for(int i=0;i<4;++i){low=nextafter(low,-INFINITY);high=nextafter(high,INFINITY);}CHECK(a>=low && a<=high);}
static void linear_case(size_t lanes,size_t cells,chem_lu_algorithm algorithm,int inplace,int solve,int kind,double initial)
{
    size_t n=kind==0 ? 3 : kind==2 ? 6 : kind==3 && solve ? 30 : kind==3 && inplace ? 5 : 10;
    random_state bits,numbers;int mt=kind==3 && solve;
    random_init(&bits,mt ? 12345u : 1u,mt);random_init(&numbers,mt ? 12345u : 1u,mt);
    chem_sparse_builder builder;REQUIRE(chem_sparse_builder_create(&builder,n)==CHEM_OK);builder.blocks=cells;builder.initial=initial;
    for(size_t i=0;i<n;++i)for(size_t j=0;j<n;++j)
        if(kind==0 || i==j || (kind!=2 && coin(&bits))) CHECK(chem_sparse_builder_add(&builder,i,j)==CHEM_OK);
    chem_sparse a;REQUIRE(chem_sparse_create(&a,&builder,CHEM_CSR,lanes,CHEM_DOUBLE)==CHEM_OK);chem_sparse_builder_destroy(&builder);
    const double dense[]={2,-1,-2,-4,6,3,-4,-2,8};
    for(size_t i=0;i<n;++i)for(size_t j=0;j<n;++j)if(a.map[i*n+j]!=CHEM_NONE) {
        double first=0;
        for(size_t b=0;b<cells;++b) {
            double x=kind==0 ? dense[i*n+j] : b && !solve && (kind==3 || (inplace && kind==1)) ? first : lognormal(&numbers,kind==2 ? 4 : 2);
            if(!b) first=x;
            CHECK(chem_sparse_set(&a,b,i,j,x)==CHEM_OK);
        }
    }
    chem_array rhs={0},solution={0};
    if(solve) {
        REQUIRE(chem_array_create(&rhs,cells,n,lanes,CHEM_DOUBLE,0)==CHEM_OK);
        for(size_t i=0;i<n;++i)for(size_t b=0;b<cells;++b) CHECK(chem_array_set(&rhs,b,i,kind==0 ? (double[]){23,42,9}[i] : lognormal(&numbers,kind==2 ? 4 : 2))==CHEM_OK);
        REQUIRE(chem_array_clone(&solution,&rhs)==CHEM_OK);
    }
    chem_lu_plan plan;REQUIRE(chem_sparse_lu_create(&plan,&a,algorithm)==CHEM_OK);chem_sparse l={0},u={0};
    if(inplace) {
        REQUIRE(chem_sparse_lu_matrix(&plan,&a,&l,initial)==CHEM_OK);
        for(size_t i=0;i<n;++i)for(size_t j=0;j<n;++j)if(l.map[i*n+j]!=CHEM_NONE)
            for(size_t b=0;b<cells;++b) CHECK(chem_sparse_set(&l,b,i,j,entry(&a,b,i,j))==CHEM_OK);
        REQUIRE(chem_sparse_lu_factor_inplace(&plan,&l)==CHEM_OK);
    } else {
        REQUIRE(chem_sparse_lu_matrices(&plan,&a,&l,&u,initial)==CHEM_OK);
        REQUIRE(chem_sparse_lu_factor(&plan,&a,&l,&u)==CHEM_OK);
    }
    CHECK(l.values.rows==a.values.rows);if(!inplace)CHECK(u.values.rows==a.values.rows);
    if(solve) {
        REQUIRE((inplace ? chem_sparse_lu_solve_inplace(&plan,&l,&solution) : chem_sparse_lu_solve(&plan,&l,&u,&solution))==CHEM_OK);
        CHECK(solution.rows==a.values.rows);CHECK(rhs.rows==a.values.rows);
        for(size_t b=0;b<cells;++b)for(size_t i=0;i<n;++i) {
            double sum=0;for(size_t j=0;j<n;++j) if(a.map[i*n+j]!=CHEM_NONE)sum+=entry(&a,b,i,j)*chem_array_get(&solution,b,j);
            if(inplace && kind==2) ulp_equal(sum,chem_array_get(&rhs,b,i));
            else NEAR(sum,chem_array_get(&rhs,b,i),kind==0 || kind==2 ? 1e-5 : kind==3 ? 2e-6 : 1e-6);
        }
    } else {
        for(size_t b=0;b<cells;++b)for(size_t i=0;i<n;++i)for(size_t j=0;j<n;++j) {
            double sum=0;
            if(inplace) {
                for(size_t k=0;k<=j;++k) {
                    if(i==k && l.map[k*n+j]!=CHEM_NONE)sum+=entry(&l,b,k,j);
                    else if(k<i && l.map[i*n+k]!=CHEM_NONE && l.map[k*n+j]!=CHEM_NONE)sum+=entry(&l,b,i,k)*entry(&l,b,k,j);
                }
            } else {
                CHECK(i>=j || l.map[i*n+j]==CHEM_NONE);CHECK(j>=i || u.map[i*n+j]==CHEM_NONE);
                for(size_t k=0;k<n;++k)if(l.map[i*n+k]!=CHEM_NONE && u.map[k*n+j]!=CHEM_NONE)sum+=entry(&l,b,i,k)*entry(&u,b,k,j);
            }
            NEAR(sum,entry(&a,b,i,j),kind==0 || kind==2 ? 1e-10 : inplace ? 1e-12 : 1e-9);
        }
    }
    chem_array_destroy(&rhs);chem_array_destroy(&solution);chem_sparse_destroy(&a);chem_sparse_destroy(&l);chem_sparse_destroy(&u);chem_lu_destroy(&plan);
}
static void extreme(size_t lanes,size_t cells,chem_lu_algorithm algorithm,int inplace,int solve)
{const double values[]={-INFINITY,-1,0,1,INFINITY};for(size_t i=0;i<5;++i)linear_case(lanes,cells,algorithm,inplace,solve,3,values[i]);}
static void ordering(size_t lanes,int inplace)
{
    random_state bits;random_init(&bits,1,0);unsigned char pattern[2500],reordered[2500];size_t permutation[50],sum=0,total=0;unsigned char seen[50]={0};
    chem_array original;REQUIRE(chem_array_create(&original,50,50,lanes,CHEM_INT,0)==CHEM_OK);
    for(size_t i=0;i<50;++i)for(size_t j=0;j<50;++j){pattern[i*50+j]=(i==j || coin(&bits));CHECK(chem_array_set(&original,i,j,pattern[i*50+j])==CHEM_OK);}
    REQUIRE(chem_markowitz(50,pattern,permutation)==CHEM_OK);
    for(size_t i=0;i<50;++i){CHECK(permutation[i]<50);CHECK(!seen[permutation[i]]);seen[permutation[i]]=1;sum+=i;total+=permutation[i];}
    CHECK(sum==total);
    for(size_t i=0;i<50;++i)for(size_t j=0;j<50;++j)reordered[i*50+j]=(unsigned char)chem_array_get(&original,permutation[i],permutation[j]);
    chem_lu_plan a,b;REQUIRE(chem_lu_create(&a,50,pattern,CHEM_LU_MOZART)==CHEM_OK);REQUIRE(chem_lu_create(&b,50,reordered,CHEM_LU_MOZART)==CHEM_OK);
    CHECK(a.nnz+(inplace ? 0 : 50)>b.nnz+(inplace ? 0 : 50));chem_lu_destroy(&a);chem_lu_destroy(&b);chem_array_destroy(&original);
}
TEST(linear_solver,LinearSolver,DenseMatrixStandardOrdering) { linear_case(0,1,CHEM_LU_MOZART,0,1,0,0); }
TEST(linear_solver,LinearSolver,RandomMatrixStandardOrdering) { linear_case(0,5,CHEM_LU_MOZART,0,1,1,0); }
TEST(linear_solver,LinearSolver,DiagonalMatrixStandardOrdering) { linear_case(0,5,CHEM_LU_MOZART,0,1,2,0); }
TEST(linear_solver,LinearSolver,DiagonalMarkowitzReorder) { ordering(0,0); }
TEST(linear_solver,LinearSolver,StandardOrderingAgnosticToInitialValue) { extreme(0,5,CHEM_LU_MOZART,0,1); }
TEST(linear_solver,LinearSolver,DenseMatrixVectorOrdering) { linear_case(1,1,CHEM_LU_MOZART,0,1,0,0); linear_case(2,1,CHEM_LU_MOZART,0,1,0,0); linear_case(3,1,CHEM_LU_MOZART,0,1,0,0); linear_case(4,1,CHEM_LU_MOZART,0,1,0,0); }
TEST(linear_solver,LinearSolver,RandomMatrixVectorOrdering) { linear_case(1,5,CHEM_LU_MOZART,0,1,1,0); linear_case(2,5,CHEM_LU_MOZART,0,1,1,0); linear_case(3,5,CHEM_LU_MOZART,0,1,1,0); linear_case(4,5,CHEM_LU_MOZART,0,1,1,0); }
TEST(linear_solver,LinearSolver,VectorOrderingAgnosticToInitialValue) { extreme(1,1,CHEM_LU_MOZART,0,1); extreme(2,2,CHEM_LU_MOZART,0,1); extreme(3,5,CHEM_LU_MOZART,0,1); extreme(4,5,CHEM_LU_MOZART,0,1); }
TEST(linear_solver,LinearSolver,DiagonalMatrixVectorOrdering) { linear_case(1,5,CHEM_LU_MOZART,0,1,2,0); linear_case(2,5,CHEM_LU_MOZART,0,1,2,0); linear_case(3,5,CHEM_LU_MOZART,0,1,2,0); linear_case(4,5,CHEM_LU_MOZART,0,1,2,0); }
TEST(linear_solver,LinearSolver,VectorDiagonalMarkowitzReordering) { ordering(1,0); ordering(2,0); ordering(3,0); ordering(4,0); }
TEST(linear_solver_in_place,LinearSolverInPlace,DenseMatrixStandardOrdering) { linear_case(0,1,CHEM_LU_MOZART,1,1,0,0); }
TEST(linear_solver_in_place,LinearSolverInPlace,RandomMatrixStandardOrdering) { linear_case(0,5,CHEM_LU_MOZART,1,1,1,0); }
TEST(linear_solver_in_place,LinearSolverInPlace,DiagonalMatrixStandardOrdering) { linear_case(0,5,CHEM_LU_MOZART,1,1,2,0); }
TEST(linear_solver_in_place,LinearSolverInPlace,DiagonalMarkowitzReorder) { ordering(0,1); }
TEST(linear_solver_in_place,LinearSolverInPlace,StandardOrderingAgnosticToInitialValue) { extreme(0,5,CHEM_LU_MOZART,1,1); }
TEST(linear_solver_in_place,LinearSolverInPlace,DenseMatrixVectorOrdering) { linear_case(1,1,CHEM_LU_MOZART,1,1,0,0); linear_case(2,1,CHEM_LU_MOZART,1,1,0,0); linear_case(3,1,CHEM_LU_MOZART,1,1,0,0); linear_case(4,1,CHEM_LU_MOZART,1,1,0,0); }
TEST(linear_solver_in_place,LinearSolverInPlace,RandomMatrixVectorOrdering) { linear_case(1,5,CHEM_LU_MOZART,1,1,1,0); linear_case(2,5,CHEM_LU_MOZART,1,1,1,0); linear_case(3,5,CHEM_LU_MOZART,1,1,1,0); linear_case(4,5,CHEM_LU_MOZART,1,1,1,0); }
TEST(linear_solver_in_place,LinearSolverInPlace,VectorOrderingAgnosticToInitialValue) { extreme(1,1,CHEM_LU_MOZART,1,1); extreme(2,2,CHEM_LU_MOZART,1,1); extreme(3,5,CHEM_LU_MOZART,1,1); extreme(4,5,CHEM_LU_MOZART,1,1); }
TEST(linear_solver_in_place,LinearSolverInPlace,DiagonalMatrixVectorOrdering) { linear_case(1,5,CHEM_LU_MOZART,1,1,2,0); linear_case(2,5,CHEM_LU_MOZART,1,1,2,0); linear_case(3,5,CHEM_LU_MOZART,1,1,2,0); linear_case(4,5,CHEM_LU_MOZART,1,1,2,0); }
TEST(linear_solver_in_place,LinearSolverInPlace,VectorDiagonalMarkowitzReordering) { ordering(1,1); ordering(2,1); ordering(3,1); ordering(4,1); }
TEST(lu_decomposition_doolittle,LuDecompositionDoolittle,DenseMatrixStandardOrdering) { linear_case(0,1,CHEM_LU_DOOLITTLE,0,0,0,0); }
TEST(lu_decomposition_doolittle,LuDecompositionDoolittle,RandomMatrixStandardOrdering) { linear_case(0,1,CHEM_LU_DOOLITTLE,0,0,1,0); linear_case(0,5,CHEM_LU_DOOLITTLE,0,0,1,0); }
TEST(lu_decomposition_doolittle,LuDecompositionDoolittle,DiagonalMatrixStandardOrdering) { linear_case(0,5,CHEM_LU_DOOLITTLE,0,0,2,0); }
TEST(lu_decomposition_doolittle,LuDecompositionDoolittle,AgnosticToInitialValueStandardOrdering) { extreme(0,5,CHEM_LU_DOOLITTLE,0,0); }
TEST(lu_decomposition_doolittle,LuDecompositionDoolittle,DenseMatrixVectorOrdering) { linear_case(1,1,CHEM_LU_DOOLITTLE,0,0,0,0); linear_case(2,1,CHEM_LU_DOOLITTLE,0,0,0,0); linear_case(3,1,CHEM_LU_DOOLITTLE,0,0,0,0); linear_case(4,1,CHEM_LU_DOOLITTLE,0,0,0,0); }
TEST(lu_decomposition_doolittle,LuDecompositionDoolittle,RandomMatrixVectorOrdering) { linear_case(1,5,CHEM_LU_DOOLITTLE,0,0,1,0); linear_case(2,5,CHEM_LU_DOOLITTLE,0,0,1,0); linear_case(3,5,CHEM_LU_DOOLITTLE,0,0,1,0); linear_case(4,5,CHEM_LU_DOOLITTLE,0,0,1,0); }
TEST(lu_decomposition_doolittle,LuDecompositionDoolittle,DiagonalMatrixVectorOrdering) { linear_case(1,5,CHEM_LU_DOOLITTLE,0,0,2,0); linear_case(2,5,CHEM_LU_DOOLITTLE,0,0,2,0); linear_case(3,5,CHEM_LU_DOOLITTLE,0,0,2,0); linear_case(4,5,CHEM_LU_DOOLITTLE,0,0,2,0); }
TEST(lu_decomposition_doolittle,LuDecompositionDoolittle,AgnosticToInitialValueVectorOrdering) { extreme(1,5,CHEM_LU_DOOLITTLE,0,0); extreme(2,5,CHEM_LU_DOOLITTLE,0,0); extreme(3,5,CHEM_LU_DOOLITTLE,0,0); extreme(4,5,CHEM_LU_DOOLITTLE,0,0); }
TEST(lu_decomposition_mozart,LuDecompositionMozart,DenseMatrixStandardOrdering) { linear_case(0,1,CHEM_LU_MOZART,0,0,0,0); }
TEST(lu_decomposition_mozart,LuDecompositionMozart,RandomMatrixStandardOrdering) { linear_case(0,1,CHEM_LU_MOZART,0,0,1,0); linear_case(0,5,CHEM_LU_MOZART,0,0,1,0); }
TEST(lu_decomposition_mozart,LuDecompositionMozart,DiagonalMatrixStandardOrdering) { linear_case(0,5,CHEM_LU_MOZART,0,0,2,0); }
TEST(lu_decomposition_mozart,LuDecompositionMozart,AgnosticToInitialValueStandardOrdering) { extreme(0,5,CHEM_LU_MOZART,0,0); }
TEST(lu_decomposition_mozart,LuDecompositionMozart,DenseMatrixVectorOrdering) { linear_case(1,1,CHEM_LU_MOZART,0,0,0,0); linear_case(2,1,CHEM_LU_MOZART,0,0,0,0); linear_case(3,1,CHEM_LU_MOZART,0,0,0,0); linear_case(4,1,CHEM_LU_MOZART,0,0,0,0); }
TEST(lu_decomposition_mozart,LuDecompositionMozart,RandomMatrixVectorOrdering) { linear_case(1,5,CHEM_LU_MOZART,0,0,1,0); linear_case(2,5,CHEM_LU_MOZART,0,0,1,0); linear_case(3,5,CHEM_LU_MOZART,0,0,1,0); linear_case(4,5,CHEM_LU_MOZART,0,0,1,0); }
TEST(lu_decomposition_mozart,LuDecompositionMozart,DiagonalMatrixVectorOrdering) { linear_case(1,5,CHEM_LU_MOZART,0,0,2,0); linear_case(2,5,CHEM_LU_MOZART,0,0,2,0); linear_case(3,5,CHEM_LU_MOZART,0,0,2,0); linear_case(4,5,CHEM_LU_MOZART,0,0,2,0); }
TEST(lu_decomposition_mozart,LuDecompositionMozart,AgnosticToInitialValueVectorOrdering) { extreme(1,5,CHEM_LU_MOZART,0,0); extreme(2,5,CHEM_LU_MOZART,0,0); extreme(3,5,CHEM_LU_MOZART,0,0); extreme(4,5,CHEM_LU_MOZART,0,0); }
TEST(lu_decomposition_doolittle_in_place,LuDecompositionDoolittleInPlace,DenseMatrixStandardOrdering) { linear_case(0,1,CHEM_LU_DOOLITTLE,1,0,0,0); }
TEST(lu_decomposition_doolittle_in_place,LuDecompositionDoolittleInPlace,RandomMatrixStandardOrdering) { linear_case(0,1,CHEM_LU_DOOLITTLE,1,0,1,0); linear_case(0,5,CHEM_LU_DOOLITTLE,1,0,1,0); }
TEST(lu_decomposition_doolittle_in_place,LuDecompositionDoolittleInPlace,DiagonalMatrixStandardOrdering) { linear_case(0,5,CHEM_LU_DOOLITTLE,1,0,2,0); }
TEST(lu_decomposition_doolittle_in_place,LuDecompositionDoolittleInPlace,AgnosticToInitialValueStandardOrdering) { extreme(0,5,CHEM_LU_DOOLITTLE,1,0); }
TEST(lu_decomposition_doolittle_in_place,LuDecompositionDoolittleInPlace,DenseMatrixVectorOrdering) { linear_case(1,1,CHEM_LU_DOOLITTLE,1,0,0,0); linear_case(2,1,CHEM_LU_DOOLITTLE,1,0,0,0); linear_case(3,1,CHEM_LU_DOOLITTLE,1,0,0,0); linear_case(4,1,CHEM_LU_DOOLITTLE,1,0,0,0); }
TEST(lu_decomposition_doolittle_in_place,LuDecompositionDoolittleInPlace,RandomMatrixVectorOrdering) { linear_case(1,5,CHEM_LU_DOOLITTLE,1,0,1,0); linear_case(2,5,CHEM_LU_DOOLITTLE,1,0,1,0); linear_case(3,5,CHEM_LU_DOOLITTLE,1,0,1,0); linear_case(4,5,CHEM_LU_DOOLITTLE,1,0,1,0); }
TEST(lu_decomposition_doolittle_in_place,LuDecompositionDoolittleInPlace,DiagonalMatrixVectorOrdering) { linear_case(1,5,CHEM_LU_DOOLITTLE,1,0,2,0); linear_case(2,5,CHEM_LU_DOOLITTLE,1,0,2,0); linear_case(3,5,CHEM_LU_DOOLITTLE,1,0,2,0); linear_case(4,5,CHEM_LU_DOOLITTLE,1,0,2,0); }
TEST(lu_decomposition_doolittle_in_place,LuDecompositionDoolittleInPlace,AgnosticToInitialValueVectorOrdering) { extreme(1,5,CHEM_LU_DOOLITTLE,1,0); extreme(2,5,CHEM_LU_DOOLITTLE,1,0); extreme(3,5,CHEM_LU_DOOLITTLE,1,0); extreme(4,5,CHEM_LU_DOOLITTLE,1,0); }
TEST(lu_decomposition_mozart_in_place,LuDecompositionMozartInPlace,DenseMatrixStandardOrdering) { linear_case(0,1,CHEM_LU_MOZART,1,0,0,0); }
TEST(lu_decomposition_mozart_in_place,LuDecompositionMozartInPlace,RandomMatrixStandardOrdering) { linear_case(0,1,CHEM_LU_MOZART,1,0,1,0); linear_case(0,5,CHEM_LU_MOZART,1,0,1,0); }
TEST(lu_decomposition_mozart_in_place,LuDecompositionMozartInPlace,DiagonalMatrixStandardOrdering) { linear_case(0,5,CHEM_LU_MOZART,1,0,2,0); }
TEST(lu_decomposition_mozart_in_place,LuDecompositionMozartInPlace,AgnosticToInitialValueStandardOrdering) { extreme(0,5,CHEM_LU_MOZART,1,0); }
TEST(lu_decomposition_mozart_in_place,LuDecompositionMozartInPlace,DenseMatrixVectorOrdering) { linear_case(1,1,CHEM_LU_MOZART,1,0,0,0); linear_case(2,1,CHEM_LU_MOZART,1,0,0,0); linear_case(3,1,CHEM_LU_MOZART,1,0,0,0); linear_case(4,1,CHEM_LU_MOZART,1,0,0,0); }
TEST(lu_decomposition_mozart_in_place,LuDecompositionMozartInPlace,RandomMatrixVectorOrdering) { linear_case(1,5,CHEM_LU_MOZART,1,0,1,0); linear_case(2,5,CHEM_LU_MOZART,1,0,1,0); linear_case(3,5,CHEM_LU_MOZART,1,0,1,0); linear_case(4,5,CHEM_LU_MOZART,1,0,1,0); }
TEST(lu_decomposition_mozart_in_place,LuDecompositionMozartInPlace,DiagonalMatrixVectorOrdering) { linear_case(1,5,CHEM_LU_MOZART,1,0,2,0); linear_case(2,5,CHEM_LU_MOZART,1,0,2,0); linear_case(3,5,CHEM_LU_MOZART,1,0,2,0); linear_case(4,5,CHEM_LU_MOZART,1,0,2,0); }
TEST(lu_decomposition_mozart_in_place,LuDecompositionMozartInPlace,AgnosticToInitialValueVectorOrdering) { extreme(1,5,CHEM_LU_MOZART,1,0); extreme(2,5,CHEM_LU_MOZART,1,0); extreme(3,5,CHEM_LU_MOZART,1,0); extreme(4,5,CHEM_LU_MOZART,1,0); }
