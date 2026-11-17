#include "chem.h"
#include "sparse.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static uint32_t rng;
static double uniform(void)
{
    rng = rng * 1664525u + 1013904223u;
    return ((rng >> 8) + 0.5) / 16777216.0;
}
static double lognormal(double sigma)
{
    return exp(-2 + sigma * sqrt(-2 * log(uniform())) * cos(6.2831853071795864769 * uniform()));
}

/* Exercise every storage width and block via the integration solver's LU
 * kernels, checking the full L*U reconstruction, including structural zeros. */
static int run_case(int kind, int inplace, ChemLUAlgorithm algorithm, size_t width, size_t blocks, double initial)
{
    size_t n = kind == 0 ? 3 : kind == 2 ? 6 : kind == 3 && inplace ? 5 : 10;
    unsigned char pattern[100] = {0};
    double a[100], actual[100], l[100], u[100], first[100] = {0};
    const double dense[] = {2,-1,-2,-4,6,3,-4,-2,8};
    double tol = kind == 0 || kind == 2 ? 1e-10 : inplace ? 1e-12 : 1e-9;
    ChemSparseBuilder builder = {0}; ChemSparseMatrix matrix = {0};
    rng = 739391u;
    if (chem_sparse_builder_init(&builder,n)) return 0;
    for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) {
        int present = kind == 0 || (kind == 2 ? i==j : i==j || uniform() < .5);
        pattern[i*n+j] = present;
        if (present && chem_sparse_builder_add(&builder,i,j)) return 0;
    }
    if (chem_sparse_init(&matrix,&builder,blocks,width,0,initial)) return 0;
    chem_sparse_builder_free(&builder);
    rng = 15339u;
    for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) if (pattern[i*n+j]) {
        for (size_t b=0;b<blocks;++b) {
            double value = kind==0 ? dense[i*n+j] : ((inplace || kind==3) && b>0) ? first[i*n+j] : lognormal(kind==2 ? 4 : 2);
            if (!b) first[i*n+j] = value;
            *chem_sparse_at(&matrix,b,i,j) = value;
        }
    }
    ChemWorkspace *w=chem_linear_workspace_create_with_initial(n,pattern,initial);
    if (!w) { chem_sparse_free(&matrix); return 0; }
    int ok=1;
    for (size_t b=0;b<blocks && ok;++b) {
        for(size_t i=0;i<n;++i) for(size_t j=0;j<n;++j) a[i*n+j]=pattern[i*n+j] ? *chem_sparse_at(&matrix,b,i,j) : 0;
        memcpy(actual,a,n*n*sizeof(double));
        for(size_t i=0;i<n*n;++i) l[i]=u[i]=initial;
        if (inplace) {
            if (chem_linear_factor_in_place(w,actual,algorithm)) { ok=0; break; }
            for(size_t i=0;i<n;++i) for(size_t j=0;j<n;++j) {
                l[i*n+j]=i==j ? 1 : i>j ? actual[i*n+j] : 0;
                u[i*n+j]=i<=j ? actual[i*n+j] : 0;
            }
        } else if (chem_linear_factor(w,a,algorithm) || chem_linear_factors(w,l,u)) { ok=0; break; }
        for(size_t i=0;i<n;++i) for(size_t j=0;j<n;++j) {
            double value=0;
            for(size_t k=0;k<n;++k) value+=l[i*n+k]*u[k*n+j];
            if (!isfinite(value) || fabs(value-a[i*n+j]) > tol || (i<j && l[i*n+j]!=0) || (i>j && u[i*n+j]!=0)) ok=0;
        }
    }
    chem_workspace_free(w); chem_sparse_free(&matrix); return ok;
}

/* The linear solver cases check A*x against b, not the answer x itself. */
static int solve_case(int kind, int inplace, size_t width, size_t blocks, double initial)
{
    size_t n=kind==0 ? 3 : kind==2 ? 6 : kind==3 ? 30 : 10;
    unsigned char pattern[900]={0}; double a[900], factored[900], x[30], b[30];
    const double dense[]={2,-1,-2,-4,6,3,-4,-2,8}, dense_rhs[]={23,42,9};
    double tol=kind==0 || kind==2 ? 1e-5 : kind==3 ? 2e-6 : 1e-6;
    ChemSparseBuilder builder={0}; ChemSparseMatrix matrix={0}; ChemMatrix rhs={0};
    rng=12345u; if(chem_sparse_builder_init(&builder,n)) return 0;
    for(size_t i=0;i<n;++i) for(size_t j=0;j<n;++j) {
        pattern[i*n+j]=kind==0 || (kind==2 ? i==j : i==j || uniform()<.5);
        if(pattern[i*n+j] && chem_sparse_builder_add(&builder,i,j)) return 0;
    }
    if(chem_sparse_init(&matrix,&builder,blocks,width,0,0) || chem_matrix_init(&rhs,blocks,n,width,width!=1,0)) return 0;
    chem_sparse_builder_free(&builder); rng=52141u;
    for(size_t i=0;i<n;++i) for(size_t j=0;j<n;++j) if(pattern[i*n+j]) for(size_t k=0;k<blocks;++k)
        *chem_sparse_at(&matrix,k,i,j)=kind==0 ? dense[i*n+j] : lognormal(kind==2 ? 4 : 2);
    for(size_t i=0;i<n;++i) for(size_t k=0;k<blocks;++k)
        *chem_matrix_at(&rhs,k,i)=kind==0 ? dense_rhs[i] : lognormal(kind==2 ? 4 : 2);
    ChemWorkspace *w=chem_linear_workspace_create_with_initial(n,pattern,initial); int ok=w!=NULL;
    for(size_t k=0;k<blocks && ok;++k) {
        for(size_t i=0;i<n;++i) {
            b[i]=x[i]=*chem_matrix_at(&rhs,k,i);
            for(size_t j=0;j<n;++j) a[i*n+j]=pattern[i*n+j] ? *chem_sparse_at(&matrix,k,i,j) : 0;
        }
        memcpy(factored,a,n*n*sizeof(double));
        if((inplace ? chem_linear_factor_in_place(w,factored,CHEM_LU_MOZART) : chem_linear_factor(w,a,CHEM_LU_MOZART)) || chem_linear_solve(w,x)) { ok=0; break; }
        for(size_t i=0;i<n;++i) {
            double result=0; for(size_t j=0;j<n;++j) result+=a[i*n+j]*x[j];
            if(!isfinite(result) || fabs(result-b[i])>tol) ok=0;
        }
    }
    chem_workspace_free(w); chem_sparse_free(&matrix); chem_matrix_free(&rhs); return ok;
}

static int markowitz_case(size_t width)
{
    enum { N=50 }; unsigned char pattern[N*N], reordered[N*N]; size_t order[N];
    ChemIntMatrix input={0}; rng=985137u;
    if(chem_int_matrix_init(&input,N,N,width,0)) return 0;
    for(size_t i=0;i<N;++i) for(size_t j=0;j<N;++j) {
        pattern[i*N+j]=(i==j || uniform()<.5);
        *chem_int_matrix_at(&input,i,j)=pattern[i*N+j];
    }
    if(chem_diagonal_markowitz_order(N,pattern,order)) { chem_int_matrix_free(&input); return 0; }
    size_t sum=0; unsigned char seen[N]={0}; int ok=1;
    for(size_t i=0;i<N;++i) {
        if(order[i]>=N || seen[order[i]]) { ok=0; break; }
        seen[order[i]]=1; sum+=order[i];
        for(size_t j=0;j<N;++j) reordered[i*N+j]=*chem_int_matrix_at(&input,order[i],order[j]);
    }
    if(sum!=N*(N-1)/2) ok=0;
    ChemWorkspace *a=chem_linear_workspace_create(N,pattern), *b=chem_linear_workspace_create(N,reordered);
    if(!a || !b || chem_linear_symbolic_nonzeros(a)<=chem_linear_symbolic_nonzeros(b)) ok=0;
    chem_workspace_free(a); chem_workspace_free(b); chem_int_matrix_free(&input); return ok;
}

int main(void)
{
    const char *prefixes[]={"lu_decomposition_doolittle:LuDecompositionDoolittle", "lu_decomposition_doolittle_in_place:LuDecompositionDoolittleInPlace", "lu_decomposition_mozart:LuDecompositionMozart", "lu_decomposition_mozart_in_place:LuDecompositionMozartInPlace"};
    const char *names[]={"DenseMatrix", "RandomMatrix", "DiagonalMatrix", "AgnosticToInitialValue"};
    const double values[]={-INFINITY,-1,0,1,INFINITY};
    for(int impl=0;impl<4;++impl) for(int vector=0;vector<2;++vector) for(int kind=0;kind<4;++kind) {
        int ok=1, inplace=impl%2; ChemLUAlgorithm alg=impl<2 ? CHEM_LU_DOOLITTLE : CHEM_LU_MOZART;
        for(size_t width=1;width<=(vector?4u:1u);++width) {
            int count=kind==3?5:1;
            for(int i=0;i<count;++i) {
                size_t blocks=kind==0 ? 1 : 5;
                ok &= run_case(kind,inplace,alg,width,blocks,values[kind==3?i:2]);
                if (kind==1 && !vector && !inplace) ok &= run_case(kind,inplace,alg,width,1,0);
            }
        }
        printf("%s.%s%sOrdering %d\n",prefixes[impl],names[kind],vector?"Vector":"Standard",ok);
    }
    const char *linear_prefix[]={"linear_solver:LinearSolver", "linear_solver_in_place:LinearSolverInPlace"};
    const char *linear_names[]={"DenseMatrix", "RandomMatrix", "DiagonalMatrix", "AgnosticToInitialValue"};
    for(int inplace=0;inplace<2;++inplace) for(int vector=0;vector<2;++vector) for(int kind=0;kind<4;++kind) {
        int ok=1;
        for(size_t width=1;width<=(vector?4u:1u);++width) for(int i=0;i<(kind==3?5:1);++i) {
            size_t blocks=kind==0 ? 1 : kind==3 && vector && width<3 ? width : 5;
            ok &= solve_case(kind,inplace,width,blocks,values[kind==3?i:2]);
        }
        if(kind==3) printf("%s.%sOrdering%s %d\n",linear_prefix[inplace],vector?"Vector":"Standard",linear_names[kind],ok);
        else printf("%s.%s%sOrdering %d\n",linear_prefix[inplace],linear_names[kind],vector?"Vector":"Standard",ok);
    }
    for(int inplace=0;inplace<2;++inplace) {
        int standard=markowitz_case(1), vector=1;
        for(size_t width=1;width<=4;++width) vector &= markowitz_case(width);
        printf("%s.DiagonalMarkowitzReorder %d\n",linear_prefix[inplace],standard);
        printf("%s.VectorDiagonalMarkowitzReordering %d\n",linear_prefix[inplace],vector);
    }
    return 0;
}
