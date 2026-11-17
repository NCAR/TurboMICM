#include "chem.h"

#include <math.h>
#include <stdio.h>

#define CHECK(condition) do { if (!(condition)) { \
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); return 1; \
} } while (0)

int main(void)
{
    enum { N=4, CELLS=3 };
    const double base[N][N]={{4,1,0,2},{2,5,1,0},{0,2,6,1},{1,0,3,7}};
    const size_t rows[]={0,4,8,12,16};
    int columns[16], lookup[16];
    const unsigned char algebraic[]={0,1,0,0};
    double matrix[16*CELLS], original[16*CELLS], rhs[N*CELLS];
    for (size_t i=0;i<N;++i) for (size_t j=0;j<N;++j) {
        size_t entry=i*N+j;
        columns[entry]=(int)j; lookup[entry]=(int)entry;
        for (size_t cell=0;cell<CELLS;++cell) {
            double value=base[i][j]+(i==j ? .1*cell : 0);
            matrix[entry*CELLS+cell]=value;
            original[entry*CELLS+cell]=value+(i==j && !algebraic[i] ? .7 : 0);
        }
    }
    ChemLinearSystem system={N,CELLS,16,rows,columns,lookup,algebraic,matrix};
    ChemLinearSolver solver=chem_doolittle_linear_solver();
    CHECK(solver.factor(&system,.7,NULL)==CHEM_OK);
    for (size_t cell=0;cell<CELLS;++cell) for (size_t i=0;i<N;++i) for (size_t j=0;j<N;++j) {
        double product=0;
        for (size_t k=0;k<N;++k) {
            double l=i==k ? 1 : i>k ? matrix[(i*N+k)*CELLS+cell] : 0;
            double u=k<=j ? matrix[(k*N+j)*CELLS+cell] : 0;
            product+=l*u;
        }
        CHECK(fabs(product-original[(i*N+j)*CELLS+cell])<1e-13);
    }
    for (size_t cell=0;cell<CELLS;++cell) for (size_t i=0;i<N;++i) {
        rhs[i*CELLS+cell]=0;
        for (size_t j=0;j<N;++j) rhs[i*CELLS+cell]+=original[(i*N+j)*CELLS+cell]*(j+1+cell);
    }
    CHECK(solver.solve(&system,rhs,NULL)==CHEM_OK);
    for (size_t cell=0;cell<CELLS;++cell) for (size_t i=0;i<N;++i)
        CHECK(fabs(rhs[i*CELLS+cell]-(i+1+cell))<1e-13);
    puts("Doolittle factors and multi-cell solves: PASS");
    return 0;
}
