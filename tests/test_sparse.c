#include "sparse.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void)
{
    /* A cycle creates fill during elimination even though the input is sparse. */
    const size_t reactants[] = {0,1,2,3};
    const ChemProduct products[] = {{1,1},{2,1},{3,1},{0,1}};
    const ChemReaction reactions[] = {{0,1,0,1},{1,1,1,1},{2,1,2,1},{3,1,3,1}};
    const ChemMechanism m = {4,4,reactions,reactants,products};
    ChemSparsePlan *p = chem_sparse_plan_create(&m);
    assert(p && chem_sparse_plan_size(p) == 4);
    double markers[16];
    const unsigned char filled[]={1,0,0,1, 1,1,0,1, 0,1,1,1, 0,0,1,1};
    for (size_t i=0;i<16;++i) markers[i]=(double)(i+1);
    chem_sparse_negate(p,markers);
    for (size_t i=0;i<16;++i) assert(markers[i]==(filled[i] ? -(double)(i+1) : (double)(i+1)));
    chem_sparse_clear(p,markers);
    for (size_t i=0;i<16;++i) assert(markers[i]==(filled[i] ? 0 : (double)(i+1)));
    double a[] = {4,0,0,2, 3,5,0,0, 0,2,6,0, 0,0,1,7};
    double b[] = {12,13,22,31}; /* A * [1,2,3,4] */
    assert(chem_sparse_factor(p,a));
    chem_sparse_solve(p,a,b);
    for (size_t i=0; i<4; ++i) assert(fabs(b[i]-(double)(i+1)) < 1e-14);
    double singular[] = {0,0,0,2, 3,5,0,0, 0,2,6,0, 0,0,1,7};
    assert(!chem_sparse_factor(p,singular));
    chem_sparse_plan_destroy(p);
    puts("sparse fill and solve PASS");
    return 0;
}
