#include "chem.h"
#include "rosenbrock.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
extern const ChemMechanism generated_mechanism;
extern const char *const generated_species[];
typedef struct { const double *rates; } Context;
static void rhs(void *c,const double *y,double *f) { chem_rhs(&generated_mechanism,y,((Context*)c)->rates,f); }
static void jac(void *c,const double *y,double *j) { chem_jacobian(&generated_mechanism,y,((Context*)c)->rates,j); }
int main(void)
{
    size_t index[4] = {4,4,4,4};
    for (size_t i=0; i<4; ++i)
        for (size_t j=0; j<4; ++j)
            if (generated_species[i][0] == 'A'+(int)j && generated_species[i][1] == '\0') index[j]=i;
    for (size_t i=0; i<4; ++i) assert(index[i]<4);
    assert(index[0] != 0); /* The fixture must actually be permuted. */
    double y[4], dense_y[4], atol[] = {1e-10,1e-10,1e-10,1e-10}, f[4];
    for (size_t i=0; i<4; ++i) y[index[i]] = (double)(i+1);
    double rates[] = {2}; Context context = {rates};
    rhs(&context,y,f);
    assert(f[index[0]] == -4 && f[index[1]] == -4 && f[index[2]] == 4 && f[index[3]] == 4);
    memcpy(dense_y,y,sizeof(y));
    ChemRosenbrockMethod method = chem_rosenbrock_three_stage();
    ChemRosenbrockOptions options = chem_rosenbrock_defaults(atol);
    assert(chem_rosenbrock_solve(4,dense_y,.1,rhs,jac,&context,&method,&options,NULL) == CHEM_SOLVE_OK);
    ChemSparsePlan *plan = chem_sparse_plan_create(&generated_mechanism);
    assert(plan); options.sparse_plan = plan;
    assert(chem_rosenbrock_solve(4,y,.1,rhs,jac,&context,&method,&options,NULL) == CHEM_SOLVE_OK);
    for (size_t i=0; i<4; ++i) assert(fabs(y[i]-dense_y[i]) < 1e-13);
    chem_sparse_plan_destroy(plan);
    puts("reordered generated mechanism PASS");
    return 0;
}
