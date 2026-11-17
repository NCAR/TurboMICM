#include "chem.h"
#include "chem_sparse.h"

#include <math.h>
#include <stdio.h>
#include <stdint.h>

static int failures;
static void report(const char *name, int ok)
{
    printf("rosenbrock:RosenbrockSolver.%s %s %s\n",name,ok ? "PASS" : "FAIL",ok ? "reference quantities matched" : "reference quantity differed");
    failures+=!ok;
}
static int alpha_minus(size_t cells)
{
    const ChemSparseCoordinate coords[]={{0,0},{0,1},{0,2},{1,0},{1,1},{1,2},{2,0},{2,2},{3,1},{3,3},{4,0},{4,2},{4,4}};
    const double base[]={-12.2,-24.3,-42.3,-.43,-23.4,-83.4,-4.74,-6.91,-59.1,-83.4,-78.5,-53.6,-1};
    ChemSparseMatrix m={0}; double z=0;
    if (chem_sparse_init(&m,5,cells,coords,13,CHEM_SPARSE_ROW,cells,sizeof z,&z)!=CHEM_OK) return 0;
    int lookup[25]; unsigned char algebraic[5]={0};
    for (size_t i=0;i<25;++i) lookup[i]=-1;
    for (size_t i=0;i<13;++i) {
        size_t r=coords[i].row,c=coords[i].column;
        lookup[r*5+c]=(int)chem_sparse_entry(&m,r,c);
        for (size_t cell=0;cell<cells;++cell) {
            double value=base[i]; if (i==1) value*=cell+2; if (i==5) value/=cell+3;
            *(double *)chem_sparse_at(&m,cell,r,c)=value;
        }
    }
    ChemLinearSystem view={5,cells,13,NULL,NULL,lookup,algebraic,(double *)(void *)m.storage.data};
    int ok=chem_linear_add_alpha(&view,42.042)==CHEM_OK;
    for (size_t i=0;i<13;++i) for (size_t cell=0;cell<cells;++cell) {
        double expected=base[i]; if (i==1) expected*=cell+2; if (i==5) expected/=cell+3;
        if (coords[i].row==coords[i].column) expected+=42.042;
        ok &= *(double *)chem_sparse_at(&m,cell,coords[i].row,coords[i].column)==expected;
    }
    chem_sparse_destroy(&m); return ok;
}
static int tolerances(void)
{
    const ChemSpecies sp[]={{"foo",0,1e-7,0,0},{"bar",0,1e-8,0,0}};
    const int reactants[]={0}; const ChemProduct products[]={{1,1}};
    const ChemReaction reaction={0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=2e-11,.c=110,.d=300}};
    const ChemMechanism m={2,1,0,sp,&reaction,reactants,products,NULL,NULL,298.15,101325}; int ok=1;
    for (size_t n=1;n<=10;++n) {
        ChemState *s=chem_state_create(&m,n);
        ok &= s && chem_absolute_tolerances(s)[0]==1e-7 && chem_absolute_tolerances(s)[1]==1e-8;
        chem_state_free(s);
    }
    return ok;
}
static int norm(size_t cells, int constrained)
{
    size_t species=constrained ? 3 : 5; double old[25], next[25], errors[25], atol[5]={1e-3,1e-3,1e-3,1e-3,1e-3};
    double expected=0;
    for (size_t i=0;i<cells;++i) for (size_t j=0;j<species;++j) {
        size_t k=j*cells+i;
        old[k]=constrained ? 1+i+.1*j : 7.7*i+j;
        next[k]=constrained ? .8+.5*i+.2*j : -13.9/(j+1)-i;
        errors[k]=constrained ? .01*(1+i+j) : 81.57/(i+7)/(j+3);
        double scale=atol[j]+1e-6*fmax(fabs(old[k]),fabs(next[k]));
        expected+=errors[k]*errors[k]/(scale*scale);
    }
    expected=fmax(sqrt(expected/(cells*species)),1e-10);
    double actual;
    int ok=chem_rosenbrock_normalized_error(species,cells,old,next,errors,atol,1e-6,&actual)==CHEM_OK;
    return ok && (constrained ? fabs(actual-expected)<=1e-12 : fabs(actual-expected)<=1e-11*fmax(fabs(actual),fabs(expected)));
}
static int rejected(void)
{
    const ChemSpecies sp[]={{"foo",1e12,1e-3,0,0},{"bar",1e12,1e-3,0,0},{"baz",1e12,1e-3,0,0},
                            {"quz",1e12,1e-3,0,0},{"quuz",1e12,1e-3,0,0}};
    const int reactants[]={0,2,1,3};
    const ChemProduct products[]={{1,1},{4,2.4},{0,1},{3,1.4}};
    const ChemReaction reactions[]={
        {0,2,0,2,0,{.kind=CHEM_ARRHENIUS,.a=2e-11,.c=110,.d=300}},
        {2,1,2,2,0,{.kind=CHEM_ARRHENIUS,.a=1e-6,.d=300}},
        {3,1,4,0,0,{.kind=CHEM_ARRHENIUS,.a=3.5e-6,.d=300}}};
    const ChemMechanism m={5,3,0,sp,reactions,reactants,products,NULL,NULL,298.15,101325};
    ChemState *a=chem_state_create(&m,1),*b=chem_state_create(&m,1);
    if (!a || !b) { chem_state_free(a); chem_state_free(b); return 0; }
    ChemLinearSolver alternate=chem_doolittle_linear_solver(); ChemRosenbrock p=chem_rosenbrock_three_stage();
    p.h_start=1000; ChemStats sa={0},sb={0};
    int ok=chem_set_linear_solver(b,&alternate,NULL)==CHEM_OK && chem_update_rates(a)==CHEM_OK && chem_update_rates(b)==CHEM_OK &&
           chem_solve(a,&p,1000,&sa)==CHEM_OK && chem_solve(b,&p,1000,&sb)==CHEM_OK;
    ok &= sa.internal_steps-sa.accepted>=2 && sa.internal_steps==sb.internal_steps && sa.accepted==sb.accepted && sa.rejected==sb.rejected;
    chem_state_free(a); chem_state_free(b); return ok;
}
int main(void)
{
    int standard=1,vector=1;
    for (size_t i=1;i<=4;++i) standard &= alpha_minus(i);
    const size_t vector_cells[]={1,4,3,2};
    for (size_t i=0;i<4;++i) vector &= alpha_minus(vector_cells[i]);
    report("StandardAlphaMinusJacobian",standard); report("VectorAlphaMinusJacobian",vector); report("CanSetTolerances",tolerances());
    standard=1; for (size_t i=1;i<=4;++i) standard &= norm(i,0);
    vector=standard && norm(5,0);
    report("StandardNormalizedError",standard); report("VectorNormalizedError",vector);
    report("StandardNormalizedErrorWithConstraints",norm(2,1) && norm(5,1));
    report("VectorNormalizedErrorWithConstraints",norm(3,1) && norm(5,1));
    report("RejectedStepAlphaMatchesInPlaceSolver",rejected());
    return failures ? 1 : 0;
}
