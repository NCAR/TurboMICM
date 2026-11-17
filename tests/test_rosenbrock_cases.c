#include "chem.h"
#include "matrix.h"
#include "sparse.h"
#include "state.h"
#include <math.h>
#include <stdio.h>

static int alpha_case(size_t cells,size_t width)
{
    const size_t points[][2]={{0,0},{0,1},{0,2},{1,0},{1,1},{1,2},{2,0},{2,2},{3,1},{3,3},{4,0},{4,2},{4,4}};
    ChemSparseBuilder b={0}; ChemSparseMatrix m={0};
    if(chem_sparse_builder_init(&b,5)) return 0;
    for(size_t i=0;i<13;++i) if(chem_sparse_builder_add(&b,points[i][0],points[i][1])) return 0;
    if(chem_sparse_init(&m,&b,cells,width,0,0)) return 0;
    chem_sparse_builder_free(&b);
    for(size_t i=0;i<cells;++i) {
        double values[]={-12.2,-24.3*(i+2),-42.3,-.43,-23.4,-83.4/(i+3),-4.74,-6.91,-59.1,-83.4,-78.5,-53.6,-1};
        for(size_t j=0;j<13;++j) *chem_sparse_at(&m,i,points[j][0],points[j][1])=values[j];
    }
    int ok=!chem_rosenbrock_alpha_minus_jacobian(&m,42.042,NULL);
    for(size_t i=0;i<cells;++i) {
        double expected[]={42.042-12.2,-24.3*(i+2),-42.3,-.43,42.042-23.4,-83.4/(i+3),-4.74,42.042-6.91,-59.1,42.042-83.4,-78.5,-53.6,42.042-1};
        for(size_t j=0;j<13;++j) if(*chem_sparse_at(&m,i,points[j][0],points[j][1])!=expected[j]) ok=0;
    }
    chem_sparse_free(&m);return ok;
}
static int norm_case(size_t cells,size_t width,int constraints)
{
    size_t n=constraints?3:5; ChemMatrix old={0},next={0},err={0};
    if(chem_matrix_init(&old,cells,n,width,width!=1,0) || chem_matrix_init(&next,cells,n,width,width!=1,0) || chem_matrix_init(&err,cells,n,width,width!=1,0)) return 0;
    double atol[]={1e-3,1e-3,1e-3,1e-3,1e-3}, y0[50],y1[50],e[50],sum=0;
    for(size_t i=0;i<cells;++i) for(size_t j=0;j<n;++j) {
        double a=constraints?1.0+i+.1*j:7.7*i+j;
        double b=constraints?.8+.5*i+.2*j:-13.9/(j+1)-i;
        double q=constraints?.01*(1+i+j):81.57/(i+7)/(j+3);
        *chem_matrix_at(&old,i,j)=a;*chem_matrix_at(&next,i,j)=b;*chem_matrix_at(&err,i,j)=q;
        double scale=atol[j]+1e-6*fmax(fabs(a),fabs(b));sum+=q*q/(scale*scale);
    }
    for(size_t i=0;i<cells;++i) {chem_matrix_get_row(&old,i,y0+i*n,n);chem_matrix_get_row(&next,i,y1+i*n,n);chem_matrix_get_row(&err,i,e+i*n,n);}
    double expected=fmax(sqrt(sum/(cells*n)),1e-10);
    double computed=chem_rosenbrock_normalized_error(cells,n,y0,y1,e,1e-6,0,atol);
    int ok=constraints?fabs(computed-expected)<=1e-12:fabs(computed-expected)/fmax(fabs(computed),fabs(expected))<=1e-11;
    chem_matrix_free(&old);chem_matrix_free(&next);chem_matrix_free(&err);return ok;
}
static int tolerances(void)
{
    ChemSpecies foo,bar; chem_species_init(&foo,"foo");chem_species_init(&bar,"bar");
    ChemPropertyValue a={.real=1e-7},b={.real=1e-8};
    if(chem_species_set(&foo,"absolute tolerance",CHEM_PROPERTY_REAL,a) || chem_species_set(&bar,"absolute tolerance",CHEM_PROPERTY_REAL,b)) return 0;
    ChemPhaseSpecies members[]={{&foo,0,0},{&bar,0,0}}; ChemPhase phase={"gas",members,2}; const char *names[]={"foo","bar"};
    int ok=1;
    for(size_t cells=1;cells<=10;++cells) {
        ChemState state={0}; ChemStateParameters p={1,2,0,names,NULL,1e-6,NULL};
        if(chem_state_init(&state,&p,cells,1) || chem_state_apply_phase_tolerances(&state,phase)) ok=0;
        if(state.parameters.variable_count!=2 || state.absolute_tolerances[0]!=1e-7 || state.absolute_tolerances[1]!=1e-8) ok=0;
        chem_state_free(&state);
    }
    chem_species_free(&foo);chem_species_free(&bar);return ok;
}
int main(void)
{
    int standard=1,vector=1,standard_norm=1,vector_norm=1,standard_constraint=1,vector_constraint=1;
    for(size_t cells=1;cells<=4;++cells) {standard &= alpha_case(cells,1);standard_norm &= norm_case(cells,1,0);}
    size_t acells[]={1,4,3,2};for(size_t width=1;width<=4;++width) vector &= alpha_case(acells[width-1],width);
    const size_t widths[]={1,2,3,4,2,3,4,8,10},counts[]={1,2,3,4,1,2,3,5,3};
    for(size_t i=0;i<9;++i) vector_norm &= norm_case(counts[i],widths[i],0);
    for(int method=0;method<5;++method) {standard_constraint &= norm_case(2,1,1);vector_constraint &= norm_case(3,4,1);}
    standard_constraint &= norm_case(5,1,1);vector_constraint &= norm_case(5,4,1);
    printf("rosenbrock:RosenbrockSolver.StandardAlphaMinusJacobian %d\n",standard);
    printf("rosenbrock:RosenbrockSolver.VectorAlphaMinusJacobian %d\n",vector);
    printf("rosenbrock:RosenbrockSolver.CanSetTolerances %d\n",tolerances());
    printf("rosenbrock:RosenbrockSolver.StandardNormalizedError %d\n",standard_norm);
    printf("rosenbrock:RosenbrockSolver.VectorNormalizedError %d\n",vector_norm);
    printf("rosenbrock:RosenbrockSolver.StandardNormalizedErrorWithConstraints %d\n",standard_constraint);
    printf("rosenbrock:RosenbrockSolver.VectorNormalizedErrorWithConstraints %d\n",vector_constraint);
    return 0;
}
