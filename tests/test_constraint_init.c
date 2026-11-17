#include "rosenbrock.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static void rhs(void *context,const double *y,double *f)
{ (void)context; f[0]=-.5*y[0]; f[1]=.5*y[0]; f[2]=0; }
static void rhs2(void *context,const double *y,double *f)
{ (void)context; f[0]=-.5*y[0]; f[1]=.5*y[0]; }
static void jac2(void *context,const double *y,double *j)
{ (void)context; (void)y; j[0]=-.5; j[1]=0; j[2]=.5; j[3]=0; }
static void jac(void *context,const double *y,double *j)
{ (void)context; (void)y; for (size_t i=0;i<9;++i) j[i]=0; j[0]=-.5; j[3]=.5; }
int main(void)
{
    const ChemConstraintTerm b[]={{1,1}}, c[]={{2,1}};
    ChemConstraint constraint={"eq",CHEM_CONSTRAINT_EQUILIBRIUM,2,b,c,1,1,2,0,8.31446261815324,298.15,NULL};
    ChemConstraintSet set={3,1,&constraint};
    const double atol[]={1e-3,1e-3,1e-3};
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol); options.constraints=&set;
    ChemRosenbrockMethod method=chem_rosenbrock_three_stage();
    ChemSolverStats stats;
    double y[]={1,.5,1};
    assert(chem_rosenbrock_solve(3,y,.001,rhs,jac,NULL,&method,&options,&stats)==CHEM_SOLVE_OK);
    assert(stats.constraint_init_iterations<=1 && fabs(2*y[1]-y[2])<=1e-6);
    const double initial_c[]={5,1000};
    for (size_t i=0;i<2;++i) {
        y[0]=1; y[1]=.5; y[2]=initial_c[i];
        assert(chem_rosenbrock_solve(3,y,.001,rhs,jac,NULL,&method,&options,&stats)==CHEM_SOLVE_OK);
        assert(stats.constraint_init_iterations>0 && fabs(2*y[1]-y[2])<=1e-6);
    }
    for (size_t cell=0;cell<4;++cell) {
        y[0]=1; y[1]=.1*(cell+1); y[2]=20;
        assert(chem_rosenbrock_solve(3,y,.001,rhs,jac,NULL,&method,&options,&stats)==CHEM_SOLVE_OK);
        assert(fabs(2*y[1]-y[2])<=1e-5);
        y[2]=999;
        assert(chem_rosenbrock_solve(3,y,.001,rhs,jac,NULL,&method,&options,&stats)==CHEM_SOLVE_OK);
        assert(fabs(2*y[1]-y[2])<=1e-6);
    }
    const double multi[][3]={{1,.5,1},{2,.3,5},{.5,.1,100}};
    for (size_t i=0;i<3;++i) {
        for(size_t j=0;j<3;++j) y[j]=multi[i][j];
        assert(chem_rosenbrock_solve(3,y,.001,rhs,jac,NULL,&method,&options,&stats)==CHEM_SOLVE_OK);
        assert(fabs(2*y[1]-y[2])<1e-5);
    }
    ChemConstraintTerm mass[]={{0,1},{1,1},{2,1}};
    ChemConstraint saved=constraint;
    method=chem_rosenbrock_four_stage_da();
    const double row_scales[]={1,1e-12,1e12};
    for(size_t i=0;i<3;++i) {
        for(size_t j=0;j<3;++j) mass[j].coefficient=row_scales[i];
        constraint=(ChemConstraint){"mass",CHEM_CONSTRAINT_LINEAR,2,mass,NULL,3,0,1.3*row_scales[i],0,0,0,NULL};
        y[0]=.7; y[1]=.15; y[2]=5;
        assert(chem_rosenbrock_solve(3,y,.001,rhs,jac,NULL,&method,&options,&stats)==CHEM_SOLVE_OK);
        assert(fabs(y[2]-.45) <= 4.5e-6);
    }
    constraint=saved; method=chem_rosenbrock_three_stage();
    const double atol2[]={1e-3,1e-3};
    ChemRosenbrockOptions ode=chem_rosenbrock_defaults(atol2);
    double ode_y[]={1,0};
    assert(chem_rosenbrock_solve(2,ode_y,.01,rhs2,jac2,NULL,&method,&ode,&stats)==CHEM_SOLVE_OK);
    assert(stats.constraint_init_iterations==0);
    ChemConstraintTerm square[]={{2,2}}; constraint.products=square;
    y[0]=1; y[1]=.5; y[2]=100; options.constraint_init_max_iterations=1;
    stats=(ChemSolverStats){0};
    assert(chem_rosenbrock_initialize_constraints(3,y,&options,&stats)==CHEM_SOLVE_CONSTRAINT_INIT);
    assert(y[0]==1 && y[1]==.5 && y[2]==100);
    options.constraint_init_max_iterations=10; y[2]=1e-6; stats=(ChemSolverStats){0};
    assert(chem_rosenbrock_initialize_constraints(3,y,&options,&stats)==CHEM_SOLVE_OK);
    assert(fabs(y[2]-1)<1e-3);
    for (size_t backtracks=0;backtracks<=2;backtracks+=2) {
        options.constraint_init_max_backtracks=backtracks; y[2]=1e-6; stats=(ChemSolverStats){0};
        assert(chem_rosenbrock_initialize_constraints(3,y,&options,&stats)==CHEM_SOLVE_CONSTRAINT_INIT);
        assert(y[1]==.5 && y[2]==1e-6);
    }
    options.constraint_init_max_backtracks=24;
    ChemConstraintTerm linear[]={{0,1},{1,-1}};
    constraint=(ChemConstraint){"scaled",CHEM_CONSTRAINT_LINEAR,1,linear,NULL,2,0,0,0,0,0,NULL};
    set=(ChemConstraintSet){2,1,&constraint};
    const double tight_atol[]={1e-12,1e-12}; options.absolute_tolerance=tight_atol; options.relative_tolerance=1e-8;
    const double scales[]={1e-12,1,1e12};
    for (size_t i=0;i<3;++i) {
        linear[0].coefficient=scales[i]; linear[1].coefficient=-scales[i];
        double x[]={1,4}; stats=(ChemSolverStats){0};
        assert(chem_rosenbrock_initialize_constraints(2,x,&options,&stats)==CHEM_SOLVE_OK);
        assert(fabs(x[1]-1)<1e-12);
    }
    linear[0].coefficient=1; linear[1].coefficient=-1;
    const double loose_atol[]={1e-8,1e-8}; options.absolute_tolerance=loose_atol; options.relative_tolerance=1e-6;
    double x[]={1,1+5e-8}, original=x[1]; stats=(ChemSolverStats){0};
    assert(chem_rosenbrock_initialize_constraints(2,x,&options,&stats)==CHEM_SOLVE_OK);
    assert(x[1]==original && stats.constraint_init_iterations==1);
    options.absolute_tolerance=tight_atol; options.relative_tolerance=1e-8; stats=(ChemSolverStats){0};
    assert(chem_rosenbrock_initialize_constraints(2,x,&options,&stats)==CHEM_SOLVE_OK);
    assert(fabs(x[1]-1)<1e-12 && stats.constraint_init_iterations==2);
    assert(!strcmp(chem_solve_status_string(CHEM_SOLVE_CONSTRAINT_INIT),"Constraint Initialization Failed"));
    puts("constraint initialization PASS");
    return 0;
}
