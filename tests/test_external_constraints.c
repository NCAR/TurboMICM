#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <string.h>
#include "rosenbrock.h"
#include "jacobian_check.h"

typedef struct { double k; size_t n; } Process;
static void rhs(void *v, const double *y, double *f)
{
    Process *p=v; memset(f,0,p->n*sizeof(*f)); f[0]-=p->k*y[0]; f[1]+=p->k*y[0];
}
static void jac(void *v,const double *y,double *j)
{
    Process *p=v; (void)y; memset(j,0,p->n*p->n*sizeof(*j)); j[0]-=p->k; j[p->n]+=p->k;
}
typedef struct { int equilibrium; double value; } ConstraintContext;
static double residual(void *v, const double *y, double temp, const double *parameter)
{
    ConstraintContext *c=v; (void)temp;
    double x=parameter ? *parameter : c->value;
    return c->equilibrium ? x*y[0]-y[1] : y[0]+y[1]-x;
}
static void subtract(void *v, const double *y, double temp, const double *parameter,double *row)
{
    ConstraintContext *c=v; (void)y; (void)temp;
    double x=parameter ? *parameter : c->value;
    row[0]-=c->equilibrium ? x : 1; row[1]-=c->equilibrium ? -1 : 1;
}
static double nonlinear_residual(void *ctx,const double *y,double temp,const double *parameter)
{ (void)ctx; (void)temp; (void)parameter; return y[1]*y[1]-y[0]; }
static void nonlinear_jacobian(void *ctx,const double *y,double temp,const double *parameter,double *row)
{ (void)ctx; (void)temp; (void)parameter; row[0]+=1; row[1]-=2*y[1]; }
typedef struct { const ChemConstraintSet *set; const double *temperature,*parameters; } FDContext;
static void constraint_forcing(void *v,size_t cell,const double *y,double *f)
{
    FDContext *c=v;
    chem_constraint_residuals(c->set,1,y,c->temperature+cell,c->parameters ? c->parameters+cell : NULL,f);
}
static void process_forcing(void *v,size_t cell,const double *y,double *f)
{
    (void)cell; Process *p=v; f[0]-=p->k*y[0]; f[1]+=p->k*y[0];
}
int main(void)
{
    const size_t dependencies[]={0,1}; ConstraintContext context={0,1};
    ChemCustomConstraint external={&context,dependencies,2,residual,subtract};
    ChemConstraint constraint={.name="mass",.type=CHEM_CONSTRAINT_CUSTOM,.row=1,.custom=&external};
    ChemConstraintSet set={2,1,&constraint}; assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_OK);
    unsigned char pattern[4]={0}; chem_constraint_sparsity(&set,pattern);
    assert(!pattern[0] && !pattern[1] && pattern[2] && pattern[3]);
    double states[]={.6,.4,.2,.8}, temperatures[]={298,300}, fd[8], analytical[8]={0};
    FDContext fdc={&set,temperatures,NULL};
    chem_constraint_subtract_jacobian(&set,2,states,temperatures,NULL,analytical);
    assert(!chem_finite_difference_jacobian(2,2,states,constraint_forcing,&fdc,1e-8,fd));
    assert(chem_compare_jacobian(2,2,pattern,analytical,fd,1e-7,1e-7).passed);
    assert(chem_check_jacobian_sparsity(2,2,pattern,fd,1e-6).passed);
    context=(ConstraintContext){1,2.5};
    double state_eq[]={3,5}, analytic_eq[4]={0}, fd_eq[4];
    chem_constraint_subtract_jacobian(&set,1,state_eq,temperatures,NULL,analytic_eq);
    assert(!chem_finite_difference_jacobian(1,2,state_eq,constraint_forcing,&fdc,1e-8,fd_eq));
    assert(chem_compare_jacobian(1,2,pattern,analytic_eq,fd_eq,1e-7,1e-7).passed);
    assert(chem_check_jacobian_sparsity(1,2,pattern,fd_eq,1e-6).passed);
    double parameter[]={3,4}, outputs[4]={0};
    chem_constraint_residuals(&set,2,states,temperatures,parameter,outputs);
    assert(fabs(outputs[1]-(3*.6-.4))<1e-15 && fabs(outputs[3]-(4*.2-.8))<1e-15);
    Process process={3.5,2}; const double process_states[]={.8,.2,.3,.7};
    double process_fd[8], process_analytic[8]; unsigned char process_pattern[]={1,0,1,0};
    for (size_t cell=0;cell<2;++cell) {
        jac(&process,process_states+2*cell,process_analytic+4*cell);
        for (size_t i=0;i<4;++i) process_analytic[4*cell+i]*=-1;
    }
    assert(!chem_finite_difference_jacobian(2,2,process_states,process_forcing,&process,1e-8,process_fd));
    assert(chem_compare_jacobian(2,2,process_pattern,process_analytic,process_fd,1e-7,1e-7).passed);
    assert(chem_check_jacobian_sparsity(2,2,process_pattern,process_fd,1e-6).passed);
    context=(ConstraintContext){0,1}; process.k=10;
    double state[]={.7,.2}, f[2], matrix[4];
    rhs(&process,state,f); chem_constraint_residuals(&set,1,state,temperatures,NULL,f);
    assert(fabs(f[0]+7)<1e-12 && fabs(f[1]+.1)<1e-12);
    jac(&process,state,matrix); for (size_t i=0;i<4;++i) matrix[i]*=-1;
    memset(matrix+2,0,2*sizeof(double));
    chem_constraint_subtract_jacobian(&set,1,state,temperatures,NULL,matrix);
    assert(fabs(matrix[0]-10)<1e-12 && fabs(matrix[2]+1)<1e-12 && fabs(matrix[3]+1)<1e-12);
    const double atol[]={1e-3,1e-3}; ChemRosenbrockMethod method=chem_rosenbrock_four_stage_da();
    pattern[0]=1; ChemSparsePlan *plan=chem_sparse_plan_from_pattern(2,pattern); assert(plan);
    process.k=.1;
    for (size_t cell=0;cell<3;++cell) for (size_t sparse=0;sparse<2;++sparse) {
        ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
        options.constraints=&set; options.sparse_plan=sparse ? plan : NULL;
        double y[]={.9-.2*cell,.1+.2*cell};
        assert(chem_rosenbrock_solve(2,y,50,rhs,jac,&process,&method,&options,NULL)==CHEM_SOLVE_OK);
        assert(fabs(y[0]+y[1]-1)<1e-4);
        y[0]=.9-.2*cell; y[1]=.1+.2*cell;
        for (size_t step=0;step<20;++step) {
            assert(chem_rosenbrock_solve(2,y,10,rhs,jac,&process,&method,&options,NULL)==CHEM_SOLVE_OK);
            assert(fabs(y[0]+y[1]-1)<1e-4);
        }
    }
    chem_sparse_plan_destroy(plan);
    /* A separate built-in equilibrium and external mass row can coexist. */
    const ChemConstraintTerm b[]={{2,1}}, c[]={{3,1}};
    ChemConstraint combined[]={constraint,
        {.name="eq",.type=CHEM_CONSTRAINT_EQUILIBRIUM,.row=3,.reactants=b,.products=c,
         .reactant_count=1,.product_count=1,.constant=5,.gas_constant=8.31446261815324,.reference_temperature=298.15}};
    ChemConstraintSet combined_set={4,2,combined};
    assert(chem_constraint_validate(&combined_set)==CHEM_CONSTRAINT_OK);
    unsigned char combined_pattern[16]={0}; chem_constraint_sparsity(&combined_set,combined_pattern);
    assert(combined_pattern[4] && combined_pattern[5] && combined_pattern[14] && combined_pattern[15]);
    double combined_state[]={.8,.2,.2,1}, combined_output[4]={0}, combined_jac[16]={0};
    chem_constraint_residuals(&combined_set,1,combined_state,temperatures,NULL,combined_output);
    chem_constraint_subtract_jacobian(&combined_set,1,combined_state,temperatures,NULL,combined_jac);
    assert(fabs(combined_output[1])<1e-15 && fabs(combined_output[3])<1e-15);
    assert(combined_jac[4]==-1 && combined_jac[5]==-1 && combined_jac[14]==-5 && combined_jac[15]==1);
    external.residual=nonlinear_residual; external.subtract_jacobian=nonlinear_jacobian;
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol); options.constraints=&set;
    double y[]={.64,.2};
    assert(chem_rosenbrock_initialize_constraints(2,y,&options,NULL)==CHEM_SOLVE_OK);
    assert(fabs(y[1]-.8)<1e-3 && y[0]==.64);
    assert(chem_rosenbrock_solve(2,y,1,rhs,jac,&process,&method,&options,NULL)==CHEM_SOLVE_OK);
    assert(fabs(y[1]*y[1]-y[0])<1e-3);
    size_t bad_dependency[]={2}; external.dependencies=bad_dependency; external.dependency_count=1;
    assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_UNKNOWN_SPECIES);
    external.dependencies=NULL; assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_INVALID);
    return 0;
}
