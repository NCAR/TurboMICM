#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <string.h>
#include "rosenbrock.h"

typedef struct { const ChemMechanism *m; const double *k; } System;
static void rhs(void *v,const double *y,double *f) { System *s=v; chem_rhs(s->m,y,s->k,f); }
static void jac(void *v,const double *y,double *j) { System *s=v; chem_jacobian(s->m,y,s->k,j); }
typedef struct { size_t r,p; double k; int temperature; } Eq;
static double equilibrium(void *v,const double *y,double t,const double *parameter)
{
    Eq *c=v; double k=parameter ? *parameter : c->k;
    if (c->temperature && !parameter) k*=exp(3000*(1/298.15-1/t));
    return k*y[c->r]-y[c->p];
}
static void equilibrium_jac(void *v,const double *y,double t,const double *parameter,double *row)
{
    Eq *c=v; (void)y; double k=parameter ? *parameter : c->k;
    if (c->temperature && !parameter) k*=exp(3000*(1/298.15-1/t));
    row[c->r]-=k; row[c->p]+=1;
}
typedef struct { size_t n; double total; } Mass;
static double mass(void *v,const double *y,double t,const double *parameter)
{
    Mass *m=v; (void)t; double total=parameter ? *parameter : m->total;
    for (size_t i=0;i<m->n;++i) total-=y[i];
    return -total;
}
static void mass_jac(void *v,const double *y,double t,const double *parameter,double *row)
{
    Mass *m=v; (void)y; (void)t; (void)parameter;
    for (size_t i=0;i<m->n;++i) row[i]-=1;
}
static void checked_solve(System *sys,double *y,double dt,ChemRosenbrockMethod *method,ChemRosenbrockOptions *options)
{
    assert(chem_rosenbrock_solve(sys->m->species_count,y,dt,rhs,jac,sys,method,options,NULL)==CHEM_SOLVE_OK);
}
int main(void)
{
    const size_t driver_react[]={0}, kin_react[]={0,1,2,1,3};
    const ChemProduct driver_prod[]={{1,1}}, kin_prod[]={{1,1},{2,1},{1,1},{3,1},{1,1}};
    const ChemReaction driver_rx[]={{0,1,0,1}}, kin_rx[]={{0,1,0,1},{1,1,1,1},{2,1,2,1},{3,1,3,1},{4,1,4,1}};
    const ChemMechanism drive3={3,1,driver_rx,driver_react,driver_prod};
    const ChemMechanism drive4={4,1,driver_rx,driver_react,driver_prod};
    const ChemMechanism kin3={3,3,kin_rx,kin_react,kin_prod};
    const ChemMechanism kin4={4,5,kin_rx,kin_react,kin_prod};
    const double driver_k[]={.1}, kin_k[]={.1,10,2,8,8.0/3};
    System driver={&drive3,driver_k}, kinetic={&kin3,kin_k};
    const double atol[]={1e-3,1e-3,1e-3,1e-3};
    Eq e1={1,2,5,0},e2={1,3,3,0}; Mass m={3,1};
    const size_t dep_e1[]={1,2}, dep_e2[]={1,3}, dep_mass3[]={0,1,2}, dep_mass4[]={0,1,2,3};
    ChemCustomConstraint custom_e1={&e1,dep_e1,2,equilibrium,equilibrium_jac};
    ChemCustomConstraint custom_e2={&e2,dep_e2,2,equilibrium,equilibrium_jac};
    ChemCustomConstraint custom_mass={&m,dep_mass3,3,mass,mass_jac};
    ChemConstraint rows[]={
        {.name="mass",.type=CHEM_CONSTRAINT_CUSTOM,.row=1,.custom=&custom_mass},
        {.name="eq1",.type=CHEM_CONSTRAINT_CUSTOM,.row=2,.custom=&custom_e1},
        {.name="eq2",.type=CHEM_CONSTRAINT_CUSTOM,.row=3,.custom=&custom_e2}};
    ChemConstraintSet set3={3,2,rows}, simple={3,1,rows+1}, set4={4,3,rows};
    ChemRosenbrockMethod ode_method=chem_rosenbrock_three_stage(), dae_method=chem_rosenbrock_four_stage_da();
    ChemRosenbrockOptions ode=chem_rosenbrock_defaults(atol), dae=chem_rosenbrock_defaults(atol);
    double ky[]={1,0,0},cy[]={1,0,0},sy[]={1,0,0};
    for (size_t step=0;step<200;++step) {
        checked_solve(&kinetic,ky,1,&ode_method,&ode);
        dae.constraints=&set3; checked_solve(&driver,cy,1,&dae_method,&dae);
        dae.constraints=&simple; checked_solve(&driver,sy,1,&dae_method,&dae);
    }
    for (size_t i=0;i<3;++i) {
        const double expected[]={0,1.0/6,5.0/6}, tol[]={1e-6,1e-3,1e-3};
        assert(fabs(ky[i]-expected[i])<tol[i] && fabs(cy[i]-expected[i])<tol[i]);
        assert(fabs(ky[i]-cy[i])<tol[i]);
    }
    assert(fabs(ky[0]+ky[1]+ky[2]-1)<1e-6 && fabs(cy[0]+cy[1]+cy[2]-1)<1e-3);
    assert(fabs(ky[2]/ky[1]-5)<1e-3 && fabs(sy[2]/sy[1]-5)<1e-3);
    assert(fabs(sy[0]+sy[1]+sy[2]-6)<.1);
    /* Compare a built-in and callback equilibrium at every reference step. */
    const ChemConstraintTerm b[]={{1,1}},c[]={{2,1}};
    ChemConstraint builtin={.name="built-in",.type=CHEM_CONSTRAINT_EQUILIBRIUM,.row=2,
        .reactants=b,.products=c,.reactant_count=1,.product_count=1,.constant=5,.gas_constant=8.31446261815324,.reference_temperature=298.15};
    ChemConstraintSet builtin_set={3,1,&builtin};
    double ext[]={1,0,0},bi[]={1,0,0};
    for (size_t step=0;step<100;++step) {
        dae.constraints=&simple; checked_solve(&driver,ext,.5,&dae_method,&dae);
        dae.constraints=&builtin_set; checked_solve(&driver,bi,.5,&dae_method,&dae);
        for (size_t i=0;i<3;++i) assert(fabs(ext[i]-bi[i])<1e-8);
        assert(fabs(5*ext[1]-ext[2])<1e-6 && fabs(5*bi[1]-bi[2])<1e-6);
    }
    /* Three composed callbacks, including a conservation row. */
    driver.m=&drive4; kinetic.m=&kin4; m.n=4; custom_mass.dependencies=dep_mass4; custom_mass.dependency_count=4;
    double k4[]={1,0,0,0},e4[]={1,0,0,0}; dae.constraints=&set4;
    assert(chem_constraint_validate(&set4)==CHEM_CONSTRAINT_OK);
    for (size_t step=0;step<200;++step) {
        checked_solve(&kinetic,k4,1,&ode_method,&ode); checked_solve(&driver,e4,1,&dae_method,&dae);
    }
    const double steady[]={0,1.0/9,5.0/9,3.0/9};
    for (size_t i=0;i<4;++i) {
        double tolerance=i ? 1e-3 : 1e-6;
        assert(fabs(k4[i]-steady[i])<tolerance && fabs(e4[i]-steady[i])<tolerance);
        if (i) assert(fabs(k4[i]-e4[i])<1e-3);
    }
    assert(fabs(k4[0]+k4[1]+k4[2]+k4[3]-1)<1e-6);
    assert(fabs(e4[0]+e4[1]+e4[2]+e4[3]-1)<1e-3);
    assert(fabs(5*e4[1]-e4[2])<1e-6 && fabs(3*e4[1]-e4[3])<1e-6);
    driver.m=&drive3; e1=(Eq){1,2,2,1};
    for (size_t j=0;j<2;++j) {
        double temp=j ? 350 : 298.15, y[]={1,0,0}; dae.constraints=&simple; dae.constraint_temperature=temp;
        double expected=2*exp(3000*(1/298.15-1/temp)), updated_parameter=0;
        dae.constraint_parameters=&updated_parameter;
        for (size_t step=0;step<200;++step) {
            updated_parameter=2*exp(3000*(1/298.15-1/temp));
            checked_solve(&driver,y,1,&dae_method,&dae);
        }
        assert(y[1]>0 && fabs(y[2]/y[1]-expected)<1e-4);
        if (j) assert(expected>2);
    }
    return 0;
}
