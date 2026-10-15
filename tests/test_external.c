/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Port of integration/test_external_model_constraints.cpp and its stub models. */
#include "test.h"
#include "chem_external.h"
#include "chem_constraints.h"
#include "chem_jacobian.h"
#include "../src/solver_internal.h"
#include <stdlib.h>

typedef enum { AEROSOL, EQUILIBRIUM, CONSERVATIVE, MASS } ModelKind;
typedef struct { ModelKind kind;double k,total;int solvent,temperature; } ModelData;
typedef struct {
    ChemExternalModel model;
    ModelData data;
    const char *variables[4],*algebraic[2],*parameters[1];
    ChemNamedJacobianElement process[4],constraint[5];
    ChemSpecies added;
} Model;
static void update_parameters(const ChemModelBinding *b,const void *context,ChemState *s)
{
    const ModelData *d=context;
    for (size_t c=0;c<s->cells;++c)
        s->custom[b->parameters[0]*s->cells+c]=d->k*exp(3000*(1/298.15-1/s->conditions[c].temperature));
}
static void process(const ChemModelBinding *b,const void *context,const ChemState *s,const double *y,double *f)
{
    const ModelData *d=context;size_t n=s->cells;
    for (size_t c=0;c<n;++c) {
        double rate=d->k*y[b->variables[0]*n+c];
        if (d->solvent) rate*=y[b->variables[2]*n+c];
        f[b->variables[0]*n+c]-=rate;f[b->variables[1]*n+c]+=rate;
    }
}
static void process_jac(const ChemModelBinding *b,const void *context,const ChemState *s,const double *y,ChemModelJacobian *j)
{
    const ModelData *d=context;size_t n=s->cells,g=b->variables[0],a=b->variables[1];
    for (size_t c=0;c<n;++c) {
        double rate=d->k*(d->solvent?y[b->variables[2]*n+c]:1);
        CHECK(chem_model_jacobian_add(j,g,g,c,-rate)==CHEM_OK);
        CHECK(chem_model_jacobian_add(j,a,g,c,rate)==CHEM_OK);
        if (d->solvent) {
            CHECK(chem_model_jacobian_add(j,g,b->variables[2],c,-d->k*y[g*n+c])==CHEM_OK);
            CHECK(chem_model_jacobian_add(j,a,b->variables[2],c,d->k*y[g*n+c])==CHEM_OK);
        }
    }
}
static void residual(const ChemModelBinding *b,const void *context,const ChemState *s,const double *y,double *f)
{
    const ModelData *d=context;size_t n=s->cells;
    for (size_t c=0;c<n;++c) {
        if (d->kind==EQUILIBRIUM) {
            double k=d->temperature?s->custom[b->parameters[0]*n+c]:d->k;
            f[b->algebraic[0]*n+c]=k*y[b->variables[0]*n+c]-y[b->variables[1]*n+c];
        } else {
            double total=-d->total;
            size_t count=d->kind==AEROSOL?2:b->variable_count;
            for (size_t v=0;v<count;++v) total+=y[b->variables[v]*n+c];
            f[b->algebraic[0]*n+c]=total;
            if (d->kind==CONSERVATIVE)
                f[b->algebraic[1]*n+c]=d->k*y[b->variables[1]*n+c]-y[b->variables[2]*n+c];
        }
    }
}
static void constraint_jac(const ChemModelBinding *b,const void *context,const ChemState *s,const double *y,ChemModelJacobian *j)
{
    (void)y;const ModelData *d=context;
    for (size_t c=0;c<s->cells;++c) {
        if (d->kind==EQUILIBRIUM) {
            double k=d->temperature?s->custom[b->parameters[0]*s->cells+c]:d->k;
            CHECK(chem_model_jacobian_add(j,b->algebraic[0],b->variables[0],c,k)==CHEM_OK);
            CHECK(chem_model_jacobian_add(j,b->algebraic[0],b->variables[1],c,-1)==CHEM_OK);
        } else {
            size_t count=d->kind==AEROSOL?2:b->variable_count;
            for (size_t v=0;v<count;++v)
                CHECK(chem_model_jacobian_add(j,b->algebraic[0],b->variables[v],c,1)==CHEM_OK);
            if (d->kind==CONSERVATIVE) {
                CHECK(chem_model_jacobian_add(j,b->algebraic[1],b->variables[1],c,d->k)==CHEM_OK);
                CHECK(chem_model_jacobian_add(j,b->algebraic[1],b->variables[2],c,-1)==CHEM_OK);
            }
        }
    }
}
static void model_init(Model *m,ModelKind kind,double k,int constraint,int solvent,int temperature,
                       size_t count,const char *const *variables,const char *row)
{
    memset(m,0,sizeof(*m));m->data=(ModelData){kind,k,1,solvent,temperature};
    for (size_t i=0;i<count;++i) m->variables[i]=variables[i];
    m->model=(ChemExternalModel){.variable_count=count,.variables=m->variables,
        .context=&m->data,.context_size=sizeof(m->data)};
    if (kind==AEROSOL) {
        m->added=(ChemSpecies){.name="AEROSOL.A_AQ"};
        m->model.additional_species_count=1;m->model.additional_species=&m->added;
        m->process[0]=(ChemNamedJacobianElement){variables[0],variables[0]};
        m->process[1]=(ChemNamedJacobianElement){variables[1],variables[0]};
        if (solvent) {
            m->process[2]=(ChemNamedJacobianElement){variables[0],variables[2]};
            m->process[3]=(ChemNamedJacobianElement){variables[1],variables[2]};
        }
        m->model.process_element_count=solvent?4:2;m->model.process_elements=m->process;
        m->model.add_forcing=process;m->model.process_jacobian=process_jac;
    }
    if (constraint) {
        m->algebraic[0]=row;m->model.algebraic=m->algebraic;m->model.algebraic_count=1;
        size_t nc=kind==AEROSOL?2:count;
        for (size_t i=0;i<nc;++i) m->constraint[i]=(ChemNamedJacobianElement){row,variables[i]};
        if (kind==CONSERVATIVE) {
            m->algebraic[1]=variables[2];m->model.algebraic_count=2;
            m->constraint[nc++]=(ChemNamedJacobianElement){variables[2],variables[1]};
            m->constraint[nc++]=(ChemNamedJacobianElement){variables[2],variables[2]};
        }
        m->model.constraint_element_count=nc;m->model.constraint_elements=m->constraint;
        m->model.constraint_residual=residual;m->model.constraint_jacobian=constraint_jac;
        if (temperature) {
            m->parameters[0]="C_K_eq";m->model.parameters=m->parameters;m->model.parameter_count=1;
            m->model.update_constraint_parameters=update_parameters;
        }
    }
}
static const char *const aerosol_names[]={"A_GAS","AEROSOL.A_AQ","S"};
static const char *const abc[]={"A","B","C","D"};
typedef struct { ChemMechanism *m;ChemState *s;ChemSolver *solver;ChemRosenbrockParameters options; } Fixture;
static Fixture build(size_t count,const char *const *names,size_t nr,const ChemReaction *reactions,
                      size_t nm,const ChemExternalModel *models,size_t nc,const ChemConstraint *constraints,size_t cells,int reorder)
{
    Fixture f={0};ChemSpecies species[8]={0};CHECK(count<=8);
    for (size_t i=0;i<count;++i) species[i].name=names[i];
    CHECK(chem_mechanism_create_with_models(count,species,nr,reactions,nm,models,&f.m)==CHEM_OK);
    CHECK(chem_mechanism_set_constraints(f.m,nc,constraints)==CHEM_OK);
    if (reorder) CHECK(chem_mechanism_reorder(f.m)==CHEM_OK);
    CHECK(chem_state_create(f.m,cells,&f.s)==CHEM_OK);
    CHECK(chem_solver_create(f.m,cells,&f.solver)==CHEM_OK);
    f.options=chem_rosenbrock_parameters(chem_constraint_count(f.m)?CHEM_ROS4_DA:CHEM_ROS3);
    for (size_t c=0;c<cells;++c) f.s->conditions[c]=(ChemConditions){298,101325,0};
    return f;
}
static void release(Fixture f) { chem_solver_free(f.solver);chem_state_free(f.s);chem_mechanism_free(f.m); }
static size_t index_of(Fixture f,const char *name)
{ size_t index=SIZE_MAX;CHECK(chem_species_index(f.m,name,&index)==CHEM_OK);return index; }
static double get(Fixture f,const char *name,size_t cell) { return f.s->values[index_of(f,name)*f.s->cells+cell]; }
static void set(Fixture f,const char *name,double value) { CHECK(chem_state_fill_concentration(f.s,name,value)==CHEM_OK); }
static void step(Fixture f,double dt)
{ ChemSolverStats stats;chem_update_rates(f.s);CHECK(chem_solve(f.solver,f.s,dt,&f.options,&stats)==CHEM_OK); }
static ChemReaction reaction(size_t *from,ChemProduct *to,double rate)
{
    ChemRate r=chem_rate_default(CHEM_ARRHENIUS);r.data.arrhenius.A=rate;
    return (ChemReaction){1,1,from,to,r};
}
static ChemConstraint builtin(const char *name,const char *row,ChemConstraintTerm *r,ChemConstraintTerm *p,double k)
{
    return (ChemConstraint){.type=CHEM_EQUILIBRIUM_CONSTRAINT,.name=name,.algebraic_species=row,
        .reactant_count=1,.product_count=1,.reactants=r,.products=p,.equilibrium=chem_vant_hoff(k,0)};
}
static Fixture aerosol_fixture(double k,int constraint,size_t cells,int solvent,int reorder)
{
    Model a;model_init(&a,AEROSOL,k,constraint,solvent,0,solvent?3:2,aerosol_names,aerosol_names[1]);
    const char *names[]={"A_GAS","S"};
    return build(solvent?2:1,names,0,NULL,1,&a.model,0,NULL,cells,reorder);
}
static void metadata(int constraints)
{
    Fixture f=aerosol_fixture(.01,constraints,1,0,0);
    CHECK(chem_species_count(f.m)==2);CHECK(chem_constraint_count(f.m)==(size_t)constraints);
    CHECK(!chem_is_algebraic(f.m,index_of(f,"A_GAS")));
    CHECK(chem_is_algebraic(f.m,index_of(f,"AEROSOL.A_AQ"))==constraints);
    if (constraints) CHECK(chem_constraint_row(f.m,0)==index_of(f,"AEROSOL.A_AQ"));
    double mass[4]={0};chem_add_mass_diagonal(f.m,f.solver->sparse,1,1,mass);
    CHECK(mass[chem_sparse_index(f.solver->sparse,0,0)]==1);
    CHECK(mass[chem_sparse_index(f.solver->sparse,1,1)]==(constraints?0:1));
    release(f);
}
static void conservation(int mode)
{
    size_t cells=mode==2?3:1;double k=mode==1?.01:.1;
    Fixture f=aerosol_fixture(k,1,cells,mode==3,0);
    size_t g=index_of(f,"A_GAS"),a=index_of(f,"AEROSOL.A_AQ");
    CHECK(chem_constraint_count(f.m)==1);CHECK(chem_is_algebraic(f.m,a));
    for (size_t c=0;c<cells;++c) {
        f.s->values[g*cells+c]=mode==1?.8:.9-.2*c;
        f.s->values[a*cells+c]=1-f.s->values[g*cells+c];
    }
    if (mode==3) set(f,"S",1);
    size_t steps=mode==0?20:mode==3?10:1;double dt=mode==0 || mode==3?10:50;
    for (size_t i=0;i<steps;++i) {
        step(f,dt);
        for (size_t c=0;c<cells;++c) NEAR(f.s->values[g*cells+c]+f.s->values[a*cells+c],1,1e-4);
    }
    if (mode==3) {
        size_t col=index_of(f,"S"),entry=chem_sparse_index(f.solver->sparse,a,col);
        CHECK(entry!=SIZE_MAX);chem_solver_negative_jacobian(f.solver,f.s,f.s->values);
        CHECK(f.solver->matrix[entry]==0);
    }
    release(f);
}
static void algebraic_rows(void)
{
    Fixture f=aerosol_fixture(10,1,1,0,0);set(f,"A_GAS",.7);set(f,"AEROSOL.A_AQ",.2);
    double forcing[2],jac[4];chem_forcing(f.s,f.s->values,forcing);chem_jacobian(f.s,f.s->values,jac);
    NEAR(forcing[0],-7,1e-12);NEAR(forcing[1],.7+.2-1,1e-12);
    NEAR(-jac[0],10,1e-12);NEAR(-jac[2],-1,1e-12);NEAR(-jac[3],-1,1e-12);
    chem_solver_negative_jacobian(f.solver,f.s,f.s->values);
    NEAR(f.solver->matrix[chem_sparse_index(f.solver->sparse,0,0)],10,1e-12);
    NEAR(f.solver->matrix[chem_sparse_index(f.solver->sparse,1,0)],-1,1e-12);
    NEAR(f.solver->matrix[chem_sparse_index(f.solver->sparse,1,1)],-1,1e-12);
    release(f);
}
static void combined(int filtering)
{
    const char *names[]={"A_GAS","B","C"};Model a;
    model_init(&a,AEROSOL,filtering?.1:.01,!filtering,0,0,2,aerosol_names,aerosol_names[1]);
    ChemConstraintTerm r={"B",1},p={filtering?"AEROSOL.A_AQ":"C",1};
    ChemConstraint c=builtin("B_C_eq",p.species,&r,&p,5);
    size_t from=0;ChemProduct to={1,1};ChemReaction rxn=reaction(&from,&to,.05);
    Fixture f=build(filtering?2:3,names,filtering?0:1,&rxn,1,&a.model,1,&c,1,0);
    CHECK(chem_constraint_count(f.m)==(filtering?1:2));CHECK(chem_is_algebraic(f.m,index_of(f,"AEROSOL.A_AQ")));
    if (filtering) CHECK(chem_sparse_index(f.solver->sparse,index_of(f,"AEROSOL.A_AQ"),index_of(f,"A_GAS"))!=SIZE_MAX);
    else CHECK(chem_is_algebraic(f.m,index_of(f,"C")));
    release(f);
}
static Fixture chemical(int mode,size_t n,int reorder)
{
    /* mode 0 is kinetic, 1 is conservative DAE, 2 is ratio-only DAE. */
    size_t from[]={0,1,2,1,3};ChemProduct to[]={{1,1},{2,1},{1,1},{3,1},{1,1}};
    double rates[]={.1,10,2,8,8.0/3};ChemReaction reactions[5];
    for (size_t i=0;i<5;++i) reactions[i]=reaction(from+i,to+i,rates[i]);
    Model models[3];ChemExternalModel definitions[3];size_t nm=0;
    if (mode && n==3) {
        model_init(models,mode==1?CONSERVATIVE:EQUILIBRIUM,5,1,0,0,mode==1?3:2,mode==1?abc:abc+1,mode==1?"B":"C");
        definitions[nm++]=models[0].model;
    } else if (mode) {
        const char *bd[]={"B","D"};
        model_init(models,EQUILIBRIUM,5,1,0,0,2,abc+1,"C");
        model_init(models+1,EQUILIBRIUM,3,1,0,0,2,bd,"D");
        model_init(models+2,MASS,0,1,0,0,4,abc,"B");
        for (size_t i=0;i<3;++i) definitions[nm++]=models[i].model;
    }
    Fixture f=build(n,abc,mode?1:n==3?3:5,reactions,nm,definitions,0,NULL,1,reorder);
    set(f,"A",1);return f;
}
static void convergence(int ratio,int multi)
{
    size_t n=multi?4:3;Fixture k=chemical(0,n,0),d=chemical(ratio?2:1,n,0);
    for (size_t i=0;i<200;++i) { step(k,1);step(d,1); }
    double kin[4],dae[4],ksum=0,dsum=0;
    for (size_t i=0;i<n;++i) { kin[i]=get(k,abc[i],0);dae[i]=get(d,abc[i],0);ksum+=kin[i];dsum+=dae[i]; }
    if (ratio) { NEAR(kin[2]/kin[1],5,1e-3);NEAR(dae[2]/dae[1],5,1e-3);NEAR(dsum,6,.1); }
    else {
        double expected[]={0,multi?1.0/9:1.0/6,multi?5.0/9:5.0/6,1.0/3};
        for (size_t i=0;i<n;++i) {
            double tol=i?1e-3:1e-6;
            NEAR(kin[i],expected[i],tol);NEAR(dae[i],expected[i],tol);NEAR(kin[i],dae[i],tol);
        }
        NEAR(dsum,1,1e-3);
        if (multi) { NEAR(5*dae[1]-dae[2],0,1e-6);NEAR(3*dae[1]-dae[3],0,1e-6); }
    }
    NEAR(ksum,1,1e-6);release(k);release(d);
}
static void stepwise(void)
{
    Model model;model_init(&model,EQUILIBRIUM,5,1,0,0,2,abc+1,"C");
    size_t from=0;ChemProduct to={1,1};ChemReaction rxn=reaction(&from,&to,.1);
    ChemConstraintTerm r={"B",1},p={"C",1};ChemConstraint constraint=builtin("B_C_eq","C",&r,&p,5);
    Fixture bi=build(3,abc,1,&rxn,0,NULL,1,&constraint,1,0),ext=build(3,abc,1,&rxn,1,&model.model,0,NULL,1,0);
    set(bi,"A",1);set(ext,"A",1);
    for (size_t i=0;i<100;++i) {
        step(bi,.5);step(ext,.5);
        for (size_t v=0;v<3;++v) NEAR(get(bi,abc[v],0),get(ext,abc[v],0),1e-8);
        NEAR(5*get(bi,"B",0)-get(bi,"C",0),0,1e-6);NEAR(5*get(ext,"B",0)-get(ext,"C",0),0,1e-6);
    }
    release(bi);release(ext);
}
static void fd_process(const double *y,double *f,void *context) { chem_external_add_forcing(context,y,f); }
static void fd_constraint(const double *y,double *f,void *context) { chem_constraints_forcing(context,y,f); }
static void finite_difference(int which)
{
    Fixture f;
    if (which<2) f=aerosol_fixture(which?.1:3.5,which,2,0,0);
    else {
        Model model;model_init(&model,EQUILIBRIUM,2.5,1,0,0,2,abc,"B");
        f=build(2,abc,0,NULL,1,&model.model,0,NULL,1,0);
    }
    size_t cells=f.s->cells;
    if (which==0) { f.s->values[0]=.8;f.s->values[1]=.3;f.s->values[2]=.2;f.s->values[3]=.7; }
    if (which==1) { f.s->values[0]=.6;f.s->values[1]=.2;f.s->values[2]=.4;f.s->values[3]=.8; }
    if (which==2) { f.s->values[0]=3;f.s->values[1]=5; }
    unsigned char pattern[4];
    for (size_t r=0;r<2;++r) for (size_t c=0;c<2;++c) pattern[r*2+c]=(unsigned char)chem_external_has_element(f.m,r,c,which!=0);
    ChemSparse *sparse=NULL;CHECK(chem_sparse_create(2,pattern,&sparse)==CHEM_OK);
    double dense[8]={0},fd[8]={0},packed[8]={0};
    if (which) chem_constraints_subtract_jacobian(f.s,f.s->values,dense);
    else chem_external_subtract_jacobian(f.s,f.s->values,dense);
    for (size_t r=0;r<2;++r) for (size_t c=0;c<2;++c) {
        size_t index=chem_sparse_index(sparse,r,c);
        if (index!=SIZE_MAX) for (size_t cell=0;cell<cells;++cell) packed[index*cells+cell]=dense[(r*2+c)*cells+cell];
    }
    CHECK(chem_finite_difference(2,cells,which?fd_constraint:fd_process,f.s,f.s->values,1e-8,fd)==CHEM_OK);
    CHECK(chem_compare_jacobian(sparse,cells,packed,fd,1e-7,1e-7).passed);
    if (which!=1) CHECK(chem_check_sparsity(sparse,cells,fd,1e-6).passed);
    chem_sparse_free(sparse);release(f);
}
static void temperature(void)
{
    Model model;model_init(&model,EQUILIBRIUM,2,1,0,1,2,abc+1,"C");
    size_t from=0;ChemProduct to={1,1};ChemReaction rxn=reaction(&from,&to,.1);
    Fixture f=build(3,abc,1,&rxn,1,&model.model,0,NULL,1,0);
    for (size_t i=0;i<2;++i) {
        set(f,"A",1);set(f,"B",0);set(f,"C",0);f.s->conditions[0].temperature=i?350:298.15;
        for (size_t j=0;j<200;++j) step(f,1);
        double k=2*exp(3000*(1/298.15-1/f.s->conditions[0].temperature));
        CHECK(get(f,"B",0)>0);NEAR(get(f,"C",0)/get(f,"B",0),k,1e-4);
        if (i) CHECK(k>2);
    }
    release(f);
}
static void no_op(const ChemModelBinding *b,const void *c,const ChemState *s,const double *y,double *f)
{ (void)b;(void)c;(void)s;(void)y;(void)f; }
static void model_order(void)
{
    Model model;model_init(&model,EQUILIBRIUM,5,1,0,0,2,abc+1,"C");
    ChemExternalModel models[]={model.model,{.add_forcing=no_op}};
    size_t from=0;ChemProduct to={1,1};ChemReaction rxn=reaction(&from,&to,.1);
    double values[2][2];
    for (size_t order=0;order<2;++order) {
        Fixture f=build(3,abc,1,&rxn,2,models,0,NULL,1,0);
        CHECK(chem_constraint_count(f.m)==1);CHECK(chem_is_algebraic(f.m,index_of(f,"C")));set(f,"A",1);
        for (size_t i=0;i<20;++i) step(f,1);
        values[order][0]=get(f,"B",0);values[order][1]=get(f,"C",0);
        CHECK(values[order][0]>0);NEAR(5*values[order][0]-values[order][1],0,1e-6);
        ChemExternalModel swap=models[0];models[0]=models[1];models[1]=swap;release(f);
    }
    CHECK(values[0][0]==values[1][0]);CHECK(values[0][1]==values[1][1]);
}
static void multiple_builds(void)
{
    Model model;model_init(&model,EQUILIBRIUM,5,1,0,0,2,abc+1,"C");
    size_t from=0;ChemProduct to={1,1};ChemReaction rxn=reaction(&from,&to,.1);
    ChemConstraintTerm r={"A",1},p={"D",1};ChemConstraint c=builtin("A_D_eq","D",&r,&p,2);
    Fixture a=build(4,abc,1,&rxn,1,&model.model,1,&c,1,0),b=build(4,abc,1,&rxn,1,&model.model,1,&c,1,0);
    Fixture f[]={a,b};
    for (size_t j=0;j<2;++j) {
        CHECK(chem_constraint_count(f[j].m)==2);set(f[j],"A",1);set(f[j],"D",2);
        for (size_t i=0;i<20;++i) step(f[j],1);
    }
    NEAR(2*get(a,"A",0)-get(a,"D",0),0,1e-6);NEAR(5*get(a,"B",0)-get(a,"C",0),0,1e-6);
    for (size_t i=0;i<4;++i) CHECK(get(a,abc[i],0)==get(b,abc[i],0));
    release(a);release(b);
}
static void ownership_and_reorder(void)
{
    /* Additional C ownership checks: source data changes after a build cannot
     * change a model, and every binding is remapped with the mechanism. */
    Fixture a=chemical(1,4,0),b=chemical(1,4,1);
    for (size_t i=0;i<20;++i) { step(a,1);step(b,1); }
    for (size_t i=0;i<4;++i) NEAR(get(a,abc[i],0),get(b,abc[i],0),1e-10);
    release(a);release(b);
    Model model;model_init(&model,EQUILIBRIUM,5,1,0,0,2,abc,"B");
    Fixture f=build(2,abc,0,NULL,1,&model.model,0,NULL,1,0);
    model.data.k=99;set(f,"A",1);double residuals[2]={0};chem_constraints_forcing(f.s,f.s->values,residuals);
    CHECK(residuals[1]==5);release(f);
}
void test_external(void)
{
#define CASE(name,call) RUN("external_model_constraints:ExternalModelConstraints." name,call)
    CASE("AddExternalModelWithConstraints",metadata(1));
    CASE("AddExternalModelProcessOnly",metadata(0));
    CASE("DAESolveEnforcesConservation",conservation(0));
    CASE("AlgebraicRowHoldsOnlyConstraintTerms",algebraic_rows());
    CASE("CombinedBuiltInAndExternalConstraints",combined(0));
    CASE("AddExternalModelOnlyStandardRosenbrock",metadata(0));
    CASE("AddExternalModelConstraintsOnly",conservation(1));
    CASE("MultiGridCell",conservation(2));
    CASE("KineticVsConservativeConstraintConvergence",convergence(0,0));
    CASE("SimpleConstraintPreservesRatioNotMass",convergence(1,0));
    CASE("BuiltInVsExternalModelConstraintStepByStep",stepwise());
    CASE("MultiEquilibriumKineticVsComposedConstraints",convergence(0,1));
    CASE("ProcessJacobianElementInAlgebraicRowSurvivesFiltering",conservation(3));
    CASE("ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering",combined(1));
    CASE("TemperatureDependentConstraintParameter",temperature());
    CASE("ConstraintActiveWhenProcessOnlyModelAddedFirst",model_order());
    CASE("BuilderCanBuildMultipleSolvers",{multiple_builds();ownership_and_reorder();});
#undef CASE
    RUN("external_model_constraints:ExternalModelFiniteDifferenceJacobian.ProcessForcingJacobian",finite_difference(0));
    RUN("external_model_constraints:ExternalModelFiniteDifferenceJacobian.ConstraintResidualJacobian",finite_difference(1));
    RUN("external_model_constraints:ExternalModelFiniteDifferenceJacobian.EquilibriumConstraintModelJacobian",finite_difference(2));
}
