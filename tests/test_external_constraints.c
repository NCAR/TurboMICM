#include "test.h"
#include "chem_external.h"
#include "chem_process.h"
#include "chem_verify.h"
#include "test_backend.h"

typedef struct {
    const char *name; size_t n; const char *terms[4]; double coefficients[4],constant;
    int constant_first; size_t row,ids[4],flat[4];
} external_row;
typedef struct {
    size_t n_rows; external_row rows[3];
    const char *reactant,*product,*solvent; double k;
    size_t r,p,v,rr,pr,rv,pv;
    int temperature_dependent,finalized; size_t parameter;
    double K,dh_over_R,reference_temperature;
} external_context;
static chem_status external_finalize(const chem_model *m,void *context)
{
    external_context *e=context; ++e->finalized;
    for (size_t i=0;i<e->n_rows;++i) {
        external_row *row=e->rows+i; row->row=chem_species_index(m,row->name); if (row->row==CHEM_NONE) return CHEM_NOT_FOUND;
        for (size_t j=0;j<row->n;++j) {
            row->ids[j]=chem_species_index(m,row->terms[j]); if (row->ids[j]==CHEM_NONE) return CHEM_NOT_FOUND;
            row->flat[j]=m->lu.map[row->row*m->n_species+row->ids[j]]; if (row->flat[j]==CHEM_NONE) return CHEM_NOT_FOUND;
        }
    }
    if (e->reactant) {
        e->r=chem_species_index(m,e->reactant); e->p=chem_species_index(m,e->product);
        if (e->r==CHEM_NONE || e->p==CHEM_NONE) return CHEM_NOT_FOUND;
        e->rr=m->lu.map[e->r*m->n_species+e->r]; e->pr=m->lu.map[e->p*m->n_species+e->r];
        if (e->rr==CHEM_NONE || e->pr==CHEM_NONE) return CHEM_NOT_FOUND;
        if (e->solvent) {
            e->v=chem_species_index(m,e->solvent); if (e->v==CHEM_NONE) return CHEM_NOT_FOUND;
            e->rv=m->lu.map[e->r*m->n_species+e->v]; e->pv=m->lu.map[e->p*m->n_species+e->v];
            if (e->rv==CHEM_NONE || e->pv==CHEM_NONE) return CHEM_NOT_FOUND;
        }
    }
    if (e->temperature_dependent) { e->parameter=chem_custom_index(m,"C_K_eq"); if (e->parameter==CHEM_NONE) return CHEM_NOT_FOUND; }
    return CHEM_OK;
}
static void external_process(const chem_state *s,const double *y,double *out,void *context)
{
    const external_context *e=context; if (!e->reactant) return;
    for (size_t c=0;c<s->cells;++c) {
        double rate=e->k*y[e->r*s->cells+c]; if (e->solvent) rate*=y[e->v*s->cells+c];
        out[e->r*s->cells+c]-=rate; out[e->p*s->cells+c]+=rate;
    }
}
static void external_process_jacobian(const chem_state *s,const double *y,double *out,void *context)
{
    const external_context *e=context; if (!e->reactant) return;
    for (size_t c=0;c<s->cells;++c) {
        double rate=e->solvent ? e->k*y[e->v*s->cells+c] : e->k;
        out[e->rr*s->cells+c]+=rate; out[e->pr*s->cells+c]-=rate;
        if (e->solvent) { rate=e->k*y[e->r*s->cells+c]; out[e->rv*s->cells+c]+=rate; out[e->pv*s->cells+c]-=rate; }
    }
}
static void external_residual(const chem_state *s,const double *y,double *out,void *context)
{
    const external_context *e=context;
    for (size_t i=0;i<e->n_rows;++i) for (size_t c=0;c<s->cells;++c) {
        const external_row *r=e->rows+i; double value=r->constant_first ? r->constant : 0;
        for (size_t j=0;j<r->n;++j) {
            double coefficient=e->temperature_dependent && i==0 && j==0 ? s->custom[e->parameter*s->cells+c] : r->coefficients[j];
            value+=coefficient*y[r->ids[j]*s->cells+c];
        }
        if (!r->constant_first) value+=r->constant;
        out[r->row*s->cells+c]=value;
    }
}
static void external_constraint_jacobian(const chem_state *s,const double *y,double *out,void *context)
{
    (void)y; const external_context *e=context;
    for (size_t i=0;i<e->n_rows;++i) for (size_t j=0;j<e->rows[i].n;++j) for (size_t c=0;c<s->cells;++c) {
        double coefficient=e->temperature_dependent && i==0 && j==0 ? s->custom[e->parameter*s->cells+c] : e->rows[i].coefficients[j];
        out[e->rows[i].flat[j]*s->cells+c]-=coefficient;
    }
}
static void external_update(chem_state *s,void *context)
{
    const external_context *e=context; if (!e->temperature_dependent) return;
    for (size_t c=0;c<s->cells;++c) s->custom[e->parameter*s->cells+c]=e->K*exp(e->dh_over_R*(1/e->reference_temperature-1/s->conditions[c].temperature));
}
typedef struct {
    external_context context; chem_species species; chem_named_value parameter;
    chem_named_element process[4],constraints[12]; const char *algebraic[3];
} external_spec;
typedef struct {
    chem_species species[4]; chem_reaction reactions[5]; size_t reactants[5]; chem_product products[5];
    chem_constraint constraint; chem_product left,right;
    external_spec specs[4]; chem_external_model external[4];
    chem_mechanism mechanism; chem_model model; chem_state state; integration_backend backend;
} external_fixture;
static void external_init(external_fixture *f,size_t n,const char *const *names)
{
    *f=(external_fixture){0}; for (size_t i=0;i<n;++i) f->species[i].name=names[i];
    f->mechanism=(chem_mechanism){.n_species=n,.species=f->species,.reactions=f->reactions,.preserve_order=1,.external=f->external,.constraints=&f->constraint};
}
static void external_reaction(external_fixture *f,size_t reactant,size_t product,double k)
{
    size_t i=f->mechanism.n_reactions++; f->reactants[i]=reactant; f->products[i]=(chem_product){product,1};
    chem_rate rate=chem_rate_default(CHEM_ARRHENIUS); rate.p.arrhenius.A=k;
    f->reactions[i]=(chem_reaction){1,1,f->reactants+i,f->products+i,rate};
}
static external_context *external_add(external_fixture *f)
{ return &f->specs[f->mechanism.n_external++].context; }
static void external_equilibrium(external_context *e,const char *r,const char *p,double K)
{ e->rows[e->n_rows++]=(external_row){.name=p,.n=2,.terms={r,p},.coefficients={K,-1}}; }
static void external_mass(external_context *e,const char *row,size_t n,const char *const *names,int first)
{
    external_row *r=e->rows+e->n_rows++; *r=(external_row){.name=row,.n=n,.constant=-1,.constant_first=first};
    for (size_t i=0;i<n;++i) { r->terms[i]=names[i]; r->coefficients[i]=1; }
}
static void external_aerosol(external_fixture *f,double k,int mass,int solvent)
{
    size_t i=f->mechanism.n_external; external_context *e=external_add(f);
    e->reactant="A_GAS"; e->product="AEROSOL.A_AQ"; e->k=k; e->solvent=solvent ? "S" : NULL;
    f->specs[i].species.name="AEROSOL.A_AQ";
    if (mass) external_mass(e,"AEROSOL.A_AQ",2,(const char *const[]){"A_GAS","AEROSOL.A_AQ"},0);
}
static void external_builtin(external_fixture *f,size_t r,size_t p,double K,const char *name)
{
    f->left=(chem_product){r,1}; f->right=(chem_product){p,1}; f->mechanism.n_constraints=1;
    f->constraint=(chem_constraint){.name=name,.type=CHEM_EQUILIBRIUM_CONSTRAINT,.algebraic_species=p,.n_reactants=1,.n_products=1,
        .reactants=&f->left,.products=&f->right,.equilibrium={.K=K}};
}
static void external_declare(external_fixture *f)
{
    for (size_t i=0;i<f->mechanism.n_external;++i) {
        external_spec *s=f->specs+i; external_context *e=&s->context;
        chem_external_model *d=f->external+i;
        *d=(chem_external_model){.n_species=s->species.name ? 1 : 0,.species=&s->species,.n_algebraic=e->n_rows,.algebraic=s->algebraic,
            .process_elements=s->process,.constraint_elements=s->constraints,.finalize=external_finalize,.update=external_update,
            .context=e,.context_size=sizeof(*e),.forcing=external_process,.jacobian=external_process_jacobian,
            .residual=external_residual,.constraint_jacobian=external_constraint_jacobian};
        if (e->reactant) {
            s->process[d->n_process_elements++]=(chem_named_element){e->reactant,e->reactant};
            s->process[d->n_process_elements++]=(chem_named_element){e->product,e->reactant};
            if (e->solvent) {
                s->process[d->n_process_elements++]=(chem_named_element){e->reactant,e->solvent};
                s->process[d->n_process_elements++]=(chem_named_element){e->product,e->solvent};
            }
        }
        for (size_t j=0;j<e->n_rows;++j) {
            s->algebraic[j]=e->rows[j].name;
            for (size_t k=0;k<e->rows[j].n;++k) s->constraints[d->n_constraint_elements++]=(chem_named_element){e->rows[j].name,e->rows[j].terms[k]};
        }
        if (e->temperature_dependent) { s->parameter=(chem_named_value){"C_K_eq",0}; d->n_parameters=1; d->parameters=&s->parameter; }
    }
}
static chem_status external_create(external_fixture *f,size_t cells)
{
    external_declare(f); chem_status status=chem_model_create(&f->model,&f->mechanism,CHEM_LU_DOOLITTLE);
    if (status==CHEM_OK) status=chem_state_create(&f->state,&f->model,cells);
    if (status==CHEM_OK) status=integration_backend_create(&f->backend,&f->state,0,CHEM_CSR,0);
    if (status==CHEM_OK) for (size_t c=0;c<cells;++c) f->state.conditions[c]=(chem_conditions){298,101325,0};
    return status;
}
static void external_destroy(external_fixture *f)
{ integration_backend_destroy(&f->backend); chem_state_destroy(&f->state); chem_model_destroy(&f->model); }
static chem_result external_solve(external_fixture *f,double dt,chem_ros_method method)
{
    chem_calculate_rates(&f->state); chem_ros_parameters p=chem_ros_parameters_default(method);
    chem_result result=chem_ros_solve(&f->state,dt,&p,NULL,NULL); chem_state_clamp(&f->state); return result;
}
static void aerosol_metadata(int mass)
{
    external_fixture f; external_init(&f,1,(const char *const[]){"A_GAS"}); external_aerosol(&f,.01,mass,0);
    REQUIRE(external_create(&f,1)==CHEM_OK); CHECK(f.model.n_species==2); CHECK(f.model.mechanism->n_constraints==(size_t)mass);
    CHECK(f.model.algebraic[0]==0); CHECK(f.model.algebraic[1]==mass); external_destroy(&f);
}
TEST(external_model_constraints,ExternalModelConstraints,AddExternalModelWithConstraints) { aerosol_metadata(1); }
TEST(external_model_constraints,ExternalModelConstraints,AddExternalModelProcessOnly) { aerosol_metadata(0); }
TEST(external_model_constraints,ExternalModelConstraints,AddExternalModelOnlyStandardRosenbrock) { aerosol_metadata(0); }
static void aerosol_conservation(size_t cells,size_t steps,double dt,double k,double initial,int solvent)
{
    external_fixture f; external_init(&f,solvent ? 2 : 1,(const char *const[]){"A_GAS","S"}); external_aerosol(&f,k,1,solvent);
    REQUIRE(external_create(&f,cells)==CHEM_OK); size_t aq=chem_species_index(&f.model,"AEROSOL.A_AQ"); CHECK(f.model.algebraic[aq]==1);
    for (size_t c=0;c<cells;++c) { f.state.y[c]=initial-.2*c; f.state.y[aq*cells+c]=1-f.state.y[c]; if (solvent) f.state.y[cells+c]=1; }
    for (size_t step=0;step<steps;++step) {
        REQUIRE(external_solve(&f,dt,CHEM_RODA3).status==CHEM_OK);
        for (size_t c=0;c<cells;++c) NEAR(f.state.y[c]+f.state.y[aq*cells+c],1,1e-4);
    }
    external_destroy(&f);
}
TEST(external_model_constraints,ExternalModelConstraints,DAESolveEnforcesConservation) { aerosol_conservation(1,20,10,.1,.9,0); }
TEST(external_model_constraints,ExternalModelConstraints,AddExternalModelConstraintsOnly) { aerosol_conservation(1,1,50,.01,.8,0); }
TEST(external_model_constraints,ExternalModelConstraints,MultiGridCell) { aerosol_conservation(3,1,50,.1,.9,0); }
TEST(external_model_constraints,ExternalModelConstraints,ProcessJacobianElementInAlgebraicRowSurvivesFiltering) { aerosol_conservation(1,10,10,.1,.9,1); }
TEST(external_model_constraints,ExternalModelConstraints,AlgebraicRowHoldsOnlyConstraintTerms)
{
    external_fixture f; external_init(&f,1,(const char *const[]){"A_GAS"}); external_aerosol(&f,10,1,0); REQUIRE(external_create(&f,1)==CHEM_OK);
    f.state.y[0]=.7; f.state.y[1]=.2; double force[2]={0};
    chem_state_forcing(&f.state,f.state.y,force); chem_constraint_forcing(&f.state,f.state.y,force);
    NEAR(force[0],-10*.7,1e-12); NEAR(force[1],.7+.2-1,1e-12);
    chem_state_jacobian(&f.state,f.state.y,f.state.jacobian); chem_constraint_jacobian(&f.state,f.state.y,f.state.jacobian);
    NEAR(f.state.jacobian[f.model.lu.map[0]],10,1e-12); NEAR(f.state.jacobian[f.model.lu.map[2]],-1,1e-12); NEAR(f.state.jacobian[f.model.lu.map[3]],-1,1e-12);
    /* Also preserve an existing value when evaluating only kinetic terms. */
    force[0]=force[1]=123; chem_state_forcing(&f.state,f.state.y,force); NEAR(force[1],123,0);
    external_destroy(&f);
}
TEST(external_model_constraints,ExternalModelConstraints,CombinedBuiltInAndExternalConstraints)
{
    external_fixture f; external_init(&f,3,(const char *const[]){"A_GAS","B","C"}); external_reaction(&f,0,1,.05);
    external_aerosol(&f,.01,1,0); external_builtin(&f,1,2,5,"B_C_eq"); REQUIRE(external_create(&f,1)==CHEM_OK);
    CHECK(f.model.mechanism->n_constraints==2); CHECK(f.model.algebraic[2]==1); CHECK(f.model.algebraic[3]==1); external_destroy(&f);
}
TEST(external_model_constraints,ExternalModelConstraints,ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering)
{
    external_fixture f; external_init(&f,2,(const char *const[]){"A_GAS","B"}); external_aerosol(&f,.1,0,0);
    external_builtin(&f,1,2,5,"B_AQ_eq"); REQUIRE(external_create(&f,1)==CHEM_OK);
    CHECK(f.model.mechanism->n_constraints==1); CHECK(f.model.algebraic[2]==1); CHECK(f.model.lu.map[2*3]!=CHEM_NONE); external_destroy(&f);
}
static void convergence_system(int kind,double *out)
{
    external_fixture f; external_init(&f,3,(const char *const[]){"A","B","C"}); external_reaction(&f,0,1,.1);
    if (!kind) { external_reaction(&f,1,2,10); external_reaction(&f,2,1,2); }
    else {
        external_context *e=external_add(&f);
        if (kind==1) external_mass(e,"B",3,(const char *const[]){"A","B","C"},0);
        external_equilibrium(e,"B","C",5);
    }
    REQUIRE(external_create(&f,1)==CHEM_OK); f.state.y[0]=1;
    for (size_t step=0;step<200;++step) REQUIRE(external_solve(&f,1,kind ? CHEM_RODA3 : CHEM_ROS3).status==CHEM_OK);
    memcpy(out,f.state.y,3*sizeof(double)); external_destroy(&f);
}
TEST(external_model_constraints,ExternalModelConstraints,KineticVsConservativeConstraintConvergence)
{
    double kin[3]={0},con[3]={0},expected[]={0,1.0/6,5.0/6}; convergence_system(0,kin); convergence_system(1,con);
    for (size_t i=0;i<3;++i) { NEAR(kin[i],expected[i],i ? 1e-3 : 1e-6); NEAR(con[i],expected[i],i ? 1e-3 : 1e-6); NEAR(kin[i],con[i],i ? 1e-3 : 1e-6); }
    NEAR(kin[0]+kin[1]+kin[2],1,1e-6); NEAR(con[0]+con[1]+con[2],1,1e-3);
}
TEST(external_model_constraints,ExternalModelConstraints,SimpleConstraintPreservesRatioNotMass)
{
    double kin[3]={0},simple[3]={0}; convergence_system(0,kin); convergence_system(2,simple);
    NEAR(kin[2]/kin[1],5,1e-3); NEAR(simple[2]/simple[1],5,1e-3);
    NEAR(kin[0]+kin[1]+kin[2],1,1e-6); NEAR(simple[0]+simple[1]+simple[2],6,.1);
}
TEST(external_model_constraints,ExternalModelConstraints,BuiltInVsExternalModelConstraintStepByStep)
{
    external_fixture f[2];
    for (size_t i=0;i<2;++i) {
        external_init(f+i,3,(const char *const[]){"A","B","C"}); external_reaction(f+i,0,1,.1);
        if (i) external_equilibrium(external_add(f+i),"B","C",5); else external_builtin(f+i,1,2,5,"B_C_eq");
        REQUIRE(external_create(f+i,1)==CHEM_OK); f[i].state.y[0]=1;
    }
    for (size_t step=0;step<100;++step) {
        for (size_t i=0;i<2;++i) { REQUIRE(external_solve(f+i,.5,CHEM_RODA3).status==CHEM_OK); NEAR(5*f[i].state.y[1]-f[i].state.y[2],0,1e-6); }
        for (size_t i=0;i<3;++i) NEAR(f[0].state.y[i],f[1].state.y[i],1e-8);
    }
    for (size_t i=0;i<2;++i) external_destroy(f+i);
}
TEST(external_model_constraints,ExternalModelConstraints,MultiEquilibriumKineticVsComposedConstraints)
{
    external_fixture f[2]; const char *names[]={"A","B","C","D"};
    for (size_t i=0;i<2;++i) {
        external_init(f+i,4,names); external_reaction(f+i,0,1,.1);
        if (!i) { external_reaction(f+i,1,2,10); external_reaction(f+i,2,1,2); external_reaction(f+i,1,3,8); external_reaction(f+i,3,1,8.0/3); }
        else { external_equilibrium(external_add(f+i),"B","C",5); external_equilibrium(external_add(f+i),"B","D",3); external_mass(external_add(f+i),"B",4,names,1); }
        REQUIRE(external_create(f+i,1)==CHEM_OK); f[i].state.y[0]=1;
    }
    for (size_t step=0;step<200;++step) for (size_t i=0;i<2;++i) REQUIRE(external_solve(f+i,1,i ? CHEM_RODA3 : CHEM_ROS3).status==CHEM_OK);
    for (size_t i=0;i<2;++i) {
        double *y=f[i].state.y; NEAR(y[0],0,1e-6); NEAR(y[1],1.0/9,1e-3); NEAR(y[2],5.0/9,1e-3); NEAR(y[3],3.0/9,1e-3);
        NEAR(y[0]+y[1]+y[2]+y[3],1,i ? 1e-3 : 1e-6);
    }
    for (size_t i=1;i<4;++i) NEAR(f[0].state.y[i],f[1].state.y[i],1e-3);
    NEAR(5*f[1].state.y[1]-f[1].state.y[2],0,1e-6); NEAR(3*f[1].state.y[1]-f[1].state.y[3],0,1e-6);
    for (size_t i=0;i<2;++i) external_destroy(f+i);
}
TEST(external_model_constraints,ExternalModelConstraints,TemperatureDependentConstraintParameter)
{
    external_fixture f; external_init(&f,3,(const char *const[]){"A","B","C"}); external_reaction(&f,0,1,.1);
    external_context *e=external_add(&f); external_equilibrium(e,"B","C",2);
    e->temperature_dependent=1; e->K=2; e->dh_over_R=3000; e->reference_temperature=298.15;
    REQUIRE(external_create(&f,1)==CHEM_OK);
    for (size_t run=0;run<2;++run) {
        if (run) {
            integration_backend_destroy(&f.backend); chem_state_destroy(&f.state);
            REQUIRE(chem_state_create(&f.state,&f.model,1)==CHEM_OK); REQUIRE(integration_backend_create(&f.backend,&f.state,0,CHEM_CSR,0)==CHEM_OK);
        }
        double temperature=run ? 350 : 298.15; f.state.conditions[0]=(chem_conditions){temperature,101325,0}; f.state.y[0]=1;
        for (size_t step=0;step<200;++step) REQUIRE(external_solve(&f,1,CHEM_RODA3).status==CHEM_OK);
        double K=2*exp(3000*(1/298.15-1/temperature)); CHECK(f.state.y[1]>0); if (run) CHECK(K>2);
        NEAR(f.state.y[2]/f.state.y[1],K,1e-4);
        NEAR(f.state.custom[chem_custom_index(&f.model,"C_K_eq")],K,0);
    }
    external_destroy(&f);
}
TEST(external_model_constraints,ExternalModelConstraints,ConstraintActiveWhenProcessOnlyModelAddedFirst)
{
    double results[2][2];
    for (size_t first=0;first<2;++first) {
        external_fixture f; external_init(&f,3,(const char *const[]){"A","B","C"}); external_reaction(&f,0,1,.1);
        if (first) external_add(&f);
        external_equilibrium(external_add(&f),"B","C",5);
        if (!first) external_add(&f);
        REQUIRE(external_create(&f,1)==CHEM_OK); CHECK(f.model.mechanism->n_constraints==1); CHECK(f.model.algebraic[2]==1); f.state.y[0]=1;
        for (size_t step=0;step<20;++step) REQUIRE(external_solve(&f,1,CHEM_RODA3).status==CHEM_OK);
        results[first][0]=f.state.y[1]; results[first][1]=f.state.y[2]; CHECK(f.state.y[1]>0); NEAR(5*f.state.y[1]-f.state.y[2],0,1e-6);
        external_destroy(&f);
    }
    for (size_t i=0;i<2;++i) NEAR(results[0][i],results[1][i],0);
}
TEST(external_model_constraints,ExternalModelConstraints,BuilderCanBuildMultipleSolvers)
{
    external_fixture f; external_init(&f,4,(const char *const[]){"A","B","C","D"}); external_reaction(&f,0,1,.1);
    external_builtin(&f,0,3,2,"A_D_eq"); external_equilibrium(external_add(&f),"B","C",5);
    REQUIRE(external_create(&f,1)==CHEM_OK); chem_model second={0}; chem_state state={0}; integration_backend backend={0};
    REQUIRE(chem_model_create(&second,&f.mechanism,CHEM_LU_DOOLITTLE)==CHEM_OK);
    REQUIRE(chem_state_create(&state,&second,1)==CHEM_OK); REQUIRE(integration_backend_create(&backend,&state,0,CHEM_CSR,0)==CHEM_OK);
    CHECK(f.specs[0].context.finalized==0);
    const external_context *c1=chem_model_external_context(&f.model,0),*c2=chem_model_external_context(&second,0);
    REQUIRE(c1 && c2); CHECK(c1!=c2); CHECK(c1->finalized==1); CHECK(c2->finalized==1);
    chem_state *states[]={&f.state,&state}; chem_ros_parameters p=chem_ros_parameters_default(CHEM_RODA3);
    for (size_t i=0;i<2;++i) {
        chem_state *s=states[i]; CHECK(s->model->mechanism->n_constraints==2); s->y[0]=1; s->y[3]=2; s->conditions[0]=(chem_conditions){298,101325,0};
        for (size_t step=0;step<20;++step) { chem_calculate_rates(s); chem_result r=chem_ros_solve(s,1,&p,NULL,NULL); chem_state_clamp(s); REQUIRE(r.status==CHEM_OK); }
    }
    NEAR(2*f.state.y[0]-f.state.y[3],0,1e-6); NEAR(5*f.state.y[1]-f.state.y[2],0,1e-6);
    for (size_t i=0;i<4;++i) NEAR(f.state.y[i],state.y[i],0);
    integration_backend_destroy(&backend); chem_state_destroy(&state); chem_model_destroy(&second); external_destroy(&f);
}
typedef struct { const chem_state *state; const external_context *context; int constraint; } external_fd_context;
static void external_fd(void *context,const double *y,double *forcing)
{
    const external_fd_context *c=context;
    if (c->constraint) external_residual(c->state,y,forcing,(void *)c->context);
    else external_process(c->state,y,forcing,(void *)c->context);
}
static void external_finite_difference(int kind)
{
    external_fixture f;
    if (kind==2) { external_init(&f,2,(const char *const[]){"A","B"}); external_equilibrium(external_add(&f),"A","B",2.5); }
    else { external_init(&f,1,(const char *const[]){"A_GAS"}); external_aerosol(&f,kind ? .1 : 3.5,kind!=0,0); }
    size_t cells=kind==2 ? 1 : 2; REQUIRE(external_create(&f,cells)==CHEM_OK);
    if (kind==2) memcpy(f.state.y,(double[]){3,5},2*sizeof(double));
    else if (kind) memcpy(f.state.y,(double[]){.6,.2,.4,.8},4*sizeof(double));
    else memcpy(f.state.y,(double[]){.8,.3,.2,.7},4*sizeof(double));
    const external_context *context=chem_model_external_context(&f.model,0);
    if (kind) external_constraint_jacobian(&f.state,f.state.y,f.state.jacobian,(void *)context);
    else external_process_jacobian(&f.state,f.state.y,f.state.jacobian,(void *)context);
    external_fd_context wrapper={&f.state,context,kind!=0}; double fd[8];
    REQUIRE(chem_finite_difference(2,cells,f.state.y,external_fd,&wrapper,1e-8,fd)==CHEM_OK);
    /* Compare and check completeness against the external model's declared
     * sparsity, without solver-added diagonal entries or LU fill. */
    size_t map[4]={CHEM_NONE,CHEM_NONE,CHEM_NONE,CHEM_NONE};
    const chem_external_model *d=&f.external[0]; size_t count=kind ? d->n_constraint_elements : d->n_process_elements;
    const chem_named_element *elements=kind ? d->constraint_elements : d->process_elements;
    for (size_t i=0;i<count;++i) {
        size_t row=chem_species_index(&f.model,elements[i].row),col=chem_species_index(&f.model,elements[i].column);
        map[row*2+col]=f.model.lu.map[row*2+col];
    }
    CHECK(chem_compare_jacobian(2,cells,map,f.state.jacobian,fd,1e-7,1e-7).passed);
    if (kind!=1) CHECK(chem_check_sparsity(2,cells,map,fd,1e-6).passed);
    external_destroy(&f);
}
TEST(external_model_constraints,ExternalModelFiniteDifferenceJacobian,ProcessForcingJacobian) { external_finite_difference(0); }
TEST(external_model_constraints,ExternalModelFiniteDifferenceJacobian,ConstraintResidualJacobian) { external_finite_difference(1); }
TEST(external_model_constraints,ExternalModelFiniteDifferenceJacobian,EquilibriumConstraintModelJacobian) { external_finite_difference(2); }
