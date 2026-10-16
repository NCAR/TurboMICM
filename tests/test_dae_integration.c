#include "test.h"
#include "chem_process.h"
#include "test_backend.h"
#include <float.h>
typedef struct {
    chem_species species[6]; chem_reaction reactions[5]; size_t reactants[5][2]; chem_product products[5][2];
    chem_constraint constraints[3]; chem_product left[3][6],right[3][2];
    chem_mechanism mechanism; chem_model model; chem_state state; integration_backend backend;
} dae_fixture;
static void dae_init(dae_fixture *f,size_t n,int reorder)
{
    *f=(dae_fixture){0}; const char *names[]={"A","B","C","D","E","F"};
    for (size_t i=0;i<n;++i) f->species[i].name=names[i];
    f->mechanism=(chem_mechanism){.n_species=n,.species=f->species,.reactions=f->reactions,.constraints=f->constraints,.preserve_order=!reorder};
}
static void dae_reaction(dae_fixture *f,size_t nr,const size_t *reactants,size_t np,const chem_product *products,double k,const char *label)
{
    size_t i=f->mechanism.n_reactions++;
    memcpy(f->reactants[i],reactants,nr*sizeof(size_t)); memcpy(f->products[i],products,np*sizeof(chem_product));
    chem_rate rate=chem_rate_default(label ? CHEM_USER : CHEM_ARRHENIUS);
    if (label) rate.label=label; else rate.p.arrhenius.A=k;
    f->reactions[i]=(chem_reaction){nr,np,f->reactants[i],f->products[i],rate};
}
static void dae_equilibrium(dae_fixture *f,const char *name,size_t reactant,double power,size_t product,double K,double dh)
{
    size_t i=f->mechanism.n_constraints++; f->left[i][0]=(chem_product){reactant,power}; f->right[i][0]=(chem_product){product,1};
    f->constraints[i]=(chem_constraint){.name=name,.type=CHEM_EQUILIBRIUM_CONSTRAINT,.algebraic_species=product,
        .n_reactants=1,.n_products=1,.reactants=f->left[i],.products=f->right[i],.equilibrium={.K=K,.delta_H=dh}};
}
static void dae_linear(dae_fixture *f,size_t algebraic,size_t count,const chem_product *terms,double total)
{
    size_t i=f->mechanism.n_constraints++; memcpy(f->left[i],terms,count*sizeof(*terms));
    f->constraints[i]=(chem_constraint){.name="mass_conservation",.type=CHEM_LINEAR_CONSTRAINT,.algebraic_species=algebraic,
        .n_reactants=count,.reactants=f->left[i],.constant=total};
}
static chem_status dae_create(dae_fixture *f,size_t width)
{
    chem_reactions_sort(f->reactions,f->mechanism.n_reactions);
    chem_status status=chem_model_create(&f->model,&f->mechanism,CHEM_LU_DOOLITTLE);
    if (status==CHEM_OK) status=chem_state_create(&f->state,&f->model,1);
    if (status==CHEM_OK) status=integration_backend_create(&f->backend,&f->state,width,CHEM_CSR,0);
    return status;
}
static void dae_destroy(dae_fixture *f)
{ integration_backend_destroy(&f->backend); chem_state_destroy(&f->state); chem_model_destroy(&f->model); }
static double *dae_y(dae_fixture *f,size_t species) { return f->state.y+f->model.spec_to_var[species]; }
static chem_result dae_solve(dae_fixture *f,double dt,chem_ros_method method)
{
    chem_ros_parameters p=chem_ros_parameters_default(method); chem_result r=chem_ros_solve(&f->state,dt,&p,NULL,NULL);
    chem_state_clamp(&f->state); return r;
}
static double eq_constant(double K,double dh,double T)
{ return K*exp((dh/CHEM_GAS_CONSTANT)*(1/T-1/298.15)); }
static void equilibrium_metadata(int multiple)
{
    dae_fixture f; dae_init(&f,multiple ? 6 : 3,0);
    dae_reaction(&f,1,(size_t[]){0},1,(chem_product[]){{1,1}},multiple ? .5 : .1,NULL);
    dae_equilibrium(&f,"B_C_eq",1,1,2,.034,-2400);
    if (multiple) {
        dae_reaction(&f,1,(size_t[]){3},1,(chem_product[]){{4,1}},.2,NULL);
        dae_equilibrium(&f,"E_F_eq",4,1,5,.012,-2000);
    }
    REQUIRE(dae_create(&f,4)==CHEM_OK); size_t n=multiple ? 6 : 3; CHECK(f.model.n_species==n);
    CHECK(f.mechanism.n_constraints==(multiple ? 2u : 1u));
    for (size_t i=0;i<n;++i) { CHECK(chem_species_index(&f.model,f.species[i].name)==i); CHECK(f.model.algebraic[i]==(i==2 || i==5)); }
    CHECK(chem_custom_index(&f.model,"B_C_eq")!=CHEM_NONE);
    memcpy(f.state.y,(double[]){1,.1,.034*.1,.5,.05,.012*.05},n*sizeof(double));
    double T=multiple ? 310 : 300; f.state.conditions[0]=(chem_conditions){T,101325,0}; chem_calculate_rates(&f.state);
    double K=eq_constant(.034,-2400,T); NEAR(f.state.custom[chem_custom_index(&f.model,"B_C_eq")],K,multiple ? 1e-10 : fabs(K)*100*DBL_EPSILON);
    if (multiple) { K=eq_constant(.012,-2000,T); NEAR(f.state.custom[chem_custom_index(&f.model,"E_F_eq")],K,fabs(K)*100*DBL_EPSILON); }
    dae_destroy(&f);
}
TEST(equilibrium_constraint,EquilibriumIntegration,SetConstraintsAPIWorks) { equilibrium_metadata(0); }
TEST(equilibrium_constraint,EquilibriumIntegration,SetConstraintsAPIMultipleConstraints) { equilibrium_metadata(1); }
static void equilibrium_solve(int kind)
{
    dae_fixture f; dae_init(&f,kind==2 ? 4 : 3,kind==1);
    dae_reaction(&f,1,(size_t[]){kind==3 ? 2 : 0},1,(chem_product[]){{kind==3 ? 0 : 1,1}},kind==3 ? .5 : 1,NULL);
    dae_equilibrium(&f,kind==3 ? "A2_B_eq" : "B_C_eq",kind==3 ? 0 : 1,kind==3 ? 2 : 1,kind==3 ? 1 : 2,kind==3 ? 10 : kind==2 ? 3 : 2,-2400);
    if (kind==2) dae_equilibrium(&f,"B_D_eq",1,1,3,5,-2400);
    REQUIRE(dae_create(&f,4)==CHEM_OK);
    if (kind==3) { *dae_y(&f,0)=.1; *dae_y(&f,1)=10*.1*.1; *dae_y(&f,2)=1; }
    else *dae_y(&f,0)=1;
    double T=kind==0 ? 270 : kind==1 ? 400 : kind==2 ? 300 : 298;
    f.state.conditions[0]=(chem_conditions){T,101325,0}; chem_calculate_rates(&f.state);
    for (size_t i=0;i<f.mechanism.n_constraints;++i) {
        double K=eq_constant(f.constraints[i].equilibrium.K,-2400,T);
        NEAR(f.state.custom[f.model.constraint_offsets[i]],K,kind==0 ? fabs(K)*100*DBL_EPSILON : 1e-10);
    }
    double time=0,total=kind==3 ? .05 : .1;
    while (time<total) {
        chem_calculate_rates(&f.state); chem_result r=dae_solve(&f,.001,CHEM_ROS3); REQUIRE(r.status==CHEM_OK);
        for (size_t i=0;i<f.mechanism.n_constraints;++i) {
            const chem_constraint *con=f.constraints+i; double value=*dae_y(&f,con->reactants[0].species);
            double residual=f.state.custom[f.model.constraint_offsets[i]]*value*(kind==3 ? value : 1)-*dae_y(&f,con->products[0].species);
            NEAR(residual,0,kind==3 ? 1e-5 : 1e-6);
        }
        time+=.001;
    }
    if (kind!=3) {
        NEAR(*dae_y(&f,0)+*dae_y(&f,1),1,.01);
        if (kind!=1) for (size_t i=0;i<f.mechanism.n_constraints;++i)
            NEAR(f.state.custom[f.model.constraint_offsets[i]] * *dae_y(&f,1),*dae_y(&f,i+2),1e-6);
    }
    dae_destroy(&f);
}
TEST(equilibrium_constraint,EquilibriumIntegration,DAESolveWithConstraint) { equilibrium_solve(0); }
TEST(equilibrium_constraint,EquilibriumIntegration,DAESolveWithConstraintAndReorderState) { equilibrium_solve(1); }
TEST(equilibrium_constraint,EquilibriumIntegration,DAESolveWithTwoCoupledConstraints) { equilibrium_solve(2); }
TEST(equilibrium_constraint,EquilibriumIntegration,DAESolveWithNonUnitStoichiometry) { equilibrium_solve(3); }
static void overshoot_case(chem_ros_method method,int equilibrium)
{
    dae_fixture f; dae_init(&f,3,0);
    dae_reaction(&f,1,(size_t[]){equilibrium ? 1 : 0},1,(chem_product[]){{equilibrium ? 2 : 1,1}},equilibrium ? 1e3 : 1e4,NULL);
    if (equilibrium) {
        dae_equilibrium(&f,"gas_aq_eq",0,1,1,10,0);
        dae_linear(&f,0,3,(chem_product[]){{1,1},{2,1},{0,1}},1e-6);
    } else dae_linear(&f,2,3,(chem_product[]){{0,1},{1,1},{2,1}},1e-6);
    REQUIRE(dae_create(&f,4)==CHEM_OK); chem_state *s=&f.state; s->rtol=1e-6;
    for (size_t i=0;i<3;++i) s->atol[i]=equilibrium && i<2 ? 1e-8 : 1e-12;
    if (equilibrium) { s->y[0]=1e-6/11; s->y[1]=10*s->y[0]; }
    else { s->y[0]=.9e-6; s->y[2]=.1e-6; }
    s->conditions[0]=(chem_conditions){298,101325,0}; chem_calculate_rates(s); double advanced=0;
    while (advanced<30) {
        chem_result r=dae_solve(&f,30-advanced,method); REQUIRE(r.status==CHEM_OK); REQUIRE(advanced+r.stats.final_time>advanced); advanced+=r.stats.final_time;
        NEAR(s->y[0]+s->y[1]+s->y[2],1e-6,1e-12);
        CHECK(s->y[equilibrium ? 0 : 2]>=-1e-18); if (equilibrium) CHECK(s->y[1]>=-1e-18);
    }
    if (equilibrium) { NEAR(s->y[2],1e-6,1e-8); CHECK(s->y[0]>=-1e-18); }
    else { CHECK(s->y[0]<1e-12); NEAR(s->y[1],.9e-6,1e-10); NEAR(s->y[2],.1e-6,1e-10); }
    dae_destroy(&f);
}
TEST(dae_constraint_overshoot,DAEConstraintOvershoot,AlgebraicVariableStaysNonNegative) { overshoot_case(CHEM_RODA3,0); }
TEST(dae_constraint_overshoot,DAEConstraintOvershoot,EquilibriumPlusConservation) { overshoot_case(CHEM_RODA3,1); }
TEST(dae_constraint_overshoot,DAEConstraintOvershoot,AllRosenbrockOrdersConstrained)
{ for (int method=CHEM_ROS3;method<=CHEM_RODA4;++method) overshoot_case((chem_ros_method)method,0); }
typedef struct { double gas,product,minimum; uint64_t accepted,rejected; } cascade_result;
static cascade_result cascade(double balance_atol,double k,double K1,double K2)
{
    cascade_result out={.minimum=1e99}; dae_fixture f; dae_init(&f,4,0);
    dae_reaction(&f,1,(size_t[]){2},1,(chem_product[]){{3,1}},k,NULL);
    dae_equilibrium(&f,"eq1",0,1,1,K1,0); dae_equilibrium(&f,"eq2",1,1,2,K2,0);
    dae_linear(&f,0,4,(chem_product[]){{1,1},{2,1},{3,1},{0,1}},1e-6);
    chem_status status=dae_create(&f,0); CHECK(status==CHEM_OK); if (status!=CHEM_OK) return out;
    chem_state *s=&f.state; s->rtol=1e-6; memcpy(s->atol,(double[]){balance_atol,1e-8,1e-8,1e-12},4*sizeof(double));
    s->y[0]=1e-6/(1+K1+K1*K2); s->y[1]=K1*s->y[0]; s->y[2]=K2*K1*s->y[0];
    s->conditions[0]=(chem_conditions){298,101325,0}; chem_calculate_rates(s);
    double advanced=0;
    while (advanced<30) {
        chem_result r=dae_solve(&f,30-advanced,CHEM_RODA3); CHECK(r.status==CHEM_OK);
        if (r.status!=CHEM_OK || advanced+r.stats.final_time<=advanced) break;
        advanced+=r.stats.final_time; out.accepted+=r.stats.accepted; out.rejected+=r.stats.rejected;
        out.minimum=fmin(out.minimum,s->y[0]);
    }
    out.gas=s->y[0]; out.product=s->y[3]; dae_destroy(&f); return out;
}
TEST(dae_algebraic_error_insensitivity,DAEAlgebraicError,ErrorSensitiveToBalanceAtol)
{
    cascade_result loose=cascade(1e-3,100,2,2),tight=cascade(1e-8,100,2,2); CHECK(tight.accepted>loose.accepted);
}
TEST(dae_algebraic_error_insensitivity,DAEAlgebraicError,AlgebraicVariableDoesNotOvershootDeeply)
{
    cascade_result r=cascade(1e-8,1e4,100,50); CHECK(r.minimum>=-1e-8); NEAR(r.product,1e-6,1e-8); NEAR(r.gas,0,1e-10);
}
TEST(linear_constraint,DAESolveWithConstraint,TerminatorAndRobertson)
{
    dae_fixture f; dae_init(&f,5,0); f.species[0].absolute_tolerance=f.species[1].absolute_tolerance=1e-20;
    dae_reaction(&f,1,(size_t[]){0},1,(chem_product[]){{1,2}},0,"terminator_k1");
    dae_reaction(&f,2,(size_t[]){1,1},1,(chem_product[]){{0,1}},1,NULL);
    dae_reaction(&f,1,(size_t[]){2},1,(chem_product[]){{3,1}},0,"robertson_r1");
    dae_reaction(&f,2,(size_t[]){3,3},2,(chem_product[]){{3,1},{4,1}},0,"robertson_r2");
    dae_reaction(&f,2,(size_t[]){3,4},2,(chem_product[]){{2,1},{4,1}},0,"robertson_r3");
    dae_linear(&f,4,3,(chem_product[]){{2,1},{3,1},{4,1}},1);
    REQUIRE(dae_create(&f,4)==CHEM_OK); chem_state *s=&f.state; s->rtol=1e-8;
    REQUIRE(chem_state_set_custom(s,"robertson_r1",(double[]){.04},1)==CHEM_OK);
    REQUIRE(chem_state_set_custom(s,"robertson_r2",(double[]){3e7},1)==CHEM_OK);
    REQUIRE(chem_state_set_custom(s,"robertson_r3",(double[]){1e4},1)==CHEM_OK);
    memcpy(s->y,(double[]){1.8e-10,1.2e-6,1,0,0},5*sizeof(double)); s->conditions[0]=(chem_conditions){298,101300,42}; chem_calculate_rates(s);
    double chlorine=2*1.8e-10+1.2e-6,dt=1;
    for (size_t step=0;step<12;++step) {
        double advanced=0;
        while (advanced<dt) {
            chem_result r=dae_solve(&f,dt-advanced,CHEM_RODA3); REQUIRE(r.status==CHEM_OK || r.status==CHEM_MAX_STEPS);
            REQUIRE(advanced+r.stats.final_time>advanced); advanced+=r.stats.final_time;
        }
        NEAR(s->y[2]+s->y[3]+s->y[4],1,1e-10); REL(2*s->y[0]+s->y[1],chlorine,1e-6);
        for (size_t i=0;i<5;++i) { CHECK(isfinite(s->y[i])); CHECK(s->y[i]>=0); }
        if (!step) { CHECK(s->y[2]<.99); CHECK(s->y[3]>0); CHECK(s->y[4]>1e-3); }
        dt*=10;
    }
    CHECK(s->y[4]>.95); CHECK(s->y[2]<.05); CHECK(s->y[0]>.4*chlorine); CHECK(s->y[1]<.1*chlorine); dae_destroy(&f);
}
