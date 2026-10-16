/* MICM test/unit/solver/test_constraint_initialization.cpp, C equivalents. */
#include "test.h"

typedef struct {
    chem_species species[3]; size_t reactants[1]; chem_product products[1], left[3], right[1];
    chem_reaction reaction; chem_constraint constraint; chem_mechanism mechanism;
    chem_model model; chem_state state;
} init_fixture;

/* kind: 0=A->B, C=2B; 1=2B-C^2=0; 2=scale*(X-Z)=0; 3=A+B+C=1.3; 4=pure ODE. */
static chem_status init_create(init_fixture *f, int kind, double scale, size_t cells)
{
    memset(f,0,sizeof(*f));
    f->species[0].name="A"; f->species[1].name="B"; f->species[2].name="C";
    f->reactants[0]=0; f->products[0]=(chem_product){1,1};
    f->reaction=(chem_reaction){.n_reactants=1,.n_products=1,.reactants=f->reactants,.products=f->products};
    f->reaction.rate=chem_rate_default(CHEM_ARRHENIUS); f->reaction.rate.p.arrhenius.A=.5;
    f->left[0]=(chem_product){1,1}; f->right[0]=(chem_product){2,1};
    f->constraint=(chem_constraint){.name="B_C_eq",.type=CHEM_EQUILIBRIUM_CONSTRAINT,.algebraic_species=2,
        .n_reactants=1,.n_products=1,.reactants=f->left,.products=f->right,.equilibrium={.K=2}};
    f->mechanism=(chem_mechanism){.n_species=3,.species=f->species,.n_reactions=1,.reactions=&f->reaction,
        .preserve_order=1,.n_constraints=1,.constraints=&f->constraint};
    if (kind==1) {
        f->mechanism.n_species=2; f->mechanism.n_reactions=0;
        f->left[0]=(chem_product){0,1}; f->right[0]=(chem_product){1,2}; f->constraint.algebraic_species=1;
    } else if (kind==2) {
        f->mechanism.n_species=2; f->mechanism.n_reactions=0;
        f->left[0]=(chem_product){0,scale}; f->left[1]=(chem_product){1,-scale};
        f->constraint.type=CHEM_LINEAR_CONSTRAINT; f->constraint.algebraic_species=1;
        f->constraint.n_reactants=2; f->constraint.n_products=0;
    } else if (kind==3) {
        for (size_t i=0; i<3; ++i) f->left[i]=(chem_product){i,scale};
        f->constraint.type=CHEM_LINEAR_CONSTRAINT; f->constraint.n_reactants=3;
        f->constraint.n_products=0; f->constraint.constant=scale*1.3;
    } else if (kind==4) { f->mechanism.n_species=2; f->mechanism.n_constraints=0; }
    chem_status rc=chem_model_create(&f->model,&f->mechanism,CHEM_LU_MOZART);
    if (rc!=CHEM_OK) return rc;
    rc=chem_state_create(&f->state,&f->model,cells);
    if (rc!=CHEM_OK) { chem_model_destroy(&f->model); return rc; }
    for (size_t c=0; c<cells; ++c) f->state.conditions[c]=(chem_conditions){298.15,101325,0};
    chem_calculate_rates(&f->state);
    return CHEM_OK;
}
static void init_destroy(init_fixture *f) { chem_state_destroy(&f->state); chem_model_destroy(&f->model); }

static void simple_case(double C,int consistent)
{
    init_fixture f; REQUIRE(init_create(&f,0,1,1)==CHEM_OK);
    f.state.y[0]=1; f.state.y[1]=.5; f.state.y[2]=C;
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3);
    chem_result r=chem_ros_solve(&f.state,.001,&p,NULL,NULL);
    CHECK(r.status==CHEM_OK);
    if (consistent) CHECK(r.stats.constraint_init_iterations<=1);
    else { CHECK(r.stats.constraint_init_iterations>0); NEAR(2*f.state.y[1]-f.state.y[2],0,1e-6); }
    init_destroy(&f);
}
TEST(constraint_initialization,ConstraintInitialization,ConsistentICsUnchanged) { simple_case(1,1); }
TEST(constraint_initialization,ConstraintInitialization,MildlyInconsistentICsCorrected) { simple_case(5,0); }
TEST(constraint_initialization,ConstraintInitialization,SeverelyInconsistentICsConverge) { simple_case(1000,0); }
TEST(constraint_initialization,ConstraintInitialization,ScaledConstraintRowsSolveIdenticallyThroughSolve)
{
    const double scales[]={1,1e-12,1e12}; double reference=0,consistent=1.3-.7-.15;
    for (size_t i=0; i<3; ++i) {
        init_fixture f; REQUIRE(init_create(&f,3,scales[i],1)==CHEM_OK);
        f.state.y[0]=.7; f.state.y[1]=.15; f.state.y[2]=5;
        chem_ros_parameters p=chem_ros_parameters_default(CHEM_RODA3);
        chem_result r=chem_ros_solve(&f.state,.001,&p,NULL,NULL);
        CHECK(r.status==CHEM_OK);
        if (!i) { reference=f.state.y[2]; NEAR(reference,consistent,1e-5*consistent); }
        else NEAR(f.state.y[2],reference,1e-5*consistent);
        init_destroy(&f);
    }
}
TEST(constraint_initialization,ConstraintInitialization,PureODESystemUnaffected)
{
    init_fixture f; REQUIRE(init_create(&f,4,1,1)==CHEM_OK);
    f.state.y[0]=1; f.state.y[1]=0;
    /* Reference does not update rate constants in this particular case. */
    f.state.rates[0]=0;
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3);
    chem_result r=chem_ros_solve(&f.state,.01,&p,NULL,NULL);
    CHECK(r.status==CHEM_OK); CHECK(r.stats.constraint_init_iterations==0);
    init_destroy(&f);
}
TEST(constraint_initialization,ConstraintInitialization,MultiCellSystems)
{
    init_fixture f; REQUIRE(init_create(&f,0,1,3)==CHEM_OK);
    const double values[]={1,2,.5, .5,.3,.1, 1,5,100}; memcpy(f.state.y,values,sizeof(values));
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3);
    chem_result r=chem_ros_solve(&f.state,.001,&p,NULL,NULL); CHECK(r.status==CHEM_OK);
    for (size_t c=0; c<3; ++c) NEAR(2*f.state.y[3+c]-f.state.y[6+c],0,1e-5);
    init_destroy(&f);
}
TEST(constraint_initialization,ConstraintInitialization,SubsequentSolveCallsReinitialize)
{
    init_fixture f; REQUIRE(init_create(&f,0,1,1)==CHEM_OK);
    f.state.y[0]=1; f.state.y[1]=.5; f.state.y[2]=1;
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3);
    chem_result r=chem_ros_solve(&f.state,.001,&p,NULL,NULL); CHECK(r.status==CHEM_OK);
    f.state.y[2]=999; chem_calculate_rates(&f.state);
    r=chem_ros_solve(&f.state,.001,&p,NULL,NULL); CHECK(r.status==CHEM_OK);
    NEAR(2*f.state.y[1]-f.state.y[2],0,1e-6); init_destroy(&f);
}
TEST(constraint_initialization,ConstraintInitialization,SolverStateToStringNewValue)
{ CHECK(!strcmp(chem_status_string(CHEM_CONSTRAINT_INIT),"Constraint Initialization Failed")); }
TEST(constraint_initialization,ConstraintInitialization,FailedInitializationRestoresCallerState)
{
    init_fixture f; REQUIRE(init_create(&f,1,1,1)==CHEM_OK);
    f.state.y[0]=.5; f.state.y[1]=100;
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3); p.init_max_iterations=1;
    chem_stats stats={0};
    CHECK(chem_initialize_constraints(&f.state,&p,&stats)==CHEM_CONSTRAINT_INIT);
    NEAR(f.state.y[0],.5,0); NEAR(f.state.y[1],100,0); init_destroy(&f);
}
TEST(constraint_initialization,ConstraintInitialization,SuccessfulInitializationKeepsTheCorrection)
{
    init_fixture f; REQUIRE(init_create(&f,1,1,1)==CHEM_OK);
    f.state.y[0]=.5; f.state.y[1]=1e-6;
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3); chem_stats stats={0};
    CHECK(chem_initialize_constraints(&f.state,&p,&stats)==CHEM_OK); NEAR(f.state.y[1],1,1e-3);
    for (size_t backtracks=0; backtracks<=2; backtracks+=2) {
        p.init_max_backtracks=backtracks; f.state.y[1]=1e-6; memset(&stats,0,sizeof(stats));
        CHECK(chem_initialize_constraints(&f.state,&p,&stats)==CHEM_CONSTRAINT_INIT);
        NEAR(f.state.y[0],.5,0); NEAR(f.state.y[1],1e-6,0);
    }
    init_destroy(&f);
}
TEST(constraint_initialization,ConstraintInitialization,WeightedCorrectionIsInvariantToConstraintRowScaling)
{
    const double scales[]={1e-12,1,1e12};
    for (size_t i=0; i<3; ++i) {
        init_fixture f; REQUIRE(init_create(&f,2,scales[i],1)==CHEM_OK);
        f.state.rtol=1e-8; f.state.atol[0]=f.state.atol[1]=1e-12;
        f.state.y[0]=1; f.state.y[1]=4;
        chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3); chem_stats stats={0};
        CHECK(chem_initialize_constraints(&f.state,&p,&stats)==CHEM_OK); NEAR(f.state.y[1],1,1e-12);
        init_destroy(&f);
    }
}
TEST(constraint_initialization,ConstraintInitialization,WeightedCorrectionUsesStateTolerances)
{
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3); NEAR(p.init_tolerance,.1,0);
    for (int tight=0; tight<2; ++tight) {
        init_fixture f; REQUIRE(init_create(&f,2,1,1)==CHEM_OK);
        f.state.rtol=tight ? 1e-8 : 1e-6; f.state.atol[0]=f.state.atol[1]=tight ? 1e-12 : 1e-8;
        f.state.y[0]=1; f.state.y[1]=1+5e-8; double initial=f.state.y[1]; chem_stats stats={0};
        CHECK(chem_initialize_constraints(&f.state,&p,&stats)==CHEM_OK);
        if (tight) { NEAR(f.state.y[1],1,1e-12); CHECK(stats.constraint_init_iterations==2); }
        else { NEAR(f.state.y[1],initial,0); CHECK(stats.constraint_init_iterations==1); }
        init_destroy(&f);
    }
}
