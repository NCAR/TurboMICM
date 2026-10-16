/* Linear and equilibrium constraint cases from MICM, Apache-2.0. */
#include "test.h"
#include "chem_verify.h"

typedef struct {
    chem_species species[3]; chem_mechanism mechanism; chem_model model; chem_state state;
} constraint_fixture;
static chem_status constraint_create(constraint_fixture *f,chem_constraint *constraint,size_t n,size_t cells)
{
    memset(f,0,sizeof(*f));
    f->species[0].name="A"; f->species[1].name="B"; f->species[2].name="C";
    f->mechanism=(chem_mechanism){.n_species=n,.species=f->species,.n_constraints=1,.constraints=constraint,.preserve_order=1};
    chem_status rc=chem_model_create(&f->model,&f->mechanism,CHEM_LU_MOZART);
    if (rc!=CHEM_OK) return rc;
    rc=chem_state_create(&f->state,&f->model,cells);
    if (rc!=CHEM_OK) { chem_model_destroy(&f->model); return rc; }
    for (size_t c=0; c<cells; ++c) f->state.conditions[c].temperature=298.15;
    chem_update_constraint_parameters(&f->state);
    return CHEM_OK;
}
static void constraint_destroy(constraint_fixture *f) { chem_state_destroy(&f->state); chem_model_destroy(&f->model); }
static double jac(constraint_fixture *f,size_t row,size_t col,size_t cell)
{
    size_t flat=f->model.lu.map[row*f->model.n_species+col];
    return flat==CHEM_NONE ? 0 : f->state.jacobian[flat*f->state.cells+cell];
}
static chem_constraint linear(const char *name,size_t row,const chem_product *terms,size_t count,double constant)
{
    return (chem_constraint){.name=name,.type=CHEM_LINEAR_CONSTRAINT,.algebraic_species=row,
        .n_reactants=count,.reactants=terms,.constant=constant};
}
static chem_constraint equilibrium(const char *name,size_t row,const chem_product *left,size_t nl,
                                    const chem_product *right,size_t nr,double K)
{
    return (chem_constraint){.name=name,.type=CHEM_EQUILIBRIUM_CONSTRAINT,.algebraic_species=row,
        .n_reactants=nl,.n_products=nr,.reactants=left,.products=right,.equilibrium={.K=K,.delta_H=-2400}};
}
static void verify_fd(constraint_fixture *f);

static void linear_metadata(int kind)
{
    const char *names[]={"A_B_conservation","A_B_C_conservation","weighted_sum","A_equals_B","fractional_conservation"};
    const size_t count[]={2,3,3,2,2};
    const double coefficients[][3]={{1,1,0},{1,1,1},{2,3,-1},{1,-1,0},{.5,1.5,0}};
    const double constants[]={1,10,5,0,2};
    chem_product terms[3];
    for (size_t i=0; i<count[kind]; ++i) terms[i]=(chem_product){i,coefficients[kind][i]};
    chem_constraint con=linear(names[kind],count[kind]-1,terms,count[kind],constants[kind]);
    constraint_fixture f; REQUIRE(constraint_create(&f,&con,count[kind],1)==CHEM_OK);
    CHECK(!strcmp(f.mechanism.constraints[0].name,names[kind]));
    NEAR(f.mechanism.constraints[0].constant,constants[kind],0);
    CHECK(f.model.constraint_rows[0]==count[kind]-1); CHECK(f.model.n_custom==0);
    CHECK(con.n_reactants==count[kind]);
    for (size_t i=0; i<count[kind]; ++i) {
        CHECK(con.reactants[i].species==i); NEAR(con.reactants[i].yield,coefficients[kind][i],0);
        CHECK(f.model.lu.map[(count[kind]-1)*count[kind]+i]!=CHEM_NONE);
    }
    constraint_destroy(&f);
}
TEST(linear,LinearConstraint,Construction) { linear_metadata(0); }
TEST(linear,LinearConstraint,AlgebraicSpecies) { linear_metadata(1); }
TEST(linear,LinearConstraint,WeightedTerms) { linear_metadata(2); }
TEST(linear,LinearConstraint,ZeroConstant) { linear_metadata(3); }
TEST(linear,LinearConstraint,FractionalCoefficients) { linear_metadata(4); }

static void linear_numeric(int kind,int finite_difference)
{
    const double coeff[][3]={{1,1,0},{2,3,-1},{1,1,1},{1,-1,0},{.5,1.5,0},{2,3,0}};
    const double constants[]={1,5,10,0,2,1};
    const size_t count[]={2,3,3,2,2,2};
    const double initial[][3]={{.3,.7,0},{1,2,3},{2,3,5},{.5,.5,0},{1,1,0},{.1,.2,0}};
    chem_product terms[3]; size_t n=count[kind],cells=finite_difference && kind==0 ? 2 : 1;
    for (size_t i=0; i<n; ++i) terms[i]=(chem_product){i,coeff[kind][i]};
    chem_constraint con=linear("test",n-1,terms,n,constants[kind]);
    constraint_fixture f; REQUIRE(constraint_create(&f,&con,n,cells)==CHEM_OK);
    for (size_t i=0; i<n; ++i) for (size_t c=0; c<cells; ++c) f.state.y[i*cells+c]=initial[kind][i];
    if (finite_difference) {
        if (kind==0) { f.state.y[1]=.8; f.state.y[3]=.2; }
        else { f.state.y[0]=1.5; f.state.y[1]=.8; f.state.y[2]=2; }
        verify_fd(&f); constraint_destroy(&f); return;
    }
    double forcing[3]={0}; chem_constraint_forcing(&f.state,f.state.y,forcing);
    if (kind!=5) NEAR(forcing[n-1],0,1e-10);
    chem_constraint_jacobian(&f.state,f.state.y,f.state.jacobian);
    for (size_t i=0; i<n; ++i) NEAR(jac(&f,n-1,i,0),-coeff[kind][i],1e-10);
    if (kind==0 || kind==2 || kind==3 || kind==5) {
        double previous[3]; for (size_t i=0; i<n; ++i) previous[i]=jac(&f,n-1,i,0);
        if (kind==0) { f.state.y[0]=.5; f.state.y[1]=.6; }
        if (kind==2) f.state.y[2]=6;
        if (kind==3) f.state.y[0]=.6;
        if (kind==5) { f.state.y[0]=10; f.state.y[1]=20; }
        memset(forcing,0,sizeof(forcing)); chem_constraint_forcing(&f.state,f.state.y,forcing);
        if (kind!=5) NEAR(forcing[n-1],kind==2 ? 1 : .1,1e-10);
        memset(f.state.jacobian,0,f.model.lu.nnz*sizeof(double));
        chem_constraint_jacobian(&f.state,f.state.y,f.state.jacobian);
        for (size_t i=0; i<n; ++i) NEAR(jac(&f,n-1,i,0),previous[i],1e-10);
    }
    constraint_destroy(&f);
}
TEST(linear,LinearConstraint,ResidualComputationThroughConstraintSet) { linear_numeric(0,0); }
TEST(linear,LinearConstraint,JacobianComputationThroughConstraintSet) { linear_numeric(0,0); }
TEST(linear,LinearConstraint,WeightedSumResidualAndJacobian) { linear_numeric(1,0); }
TEST(linear,LinearConstraint,ThreeSpeciesConservationResidual) { linear_numeric(2,0); }
TEST(linear,LinearConstraint,ZeroConstantResidual) { linear_numeric(3,0); }
TEST(linear,LinearConstraint,FractionalCoefficientsResidualAndJacobian) { linear_numeric(4,0); }
TEST(linear,LinearConstraint,JacobianIndependentOfConcentrations) { linear_numeric(5,0); }
TEST(linear,LinearConstraint,FiniteDifferenceJacobianSimpleConservation) { linear_numeric(0,1); }
TEST(linear,LinearConstraint,FiniteDifferenceJacobianWeightedSum) { linear_numeric(1,1); }

static void equilibrium_metadata(int kind)
{
    chem_product left[]={{0,1},{1,1}}, right[]={{2,1},{2,1}};
    size_t nl=kind==1 ? 1 : 2,nr=1,row=2,n=3;
    const char *name="A_B_equilibrium";
    if (kind==1) { name="A_B_simple"; right[0].species=1; row=1; n=2; }
    if (kind==2) { name="dissociation"; nl=1; nr=2; left[0].yield=2; right[0].species=1; row=1; }
    chem_constraint con=equilibrium(name,row,left,nl,right,nr,kind==1 ? 10 : kind==2 ? 100 : 1000);
    constraint_fixture f; REQUIRE(constraint_create(&f,&con,n,1)==CHEM_OK);
    CHECK(!strcmp(f.mechanism.constraints[0].name,name)); CHECK(con.n_reactants==nl); CHECK(con.n_products==nr);
    CHECK(f.model.constraint_rows[0]==row); CHECK(f.model.n_custom==1);
    CHECK(chem_custom_index(&f.model,name)==0);
    size_t position=0;
    for (size_t i=0; i<nl; ++i) { CHECK(con.reactants[i].species==position++); }
    for (size_t i=0; i<nr; ++i) { CHECK(con.products[i].species==position++); }
    if (kind==2) NEAR(con.reactants[0].yield,2,0);
    constraint_destroy(&f);
}
TEST(equilibrium,EquilibriumConstraint,Construction) { equilibrium_metadata(0); }
TEST(equilibrium,EquilibriumConstraint,AlgebraicSpecies) { equilibrium_metadata(0); }
TEST(equilibrium,EquilibriumConstraint,SingleReactantSingleProduct) { equilibrium_metadata(1); }
TEST(equilibrium,EquilibriumConstraint,MultipleReactantsAndProducts) { equilibrium_metadata(2); }
TEST(equilibrium,EquilibriumConstraint,InvalidEquilibriumConstant)
{
    chem_product left[]={{0,1}},right[]={{1,1}};
    for (int K=-1; K<=0; ++K) {
        chem_constraint con=equilibrium("invalid",1,left,1,right,1,K);
        CHECK(chem_constraint_validate(&con,2)==CHEM_INVALID);
    }
}
TEST(equilibrium,EquilibriumConstraint,EmptyReactantsThrows)
{
    chem_product right[]={{1,1}}; chem_constraint con=equilibrium("invalid",1,NULL,0,right,1,1);
    CHECK(chem_constraint_validate(&con,2)==CHEM_INVALID);
}
TEST(equilibrium,EquilibriumConstraint,EmptyProductsThrows)
{
    chem_product left[]={{0,1}}; chem_constraint con=equilibrium("invalid",1,left,1,NULL,0,1);
    CHECK(chem_constraint_validate(&con,2)==CHEM_INVALID);
}
TEST(equilibrium,EquilibriumConstraint,InvalidStoichiometryThrows)
{
    chem_product left[]={{0,1}},right[]={{1,1}};
    chem_constraint con=equilibrium("invalid",1,left,1,right,1,1);
    left[0].yield=0; CHECK(chem_constraint_validate(&con,2)==CHEM_INVALID);
    left[0].yield=-1; CHECK(chem_constraint_validate(&con,2)==CHEM_INVALID);
    left[0].yield=1; right[0].yield=0; CHECK(chem_constraint_validate(&con,2)==CHEM_INVALID);
    right[0].yield=-2; CHECK(chem_constraint_validate(&con,2)==CHEM_INVALID);
}
static void equilibrium_numeric(int kind)
{
    chem_product left[]={{0,1},{1,1}},right[]={{2,1},{2,1}};
    size_t nl=2,nr=1,row=2,cells=kind==3 ? 2 : 1;
    int complex=kind==2 || kind==4;
    if (complex) { nl=1; nr=2; left[0].yield=2; right[0].species=1; row=1; }
    chem_constraint con=equilibrium("eq",row,left,nl,right,nr,complex ? 100 : 1000);
    constraint_fixture f; REQUIRE(constraint_create(&f,&con,3,cells)==CHEM_OK);
    if (kind==0) {
        double force[3]={0}; f.state.y[0]=f.state.y[1]=f.state.y[2]=.001;
        chem_constraint_forcing(&f.state,f.state.y,force); NEAR(force[2],0,1e-10);
        f.state.y[0]=.02; f.state.y[1]=.03; f.state.y[2]=.05;
        chem_constraint_forcing(&f.state,f.state.y,force); NEAR(force[2],.55,1e-10);
    } else if (kind==1) {
        f.state.y[0]=.01; f.state.y[1]=.02; f.state.y[2]=.05;
        chem_constraint_jacobian(&f.state,f.state.y,f.state.jacobian);
        NEAR(jac(&f,2,0,0),-20,1e-8); NEAR(jac(&f,2,1,0),-10,1e-8); NEAR(jac(&f,2,2,0),1,1e-8);
    } else if (kind==2) {
        double force[3]={0}; f.state.y[0]=.1; f.state.y[1]=.5; f.state.y[2]=2;
        chem_constraint_forcing(&f.state,f.state.y,force); NEAR(force[1],0,1e-10);
    } else {
        if (kind==3) { const double values[]={.02,.1,.03,.2,.05,5}; memcpy(f.state.y,values,sizeof(values)); }
        else { f.state.y[0]=.15; f.state.y[1]=.8; f.state.y[2]=1.5; }
        verify_fd(&f);
    }
    constraint_destroy(&f);
}
TEST(equilibrium,EquilibriumConstraint,ResidualComputationThroughConstraintSet) { equilibrium_numeric(0); }
TEST(equilibrium,EquilibriumConstraint,JacobianComputationThroughConstraintSet) { equilibrium_numeric(1); }
TEST(equilibrium,EquilibriumConstraint,ComplexStoichiometryResidual) { equilibrium_numeric(2); }
TEST(equilibrium,EquilibriumConstraint,FiniteDifferenceJacobianSimple) { equilibrium_numeric(3); }
TEST(equilibrium,EquilibriumConstraint,FiniteDifferenceJacobianComplexStoichiometry) { equilibrium_numeric(4); }

static void constraint_force_callback(void *context,const double *y,double *f)
{ chem_constraint_forcing((const chem_state *)context,y,f); }
static void verify_fd(constraint_fixture *f)
{
    size_t n=f->model.n_species,cells=f->state.cells;
    double *fd=calloc(n*n*cells,sizeof(double)); REQUIRE(fd!=NULL);
    chem_constraint_jacobian(&f->state,f->state.y,f->state.jacobian);
    CHECK(chem_finite_difference(n,cells,f->state.y,constraint_force_callback,&f->state,1e-8,fd)==CHEM_OK);
    chem_jacobian_comparison result=chem_compare_jacobian(n,cells,f->model.lu.map,f->state.jacobian,fd,1e-7,1e-7);
    CHECK(result.passed);
    result=chem_check_sparsity(n,cells,f->model.lu.map,fd,1e-6); CHECK(result.passed);
    free(fd);
}
