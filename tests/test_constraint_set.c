/* ConstraintSet cases ported from MICM, Apache-2.0. */
#include "test.h"

typedef struct {
    chem_species species[5]; chem_product terms[8]; chem_constraint constraints[2];
    chem_mechanism mechanism; chem_model model; chem_state state;
} set_fixture;
static chem_constraint make_constraint(const char *name,size_t row,const chem_product *left,size_t count,const chem_product *right)
{
    return (chem_constraint){.name=name,.type=CHEM_EQUILIBRIUM_CONSTRAINT,.algebraic_species=row,
        .reactants=left,.n_reactants=count,.products=right,.n_products=1,.equilibrium={.K=.033,.delta_H=-24000}};
}
/* 0: A+B=AB; 1: B=C; 2: independent A+B=AB,C=D; 3: X=Y;
   4: A=B,C*D=A; 5: A=B,A=C; 6: empty; 7: unknown species */
static chem_status set_create(set_fixture *f,int kind,size_t cells)
{
    memset(f,0,sizeof(*f));
    const char *names[]={"A","B","AB","C","D"};
    for (size_t i=0; i<5; ++i) { f->species[i].name=names[i]; f->terms[i]=(chem_product){i,1}; }
    f->constraints[0]=make_constraint("A_B_eq",2,&f->terms[0],2,&f->terms[2]);
    f->mechanism=(chem_mechanism){.n_species=3,.species=f->species,.n_constraints=1,.constraints=f->constraints,.preserve_order=1};
    if (kind==1) f->constraints[0]=make_constraint("B_C_eq",2,&f->terms[1],1,&f->terms[2]);
    if (kind==2) {
        f->mechanism.n_species=5; f->mechanism.n_constraints=2;
        f->constraints[1]=make_constraint("C_D_eq",4,&f->terms[3],1,&f->terms[4]);
    }
    if (kind==3 || kind==4 || kind==5) f->constraints[0]=make_constraint("A_B_eq",1,&f->terms[0],1,&f->terms[1]);
    if (kind==4) {
        f->mechanism.n_species=4; f->mechanism.n_constraints=2;
        f->constraints[1]=make_constraint("CD_A_eq",0,&f->terms[2],2,&f->terms[0]);
    }
    if (kind==5) {
        f->mechanism.n_constraints=2; f->constraints[1]=make_constraint("A_C_eq",2,&f->terms[0],1,&f->terms[2]);
    }
    if (kind==6) f->mechanism.n_constraints=0;
    if (kind==7) f->constraints[0].algebraic_species=5;
    chem_status status=chem_model_create(&f->model,&f->mechanism,CHEM_LU_MOZART);
    if (status!=CHEM_OK) return status;
    status=chem_state_create(&f->state,&f->model,cells);
    if (status!=CHEM_OK) { chem_model_destroy(&f->model); return status; }
    for (size_t c=0; c<cells; ++c) f->state.conditions[c].temperature=298.15;
    chem_update_constraint_parameters(&f->state); return CHEM_OK;
}
static void set_destroy(set_fixture *f) { chem_state_destroy(&f->state); chem_model_destroy(&f->model); }
static size_t nonzeros(const set_fixture *f)
{
    size_t count=0,n=f->model.n_species;
    for (size_t row=0; row<n; ++row) if (f->model.algebraic[row])
        for (size_t col=0; col<n; ++col) count+=f->model.jacobian_pattern[row*n+col]!=0;
    return count;
}
static void structure(int kind)
{
    set_fixture f; REQUIRE(set_create(&f,kind,1)==CHEM_OK);
    size_t expected_count=kind==2 || kind==4 || kind==5 ? 2 : 1;
    CHECK(f.mechanism.n_constraints==expected_count);
    const size_t expected_nnz[]={3,2,5,2,5,4}; CHECK(nonzeros(&f)==expected_nnz[kind]);
    for (size_t k=0; k<expected_count; ++k) {
        chem_constraint *con=&f.constraints[k]; size_t row=con->algebraic_species;
        CHECK(f.model.constraint_rows[k]==row); CHECK(f.model.algebraic[row]);
        for (size_t j=0; j<con->n_reactants; ++j) CHECK(f.model.jacobian_pattern[row*f.model.n_species+con->reactants[j].species]);
        CHECK(f.model.jacobian_pattern[row*f.model.n_species+con->products[0].species]);
    }
    set_destroy(&f);
}
TEST(constraint_set,ConstraintSet,Construction) { structure(0); }
TEST(constraint_set,ConstraintSet,ReplaceStateRowsMapsToAlgebraicSpecies) { structure(1); }
TEST(constraint_set,ConstraintSet,NonZeroJacobianElements) { structure(0); }
TEST(constraint_set,ConstraintSet,MultipleConstraints) { structure(2); }
TEST(constraint_set,ConstraintSet,ThreeDStateOneConstraint) { structure(3); }

static void basic_forcing(size_t group_size)
{
    const double data[2][3]={{.2,.4,.6},{.1,.3,.7}},expected[]={-.59736,-.69901};
    for (size_t begin=0; begin<2; begin+=group_size) {
        size_t cells=2-begin<group_size ? 2-begin : group_size;
        set_fixture f; REQUIRE(set_create(&f,0,cells)==CHEM_OK);
        for (size_t c=0; c<cells; ++c) for (size_t i=0; i<3; ++i) f.state.y[i*cells+c]=data[begin+c][i];
        double force[6]={0}; chem_constraint_forcing(&f.state,f.state.y,force);
        for (size_t c=0; c<cells; ++c) NEAR(force[2*cells+c],expected[begin+c],1e-5);
        set_destroy(&f);
    }
}
static double value(const set_fixture *f,size_t row,size_t col,size_t cell)
{ return f->state.jacobian[f->model.lu.map[row*f->model.n_species+col]*f->state.cells+cell]; }
static void basic_jacobian(void)
{
    set_fixture f; REQUIRE(set_create(&f,0,1)==CHEM_OK);
    f.state.y[0]=.01; f.state.y[1]=.02; f.state.y[2]=.05;
    chem_constraint_jacobian(&f.state,f.state.y,f.state.jacobian);
    NEAR(value(&f,2,0,0),-.00066,1e-10); NEAR(value(&f,2,1,0),-.00033,1e-10); NEAR(value(&f,2,2,0),1,1e-10);
    set_destroy(&f);
}
TEST(constraint_set,ConstraintSet,AddForcingTerms) { basic_forcing(2); }
TEST(constraint_set,ConstraintSet,SubtractJacobianTerms) { basic_jacobian(); }
TEST(constraint_set,ConstraintSet,EmptyConstraintSet)
{
    set_fixture f; REQUIRE(set_create(&f,6,1)==CHEM_OK);
    CHECK(f.mechanism.n_constraints==0); CHECK(nonzeros(&f)==0);
    f.state.y[0]=.1; f.state.y[1]=.2; f.state.y[2]=.3;
    double force[]={1,1,1,1}; chem_constraint_forcing(&f.state,f.state.y,force);
    for (size_t i=0; i<4; ++i) NEAR(force[i],1,0);
    set_destroy(&f);
}
TEST(constraint_set,ConstraintSet,UnknownSpeciesThrows)
{ set_fixture f; CHECK(set_create(&f,7,1)==CHEM_NOT_FOUND); }
TEST(constraint_set,ConstraintSet,FourDStateTwoConstraints)
{
    structure(4);
    set_fixture f; REQUIRE(set_create(&f,4,3)==CHEM_OK);
    const double values[]={.33,.33,.33, .01089,.01089,.02, 10,5,5, 1,1,1};
    memcpy(f.state.y,values,sizeof(values)); double force[12]={0};
    chem_constraint_forcing(&f.state,f.state.y,force);
    NEAR(force[3],0,1e-5); NEAR(force[0],0,1e-10);
    NEAR(force[4],0,1e-5); NEAR(force[1],-.165,1e-10);
    NEAR(force[5],-.00911,1e-5); NEAR(force[2],-.165,1e-10);
    chem_constraint_jacobian(&f.state,f.state.y,f.state.jacobian);
    for (size_t c=0; c<3; ++c) {
        NEAR(value(&f,1,0,c),-.033,1e-10); NEAR(value(&f,1,1,c),1,1e-10);
        NEAR(value(&f,0,2,c),-.033*f.state.y[9+c],1e-10);
        NEAR(value(&f,0,3,c),-.033*f.state.y[6+c],1e-10); NEAR(value(&f,0,0,c),1,1e-10);
    }
    set_destroy(&f);
}
TEST(constraint_set,ConstraintSet,CoupledConstraintsSharedSpecies)
{
    structure(5);
    set_fixture f; REQUIRE(set_create(&f,5,1)==CHEM_OK);
    f.state.y[0]=.1; f.state.y[1]=f.state.y[2]=.0033; double force[3]={0};
    chem_constraint_forcing(&f.state,f.state.y,force); NEAR(force[1],0,1e-10); NEAR(force[2],0,1e-10);
    chem_constraint_jacobian(&f.state,f.state.y,f.state.jacobian);
    for (size_t row=1; row<3; ++row) { NEAR(value(&f,row,0,0),-.033,1e-10); NEAR(value(&f,row,row,0),1,1e-10); }
    set_destroy(&f);
}
TEST(constraint_set,ConstraintSet,VectorizedMatricesRespectGridCellIndexing)
{
    set_fixture f; REQUIRE(set_create(&f,0,3)==CHEM_OK);
    const double values[]={.01,.03,.001, .02,.01,.002, .05,.2,.004}; memcpy(f.state.y,values,sizeof(values));
    double force[9]={0}; chem_constraint_forcing(&f.state,f.state.y,force);
    chem_constraint_jacobian(&f.state,f.state.y,f.state.jacobian);
    for (size_t c=0; c<3; ++c) {
        NEAR(force[6+c],.033*values[c]*values[3+c]-values[6+c],1e-9);
        NEAR(value(&f,2,0,c),-.033*values[3+c],1e-12); NEAR(value(&f,2,1,c),-.033*values[c],1e-12);
        NEAR(value(&f,2,2,c),1,1e-12);
    }
    set_destroy(&f);
}
TEST(constraint_set,ConstraintSet,VectorMatrix1) { structure(0); basic_forcing(1); basic_jacobian(); }
TEST(constraint_set,ConstraintSet,VectorMatrix2) { structure(0); basic_forcing(2); basic_jacobian(); }
TEST(constraint_set,ConstraintSet,VectorMatrix3) { structure(0); basic_forcing(3); basic_jacobian(); }
