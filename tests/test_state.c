#include "test.h"
#include <float.h>

typedef struct { chem_species species[4]; chem_reaction reactions[10]; chem_mechanism mechanism; chem_model model; chem_state state; } state_fixture;
static const char *state_names[]={"foo","bar","baz","quz"};
static const char *custom_names[]={"O1","O2","O3","AAA","BBB"};
static chem_status state_fixture_create(state_fixture *f,size_t cells,int full_custom)
{
    *f=(state_fixture){0};
    for (size_t i=0;i<4;++i) f->species[i]=(chem_species){.name=state_names[i],.absolute_tolerance=1e-10};
    for (size_t i=0;i<10;++i) {
        f->reactions[i].rate=chem_rate_default(i<(full_custom ? 5u : 2u) ? CHEM_USER : CHEM_ARRHENIUS);
        if (i<(full_custom ? 5u : 2u)) f->reactions[i].rate.label=full_custom ? custom_names[i] : i==0 ? "quux" : "corge";
    }
    f->mechanism=(chem_mechanism){.n_species=4,.n_reactions=10,.species=f->species,.reactions=f->reactions,.preserve_order=1};
    chem_status status=chem_model_create(&f->model,&f->mechanism,CHEM_LU_MOZART);
    if (status==CHEM_OK) status=chem_state_create(&f->state,&f->model,cells);
    f->state.rtol=1e-5; return status;
}
static void state_fixture_destroy(state_fixture *f) { chem_state_destroy(&f->state); chem_model_destroy(&f->model); }
TEST(state,State,DefaultConstructor)
{
    chem_state s={0},states[10]={0}; CHECK(s.cells==0); chem_state_destroy(&s);
    for (size_t i=0;i<10;++i) { CHECK(states[i].model==NULL); chem_state_destroy(states+i); }
}
TEST(state,State,Constructor)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,3,0)==CHEM_OK); chem_state *s=&f.state;
    CHECK(s->cells==3); CHECK(s->conditions!=NULL); CHECK(s->model->n_species==4);
    CHECK(s->model->n_reactions==10); CHECK(s->model->n_custom==2);
    for (size_t i=0;i<4;++i) CHECK(chem_species_index(s->model,state_names[i])==i);
    CHECK(chem_custom_index(s->model,"quux")==0); CHECK(chem_custom_index(s->model,"corge")==1);
    for (size_t c=0;c<3;++c) {
        for (size_t i=0;i<4;++i) NEAR(s->y[i*3+c],0,0);
        for (size_t i=0;i<2;++i) NEAR(s->custom[i*3+c],0,0);
        for (size_t i=0;i<10;++i) NEAR(s->rates[i*3+c],0,0);
    }
    state_fixture_destroy(&f);
}
static void state_copy_test(int move,int assigned)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,3,0)==CHEM_OK); chem_state copy={0};
    if (assigned) REQUIRE(chem_state_create(&copy,&f.model,1)==CHEM_OK);
    /* The source tests copy zeros; also exercise ownership with nonzero values. */
    for (size_t i=0;i<12;++i) f.state.y[i]=i+.25;
    for (size_t i=0;i<6;++i) f.state.custom[i]=i+.5;
    for (size_t i=0;i<30;++i) f.state.rates[i]=i+.75;
    f.state.conditions[2]=(chem_conditions){280,101325,42}; f.state.work[3]=123; f.state.jacobian[0]=456;
    double *original_y=f.state.y;
    REQUIRE((move ? chem_state_move(&copy,&f.state) : chem_state_copy(&copy,&f.state))==CHEM_OK);
    CHECK(copy.cells==3); CHECK(copy.model==&f.model); CHECK(copy.model->n_species==4);
    CHECK(copy.model->n_reactions==10); CHECK(copy.model->n_custom==2);
    NEAR(copy.rtol,1e-5,1e-5*10*DBL_EPSILON);
    for (size_t i=0;i<4;++i) { NEAR(copy.atol[i],1e-10,0); CHECK(chem_species_index(copy.model,state_names[i])==i); }
    CHECK(chem_custom_index(copy.model,"quux")==0); CHECK(chem_custom_index(copy.model,"corge")==1);
    if (move) { CHECK(copy.y==original_y); CHECK(f.state.y==NULL); CHECK(f.state.cells==0); }
    else { CHECK(copy.y!=original_y); chem_state_destroy(&f.state); }
    for (size_t i=0;i<12;++i) NEAR(copy.y[i],i+.25,0);
    for (size_t i=0;i<6;++i) NEAR(copy.custom[i],i+.5,0);
    for (size_t i=0;i<30;++i) NEAR(copy.rates[i],i+.75,0);
    NEAR(copy.conditions[2].temperature,280,0); NEAR(copy.work[3],123,0); NEAR(copy.jacobian[0],456,0);
    chem_state_destroy(&copy); state_fixture_destroy(&f);
}
TEST(state,State,CopyConstructor) { state_copy_test(0,0); }
TEST(state,State,CopyAssignmentOperator) { state_copy_test(0,1); }
TEST(state,State,MoveConstructor) { state_copy_test(1,0); }
TEST(state,State,MoveAssignmentOperator) { state_copy_test(1,1); }
TEST(state,State,SettingSingleConcentrationWithInvalidArgumentsThowsException)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,3,0)==CHEM_OK); double values[]={1,2,3,4};
    CHECK(chem_state_assign_species(&f.state,"foo",values,1)==CHEM_SIZE);
    CHECK(chem_state_assign_species(&f.state,"foo",values,2)==CHEM_SIZE);
    CHECK(chem_state_assign_species(&f.state,"not foo",values,1)==CHEM_NOT_FOUND);
    CHECK(chem_state_assign_species(&f.state,"not foo",values,4)==CHEM_NOT_FOUND);
    state_fixture_destroy(&f);
}
TEST(state,State,SetSingleConcentration)
{
    for (size_t cells=1;cells<=3;cells+=2) {
        state_fixture f; REQUIRE(state_fixture_create(&f,cells,0)==CHEM_OK);
        double values[]={12,42,35.2}; if (cells==1) values[0]=324.2;
        REQUIRE(chem_state_assign_species(&f.state,f.species[1].name,values,cells)==CHEM_OK);
        for (size_t c=0;c<cells;++c) NEAR(f.state.y[cells+c],values[c],0);
        state_fixture_destroy(&f);
    }
}
TEST(state,State,SetConcentrationByElementSingleValue)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,1,0)==CHEM_OK);
    REQUIRE(chem_state_assign_species(&f.state,"baz",(double[]){42},1)==CHEM_OK);
    REQUIRE(chem_state_assign_species(&f.state,"quz",(double[]){12},1)==CHEM_OK);
    NEAR(f.state.y[2],42,0); NEAR(f.state.y[3],12,0); state_fixture_destroy(&f);
}
TEST(state,State,SetConcentrationByElementVector)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,3,0)==CHEM_OK); double values[]={12,42,35.2};
    REQUIRE(chem_state_assign_species(&f.state,"foo",values,3)==CHEM_OK);
    for (size_t c=0;c<3;++c) NEAR(f.state.y[c],values[c],0);
    state_fixture_destroy(&f);
}
TEST(state,State,SettingConcentrationsWithInvalidArguementsThrowsException)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,3,0)==CHEM_OK);
    chem_named_column columns[]={{"FUU",(double[]){.1},1},{"bar",(double[]){.2},1},{"baz",(double[]){.3},1},{"quz",(double[]){.4},1}};
    CHECK(chem_state_assign_species_many(&f.state,4,columns)!=CHEM_OK); state_fixture_destroy(&f);
}
TEST(state,State,SetConcentrations)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,3,0)==CHEM_OK);
    chem_named_column columns[]={{"bar",(double[]){.2,.22,.222},3},{"baz",(double[]){.3,.33,.333},3},
        {"foo",(double[]){.9,.99,.999},3},{"quz",(double[]){.4,.44,.444},3}};
    REQUIRE(chem_state_assign_species_many(&f.state,4,columns)==CHEM_OK);
    double expected[]={.9,.2,.3,.4,.99,.22,.33,.44,.999,.222,.333,.444};
    for (size_t c=0;c<3;++c) for (size_t i=0;i<4;++i) NEAR(f.state.y[i*3+c],expected[c*4+i],0);
    state_fixture_destroy(&f);
}
TEST(state,State,SetStateWithSquareBracketOperator)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,1,0)==CHEM_OK); double values[]={42,12,35.2,24.2};
    for (size_t i=0;i<4;++i) { double *column=chem_state_species_column(&f.state,state_names[i]); REQUIRE(column); *column=values[i]; }
    const chem_state *cs=&f.state;
    for (size_t i=0;i<4;++i) { const double *column=chem_state_species_column_const(cs,state_names[i]); REQUIRE(column); NEAR(*column,values[i],fabs(values[i])*10*DBL_EPSILON); }
    double *foo=chem_state_species_column(&f.state,"foo"); *foo+=8; NEAR(*foo,50,0); *foo-=20; NEAR(*foo,30,0);
    *foo*=2; NEAR(*foo,60,0); *foo/=3; NEAR(*foo,20,0); state_fixture_destroy(&f);
}
TEST(state,State,SetMultiCellStateWithSquareBrackets)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,3,0)==CHEM_OK);
    for (size_t i=0;i<2;++i) {
        double *column=chem_state_species_column(&f.state,state_names[i]); REQUIRE(column);
        for (size_t c=0;c<3;++c) column[c]=(double)(i*3+c+1);
    }
    const chem_state *cs=&f.state;
    for (size_t i=0;i<2;++i) for (size_t c=0;c<3;++c) {
        NEAR(f.state.y[i*3+c],i*3+c+1,0);
        NEAR(chem_state_species_column(&f.state,state_names[i])[c],i*3+c+1,0);
        NEAR(chem_state_species_column_const(cs,state_names[i])[c],i*3+c+1,0);
    }
    state_fixture_destroy(&f);
}
static void invalid_custom_test(int bad_label)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,bad_label ? 1 : 3,1)==CHEM_OK);
    chem_named_column columns[]={{"O3",(double[]){.3},1},{"O1",(double[]){.1},1},{"O2",(double[]){.5},1},
        {bad_label ? "CCC" : "BBB",(double[]){.7},1},{"AAA",(double[]){.5},1}};
    CHECK(chem_state_assign_custom_many(&f.state,5,columns)==(bad_label ? CHEM_NOT_FOUND : CHEM_SIZE));
    for (size_t i=0;i<5*f.state.cells;++i) NEAR(f.state.custom[i],0,0);
    state_fixture_destroy(&f);
}
TEST(state,State,SettingCustomRateParameterWithInvalidVectorSizeThrowsException) { invalid_custom_test(0); }
TEST(state,State,SettingCustomRateParameterWithInvalidLabelThrowsException) { invalid_custom_test(1); }
TEST(state,State,SetCustomRateParameter)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,1,1)==CHEM_OK);
    REQUIRE(chem_state_assign_custom(&f.state,"O2",(double[]){42.3},1)==CHEM_OK);
    NEAR(f.state.custom[1],42.3,42.3*10*DBL_EPSILON); state_fixture_destroy(&f);
}
TEST(state,State,SetCustomRateParameters)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,3,1)==CHEM_OK);
    chem_named_column columns[]={{"O3",(double[]){.3,.33,.333},3},{"O1",(double[]){.1,.11,.111},3},
        {"O2",(double[]){.5,.55,.555},3},{"BBB",(double[]){.7,.77,.777},3},{"AAA",(double[]){.5,.55,.555},3}};
    REQUIRE(chem_state_assign_custom_many(&f.state,5,columns)==CHEM_OK);
    double expected[]={.1,.5,.3,.5,.7,.11,.55,.33,.55,.77,.111,.555,.333,.555,.777};
    for (size_t c=0;c<3;++c) for (size_t i=0;i<5;++i) NEAR(f.state.custom[i*3+c],expected[c*5+i],0);
    state_fixture_destroy(&f);
}
static void raw_custom_test(size_t cells)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,cells,1)==CHEM_OK);
    double values[]={.1,.2,.3,.4,.5,.1,.2,.3,.4,.5,.1,.2,.3,.4,.5};
    REQUIRE(chem_state_import_custom(&f.state,values,cells,5)==CHEM_OK);
    for (size_t c=0;c<cells;++c) for (size_t i=0;i<5;++i) NEAR(f.state.custom[i*cells+c],values[c*5+i],fabs(values[c*5+i])*10*DBL_EPSILON);
    state_fixture_destroy(&f);
}
TEST(state,State,UnsafelySetCustomRateParameterOneCell) { raw_custom_test(1); }
TEST(state,State,UnsafelySetCustomRateParameterMultiCell) { raw_custom_test(3); }
TEST(state,State,UnsafelySetCustomRateParameterCatchesTooFewGridCells)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,2,1)==CHEM_OK);
    CHECK(chem_state_import_custom(&f.state,(double[]){.1,.2,.3,.4,.5},1,5)==CHEM_SIZE); state_fixture_destroy(&f);
}
TEST(state,State,UnsafelySetCustomRateParameterCatchesTooParameters)
{
    state_fixture f; REQUIRE(state_fixture_create(&f,2,1)==CHEM_OK);
    CHECK(chem_state_import_custom(&f.state,(double[]){.1,.2,.3},1,3)==CHEM_SIZE);
    CHECK(chem_state_import_custom(&f.state,(double[]){.1,.2,.3,.1,.2,.3},2,3)==CHEM_SIZE); state_fixture_destroy(&f);
}
