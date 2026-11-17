#undef NDEBUG
#include <assert.h>
#include <string.h>
#include "state.h"

static void check_copy(const ChemState *a, const ChemState *b)
{
    assert(a->cells==b->cells && a->species_count==b->species_count);
    assert(a->parameter_count==b->parameter_count && a->rate_count==b->rate_count);
    assert(a->relative_tolerance==b->relative_tolerance);
    for (size_t i=0;i<a->species_count;++i) {
        assert(!strcmp(a->species_names[i],b->species_names[i]));
        assert(a->absolute_tolerances[i]==b->absolute_tolerances[i]);
    }
    for (size_t i=0;i<a->parameter_count;++i) assert(!strcmp(a->parameter_names[i],b->parameter_names[i]));
    for (size_t i=0;i<a->cells*a->species_count;++i) assert(a->concentrations[i]==b->concentrations[i]);
    for (size_t i=0;i<a->cells*a->parameter_count;++i) assert(a->parameters[i]==b->parameters[i]);
    for (size_t i=0;i<a->cells*a->rate_count;++i) assert(a->rates[i]==b->rates[i]);
}
int main(void)
{
    const char *species[] = {"foo","bar","baz","quz"};
    const char *labels[] = {"O1","O2","O3","AAA","BBB"};
    const double atol[] = {1e-10,1e-10,1e-10,1e-10};
    ChemState s={0}, copy={0}, moved={0};
    assert(!s.cells && !s.species_count && !s.concentrations);
    assert(chem_state_clone(&copy,&s)==CHEM_STATE_OK);
    chem_state_destroy(&copy);
    ChemStateConfig c={3,4,5,10,species,labels,1e-5,atol};
    assert(chem_state_init(&s,&c)==CHEM_STATE_OK);
    assert(s.cells==3 && s.species_count==4 && s.parameter_count==5 && s.rate_count==10);
    assert(!strcmp(s.species_names[0],"foo") && !strcmp(s.species_names[3],"quz"));
    assert(!strcmp(s.parameter_names[0],"O1") && !strcmp(s.parameter_names[4],"BBB"));
    double bar[]={12,42,35.2}, foo[]={.9,.99,.999}, b[]={.2,.22,.222};
    double baz[]={.3,.33,.333}, quz[]={.4,.44,.444};
    assert(chem_state_set_concentration(&s,"foo",foo,1)==CHEM_STATE_SIZE);
    assert(chem_state_set_concentration(&s,"foo",foo,2)==CHEM_STATE_SIZE);
    assert(chem_state_set_concentration(&s,"not foo",foo,3)==CHEM_STATE_UNKNOWN_NAME);
    assert(chem_state_set_concentration(&s,"not foo",foo,4)==CHEM_STATE_UNKNOWN_NAME);
    double abc[]={1,2,3}, def[]={4,5,6};
    assert(chem_state_set_concentration(&s,"foo",abc,3)==CHEM_STATE_OK);
    assert(chem_state_set_concentration(&s,"bar",def,3)==CHEM_STATE_OK);
    for (size_t i=0;i<3;++i) {
        assert(*chem_state_concentration(&s,i,"foo")==abc[i]);
        assert(*chem_state_concentration(&s,i,"bar")==def[i]);
    }
    assert(chem_state_set_concentration(&s,"bar",bar,3)==CHEM_STATE_OK);
    for (size_t i=0;i<3;++i) assert(*chem_state_concentration(&s,i,"bar")==bar[i]);
    ChemStateColumn concentration_columns[]={{"bar",b,3},{"baz",baz,3},{"foo",foo,3},{"quz",quz,3}};
    ChemStateColumn invalid_concentrations[]={{"FUU",foo,1},{"bar",b,1},{"baz",baz,1},{"quz",quz,1}};
    assert(chem_state_set_concentrations(&s,invalid_concentrations,4)!=CHEM_STATE_OK);
    assert(chem_state_set_concentrations(&s,concentration_columns,4)==CHEM_STATE_OK);
    double ordered[]={.9,.2,.3,.4,.99,.22,.33,.44,.999,.222,.333,.444};
    for (size_t i=0;i<12;++i) assert(s.concentrations[i]==ordered[i]);
    const double param[5][3]={{.1,.11,.111},{.5,.55,.555},{.3,.33,.333},{.5,.55,.555},{.7,.77,.777}};
    assert(chem_state_set_parameter(&s,"O3",param[2],1)==CHEM_STATE_SIZE);
    assert(chem_state_set_parameter(&s,"CCC",param[0],3)==CHEM_STATE_UNKNOWN_NAME);
    ChemStateColumn parameters[]={{"O3",param[2],3},{"O1",param[0],3},{"O2",param[1],3},
                                  {"BBB",param[4],3},{"AAA",param[3],3}};
    ChemStateColumn invalid_parameters[]={{"O3",param[2],1},{"O1",param[0],1}};
    assert(chem_state_set_parameters(&s,invalid_parameters,2)==CHEM_STATE_SIZE);
    invalid_parameters[0]=(ChemStateColumn){"CCC",param[2],3};
    assert(chem_state_set_parameters(&s,invalid_parameters,2)==CHEM_STATE_UNKNOWN_NAME);
    assert(chem_state_set_parameters(&s,parameters,5)==CHEM_STATE_OK);
    for (size_t cell=0;cell<3;++cell) for (size_t col=0;col<5;++col)
        assert(s.parameters[cell*5+col]==param[col][cell]);
    const double table[]={.1,.2,.3,.4,.5,.1,.2,.3,.4,.5,.1,.2,.3,.4,.5};
    assert(chem_state_set_parameter_table(&s,table,1,5)==CHEM_STATE_SIZE);
    assert(chem_state_set_parameter_table(&s,table,3,3)==CHEM_STATE_SIZE);
    assert(chem_state_set_parameter_table(&s,table,3,5)==CHEM_STATE_OK);
    for (size_t i=0;i<15;++i) assert(s.parameters[i]==table[i]);
    s.rates[12]=4.2; s.conditions[1].temperature=275;
    assert(chem_state_clone(&copy,&s)==CHEM_STATE_OK);
    check_copy(&s,&copy);
    assert(copy.conditions[1].temperature==275);
    assert(copy.concentrations!=s.concentrations && copy.species_names[0]!=s.species_names[0]);
    chem_state_move(&moved,&copy);
    check_copy(&s,&moved);
    assert(!copy.cells && !copy.concentrations);
    assert(chem_state_clone(&copy,&s)==CHEM_STATE_OK);
    assert(chem_state_clone(&copy,&s)==CHEM_STATE_OK);
    check_copy(&s,&copy);
    chem_state_move(&copy,&moved);
    check_copy(&s,&copy);
    chem_state_destroy(&s); chem_state_destroy(&copy); chem_state_destroy(&moved);
    c.cells=1;
    assert(chem_state_init(&s,&c)==CHEM_STATE_OK);
    double v=42;
    assert(chem_state_set_concentration(&s,"baz",&v,1)==CHEM_STATE_OK);
    v=12; assert(chem_state_set_concentration(&s,"quz",&v,1)==CHEM_STATE_OK);
    assert(*chem_state_concentration(&s,0,"baz")==42 && *chem_state_concentration(&s,0,"quz")==12);
    v=324.2;
    assert(chem_state_set_concentration(&s,"bar",&v,1)==CHEM_STATE_OK);
    assert(*chem_state_concentration(&s,0,"bar")==324.2);
    *chem_state_concentration(&s,0,"foo")=42;
    *chem_state_concentration(&s,0,"bar")=12;
    *chem_state_concentration(&s,0,"baz")=35.2;
    *chem_state_concentration(&s,0,"quz")=24.2;
    assert(s.concentrations[0]==42 && s.concentrations[1]==12 && s.concentrations[2]==35.2 && s.concentrations[3]==24.2);
    double *p=chem_state_concentration(&s,0,"foo");
    *p+=8; assert(*p==50); *p-=20; assert(*p==30); *p*=2; assert(*p==60); *p/=3; assert(*p==20);
    v=42.3; assert(chem_state_set_parameter(&s,"O2",&v,1)==CHEM_STATE_OK && s.parameters[1]==42.3);
    assert(chem_state_set_parameter_table(&s,table,1,5)==CHEM_STATE_OK);
    for (size_t i=0;i<5;++i) assert(s.parameters[i]==table[i]);
    chem_state_destroy(&s);
    return 0;
}
