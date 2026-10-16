#include "test.h"
#include "chem_system.h"

static const chem_property initial_properties[]={
    {"name [units]",CHEM_PROPERTY_DOUBLE,{.real=1}},
    {"name2 [units2]",CHEM_PROPERTY_DOUBLE,{.real=2}}
};
static void check_properties(const chem_species_info *s,int extended)
{
    CHECK(!strcmp(s->species.name,"thing"));
    CHECK(chem_species_property_count(s,CHEM_PROPERTY_DOUBLE)==2);
    chem_property_value value;
    REQUIRE(chem_species_property_get(s,"name [units]",CHEM_PROPERTY_DOUBLE,&value)==CHEM_OK); NEAR(value.real,1,0);
    REQUIRE(chem_species_property_get(s,"name2 [units2]",CHEM_PROPERTY_DOUBLE,&value)==CHEM_OK); NEAR(value.real,2,0);
    if (extended) {
        REQUIRE(chem_species_property_get(s,"foo",CHEM_PROPERTY_STRING,&value)==CHEM_OK); CHECK(!strcmp(value.string,"bar"));
        REQUIRE(chem_species_property_get(s,"baz",CHEM_PROPERTY_INT,&value)==CHEM_OK); CHECK(value.integer==42);
        REQUIRE(chem_species_property_get(s,"qux",CHEM_PROPERTY_BOOL,&value)==CHEM_OK); CHECK(value.boolean);
    }
}
TEST(species,Species,StringConstructor)
{
    chem_species_info s={0}; REQUIRE(chem_species_info_create(&s,"thing",0,NULL)==CHEM_OK);
    CHECK(!strcmp(s.species.name,"thing")); CHECK(s.species.type==CHEM_SPECIES_VARIABLE);
    chem_species_parameterize(&s,(chem_parameterization){.temperature=100}); CHECK(s.species.type==CHEM_SPECIES_PARAMETER);
    NEAR(chem_species_evaluate(&s.species,&(chem_conditions){.temperature=12.5}),1250,0);
    chem_species_info_destroy(&s);
}
TEST(species,Species,StringAndVectorConstructor)
{
    chem_species_info s={0}; REQUIRE(chem_species_info_create(&s,"thing",2,initial_properties)==CHEM_OK);
    check_properties(&s,0); chem_species_info_destroy(&s);
}
TEST(species,Species,GetProperty)
{
    chem_species_info s={0}; REQUIRE(chem_species_info_create(&s,"thing",2,initial_properties)==CHEM_OK);
    check_properties(&s,0); chem_property_value value;
    for (int type=CHEM_PROPERTY_DOUBLE;type<=CHEM_PROPERTY_BOOL;++type)
        CHECK(chem_species_property_get(&s,"not there",(chem_property_type)type,&value)==CHEM_NOT_FOUND);
    CHECK(chem_species_property_get(&s,"name [units]",(chem_property_type)99,&value)==CHEM_INVALID);
    CHECK(chem_species_property_set(&s,"bad",(chem_property_type)99,(chem_property_value){0})==CHEM_INVALID);
    CHECK(chem_species_has_property(&s,"name [units]")); CHECK(!chem_species_has_property(&s,"not there"));
    chem_species_info_destroy(&s);
}
TEST(species,Species,SetThirdBody)
{
    chem_species_info s={0}; REQUIRE(chem_species_info_create(&s,"M",0,NULL)==CHEM_OK);
    chem_species_third_body(&s); CHECK(!strcmp(s.species.name,"M")); CHECK(s.species.type==CHEM_SPECIES_THIRD_BODY);
    NEAR(chem_species_evaluate(&s.species,&(chem_conditions){.air_density=42.4}),42.4,0); chem_species_info_destroy(&s);
}
static void check_species_copy(int replace)
{
    for (int parameterized=0;parameterized<2;++parameterized) {
        chem_species_info s={0},copy={0}; REQUIRE(chem_species_info_create(&s,"thing",2,initial_properties)==CHEM_OK);
        REQUIRE(chem_species_property_set(&s,"foo",CHEM_PROPERTY_STRING,(chem_property_value){.string="bar"})==CHEM_OK);
        REQUIRE(chem_species_property_set(&s,"baz",CHEM_PROPERTY_INT,(chem_property_value){.integer=42})==CHEM_OK);
        REQUIRE(chem_species_property_set(&s,"qux",CHEM_PROPERTY_BOOL,(chem_property_value){.boolean=1})==CHEM_OK);
        if (parameterized) chem_species_parameterize(&s,(chem_parameterization){.constant=15.4});
        if (replace) {
            REQUIRE(chem_species_info_create(&copy,"old",0,NULL)==CHEM_OK);
            REQUIRE(chem_species_info_replace(&copy,&s)==CHEM_OK);
            REQUIRE(chem_species_info_replace(&copy,&copy)==CHEM_OK);
        } else REQUIRE(chem_species_info_clone(&copy,&s)==CHEM_OK);
        /* Destroying the original also proves names and strings were deep copied. */
        chem_species_info_destroy(&s); check_properties(&copy,1);
        CHECK(copy.species.type==(parameterized ? CHEM_SPECIES_PARAMETER : CHEM_SPECIES_VARIABLE));
        if (parameterized) NEAR(chem_species_evaluate(&copy.species,&(chem_conditions){0}),15.4,0);
        chem_species_info_destroy(&copy);
    }
}
TEST(species,Species,CopyConstructor) { check_species_copy(0); }
TEST(species,Species,CopyAssignment) { check_species_copy(1); }
static void phase_species_test(int mode)
{
    chem_species_info s={0}; chem_phase_species ps={0}; double diffusion=1.23e-5;
    REQUIRE(chem_species_info_create(&s,"CO2",0,NULL)==CHEM_OK);
    REQUIRE(chem_phase_species_create(&ps,&s,mode==1 ? &diffusion : NULL,NULL)==CHEM_OK);
    chem_species_info_destroy(&s); CHECK(!strcmp(ps.info.species.name,"CO2"));
    CHECK(ps.has_diffusion==(mode==1));
    if (mode==1) NEAR(ps.diffusion,diffusion,0);
    if (mode==2) { chem_phase_species_set_diffusion(&ps,2.5e-6); CHECK(ps.has_diffusion); NEAR(ps.diffusion,2.5e-6,0); }
    chem_phase_species_destroy(&ps);
}
TEST(phase,PhaseSpecies,ConstructsWithoutDiffusionCoefficient) { phase_species_test(0); }
TEST(phase,PhaseSpecies,ConstructsWithDiffusionCoefficient) { phase_species_test(1); }
TEST(phase,PhaseSpecies,SetDiffusionCoefficient) { phase_species_test(2); }
static chem_status phase_fixture(chem_phase *p,int mode)
{
    const char *names[3]={mode==0 ? "species1" : "foo", mode==0 ? "species2" : "bar",mode==3 ? "paramSpecies" : "baz"};
    chem_phase_species ps[3]={0}; chem_status status=CHEM_OK; size_t n=mode==0 || mode==2 ? 2 : 3;
    for (size_t i=0;i<n;++i) {
        status=chem_species_info_create(&ps[i].info,names[i],0,NULL); if (status!=CHEM_OK) break;
        if ((mode==1 && i==1) || (mode==3 && i==2)) chem_species_parameterize(&ps[i].info,(chem_parameterization){.constant=mode==1 ? 42 : 64.2});
    }
    if (status==CHEM_OK) status=chem_phase_create(p,mode==0 ? "test_phase" : "gas",n,ps);
    for (size_t i=0;i<n;++i) chem_phase_species_destroy(ps+i);
    return status;
}
TEST(phase,Phase,Constructor)
{
    chem_phase p={0}; REQUIRE(phase_fixture(&p,0)==CHEM_OK);
    CHECK(!strcmp(p.name,"test_phase")); CHECK(p.n_species==2); CHECK(chem_phase_state_size(&p)==2);
    chem_names names={0}; REQUIRE(chem_phase_names(&p,1,&names)==CHEM_OK); REQUIRE(names.count==2);
    CHECK(!strcmp(names.names[0],"test_phase.species1")); CHECK(!strcmp(names.names[1],"test_phase.species2"));
    chem_names_destroy(&names); chem_phase_destroy(&p);
}
TEST(phase,Phase,ConstructorWithParameterizedSpecies)
{
    chem_phase p={0}; REQUIRE(phase_fixture(&p,1)==CHEM_OK);
    CHECK(p.n_species==3); CHECK(chem_phase_state_size(&p)==2); chem_phase_destroy(&p);
}
TEST(phase,Phase,UniqueNamesWithParameterizedSpecies)
{
    chem_phase p={0}; REQUIRE(phase_fixture(&p,1)==CHEM_OK);
    chem_names names={0}; REQUIRE(chem_phase_names(&p,1,&names)==CHEM_OK); REQUIRE(names.count==2);
    CHECK(!strcmp(names.names[0],"gas.foo")); CHECK(!strcmp(names.names[1],"gas.baz"));
    chem_names_destroy(&names); chem_phase_destroy(&p);
}
static void system_test(int parameterized)
{
    chem_phase p={0}; REQUIRE(phase_fixture(&p,parameterized ? 3 : 2)==CHEM_OK);
    chem_system s={0}; REQUIRE(chem_system_create(&s,&p)==CHEM_OK); chem_phase_destroy(&p);
    CHECK(s.gas.n_species==(parameterized ? 3u : 2u)); CHECK(chem_system_state_size(&s)==2);
    chem_names names={0}; REQUIRE(chem_system_names(&s,NULL,&names)==CHEM_OK); REQUIRE(names.count==2);
    CHECK(!strcmp(names.names[0],"foo")); CHECK(!strcmp(names.names[1],"bar")); chem_names_destroy(&names);
    if (parameterized) {
        REQUIRE(chem_system_names(&s,(size_t[]){1,0},&names)==CHEM_OK); REQUIRE(names.count==2);
        CHECK(!strcmp(names.names[0],"bar")); CHECK(!strcmp(names.names[1],"foo")); chem_names_destroy(&names);
        CHECK(chem_system_names(&s,(size_t[]){0,0},&names)==CHEM_INVALID);
    }
    chem_system_destroy(&s);
}
TEST(system,System,ConstructorWithAllParameters) { system_test(0); }
TEST(system,System,ConstructorWithParameterizedSpecies) { system_test(1); }
