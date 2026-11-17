#undef NDEBUG
#include <assert.h>
#include <string.h>
#include "species.h"

static void property(ChemSpecies *s,const char *key,ChemPropertyType type,ChemPropertyValue v)
{ assert(chem_species_set_property(s,key,type,v)==CHEM_SPECIES_OK); }
static void check_properties(const ChemSpecies *s,int parameterized)
{
    ChemPropertyValue v;
    assert(!strcmp(s->name,"thing") && s->property_count==5 && s->is_parameterized==parameterized);
    assert(chem_species_get_property(s,"name [units]",CHEM_PROPERTY_DOUBLE,&v)==CHEM_SPECIES_OK && v.real==1);
    assert(chem_species_get_property(s,"name2 [units2]",CHEM_PROPERTY_DOUBLE,&v)==CHEM_SPECIES_OK && v.real==2);
    assert(chem_species_get_property(s,"foo",CHEM_PROPERTY_STRING,&v)==CHEM_SPECIES_OK && !strcmp(v.string,"bar"));
    assert(chem_species_get_property(s,"baz",CHEM_PROPERTY_INT,&v)==CHEM_SPECIES_OK && v.integer==42);
    assert(chem_species_get_property(s,"qux",CHEM_PROPERTY_BOOL,&v)==CHEM_SPECIES_OK && v.boolean==1);
}
int main(void)
{
    ChemSpecies thing={0},copy={0},assigned={0};
    assert(chem_species_init(&thing,"thing")==CHEM_SPECIES_OK && !strcmp(thing.name,"thing") && !thing.is_parameterized);
    chem_species_set_parameterization(&thing,(ChemParameterizedSpecies){0,100,0,0});
    assert(thing.is_parameterized && chem_species_evaluate(&thing,(ChemConditions){12.5,0,0})==1250);
    thing.is_parameterized=0;
    property(&thing,"name [units]",CHEM_PROPERTY_DOUBLE,(ChemPropertyValue){.real=1});
    property(&thing,"name2 [units2]",CHEM_PROPERTY_DOUBLE,(ChemPropertyValue){.real=2});
    ChemPropertyValue v;
    assert(thing.property_count==2 && chem_species_get_property(&thing,"name [units]",CHEM_PROPERTY_DOUBLE,&v)==CHEM_SPECIES_OK && v.real==1);
    assert(chem_species_get_property(&thing,"name2 [units2]",CHEM_PROPERTY_DOUBLE,&v)==CHEM_SPECIES_OK && v.real==2);
    for (int t=CHEM_PROPERTY_DOUBLE;t<=CHEM_PROPERTY_BOOL;++t)
        assert(chem_species_get_property(&thing,"not there",(ChemPropertyType)t,&v)==CHEM_SPECIES_NOT_FOUND);
    assert(chem_species_get_property(&thing,"name [units]",(ChemPropertyType)99,&v)==CHEM_SPECIES_WRONG_TYPE);
    property(&thing,"foo",CHEM_PROPERTY_STRING,(ChemPropertyValue){.string="bar"});
    property(&thing,"baz",CHEM_PROPERTY_INT,(ChemPropertyValue){.integer=42});
    property(&thing,"qux",CHEM_PROPERTY_BOOL,(ChemPropertyValue){.boolean=1});
    for (int parameterized=0;parameterized<2;++parameterized) {
        thing.is_parameterized=0;
        if (parameterized) chem_species_set_parameterization(&thing,(ChemParameterizedSpecies){15.4,0,0,0});
        assert(chem_species_clone(&copy,&thing)==CHEM_SPECIES_OK);
        assert(chem_species_clone(&assigned,&thing)==CHEM_SPECIES_OK);
        check_properties(&copy,parameterized); check_properties(&assigned,parameterized);
        assert(copy.name!=thing.name && copy.properties!=thing.properties);
        if (parameterized) assert(chem_species_evaluate(&copy,(ChemConditions){0})==15.4);
    }
    property(&thing,"foo",CHEM_PROPERTY_STRING,(ChemPropertyValue){.string="changed"});
    assert(chem_species_get_property(&copy,"foo",CHEM_PROPERTY_STRING,&v)==CHEM_SPECIES_OK && !strcmp(v.string,"bar"));
    ChemSpecies third={0}; assert(chem_species_init(&third,"M")==CHEM_SPECIES_OK);
    chem_species_set_third_body(&third);
    assert(!strcmp(third.name,"M") && third.is_parameterized && chem_species_evaluate(&third,(ChemConditions){0,0,42.4})==42.4);

    {
        ChemSpecies co2={0}; ChemPhaseSpecies no_diff={0},with_diff={0};
        assert(chem_species_init(&co2,"CO2")==CHEM_SPECIES_OK);
        assert(chem_phase_species_init(&no_diff,&co2)==CHEM_SPECIES_OK);
        assert(!strcmp(no_diff.species.name,"CO2") && !no_diff.has_diffusion);
        assert(chem_phase_species_init(&with_diff,&co2)==CHEM_SPECIES_OK);
        chem_phase_species_set_diffusion(&with_diff,1.23e-5);
        assert(!strcmp(with_diff.species.name,"CO2") && with_diff.has_diffusion && with_diff.diffusion==1.23e-5);
        chem_phase_species_set_diffusion(&no_diff,2.5e-6);
        assert(no_diff.has_diffusion && no_diff.diffusion==2.5e-6);
        chem_phase_species_destroy(&no_diff); chem_phase_species_destroy(&with_diff); chem_species_destroy(&co2);
    }
    {
        ChemSpecies items[2]={{0}}; ChemPhaseSpecies members[2]={{0}}; ChemPhase p={0}; char n[64];
        assert(chem_species_init(&items[0],"species1")==CHEM_SPECIES_OK);
        assert(chem_species_init(&items[1],"species2")==CHEM_SPECIES_OK);
        for (size_t i=0;i<2;++i) assert(chem_phase_species_init(members+i,items+i)==CHEM_SPECIES_OK);
        assert(chem_phase_init(&p,"test_phase",2,members)==CHEM_SPECIES_OK && p.count==2 && chem_phase_state_size(&p)==2);
        assert(chem_phase_state_name(&p,0,1,n,sizeof(n))==CHEM_SPECIES_OK && !strcmp(n,"test_phase.species1"));
        assert(chem_phase_state_name(&p,1,1,n,sizeof(n))==CHEM_SPECIES_OK && !strcmp(n,"test_phase.species2"));
        chem_phase_destroy(&p);
        for (size_t i=0;i<2;++i) { chem_phase_species_destroy(members+i); chem_species_destroy(items+i); }
    }

    ChemSpecies species[3]={{0}}; const char *names[]={"foo","bar","baz"};
    ChemPhaseSpecies entries[3]={{0}};
    for (size_t i=0;i<3;++i) {
        assert(chem_species_init(species+i,names[i])==CHEM_SPECIES_OK);
        assert(chem_phase_species_init(entries+i,species+i)==CHEM_SPECIES_OK);
    }
    assert(!strcmp(entries[0].species.name,"foo") && !entries[0].has_diffusion);
    chem_phase_species_set_diffusion(entries,1.23e-5);
    assert(entries[0].has_diffusion && entries[0].diffusion==1.23e-5);
    chem_phase_species_set_diffusion(entries,2.5e-6);
    assert(entries[0].has_diffusion && entries[0].diffusion==2.5e-6);
    ChemPhase phase={0};
    assert(chem_phase_init(&phase,"test_phase",2,entries)==CHEM_SPECIES_OK);
    assert(!strcmp(phase.name,"test_phase") && phase.count==2 && chem_phase_state_size(&phase)==2);
    char name[64];
    assert(chem_phase_state_name(&phase,0,1,name,sizeof(name))==CHEM_SPECIES_OK && !strcmp(name,"test_phase.foo"));
    assert(chem_phase_state_name(&phase,1,1,name,sizeof(name))==CHEM_SPECIES_OK && !strcmp(name,"test_phase.bar"));
    ChemSystem system={0}; assert(chem_system_init(&system,&phase)==CHEM_SPECIES_OK);
    assert(system.gas.count==2 && chem_phase_state_name(&system.gas,0,0,name,sizeof(name))==CHEM_SPECIES_OK && !strcmp(name,"foo"));
    chem_species_set_parameterization(&entries[1].species,(ChemParameterizedSpecies){42,0,0,0});
    assert(chem_phase_init(&phase,"gas",3,entries)==CHEM_SPECIES_OK && phase.count==3 && chem_phase_state_size(&phase)==2);
    assert(chem_phase_state_name(&phase,0,1,name,sizeof(name))==CHEM_SPECIES_OK && !strcmp(name,"gas.foo"));
    assert(chem_phase_state_name(&phase,1,1,name,sizeof(name))==CHEM_SPECIES_OK && !strcmp(name,"gas.baz"));
    assert(chem_phase_state_name(&phase,2,1,name,sizeof(name))==CHEM_SPECIES_RANGE);
    assert(chem_phase_state_name(&phase,0,1,name,3)==CHEM_SPECIES_RANGE);
    assert(chem_system_init(&system,&phase)==CHEM_SPECIES_OK && system.gas.count==3 && chem_phase_state_size(&system.gas)==2);
    const size_t reordered[]={1,0};
    assert(chem_phase_state_name(&system.gas,reordered[0],0,name,sizeof(name))==CHEM_SPECIES_OK && !strcmp(name,"baz"));
    assert(chem_phase_state_name(&system.gas,reordered[1],0,name,sizeof(name))==CHEM_SPECIES_OK && !strcmp(name,"foo"));
    entries[1].species.is_parameterized=0;
    assert(chem_species_init(species+2,"paramSpecies")==CHEM_SPECIES_OK);
    chem_species_set_parameterization(species+2,(ChemParameterizedSpecies){64.2,0,0,0});
    assert(chem_phase_species_init(entries+2,species+2)==CHEM_SPECIES_OK);
    assert(chem_phase_init(&phase,"gas",3,entries)==CHEM_SPECIES_OK);
    assert(chem_system_init(&system,&phase)==CHEM_SPECIES_OK && system.gas.count==3 && chem_phase_state_size(&system.gas)==2);
    assert(chem_phase_state_name(&system.gas,1,0,name,sizeof(name))==CHEM_SPECIES_OK && !strcmp(name,"bar"));
    assert(chem_phase_state_name(&system.gas,reordered[0],0,name,sizeof(name))==CHEM_SPECIES_OK && !strcmp(name,"bar"));
    assert(chem_phase_state_name(&system.gas,reordered[1],0,name,sizeof(name))==CHEM_SPECIES_OK && !strcmp(name,"foo"));
    chem_system_destroy(&system); chem_phase_destroy(&phase);
    for (size_t i=0;i<3;++i) { chem_phase_species_destroy(entries+i); chem_species_destroy(species+i); }
    chem_species_destroy(&third); chem_species_destroy(&assigned); chem_species_destroy(&copy); chem_species_destroy(&thing);
    return 0;
}
