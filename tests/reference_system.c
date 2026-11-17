#include "chem_system.h"

#include <stdio.h>
#include <string.h>

static int constructor(void)
{
    ChemSpeciesInfo species={0}; int ok=chem_species_info_init(&species,"thing")==CHEM_OK;
    ok &= !strcmp(species.name,"thing") && !species.parameterized;
    chem_species_info_set_parameterization(&species,(ChemParameterizedFunction){.temperature=100});
    ok &= species.parameterized && chem_parameterized_value(species.parameterization,12.5,0,0)==1250;
    chem_species_info_destroy(&species); return ok;
}
static int properties(int missing)
{
    ChemSpeciesInfo species={0}; ChemPropertyValue value;
    int ok=chem_species_info_init(&species,"thing")==CHEM_OK;
    ok &= chem_species_info_set(&species,"name [units]",CHEM_PROPERTY_REAL,(ChemPropertyValue){.real=1})==CHEM_OK;
    ok &= chem_species_info_set(&species,"name2 [units2]",CHEM_PROPERTY_REAL,(ChemPropertyValue){.real=2})==CHEM_OK;
    ok &= !strcmp(species.name,"thing") && species.property_count==2;
    ok &= chem_species_info_get(&species,"name [units]",CHEM_PROPERTY_REAL,&value)==CHEM_OK && value.real==1;
    ok &= chem_species_info_get(&species,"name2 [units2]",CHEM_PROPERTY_REAL,&value)==CHEM_OK && value.real==2;
    if (missing) {
        for (int type=CHEM_PROPERTY_STRING;type<=CHEM_PROPERTY_BOOL;++type)
            ok &= chem_species_info_get(&species,"not there",(ChemPropertyType)type,&value)==CHEM_PROPERTY_NOT_FOUND;
        ok &= chem_species_info_get(&species,"name [units]",(ChemPropertyType)99,&value)==CHEM_INVALID_PROPERTY_TYPE;
        ok &= chem_species_info_get(&species,"name [units]",CHEM_PROPERTY_STRING,&value)==CHEM_PROPERTY_NOT_FOUND;
        ok &= chem_species_info_has(&species,"name [units]") && !chem_species_info_has(&species,"not there");
    }
    chem_species_info_destroy(&species); return ok;
}
static int third_body(void)
{
    ChemSpeciesInfo species={0}; int ok=chem_species_info_init(&species,"M")==CHEM_OK;
    chem_species_info_set_third_body(&species);
    ok &= !strcmp(species.name,"M") && species.parameterized &&
          chem_parameterized_value(species.parameterization,0,0,42.4)==42.4;
    chem_species_info_destroy(&species); return ok;
}
static int copy_species(int assignment)
{
    int ok=1;
    for (int parameterized=0;parameterized<2;++parameterized) {
        ChemSpeciesInfo species={0},copy={0}; ChemPropertyValue value;
        ok &= chem_species_info_init(&species,"thing")==CHEM_OK;
        ok &= chem_species_info_set(&species,"name [units]",CHEM_PROPERTY_REAL,(ChemPropertyValue){.real=1})==CHEM_OK;
        ok &= chem_species_info_set(&species,"name2 [units2]",CHEM_PROPERTY_REAL,(ChemPropertyValue){.real=2})==CHEM_OK;
        ok &= chem_species_info_set(&species,"foo",CHEM_PROPERTY_STRING,(ChemPropertyValue){.string="bar"})==CHEM_OK;
        ok &= chem_species_info_set(&species,"baz",CHEM_PROPERTY_INT,(ChemPropertyValue){.integer=42})==CHEM_OK;
        ok &= chem_species_info_set(&species,"qux",CHEM_PROPERTY_BOOL,(ChemPropertyValue){.boolean=1})==CHEM_OK;
        if (parameterized) chem_species_info_set_parameterization(&species,(ChemParameterizedFunction){.constant=15.4});
        if (assignment) ok &= chem_species_info_init(&copy,"old") == CHEM_OK;
        ok &= chem_species_info_copy(&copy,&species)==CHEM_OK && copy.parameterized==parameterized;
        ok &= !strcmp(copy.name,"thing") && copy.property_count==5;
        ok &= chem_species_info_get(&copy,"name [units]",CHEM_PROPERTY_REAL,&value)==CHEM_OK && value.real==1;
        ok &= chem_species_info_get(&copy,"name2 [units2]",CHEM_PROPERTY_REAL,&value)==CHEM_OK && value.real==2;
        ok &= chem_species_info_get(&copy,"foo",CHEM_PROPERTY_STRING,&value)==CHEM_OK && !strcmp(value.string,"bar");
        ok &= chem_species_info_get(&copy,"baz",CHEM_PROPERTY_INT,&value)==CHEM_OK && value.integer==42;
        ok &= chem_species_info_get(&copy,"qux",CHEM_PROPERTY_BOOL,&value)==CHEM_OK && value.boolean==1;
        if (parameterized) ok &= chem_parameterized_value(copy.parameterization,0,0,0)==15.4;
        ok &= chem_species_info_set(&species,"foo",CHEM_PROPERTY_STRING,(ChemPropertyValue){.string="changed"})==CHEM_OK;
        chem_species_info_destroy(&species);
        ok &= chem_species_info_get(&copy,"foo",CHEM_PROPERTY_STRING,&value)==CHEM_OK && !strcmp(value.string,"bar");
        ok &= chem_species_info_copy(&copy,&copy)==CHEM_OK;
        chem_species_info_destroy(&copy);
    }
    return ok;
}
static int diffusion(int mode)
{
    ChemSpeciesInfo species={0}; int ok=chem_species_info_init(&species,"CO2")==CHEM_OK;
    ChemPhaseSpecies gas={.species=&species};
    if (mode==1) { gas.has_diffusion=1; gas.diffusion=1.23e-5; }
    ok &= !strcmp(gas.species->name,"CO2") && gas.has_diffusion==(mode==1);
    if (mode==1) ok &= gas.diffusion==1.23e-5;
    if (mode==2) { gas.has_diffusion=1; gas.diffusion=2.5e-6; ok &= gas.has_diffusion && gas.diffusion==2.5e-6; }
    chem_species_info_destroy(&species); return ok;
}
static int phase_checks(int parameterized,int system_names)
{
    ChemSpeciesInfo species[3]={{0}}; ChemPhaseSpecies entries[3];
    const char *names[]={parameterized ? "foo" : "species1",parameterized ? "bar" : "species2","baz"};
    size_t count=parameterized ? 3 : 2; int ok=1;
    for (size_t i=0;i<count;++i) {
        ok &= chem_species_info_init(&species[i],names[i])==CHEM_OK;
        entries[i]=(ChemPhaseSpecies){.species=&species[i]};
    }
    ok &= chem_species_info_set(&species[0],"absolute tolerance",CHEM_PROPERTY_REAL,(ChemPropertyValue){.real=1e-8})==CHEM_OK;
    ok &= chem_species_info_set(&species[0],"initial concentration [mol m-3]",CHEM_PROPERTY_REAL,(ChemPropertyValue){.real=2})==CHEM_OK;
    entries[0].has_diffusion=1; entries[0].diffusion=3e-5;
    if (parameterized) chem_species_info_set_parameterization(&species[1],(ChemParameterizedFunction){.constant=42});
    ChemPhase phase={parameterized ? "gas" : "test_phase",entries,count}; ChemSystem system={&phase};
    ok &= chem_phase_state_size(&phase)==2 && chem_system_state_size(&system)==2;
    ChemSpecies runtime[2];
    ok &= chem_phase_export_species(&phase,runtime,2)==CHEM_OK && runtime[0].initial==2 &&
          runtime[0].absolute_tolerance==1e-8 && runtime[0].diffusion==3e-5;
    ChemMechanism mechanism={2,0,0,runtime,NULL,NULL,NULL,NULL,NULL,298,101325};
    ChemState *state=chem_state_create(&mechanism,1);
    ok &= state && chem_species_values(state,0)[0]==2 && chem_find_species(&mechanism,names[0])==0;
    chem_state_free(state);
    char first[64],second[64],reverse[2][64];
    if (system_names) {
        ok &= chem_system_state_name(&system,0,first,sizeof(first))==CHEM_OK &&
              chem_system_state_name(&system,1,second,sizeof(second))==CHEM_OK;
        ok &= !strcmp(first,names[0]) && !strcmp(second,names[parameterized ? 2 : 1]);
        size_t reorder[]={1,0};
        for (size_t i=0;i<2;++i) ok &= chem_system_state_name(&system,reorder[i],reverse[i],sizeof(reverse[i]))==CHEM_OK;
        ok &= !strcmp(reverse[0],second) && !strcmp(reverse[1],first);
    } else {
        ok &= chem_phase_state_name(&phase,0,first,sizeof(first))==CHEM_OK &&
              chem_phase_state_name(&phase,1,second,sizeof(second))==CHEM_OK;
        ok &= !strcmp(first,parameterized ? "gas.foo" : "test_phase.species1") &&
              !strcmp(second,parameterized ? "gas.baz" : "test_phase.species2");
    }
    ok &= chem_phase_state_name(&phase,2,first,sizeof(first))==CHEM_INVALID;
    for (size_t i=0;i<count;++i) chem_species_info_destroy(&species[i]);
    return ok;
}
static int system_checks(int parameterized)
{
    ChemSpeciesInfo species[3]={{0}}; ChemPhaseSpecies entries[3];
    const char *names[]={"foo","bar","paramSpecies"}; size_t count=parameterized ? 3 : 2; int ok=1;
    for (size_t i=0;i<count;++i) {
        ok &= chem_species_info_init(&species[i],names[i])==CHEM_OK; entries[i]=(ChemPhaseSpecies){.species=&species[i]};
    }
    if (parameterized) chem_species_info_set_parameterization(&species[2],(ChemParameterizedFunction){.constant=64.2});
    ChemPhase phase={"gas",entries,count}; ChemSystem system={&phase};
    char first[16],second[16],reversed[2][16]; size_t order[]={1,0};
    ok &= phase.species_count==count && chem_system_state_size(&system)==2 &&
          chem_system_state_name(&system,0,first,sizeof(first))==CHEM_OK && !strcmp(first,"foo") &&
          chem_system_state_name(&system,1,second,sizeof(second))==CHEM_OK && !strcmp(second,"bar");
    for (size_t i=0;i<2;++i) ok &= chem_system_state_name(&system,order[i],reversed[i],sizeof(reversed[i]))==CHEM_OK;
    ok &= !strcmp(reversed[0],second) && !strcmp(reversed[1],first);
    for (size_t i=0;i<count;++i) chem_species_info_destroy(&species[i]);
    return ok;
}
static int failures;
static void report(const char *prefix,const char *name,int ok)
{
    printf("%s:%s %s %s\n",prefix,name,ok ? "PASS" : "FAIL",ok ? "reference metadata matched" : "metadata differed");
    failures+=!ok;
}
int main(void)
{
    report("species","Species.StringConstructor",constructor());
    report("species","Species.StringAndVectorConstructor",properties(0));
    report("species","Species.GetProperty",properties(1));
    report("species","Species.SetThirdBody",third_body());
    report("species","Species.CopyConstructor",copy_species(0));
    report("species","Species.CopyAssignment",copy_species(1));
    report("phase","PhaseSpecies.ConstructsWithoutDiffusionCoefficient",diffusion(0));
    report("phase","PhaseSpecies.ConstructsWithDiffusionCoefficient",diffusion(1));
    report("phase","PhaseSpecies.SetDiffusionCoefficient",diffusion(2));
    report("phase","Phase.Constructor",phase_checks(0,0));
    report("phase","Phase.ConstructorWithParameterizedSpecies",phase_checks(1,0));
    report("phase","Phase.UniqueNamesWithParameterizedSpecies",phase_checks(1,0));
    report("system","System.ConstructorWithAllParameters",system_checks(0));
    report("system","System.ConstructorWithParameterizedSpecies",system_checks(1));
    return failures ? 1 : 0;
}
