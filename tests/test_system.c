/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "test.h"
#include "chem_system.h"

static ChemSpeciesInfo *species(const char *name,int properties)
{
    ChemProperty p[]={{"name [units]",CHEM_PROPERTY_REAL,{.real=1}},{"name2 [units2]",CHEM_PROPERTY_REAL,{.real=2}}};
    ChemSpeciesInfo *s=NULL;CHECK(chem_species_info_create(name,properties?2:0,p,&s)==CHEM_OK);return s;
}
static ChemPropertyValue property(const ChemSpeciesInfo *s,const char *name,ChemPropertyType type)
{ChemPropertyValue v={0};CHECK(chem_species_get_property(s,name,type,&v)==CHEM_OK);return v;}
static void string_constructor(void)
{
    ChemSpeciesInfo *s=species("thing",0);CHECK(!strcmp(chem_species_info_name(s),"thing"));CHECK(!chem_species_is_parameterized(s));
    chem_species_parameterize(s,(ChemParameterization){.temperature=100,.enabled=1});
    CHECK(chem_species_is_parameterized(s));CHECK(chem_species_parameter_value(s,&(ChemConditions){.temperature=12.5})==1250);
    chem_species_info_free(s);
}
static void vector_constructor(void)
{
    ChemSpeciesInfo *s=species("thing",1);CHECK(!strcmp(chem_species_info_name(s),"thing"));CHECK(chem_species_property_count(s,CHEM_PROPERTY_REAL)==2);
    CHECK(property(s,"name [units]",CHEM_PROPERTY_REAL).real==1);CHECK(property(s,"name2 [units2]",CHEM_PROPERTY_REAL).real==2);
    chem_species_info_free(s);
}
static void get_property(void)
{
    ChemSpeciesInfo *s=species("thing",1);ChemPropertyValue v={0};
    CHECK(property(s,"name [units]",CHEM_PROPERTY_REAL).real==1);CHECK(property(s,"name2 [units2]",CHEM_PROPERTY_REAL).real==2);
    for(int t=CHEM_PROPERTY_REAL;t<=CHEM_PROPERTY_BOOL;++t)CHECK(chem_species_get_property(s,"not there",(ChemPropertyType)t,&v)==CHEM_NOT_FOUND);
    /* An unsupported type tag replaces the forbidden extended floating type. */
    CHECK(chem_species_get_property(s,"name [units]",(ChemPropertyType)99,&v)==CHEM_INVALID);
    chem_species_info_free(s);
}
static void third_body(void)
{
    ChemSpeciesInfo *s=species("M",0);chem_species_set_third_body(s);
    CHECK(!strcmp(chem_species_info_name(s),"M"));CHECK(chem_species_is_parameterized(s));
    CHECK(chem_species_parameter_value(s,&(ChemConditions){.air_density=42.4})==42.4);chem_species_info_free(s);
}
static void copy_species(int assignment)
{
    for(int parameterized=0;parameterized<=1;++parameterized){
        ChemSpeciesInfo *s=species("thing",1),*copy=NULL;
        if(parameterized)chem_species_parameterize(s,(ChemParameterization){.constant=15.4,.enabled=1});
        CHECK(chem_species_set_property(s,(ChemProperty){"foo",CHEM_PROPERTY_STRING,{.string="bar"}})==CHEM_OK);
        CHECK(chem_species_set_property(s,(ChemProperty){"baz",CHEM_PROPERTY_INT,{.integer=42}})==CHEM_OK);
        CHECK(chem_species_set_property(s,(ChemProperty){"qux",CHEM_PROPERTY_BOOL,{.boolean=1}})==CHEM_OK);
        if(assignment){copy=species("old",0);CHECK(chem_species_info_assign(&copy,s)==CHEM_OK);}
        else CHECK(chem_species_info_copy(s,&copy)==CHEM_OK);
        CHECK(!strcmp(chem_species_info_name(copy),"thing"));CHECK(chem_species_property_count(copy,CHEM_PROPERTY_REAL)==2);
        CHECK(property(copy,"name [units]",CHEM_PROPERTY_REAL).real==1);CHECK(property(copy,"name2 [units2]",CHEM_PROPERTY_REAL).real==2);
        CHECK(!strcmp(property(copy,"foo",CHEM_PROPERTY_STRING).string,"bar"));
        CHECK(property(copy,"baz",CHEM_PROPERTY_INT).integer==42);CHECK(property(copy,"qux",CHEM_PROPERTY_BOOL).boolean);
        CHECK(chem_species_is_parameterized(copy)==parameterized);
        if(parameterized){CHECK(chem_species_parameter_value(s,&(ChemConditions){0})==15.4);CHECK(chem_species_parameter_value(copy,&(ChemConditions){0})==15.4);}
        chem_species_info_free(s);
        CHECK(!strcmp(property(copy,"foo",CHEM_PROPERTY_STRING).string,"bar"));chem_species_info_free(copy);
    }
}
static void phase_species(int mode)
{
    ChemSpeciesInfo *s=species("CO2",0);ChemPhaseSpecies p={0};CHECK(chem_phase_species_create(s,&p)==CHEM_OK);
    chem_species_info_free(s);
    CHECK(!strcmp(chem_species_info_name(p.species),"CO2"));CHECK(!p.has_diffusion);
    if(mode){double d=mode==1?1.23e-5:2.5e-6;chem_phase_species_diffusion(&p,d);CHECK(p.has_diffusion);CHECK(p.diffusion==d);}
    chem_phase_species_free(&p);
}
static ChemPhase *phase(const char *name,size_t count,const char *const *names,int parameterized)
{
    ChemPhaseSpecies entries[4]={{0}};CHECK(count<=4);
    for(size_t i=0;i<count;++i){
        ChemSpeciesInfo *s=species(names[i],0);CHECK(chem_phase_species_create(s,entries+i)==CHEM_OK);chem_species_info_free(s);
        if((int)i==parameterized)chem_species_parameterize(entries[i].species,(ChemParameterization){.constant=42,.enabled=1});
    }
    ChemPhase *p=NULL;CHECK(chem_phase_create(name,count,entries,&p)==CHEM_OK);
    for(size_t i=0;i<count;++i)chem_phase_species_free(entries+i);
    return p;
}
static void phase_constructor(void)
{
    const char *names[]={"species1","species2"};ChemPhase *p=phase("test_phase",2,names,-1);
    CHECK(!strcmp(chem_phase_name(p),"test_phase"));CHECK(chem_phase_count(p)==2);CHECK(chem_phase_state_size(p)==2);
    ChemNames unique={0};CHECK(chem_phase_names(p,1,&unique)==CHEM_OK);
    CHECK(unique.count==2);CHECK(!strcmp(unique.values[0],"test_phase.species1"));CHECK(!strcmp(unique.values[1],"test_phase.species2"));
    chem_names_free(&unique);chem_phase_free(p);
}
static void parameterized_phase(int unique_names)
{
    const char *names[]={"foo","bar","baz"};ChemPhase *p=phase("gas",3,names,1);
    if(!unique_names){CHECK(chem_phase_count(p)==3);CHECK(chem_phase_state_size(p)==2);}
    else {ChemNames unique={0};CHECK(chem_phase_names(p,1,&unique)==CHEM_OK);CHECK(unique.count==2);
        CHECK(!strcmp(unique.values[0],"gas.foo"));CHECK(!strcmp(unique.values[1],"gas.baz"));chem_names_free(&unique);}
    chem_phase_free(p);
}
static const char *reorder(const char *const *names,size_t count,size_t i,void *ctx)
{(void)count;return names[((const size_t *)ctx)[i]];}
static void system_constructor(int parameterized)
{
    const char *names[]={"foo","bar","paramSpecies"};ChemPhase *p=phase("gas",parameterized?3:2,names,parameterized?2:-1);
    /* Set the original system fixture's distinct parameter value before copying. */
    if(parameterized)chem_species_parameterize(chem_phase_species(p,2)->species,(ChemParameterization){.constant=64.2,.enabled=1});
    ChemSystem *s=NULL;CHECK(chem_system_create(p,&s)==CHEM_OK);chem_phase_free(p);
    CHECK(chem_phase_count(chem_system_gas_phase(s))==(size_t)(parameterized?3:2));CHECK(chem_system_state_size(s)==2);
    ChemNames unique={0};CHECK(chem_system_names(s,NULL,NULL,&unique)==CHEM_OK);CHECK(unique.count==2);
    CHECK((!strcmp(unique.values[0],"foo")&&!strcmp(unique.values[1],"bar"))||(!strcmp(unique.values[1],"foo")&&!strcmp(unique.values[0],"bar")));
    if(parameterized){
        size_t order[]={1,0};ChemNames sorted={0};CHECK(chem_system_names(s,reorder,order,&sorted)==CHEM_OK);CHECK(sorted.count==2);
        CHECK(!strcmp(sorted.values[0],unique.values[1]));CHECK(!strcmp(sorted.values[1],unique.values[0]));chem_names_free(&sorted);
    }
    chem_names_free(&unique);chem_system_free(s);
}
static void reaction_build(void)
{
    const char *names[]={"O3","NO","NO2","O2"};ChemPhase *p=phase("gas",4,names,-1);
    const char *reactants[]={"O3","NO"};ChemNamedProduct products[]={{"NO2",1},{"O2",1}};ChemReactionInfo *r=NULL;
    CHECK(chem_reaction_info_create(p,2,reactants,2,products,chem_rate_default(CHEM_ARRHENIUS),&r)==CHEM_OK);
    CHECK(chem_reaction_definition(r)->reactant_count==2);CHECK(!strcmp(chem_reaction_product(r,0),"NO2"));
    CHECK(!strcmp(chem_phase_name(chem_reaction_phase(r)),"gas"));
    /* Exercise the metadata-to-runtime path, including name resolution. */
    ChemSystem *s=NULL;CHECK(chem_system_create(p,&s)==CHEM_OK);ChemMechanism *m=NULL;
    const ChemReactionInfo *rs[]={r};CHECK(chem_system_mechanism(s,1,rs,&m)==CHEM_OK);
    chem_system_free(s);chem_reaction_info_free(r);chem_phase_free(p);
    CHECK(chem_species_count(m)==4);CHECK(chem_reaction_count(m)==1);
    CHECK(!strcmp(chem_phase_name(chem_system_gas_phase(chem_mechanism_system(m))),"gas"));
    chem_mechanism_free(m);
}
static void reaction_copy(void)
{
    const char *names[]={"foo","bar"};ChemPhase *p=phase("gas",2,names,-1);
    ChemSpeciesInfo *foo=chem_phase_species(p,0)->species,*bar=chem_phase_species(p,1)->species;
    CHECK(chem_species_set_property(foo,(ChemProperty){"molecular weight [kg mol-1]",CHEM_PROPERTY_REAL,{.real=.025}})==CHEM_OK);
    CHECK(chem_species_set_property(foo,(ChemProperty){"diffusion coefficient [m2 s-1]",CHEM_PROPERTY_REAL,{.real=2.3e2}})==CHEM_OK);
    chem_species_parameterize(bar,(ChemParameterization){.air_density=.82,.enabled=1});
    ChemRate rate=chem_rate_default(CHEM_ARRHENIUS);rate.data.arrhenius.A=12.2;rate.data.arrhenius.C=300;
    ChemReactionInfo *r=NULL,*copy=NULL;CHECK(chem_reaction_info_create(p,2,names,0,NULL,rate,&r)==CHEM_OK);
    CHECK(chem_reaction_info_assign(&copy,r)==CHEM_OK);
    CHECK(!strcmp(chem_reaction_reactant(copy,0),chem_reaction_reactant(r,0)));
    const ChemReaction *original=chem_reaction_definition(r),*copied=chem_reaction_definition(copy);
    CHECK(copied->product_count==original->product_count);
    CHECK(!strcmp(chem_phase_name(chem_reaction_phase(copy)),chem_phase_name(chem_reaction_phase(r))));
    CHECK(copied->rate.type==CHEM_ARRHENIUS);CHECK(original->rate.type==CHEM_ARRHENIUS);
    CHECK(copied->rate.data.arrhenius.A==original->rate.data.arrhenius.A);
    /* A built mechanism retains coefficients independently of all builder input. */
    ChemSystem *s=NULL;CHECK(chem_system_create(p,&s)==CHEM_OK);ChemMechanism *m=NULL;
    const ChemReactionInfo *rs[]={copy};CHECK(chem_system_mechanism(s,1,rs,&m)==CHEM_OK);
    chem_system_free(s);chem_reaction_info_free(r);chem_reaction_info_free(copy);chem_phase_free(p);
    ChemState *state=NULL;CHECK(chem_state_create(m,1,&state)==CHEM_OK);state->conditions[0]=(ChemConditions){.temperature=300,.air_density=10};
    chem_update_rates(state);REL(state->rates[0],12.2*exp(1)*8.2,1e-15);
    chem_state_free(state);chem_mechanism_free(m);
}
void test_system(void)
{
    RUN("species:Species.StringConstructor",string_constructor());
    RUN("species:Species.StringAndVectorConstructor",vector_constructor());
    RUN("species:Species.GetProperty",get_property());RUN("species:Species.SetThirdBody",third_body());
    RUN("species:Species.CopyConstructor",copy_species(0));RUN("species:Species.CopyAssignment",copy_species(1));
    RUN("phase:PhaseSpecies.ConstructsWithoutDiffusionCoefficient",phase_species(0));
    RUN("phase:PhaseSpecies.ConstructsWithDiffusionCoefficient",phase_species(1));
    RUN("phase:PhaseSpecies.SetDiffusionCoefficient",phase_species(2));
    RUN("phase:Phase.Constructor",phase_constructor());
    RUN("phase:Phase.ConstructorWithParameterizedSpecies",parameterized_phase(0));
    RUN("phase:Phase.UniqueNamesWithParameterizedSpecies",parameterized_phase(1));
    RUN("system:System.ConstructorWithAllParameters",system_constructor(0));
    RUN("system:System.ConstructorWithParameterizedSpecies",system_constructor(1));
    RUN("process:Process.BuildsChemicalReaction",reaction_build());
    RUN("process:Process.ChemicalReactionCopyAssignmentSucceeds",reaction_copy());
}
