/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * C ownership, named-column and checked-setter equivalents of MICM State. */
#include "test.h"
#include <float.h>
static const char *const names[]={"foo","bar","baz","quz"};
static const char *const two_labels[]={"quux","corge"};
static const char *const five_labels[]={"O1","O2","O3","AAA","BBB"};
static ChemMechanism *state_mechanism(size_t custom)
{
    ChemSpecies species[4]={{0}};for (size_t i=0;i<4;++i) species[i].name=names[i];
    ChemReaction reactions[10]={{0}};
    for (size_t i=0;i<10;++i) {
        reactions[i].rate=chem_rate_default(i<10-custom?CHEM_ARRHENIUS:CHEM_USER);
        if (i>=10-custom) reactions[i].rate.label=(custom==2?two_labels:five_labels)[i-(10-custom)];
    }
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create(4,species,10,reactions,&m)==CHEM_OK);return m;
}
static ChemState *create(ChemMechanism *m,size_t cells)
{ ChemState *s=NULL;CHECK(chem_state_create(m,cells,&s)==CHEM_OK);return s; }
static void default_constructor(void)
{
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create(0,NULL,0,NULL,&m)==CHEM_OK);
    ChemState *state=create(m,0);chem_state_free(state);
    ChemState *states[10];for (size_t i=0;i<10;++i) states[i]=create(m,0);
    for (size_t i=0;i<10;++i) chem_state_free(states[i]);
    chem_mechanism_free(m);
}
static void dimensions(const ChemState *s)
{
    CHECK(s->cells==3);CHECK(chem_species_count(s->mechanism)==4);CHECK(chem_custom_count(s->mechanism)==2);
    CHECK(chem_reaction_count(s->mechanism)==10);CHECK(s->conditions && s->values && s->custom && s->rates);
    for (size_t i=0;i<4;++i) { size_t index=99;CHECK(chem_species_index(s->mechanism,names[i],&index)==CHEM_OK);CHECK(index==i); }
    for (size_t i=0;i<2;++i) { size_t index=99;CHECK(chem_custom_index(s->mechanism,two_labels[i],&index)==CHEM_OK);CHECK(index==i); }
}
static void constructor(void)
{ ChemMechanism *m=state_mechanism(2);ChemState *s=create(m,3);dimensions(s);chem_state_free(s);chem_mechanism_free(m); }
static void ownership(int operation)
{
    ChemMechanism *m=state_mechanism(2);ChemState *original=create(m,3),*copy=NULL;
    original->relative_tolerance=1e-5;for (size_t i=0;i<4;++i) original->absolute_tolerance[i]=1e-10;
    double expected_values[12],expected_custom[6],expected_rates[30];
    memcpy(expected_values,original->values,sizeof(expected_values));memcpy(expected_custom,original->custom,sizeof(expected_custom));
    memcpy(expected_rates,original->rates,sizeof(expected_rates));
    if (operation==0) CHECK(chem_state_copy(original,&copy)==CHEM_OK);
    else if (operation==1) { copy=create(m,0);CHECK(chem_state_assign(&copy,original)==CHEM_OK); }
    else {
        if (operation==3) copy=create(m,0);
        chem_state_move(&copy,&original);CHECK(original==NULL);
    }
    dimensions(copy);NEAR(copy->relative_tolerance,1e-5,1e-5*10*DBL_EPSILON);
    for (size_t i=0;i<4;++i) CHECK(copy->absolute_tolerance[i]==1e-10);
    for (size_t i=0;i<12;++i) CHECK(copy->values[i]==expected_values[i]);
    for (size_t i=0;i<6;++i) CHECK(copy->custom[i]==expected_custom[i]);
    for (size_t i=0;i<30;++i) CHECK(copy->rates[i]==expected_rates[i]);
    if (original) { copy->values[0]=42;CHECK(original->values[0]==0); }
    chem_state_free(original);chem_state_free(copy);chem_mechanism_free(m);
}
static void invalid_concentration(void)
{
    ChemMechanism *m=state_mechanism(2);ChemState *s=create(m,3);double values[]={1,2,3,4};
    CHECK(chem_state_set_concentration(s,"foo",values,1)==CHEM_DIMENSION);
    CHECK(chem_state_set_concentration(s,"foo",values,2)==CHEM_DIMENSION);
    CHECK(chem_state_set_concentration(s,"not foo",values,1)==CHEM_NOT_FOUND);
    CHECK(chem_state_set_concentration(s,"not foo",values,4)==CHEM_NOT_FOUND);
    chem_state_free(s);chem_mechanism_free(m);
}
static void concentration(int which)
{
    ChemMechanism *m=state_mechanism(2);ChemState *s=create(m,which==1?1:3);double values[]={12,42,35.2};
    if (which==1) {
        double baz=42,quz=12;CHECK(chem_state_set_concentration(s,"baz",&baz,1)==CHEM_OK);
        CHECK(chem_state_set_concentration(s,"quz",&quz,1)==CHEM_OK);CHECK(s->values[2]==42);CHECK(s->values[3]==12);
    } else {
        size_t index=which==0?1:0;CHECK(chem_state_set_concentration(s,names[index],values,3)==CHEM_OK);
        for (size_t c=0;c<3;++c) CHECK(s->values[index*3+c]==values[c]);
        if (which==0) {
            chem_state_free(s);s=create(m,1);double value=324.2;CHECK(chem_state_set_concentration(s,"bar",&value,1)==CHEM_OK);
            NEAR(s->values[1],324.2,324.2*10*DBL_EPSILON);
        }
    }
    chem_state_free(s);chem_mechanism_free(m);
}
static void concentrations(int invalid)
{
    ChemMechanism *m=state_mechanism(2);ChemState *s=create(m,3);
    double values[][3]={{.2,.22,.222},{.3,.33,.333},{.9,.99,.999},{.4,.44,.444}};
    ChemNamedValues fields[]={{"bar",values[0],3},{"baz",values[1],3},{"foo",values[2],3},{"quz",values[3],3}};
    if (invalid) {
        double input[]={.1,.2,.3,.4};const char *bad_names[]={"FUU","bar","baz","quz"};
        for (size_t i=0;i<4;++i) fields[i]=(ChemNamedValues){bad_names[i],input+i,1};
        CHECK(chem_state_set_concentrations(s,4,fields)!=CHEM_OK);
    } else {
        CHECK(chem_state_set_concentrations(s,4,fields)==CHEM_OK);
        double expected[]={.9,.2,.3,.4,.99,.22,.33,.44,.999,.222,.333,.444};
        for (size_t c=0;c<3;++c) for (size_t v=0;v<4;++v) CHECK(s->values[v*3+c]==expected[c*4+v]);
    }
    chem_state_free(s);chem_mechanism_free(m);
}
static void named_columns(int multi)
{
    ChemMechanism *m=state_mechanism(2);ChemState *s=create(m,multi?3:1);const ChemState *read_only=s;
    if (multi) {
        double foo[]={1,2,3},bar[]={4,5,6};
        CHECK(chem_state_set_concentration(s,"foo",foo,3)==CHEM_OK);CHECK(chem_state_set_concentration(s,"bar",bar,3)==CHEM_OK);
        for (size_t c=0;c<3;++c) {
            CHECK(s->values[c]==foo[c]);CHECK(s->values[3+c]==bar[c]);
            CHECK(chem_state_concentration_data(s,"foo")[c]==foo[c]);CHECK(chem_state_concentration_data(s,"bar")[c]==bar[c]);
            CHECK(chem_state_concentration_const_data(read_only,"foo")[c]==foo[c]);CHECK(chem_state_concentration_const_data(read_only,"bar")[c]==bar[c]);
        }
    } else {
        double values[]={42,12,35.2,24.2};for (size_t i=0;i<4;++i) *chem_state_concentration_data(s,names[i])=values[i];
        for (size_t i=0;i<4;++i) {
            CHECK(s->values[i]==values[i]);NEAR(*chem_state_concentration_const_data(read_only,names[i]),values[i],fabs(values[i])*10*DBL_EPSILON);
        }
        double *foo=chem_state_concentration_data(s,"foo");*foo+=8;CHECK(*foo==50);*foo-=20;CHECK(*foo==30);
        *foo*=2;CHECK(*foo==60);*foo/=3;CHECK(*foo==20);
    }
    chem_state_free(s);chem_mechanism_free(m);
}
static void custom(int which)
{
    ChemMechanism *m=state_mechanism(5);ChemState *s=create(m,which==1 || which==2?1:3);
    double input[][3]={{.3,.33,.333},{.1,.11,.111},{.5,.55,.555},{.7,.77,.777},{.5,.55,.555}};
    const char *labels[]={"O3","O1","O2","BBB","AAA"};ChemNamedValues fields[5];
    for (size_t i=0;i<5;++i) fields[i]=(ChemNamedValues){labels[i],input[i],which<=1?1:3};
    if (which==0) CHECK(chem_state_set_custom_parameters(s,5,fields)==CHEM_DIMENSION);
    else if (which==1) { fields[3].name="CCC";CHECK(chem_state_set_custom_parameters(s,5,fields)==CHEM_NOT_FOUND); }
    else if (which==2) { double value=42.3;CHECK(chem_state_set_custom(s,"O2",&value,1)==CHEM_OK);NEAR(s->custom[1],42.3,42.3*10*DBL_EPSILON); }
    else {
        CHECK(chem_state_set_custom_parameters(s,5,fields)==CHEM_OK);
        double expected[]={.1,.5,.3,.5,.7,.11,.55,.33,.55,.77,.111,.555,.333,.555,.777};
        for (size_t c=0;c<3;++c) for (size_t p=0;p<5;++p) CHECK(s->custom[p*3+c]==expected[c*5+p]);
    }
    chem_state_free(s);chem_mechanism_free(m);
}
static void custom_matrix(int which)
{
    ChemMechanism *m=state_mechanism(5);size_t cells=which==0?1:which==1?3:2;ChemState *s=create(m,cells);
    double matrix[]={.1,.2,.3,.4,.5,.1,.2,.3,.4,.5,.1,.2,.3,.4,.5};
    if (which>=2) CHECK(chem_state_set_custom_matrix(s,1,which==3?3:5,matrix)==CHEM_DIMENSION);
    else {
        CHECK(chem_state_set_custom_matrix(s,cells,5,matrix)==CHEM_OK);
        for (size_t c=0;c<cells;++c) for (size_t p=0;p<5;++p)
            NEAR(s->custom[p*cells+c],matrix[p],fabs(matrix[p])*10*DBL_EPSILON);
    }
    chem_state_free(s);chem_mechanism_free(m);
}
void test_state(void)
{
#define S(name,call) RUN("state:State." name,call)
    S("DefaultConstructor",default_constructor());S("Constructor",constructor());S("CopyConstructor",ownership(0));
    S("CopyAssignmentOperator",ownership(1));S("MoveConstructor",ownership(2));S("MoveAssignmentOperator",ownership(3));
    S("SettingSingleConcentrationWithInvalidArgumentsThowsException",invalid_concentration());
    S("SetSingleConcentration",concentration(0));S("SetConcentrationByElementSingleValue",concentration(1));S("SetConcentrationByElementVector",concentration(2));
    S("SettingConcentrationsWithInvalidArguementsThrowsException",concentrations(1));S("SetConcentrations",concentrations(0));
    S("SetStateWithSquareBracketOperator",named_columns(0));S("SetMultiCellStateWithSquareBrackets",named_columns(1));
    S("SettingCustomRateParameterWithInvalidVectorSizeThrowsException",custom(0));S("SettingCustomRateParameterWithInvalidLabelThrowsException",custom(1));
    S("SetCustomRateParameter",custom(2));S("SetCustomRateParameters",custom(3));
    S("UnsafelySetCustomRateParameterOneCell",custom_matrix(0));S("UnsafelySetCustomRateParameterMultiCell",custom_matrix(1));
    S("UnsafelySetCustomRateParameterCatchesTooFewGridCells",custom_matrix(2));S("UnsafelySetCustomRateParameterCatchesTooParameters",custom_matrix(3));
#undef S
}
