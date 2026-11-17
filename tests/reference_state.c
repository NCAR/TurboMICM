#include "chem.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures;
static void report(const char *name, int ok)
{
    printf("state:State.%s %s %s\n",name,ok ? "PASS" : "FAIL",ok ? "reference quantities matched" : "reference quantity differed");
    failures+=!ok;
}
static const ChemSpecies species[]={{"foo",0,1e-10,0,0},{"bar",0,1e-10,0,0},
                                    {"baz",0,1e-10,0,0},{"quz",0,1e-10,0,0}};
#define EMPTY_REACTION {0,0,0,0,0,{.kind=CHEM_ARRHENIUS,.a=1,.d=300}}
static const ChemReaction reactions[]={EMPTY_REACTION,EMPTY_REACTION,EMPTY_REACTION,EMPTY_REACTION,EMPTY_REACTION,
                                       EMPTY_REACTION,EMPTY_REACTION,EMPTY_REACTION,EMPTY_REACTION,EMPTY_REACTION};
#undef EMPTY_REACTION
static const char *const labels[]={"O1","O2","O3","AAA","BBB"};
static const char *const two_labels[]={"quux","corge"};
static const double defaults[]={0,0,0,0,0};
static const ChemMechanism five={4,10,5,species,reactions,NULL,NULL,labels,defaults,298,101325};
static const ChemMechanism two={4,10,2,species,reactions,NULL,NULL,two_labels,defaults,298,101325};

static int constructor(void)
{
    ChemState *s=chem_state_create(&two,3); if (!s) return 0;
    int ok=chem_cell_count(s)==3;
    for (size_t i=0;i<4;++i) ok &= chem_find_species(&two,species[i].name)==(int)i && chem_species_values(s,i)!=NULL;
    for (size_t i=0;i<2;++i) ok &= chem_find_parameter(&two,two_labels[i])==(int)i && chem_parameter_values(s,i)!=NULL;
    for (size_t i=0;i<10;++i) ok &= chem_rate_constants(s,i)!=NULL;
    ok &= !chem_species_values(s,4) && !chem_parameter_values(s,2) && !chem_rate_constants(s,10);
    chem_state_free(s); return ok;
}

static int copy_or_move(int move, int assignment)
{
    ChemState *original=chem_state_create(&two,3), *target=assignment ? chem_state_create(&five,1) : NULL;
    if (!original || (assignment && !target)) { chem_state_free(original); chem_state_free(target); return 0; }
    for (size_t i=0;i<4;++i) for (size_t c=0;c<3;++c) chem_species_values(original,i)[c]=10*i+c;
    for (size_t i=0;i<2;++i) for (size_t c=0;c<3;++c) chem_parameter_values(original,i)[c]=100*i+c;
    int ok=chem_update_rates(original)==CHEM_OK;
    chem_state_free(target);
    target=move ? original : chem_state_clone(original);
    if (move) original=NULL;
    if (!target) { chem_state_free(original); return 0; }
    ok &= chem_cell_count(target)==3;
    for (size_t i=0;i<4;++i) {
        ok &= chem_find_species(&two,species[i].name)==(int)i && chem_absolute_tolerances(target)[i]==1e-10;
        for (size_t c=0;c<3;++c) ok &= chem_species_values(target,i)[c]==10*i+c;
    }
    for (size_t i=0;i<2;++i) for (size_t c=0;c<3;++c) ok &= chem_parameter_values(target,i)[c]==100*i+c;
    for (size_t i=0;i<10;++i) for (size_t c=0;c<3;++c) ok &= chem_rate_constants(target,i)[c]==1;
    if (!move) {
        chem_species_values(original,0)[0]=-1;
        chem_parameter_values(original,0)[0]=-1;
        ok &= chem_species_values(target,0)[0]==0 && chem_parameter_values(target,0)[0]==0;
    }
    chem_state_free(original); chem_state_free(target); return ok;
}

static int invalid_concentrations(void)
{
    ChemState *s=chem_state_create(&two,3); if (!s) return 0;
    const double vals[]={1,2,3,4};
    ChemNamedValues item={"foo",vals,1};
    int ok=chem_set_concentrations(s,&item,1)==CHEM_INVALID;
    item.count=2; ok &= chem_set_concentrations(s,&item,1)==CHEM_INVALID;
    item.name="not foo"; item.count=1; ok &= chem_set_concentrations(s,&item,1)==CHEM_INVALID;
    item.count=4; ok &= chem_set_concentrations(s,&item,1)==CHEM_INVALID;
    chem_state_free(s); return ok;
}

static int set_one(int cells, const char *name, const double *values, size_t count)
{
    ChemState *s=chem_state_create(&two,(size_t)cells); if (!s) return 0;
    ChemNamedValues item={name,values,count};
    int index=chem_find_species(&two,name), ok=chem_set_concentrations(s,&item,1)==CHEM_OK && index>=0;
    for (int c=0;c<cells && ok;++c) ok &= chem_species_values(s,(size_t)index)[c]==values[c];
    chem_state_free(s); return ok;
}

static int named_bulk(int parameters, int invalid_name, int invalid_size)
{
    static const double concentrations[][3]={{.2,.22,.222},{.3,.33,.333},{.9,.99,.999},{.4,.44,.444}};
    static const double param_values[][3]={{.3,.33,.333},{.1,.11,.111},{.5,.55,.555},{.7,.77,.777},{.5,.55,.555}};
    static const char *const concentration_names[]={"bar","baz","foo","quz"};
    static const char *const parameter_names[]={"O3","O1","O2","BBB","AAA"};
    ChemState *s=chem_state_create(&five,3); if (!s) return 0;
    ChemNamedValues items[5]; size_t count=parameters ? 5 : 4;
    for (size_t i=0;i<count;++i) items[i]=(ChemNamedValues){parameters ? parameter_names[i] : concentration_names[i],
        parameters ? param_values[i] : concentrations[i],invalid_size ? 1 : 3};
    if (invalid_name) items[0].name="not a label";
    int status=parameters ? chem_set_parameters(s,items,count) : chem_set_concentrations(s,items,count);
    int ok=status==(invalid_name || invalid_size ? CHEM_INVALID : CHEM_OK);
    for (size_t i=0;i<count && ok;++i) {
        if (invalid_name && i==0) continue;
        int index=parameters ? chem_find_parameter(&five,items[i].name) : chem_find_species(&five,items[i].name);
        double *values=parameters ? chem_parameter_values(s,(size_t)index) : chem_species_values(s,(size_t)index);
        for (size_t c=0;c<3;++c) ok &= values[c]==(invalid_name || invalid_size ? 0 : items[i].values[c]);
    }
    chem_state_free(s); return ok;
}

static int single_proxy(void)
{
    ChemState *s=chem_state_create(&two,1); if (!s) return 0;
    const double values[]={42,12,35.2,24.2};
    for (size_t i=0;i<4;++i) chem_species_values(s,i)[0]=values[i];
    const ChemState *read=s; int ok=1;
    for (size_t i=0;i<4;++i) ok &= fabs(chem_species_values_const(read,i)[0]-values[i])<=10*DBL_EPSILON*fabs(values[i]);
    double *foo=chem_species_values(s,0); *foo+=8; ok &= *foo==50;
    *foo-=20; ok &= *foo==30; *foo*=2; ok &= *foo==60; *foo/=3; ok &= *foo==20;
    chem_state_free(s); return ok;
}

static int multicell_proxy(void)
{
    ChemState *s=chem_state_create(&two,3); if (!s) return 0;
    const double foo[]={1,2,3},bar[]={4,5,6};
    memcpy(chem_species_values(s,0),foo,sizeof foo); memcpy(chem_species_values(s,1),bar,sizeof bar);
    const ChemState *read=s; int ok=1;
    for (size_t c=0;c<3;++c) ok &= chem_species_values_const(read,0)[c]==foo[c] && chem_species_values_const(read,1)[c]==bar[c];
    chem_state_free(s); return ok;
}

static int one_parameter(void)
{
    ChemState *s=chem_state_create(&five,1); if (!s) return 0;
    const double value=42.3; ChemNamedValues item={"O2",&value,1};
    int ok=chem_set_parameters(s,&item,1)==CHEM_OK && chem_parameter_values(s,1)[0]==value;
    chem_state_free(s); return ok;
}

static int row_major(size_t cells, size_t supplied_cells, size_t supplied_parameters, int valid)
{
    ChemState *s=chem_state_create(&five,cells); if (!s) return 0;
    const double values[]={.1,.2,.3,.4,.5,.1,.2,.3,.4,.5,.1,.2,.3,.4,.5};
    int result=chem_set_parameters_row_major(s,values,supplied_cells,supplied_parameters);
    int ok=result==(valid ? CHEM_OK : CHEM_INVALID);
    if (valid) for (size_t c=0;c<cells;++c) for (size_t i=0;i<5;++i)
        ok &= chem_parameter_values(s,i)[c]==values[c*5+i];
    chem_state_free(s); return ok;
}

int main(void)
{
    ChemState *states[10]={0}; int ok=1;
    for (size_t i=0;i<10;++i) { ok &= states[i]==NULL; chem_state_free(states[i]); }
    report("DefaultConstructor",ok);
    report("Constructor",constructor());
    report("CopyConstructor",copy_or_move(0,0)); report("CopyAssignmentOperator",copy_or_move(0,1));
    report("MoveConstructor",copy_or_move(1,0)); report("MoveAssignmentOperator",copy_or_move(1,1));
    report("SettingSingleConcentrationWithInvalidArgumentsThowsException",invalid_concentrations());
    const double three[]={12,42,35.2}, one[]={324.2}, baz[]={42}, quz[]={12};
    report("SetSingleConcentration",set_one(3,"bar",three,3) && set_one(1,"bar",one,1));
    report("SetConcentrationByElementSingleValue",set_one(1,"baz",baz,1) && set_one(1,"quz",quz,1));
    report("SetConcentrationByElementVector",set_one(3,"foo",three,3));
    report("SettingConcentrationsWithInvalidArguementsThrowsException",named_bulk(0,1,0));
    report("SetConcentrations",named_bulk(0,0,0));
    report("SetStateWithSquareBracketOperator",single_proxy()); report("SetMultiCellStateWithSquareBrackets",multicell_proxy());
    report("SettingCustomRateParameterWithInvalidVectorSizeThrowsException",named_bulk(1,0,1));
    report("SettingCustomRateParameterWithInvalidLabelThrowsException",named_bulk(1,1,0));
    report("SetCustomRateParameter",one_parameter()); report("SetCustomRateParameters",named_bulk(1,0,0));
    report("UnsafelySetCustomRateParameterOneCell",row_major(1,1,5,1));
    report("UnsafelySetCustomRateParameterMultiCell",row_major(3,3,5,1));
    report("UnsafelySetCustomRateParameterCatchesTooFewGridCells",row_major(2,1,5,0));
    report("UnsafelySetCustomRateParameterCatchesTooParameters",row_major(2,1,3,0));
    return failures ? 1 : 0;
}
