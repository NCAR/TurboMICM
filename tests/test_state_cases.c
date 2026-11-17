#include "state.h"
#include <stdio.h>
#include <string.h>

static const char *variables[]={"foo","bar","baz","quz"};
static const char *params2[]={"quux","corge"};
static const char *params5[]={"O1","O2","O3","AAA","BBB"};
static const double atol[]={1e-10,1e-10,1e-10,1e-10};
static int setup(ChemState *s,size_t cells,int many,size_t width)
{
    ChemStateParameters p={10,4,many?5u:2u,variables,many?params5:params2,1e-5,atol};
    return chem_state_init(s,&p,cells,width);
}
static int compare(const ChemState *a,const ChemState *b)
{
    if(a->cells!=b->cells || a->parameters.variable_count!=b->parameters.variable_count ||
       a->parameters.parameter_count!=b->parameters.parameter_count || a->parameters.rate_count!=b->parameters.rate_count ||
       a->parameters.relative_tolerance!=b->parameters.relative_tolerance) return 0;
    for(size_t i=0;i<a->parameters.variable_count;++i) if(strcmp(a->parameters.variable_names[i],b->parameters.variable_names[i]) || a->absolute_tolerances[i]!=b->absolute_tolerances[i]) return 0;
    for(size_t i=0;i<a->parameters.parameter_count;++i) if(strcmp(a->parameters.parameter_names[i],b->parameters.parameter_names[i])) return 0;
    for(size_t i=0;i<a->cells;++i) {
        if(a->conditions[i].temperature!=b->conditions[i].temperature || a->conditions[i].pressure!=b->conditions[i].pressure) return 0;
        for(size_t j=0;j<a->parameters.variable_count;++j) if(*chem_matrix_const_at(&a->variables,i,j)!=*chem_matrix_const_at(&b->variables,i,j)) return 0;
        for(size_t j=0;j<a->parameters.parameter_count;++j) if(*chem_matrix_const_at(&a->custom_parameters,i,j)!=*chem_matrix_const_at(&b->custom_parameters,i,j)) return 0;
        for(size_t j=0;j<a->parameters.rate_count;++j) if(*chem_matrix_const_at(&a->rate_constants,i,j)!=*chem_matrix_const_at(&b->rate_constants,i,j)) return 0;
    }
    return 1;
}
static int run(int kind,size_t width)
{
    ChemState s={0},other={0},expected={0}; int ok=1;
    if(kind==0) { ChemState list[10]={{0}}; for(int i=0;i<10;++i) chem_state_free(&list[i]); return !s.cells; }
    size_t cells=(kind==7 || kind==8 || kind==12 || kind==15 || kind==16 || kind==18)?1:(kind==20 || kind==21)?2:3;
    int many=kind>=14;
    if(setup(&s,cells,many,width)) return 0;
    if(kind==1) {
        ok=s.cells==3 && s.variables.rows==3 && s.variables.columns==4 && s.custom_parameters.rows==3 && s.custom_parameters.columns==2 && s.rate_constants.rows==3 && s.rate_constants.columns==10;
        for(size_t i=0;i<4;++i) ok &= chem_state_variable_index(&s,variables[i])==(ptrdiff_t)i;
        for(size_t i=0;i<2;++i) ok &= chem_state_parameter_index(&s,params2[i])==(ptrdiff_t)i;
    } else if(kind>=2 && kind<=5) {
        for(size_t i=0;i<cells;++i) {
            s.conditions[i].temperature=273+i; s.conditions[i].pressure=1000+i;
            for(size_t j=0;j<4;++j) *chem_matrix_at(&s.variables,i,j)=10*i+j;
            for(size_t j=0;j<2;++j) *chem_matrix_at(&s.custom_parameters,i,j)=20*i+j;
            for(size_t j=0;j<10;++j) *chem_matrix_at(&s.rate_constants,i,j)=30*i+j;
        }
        if(chem_state_copy(&expected,&s)) ok=0;
        if(kind==3 || kind==5) { if(setup(&other,1,1,width)) ok=0; }
        if(kind==2 || kind==3) { if(chem_state_copy(&other,&s) || !compare(&other,&s)) ok=0; *chem_matrix_at(&s.variables,0,0)=999; if(*chem_matrix_at(&other.variables,0,0)==999) ok=0; }
        else { chem_state_move(&other,&s); if(s.cells || !compare(&other,&expected)) ok=0; }
        if(!compare(&other,&expected)) ok=0;
    } else if(kind==6) {
        double v[]={1,2,3,4};
        ok=chem_state_set_concentration(&s,"foo",v,1) && chem_state_set_concentration(&s,"foo",v,2) && chem_state_set_concentration(&s,"not foo",v,1) && chem_state_set_concentration(&s,"not foo",v,4);
    } else if(kind==7) {
        double v=324.2; ok=!chem_state_set_concentration(&s,"bar",&v,1) && *chem_state_variable(&s,0,"bar")==v;
        ChemState multi={0}; double q[]={12,42,35.2}; if(setup(&multi,3,0,width) || chem_state_set_concentration(&multi,"bar",q,3)) ok=0;
        for(size_t i=0;i<3;++i) if(*chem_state_variable(&multi,i,"bar")!=q[i]) ok=0;
        chem_state_free(&multi);
    } else if(kind==8) {
        double a=42,b=12; ok=!chem_state_set_concentration(&s,"baz",&a,1) && !chem_state_set_concentration(&s,"quz",&b,1) && *chem_state_variable(&s,0,"baz")==a && *chem_state_variable(&s,0,"quz")==b;
    } else if(kind==9) {
        double a[]={12,42,35.2}; ok=!chem_state_set_concentration(&s,"foo",a,3);
        for(size_t i=0;i<3;++i) if(*chem_state_variable(&s,i,"foo")!=a[i]) ok=0;
    } else if(kind==10) {
        double a[]={.1},b[]={.2}; ChemNamedValues v[]={{"FUU",a,1},{"bar",b,1}}; ok=chem_state_set_concentrations(&s,v,2)!=0;
    } else if(kind==11) {
        double a[]={.2,.22,.222},b[]={.3,.33,.333},c[]={.9,.99,.999},d[]={.4,.44,.444};
        ChemNamedValues v[]={{"bar",a,3},{"baz",b,3},{"foo",c,3},{"quz",d,3}};
        double expected_values[]={.9,.2,.3,.4,.99,.22,.33,.44,.999,.222,.333,.444};
        ok=!chem_state_set_concentrations(&s,v,4);
        for(size_t i=0;i<3;++i) for(size_t j=0;j<4;++j) if(*chem_matrix_at(&s.variables,i,j)!=expected_values[4*i+j]) ok=0;
    } else if(kind==12) {
        double *foo=chem_state_variable(&s,0,"foo"),*bar=chem_state_variable(&s,0,"bar"),*baz=chem_state_variable(&s,0,"baz"),*quz=chem_state_variable(&s,0,"quz");
        *foo=42;*bar=12;*baz=35.2;*quz=24.2;
        ok=*chem_state_const_variable(&s,0,"foo")==42 && *chem_state_const_variable(&s,0,"bar")==12 && *chem_state_const_variable(&s,0,"baz")==35.2 && *chem_state_const_variable(&s,0,"quz")==24.2;
        *foo+=8; *foo-=20; *foo*=2; *foo/=3; if(*foo!=20) ok=0;
    } else if(kind==13) {
        double a[]={1,2,3},b[]={4,5,6}; ok=!chem_state_set_concentration(&s,"foo",a,3) && !chem_state_set_concentration(&s,"bar",b,3);
        for(size_t i=0;i<3;++i) if(*chem_state_const_variable(&s,i,"foo")!=a[i] || *chem_state_const_variable(&s,i,"bar")!=b[i]) ok=0;
    } else if(kind==14 || kind==15) {
        double a[]={.3},b[]={.1}; ChemNamedValues v[]={{"O3",a,1},{kind==15?"CCC":"O1",b,1}};
        ok=chem_state_set_parameters(&s,v,2)!=0;
    } else if(kind==16) {
        double a=42.3; ok=!chem_state_set_parameter(&s,"O2",&a,1) && *chem_matrix_at(&s.custom_parameters,0,1)==a;
    } else if(kind==17) {
        double a[]={.3,.33,.333},b[]={.1,.11,.111},c[]={.5,.55,.555},d[]={.7,.77,.777};
        ChemNamedValues v[]={{"O3",a,3},{"O1",b,3},{"O2",c,3},{"BBB",d,3},{"AAA",c,3}};
        double expected_values[]={.1,.5,.3,.5,.7,.11,.55,.33,.55,.77,.111,.555,.333,.555,.777};
        ok=!chem_state_set_parameters(&s,v,5);
        for(size_t i=0;i<3;++i) for(size_t j=0;j<5;++j) if(*chem_matrix_at(&s.custom_parameters,i,j)!=expected_values[5*i+j]) ok=0;
    } else {
        double a[]={.1,.2,.3,.4,.5}; const double *rows[]={a,a,a}; size_t sizes[]={5,5,5};
        if(kind==20) ok=chem_state_set_parameter_rows(&s,rows,sizes,1)!=0;
        else if(kind==21) { sizes[0]=3; ok=chem_state_set_parameter_rows(&s,rows,sizes,1)!=0 && chem_state_set_parameter_rows(&s,rows,sizes,2)!=0; }
        else { ok=!chem_state_set_parameter_rows(&s,rows,sizes,cells); for(size_t i=0;i<cells;++i) for(size_t j=0;j<5;++j) if(*chem_matrix_at(&s.custom_parameters,i,j)!=a[j]) ok=0; }
    }
    chem_state_free(&expected); chem_state_free(&other); chem_state_free(&s); return ok;
}
int main(void)
{
    const char *names[]={"DefaultConstructor","Constructor","CopyConstructor","CopyAssignmentOperator","MoveConstructor","MoveAssignmentOperator","SettingSingleConcentrationWithInvalidArgumentsThowsException","SetSingleConcentration","SetConcentrationByElementSingleValue","SetConcentrationByElementVector","SettingConcentrationsWithInvalidArguementsThrowsException","SetConcentrations","SetStateWithSquareBracketOperator","SetMultiCellStateWithSquareBrackets","SettingCustomRateParameterWithInvalidVectorSizeThrowsException","SettingCustomRateParameterWithInvalidLabelThrowsException","SetCustomRateParameter","SetCustomRateParameters","UnsafelySetCustomRateParameterOneCell","UnsafelySetCustomRateParameterMultiCell","UnsafelySetCustomRateParameterCatchesTooFewGridCells","UnsafelySetCustomRateParameterCatchesTooParameters"};
    for(int i=0;i<22;++i) printf("state:State.%s %d\n",names[i],run(i,1) && run(i,4));
    return 0;
}
