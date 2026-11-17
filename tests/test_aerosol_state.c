#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "composite.h"
#include "state.h"

static const char *const gas[]={"FO2","BAR"};
static const char *const model1[]={"STUB1.MODE1.QUUX.BAZ","STUB1.MODE1.QUUX.QUX",
    "STUB1.MODE2.QUUX.BAZ","STUB1.MODE2.QUUX.QUX","STUB1.MODE2.CORGE.FO2",
    "STUB1.MODE2.CORGE.BAZ","STUB1.MODE2.CORGE.QUX"};
static const char *const model2[]={"STUB2.MODE1.NUMBER","STUB2.MODE1.QUUX.BAZ","STUB2.MODE1.QUUX.QUX",
    "STUB2.MODE2.NUMBER","STUB2.MODE2.CORGE.FO2","STUB2.MODE2.CORGE.BAZ","STUB2.MODE2.CORGE.QUX",
    "STUB2.MODE3.NUMBER","STUB2.MODE3.QUUX.BAZ","STUB2.MODE3.QUUX.QUX",
    "STUB2.MODE3.CORGE.FO2","STUB2.MODE3.CORGE.BAZ","STUB2.MODE3.CORGE.QUX"};
static const char *const params[]={"STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT",
    "STUB2.PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT"};
typedef struct { size_t n,from[2],to[2]; } Transfer;
static void forcing(void *v,const double *y,double *f)
{
    Transfer *t=v; const double k[]={.001,.002};
    for (size_t i=0;i<2;++i) { double r=k[i]*y[t->from[i]]; f[t->from[i]]-=r; f[t->to[i]]+=r; }
}
static void jacobian(void *v,const double *y,double *j)
{
    Transfer *t=v; const double k[]={.001,.002}; (void)y;
    for (size_t i=0;i<2;++i) { j[t->from[i]*t->n+t->from[i]]-=k[i]; j[t->to[i]*t->n+t->from[i]]+=k[i]; }
}
static void run(size_t cells)
{
    const double gas_atol[]={.001,.001},m1_atol[]={.001,.001,.001,.001,.001,.001,.001};
    double m2_atol[13]; for (size_t i=0;i<13;++i) m2_atol[i]=.001;
    ChemStateConfig base={cells,2,0,0,gas,NULL,1e-6,gas_atol};
    ChemStateFragment fragments[]={{7,0,model1,NULL,m1_atol},{13,2,model2,params,m2_atol}};
    ChemState state={0};
    assert(chem_state_init_composed(&state,&base,2,fragments)==CHEM_STATE_OK);
    assert(state.species_count==2+7+13 && state.parameter_count==2);
    for (size_t i=0;i<2;++i) assert(chem_state_species_index(&state,gas[i])<state.species_count);
    for (size_t i=0;i<7;++i) assert(chem_state_species_index(&state,model1[i])<state.species_count);
    for (size_t i=0;i<13;++i) assert(chem_state_species_index(&state,model2[i])<state.species_count);
    for (size_t i=0;i<2;++i) assert(chem_state_parameter_index(&state,params[i])==i);
    for (size_t i=0;i<state.species_count;++i) assert(state.absolute_tolerances[i]==.001);
    const char *keys[]={"FO2","BAR","STUB1.MODE1.QUUX.BAZ","STUB1.MODE2.CORGE.FO2",
        "STUB2.MODE3.CORGE.QUX","STUB2.MODE3.CORGE.BAZ","STUB2.MODE1.NUMBER","STUB2.MODE2.NUMBER"};
    const double values[8][3]={{1.34,1.35,1.36},{2.53,2.54,2.55},{.75,.76,.77},{1.23,1.24,1.25},
        {.42,.43,.44},{.33,.34,.35},{1000,1001,1002},{500,501,502}};
    for (size_t i=0;i<8;++i) {
        assert(chem_state_set_concentration(&state,keys[i],values[i],cells)==CHEM_STATE_OK);
        size_t index=chem_state_species_index(&state,keys[i]);
        assert(index<state.species_count);
        for (size_t c=0;c<cells;++c) {
            assert(*chem_state_concentration(&state,c,keys[i])==values[i][c]);
            assert(state.concentrations[c*state.species_count+index]==values[i][c]);
        }
    }
    char composed[64];
    snprintf(composed,sizeof(composed),"%s.MODE%u.%s.%s","STUB2",3u,"CORGE","QUX");
    assert(!strcmp(composed,keys[4]));
    snprintf(composed,sizeof(composed),"%s.MODE%u.NUMBER","STUB2",2u);
    assert(!strcmp(composed,keys[7]));

    Transfer t={state.species_count,{chem_state_species_index(&state,"FO2"),
        chem_state_species_index(&state,"STUB1.MODE1.QUUX.BAZ")},
        {chem_state_species_index(&state,"STUB1.MODE2.CORGE.FO2"),
         chem_state_species_index(&state,"STUB1.MODE2.QUUX.BAZ")}};
    ChemOdeComponent component={&t,NULL,forcing,jacobian}; ChemCompositeOde system={t.n,1,&component};
    double f[22],j[22*22];
    chem_composite_rhs(&system,state.concentrations,f);
    chem_composite_jacobian(&system,state.concentrations,j);
    for (size_t i=0;i<2;++i) {
        double k=i ? .002 : .001, expected=k*state.concentrations[t.from[i]];
        assert(f[t.from[i]]==-expected && f[t.to[i]]==expected);
        /* Reference external Jacobians use the negative-derivative convention. */
        assert(-j[t.from[i]*t.n+t.from[i]]==k && -j[t.to[i]*t.n+t.from[i]]==-k);
    }
    const char *duplicate[]={"FO2"};
    ChemStateFragment bad={1,0,duplicate,NULL,NULL};
    double *old=state.concentrations;
    assert(chem_state_init_composed(&state,&base,1,&bad)==CHEM_STATE_DUPLICATE_NAME);
    assert(state.concentrations==old && state.species_count==22 && *chem_state_concentration(&state,0,"FO2")==1.34);
    chem_state_destroy(&state);
}
int main(void) { run(1); run(3); return 0; }
