#include "rate_store.h"
#include <math.h>
#include <stdio.h>

static double callback(ChemConditions c,void *context)
{(void)context;return c.temperature*1e-3;}
static int run(int kind)
{
    const char *names[]={"a","b","c"},*labels[]={"rate1","rate2","radius","count","lambda"};
    ChemTerm reactants[]={{0,1},{2,1}};ChemReaction r[6]={0};size_t count=0,parameters=5;
    switch(kind) {
    case 0:
        r[0].type=r[1].type=CHEM_ARRHENIUS;r[2].type=CHEM_TROE;r[3].type=CHEM_TUNNELING;r[4].type=CHEM_USER;r[5].type=CHEM_SURFACE;count=6;break;
    case 1:r[0].type=CHEM_ARRHENIUS;r[0].p[0]=2.15e-4;r[0].p[1]=1.2;r[0].p[2]=110;r[0].p[3]=300;count=1;break;
    case 2:r[0].type=CHEM_BRANCH_ALKOXY;r[0].p[0]=1;r[0].p[2]=.5;r[0].p[3]=3;count=1;break;
    case 3:r[0].type=r[1].type=CHEM_USER;r[0].parameter=0;r[1].parameter=1;r[0].p[0]=2;r[1].p[0]=.5;count=2;break;
    case 4:r[0].type=CHEM_SURFACE;r[0].parameter=0;r[0].p[0]=230;r[0].p[1]=.025;r[0].p[2]=.74;count=1;break;
    case 5:r[0].type=CHEM_USER;r[0].parameter=0;r[1].type=CHEM_SURFACE;r[1].parameter=1;count=2;break;
    case 6:r[0].type=CHEM_ARRHENIUS;r[1].type=CHEM_CALLBACK;r[1].parameter=4;count=2;break;
    case 7:r[0].type=CHEM_ARRHENIUS;r[1].type=CHEM_TROE;r[2].type=CHEM_USER;r[3].type=CHEM_SURFACE;r[4].type=CHEM_CALLBACK;r[4].parameter=4;count=5;break;
    }
    for(size_t i=0;i<count;++i) {
        r[i].reactant_count=1;
        if(r[i].type==CHEM_SURFACE) {r[i].reactant_begin=1;if(!r[i].p[0])r[i].p[0]=1e-5;if(!r[i].p[1])r[i].p[1]=.025;if(!r[i].p[2])r[i].p[2]=1;}
    }
    ChemMechanism m={.species_count=3,.reaction_count=count,.parameter_count=parameters,.species_names=names,.parameter_names=labels,.reactions=r,.reactants=reactants};
    ChemCallback callbacks[5]={{0}};callbacks[4]=(ChemCallback){callback,NULL};ChemRateStore store={0};
    if(chem_rate_store_build(&store,&m,callbacks))return 0;
    int ok=1;
    if(kind==0) {
        size_t expected[]={0,2,3,3,3,3,4,4,4,5,6,6};
        for(size_t i=0;i<CHEM_RATE_GROUPS+1;++i)if(store.offsets[i]!=expected[i])ok=0;
        if(store.count!=6)ok=0;
    } else if(kind==1) {
        ChemReaction a=store.entries[0].reaction;
        ok=store.offsets[1]-store.offsets[0]==1 && a.p[0]==2.15e-4 && a.p[1]==1.2 && a.p[2]==110 && a.p[3]==300;
    } else if(kind==2) {
        double na=6.02214076e23,k0=2e-22*na*1e-6*exp(3),air=2.45e19/na*1e6;
        double a=k0*air,b=.43*pow(293.0/298,-8),q=a/b,z=a/(1+q)*pow(.41,1/(1+pow(log10(q),2)));
        ok=store.offsets[CHEM_RATE_BRANCHED+1]-store.offsets[CHEM_RATE_BRANCHED]==1 &&
           fabs(store.entries[0].branch_k0-k0)<=1e-10*k0 && fabs(store.entries[0].branch_z-z)<=1e-10*fabs(z);
    } else if(kind==3) ok=store.entries[0].reaction.parameter==0 && store.entries[0].reaction.p[0]==2 && store.entries[1].reaction.parameter==1 && store.entries[1].reaction.p[0]==.5;
    else if(kind==4) {
        ChemRateEntry a=store.entries[0];double expected=8*(1.380649e-23*6.02214076e23)/(3.14159265358979323846*.025);
        ok=a.reaction.p[0]==230 && fabs(a.surface_speed_factor-expected)<=1e-14 && a.reaction.p[2]==.74 && a.reaction.parameter==0;
    } else if(kind==5) ok=store.entries[0].reaction.parameter==0 && store.entries[1].reaction.parameter==1;
    else if(kind==6) {
        ChemRateEntry a=store.entries[store.offsets[CHEM_RATE_CALLBACK]];ChemConditions cond={300,0,0,0};
        ok=store.offsets[CHEM_RATE_TROE]==1 && a.reaction_index==1 && a.callback.function && fabs(a.callback.function(cond,a.callback.context)-.3)<=1e-14;
    } else ok=store.count==5 && store.offsets[CHEM_RATE_GROUPS]==count && store.offsets[CHEM_RATE_CALLBACK]==4;
    chem_rate_store_free(&store);return ok;
}
int main(void)
{
    const char *names[]={"OffsetsAreContiguousCumulativeSizes","ArrheniusParametersPreserved","BranchedDerivedFieldsComputed","UserDefinedCustomParamIndex","SurfaceDataFieldsAndCustomParamIndex","SurfaceCustomParamIndexAfterUserDefined","LambdaEntriesRcIndex","TotalRateConstantsMatchesProcessCount"};
    for(int i=0;i<8;++i)printf("reaction_rate_store:ReactionRateConstantStore.%s %d\n",names[i],run(i));
    return 0;
}
