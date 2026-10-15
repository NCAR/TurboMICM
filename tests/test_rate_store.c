/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "test.h"
static double lambda(const ChemConditions *c,void *context) { (void)context;return c->temperature*1e-3; }
static double density(const ChemConditions *c,void *context) { (void)context;return 2*c->air_density; }
static ChemRate surface(double diffusion,double weight,double probability)
{
    ChemRate r=chem_rate_default(CHEM_SURFACE);r.label="surf";
    r.data.surface.diffusion=diffusion;r.data.surface.molecular_weight=weight;r.data.surface.probability=probability;return r;
}
static ChemMechanism *build(size_t count,const ChemRate *rates)
{
    ChemSpecies species[]={{.name="a"},{.name="b"},{.name="c",.molecular_weight=.025,.diffusion=1e-5}};
    size_t a=0,c=2;ChemProduct product={1,1};ChemReaction reactions[6];
    for (size_t i=0;i<count;++i) reactions[i]=(ChemReaction){1,1,rates[i].type==CHEM_SURFACE?&c:&a,&product,rates[i]};
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create(3,species,count,reactions,&m)==CHEM_OK);return m;
}
static size_t rate_count(const ChemMechanism *m,ChemRateType type)
{ return chem_rate_offset(m,(ChemRateType)(type+1))-chem_rate_offset(m,type); }
static void offsets(int total_only)
{
    ChemRate rates[]={chem_rate_default(CHEM_ARRHENIUS),chem_rate_default(CHEM_ARRHENIUS),chem_rate_default(CHEM_TROE),
        chem_rate_default(CHEM_TUNNELING),chem_rate_default(CHEM_USER),surface(1e-5,.025,1)};
    rates[4].label="ud1";
    if (total_only) {
        rates[1]=rates[2];rates[2]=rates[4];rates[3]=rates[5];rates[4]=chem_rate_default(CHEM_CALLBACK);
        rates[4].label="lam";rates[4].data.callback.function=lambda;
    }
    ChemMechanism *m=build(total_only?5:6,rates);
    if (!total_only) {
        size_t counts[]={2,1,0,0,0,1,0,0,1,1,0},offsets[]={0,2,3,3,3,3,4,4,4,5,6};
        for (size_t i=0;i<CHEM_RATE_TYPE_COUNT;++i) { CHECK(rate_count(m,(ChemRateType)i)==counts[i]);CHECK(chem_rate_offset(m,(ChemRateType)i)==offsets[i]); }
    }
    CHECK(chem_rate_offset(m,CHEM_CALLBACK)+rate_count(m,CHEM_CALLBACK)==chem_reaction_count(m));chem_mechanism_free(m);
}
static void arrhenius(void)
{
    ChemRate rate=chem_rate_default(CHEM_ARRHENIUS);rate.data.arrhenius=(ChemArrhenius){2.15e-4,1.2,110,300,0};
    ChemMechanism *m=build(1,&rate);CHECK(rate_count(m,CHEM_ARRHENIUS)==1);
    const ChemArrhenius *p=&chem_rate_definition(m,0)->data.arrhenius;
    CHECK(p->A==2.15e-4);CHECK(p->B==1.2);CHECK(p->C==110);CHECK(p->D==300);chem_mechanism_free(m);
}
static void branched(void)
{
    ChemRate rate=chem_rate_default(CHEM_BRANCHED);rate.data.branched=(ChemBranched){.branch=0,.X=1,.Y=0,.a0=.5,.n=3};
    ChemMechanism *m=build(1,&rate);CHECK(rate_count(m,CHEM_BRANCHED)==1);const ChemBranched *p=&chem_rate_definition(m,0)->data.branched;
    double k0=2e-22*CHEM_AVOGADRO*1e-6*exp(3),air=2.45e19/CHEM_AVOGADRO*1e6,a=k0*air,b=.43*pow(293.0/298,-8);
    double value=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2))),z=value*(1-.5)/.5;
    REL(p->k0,k0,1e-10);REL(p->z,z,1e-10);chem_mechanism_free(m);
}
static void user(void)
{
    ChemRate rates[]={chem_rate_default(CHEM_USER),chem_rate_default(CHEM_USER)};
    rates[0].label="rate1";rates[0].data.user.scaling=2;rates[1].label="rate2";rates[1].data.user.scaling=.5;
    ChemMechanism *m=build(2,rates);CHECK(rate_count(m,CHEM_USER)==2);
    for (size_t i=0;i<2;++i) { const ChemRate *r=chem_rate_definition(m,i);CHECK(r->custom_index==i);CHECK(r->data.user.scaling==(i?.5:2)); }
    chem_mechanism_free(m);
}
static void surface_test(int after_user)
{
    ChemRate rates[]={chem_rate_default(CHEM_USER),surface(after_user?1e-5:2.3e2,.025,after_user?1:.74)};
    rates[0].label="ud";ChemMechanism *m=build(after_user?2:1,rates+(after_user?0:1));
    CHECK(rate_count(m,CHEM_SURFACE)==1);const ChemRate *r=chem_rate_definition(m,after_user?1:0);
    CHECK(r->custom_index==(size_t)after_user);
    if (after_user) { CHECK(rate_count(m,CHEM_USER)==1);CHECK(chem_rate_definition(m,0)->custom_index==0); }
    else { CHECK(r->data.surface.diffusion==2.3e2);CHECK(r->data.surface.probability==.74);
        NEAR(r->data.surface.speed_factor,8*CHEM_GAS_CONSTANT/(CHEM_PI*.025),1e-14); }
    chem_mechanism_free(m);
}
static void callback(void)
{
    ChemRate rates[]={chem_rate_default(CHEM_ARRHENIUS),chem_rate_default(CHEM_CALLBACK)};
    rates[1].label="lam";rates[1].data.callback.function=lambda;ChemMechanism *m=build(2,rates);
    CHECK(rate_count(m,CHEM_ARRHENIUS)==1);CHECK(rate_count(m,CHEM_CALLBACK)==1);CHECK(chem_rate_offset(m,CHEM_CALLBACK)==1);
    const ChemRate *r=chem_rate_definition(m,1);CHECK(r->data.callback.function!=NULL);
    ChemConditions conditions={.temperature=300};NEAR(r->data.callback.function(&conditions,r->data.callback.context),.3,1e-14);
    chem_mechanism_free(m);
}
static void multipliers(void)
{
    ChemSpecies species[]={{.name="a"},{.name="b",.parameter=density}};size_t a=0,b=1;ChemProduct pa={0,1},pb={1,1};
    ChemRate rate=chem_rate_default(CHEM_ARRHENIUS);ChemReaction reactions[]={{1,1,&a,&pb,rate},{1,1,&b,&pa,rate}};
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create(2,species,2,reactions,&m)==CHEM_OK);
    CHECK(chem_parameterized_reactant_count(m,0)==0);CHECK(chem_parameterized_reactant_count(m,1)==1);
    ChemConditions c={.air_density=5};NEAR(chem_reaction_multiplier(m,1,&c),10,1e-14);chem_mechanism_free(m);
}
static void missing_property(int diffusion)
{
    ChemRate rate=surface(diffusion?0:1e-5,diffusion?.025:0,1);ChemSpecies species[]={{.name="c"},{.name="b"}};
    size_t c=0;ChemProduct b={1,1};ChemReaction reaction={1,1,&c,&b,rate};ChemMechanism *m=NULL;
    CHECK(chem_mechanism_create(2,species,1,&reaction,&m)==CHEM_INVALID);CHECK(m==NULL);
}
void test_rate_store(void)
{
#define S(name,call) RUN("reaction_rate_store:ReactionRateConstantStore." name,call)
    S("OffsetsAreContiguousCumulativeSizes",offsets(0));S("ArrheniusParametersPreserved",arrhenius());
    S("BranchedDerivedFieldsComputed",branched());S("UserDefinedCustomParamIndex",user());
    S("SurfaceDataFieldsAndCustomParamIndex",surface_test(0));S("SurfaceCustomParamIndexAfterUserDefined",surface_test(1));
    S("LambdaEntriesRcIndex",callback());S("ParameterizedMultipliers",multipliers());
    S("TotalRateConstantsMatchesProcessCount",offsets(1));S("SurfaceMissingDiffusionCoefficientThrows",missing_property(1));
    S("SurfaceMissingMolecularWeightThrows",missing_property(0));
#undef S
}
