#include "chem.h"

#include <math.h>
#include <stdio.h>

static int failures;
static void report(const char *name, int ok)
{
    printf("reaction_rate_store:ReactionRateConstantStore.%s %s %s\n",name,ok ? "PASS" : "FAIL",
           ok ? "reference quantities matched" : "reference quantity differed");
    failures+=!ok;
}
static const ChemSpecies species[]={{"a",1,1e-3,0,0},{"b",0,1e-3,0,0},{"c",0,1e-3,.025,1e-5}};
static const int reactants[]={0};
static const ChemProduct products[]={{1,1}};
static const char *names[]={"rate1","rate2","radius","number"};
static const double defaults[]={3,4,1e-7,2e6};
static double callback(double temperature, double pressure, double density, void *context)
{ (void)pressure; (void)density; (void)context; return temperature*1e-3; }

static int offsets(void)
{
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=1,.d=300}}, {0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=1,.d=300}},
        {0,1,0,1,0,{.kind=CHEM_TROE}}, {0,1,0,1,0,{.kind=CHEM_TUNNELING}},
        {0,1,0,1,0,{.kind=CHEM_USER}}, {0,1,0,1,0,{.kind=CHEM_SURFACE}}};
    const ChemMechanism m={3,6,0,species,reactions,reactants,products,NULL,NULL,300,101325};
    ChemRateIndex index={0}; int ok=chem_rate_index_init(&index,&m)==CHEM_OK;
    ok &= index.offsets[CHEM_TROE]==2 && index.offsets[CHEM_TERNARY]==3 &&
          index.offsets[CHEM_TERNARY_JPL19]==3 && index.offsets[CHEM_BRANCHED]==3 &&
          index.offsets[CHEM_TUNNELING]==3 && index.offsets[CHEM_TAYLOR]==4 &&
          index.offsets[CHEM_REVERSIBLE]==4 && index.offsets[CHEM_USER]==4 &&
          index.offsets[CHEM_SURFACE]==5 && index.offsets[CHEM_CALLBACK]==6 &&
          index.offsets[CHEM_CALLBACK+1]==6;
    for (size_t i=0;i<6 && ok;++i) ok &= index.reaction_indices[i]==i;
    chem_rate_index_destroy(&index); return ok;
}
static int arrhenius(void)
{
    const ChemReaction reaction={0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=2.15e-4,.b=1.2,.c=110,.d=300,.e=0}};
    const ChemMechanism m={3,1,0,species,&reaction,reactants,products,NULL,NULL,300,101325};
    ChemRateIndex index={0}; int ok=chem_rate_index_init(&index,&m)==CHEM_OK;
    const ChemRate *rate=&m.reactions[index.reaction_indices[0]].rate;
    ok &= index.offsets[CHEM_TROE]==1 && rate->a==2.15e-4 && rate->b==1.2 && rate->c==110 && rate->d==300;
    chem_rate_index_destroy(&index); return ok;
}
static int branched(void)
{
    ChemRate rate=chem_rate_defaults(CHEM_BRANCHED);
    rate.branch_x=1; rate.branch_y=0; rate.branch_a0=.5; rate.branch_n=3;
    double k0,z, expected_k0=2e-22*6.02214076e23*1e-6*exp(3.0);
    double density=2.45e19/6.02214076e23*1e6, a=expected_k0*density, b=.43*pow(293.0/298,-8);
    double expected_z=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)));
    return chem_branched_constants(&rate,&k0,&z)==CHEM_OK && fabs(k0-expected_k0)<=1e-10*expected_k0 &&
           fabs(z-expected_z)<=1e-10*fabs(expected_z);
}
static int user_indices(void)
{
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_USER,.parameter=0,.a=2}},
        {0,1,0,1,0,{.kind=CHEM_USER,.parameter=1,.a=.5}}};
    const ChemMechanism m={3,2,2,species,reactions,reactants,products,names,defaults,300,101325};
    ChemState *s=chem_state_create(&m,1); if (!s) return 0;
    int ok=reactions[0].rate.parameter==0 && reactions[0].rate.a==2 && reactions[1].rate.parameter==1 &&
           reactions[1].rate.a==.5 && chem_update_rates(s)==CHEM_OK &&
           chem_rate_constants(s,0)[0]==6 && chem_rate_constants(s,1)[0]==2;
    chem_state_free(s); return ok;
}
static int surface_fields(void)
{
    ChemRate rate=chem_rate_defaults(CHEM_SURFACE);
    rate.parameter=0; rate.diffusion=230; rate.molecular_weight=.025; rate.probability=.74;
    double factor, expected=8*(1.380649e-23*6.02214076e23)/(3.14159265358979323846*.025);
    return rate.diffusion==230 && rate.probability==.74 && rate.parameter==0 &&
           chem_surface_speed_factor(&rate,&factor)==CHEM_OK && fabs(factor-expected)<=1e-14;
}
static int surface_after_user(void)
{
    const char *labels[]={"ud","radius","number"}; const double values[]={2,1e-7,2e6};
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_USER,.parameter=0,.a=1}},
        {0,1,0,1,0,{.kind=CHEM_SURFACE,.parameter=1,.diffusion=1e-5,.molecular_weight=.025,.probability=1}}};
    const ChemMechanism m={3,2,3,species,reactions,reactants,products,labels,values,300,101325};
    ChemState *s=chem_state_create(&m,1); if (!s) return 0;
    int ok=reactions[0].rate.parameter==0 && reactions[1].rate.parameter==1 &&
           chem_update_rates(s)==CHEM_OK && chem_rate_constants(s,0)[0]==2 && chem_rate_constants(s,1)[0]>0;
    chem_state_free(s); return ok;
}
static int lambda_index(void)
{
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=1,.d=300}},
        {0,1,0,1,0,{.kind=CHEM_CALLBACK,.label="lam",.function=callback}}};
    const ChemMechanism m={3,2,0,species,reactions,reactants,products,NULL,NULL,300,101325};
    ChemRateIndex index={0}; ChemState *s=chem_state_create(&m,1);
    int ok=s && chem_rate_index_init(&index,&m)==CHEM_OK && index.offsets[CHEM_CALLBACK]==1 &&
           index.reaction_indices[1]==1 && reactions[1].rate.function && chem_update_rates(s)==CHEM_OK &&
           fabs(chem_rate_constants(s,1)[0]-.3)<=1e-14;
    chem_rate_index_destroy(&index); chem_state_free(s); return ok;
}
static int parameterized(void)
{
    const ChemParameterizedFunction factor={0,0,0,2};
    const ChemReaction reactions[]={
        {0,1,0,0,0,{.kind=CHEM_ARRHENIUS,.a=1,.d=300}},
        {0,0,0,1,0,{.kind=CHEM_ARRHENIUS,.a=1,.d=300,.parameterization_count=1,.parameterizations=&factor}}};
    const ChemProduct source[]={{0,1}};
    const ChemMechanism m={1,2,0,species,reactions,reactants,source,NULL,NULL,300,101325};
    ChemState *s=chem_state_create(&m,1); if (!s) return 0;
    int ok=chem_set_air_density(s,0,5)==CHEM_OK && chem_update_rates(s)==CHEM_OK &&
           fabs(chem_rate_constants(s,1)[0]-10)<=1e-14 && chem_rate_constants(s,0)[0]==1;
    ChemRate invalid=reactions[1].rate; double value;
    invalid.parameterization_count=5;
    ok &= chem_rate_calculate(&invalid,300,101325,5,NULL,0,&value)==CHEM_INVALID;
    chem_state_free(s); return ok;
}
static int total(void)
{
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_CALLBACK,.function=callback}}, {0,1,0,1,0,{.kind=CHEM_SURFACE}},
        {0,1,0,1,0,{.kind=CHEM_TROE}}, {0,1,0,1,0,{.kind=CHEM_USER}},
        {0,1,0,1,0,{.kind=CHEM_ARRHENIUS}}};
    const ChemMechanism m={3,5,0,species,reactions,reactants,products,NULL,NULL,300,101325};
    ChemRateIndex index={0}; int ok=chem_rate_index_init(&index,&m)==CHEM_OK;
    ok &= index.offsets[CHEM_CALLBACK]==4 && index.offsets[CHEM_CALLBACK+1]==5;
    const size_t expected[]={4,2,3,1,0};
    for (size_t i=0;i<5 && ok;++i) ok &= index.reaction_indices[i]==expected[i];
    chem_rate_index_destroy(&index); return ok;
}
static int missing(int diffusion)
{
    ChemRate rate=chem_rate_defaults(CHEM_SURFACE); rate.parameter=2; rate.probability=1;
    rate.diffusion=diffusion ? 0 : 1e-5; rate.molecular_weight=diffusion ? .025 : 0;
    double result;
    return chem_rate_calculate(&rate,300,101325,42,defaults,4,&result)==
           (diffusion ? CHEM_MISSING_DIFFUSION : CHEM_MISSING_MOLECULAR_WEIGHT);
}
int main(void)
{
    report("OffsetsAreContiguousCumulativeSizes",offsets()); report("ArrheniusParametersPreserved",arrhenius());
    report("BranchedDerivedFieldsComputed",branched()); report("UserDefinedCustomParamIndex",user_indices());
    report("SurfaceDataFieldsAndCustomParamIndex",surface_fields()); report("SurfaceCustomParamIndexAfterUserDefined",surface_after_user());
    report("LambdaEntriesRcIndex",lambda_index()); report("ParameterizedMultipliers",parameterized());
    report("TotalRateConstantsMatchesProcessCount",total());
    report("SurfaceMissingDiffusionCoefficientThrows",missing(1)); report("SurfaceMissingMolecularWeightThrows",missing(0));
    return failures ? 1 : 0;
}
