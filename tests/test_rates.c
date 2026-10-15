/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Port of test/unit/process/rate_constant numerical checks. */
#include "test.h"

static void arrhenius_test(void)
{
    ChemConditions c = {.temperature=301.24};
    ChemRate r = chem_rate_default(CHEM_ARRHENIUS);
    REL(chem_rate_evaluate(&r,&c,NULL),1,1e-13);
    r.data.arrhenius.A=1;
    REL(chem_rate_evaluate(&r,&c,NULL),1,1e-13);
    r.data.arrhenius.A=2.2e-10;
    REL(chem_rate_evaluate(&r,&c,NULL),2.2e-10,1e-13);
    r.data.arrhenius.A=3e-11; r.data.arrhenius.C=-200;
    REL(chem_rate_evaluate(&r,&c,NULL),3e-11*exp(-200/301.24),1e-13);
    r.data.arrhenius.A=2.6e-12; r.data.arrhenius.C=-350;
    REL(chem_rate_evaluate(&r,&c,NULL),2.6e-12*exp(-350/301.24),1e-13);
}
static void branched_test(int branch)
{
    ChemConditions c = {.temperature=301.24,.air_density=42.2};
    ChemRate r=chem_rate_default(CHEM_BRANCHED);
    r.data.branched=(ChemBranched){.branch=branch,.X=1.2,.Y=204.3,.a0=1e-3,.n=2};
    CHECK(chem_rate_prepare(&r)==CHEM_OK);
    double a=2e-22*exp(2)*2.45e19, b=.43*pow(293.0/298,-8);
    double z=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)))*(1-1e-3)/1e-3;
    a=2e-22*exp(2)*(42.2*CHEM_AVOGADRO*1e-6); b=.43*pow(c.temperature/298,-8);
    double A=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)));
    double expected=1.2*exp(-204.3/c.temperature)*(branch ? A/(z+A) : z/(z+A));
    REL(chem_rate_evaluate(&r,&c,NULL),expected,1e-13);
}
static void falloff_test(ChemRateType type,int full)
{
    ChemConditions c={.temperature=301.24,.air_density=42.2};
    ChemRate r=chem_rate_default(type);
    double expected,k0=1,ki=1;
    if (full) {
        r.data.falloff=(ChemFalloff){1.2,2.3,302.3,2.6,-3.1,402.1,.9,1.2};
        k0=1.2*exp(302.3/c.temperature)*pow(c.temperature/300,2.3);
        ki=2.6*exp(402.1/c.temperature)*pow(c.temperature/300,-3.1);
        expected=(type==CHEM_TROE?42.2*k0:k0)/(1+42.2*k0/ki)*
            pow(.9,1/(1+1/1.2*pow(log10(42.2*k0/ki),2)));
    } else {
        expected=(type==CHEM_TROE?42.2*k0:k0)/(1+42.2*k0/ki)*
            pow(.6,1/(1+pow(log10(42.2*k0/ki),2)));
    }
    NEAR(chem_rate_evaluate(&r,&c,NULL),expected,0);
}
static void troe_analytical(int bc)
{
    ChemConditions c={.temperature=301.24,.air_density=42.2};
    ChemRate r=chem_rate_default(CHEM_TROE);
    double k0,ki,expected;
    if (!bc) {
        r.data.falloff.k0_A=4e-18; k0=4e-18; ki=1;
        expected=k0*42.2/(1+k0*42.2/ki)*pow(.6,1/(1+1/1.0*pow(log10(k0*42.2/ki),2)));
        NEAR(chem_rate_evaluate(&r,&c,NULL),expected,0);
    } else {
        r.data.falloff=(ChemFalloff){1.2e-12,167,3,136,5,24,.9,.8};
        k0=1.2e-12*exp(3/301.24)*pow(301.24/300,167);
        ki=136*exp(24/301.24)*pow(301.24/300,5);
        expected=k0*42.2/(1+k0*42.2/ki)*pow(.9,1/(1+(1/.8)*pow(log10(k0*42.2/ki),2)));
        double actual=chem_rate_evaluate(&r,&c,NULL);
        CHECK(fabs(actual-expected)/fmax(fabs(actual),fabs(expected))<=1e-14);
    }
}
static void jpl_test(int zero)
{
    ChemConditions c={.temperature=301.24,.air_density=42.2};
    ChemRate r=chem_rate_default(CHEM_TERNARY_JPL19);
    double k0,ki,km=0;
    if (zero) {
        c.temperature=250; c.air_density=30;
        r.data.jpl19.k0_A=3e-3; r.data.jpl19.k0_B=-2;
        r.data.jpl19.kinf_A=1.5; r.data.jpl19.kinf_B=-1; r.data.jpl19.kint_A=0;
        k0=3e-3*pow(250.0/298,-2); ki=1.5*pow(250.0/298,-1);
    } else {
        r.data.jpl19=(ChemJpl19){1.2,2.3,302.3,298,2.6,-3.1,402.1,310,.7,.5,120,290,.9,1.2};
        k0=1.2*exp(302.3/c.temperature)*pow(c.temperature/298,2.3);
        ki=2.6*exp(402.1/c.temperature)*pow(c.temperature/310,-3.1);
        km=.7*exp(120/c.temperature)*pow(c.temperature/290,.5);
    }
    double ratio=k0*c.air_density/ki;
    double kf=k0*c.air_density/(1+ratio)*(zero?pow(.6,1/(1+pow(log10(ratio),2))):pow(.9,1.2/(1.2+pow(log10(ratio),2))));
    NEAR(chem_rate_evaluate(&r,&c,NULL),kf+km*(1-kf/ki),zero?1e-30:0);
}
static void tunneling_test(int full)
{
    ChemConditions c={.temperature=301.24};
    ChemRate r=chem_rate_default(CHEM_TUNNELING);
    double expected=1;
    if (full) {
        r.data.tunneling=(ChemTunneling){1.2,2.3,302.3};
        expected=1.2*exp(-2.3/c.temperature)*exp(302.3/pow(c.temperature,3));
    }
    REL(chem_rate_evaluate(&r,&c,NULL),expected,1e-13);
}
static void taylor_test(void)
{
    ChemConditions c={.temperature=301.24};
    ChemRate r=chem_rate_default(CHEM_TAYLOR);
    ChemTaylor *p=&r.data.taylor;
    REL(chem_rate_evaluate(&r,&c,NULL),1,1e-13);
    p->arrhenius.A=1; REL(chem_rate_evaluate(&r,&c,NULL),1,1e-13);
    p->arrhenius.A=2.2e-10; REL(chem_rate_evaluate(&r,&c,NULL),2.2e-10,1e-13);
    p->arrhenius.A=3e-11; p->arrhenius.C=-200;
    p->count=3; p->coefficients[0]=12.5; p->coefficients[1]=1.3e-2; p->coefficients[2]=5.2e-4;
    double expected=3e-11*exp(-200/301.24)*(12.5+1.3e-2*301.24+5.2e-4*pow(301.24,2));
    REL(chem_rate_evaluate(&r,&c,NULL),expected,1e-13);
    p->arrhenius.A=2.6e-12; p->arrhenius.C=-350;
    p->coefficients[0]=1; p->coefficients[1]=4.3e-1; p->coefficients[2]=7.3e-3;
    expected=2.6e-12*exp(-350/301.24)*(1+4.3e-1*301.24+7.3e-3*pow(301.24,2));
    REL(chem_rate_evaluate(&r,&c,NULL),expected,1e-13);
}
static void reversible_test(double A,double C,double kr,double temperature)
{
    ChemConditions c={.temperature=temperature};
    ChemRate r=chem_rate_default(CHEM_REVERSIBLE);
    r.data.reversible=(ChemReversible){A,C,kr};
    double expected=A*exp(C/temperature)*kr;
    NEAR(chem_rate_evaluate(&r,&c,NULL),expected,1e-13*(expected?expected:1));
}
static void reversible_default(void)
{
    ChemConditions c={.temperature=298.15};
    ChemRate r=chem_rate_default(CHEM_REVERSIBLE);
    NEAR(chem_rate_evaluate(&r,&c,NULL),0,1e-13);
}
static void reversible_temperature(void)
{
    reversible_test(1.5,2000,.5,273.15);
    reversible_test(1.5,2000,.5,298.15);
    reversible_test(1.5,2000,.5,323.15);
}
static void surface_test(double probability)
{
    ChemConditions c={.temperature=273.65};
    ChemRate r=chem_rate_default(CHEM_SURFACE);
    r.data.surface.diffusion=2.3e2; r.data.surface.molecular_weight=.025;
    if (probability!=1) r.data.surface.probability=probability;
    CHECK(chem_rate_prepare(&r)==CHEM_OK);
    double custom[]={1e-7,2.5e6};
    double expected=4*2.5e6*CHEM_PI*pow(1e-7,2)/(1e-7/2.3e2+4/(sqrt(8*CHEM_GAS_CONSTANT*273.65/(CHEM_PI*.025))*probability));
    REL(chem_rate_evaluate(&r,&c,custom),expected,1e-13);
}
static void surface_missing(int weight)
{
    ChemSpecies sp={.name="foo"};
    ChemRate r=chem_rate_default(CHEM_SURFACE);
    r.label="foo";
    r.data.surface.diffusion=weight?2.3e2:0;
    r.data.surface.molecular_weight=weight?0:.025;
    const size_t reactant=0;
    ChemReaction reaction={.reactant_count=1,.reactants=&reactant,.rate=r};
    ChemMechanism *m=NULL;
    CHECK(chem_mechanism_create(1,&sp,1,&reaction,&m)==CHEM_INVALID);
    CHECK(m==NULL); chem_mechanism_free(m);
}
static void user_test(double v,double scale,int named)
{
    ChemRate r=chem_rate_default(CHEM_USER);
    ChemConditions c={0};
    if (named) r.label="test";
    r.data.user.scaling=scale;
    NEAR(chem_rate_evaluate(&r,&c,&v),v*scale,0);
}
static double lambda(const ChemConditions *c,void *context)
{ (void)context; return 2*c->temperature; }
static void lambda_test(void)
{
    ChemRate r=chem_rate_default(CHEM_CALLBACK);
    r.label="test lambda"; r.data.callback.function=lambda;
    ChemConditions c={.temperature=300};
    REL(chem_rate_evaluate(&r,&c,NULL),600,1e-13);
}
void test_rates(void)
{
#define R(file,suite,test,call) RUN(file ":" suite "." test,call)
    R("arrhenius_rate_constant","ArrheniusRateConstant","CalculateWithSystem",arrhenius_test());
    R("arrhenius_rate_constant","ArrheniusRateConstant","CalculateWithPrescribedArguments",arrhenius_test());
    R("branched_rate_constant","BranchedRateConstant","CalculateAlkoxyBranchWithAllArguments",branched_test(0));
    R("branched_rate_constant","BranchedRateConstant","CalculateNitrateBranchWithAllArguments",branched_test(1));
    R("troe_rate_constant","TroeRateConstant","CalculateWithMinimalArguments",falloff_test(CHEM_TROE,0));
    R("troe_rate_constant","TroeRateConstant","CalculateWithAllArguments",falloff_test(CHEM_TROE,1));
    R("troe_rate_constant","TroeRateConstant","AnalyticalTroeExampleAB",troe_analytical(0));
    R("troe_rate_constant","TroeRateConstant","AnalyticalTroeExampleBC",troe_analytical(1));
    R("ternary_chemical_activation_rate_constant","TernaryChemicalActivationRateConstant","CalculateWithMinimalArguments",falloff_test(CHEM_TERNARY,0));
    R("ternary_chemical_activation_rate_constant","TernaryChemicalActivationRateConstant","CalculateWithAllArguments",falloff_test(CHEM_TERNARY,1));
    R("ternary_chemical_activation_rate_constant","TernaryChemicalActivationJPL19RateConstant","CalculateWithAllArguments",jpl_test(0));
    R("ternary_chemical_activation_rate_constant","TernaryChemicalActivationJPL19RateConstant","ZeroKintReducesToTroe",jpl_test(1));
    R("tunneling_rate_constant","TunnelingRateConstant","CalculateWithMinimalArguments",tunneling_test(0));
    R("tunneling_rate_constant","TunnelingRateConstant","CalculateWithAllArguments",tunneling_test(1));
    R("taylor_series_rate_constant","TaylorSeriesRateConstant","CalculateWithSystem",taylor_test());
    R("taylor_series_rate_constant","TaylorSeriesRateConstant","CalculateWithPrescribedArguments",taylor_test());
    R("reversible_rate_constant","ReversibleRateConstant","DefaultConstructor",reversible_default());
    R("reversible_rate_constant","ReversibleRateConstant","CalculateWithSystem",reversible_test(1.14e-2,2300,.32,301.24));
    R("reversible_rate_constant","ReversibleRateConstant","CalculateWithPrescribedArguments",reversible_test(1.14e-2,2300,.32,298.15));
    R("reversible_rate_constant","ReversibleRateConstant","TemperatureDependence",reversible_temperature());
    R("reversible_rate_constant","ReversibleRateConstant","ZeroReverseRate",reversible_test(5,1000,0,298.15));
    R("reversible_rate_constant","ReversibleRateConstant","EquilibriumConstantOnly",reversible_test(3.5,1800,1,298.15));
    R("surface_rate_constant","SurfaceRateConstant","CalculateDefaultProbability",surface_test(1));
    R("surface_rate_constant","SurfaceRateConstant","CalculateSpecifiedProbability",surface_test(.74));
    R("surface_rate_constant","SurfaceRateConstant","DiffusionCoefficientIsMissing",surface_missing(0));
    R("surface_rate_constant","SurfaceRateConstant","MolecularWeightIsMissing",surface_missing(1));
    R("user_defined_rate_constant","UserDefinedRateConstant","CalculateWithSystem",user_test(.5,1,0));
    R("user_defined_rate_constant","UserDefinedRateConstant","ConstructorWithRate",user_test(1.1,1,0));
    R("user_defined_rate_constant","UserDefinedRateConstant","ConstructorWithRateAndName",user_test(1.1,1,1));
    R("user_defined_rate_constant","UserDefinedRateConstant","CustomScalingFactor",user_test(1.2,2,0));
    R("lambda_rate_constant","LambdaRateConstant","CalculateWithSystem",lambda_test());
#undef R
}
