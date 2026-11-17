#include "chem.h"

#include <math.h>
#include <stdio.h>

static int failures;
static void report(const char *id, int ok)
{
    printf("%s %s %s\n", id, ok ? "PASS" : "FAIL", ok ? "reference quantities matched" : "reference quantity differed");
    failures += !ok;
}
static int near(double actual, double expected, double tolerance)
{ return isfinite(actual) && fabs(actual-expected) <= tolerance; }
static double rate(ChemRate p, double t, double d, const double *custom, size_t n)
{
    double result;
    return chem_rate_calculate(&p, t, 0, d, custom, n, &result) == CHEM_OK ? result : NAN;
}
static double doubled_temperature(double t, double p, double d, void *context)
{ (void)p; (void)d; (void)context; return 2*t; }

static int arrhenius(void)
{
    ChemRate p = chem_rate_defaults(CHEM_ARRHENIUS);
    double t = 301.24, expected;
    if (!near(rate(p,t,0,NULL,0), 1, 1e-13)) return 0;
    p.a = 1; if (!near(rate(p,t,0,NULL,0), 1, 1e-13)) return 0;
    p.a = 2.2e-10; if (!near(rate(p,t,0,NULL,0), p.a, 1e-13*p.a)) return 0;
    p.a = 3e-11; p.c = -200; expected = p.a*exp(-200/t);
    if (!near(rate(p,t,0,NULL,0), expected, 1e-13*expected)) return 0;
    p.a = 2.6e-12; p.c = -350; expected = p.a*exp(-350/t);
    return near(rate(p,t,0,NULL,0), expected, 1e-13*expected);
}

static int branched(int nitrate)
{
    ChemRate p = chem_rate_defaults(CHEM_BRANCHED);
    p.branch_x=1.2; p.branch_y=204.3; p.branch_a0=1e-3; p.branch_n=2; p.nitrate_branch=nitrate;
    double t=301.24, density=42.2;
    double a=2e-22*exp(2)*2.45e19, b=.43*pow(293.0/298.0,-8);
    double z=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)))*(1-1e-3)/1e-3;
    a=2e-22*exp(2)*density*6.02214076e23*1e-6;
    b=.43*pow(t/298.0,-8);
    double A=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)));
    double expected=1.2*exp(-204.3/t)*(nitrate ? A/(z+A) : z/(z+A));
    return near(rate(p,t,density,NULL,0),expected,1e-13*expected);
}

static int reversible(int variant)
{
    ChemRate p = chem_rate_defaults(CHEM_REVERSIBLE);
    double t=298.15, expected;
    if (variant==0) return near(rate(p,t,0,NULL,0),0,1e-13);
    if (variant==1 || variant==2) {
        t=variant==1 ? 301.24 : 298.15; p.a=1.14e-2; p.c=2300; p.b=.32;
        expected=p.a*exp(p.c/t)*p.b;
        return near(rate(p,t,0,NULL,0),expected,1e-13*expected);
    }
    if (variant==3) {
        p.a=1.5; p.c=2000; p.b=.5;
        const double ts[]={273.15,298.15,323.15};
        for (int i=0;i<3;++i) {
            expected=1.5*exp(2000/ts[i])*.5;
            if (!near(rate(p,ts[i],0,NULL,0),expected,1e-13*expected)) return 0;
        }
        return 1;
    }
    if (variant==4) { p.a=5; p.c=1000; p.b=0; return near(rate(p,t,0,NULL,0),0,1e-13); }
    p.a=3.5; p.c=1800; p.b=1; expected=3.5*exp(1800/t);
    return near(rate(p,t,0,NULL,0),expected,1e-13*expected);
}

static int surface(double probability)
{
    const double pi=3.14159265358979323846, R=1.380649e-23*6.02214076e23;
    const double params[]={1e-7,2.5e6}, t=273.65;
    ChemRate p=chem_rate_defaults(CHEM_SURFACE);
    p.parameter=0; p.diffusion=230; p.molecular_weight=.025; p.probability=probability;
    double expected=4*2.5e6*pi*pow(1e-7,2)/(1e-7/230+4/(sqrt(8*R*t/(pi*.025))*probability));
    return near(rate(p,t,0,params,2),expected,1e-13*expected);
}

static int taylor(void)
{
    ChemRate p=chem_rate_defaults(CHEM_TAYLOR);
    double t=301.24, expected;
    if (!near(rate(p,t,0,NULL,0),1,1e-13)) return 0;
    p.a=1; if (!near(rate(p,t,0,NULL,0),1,1e-13)) return 0;
    p.a=2.2e-10; if (!near(rate(p,t,0,NULL,0),p.a,1e-13*p.a)) return 0;
    p.a=3e-11; p.c=-200; p.coefficient_count=3; p.coefficients[0]=12.5; p.coefficients[1]=1.3e-2; p.coefficients[2]=5.2e-4;
    expected=p.a*exp(-200/t)*(12.5+1.3e-2*t+5.2e-4*pow(t,2));
    if (!near(rate(p,t,0,NULL,0),expected,1e-13*expected)) return 0;
    p.a=2.6e-12; p.c=-350; p.coefficients[0]=1; p.coefficients[1]=4.3e-1; p.coefficients[2]=7.3e-3;
    expected=p.a*exp(-350/t)*(1+4.3e-1*t+7.3e-3*pow(t,2));
    return near(rate(p,t,0,NULL,0),expected,1e-13*expected);
}

static int falloff(int troe, int variant)
{
    double t=301.24, d=42.2, k0=1, ki=1, fc=.6, nn=1;
    ChemRate p=chem_rate_defaults(troe ? CHEM_TROE : CHEM_TERNARY);
    if (variant==1) {
        p.k0_a=1.2; p.k0_b=2.3; p.k0_c=302.3; p.kinf_a=2.6; p.kinf_b=-3.1; p.kinf_c=402.1; p.fc=.9; p.n=1.2;
        k0=1.2*exp(302.3/t)*pow(t/300,2.3); ki=2.6*exp(402.1/t)*pow(t/300,-3.1); fc=.9; nn=1.2;
    } else if (variant==2) { p.k0_a=k0=4e-18; }
    else if (variant==3) {
        p.k0_a=1.2e-12; p.k0_b=167; p.k0_c=3; p.kinf_a=136; p.kinf_b=5; p.kinf_c=24; p.fc=.9; p.n=.8;
        k0=1.2e-12*exp(3/t)*pow(t/300,167); ki=136*exp(24/t)*pow(t/300,5); fc=.9; nn=.8;
    }
    double ratio=k0*d/ki;
    double expected=(troe ? d*k0 : k0)/(1+ratio)*pow(fc,1/(1+(1/nn)*pow(log10(ratio),2)));
    double tolerance=variant==3 ? 1e-14*fmax(fabs(expected),fabs(rate(p,t,d,NULL,0))) : 0;
    return near(rate(p,t,d,NULL,0),expected,tolerance);
}

static int jpl(int zero)
{
    ChemRate p=chem_rate_defaults(CHEM_TERNARY_JPL19);
    double t, density;
    if (!zero) {
        t=301.24; density=42.2;
        p.k0_a=1.2; p.k0_b=2.3; p.k0_c=302.3; p.k0_d=298;
        p.kinf_a=2.6; p.kinf_b=-3.1; p.kinf_c=402.1; p.kinf_d=310;
        p.kint_a=.7; p.kint_b=.5; p.kint_c=120; p.kint_d=290; p.fc=.9; p.n=1.2;
    } else {
        t=250; density=30; p.k0_a=3e-3; p.k0_b=-2; p.kinf_a=1.5; p.kinf_b=-1; p.kint_a=0;
    }
    double k0=p.k0_a*exp(p.k0_c/t)*pow(t/p.k0_d,p.k0_b);
    double kinf=p.kinf_a*exp(p.kinf_c/t)*pow(t/p.kinf_d,p.kinf_b);
    double kint=p.kint_a*exp(p.kint_c/t)*pow(t/p.kint_d,p.kint_b);
    double ratio=k0*density/kinf;
    double kf=k0*density/(1+ratio)*pow(p.fc,p.n/(p.n+pow(log10(ratio),2)));
    double expected=kf+kint*(1-kf/kinf);
    return near(rate(p,t,density,NULL,0),expected,zero ? 1e-30 : 0);
}

int main(void)
{
    int a=arrhenius();
    report("arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithSystem",a);
    report("arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithPrescribedArguments",a);
    report("branched_rate_constant:BranchedRateConstant.CalculateAlkoxyBranchWithAllArguments",branched(0));
    report("branched_rate_constant:BranchedRateConstant.CalculateNitrateBranchWithAllArguments",branched(1));
    ChemRate callback=chem_rate_defaults(CHEM_CALLBACK); callback.function=doubled_temperature;
    report("lambda_rate_constant:LambdaRateConstant.CalculateWithSystem",near(rate(callback,300,0,NULL,0),600,6e-11));
    const char *rev[]={"DefaultConstructor","CalculateWithSystem","CalculateWithPrescribedArguments","TemperatureDependence","ZeroReverseRate","EquilibriumConstantOnly"};
    for (int i=0;i<6;++i) {
        char id[128]; snprintf(id,sizeof id,"reversible_rate_constant:ReversibleRateConstant.%s",rev[i]); report(id,reversible(i));
    }
    report("surface_rate_constant:SurfaceRateConstant.CalculateDefaultProbability",surface(1));
    report("surface_rate_constant:SurfaceRateConstant.CalculateSpecifiedProbability",surface(.74));
    ChemRate sf=chem_rate_defaults(CHEM_SURFACE); sf.parameter=0; sf.molecular_weight=.025;
    double custom[]={1e-7,2.5e6}, result;
    report("surface_rate_constant:SurfaceRateConstant.DiffusionCoefficientIsMissing",chem_rate_calculate(&sf,273.65,0,0,custom,2,&result)==CHEM_MISSING_DIFFUSION);
    sf.diffusion=230; sf.molecular_weight=0;
    report("surface_rate_constant:SurfaceRateConstant.MolecularWeightIsMissing",chem_rate_calculate(&sf,273.65,0,0,custom,2,&result)==CHEM_MISSING_MOLECULAR_WEIGHT);
    a=taylor();
    report("taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithSystem",a);
    report("taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithPrescribedArguments",a);
    report("ternary_chemical_activation_rate_constant:TernaryChemicalActivationRateConstant.CalculateWithMinimalArguments",falloff(0,0));
    report("ternary_chemical_activation_rate_constant:TernaryChemicalActivationRateConstant.CalculateWithAllArguments",falloff(0,1));
    report("ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.CalculateWithAllArguments",jpl(0));
    report("ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.ZeroKintReducesToTroe",jpl(1));
    const char *troe[]={"CalculateWithMinimalArguments","CalculateWithAllArguments","AnalyticalTroeExampleAB","AnalyticalTroeExampleBC"};
    for (int i=0;i<4;++i) {
        char id[128]; snprintf(id,sizeof id,"troe_rate_constant:TroeRateConstant.%s",troe[i]); report(id,falloff(1,i));
    }
    ChemRate tunneling=chem_rate_defaults(CHEM_TUNNELING); double t=301.24;
    report("tunneling_rate_constant:TunnelingRateConstant.CalculateWithMinimalArguments",near(rate(tunneling,t,0,NULL,0),1,1e-13));
    tunneling.a=1.2; tunneling.b=2.3; tunneling.c=302.3;
    double expected=1.2*exp(-2.3/t)*exp(302.3/pow(t,3));
    report("tunneling_rate_constant:TunnelingRateConstant.CalculateWithAllArguments",near(rate(tunneling,t,0,NULL,0),expected,1e-13*expected));
    ChemRate user=chem_rate_defaults(CHEM_USER); user.parameter=0;
    const char *ud[]={"CalculateWithSystem","ConstructorWithRate","ConstructorWithRateAndName","CustomScalingFactor"};
    for (int i=0;i<4;++i) {
        double v=i==0 ? .5 : i==3 ? 1.2 : 1.1; user.a=i==3 ? 2 : 1;
        char id[128]; snprintf(id,sizeof id,"user_defined_rate_constant:UserDefinedRateConstant.%s",ud[i]);
        report(id,near(rate(user,300,0,&v,1),v*user.a,0));
    }
    return failures ? 1 : 0;
}
