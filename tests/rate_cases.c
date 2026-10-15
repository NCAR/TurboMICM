/* Case-by-case numerical assertions ported from the reference tests.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "chem.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(a,b,tol) do { double got_=(a), want_=(b), tol_=(tol); \
    if (!isfinite(got_) || fabs(got_-want_) > tol_) { \
        fprintf(stderr,"line %d: %.17g != %.17g (tol %.17g)\n",__LINE__,got_,want_,tol_); return 1; } } while(0)
static double calc(ChemRate r, double t, double air) { return chem_rate_evaluate(&r, t, 0, air, NULL); }
static int arrhenius(int taylor)
{
    ChemRate r = chem_rate_default(taylor ? CHEM_TAYLOR : CHEM_ARRHENIUS);
    double t = 301.24;
    CHECK(calc(r,t,0),1,1e-13);
    r.p[0]=1; CHECK(calc(r,t,0),1,1e-13);
    r.p[0]=2.2e-10; CHECK(calc(r,t,0),2.2e-10,1e-13*2.2e-10);
    r.p[0]=3e-11; r.p[2]=-200;
    double expected=3e-11*exp(-200/t);
    if (taylor) {
        r.coefficient_count=3; r.p[5]=12.5; r.p[6]=1.3e-2; r.p[7]=5.2e-4;
        expected *= (12.5+1.3e-2*t+5.2e-4*pow(t,2));
    }
    CHECK(calc(r,t,0),expected,1e-13*expected);
    r.p[0]=2.6e-12; r.p[2]=-350; expected=2.6e-12*exp(-350/t);
    if (taylor) { r.p[5]=1; r.p[6]=4.3e-1; r.p[7]=7.3e-3; expected *= (1+4.3e-1*t+7.3e-3*pow(t,2)); }
    CHECK(calc(r,t,0),expected,1e-13*expected);
    return 0;
}
static int falloff(int ternary, int all, int analytical)
{
    ChemRate r=chem_rate_default(ternary ? CHEM_TERNARY : CHEM_TROE);
    double t=301.24, air=42.2;
    if(all) {
        r.p[0]=1.2; r.p[1]=2.3; r.p[2]=302.3; r.p[3]=2.6; r.p[4]=-3.1; r.p[5]=402.1; r.p[6]=.9; r.p[7]=1.2;
    }
    if(analytical==1) r.p[0]=4e-18;
    if(analytical==2) {
        r.p[0]=1.2e-12; r.p[1]=167; r.p[2]=3; r.p[3]=136; r.p[4]=5; r.p[5]=24; r.p[6]=.9; r.p[7]=.8;
    }
    double k0=r.p[0]*exp(r.p[2]/t)*pow(t/300,r.p[1]);
    double kinf=r.p[3]*exp(r.p[5]/t)*pow(t/300,r.p[4]);
    double numerator=ternary ? k0 : air*k0;
    double expected=numerator/(1+air*k0/kinf)*pow(r.p[6],1/(1+1/r.p[7]*pow(log10(air*k0/kinf),2)));
    CHECK(calc(r,t,air),expected,analytical==2 ? 1e-14*fmax(fabs(calc(r,t,air)),fabs(expected)) : 0);
    return 0;
}
static int jpl19(int zero)
{
    ChemRate r=chem_rate_default(CHEM_TERNARY_JPL19);
    double t,air,k0,kinf,kint;
    if(zero) {
        t=250; air=30;
        r.p[0]=3e-3; r.p[1]=-2; r.p[4]=1.5; r.p[5]=-1; r.p[8]=0;
        k0=3e-3*pow(t/298,-2); kinf=1.5*pow(t/298,-1); kint=0;
    } else {
        t=301.24; air=42.2;
        const double p[]={1.2,2.3,302.3,298,2.6,-3.1,402.1,310,.7,.5,120,290,.9,1.2};
        memcpy(r.p,p,sizeof(p));
        k0=1.2*exp(302.3/t)*pow(t/298,2.3);
        kinf=2.6*exp(402.1/t)*pow(t/310,-3.1);
        kint=.7*exp(120/t)*pow(t/290,.5);
    }
    double ratio=k0*air/kinf;
    double kf=k0*air/(1+ratio)*pow(r.p[12],r.p[13]/(r.p[13]+pow(log10(ratio),2)));
    double expected=kf+kint*(1-kf/kinf);
    CHECK(calc(r,t,air),expected,zero ? 1e-30 : 0);
    return 0;
}
static int branched(int nitrate)
{
    ChemRate r={.type=nitrate?CHEM_BRANCHED_NITRATE:CHEM_BRANCHED_ALKOXY,.p={1.2,204.3,1e-3,2}};
    double t=301.24;
    double a=2e-22*exp(2)*2.45e19, b=.43*pow(293.0/298,-8);
    double z=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)))*(1-1e-3)/1e-3;
    a=2e-22*exp(2)*(42.2*6.02214076e23*1e-6); b=.43*pow(t/298,-8);
    double av=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)));
    double expected=1.2*exp(-204.3/t)*(nitrate?av/(z+av):z/(z+av));
    CHECK(calc(r,t,42.2),expected,1e-13*expected);
    return 0;
}
static int reversible(int variant)
{
    ChemRate r=chem_rate_default(CHEM_REVERSIBLE);
    double t=298.15;
    if(variant==0) { CHECK(calc(r,t,0),0,1e-13); return 0; }
    if(variant==1 || variant==2) {
        r.p[0]=1.14e-2; r.p[1]=2300; r.p[2]=.32;
        if(variant==1) t=301.24;
    } else if(variant==3) {
        r.p[0]=1.5; r.p[1]=2000; r.p[2]=.5;
        double temps[]={273.15,298.15,323.15};
        for(int i=0;i<3;++i) { double expected=1.5*exp(2000/temps[i])*.5; CHECK(calc(r,temps[i],0),expected,1e-13*expected); }
        return 0;
    } else if(variant==4) { r.p[0]=5; r.p[1]=1000; r.p[2]=0; }
    else { r.p[0]=3.5; r.p[1]=1800; r.p[2]=1; }
    double expected=r.p[0]*exp(r.p[1]/t)*r.p[2];
    CHECK(calc(r,t,0),expected,expected==0?1e-13:1e-13*expected);
    return 0;
}
static int surface(int variant)
{
    ChemRate rate={.type=CHEM_SURFACE,.p={2.3e2,.025,variant==1?.74:1},.parameter=0};
    if(variant>=2) {
        const char *names[]={"foo"}, *parameters[]={"radius","number"};
        int reactant=0;
        if(variant==2) rate.p[0]=0; else rate.p[1]=0;
        ChemReaction rx={.reactant_count=1,.reactants=&reactant,.rate=rate};
        ChemMechanismData data={.species_count=1,.reaction_count=1,.parameter_count=2,
            .species_names=names,.parameter_names=parameters,.reactions=&rx};
        ChemMechanism *m=NULL;
        CHECK(chem_mechanism_create(&data,&m),CHEM_INVALID,0);
        chem_mechanism_free(m);
        return 0;
    }
    double custom[]={1e-7,2.5e6}, pi=3.14159265358979323846;
    double expected=4*2.5e6*pi*pow(1e-7,2)/(1e-7/2.3e2+4/(sqrt(8*(1.380649e-23*6.02214076e23)*273.65/(pi*.025))*rate.p[2]));
    CHECK(chem_rate_evaluate(&rate,273.65,0,0,custom),expected,1e-13*expected);
    return 0;
}
static double callback(double t,double p,double air,const double *c,void *ctx)
{ (void)p;(void)air;(void)c;(void)ctx;return 2*t; }
static int lambda(void)
{
    ChemRate r={.type=CHEM_CALLBACK,.function=callback};
    CHECK(calc(r,300,0),600,1e-13*600);
    return 0;
}
static int tunneling(int all)
{
    ChemRate r=chem_rate_default(CHEM_TUNNELING);
    double expected=1,t=301.24;
    if(all) { r.p[0]=1.2;r.p[1]=2.3;r.p[2]=302.3;expected=1.2*exp(-2.3/t)*exp(302.3/pow(t,3)); }
    CHECK(calc(r,t,0),expected,1e-13*expected);
    return 0;
}
static int user(int variant)
{
    ChemRate r=chem_rate_default(CHEM_USER);
    double custom=variant==0?.5:variant==3?1.2:1.1;
    if(variant==3) r.p[0]=2;
    CHECK(chem_rate_evaluate(&r,300,0,0,&custom),custom*r.p[0],0);
    return 0;
}

typedef struct { const char *name; int family, variant; } Case;
static const Case cases[]={
 {"ArrheniusRateConstant.CalculateWithSystem",0,0},
 {"ArrheniusRateConstant.CalculateWithPrescribedArguments",0,0},
 {"TaylorSeriesRateConstant.CalculateWithSystem",0,1},
 {"TaylorSeriesRateConstant.CalculateWithPrescribedArguments",0,1},
 {"TroeRateConstant.CalculateWithMinimalArguments",1,0},
 {"TroeRateConstant.CalculateWithAllArguments",1,1},
 {"TroeRateConstant.AnalyticalTroeExampleAB",1,2},
 {"TroeRateConstant.AnalyticalTroeExampleBC",1,3},
 {"TernaryChemicalActivationRateConstant.CalculateWithMinimalArguments",2,0},
 {"TernaryChemicalActivationRateConstant.CalculateWithAllArguments",2,1},
 {"TernaryChemicalActivationJPL19RateConstant.CalculateWithAllArguments",3,0},
 {"TernaryChemicalActivationJPL19RateConstant.ZeroKintReducesToTroe",3,1},
 {"BranchedRateConstant.CalculateAlkoxyBranchWithAllArguments",4,0},
 {"BranchedRateConstant.CalculateNitrateBranchWithAllArguments",4,1},
 {"ReversibleRateConstant.DefaultConstructor",5,0},
 {"ReversibleRateConstant.CalculateWithSystem",5,1},
 {"ReversibleRateConstant.CalculateWithPrescribedArguments",5,2},
 {"ReversibleRateConstant.TemperatureDependence",5,3},
 {"ReversibleRateConstant.ZeroReverseRate",5,4},
 {"ReversibleRateConstant.EquilibriumConstantOnly",5,5},
 {"SurfaceRateConstant.CalculateDefaultProbability",6,0},
 {"SurfaceRateConstant.CalculateSpecifiedProbability",6,1},
 {"SurfaceRateConstant.DiffusionCoefficientIsMissing",6,2},
 {"SurfaceRateConstant.MolecularWeightIsMissing",6,3},
 {"LambdaRateConstant.CalculateWithSystem",7,0},
 {"TunnelingRateConstant.CalculateWithMinimalArguments",8,0},
 {"TunnelingRateConstant.CalculateWithAllArguments",8,1},
 {"UserDefinedRateConstant.CalculateWithSystem",9,0},
 {"UserDefinedRateConstant.ConstructorWithRate",9,1},
 {"UserDefinedRateConstant.ConstructorWithRateAndName",9,2},
 {"UserDefinedRateConstant.CustomScalingFactor",9,3},
};
int main(int argc,char **argv)
{
    if(argc!=2) return 2;
    for(size_t i=0;i<sizeof(cases)/sizeof(cases[0]);++i) {
        const Case *c=&cases[i];
        if(!strcmp(argv[1],"--list")) { puts(c->name);continue; }
        if(strcmp(argv[1],c->name)) continue;
        int v=c->variant;
        switch(c->family) {
        case 0:return arrhenius(v);
        case 1:return falloff(0,v==1,v>=2?v-1:0);
        case 2:return falloff(1,v,0);
        case 3:return jpl19(v);
        case 4:return branched(v);
        case 5:return reversible(v);
        case 6:return surface(v);
        case 7:return lambda();
        case 8:return tunneling(v);
        case 9:return user(v);
        }
    }
    return strcmp(argv[1],"--list")?2:0;
}
