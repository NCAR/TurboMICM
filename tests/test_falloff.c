#include "rates.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
int main(void)
{
    ChemConditions c = {301.24,101253.3,42.2};
    ChemRate r = {CHEM_RATE_ALKOXY,{1.2,204.3,1e-3,2},0};
    double avogadro=6.02214076e23;
    double a=2e-22*exp(2)*2.45e19;
    double b=.43*pow(293.0/298.0,-8);
    double z=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)))*(1-1e-3)/1e-3;
    a=2e-22*exp(2)*42.2*avogadro*1e-6;
    b=.43*pow(c.temperature/298.0,-8);
    double branch=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)));
    double expected=1.2*exp(-204.3/c.temperature)*z/(z+branch);
    assert(fabs(chem_rate(&r,c,NULL)-expected) <= 1e-13*expected);
    r.type=CHEM_RATE_NITRATE;
    expected=1.2*exp(-204.3/c.temperature)*branch/(z+branch);
    assert(fabs(chem_rate(&r,c,NULL)-expected) <= 1e-13*expected);
    r=(ChemRate){CHEM_RATE_TERNARY,{1,0,0,1,0,0,.6,1},0};
    expected=1.0/(1+42.2)*pow(.6,1/(1+pow(log10(42.2),2)));
    assert(chem_rate(&r,c,NULL)==expected);
    r.type=CHEM_RATE_TROE;
    expected=42.2/(1+42.2)*pow(.6,1/(1+pow(log10(42.2),2)));
    assert(chem_rate(&r,c,NULL)==expected);
    r=(ChemRate){CHEM_RATE_TERNARY,{1.2,2.3,302.3,2.6,-3.1,402.1,.9,1.2},0};
    double low=1.2*exp(302.3/301.24)*pow(301.24/300.0,2.3);
    double high=2.6*exp(402.1/301.24)*pow(301.24/300.0,-3.1);
    double q=low*42.2/high;
    expected=low/(1+q)*pow(.9,1/(1+(1/1.2)*pow(log10(q),2)));
    assert(chem_rate(&r,c,NULL)==expected);
    r.type=CHEM_RATE_TROE;
    expected=42.2*low/(1+42.2*low/high)*pow(.9,1/(1+(1/1.2)*pow(log10(42.2*low/high),2)));
    assert(chem_rate(&r,c,NULL)==expected);
    r=(ChemRate){CHEM_RATE_TROE,{4e-18,0,0,1,0,0,.6,1},0};
    double k0=4e-18, ratio=k0*42.2;
    expected=ratio/(1+ratio)*pow(.6,1/(1+pow(log10(ratio),2)));
    assert(chem_rate(&r,c,NULL)==expected);
    r=(ChemRate){CHEM_RATE_TROE,{1.2e-12,167,3,136,5,24,.9,.8},0};
    k0=1.2e-12*exp(3.0/301.24)*pow(301.24/300.0,167);
    double ki=136*exp(24.0/301.24)*pow(301.24/300.0,5);
    ratio=k0*42.2/ki;
    expected=k0*42.2/(1+ratio)*pow(.9,1/(1+(1/.8)*pow(log10(ratio),2)));
    assert(fabs(chem_rate(&r,c,NULL)-expected) <= 1e-14*fabs(expected));
    c.temperature=273.65;
    double custom[] = {1e-7,2.5e6};
    r=(ChemRate){CHEM_RATE_SURFACE,{230,.025,1},0};
    const double pi=3.14159265358979323846, gas=1.380649e-23*6.02214076e23;
    expected=4*2.5e6*pi*pow(1e-7,2)/(1e-7/230+4/sqrt(8*gas*273.65/(pi*.025)));
    assert(fabs(chem_rate(&r,c,custom)-expected) <= 1e-13*expected);
    r.p[2]=.74;
    expected=4*2.5e6*pi*pow(1e-7,2)/(1e-7/230+4/(sqrt(8*gas*273.65/(pi*.025))*.74));
    assert(fabs(chem_rate(&r,c,custom)-expected) <= 1e-13*expected);
    puts("falloff and branching cases PASS");
    return 0;
}
