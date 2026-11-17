#include "rates.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>

static void near_relative(double actual, double expected, double tolerance)
{
    assert(fabs(actual-expected) <= tolerance*fabs(expected));
}
int main(void)
{
    ChemConditions c = {298.15,101253.3,42.2};
    ChemRate r = {CHEM_RATE_REVERSIBLE,{1,0,0},0};
    assert(chem_rate(&r,c,NULL) == 0);
    r.p[0]=1.14e-2; r.p[1]=2300; r.p[2]=.32; c.temperature=301.24;
    near_relative(chem_rate(&r,c,NULL),1.14e-2*exp(2300.0/301.24)*.32,1e-13);
    c.temperature=298.15;
    near_relative(chem_rate(&r,c,NULL),1.14e-2*exp(2300.0/298.15)*.32,1e-13);
    r.p[0]=1.5; r.p[1]=2000; r.p[2]=.5;
    const double temperatures[] = {273.15,298.15,323.15};
    for (size_t i=0;i<3;++i) {
        c.temperature=temperatures[i];
        near_relative(chem_rate(&r,c,NULL),1.5*exp(2000.0/c.temperature)*.5,1e-13);
    }
    r.p[0]=5; r.p[1]=1000; r.p[2]=0;
    assert(chem_rate(&r,c,NULL) == 0);
    r.p[0]=3.5; r.p[1]=1800; r.p[2]=1; c.temperature=298.15;
    near_relative(chem_rate(&r,c,NULL),3.5*exp(1800.0/298.15),1e-13);

    c.temperature=301.24;
    r=(ChemRate){CHEM_RATE_TAYLOR,{1,0,0,300,0,1},1};
    near_relative(chem_rate(&r,c,NULL),1,1e-13);
    r.p[0]=2.2e-10; near_relative(chem_rate(&r,c,NULL),2.2e-10,1e-13);
    r.p[0]=3e-11; r.p[2]=-200; r.p[5]=12.5; r.p[6]=1.3e-2; r.p[7]=5.2e-4; r.custom_index=3;
    double expected=3e-11*exp(-200/301.24)*(12.5+1.3e-2*301.24+5.2e-4*pow(301.24,2));
    near_relative(chem_rate(&r,c,NULL),expected,1e-13);
    r.p[0]=2.6e-12; r.p[2]=-350; r.p[5]=1; r.p[6]=.43; r.p[7]=.0073;
    expected=2.6e-12*exp(-350/301.24)*(1+.43*301.24+.0073*pow(301.24,2));
    near_relative(chem_rate(&r,c,NULL),expected,1e-13);

    r=(ChemRate){CHEM_RATE_TERNARY_JPL19,{1.2,2.3,302.3,298,2.6,-3.1,402.1,310,.7,.5,120,290,.9,1.2},0};
    double k0=1.2*exp(302.3/c.temperature)*pow(c.temperature/298,2.3);
    double ki=2.6*exp(402.1/c.temperature)*pow(c.temperature/310,-3.1);
    double kint=.7*exp(120.0/c.temperature)*pow(c.temperature/290,.5);
    double ratio=k0*c.air_density/ki;
    double kf=k0*c.air_density/(1+ratio)*pow(.9,1.2/(1.2+pow(log10(ratio),2)));
    expected=kf+kint*(1-kf/ki);
    assert(chem_rate(&r,c,NULL) == expected);
    c.temperature=250; c.air_density=30;
    r=(ChemRate){CHEM_RATE_TERNARY_JPL19,{.003,-2,0,298,1.5,-1,0,298,0,0,0,298,.6,1},0};
    k0=.003*pow(250.0/298,-2); ki=1.5*pow(250.0/298,-1); ratio=k0*30/ki;
    expected=k0*30/(1+ratio)*pow(.6,1/(1+pow(log10(ratio),2)));
    assert(fabs(chem_rate(&r,c,NULL)-expected) <= 1e-30);
    puts("additional rate constant cases PASS");
    return 0;
}
