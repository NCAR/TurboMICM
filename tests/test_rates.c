#include "rates.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void)
{
    const ChemConditions c = {301.24,101253.3,1e6};
    ChemRate rate = {CHEM_RATE_ARRHENIUS,{1,0,0,300,0},0};
    assert(fabs(chem_rate(&rate,c,NULL)-1) <= 1e-13);
    rate.p[0] = 2.2e-10;
    assert(fabs(chem_rate(&rate,c,NULL)-2.2e-10) <= 1e-13*2.2e-10);
    rate.p[0] = 3e-11; rate.p[2] = -200;
    double expected = 3e-11*exp(-200/301.24);
    assert(fabs(chem_rate(&rate,c,NULL)-expected) <= 1e-13*expected);
    rate.p[0] = 2.6e-12; rate.p[2] = -350;
    expected = 2.6e-12*exp(-350/301.24);
    assert(fabs(chem_rate(&rate,c,NULL)-expected) <= 1e-13*expected);
    rate = (ChemRate){CHEM_RATE_TUNNELING,{1,0,0},0};
    assert(fabs(chem_rate(&rate,c,NULL)-1) <= 1e-13);
    rate = (ChemRate){CHEM_RATE_TUNNELING,{1.2,2.3,302.3},0};
    expected = 1.2*exp(-2.3/c.temperature)*exp(302.3/pow(c.temperature,3));
    assert(fabs(chem_rate(&rate,c,NULL)-expected) <= 1e-13*expected);
    double custom[] = {.5}; rate = (ChemRate){CHEM_RATE_USER,{1},0};
    assert(chem_rate(&rate,c,custom) == .5);
    custom[0] = 1.1; assert(chem_rate(&rate,c,custom) == 1.1);
    custom[0] = 1.2; rate.p[0] = 2;
    assert(chem_rate(&rate,c,custom) == 1.2*2.0);
    puts("reference rate constant cases PASS");
    return 0;
}
