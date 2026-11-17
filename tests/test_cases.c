#include "chem.h"
#include <math.h>
#include <stdio.h>

static double rate(ChemRateType type, const double *p, size_t n)
{
    ChemReaction r = {0};
    r.type = type;
    for (size_t i = 0; i < n; ++i) r.p[i] = p[i];
    ChemMechanism m = {0};
    m.reaction_count = 1; m.reactions = &r;
    double value;
    if (chem_rate_constants(&m, (ChemConditions){301.24, 0}, NULL, &value)) return NAN;
    return value;
}

static int arrhenius(void)
{
    double p[] = {1, 0, 0, 300, 0};
    if (fabs(rate(CHEM_ARRHENIUS, p, 5) - 1) > 1e-13) return 0;
    p[0] = 1;
    if (fabs(rate(CHEM_ARRHENIUS, p, 5) - 1) > 1e-13) return 0;
    p[0] = 2.2e-10;
    if (fabs(rate(CHEM_ARRHENIUS, p, 5) - p[0]) > 1e-13*p[0]) return 0;
    p[0] = 3e-11; p[2] = -200;
    double expected = p[0]*exp(-200/301.24);
    if (fabs(rate(CHEM_ARRHENIUS, p, 5)-expected) > 1e-13*expected) return 0;
    p[0] = 2.6e-12; p[2] = -350;
    expected = p[0]*exp(-350/301.24);
    return fabs(rate(CHEM_ARRHENIUS, p, 5)-expected) <= 1e-13*expected;
}

static int tunneling(int full)
{
    double p[] = {1, 0, 0};
    if (full) { p[0] = 1.2; p[1] = 2.3; p[2] = 302.3; }
    double expected = p[0]*exp(-p[1]/301.24)*exp(p[2]/pow(301.24, 3));
    return fabs(rate(CHEM_TUNNELING, p, 3)-expected) <= 1e-13*expected;
}

int main(void)
{
    printf("arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithSystem %d\n", arrhenius());
    printf("arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithPrescribedArguments %d\n", arrhenius());
    printf("tunneling_rate_constant:TunnelingRateConstant.CalculateWithMinimalArguments %d\n", tunneling(0));
    printf("tunneling_rate_constant:TunnelingRateConstant.CalculateWithAllArguments %d\n", tunneling(1));
    return 0;
}
