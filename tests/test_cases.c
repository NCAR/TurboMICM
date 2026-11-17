#include "chem.h"
#include <math.h>
#include <stdio.h>

static double rate_at(ChemRateType type, const double *p, size_t n, double temperature)
{
    ChemReaction r = {0};
    r.type = type;
    for (size_t i = 0; i < n; ++i) r.p[i] = p[i];
    ChemMechanism m = {0};
    m.reaction_count = 1; m.reactions = &r;
    double value;
    if (chem_rate_constants(&m, (ChemConditions){temperature, 0}, NULL, &value)) return NAN;
    return value;
}

static double rate(ChemRateType type, const double *p, size_t n)
{
    return rate_at(type, p, n, 301.24);
}

static int reversible(double a, double c, double reverse, double temperature)
{
    double p[] = {a, c, reverse};
    double expected = a*exp(c/temperature)*reverse;
    return fabs(rate_at(CHEM_REVERSIBLE, p, 3, temperature)-expected) <= 1e-13*(expected == 0 ? 1 : fabs(expected));
}

static int user_rate(double value, double scale)
{
    ChemReaction reaction = {0}; reaction.type = CHEM_USER; reaction.parameter = 0; reaction.p[0] = scale;
    ChemMechanism m = {0}; m.reaction_count = 1; m.parameter_count = 1; m.reactions = &reaction;
    double result;
    return !chem_rate_constants(&m, (ChemConditions){301.24, 0}, &value, &result) && result == value*scale;
}

static int taylor(void)
{
    double p[22] = {1, 0, 0, 300, 0, 1, 1};
    if (fabs(rate(CHEM_TAYLOR, p, 22)-1) > 1e-13) return 0;
    p[0] = 1;
    if (fabs(rate(CHEM_TAYLOR, p, 22)-1) > 1e-13) return 0;
    p[0] = 2.2e-10;
    if (fabs(rate(CHEM_TAYLOR, p, 22)-p[0]) > 1e-13*p[0]) return 0;
    p[0] = 3e-11; p[2] = -200; p[5] = 3; p[6] = 12.5; p[7] = 1.3e-2; p[8] = 5.2e-4;
    double expected = p[0]*exp(p[2]/301.24)*(p[6]+p[7]*301.24+p[8]*pow(301.24,2));
    if (fabs(rate(CHEM_TAYLOR, p, 22)-expected) > 1e-13*expected) return 0;
    p[0] = 2.6e-12; p[2] = -350; p[6] = 1; p[7] = 4.3e-1; p[8] = 7.3e-3;
    expected = p[0]*exp(p[2]/301.24)*(p[6]+p[7]*301.24+p[8]*pow(301.24,2));
    return fabs(rate(CHEM_TAYLOR, p, 22)-expected) <= 1e-13*expected;
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
    printf("user_defined_rate_constant:UserDefinedRateConstant.CalculateWithSystem %d\n", user_rate(.5,1));
    printf("user_defined_rate_constant:UserDefinedRateConstant.ConstructorWithRate %d\n", user_rate(1.1,1));
    printf("user_defined_rate_constant:UserDefinedRateConstant.ConstructorWithRateAndName %d\n", user_rate(1.1,1));
    printf("user_defined_rate_constant:UserDefinedRateConstant.CustomScalingFactor %d\n", user_rate(1.2,2));
    printf("reversible_rate_constant:ReversibleRateConstant.DefaultConstructor %d\n", reversible(1,0,0,298.15));
    printf("reversible_rate_constant:ReversibleRateConstant.CalculateWithSystem %d\n", reversible(1.14e-2,2300,.32,301.24));
    printf("reversible_rate_constant:ReversibleRateConstant.CalculateWithPrescribedArguments %d\n", reversible(1.14e-2,2300,.32,298.15));
    printf("reversible_rate_constant:ReversibleRateConstant.TemperatureDependence %d\n", reversible(1.5,2000,.5,273.15) && reversible(1.5,2000,.5,298.15) && reversible(1.5,2000,.5,323.15));
    printf("reversible_rate_constant:ReversibleRateConstant.ZeroReverseRate %d\n", reversible(5,1000,0,298.15));
    printf("reversible_rate_constant:ReversibleRateConstant.EquilibriumConstantOnly %d\n", reversible(3.5,1800,1,298.15));
    printf("taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithSystem %d\n", taylor());
    printf("taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithPrescribedArguments %d\n", taylor());
    return 0;
}
