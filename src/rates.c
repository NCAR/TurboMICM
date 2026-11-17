#include "rates.h"
#include <math.h>

static double branch_a(double k0, double temperature, double air)
{
    double a = k0*air;
    double b = 0.43*pow(temperature/298.0, -8.0);
    return a/(1+a/b)*pow(0.41, 1/(1+pow(log10(a/b), 2)));
}

double chem_rate(const ChemRate *rate, ChemConditions c, const double *custom)
{
    const double *p = rate->p;
    double t = c.temperature;
    switch (rate->type) {
    case CHEM_RATE_ARRHENIUS:
        return p[0]*exp(p[2]/t)*pow(t/p[3], p[1])*(1+p[4]*c.pressure);
    case CHEM_RATE_TUNNELING:
        return p[0]*exp(-p[1]/t+p[2]/(t*t*t));
    case CHEM_RATE_TROE:
    case CHEM_RATE_TERNARY: {
        double k0 = p[0]*exp(p[2]/t)*pow(t/300, p[1]);
        double ki = p[3]*exp(p[5]/t)*pow(t/300, p[4]);
        double ratio = k0*c.air_density/ki;
        return k0*(rate->type == CHEM_RATE_TROE ? c.air_density : 1.0)/(1+ratio)*
            pow(p[6], p[7]/(p[7]+pow(log10(ratio),2)));
    }
    case CHEM_RATE_ALKOXY:
    case CHEM_RATE_NITRATE: {
        const double avogadro = 6.02214076e23;
        double k0 = 2e-22*avogadro*1e-6*exp(p[3]);
        double z = branch_a(k0, 293, 2.45e19/avogadro*1e6)*(1-p[2])/p[2];
        double a = branch_a(k0, t, c.air_density);
        double pre = p[0]*exp(-p[1]/t);
        return pre*(rate->type == CHEM_RATE_ALKOXY ? z/(z+a) : a/(a+z));
    }
    case CHEM_RATE_SURFACE: {
        const double pi = 3.14159265358979323846;
        const double gas_constant = 1.380649e-23*6.02214076e23;
        double radius = custom[rate->custom_index];
        double number = custom[rate->custom_index+1];
        double speed = sqrt(8*gas_constant*t/(pi*p[1]));
        return 4*number*pi*radius*radius/(radius/p[0]+4/(speed*p[2]));
    }
    case CHEM_RATE_USER:
        return custom[rate->custom_index]*p[0];
    }
    return NAN;
}

void chem_rates(size_t count, const ChemRate *rates, ChemConditions conditions,
                const double *custom, double *out)
{
    for (size_t i = 0; i < count; ++i) out[i] = chem_rate(&rates[i], conditions, custom);
}
