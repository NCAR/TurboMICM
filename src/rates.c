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
    case CHEM_RATE_REVERSIBLE:
        return p[0]*exp(p[1]/t)*p[2];
    case CHEM_RATE_TAYLOR: {
        double polynomial = 0, power = 1;
        if (rate->custom_index > 16) return NAN;
        for (size_t i = 0; i < rate->custom_index; ++i, power *= t)
            polynomial += p[5+i]*power;
        return polynomial*p[0]*exp(p[2]/t)*pow(t/p[3],p[1])*(1+p[4]*c.pressure);
    }
    case CHEM_RATE_TERNARY_JPL19: {
        double k0 = p[0]*exp(p[2]/t)*pow(t/p[3],p[1]);
        double ki = p[4]*exp(p[6]/t)*pow(t/p[7],p[5]);
        double kint = p[8]*exp(p[10]/t)*pow(t/p[11],p[9]);
        double k0m = k0*c.air_density, ratio = k0m/ki;
        double falloff = k0m/(1+ratio)*pow(p[12],p[13]/(p[13]+pow(log10(ratio),2)));
        return falloff+kint*(1-falloff/ki);
    }
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
    case CHEM_RATE_CALLBACK:
        return NAN; /* Callback bindings are supplied to the batch evaluator. */
    }
    return NAN;
}

void chem_rates(size_t count, const ChemRate *rates, ChemConditions conditions,
                const double *custom, double *out)
{
    for (size_t i = 0; i < count; ++i) out[i] = chem_rate(&rates[i], conditions, custom);
}

void chem_rates_with_air(size_t count, const ChemRate *rates, const unsigned *air_orders,
                         ChemConditions conditions, const double *custom, double *out)
{
    for (size_t i = 0; i < count; ++i) {
        double rate = chem_rate(&rates[i], conditions, custom);
        for (unsigned j = 0; j < air_orders[i]; ++j) rate *= conditions.air_density;
        out[i] = rate;
    }
}

void chem_rates_parameterized_with_callbacks(size_t count, const ChemRate *rates,
    const ChemRateMultiplier *multipliers, const ChemParameterizedSpecies *functions,
    const ChemRateBinding *bindings, size_t binding_count,
    ChemConditions conditions, const double *custom, double *out)
{
    for (size_t i=0;i<count;++i) {
        double rate;
        if (rates[i].type==CHEM_RATE_CALLBACK) {
            size_t index=rates[i].custom_index;
            rate=bindings && index<binding_count && bindings[index].function ?
                bindings[index].function(conditions,bindings[index].context) : NAN;
        } else rate=chem_rate(&rates[i],conditions,custom);
        double multiplier=1;
        for (size_t j=0;j<multipliers[i].count;++j) {
            const ChemParameterizedSpecies *p=functions+multipliers[i].begin+j;
            multiplier *= p->constant + p->temperature*conditions.temperature +
                          p->pressure*conditions.pressure + p->air_density*conditions.air_density;
        }
        out[i]=rate*multiplier;
    }
}

void chem_rates_parameterized(size_t count, const ChemRate *rates,
    const ChemRateMultiplier *multipliers, const ChemParameterizedSpecies *functions,
    ChemConditions conditions, const double *custom, double *out)
{
    chem_rates_parameterized_with_callbacks(count,rates,multipliers,functions,NULL,0,conditions,custom,out);
}
