/* Numerical formulas ported from MICM, Copyright UCAR 2023-2026, Apache-2.0. */
#include "chem.h"
#include <math.h>
#include <string.h>

void chem_ideal_air_density(chem_conditions *c) { c->air_density = c->pressure / (CHEM_GAS_CONSTANT * c->temperature); }

chem_rate chem_rate_default(chem_rate_type type)
{
    chem_rate r = {0};
    r.type = type;
    switch (type) {
    case CHEM_ARRHENIUS: r.p.arrhenius = (chem_arrhenius){1,0,0,300,0}; break;
    case CHEM_TROE: case CHEM_TERNARY: r.p.falloff = (chem_falloff){1,0,0,1,0,0,0.6,1}; break;
    case CHEM_TERNARY_JPL19: r.p.jpl19 = (chem_jpl19){1,0,0,298,1,0,0,298,1,0,0,298,0.6,1}; break;
    case CHEM_TUNNELING: r.p.tunneling.A = 1; break;
    case CHEM_TAYLOR: r.p.taylor.arrhenius = (chem_arrhenius){1,0,0,300,0}; r.p.taylor.n=1; r.p.taylor.coefficients[0]=1; break;
    case CHEM_REVERSIBLE: r.p.reversible.A = 1; break;
    case CHEM_USER: r.p.user.scaling = 1; break;
    case CHEM_SURFACE: r.p.surface.probability = 1; break;
    default: break;
    }
    return r;
}

chem_status chem_rate_prepare(chem_rate *r)
{
    if (r->type == CHEM_BRANCHED) {
        chem_branched *p = &r->p.branched;
        p->k0 = 2.0e-22 * 6.02214076e23 * 1.0e-6 * exp((double)p->n);
        double a = p->k0 * (2.45e19 / 6.02214076e23 * 1.0e6);
        double b = 0.43 * pow(293.0 / 298.0, -8.0);
        double A = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
        p->z = A * (1.0 - p->a0) / p->a0;
    } else if (r->type == CHEM_SURFACE) {
        chem_surface *p = &r->p.surface;
        if (p->diffusion <= 0 || p->molecular_weight <= 0) return CHEM_INVALID;
        p->speed_factor = 8.0 * CHEM_GAS_CONSTANT / (CHEM_PI * p->molecular_weight);
    } else if (r->type == CHEM_TAYLOR && (r->p.taylor.n == 0 || r->p.taylor.n > 16)) {
        return CHEM_INVALID;
    } else if (r->type == CHEM_LAMBDA && !r->p.lambda.function) return CHEM_INVALID;
    return CHEM_OK;
}

static double arrhenius(const chem_arrhenius *p, double T, double P)
{
    return p->A * exp(p->C / T) * pow(T / p->D, p->B) * (1.0 + p->E * P);
}

double chem_rate_calculate(const chem_rate *r, const chem_conditions *c, const double *custom)
{
    double T = c->temperature, M = c->air_density;
    switch (r->type) {
    case CHEM_ARRHENIUS: return arrhenius(&r->p.arrhenius, T, c->pressure);
    case CHEM_TROE: case CHEM_TERNARY: {
        const chem_falloff *p = &r->p.falloff;
        double k0 = p->k0_A * exp(p->k0_C / T) * pow(T / 300.0, p->k0_B);
        double kinf = p->kinf_A * exp(p->kinf_C / T) * pow(T / 300.0, p->kinf_B);
        double ratio = k0 * M / kinf;
        return k0 * (r->type == CHEM_TROE ? M : 1.0) / (1.0 + ratio)
            * pow(p->Fc, p->N / (p->N + pow(log10(ratio), 2.0)));
    }
    case CHEM_TERNARY_JPL19: {
        const chem_jpl19 *p = &r->p.jpl19;
        double k0 = p->k0_A * exp(p->k0_C / T) * pow(T / p->k0_D, p->k0_B);
        double kinf = p->kinf_A * exp(p->kinf_C / T) * pow(T / p->kinf_D, p->kinf_B);
        double kint = p->kint_A * exp(p->kint_C / T) * pow(T / p->kint_D, p->kint_B);
        double k0_M = k0 * M, ratio = k0_M / kinf;
        double kf = k0_M / (1.0 + ratio) * pow(p->Fc, p->N / (p->N + pow(log10(ratio), 2.0)));
        return kf + kint * (1.0 - kf / kinf);
    }
    case CHEM_TUNNELING: {
        const chem_tunneling *p = &r->p.tunneling;
        return p->A * exp(-p->B / T + p->C / (T * T * T));
    }
    case CHEM_BRANCHED: {
        const chem_branched *p = &r->p.branched;
        double a = p->k0 * M, b = 0.43 * pow(T / 298.0, -8.0);
        double A = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
        double pre = p->X * exp(-p->Y / T);
        return p->nitrate ? pre * (A / (A + p->z)) : pre * (p->z / (p->z + A));
    }
    case CHEM_TAYLOR: {
        const chem_taylor *p = &r->p.taylor;
        double poly = 0, tpow = 1;
        for (size_t j=0; j<p->n; ++j, tpow*=T) poly += p->coefficients[j] * tpow;
        const chem_arrhenius *a = &p->arrhenius;
        return poly * a->A * exp(a->C / T) * pow(T / a->D, a->B) * (1.0 + a->E * c->pressure);
    }
    case CHEM_REVERSIBLE: return r->p.reversible.A * exp(r->p.reversible.C / T) * r->p.reversible.k_r;
    case CHEM_USER: return custom[0] * r->p.user.scaling;
    case CHEM_SURFACE: {
        const chem_surface *p = &r->p.surface;
        double speed = sqrt(p->speed_factor * T), radius = custom[0], number = custom[1];
        return 4.0 * number * CHEM_PI * radius * radius / (radius / p->diffusion + 4.0 / (speed * p->probability));
    }
    case CHEM_LAMBDA: return r->p.lambda.function(c, custom, r->p.lambda.context);
    }
    return NAN;
}
