/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Rate formulas translated from MICM to C. */
#include "chem.h"
#include <math.h>
#include <string.h>

const char *chem_status_string(ChemStatus s)
{
    static const char *const names[] = {"success", "invalid argument", "allocation failed",
        "not found", "dimension mismatch", "singular matrix", "maximum steps exceeded",
        "step size too small", "constraint initialization failed"};
    return (unsigned)s < sizeof(names)/sizeof(*names) ? names[s] : "unknown status";
}
void chem_ideal_air_density(ChemConditions *c)
{ c->air_density = c->pressure / (CHEM_GAS_CONSTANT * c->temperature); }
double chem_parameter_air_density(const ChemConditions *c, void *context)
{ (void)context; return c->air_density; }

ChemRate chem_rate_default(ChemRateType type)
{
    ChemRate r;
    memset(&r, 0, sizeof(r));
    r.type = type;
    r.label = "";
    switch (type) {
    case CHEM_ARRHENIUS: r.data.arrhenius = (ChemArrhenius){1, 0, 0, 300, 0}; break;
    case CHEM_TROE: case CHEM_TERNARY:
        r.data.falloff = (ChemFalloff){1, 0, 0, 1, 0, 0, .6, 1}; break;
    case CHEM_TERNARY_JPL19:
        r.data.jpl19 = (ChemJpl19){1, 0, 0, 298, 1, 0, 0, 298, 1, 0, 0, 298, .6, 1}; break;
    case CHEM_TUNNELING: r.data.tunneling.A = 1; break;
    case CHEM_TAYLOR:
        r.data.taylor.arrhenius = (ChemArrhenius){1, 0, 0, 300, 0};
        r.data.taylor.count = 1; r.data.taylor.coefficients[0] = 1; break;
    case CHEM_REVERSIBLE: r.data.reversible.A = 1; break;
    case CHEM_USER: r.data.user.scaling = 1; break;
    case CHEM_SURFACE: r.data.surface.probability = 1; break;
    default: break;
    }
    return r;
}
ChemStatus chem_rate_prepare(ChemRate *r)
{
    if (!r || (unsigned)r->type >= CHEM_RATE_TYPE_COUNT) return CHEM_INVALID;
    if (r->type == CHEM_BRANCHED) {
        ChemBranched *p = &r->data.branched;
        p->k0 = 2e-22 * CHEM_AVOGADRO * 1e-6 * exp((double)p->n);
        double a = p->k0 * (2.45e19 / CHEM_AVOGADRO * 1e6);
        double b = .43 * pow(293.0 / 298.0, -8.0);
        double A = a / (1.0 + a / b) * pow(.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
        p->z = A * (1.0 - p->a0) / p->a0;
    } else if (r->type == CHEM_SURFACE) {
        ChemSurface *p = &r->data.surface;
        if (!(p->diffusion > 0) || !(p->molecular_weight > 0)) return CHEM_INVALID;
        p->speed_factor = 8.0 * CHEM_GAS_CONSTANT / (CHEM_PI * p->molecular_weight);
    } else if (r->type == CHEM_TAYLOR) {
        if (!r->data.taylor.count || r->data.taylor.count > 16) return CHEM_INVALID;
    } else if (r->type == CHEM_CALLBACK && !r->data.callback.function) return CHEM_INVALID;
    return CHEM_OK;
}
static double arrhenius(const ChemArrhenius *p, double t, double pressure)
{ return p->A * exp(p->C / t) * pow(t / p->D, p->B) * (1.0 + p->E * pressure); }
static double falloff(const ChemFalloff *p, double t, double density, double scale)
{
    double k0 = p->k0_A * exp(p->k0_C / t) * pow(t / 300.0, p->k0_B);
    double ki = p->kinf_A * exp(p->kinf_C / t) * pow(t / 300.0, p->kinf_B);
    double ratio = k0 * density / ki;
    return k0 * scale / (1.0 + ratio) * pow(p->Fc, p->N / (p->N + pow(log10(ratio), 2.0)));
}
double chem_rate_evaluate(const ChemRate *r, const ChemConditions *c, const double *custom)
{
    double t = c->temperature;
    switch (r->type) {
    case CHEM_ARRHENIUS: return arrhenius(&r->data.arrhenius, t, c->pressure);
    case CHEM_TROE: return falloff(&r->data.falloff, t, c->air_density, c->air_density);
    case CHEM_TERNARY: return falloff(&r->data.falloff, t, c->air_density, 1.0);
    case CHEM_TERNARY_JPL19: {
        const ChemJpl19 *p = &r->data.jpl19;
        double k0 = p->k0_A * exp(p->k0_C / t) * pow(t / p->k0_D, p->k0_B);
        double ki = p->kinf_A * exp(p->kinf_C / t) * pow(t / p->kinf_D, p->kinf_B);
        double km = p->kint_A * exp(p->kint_C / t) * pow(t / p->kint_D, p->kint_B);
        double k0m = k0 * c->air_density, ratio = k0m / ki;
        double kf = k0m / (1.0 + ratio) * pow(p->Fc, p->N / (p->N + pow(log10(ratio), 2.0)));
        return kf + km * (1.0 - kf / ki);
    }
    case CHEM_BRANCHED: {
        const ChemBranched *p = &r->data.branched;
        double a = p->k0 * c->air_density, b = .43 * pow(t / 298.0, -8.0);
        double A = a / (1.0 + a / b) * pow(.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
        double pre = p->X * exp(-p->Y / t);
        return p->branch == 0 ? pre * (p->z / (p->z + A)) : pre * (A / (A + p->z));
    }
    case CHEM_TUNNELING: {
        const ChemTunneling *p = &r->data.tunneling;
        return p->A * exp(-p->B / t + p->C / (t * t * t));
    }
    case CHEM_TAYLOR: {
        const ChemTaylor *p = &r->data.taylor;
        double poly = 0, power = 1;
        for (size_t i = 0; i < p->count; ++i, power *= t) poly += p->coefficients[i] * power;
        const ChemArrhenius *a = &p->arrhenius;
        return poly * a->A * exp(a->C / t) * pow(t / a->D, a->B) * (1.0 + a->E * c->pressure);
    }
    case CHEM_REVERSIBLE:
        return r->data.reversible.A * exp(r->data.reversible.C / t) * r->data.reversible.reverse_rate;
    case CHEM_USER: return custom[0] * r->data.user.scaling;
    case CHEM_SURFACE: {
        const ChemSurface *p = &r->data.surface;
        double speed = sqrt(p->speed_factor * t), radius = custom[0], count = custom[1];
        return 4.0 * count * CHEM_PI * radius * radius /
            (radius / p->diffusion + 4.0 / (speed * p->probability));
    }
    case CHEM_CALLBACK: return r->data.callback.function(c, r->data.callback.context);
    default: return NAN;
    }
}
