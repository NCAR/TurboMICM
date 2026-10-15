/* Numerical rate formulas ported from the specified reference revision.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "chem.h"
#include <math.h>

static const double gas_constant = 1.380649e-23 * 6.02214076e23;
static const double pi = 3.14159265358979323846;

double chem_parameterization_evaluate(const ChemParameterization *p, double t, double pressure, double air)
{ return p->constant + p->temperature * t + p->pressure * pressure + p->air_density * air; }

ChemRateDerived chem_rate_derived(const ChemRate *r)
{
    ChemRateDerived result = {0};
    if (r->type == CHEM_BRANCHED_ALKOXY || r->type == CHEM_BRANCHED_NITRATE) {
        result.branch_k0 = 2e-22 * 6.02214076e23 * 1e-6 * exp(r->p[3]);
        double ar = result.branch_k0 * (2.45e19 / 6.02214076e23 * 1e6);
        double br = .43 * pow(293.0 / 298.0, -8.0);
        double ref = ar / (1.0 + ar / br) * pow(.41, 1.0 / (1.0 + pow(log10(ar / br), 2.0)));
        result.branch_z = ref * (1.0 - r->p[2]) / r->p[2];
    }
    if (r->type == CHEM_SURFACE) result.surface_speed_factor = 8.0 * gas_constant / (pi * r->p[1]);
    return result;
}

static double arrhenius(const double *p, double t, double pressure)
{
    return p[0] * exp(p[2] / t) * pow(t / p[3], p[1]) * (1.0 + p[4] * pressure);
}

double chem_rate_evaluate(const ChemRate *r, double t, double pressure, double air, const double *custom)
{
    const double *p = r->p;
    double value = NAN;
    switch (r->type) {
    case CHEM_ARRHENIUS:
        value = arrhenius(p, t, pressure);
        break;
    case CHEM_TROE:
    case CHEM_TERNARY: {
        double k0 = p[0] * exp(p[2] / t) * pow(t / 300.0, p[1]);
        double kinf = p[3] * exp(p[5] / t) * pow(t / 300.0, p[4]);
        double ratio = k0 * air / kinf;
        double numerator = r->type == CHEM_TROE ? air : 1.0;
        value = k0 * numerator / (1.0 + ratio) * pow(p[6], p[7] / (p[7] + pow(log10(ratio), 2.0)));
        break;
    }
    case CHEM_TERNARY_JPL19: {
        double k0 = p[0] * exp(p[2] / t) * pow(t / p[3], p[1]);
        double kinf = p[4] * exp(p[6] / t) * pow(t / p[7], p[5]);
        double kint = p[8] * exp(p[10] / t) * pow(t / p[11], p[9]);
        double k0m = k0 * air, ratio = k0m / kinf;
        double kf = k0m / (1.0 + ratio) * pow(p[12], p[13] / (p[13] + pow(log10(ratio), 2.0)));
        value = kf + kint * (1.0 - kf / kinf);
        break;
    }
    case CHEM_TUNNELING:
        value = p[0] * exp(-p[1] / t + p[2] / (t * t * t));
        break;
    case CHEM_BRANCHED_ALKOXY:
    case CHEM_BRANCHED_NITRATE: {
        ChemRateDerived derived = chem_rate_derived(r);
        double k0 = derived.branch_k0, z = derived.branch_z;
        double a = k0 * air, b = .43 * pow(t / 298.0, -8.0);
        double av = a / (1.0 + a / b) * pow(.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
        double pre = p[0] * exp(-p[1] / t);
        value = r->type == CHEM_BRANCHED_ALKOXY ? pre * (z / (z + av)) : pre * (av / (av + z));
        break;
    }
    case CHEM_TAYLOR: {
        double poly = 0, power = 1;
        for (int i = 0; i < r->coefficient_count; ++i, power *= t)
            poly += p[5 + i] * power;
        value = poly * p[0] * exp(p[2] / t) * pow(t / p[3], p[1]) * (1.0 + p[4] * pressure);
        break;
    }
    case CHEM_REVERSIBLE:
        value = p[0] * exp(p[1] / t) * p[2];
        break;
    case CHEM_USER:
        value = custom[r->parameter] * p[0];
        break;
    case CHEM_SURFACE: {
        double radius = custom[r->parameter], number = custom[r->parameter + 1];
        double speed = sqrt(chem_rate_derived(r).surface_speed_factor * t);
        value = 4.0 * number * pi * radius * radius / (radius / p[0] + 4.0 / (speed * p[2]));
        break;
    }
    case CHEM_CALLBACK:
        if (r->function) value = r->function(t, pressure, air, custom, r->context);
        break;
    }
    for (int i = 0; i < r->third_body_count; ++i) value *= air;
    for (size_t i = 0; i < r->multiplier_count; ++i)
        value *= chem_parameterization_evaluate(&r->multipliers[i], t, pressure, air);
    return value;
}

ChemRate chem_rate_default(ChemRateType type)
{
    ChemRate r = {.type = type, .parameter = -1};
    switch (type) {
    case CHEM_ARRHENIUS: case CHEM_TAYLOR:
        r.p[0] = 1; r.p[3] = 300;
        if (type == CHEM_TAYLOR) { r.coefficient_count = 1; r.p[5] = 1; }
        break;
    case CHEM_TROE: case CHEM_TERNARY:
        r.p[0] = r.p[3] = r.p[7] = 1; r.p[6] = .6;
        break;
    case CHEM_TERNARY_JPL19:
        r.p[0] = r.p[4] = r.p[8] = r.p[13] = 1;
        r.p[3] = r.p[7] = r.p[11] = 298; r.p[12] = .6;
        break;
    case CHEM_TUNNELING: case CHEM_REVERSIBLE: case CHEM_USER:
        r.p[0] = 1;
        if (type == CHEM_USER) r.parameter = 0;
        break;
    case CHEM_SURFACE: r.p[2] = 1; r.parameter = 0; break;
    default: break;
    }
    return r;
}
