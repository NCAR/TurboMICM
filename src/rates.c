/* Numerical rate formulas ported from MICM.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

void ch_ideal_air_density(ch_conditions *c) {
  if (c) c->air_density = c->pressure / (CH_GAS_CONSTANT * c->temperature);
}

ch_rate ch_rate_default(int type) {
  ch_rate r = {0};
  r.type = type;
  r.custom_index = -1;
  r.p[0] = 1.0;
  if (type == CH_ARRHENIUS || type == CH_TAYLOR) r.p[3] = 300.0;
  if (type == CH_TAYLOR) { r.coefficient_count = 1; r.p[5] = 1.0; }
  if (type == CH_TROE || type == CH_TERNARY) {
    r.p[3] = 1.0; r.p[6] = 0.6; r.p[7] = 1.0;
  }
  if (type == CH_TERNARY_JPL19) {
    r.p[3] = r.p[7] = r.p[11] = 298.0;
    r.p[4] = r.p[8] = r.p[13] = 1.0;
    r.p[12] = 0.6;
  }
  if (type == CH_SURFACE) r.p[2] = 1.0;
  return r;
}

void ch_rate_prepare(ch_rate *r) {
  if (r->type == CH_BRANCHED) {
    double *p = r->p;
    p[5] = 2.0e-22 * CH_AVOGADRO * 1.0e-6 * exp(p[3]);
    double a = p[5] * (2.45e19 / CH_AVOGADRO * 1.0e6);
    double b = 0.43 * pow(293.0 / 298.0, -8.0);
    double v = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a/b), 2.0)));
    p[6] = v * (1.0 - p[2]) / p[2];
  }
}

double ch_rate_calculate(const ch_rate *r, const ch_conditions *cond, const double *custom) {
  const double *p = r->p;
  double t = cond->temperature, air = cond->air_density;
  switch (r->type) {
    case CH_ARRHENIUS:
      return p[0] * exp(p[2]/t) * pow(t/p[3], p[1]) * (1.0 + p[4]*cond->pressure);
    case CH_TROE:
    case CH_TERNARY: {
      double k0 = p[0] * exp(p[2]/t) * pow(t/300.0, p[1]);
      double ki = p[3] * exp(p[5]/t) * pow(t/300.0, p[4]);
      double ratio = k0 * air / ki;
      return k0 * (r->type == CH_TROE ? air : 1.0) / (1.0 + ratio)
             * pow(p[6], p[7] / (p[7] + pow(log10(ratio), 2.0)));
    }
    case CH_TERNARY_JPL19: {
      double k0 = p[0] * exp(p[2]/t) * pow(t/p[3], p[1]);
      double ki = p[4] * exp(p[6]/t) * pow(t/p[7], p[5]);
      double km = p[8] * exp(p[10]/t) * pow(t/p[11], p[9]);
      double ratio = k0 * air / ki;
      double kf = k0 * air / (1.0 + ratio) * pow(p[12], p[13]/(p[13]+pow(log10(ratio),2.0)));
      return kf + km * (1.0 - kf / ki);
    }
    case CH_TUNNELING: return p[0] * exp(-p[1]/t + p[2]/(t*t*t));
    case CH_BRANCHED: {
      double a = p[5]*air, b = 0.43*pow(t/298.0,-8.0);
      double v = a/(1.0+a/b) * pow(0.41,1.0/(1.0+pow(log10(a/b),2.0)));
      double pre = p[0]*exp(-p[1]/t);
      return p[4] == 0.0 ? pre*(p[6]/(p[6]+v)) : pre*(v/(v+p[6]));
    }
    case CH_SURFACE: {
      double radius = custom[r->custom_index], number = custom[r->custom_index+1];
      return 4.0 * number * 3.14159265358979323846 * radius * radius /
        (radius/p[0] + 4.0/(sqrt(p[1]*t)*p[2]));
    }
    case CH_USER: return custom[r->custom_index] * p[0];
    case CH_TAYLOR: {
      double poly = 0.0, power = 1.0;
      for (int i = 0; i < r->coefficient_count; ++i, power *= t) poly += p[5+i] * power;
      return poly * p[0] * exp(p[2]/t) * pow(t/p[3],p[1]) * (1.0+p[4]*cond->pressure);
    }
    case CH_REVERSIBLE: return p[0]*exp(p[1]/t)*p[2];
    case CH_LAMBDA: return r->function ? r->function(cond,custom,r->context) : NAN;
    default: return NAN;
  }
}
