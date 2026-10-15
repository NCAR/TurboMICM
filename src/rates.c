/* Translated formulas: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
 * SPDX-License-Identifier: Apache-2.0 */
#include "chem.h"
#include <math.h>
#include <string.h>

void chem_ideal_air_density(ChemConditions *c) {
  c->air_density = c->pressure / (CHEM_GAS_CONSTANT * c->temperature);
}

ChemRate chem_rate(ChemRateType type) {
  ChemRate r;
  memset(&r, 0, sizeof r);
  r.type = type;
  r.A = r.X = 1.0;
  r.D = 300.0;
  r.k0_A = r.kinf_A = r.kint_A = 1.0;
  r.k0_D = r.kinf_D = r.kint_D = 298.0;
  r.Fc = 0.6;
  r.N = r.probability = 1.0;
  r.custom_index = SIZE_MAX;
  r.diffusion = r.molecular_weight = NAN;
  return r;
}

double chem_rate_value(const ChemRate *r, const ChemConditions *c, const double *u) {
  double t = c->temperature, density = c->air_density;
  switch (r->type) {
    case CHEM_ARRHENIUS:
      return r->A * exp(r->C / t) * pow(t / r->D, r->B) * (1.0 + r->E * c->pressure);
    case CHEM_TROE:
    case CHEM_TERNARY:
    case CHEM_JPL19: {
      double d0 = r->type == CHEM_JPL19 ? r->k0_D : 300.0;
      double di = r->type == CHEM_JPL19 ? r->kinf_D : 300.0;
      double k0 = r->k0_A * exp(r->k0_C / t) * pow(t / d0, r->k0_B);
      double ki = r->kinf_A * exp(r->kinf_C / t) * pow(t / di, r->kinf_B);
      double ratio = k0 * density / ki;
      double scale = r->type == CHEM_TERNARY ? 1.0 : density;
      double k = k0 * scale / (1.0 + ratio) * pow(r->Fc, r->N / (r->N + pow(log10(ratio), 2.0)));
      if (r->type == CHEM_JPL19) {
        double kint = r->kint_A * exp(r->kint_C / t) * pow(t / r->kint_D, r->kint_B);
        return k + kint * (1.0 - k / ki);
      }
      return k;
    }
    case CHEM_BRANCHED: {
      double k0 = 2.0e-22 * 6.02214076e23 * 1.0e-6 * exp(r->n);
      double a_ref = k0 * (2.45e19 / 6.02214076e23 * 1.0e6);
      double b_ref = 0.43 * pow(293.0 / 298.0, -8.0);
      double v_ref = a_ref / (1.0 + a_ref / b_ref) * pow(0.41, 1.0 / (1.0 + pow(log10(a_ref / b_ref), 2.0)));
      double z = v_ref * (1.0 - r->a0) / r->a0;
      double a = k0 * density, b = 0.43 * pow(t / 298.0, -8.0);
      double v = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
      double pre = r->X * exp(-r->Y / t);
      return r->branch == 0 ? pre * (z / (z + v)) : pre * (v / (v + z));
    }
    case CHEM_TUNNELING:
      return r->A * exp(-r->B / t + r->C / (t * t * t));
    case CHEM_TAYLOR: {
      double polynomial = 0.0, power = 1.0;
      for (size_t j = 0; j < r->ncoefficients; ++j, power *= t)
        polynomial += r->coefficients[j] * power;
      return polynomial * r->A * exp(r->C / t) * pow(t / r->D, r->B) * (1.0 + r->E * c->pressure);
    }
    case CHEM_REVERSIBLE:
      return r->A * exp(r->C / t) * r->B; /* B stores k_r */
    case CHEM_USER:
      return u[r->custom_index] * r->A; /* A stores the scaling factor */
    case CHEM_SURFACE: {
      double radius = u[r->custom_index], number = u[r->custom_index + 1];
      double speed_factor = 8.0 * CHEM_GAS_CONSTANT / (CHEM_PI * r->molecular_weight);
      double speed = sqrt(speed_factor * t);
      return 4.0 * number * CHEM_PI * radius * radius / (radius / r->diffusion + 4.0 / (speed * r->probability));
    }
    case CHEM_LAMBDA:
      return r->function ? r->function(c, u, r->context) : 0.0;
  }
  return NAN;
}
