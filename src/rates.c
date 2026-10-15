#include "chem.h"
#include <math.h>
#include <string.h>

const char *ck_status_string(ck_status status)
{
  static const char *const names[] = {"Converged", "Invalid argument", "Allocation failed", "Unknown name",
    "Convergence exceeded maximum steps", "Step size too small", "NaN detected", "Inf detected",
    "Constraint Initialization Failed", "Singular matrix", "Accepting unconverged integration"};
  return (unsigned)status < sizeof names / sizeof *names ? names[status] : "Unknown status";
}
void ck_ideal_air_density(ck_conditions *c) { c->air_density = c->pressure / (CK_GAS_CONSTANT * c->temperature); }
ck_rate ck_rate_defaults(ck_rate_type type)
{
  ck_rate r = {0};
  r.type = type; r.A = 1; r.D = 300; r.scaling = 1;
  r.k0_A = r.kinf_A = 1; r.k0_D = r.kinf_D = r.kint_D = 300;
  r.Fc = 0.6; r.N = 1; r.probability = 1;
  if (type == CK_TERNARY_JPL19) r.k0_D = r.kinf_D = r.kint_D = 298;
  if (type == CK_TAYLOR) { r.ncoefficients = 1; r.coefficients[0] = 1; }
  return r;
}
ck_status ck_rate_prepare(ck_rate *r)
{
  if (!r || r->ncoefficients > 16 || (unsigned)r->type > CK_LAMBDA) return CK_INVALID;
  if (r->type == CK_SURFACE && (!(r->diffusion > 0) || !(r->molecular_weight > 0))) return CK_INVALID;
  if (r->type == CK_BRANCH_ALKOXY || r->type == CK_BRANCH_NITRATE) {
    r->k0 = 2.0e-22 * CK_AVOGADRO * 1.0e-6 * exp(r->n);
    double air_ref = 2.45e19 / CK_AVOGADRO * 1.0e6;
    double a = r->k0 * air_ref;
    double b = 0.43 * pow(293.0 / 298.0, -8.0);
    double val = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
    r->z = val * (1.0 - r->a0) / r->a0;
  }
  return CK_OK;
}
static double arrhenius(const ck_rate *p, double t, double pressure)
{ return p->A * exp(p->C / t) * pow(t / p->D, p->B) * (1.0 + p->E * pressure); }
double ck_rate_calculate(const ck_rate *p, const ck_conditions *c, const double *custom)
{
  double t = c->temperature, air = c->air_density;
  switch (p->type) {
    case CK_ARRHENIUS: return arrhenius(p, t, c->pressure);
    case CK_TROE:
    case CK_TERNARY: {
      double k0 = p->k0_A * exp(p->k0_C / t) * pow(t / 300.0, p->k0_B);
      double kinf = p->kinf_A * exp(p->kinf_C / t) * pow(t / 300.0, p->kinf_B);
      double ratio = k0 * air / kinf;
      return k0 * (p->type == CK_TROE ? air : 1.0) / (1.0 + ratio) *
        pow(p->Fc, p->N / (p->N + pow(log10(ratio), 2.0)));
    }
    case CK_TERNARY_JPL19: {
      double k0 = p->k0_A * exp(p->k0_C / t) * pow(t / p->k0_D, p->k0_B);
      double kinf = p->kinf_A * exp(p->kinf_C / t) * pow(t / p->kinf_D, p->kinf_B);
      double kint = p->kint_A * exp(p->kint_C / t) * pow(t / p->kint_D, p->kint_B);
      double k0_m = k0 * air, ratio = k0_m / kinf;
      double kf = k0_m / (1.0 + ratio) * pow(p->Fc, p->N / (p->N + pow(log10(ratio), 2.0)));
      return kf + kint * (1.0 - kf / kinf);
    }
    case CK_TUNNELING: return p->A * exp(-p->B / t + p->C / (t * t * t));
    case CK_BRANCH_ALKOXY:
    case CK_BRANCH_NITRATE: {
      double a = p->k0 * air, b = 0.43 * pow(t / 298.0, -8.0);
      double val = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
      double pre = p->X * exp(-p->Y / t);
      return p->type == CK_BRANCH_ALKOXY ? pre * (p->z / (p->z + val)) : pre * (val / (val + p->z));
    }
    case CK_TAYLOR: {
      double poly = 0, power = 1;
      for (size_t j = 0; j < p->ncoefficients; ++j, power *= t) poly += p->coefficients[j] * power;
      return poly * p->A * exp(p->C / t) * pow(t / p->D, p->B) * (1.0 + p->E * c->pressure);
    }
    case CK_REVERSIBLE: return p->A * exp(p->C / t) * p->reverse;
    case CK_USER: return custom[0] * p->scaling;
    case CK_SURFACE: {
      double speed_factor = 8.0 * CK_GAS_CONSTANT / (CK_PI * p->molecular_weight);
      double speed = sqrt(speed_factor * t), radius = custom[0], number = custom[1];
      return 4.0 * number * CK_PI * radius * radius / (radius / p->diffusion + 4.0 / (speed * p->probability));
    }
    case CK_LAMBDA: return p->function ? p->function(c, custom, p->context) : NAN;
  }
  return NAN;
}
