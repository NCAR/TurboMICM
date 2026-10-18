/* Rate constants: constructors and MICM's formulas (rate_constant_functions.hpp). */
#include "internal.h"

#include <math.h>
#include <stdio.h>

#define TMC_PI 3.14159265358979323846
#define TMC_GAS_CONSTANT 8.31446261815324
#define TMC_AVOGADRO 6.02214076e23

tmc_arrhenius_parameters tmc_arrhenius_defaults(void)
{
  tmc_arrhenius_parameters p = { 1.0, 0.0, 0.0, 300.0, 0.0 };
  return p;
}

tmc_falloff_parameters tmc_falloff_defaults(void)
{
  tmc_falloff_parameters p = { 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.6, 1.0 };
  return p;
}

tmc_ternary_jpl19_parameters tmc_ternary_jpl19_defaults(void)
{
  tmc_ternary_jpl19_parameters p = { 1.0, 0.0, 0.0, 298.0, 1.0, 0.0, 0.0, 298.0, 1.0, 0.0, 0.0, 298.0, 0.6, 1.0 };
  return p;
}

tmc_tunneling_parameters tmc_tunneling_defaults(void)
{
  tmc_tunneling_parameters p = { 1.0, 0.0, 0.0 };
  return p;
}

tmc_taylor_series_parameters tmc_taylor_series_defaults(void)
{
  tmc_taylor_series_parameters p;
  memset(&p, 0, sizeof p);
  p.A = 1.0;
  p.D = 300.0;
  p.n_coefficients = 1;
  p.coefficients[0] = 1.0;
  return p;
}

tmc_reversible_parameters tmc_reversible_defaults(void)
{
  tmc_reversible_parameters p = { 1.0, 0.0, 0.0 };
  return p;
}

static tmc_rate_constant blank(tmc_rate_type type)
{
  tmc_rate_constant rc;
  memset(&rc, 0, sizeof rc);
  rc.type = type;
  return rc;
}

tmc_rate_constant tmc_rate_arrhenius(tmc_arrhenius_parameters p)
{
  tmc_rate_constant rc = blank(TMC_RATE_ARRHENIUS);
  rc.p.arrhenius = p;
  return rc;
}

tmc_rate_constant tmc_rate_troe(tmc_falloff_parameters p)
{
  tmc_rate_constant rc = blank(TMC_RATE_TROE);
  rc.p.troe = p;
  return rc;
}

tmc_rate_constant tmc_rate_ternary_chemical_activation(tmc_falloff_parameters p)
{
  tmc_rate_constant rc = blank(TMC_RATE_TERNARY_CHEMICAL_ACTIVATION);
  rc.p.troe = p;
  return rc;
}

tmc_rate_constant tmc_rate_ternary_jpl19(tmc_ternary_jpl19_parameters p)
{
  tmc_rate_constant rc = blank(TMC_RATE_TERNARY_JPL19);
  rc.p.ternary_jpl19 = p;
  return rc;
}

tmc_rate_constant tmc_rate_branched(tmc_branched_parameters p)
{
  tmc_rate_constant rc = blank(TMC_RATE_BRANCHED);
  rc.p.branched = p;
  return rc;
}

tmc_rate_constant tmc_rate_tunneling(tmc_tunneling_parameters p)
{
  tmc_rate_constant rc = blank(TMC_RATE_TUNNELING);
  rc.p.tunneling = p;
  return rc;
}

tmc_rate_constant tmc_rate_taylor_series(tmc_taylor_series_parameters p)
{
  tmc_rate_constant rc = blank(TMC_RATE_TAYLOR_SERIES);
  rc.p.taylor = p;
  return rc;
}

tmc_rate_constant tmc_rate_reversible(tmc_reversible_parameters p)
{
  tmc_rate_constant rc = blank(TMC_RATE_REVERSIBLE);
  rc.p.reversible = p;
  return rc;
}

tmc_rate_constant tmc_rate_user_defined(const char *label, double scaling_factor)
{
  tmc_rate_constant rc = blank(TMC_RATE_USER_DEFINED);
  tmc_strlcpy(rc.p.user_defined.label, label, TMC_NAME_MAX);
  rc.p.user_defined.scaling_factor = scaling_factor;
  return rc;
}

tmc_rate_constant tmc_rate_surface(const char *label, const tmc_phase_species *species, double reaction_probability)
{
  tmc_rate_constant rc = blank(TMC_RATE_SURFACE);
  tmc_strlcpy(rc.p.surface.label, label, TMC_NAME_MAX);
  rc.p.surface.species = *species;
  rc.p.surface.reaction_probability = reaction_probability;
  return rc;
}

tmc_rate_constant tmc_rate_lambda(const char *label, tmc_lambda_function function, void *context)
{
  tmc_rate_constant rc = blank(TMC_RATE_LAMBDA);
  tmc_strlcpy(rc.p.lambda.label, label, TMC_NAME_MAX);
  rc.p.lambda.function = function;
  rc.p.lambda.context = context;
  return rc;
}

int tmc_rate_constant_custom_parameter_count(const tmc_rate_constant *rc)
{
  switch (rc->type)
  {
    case TMC_RATE_USER_DEFINED: return 1;
    case TMC_RATE_SURFACE: return 2;
    default: return 0;
  }
}

int tmc_rate_constant_custom_parameter_labels(const tmc_rate_constant *rc, char (*labels)[TMC_NAME_MAX])
{
  switch (rc->type)
  {
    case TMC_RATE_USER_DEFINED: tmc_strlcpy(labels[0], rc->p.user_defined.label, TMC_NAME_MAX); return 1;
    case TMC_RATE_SURFACE:
    {
      char buffer[2 * TMC_NAME_MAX];
      snprintf(buffer, sizeof buffer, "%s.effective radius [m]", rc->p.surface.label);
      tmc_strlcpy(labels[0], buffer, TMC_NAME_MAX);
      snprintf(buffer, sizeof buffer, "%s.particle number concentration [# m-3]", rc->p.surface.label);
      tmc_strlcpy(labels[1], buffer, TMC_NAME_MAX);
      return 2;
    }
    default: return 0;
  }
}

/* Falloff kernel shared by Troe (numerator_scale = air density) and ternary
 * chemical activation (numerator_scale = 1). */
static double falloff(const tmc_falloff_parameters *p, double temperature, double air_density, double numerator_scale)
{
  double k0 = p->k0_A * exp(p->k0_C / temperature) * pow(temperature / 300.0, p->k0_B);
  double kinf = p->kinf_A * exp(p->kinf_C / temperature) * pow(temperature / 300.0, p->kinf_B);
  double ratio = k0 * air_density / kinf;
  return k0 * numerator_scale / (1.0 + ratio) * pow(p->Fc, p->N / (p->N + pow(log10(ratio), 2.0)));
}

static double arrhenius(const tmc_arrhenius_parameters *p, double temperature, double pressure)
{
  return p->A * exp(p->C / temperature) * pow(temperature / p->D, p->B) * (1.0 + p->E * pressure);
}

static double ternary_jpl19(const tmc_ternary_jpl19_parameters *p, double temperature, double air_density)
{
  double k0 = p->k0_A * exp(p->k0_C / temperature) * pow(temperature / p->k0_D, p->k0_B);
  double kinf = p->kinf_A * exp(p->kinf_C / temperature) * pow(temperature / p->kinf_D, p->kinf_B);
  double kint = p->kint_A * exp(p->kint_C / temperature) * pow(temperature / p->kint_D, p->kint_B);
  double k0_M = k0 * air_density;
  double ratio = k0_M / kinf;
  double k_f = k0_M / (1.0 + ratio) * pow(p->Fc, p->N / (p->N + pow(log10(ratio), 2.0)));
  return k_f + kint * (1.0 - k_f / kinf);
}

static double tunneling(const tmc_tunneling_parameters *p, double temperature)
{
  return p->A * exp(-p->B / temperature + p->C / (temperature * temperature * temperature));
}

/* Branched: MICM derives k0 and z from a0 and n when the solver is built. */
static double branched(const tmc_branched_parameters *p, double temperature, double air_density)
{
  double k0 = 2.0e-22 * TMC_AVOGADRO * 1.0e-6 * exp((double)p->n);
  double air_ref = 2.45e19 / TMC_AVOGADRO * 1.0e6;
  double a_ref = k0 * air_ref;
  double b_ref = 0.43 * pow(293.0 / 298.0, -8.0);
  double A_ref = a_ref / (1.0 + a_ref / b_ref) * pow(0.41, 1.0 / (1.0 + pow(log10(a_ref / b_ref), 2.0)));
  double z = A_ref * (1.0 - p->a0) / p->a0;

  double a = k0 * air_density;
  double b = 0.43 * pow(temperature / 298.0, -8.0);
  double A_val = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
  double pre = p->X * exp(-p->Y / temperature);
  return (p->branch == TMC_BRANCH_ALKOXY) ? pre * (z / (z + A_val)) : pre * (A_val / (A_val + z));
}

static double taylor_series(const tmc_taylor_series_parameters *p, double temperature, double pressure)
{
  double poly = 0.0;
  double t_pow = 1.0;
  for (int j = 0; j < p->n_coefficients; ++j, t_pow *= temperature)
    poly += p->coefficients[j] * t_pow;
  return poly * p->A * exp(p->C / temperature) * pow(temperature / p->D, p->B) * (1.0 + p->E * pressure);
}

static double surface(const tmc_surface_parameters *p, double temperature, double radius, double number_concentration)
{
  double mw = 0.0;
  tmc_species_get_double(&p->species.species, TMC_PROPERTY_MOLECULAR_WEIGHT, &mw, NULL);
  double mean_free_speed_factor = 8.0 * TMC_GAS_CONSTANT / (TMC_PI * mw);
  double mean_free_speed = sqrt(mean_free_speed_factor * temperature);
  return 4.0 * number_concentration * TMC_PI * radius * radius /
         (radius / p->species.diffusion_coefficient + 4.0 / (mean_free_speed * p->reaction_probability));
}

double tmc_rate_constant_calculate(const tmc_rate_constant *rc, const tmc_conditions *c, const double *custom)
{
  switch (rc->type)
  {
    case TMC_RATE_ARRHENIUS: return arrhenius(&rc->p.arrhenius, c->temperature, c->pressure);
    case TMC_RATE_TROE: return falloff(&rc->p.troe, c->temperature, c->air_density, c->air_density);
    case TMC_RATE_TERNARY_CHEMICAL_ACTIVATION: return falloff(&rc->p.troe, c->temperature, c->air_density, 1.0);
    case TMC_RATE_TERNARY_JPL19: return ternary_jpl19(&rc->p.ternary_jpl19, c->temperature, c->air_density);
    case TMC_RATE_BRANCHED: return branched(&rc->p.branched, c->temperature, c->air_density);
    case TMC_RATE_TUNNELING: return tunneling(&rc->p.tunneling, c->temperature);
    case TMC_RATE_TAYLOR_SERIES: return taylor_series(&rc->p.taylor, c->temperature, c->pressure);
    case TMC_RATE_REVERSIBLE: return rc->p.reversible.A * exp(rc->p.reversible.C / c->temperature) * rc->p.reversible.k_r;
    case TMC_RATE_USER_DEFINED: return custom[0] * rc->p.user_defined.scaling_factor;
    case TMC_RATE_SURFACE: return surface(&rc->p.surface, c->temperature, custom[0], custom[1]);
    case TMC_RATE_LAMBDA: return rc->p.lambda.function ? rc->p.lambda.function(c, rc->p.lambda.context) : 0.0;
    default: return 0.0;
  }
}
