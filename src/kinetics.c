#include "kinetics.h"

#include <math.h>
#include <stddef.h>

#define CHEM_GAS_CONSTANT (1.380649e-23 * 6.02214076e23)
#define CHEM_PI 3.141592653589793238462643383279502884

static double falloff(double k0_a, double k0_b, double k0_c,
                      double kinf_a, double kinf_b, double kinf_c,
                      double fc, double n, double temperature,
                      double air_density, double numerator_scale)
{
  const double k0 = k0_a * exp(k0_c / temperature) *
                    pow(temperature / 300.0, k0_b);
  const double kinf = kinf_a * exp(kinf_c / temperature) *
                      pow(temperature / 300.0, kinf_b);
  const double ratio = k0 * air_density / kinf;
  return k0 * numerator_scale / (1.0 + ratio) *
         pow(fc, n / (n + pow(log10(ratio), 2.0)));
}

double chem_rate_constant(const chem_mechanism *m, const chem_reaction *r,
                          double temperature, double pressure,
                          const double *custom)
{
  const double air = pressure / (CHEM_GAS_CONSTANT * temperature);
  double rate = 0.0;
  switch ((chem_rate_kind)r->kind) {
  case CHEM_ARRHENIUS:
    rate = r->p[0] * exp(r->p[2] / temperature) *
           pow(temperature / r->p[3], r->p[1]) * (1.0 + r->p[4] * pressure);
    break;
  case CHEM_TROE:
    rate = falloff(r->p[0], r->p[1], r->p[2], r->p[3], r->p[4], r->p[5],
                   r->p[6], r->p[7], temperature, air, air);
    break;
  case CHEM_TERNARY:
    rate = falloff(r->p[0], r->p[1], r->p[2], r->p[3], r->p[4], r->p[5],
                   r->p[6], r->p[7], temperature, air, 1.0);
    break;
  case CHEM_TERNARY_JPL19: {
    const double k0 = r->p[0] * exp(r->p[2] / temperature) *
                      pow(temperature / r->p[3], r->p[1]);
    const double kinf = r->p[4] * exp(r->p[6] / temperature) *
                        pow(temperature / r->p[7], r->p[5]);
    const double kint = r->p[8] * exp(r->p[10] / temperature) *
                        pow(temperature / r->p[11], r->p[9]);
    const double k0_air = k0 * air;
    const double ratio = k0_air / kinf;
    const double kf = k0_air / (1.0 + ratio) *
                      pow(r->p[12], r->p[13] /
                          (r->p[13] + pow(log10(ratio), 2.0)));
    rate = kf + kint * (1.0 - kf / kinf);
    break;
  }
  case CHEM_BRANCHED: {
    const double k0 = 2.0e-22 * 6.02214076e23 * 1.0e-6 * exp(r->p[3]);
    const double air_ref = 2.45e19 / 6.02214076e23 * 1.0e6;
    const double a0 = k0 * air_ref;
    const double b0 = 0.43 * pow(293.0 / 298.0, -8.0);
    const double aval0 = a0 / (1.0 + a0 / b0) *
                         pow(0.41, 1.0 / (1.0 + pow(log10(a0 / b0), 2.0)));
    const double z = aval0 * (1.0 - r->p[2]) / r->p[2];
    const double a = k0 * air;
    const double b = 0.43 * pow(temperature / 298.0, -8.0);
    const double aval = a / (1.0 + a / b) *
                        pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
    const double pre = r->p[1] * exp(-r->p[0] / temperature);
    rate = r->p[4] == 0.0 ? pre * z / (z + aval) : pre * aval / (aval + z);
    break;
  }
  case CHEM_TUNNELING:
    rate = r->p[0] * exp(-r->p[1] / temperature +
                         r->p[2] / (temperature * temperature * temperature));
    break;
  case CHEM_TAYLOR: {
    double polynomial = 0.0, t_power = 1.0;
    size_t j;
    for (j = 0; j < r->coefficient_count; ++j) {
      polynomial += m->coefficients[r->coefficient_offset + j] * t_power;
      t_power *= temperature;
    }
    rate = polynomial * r->p[0] * exp(r->p[2] / temperature) *
           pow(temperature / r->p[3], r->p[1]) * (1.0 + r->p[4] * pressure);
    break;
  }
  case CHEM_REVERSIBLE:
    rate = r->p[0] * exp(r->p[1] / temperature) * r->p[2];
    break;
  case CHEM_USER_DEFINED:
    rate = custom[r->custom0] * r->p[0];
    break;
  case CHEM_SURFACE: {
    const chem_species *s = &m->species[r->surface_species];
    const double radius = custom[r->custom0];
    const double number = custom[r->custom1];
    const double speed_factor = 8.0 * CHEM_GAS_CONSTANT / (CHEM_PI * s->molecular_weight);
    const double mean_free_speed = sqrt(speed_factor * temperature);
    rate = 4.0 * number * CHEM_PI * radius * radius /
           (radius / s->diffusion_coefficient +
            4.0 / (mean_free_speed * r->p[0]));
    break;
  }
  case CHEM_CALLBACK:
    if (r->callback)
      rate = r->callback(temperature, pressure, air, custom, r->callback_context);
    break;
  }
  return rate;
}
