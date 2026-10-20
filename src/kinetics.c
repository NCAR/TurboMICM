#include "kinetics.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

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
  case CHEM_RATE_KIND_COUNT:
    break;
  }
  return rate;
}

double chem_effective_rate_constant(const chem_mechanism *mechanism,
                                    const chem_reaction *reaction,
                                    double temperature, double pressure,
                                    const double *custom_parameters)
{
  const double air = pressure / (CHEM_GAS_CONSTANT * temperature);
  return chem_rate_constant(mechanism, reaction, temperature, pressure,
                            custom_parameters) *
         pow(air, (double)reaction->body_order);
}

int chem_rate_store_layout(const chem_mechanism *mechanism,
                           size_t counts[CHEM_RATE_KIND_COUNT],
                           size_t offsets[CHEM_RATE_KIND_COUNT + 1],
                           size_t *ordered_reaction_indices)
{
  size_t cursor[CHEM_RATE_KIND_COUNT];
  size_t reaction, kind;
  if (!mechanism || !counts || !offsets ||
      (mechanism->reaction_count &&
       (!mechanism->reactions || !ordered_reaction_indices)))
    return -1;
  memset(counts, 0, CHEM_RATE_KIND_COUNT * sizeof(*counts));
  for (reaction = 0; reaction < mechanism->reaction_count; ++reaction) {
    const int value = mechanism->reactions[reaction].kind;
    if (value < 0 || value >= CHEM_RATE_KIND_COUNT)
      return -1;
    ++counts[value];
  }
  offsets[0] = 0;
  for (kind = 0; kind < CHEM_RATE_KIND_COUNT; ++kind)
    offsets[kind + 1] = offsets[kind] + counts[kind];
  memcpy(cursor, offsets, CHEM_RATE_KIND_COUNT * sizeof(*cursor));
  for (reaction = 0; reaction < mechanism->reaction_count; ++reaction) {
    const size_t bucket = (size_t)mechanism->reactions[reaction].kind;
    ordered_reaction_indices[cursor[bucket]++] = reaction;
  }
  return 0;
}

static int kinetic_tables_valid(const chem_mechanism *mechanism)
{
  size_t reaction, term;
  if (!mechanism || !mechanism->species || !mechanism->species_count ||
      (mechanism->reaction_count && !mechanism->reactions) ||
      (mechanism->reactant_count && !mechanism->reactants) ||
      (mechanism->product_count && (!mechanism->products || !mechanism->yields)))
    return 0;
  for (reaction = 0; reaction < mechanism->reaction_count; ++reaction) {
    const chem_reaction *r = &mechanism->reactions[reaction];
    if (r->reactant_offset > mechanism->reactant_count ||
        r->reactant_count > mechanism->reactant_count - r->reactant_offset ||
        r->product_offset > mechanism->product_count ||
        r->product_count > mechanism->product_count - r->product_offset)
      return 0;
    for (term = 0; term < r->reactant_count; ++term)
      if (mechanism->reactants[r->reactant_offset + term] >= mechanism->species_count)
        return 0;
    for (term = 0; term < r->product_count; ++term)
      if (mechanism->products[r->product_offset + term] >= mechanism->species_count)
        return 0;
  }
  return 1;
}

int chem_mechanism_forcing(const chem_mechanism *mechanism,
                           const double *concentrations,
                           const double *rate_constants, double *forcing)
{
  size_t reaction, species;
  if (!kinetic_tables_valid(mechanism) || !forcing ||
      (mechanism->reaction_count && (!concentrations || !rate_constants)))
    return -1;
  for (species = 0; species < mechanism->species_count; ++species)
    forcing[species] = 0.0;
  for (reaction = 0; reaction < mechanism->reaction_count; ++reaction) {
    const chem_reaction *r = &mechanism->reactions[reaction];
    double flux = rate_constants[reaction];
    size_t term;
    for (term = 0; term < r->reactant_count; ++term)
      flux *= concentrations[mechanism->reactants[r->reactant_offset + term]];
    for (term = 0; term < r->reactant_count; ++term)
      forcing[mechanism->reactants[r->reactant_offset + term]] -= flux;
    for (term = 0; term < r->product_count; ++term)
      forcing[mechanism->products[r->product_offset + term]] +=
          mechanism->yields[r->product_offset + term] * flux;
  }
  return 0;
}

int chem_mechanism_jacobian_dense(const chem_mechanism *mechanism,
                                  const double *concentrations,
                                  const double *rate_constants,
                                  double *jacobian)
{
  size_t reaction, term;
  const size_t n = mechanism ? mechanism->species_count : 0;
  if (!kinetic_tables_valid(mechanism) || !jacobian ||
      n > SIZE_MAX / sizeof(*jacobian) / n ||
      (mechanism->reaction_count && (!concentrations || !rate_constants)))
    return -1;
  memset(jacobian, 0, n * n * sizeof(*jacobian));
  for (reaction = 0; reaction < mechanism->reaction_count; ++reaction) {
    const chem_reaction *r = &mechanism->reactions[reaction];
    for (term = 0; term < r->reactant_count; ++term) {
      double derivative = rate_constants[reaction];
      const size_t column = mechanism->reactants[r->reactant_offset + term];
      size_t other;
      for (other = 0; other < r->reactant_count; ++other)
        if (other != term)
          derivative *= concentrations[mechanism->reactants[r->reactant_offset + other]];
      for (other = 0; other < r->reactant_count; ++other) {
        const size_t row = mechanism->reactants[r->reactant_offset + other];
        jacobian[row * n + column] -= derivative;
      }
      for (other = 0; other < r->product_count; ++other) {
        const size_t row = mechanism->products[r->product_offset + other];
        jacobian[row * n + column] +=
            mechanism->yields[r->product_offset + other] * derivative;
      }
    }
  }
  return 0;
}
