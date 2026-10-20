#include "chemistry.h"

#include <math.h>
#include <string.h>

double chem_arrhenius(double a, double b, double c, double d, double e,
                      double temperature, double air_density) {
  const double pressure = air_density * 8.31446261815324 * temperature;
  return a * exp(c / temperature) * pow(temperature / d, b) *
         (1.0 + e * pressure);
}

double chem_reaction_rate(const chem_reaction *reaction, double temperature,
                          double air_density, const double *rate_parameters) {
  if (reaction == NULL || temperature <= 0.0) return 0.0;
  (void)air_density;
  switch (reaction->kind) {
  case CHEM_ARRHENIUS:
    return chem_arrhenius(reaction->p[0], reaction->p[1], reaction->p[2],
                          reaction->p[3], reaction->p[4], temperature,
                          air_density);
  case CHEM_PHOTOLYSIS:
  case CHEM_USER_DEFINED:
  case CHEM_SURFACE:
    return rate_parameters == NULL ? 0.0 :
           rate_parameters[(size_t)reaction->p[7]] * reaction->p[0];
  case CHEM_TROE: {
    double k0 = reaction->p[0] * pow(temperature / 300.0, reaction->p[1]) *
                exp(reaction->p[2] / temperature);
    double kinf = reaction->p[3] * pow(temperature / 300.0, reaction->p[4]) *
                  exp(reaction->p[5] / temperature);
    double pr = air_density * k0 / kinf;
    double log_pr = log10(pr);
    return k0 * air_density / (1.0 + pr) *
           pow(reaction->p[6], reaction->p[7] /
                                      (reaction->p[7] + log_pr * log_pr));
  }
  case CHEM_TERNARY_ACTIVATION: {
    double k0 = reaction->p[0] * pow(temperature / 300.0, reaction->p[1]) *
                exp(reaction->p[2] / temperature);
    double kinf = reaction->p[3] * pow(temperature / 300.0, reaction->p[4]) *
                  exp(reaction->p[5] / temperature);
    double ratio = air_density * k0 / kinf;
    double lratio = log10(ratio);
    return k0 / (1.0 + ratio) *
           pow(reaction->p[6], reaction->p[7] /
                                      (reaction->p[7] + lratio * lratio));
  }
  case CHEM_TUNNELING:
    return reaction->p[0] * exp(-reaction->p[1] / temperature +
                                 reaction->p[2] / (temperature * temperature * temperature));
  case CHEM_BRANCHED_ALKOXY:
  case CHEM_BRANCHED_NITRATE: {
    const double k0 = 2.0e-22 * 6.02214076e23 * 1.0e-6 * exp(reaction->p[3]);
    const double a = k0 * air_density;
    const double b = 0.43 * pow(temperature / 298.0, -8.0);
    const double ratio = a / b;
    const double av = a / (1.0 + ratio) * pow(0.41, 1.0 / (1.0 + pow(log10(ratio), 2.0)));
    const double z = reaction->p[4];
    const double fraction = reaction->kind == CHEM_BRANCHED_ALKOXY ? z / (z + av) : av / (av + z);
    return reaction->p[0] * exp(-reaction->p[1] / temperature) * fraction;
  }
  }
  return 0.0;
}

static double chem_surface_rate(const chem_reaction *reaction,
                                const chem_species *species, double temperature,
                                const double *rate_parameters) {
  if (rate_parameters == NULL || reaction->surface_species == (size_t)-1) return 0.0;
  const chem_species *sp = &species[reaction->surface_species];
  if (sp->molecular_weight <= 0.0 || sp->diffusion <= 0.0 || reaction->p[0] <= 0.0)
    return 0.0;
  const double radius = rate_parameters[(size_t)reaction->p[7]];
  const double number = rate_parameters[(size_t)reaction->p[7] + 1];
  const double mean_speed = sqrt(8.0 * 8.31446261815324 * temperature /
                                 (3.14159265358979323846 * sp->molecular_weight));
  return 4.0 * number * 3.14159265358979323846 * radius * radius /
         (radius / sp->diffusion + 4.0 / (mean_speed * reaction->p[0]));
}

int chem_rate_constants(const chem_mechanism *mechanism, double temperature,
                        double pressure, double air_density,
                        const double *rate_parameters,
                        double *rate_constants) {
  if (!mechanism || !rate_constants ||
      (mechanism->reaction_count && !mechanism->reactions) ||
      (mechanism->species_count && !mechanism->species) ||
      !isfinite(temperature) || temperature <= 0.0 || !isfinite(pressure) ||
      !isfinite(air_density) || air_density < 0.0) return 1;
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    const chem_reaction *reaction = &mechanism->reactions[r];
    if (reaction->kind == CHEM_SURFACE) {
      rate_constants[r] = chem_surface_rate(reaction, mechanism->species,
                                            temperature, rate_parameters);
    } else if (reaction->kind == CHEM_ARRHENIUS) {
      rate_constants[r] = reaction->p[0] * exp(reaction->p[2] / temperature) *
          pow(temperature / reaction->p[3], reaction->p[1]) *
          (1.0 + reaction->p[4] * pressure);
    } else {
      rate_constants[r] = chem_reaction_rate(reaction, temperature,
                                             air_density, rate_parameters);
    }
    /* Parameterized reactants are folded into MICM's rate constant before
       ordinary reactant concentrations are multiplied into the rate. */
    for (size_t i = 0; i < reaction->reactant_count; ++i) {
      const chem_term *term = &reaction->reactants[i];
      if (mechanism->species[term->species].third_body)
        rate_constants[r] *= pow(air_density, term->coefficient);
    }
  }
  return 0;
}

static int chem_forcing_impl(const chem_mechanism *mechanism, const double *concentration,
                             double temperature, double air_density,
                             const double *rate_parameters, const double *constants,
                             double *forcing) {
  if (mechanism == NULL || concentration == NULL || forcing == NULL ||
      (mechanism->species_count != 0 && mechanism->species == NULL) ||
      (mechanism->reaction_count != 0 && mechanism->reactions == NULL)) return 1;
  memset(forcing, 0, mechanism->species_count * sizeof(*forcing));
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    const chem_reaction *reaction = &mechanism->reactions[r];
    if ((reaction->reactant_count != 0 && reaction->reactants == NULL) ||
        (reaction->product_count != 0 && reaction->products == NULL)) return 1;
    double rate_constant = constants ? constants[r] : reaction->kind == CHEM_SURFACE
                      ? chem_surface_rate(reaction, mechanism->species, temperature,
                                          rate_parameters)
                      : chem_reaction_rate(reaction, temperature, air_density,
                                           rate_parameters);
    double rate = rate_constant;
    for (size_t i = 0; i < reaction->reactant_count; ++i) {
      const chem_term *term = &reaction->reactants[i];
      if (term->species >= mechanism->species_count) return 1;
      if (mechanism->species[term->species].third_body && !constants)
        rate *= pow(air_density, term->coefficient);
      else if (!mechanism->species[term->species].third_body)
        rate *= pow(concentration[term->species], term->coefficient);
    }
    for (size_t i = 0; i < reaction->reactant_count; ++i) {
      const chem_term *term = &reaction->reactants[i];
      if (term->species >= mechanism->species_count) return 1;
      if (mechanism->species[term->species].third_body) continue;
      forcing[term->species] -= term->coefficient * rate;
    }
    for (size_t i = 0; i < reaction->product_count; ++i) {
      const chem_term *term = &reaction->products[i];
      if (term->species >= mechanism->species_count) return 1;
      if (mechanism->species[term->species].third_body) continue;
      forcing[term->species] += term->coefficient * rate;
    }
  }
  return 0;
}

int chem_forcing(const chem_mechanism *mechanism, const double *concentration,
                 double temperature, double air_density,
                 const double *rate_parameters, double *forcing) {
  return chem_forcing_impl(mechanism, concentration, temperature, air_density,
                           rate_parameters, NULL, forcing);
}

int chem_forcing_with_constants(const chem_mechanism *mechanism,
                                const double *concentration,
                                const double *rate_constants, double air_density,
                                double *forcing) {
  if (!rate_constants) return 1;
  return chem_forcing_impl(mechanism, concentration, 0.0, air_density, NULL,
                           rate_constants, forcing);
}

static int chem_jacobian_impl(const chem_mechanism *mechanism, const double *concentration,
                              double temperature, double air_density,
                              const double *rate_parameters, const double *constants,
                              double *jacobian) {
  if (mechanism == NULL || concentration == NULL || jacobian == NULL ||
      (mechanism->species_count != 0 && mechanism->species == NULL) ||
      (mechanism->reaction_count != 0 && mechanism->reactions == NULL)) return 1;
  const size_t n = mechanism->species_count;
  memset(jacobian, 0, n * n * sizeof(*jacobian));
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    const chem_reaction *reaction = &mechanism->reactions[r];
    double rate_constant = constants ? constants[r] : reaction->kind == CHEM_SURFACE
                               ? chem_surface_rate(reaction, mechanism->species,
                                                   temperature, rate_parameters)
                               : chem_reaction_rate(reaction, temperature,
                                                    air_density, rate_parameters);
    double rate = rate_constant;
    for (size_t i = 0; i < reaction->reactant_count; ++i) {
      const chem_term *term = &reaction->reactants[i];
      if (mechanism->species[term->species].third_body && !constants)
        rate *= pow(air_density, term->coefficient);
      else if (!mechanism->species[term->species].third_body)
        rate *= pow(concentration[term->species], term->coefficient);
    }
    for (size_t q = 0; q < reaction->reactant_count; ++q) {
      const chem_term *variable = &reaction->reactants[q];
      const size_t col = variable->species;
      if (mechanism->species[col].third_body) continue;
      double derivative;
      if (concentration[col] != 0.0) {
        derivative = rate * variable->coefficient / concentration[col];
      } else if (variable->coefficient == 1.0) {
        derivative = rate_constant;
        for (size_t i = 0; i < reaction->reactant_count; ++i) {
          if (i == q) continue;
          const chem_term *term = &reaction->reactants[i];
          if (mechanism->species[term->species].third_body && !constants)
            derivative *= pow(air_density, term->coefficient);
          else if (!mechanism->species[term->species].third_body)
            derivative *= pow(concentration[term->species], term->coefficient);
        }
      } else {
        derivative = 0.0;
      }
      for (size_t i = 0; i < reaction->reactant_count; ++i) {
        const chem_term *term = &reaction->reactants[i];
        if (!mechanism->species[term->species].third_body)
          jacobian[term->species * n + col] -= term->coefficient * derivative;
      }
      for (size_t i = 0; i < reaction->product_count; ++i) {
        const chem_term *term = &reaction->products[i];
        if (!mechanism->species[term->species].third_body)
          jacobian[term->species * n + col] += term->coefficient * derivative;
      }
    }
  }
  return 0;
}

int chem_jacobian(const chem_mechanism *mechanism, const double *concentration,
                  double temperature, double air_density,
                  const double *rate_parameters, double *jacobian) {
  return chem_jacobian_impl(mechanism, concentration, temperature, air_density,
                            rate_parameters, NULL, jacobian);
}

int chem_jacobian_with_constants(const chem_mechanism *mechanism,
                                 const double *concentration,
                                 const double *rate_constants,
                                 double air_density, double *jacobian) {
  if (!rate_constants) return 1;
  return chem_jacobian_impl(mechanism, concentration, 0.0, air_density, NULL,
                            rate_constants, jacobian);
}

int chem_negative_jacobian_with_constants(
    const chem_mechanism *mechanism, const double *concentration,
    const double *rate_constants, double *negative_jacobian) {
  if (!mechanism || !concentration || !rate_constants || !negative_jacobian ||
      (mechanism->species_count && !mechanism->species) ||
      (mechanism->reaction_count && !mechanism->reactions)) return 1;
  const size_t n = mechanism->species_count;
  memset(negative_jacobian, 0, n*n*sizeof(*negative_jacobian));
  /* MICM accumulates Jacobian terms by independent variable, preserving
     reaction order within each column. Match that order for solver arithmetic. */
  for (size_t col = 0; col < n; ++col) {
    if (mechanism->species[col].third_body) continue;
    for (size_t r = 0; r < mechanism->reaction_count; ++r) {
      const chem_reaction *reaction = &mechanism->reactions[r];
      for (size_t q = 0; q < reaction->reactant_count; ++q) {
        const chem_term *independent = &reaction->reactants[q];
        if (independent->species != col ||
            mechanism->species[col].third_body) continue;
        double derivative = rate_constants[r];
        int skipped_independent = 0;
        for (size_t i = 0; i < reaction->reactant_count; ++i) {
          const chem_term *term = &reaction->reactants[i];
          if (!skipped_independent && term->species == col) {
            skipped_independent = 1;
            continue;
          }
          if (!mechanism->species[term->species].third_body)
            derivative *= pow(concentration[term->species], term->coefficient);
        }
        skipped_independent = 0;
        for (size_t i = 0; i < reaction->reactant_count; ++i) {
          const chem_term *term = &reaction->reactants[i];
          if (!skipped_independent && term->species == col) {
            skipped_independent = 1;
            continue;
          }
          if (!mechanism->species[term->species].third_body)
            negative_jacobian[term->species*n + col] += derivative;
        }
        negative_jacobian[col*n + col] += derivative;
        for (size_t i = 0; i < reaction->product_count; ++i) {
          const chem_term *term = &reaction->products[i];
          if (!mechanism->species[term->species].third_body)
            negative_jacobian[term->species*n + col] -=
                term->coefficient * derivative;
        }
      }
    }
  }
  return 0;
}
