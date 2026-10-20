#include "chemistry.h"

#include <math.h>
#include <string.h>

static int chem_density_parameterized(const chem_species *species) {
  return species != NULL && (species->third_body || species->density_parameterized);
}

static double chem_density_factor(const chem_species *species) {
  return species->third_body ? 1.0 : species->density_factor;
}

int chem_species_is_parameterized(const chem_species *species) {
  return species != NULL && (species->third_body || species->density_parameterized ||
      species->temperature_parameterization_factor != 0.0);
}

double chem_species_parameterize(const chem_species *species,
                                 double air_density) {
  if (species == NULL || !isfinite(air_density) || air_density < 0.0)
    return 0.0;
  return air_density * chem_density_factor(species);
}

double chem_species_parameterize_conditions(const chem_species *species,
                                           double temperature,
                                           double air_density) {
  if (species == NULL || !isfinite(temperature) || !isfinite(air_density) ||
      temperature < 0.0 || air_density < 0.0) return 0.0;
  return species->temperature_parameterization_factor*temperature +
      ((species->third_body || species->density_parameterized)
          ? chem_species_parameterize(species,air_density) : 0.0);
}

int chem_species_set_third_body(chem_species *species) {
  if (species == NULL || species->name == NULL) return 1;
  species->third_body = 1;
  return 0;
}

static chem_property_status find_species_property(const chem_species *species,
                                                  const char *name,
                                                  const chem_property **property) {
  if (species == NULL || name == NULL || property == NULL ||
      (species->property_count && species->properties == NULL))
    return CHEM_PROPERTY_STATUS_INVALID;
  for (size_t i = 0; i < species->property_count; ++i) {
    const chem_property *candidate = &species->properties[i];
    if (candidate->name && strcmp(candidate->name,name) == 0) {
      *property = candidate;
      return CHEM_PROPERTY_STATUS_OK;
    }
  }
  return CHEM_PROPERTY_STATUS_NOT_FOUND;
}

#define DEFINE_PROPERTY_GETTER(function_name, ctype, expected_kind, member) \
chem_property_status function_name(const chem_species *species, const char *name, \
                                   ctype *value) { \
  if (value == NULL) return CHEM_PROPERTY_STATUS_INVALID; \
  const chem_property *property = NULL; \
  const chem_property_status status = find_species_property(species,name,&property); \
  if (status != CHEM_PROPERTY_STATUS_OK) return status; \
  if (property->kind != expected_kind) return CHEM_PROPERTY_STATUS_TYPE_MISMATCH; \
  *value = property->value.member; \
  return CHEM_PROPERTY_STATUS_OK; \
}

DEFINE_PROPERTY_GETTER(chem_species_get_double_property,double,
                       CHEM_PROPERTY_DOUBLE,real_value)
DEFINE_PROPERTY_GETTER(chem_species_get_integer_property,int,
                       CHEM_PROPERTY_INTEGER,integer_value)
DEFINE_PROPERTY_GETTER(chem_species_get_boolean_property,int,
                       CHEM_PROPERTY_BOOLEAN,boolean_value)
DEFINE_PROPERTY_GETTER(chem_species_get_string_property,const char *,
                       CHEM_PROPERTY_STRING,string_value)

#undef DEFINE_PROPERTY_GETTER

int chem_mechanism_absolute_tolerances(const chem_mechanism *mechanism,
                                       double *values, size_t value_count) {
  if (!mechanism || !mechanism->species || !values ||
      value_count != mechanism->species_count) return 1;
  for (size_t i = 0; i < mechanism->species_count; ++i) {
    values[i] = 1.0e-3;
    double property_value = 0.0;
    const chem_property_status status = chem_species_get_double_property(
        &mechanism->species[i], "absolute tolerance", &property_value);
    if (status == CHEM_PROPERTY_STATUS_NOT_FOUND) continue;
    if (status != CHEM_PROPERTY_STATUS_OK || !(property_value > 0.0) ||
        !isfinite(property_value)) return 1;
    values[i] = property_value;
  }
  return 0;
}

int chem_mechanism_validate_species_usage(const chem_mechanism *mechanism,
                                          unsigned char *used_species,
                                          size_t used_species_count) {
  if (!mechanism || !mechanism->species || !mechanism->species_count ||
      used_species_count != mechanism->species_count || !used_species ||
      (mechanism->reaction_count && !mechanism->reactions)) return 1;
  memset(used_species,0,used_species_count*sizeof(*used_species));
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    const chem_reaction *reaction = &mechanism->reactions[r];
    if ((reaction->reactant_count && !reaction->reactants) ||
        (reaction->product_count && !reaction->products)) return 1;
    for (size_t i = 0; i < reaction->reactant_count; ++i) {
      const size_t species = reaction->reactants[i].species;
      if (species >= mechanism->species_count) return 1;
      used_species[species] = 1;
    }
    for (size_t i = 0; i < reaction->product_count; ++i) {
      const size_t species = reaction->products[i].species;
      if (species >= mechanism->species_count) return 1;
      used_species[species] = 1;
    }
  }
  for (size_t i = 0; i < used_species_count; ++i)
    if (!used_species[i]) return 2;
  return 0;
}

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
    double k0 = reaction->p[0] * exp(reaction->p[2] / temperature) *
                pow(temperature / 300.0, reaction->p[1]);
    double kinf = reaction->p[3] * exp(reaction->p[5] / temperature) *
                  pow(temperature / 300.0, reaction->p[4]);
    double pr = k0 * air_density / kinf;
    double log_pr = log10(pr);
    return k0 * air_density / (1.0 + pr) *
           pow(reaction->p[6], reaction->p[7] /
                                      (reaction->p[7] + log_pr * log_pr));
  }
  case CHEM_TERNARY_ACTIVATION: {
    double k0 = reaction->p[0] * exp(reaction->p[2] / temperature) *
                pow(temperature / 300.0, reaction->p[1]);
    double kinf = reaction->p[3] * exp(reaction->p[5] / temperature) *
                  pow(temperature / 300.0, reaction->p[4]);
    double ratio = k0 * air_density / kinf;
    double lratio = log10(ratio);
    return k0 / (1.0 + ratio) *
           pow(reaction->p[6], reaction->p[7] /
                                      (reaction->p[7] + lratio * lratio));
  }
  case CHEM_TUNNELING:
    return reaction->p[0] * exp(-reaction->p[1] / temperature +
                                 reaction->p[2] / (temperature * temperature * temperature));
  case CHEM_TERNARY_ACTIVATION_JPL19: {
    double k0 = reaction->p[0] * exp(reaction->p[2] / temperature) *
                pow(temperature / reaction->p[3], reaction->p[1]);
    double kinf = reaction->p[4] * exp(reaction->p[6] / temperature) *
                  pow(temperature / reaction->p[7], reaction->p[5]);
    double kint = reaction->p[8] * exp(reaction->p[10] / temperature) *
                  pow(temperature / reaction->p[11], reaction->p[9]);
    double k0_m = k0 * air_density;
    double ratio = k0_m / kinf;
    double kf = k0_m / (1.0 + ratio) *
        pow(reaction->p[12], reaction->p[13] /
            (reaction->p[13] + pow(log10(ratio), 2.0)));
    return kf + kint * (1.0 - kf / kinf);
  }
  case CHEM_TAYLOR_SERIES: {
    double polynomial = 0.0, power = 1.0;
    for (size_t i = 0; i < (size_t)reaction->p[5] && i < 16; ++i) {
      polynomial += reaction->p[6+i] * power;
      power *= temperature;
    }
    double pressure = air_density * 8.31446261815324 * temperature;
    return polynomial * reaction->p[0] * exp(reaction->p[2] / temperature) *
           pow(temperature / reaction->p[3], reaction->p[1]) *
           (1.0 + reaction->p[4] * pressure);
  }
  case CHEM_REVERSIBLE:
    return reaction->p[0] * exp(reaction->p[1] / temperature) * reaction->p[2];
  case CHEM_LAMBDA:
    return reaction->rate_function
        ? reaction->rate_function(temperature,
              air_density * 8.31446261815324 * temperature, air_density,
              rate_parameters)
        : 0.0;
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
  const double mean_speed_factor = 8.0 * 8.31446261815324 /
      (3.14159265358979323846 * sp->molecular_weight);
  const double mean_speed = sqrt(mean_speed_factor * temperature);
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
      if (reaction->surface_species >= mechanism->species_count ||
          mechanism->species[reaction->surface_species].molecular_weight <= 0.0 ||
          mechanism->species[reaction->surface_species].diffusion <= 0.0)
        return 1;
      rate_constants[r] = chem_surface_rate(reaction, mechanism->species,
                                            temperature, rate_parameters);
    } else if (reaction->kind == CHEM_ARRHENIUS) {
      rate_constants[r] = reaction->p[0] * exp(reaction->p[2] / temperature) *
          pow(temperature / reaction->p[3], reaction->p[1]) *
          (1.0 + reaction->p[4] * pressure);
    } else if (reaction->kind == CHEM_TAYLOR_SERIES) {
      double polynomial = 0.0, power = 1.0;
      for (size_t i = 0; i < (size_t)reaction->p[5] && i < 16; ++i) {
        polynomial += reaction->p[6+i] * power;
        power *= temperature;
      }
      rate_constants[r] = polynomial * reaction->p[0] *
          exp(reaction->p[2] / temperature) *
          pow(temperature / reaction->p[3], reaction->p[1]) *
          (1.0 + reaction->p[4] * pressure);
    } else if (reaction->kind == CHEM_LAMBDA) {
      rate_constants[r] = reaction->rate_function
          ? reaction->rate_function(temperature, pressure, air_density,
                                    rate_parameters)
          : 0.0;
    } else {
      rate_constants[r] = chem_reaction_rate(reaction, temperature,
                                             air_density, rate_parameters);
    }
    /* Parameterized reactants are folded into MICM's rate constant before
       ordinary reactant concentrations are multiplied into the rate. */
    for (size_t i = 0; i < reaction->reactant_count; ++i) {
      const chem_term *term = &reaction->reactants[i];
      if (chem_density_parameterized(&mechanism->species[term->species]))
        rate_constants[r] *= pow(air_density *
            chem_density_factor(&mechanism->species[term->species]),
            term->coefficient);
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
      if (chem_density_parameterized(&mechanism->species[term->species]) && !constants)
        rate *= pow(air_density * chem_density_factor(&mechanism->species[term->species]),
                    term->coefficient);
      else if (!chem_density_parameterized(&mechanism->species[term->species]))
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
      if (chem_density_parameterized(&mechanism->species[term->species]) && !constants)
        rate *= pow(air_density * chem_density_factor(&mechanism->species[term->species]),
                    term->coefficient);
      else if (!chem_density_parameterized(&mechanism->species[term->species]))
        rate *= pow(concentration[term->species], term->coefficient);
    }
    for (size_t q = 0; q < reaction->reactant_count; ++q) {
      const chem_term *variable = &reaction->reactants[q];
      const size_t col = variable->species;
      if (chem_density_parameterized(&mechanism->species[col])) continue;
      double derivative;
      if (concentration[col] != 0.0) {
        derivative = rate * variable->coefficient / concentration[col];
      } else if (variable->coefficient == 1.0) {
        derivative = rate_constant;
        for (size_t i = 0; i < reaction->reactant_count; ++i) {
          if (i == q) continue;
          const chem_term *term = &reaction->reactants[i];
          if (chem_density_parameterized(&mechanism->species[term->species]) && !constants)
            derivative *= pow(air_density * chem_density_factor(&mechanism->species[term->species]),
                              term->coefficient);
          else if (!chem_density_parameterized(&mechanism->species[term->species]))
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
    if (chem_density_parameterized(&mechanism->species[col])) continue;
    for (size_t r = 0; r < mechanism->reaction_count; ++r) {
      const chem_reaction *reaction = &mechanism->reactions[r];
      for (size_t q = 0; q < reaction->reactant_count; ++q) {
        const chem_term *independent = &reaction->reactants[q];
        if (independent->species != col ||
            chem_density_parameterized(&mechanism->species[col])) continue;
        double derivative = rate_constants[r];
        for (size_t i = 0; i < reaction->reactant_count; ++i) {
          if (i == q) continue;
          const chem_term *term = &reaction->reactants[i];
          if (!chem_density_parameterized(&mechanism->species[term->species]))
            derivative *= pow(concentration[term->species], term->coefficient);
        }
        derivative *= independent->coefficient;
        if (independent->coefficient != 1.0)
          derivative *= pow(concentration[col], independent->coefficient - 1.0);
        for (size_t i = 0; i < reaction->reactant_count; ++i) {
          if (i == q) continue;
          const chem_term *term = &reaction->reactants[i];
          if (!mechanism->species[term->species].third_body)
            negative_jacobian[term->species*n + col] +=
                term->coefficient * derivative;
        }
        negative_jacobian[col*n + col] +=
            independent->coefficient * derivative;
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
