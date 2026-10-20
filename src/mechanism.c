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
                exp(-reaction->p[2] / temperature);
    double kinf = reaction->p[3] * pow(temperature / 300.0, reaction->p[4]) *
                  exp(-reaction->p[5] / temperature);
    double pr = air_density * k0 / kinf;
    double log_pr = log10(pr);
    return k0 * air_density / (1.0 + pr) *
           pow(reaction->p[6], reaction->p[7] /
                                      (reaction->p[7] + log_pr * log_pr));
  }
  }
  return 0.0;
}

int chem_forcing(const chem_mechanism *mechanism, const double *concentration,
                 double temperature, double air_density,
                 const double *rate_parameters, double *forcing) {
  if (mechanism == NULL || concentration == NULL || forcing == NULL) return 1;
  memset(forcing, 0, mechanism->species_count * sizeof(*forcing));
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    const chem_reaction *reaction = &mechanism->reactions[r];
    double rate = chem_reaction_rate(reaction, temperature, air_density,
                                     rate_parameters);
    for (size_t i = 0; i < reaction->reactant_count; ++i) {
      const chem_term *term = &reaction->reactants[i];
      rate *= pow(concentration[term->species], term->coefficient);
    }
    for (size_t i = 0; i < reaction->reactant_count; ++i) {
      const chem_term *term = &reaction->reactants[i];
      forcing[term->species] -= term->coefficient * rate;
    }
    for (size_t i = 0; i < reaction->product_count; ++i) {
      const chem_term *term = &reaction->products[i];
      forcing[term->species] += term->coefficient * rate;
    }
  }
  return 0;
}
