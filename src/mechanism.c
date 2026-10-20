#include "luna.h"

#include <math.h>
#include <string.h>

double luna_arrhenius(double a, double b, double c, double d, double e,
                      double temperature, double air_density) {
  const double pressure = air_density * 8.31446261815324 * temperature;
  return a * exp(c / temperature) * pow(temperature / d, b) *
         (1.0 + e * pressure);
}

double luna_reaction_rate(const luna_reaction *reaction, double temperature,
                          double air_density, const double *rate_parameters) {
  (void)air_density;
  switch (reaction->kind) {
  case LUNA_ARRHENIUS:
    return luna_arrhenius(reaction->p[0], reaction->p[1], reaction->p[2],
                          reaction->p[3], reaction->p[4], temperature,
                          air_density);
  case LUNA_PHOTOLYSIS:
  case LUNA_USER_DEFINED:
  case LUNA_SURFACE:
    return rate_parameters == NULL ? 0.0 :
           rate_parameters[(size_t)reaction->p[7]] * reaction->p[0];
  case LUNA_TROE: {
    double k0 = reaction->p[0] * pow(temperature / 300.0, reaction->p[1]) *
                exp(-reaction->p[2] / temperature);
    double kinf = reaction->p[3] * pow(temperature / 300.0, reaction->p[4]) *
                  exp(-reaction->p[5] / temperature);
    double pr = air_density * k0 / kinf;
    double logfc = log10(reaction->p[6]);
    double c = -0.4 - 0.67 * logfc;
    double n = 0.75 - 1.27 * logfc;
    double x = log10(pr) + c;
    double f = x / (n - 0.14 * x);
    return kinf * pr / (1.0 + pr) * pow(10.0, logfc / (1.0 + f * f));
  }
  }
  return 0.0;
}

int luna_forcing(const luna_mechanism *mechanism, const double *concentration,
                 double temperature, double air_density,
                 const double *rate_parameters, double *forcing) {
  if (mechanism == NULL || concentration == NULL || forcing == NULL) return 1;
  memset(forcing, 0, mechanism->species_count * sizeof(*forcing));
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    const luna_reaction *reaction = &mechanism->reactions[r];
    double rate = luna_reaction_rate(reaction, temperature, air_density,
                                     rate_parameters);
    for (size_t i = 0; i < reaction->reactant_count; ++i) {
      const luna_term *term = &reaction->reactants[i];
      rate *= pow(concentration[term->species], term->coefficient);
    }
    for (size_t i = 0; i < reaction->reactant_count; ++i) {
      const luna_term *term = &reaction->reactants[i];
      forcing[term->species] -= term->coefficient * rate;
    }
    for (size_t i = 0; i < reaction->product_count; ++i) {
      const luna_term *term = &reaction->products[i];
      forcing[term->species] += term->coefficient * rate;
    }
  }
  return 0;
}
