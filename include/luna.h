#ifndef LUNA_H
#define LUNA_H

#include <stddef.h>

typedef enum {
  LUNA_ARRHENIUS,
  LUNA_PHOTOLYSIS,
  LUNA_USER_DEFINED,
  LUNA_TROE,
  LUNA_SURFACE
} luna_rate_kind;

typedef struct {
  size_t species;
  double coefficient;
} luna_term;

typedef struct {
  luna_rate_kind kind;
  const char *label;
  const luna_term *reactants;
  size_t reactant_count;
  const luna_term *products;
  size_t product_count;
  double p[8];
} luna_reaction;

typedef struct {
  const char *name;
  int third_body;
  double molecular_weight;
  double diffusion;
} luna_species;

typedef struct {
  const luna_species *species;
  size_t species_count;
  const luna_reaction *reactions;
  size_t reaction_count;
} luna_mechanism;

double luna_arrhenius(double a, double b, double c, double d, double e,
                      double temperature, double air_density);
double luna_reaction_rate(const luna_reaction *reaction, double temperature,
                          double air_density, const double *rate_parameters);
int luna_forcing(const luna_mechanism *mechanism, const double *concentration,
                 double temperature, double air_density,
                 const double *rate_parameters, double *forcing);

#endif
