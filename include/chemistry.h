#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>

typedef enum {
  CHEM_ARRHENIUS,
  CHEM_PHOTOLYSIS,
  CHEM_USER_DEFINED,
  CHEM_TROE,
  CHEM_SURFACE,
  CHEM_TERNARY_ACTIVATION,
  CHEM_TUNNELING,
  CHEM_BRANCHED_ALKOXY,
  CHEM_BRANCHED_NITRATE
} chem_rate_kind;

typedef struct {
  size_t species;
  double coefficient;
} chem_term;

typedef struct {
  chem_rate_kind kind;
  const char *label;
  const chem_term *reactants;
  size_t reactant_count;
  const chem_term *products;
  size_t product_count;
  size_t surface_species;
  double p[8];
} chem_reaction;

typedef struct {
  const char *name;
  int third_body;
  double molecular_weight;
  double diffusion;
} chem_species;

typedef struct {
  const chem_species *species;
  size_t species_count;
  const chem_reaction *reactions;
  size_t reaction_count;
} chem_mechanism;

double chem_arrhenius(double a, double b, double c, double d, double e,
                      double temperature, double air_density);
double chem_reaction_rate(const chem_reaction *reaction, double temperature,
                          double air_density, const double *rate_parameters);
int chem_forcing(const chem_mechanism *mechanism, const double *concentration,
                 double temperature, double air_density,
                 const double *rate_parameters, double *forcing);

#endif
