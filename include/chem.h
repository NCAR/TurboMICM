#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  CHEM_RATE_ARRHENIUS = 0,
  CHEM_RATE_TROE,
  CHEM_RATE_TERNARY,
  CHEM_RATE_TERNARY_JPL19,
  CHEM_RATE_BRANCHED,
  CHEM_RATE_TUNNELING,
  CHEM_RATE_TAYLOR,
  CHEM_RATE_REVERSIBLE,
  CHEM_RATE_USER,
  CHEM_RATE_SURFACE
} chem_rate_kind;

enum { CHEM_RATE_PARAMETER_COUNT = 24, CHEM_SPECIES_THIRD_BODY = 1u };

typedef struct {
  const char *name;
  double molecular_weight;
  double diffusion_coefficient;
  double initial_value;
  double absolute_tolerance;
  unsigned flags;
} chem_species;

typedef struct {
  unsigned species;
  double coefficient;
} chem_term;

typedef struct {
  chem_rate_kind kind;
  const char *label;
  double p[CHEM_RATE_PARAMETER_COUNT];
  const double *coefficients;
  size_t coefficient_count;
  unsigned custom_parameter;
  unsigned surface_species;
  unsigned air_density_power;
  unsigned branch;
} chem_rate;

typedef struct {
  unsigned reactant_offset;
  unsigned reactant_count;
  unsigned product_offset;
  unsigned product_count;
  chem_rate rate;
} chem_reaction;

typedef struct {
  const char *name;
  double initial_value;
} chem_parameter;

typedef struct {
  const char *name;
  size_t species_count;
  size_t reaction_count;
  size_t reactant_term_count;
  size_t product_term_count;
  size_t parameter_count;
  const chem_species *species;
  const chem_reaction *reactions;
  const chem_term *reactant_terms;
  const chem_term *product_terms;
  const chem_parameter *parameters;
  double temperature;
  double pressure;
  double relative_tolerance;
} chem_mechanism;

typedef struct {
  unsigned long internal_steps;
  unsigned long accepted;
  unsigned long rejected;
  unsigned long function_calls;
  unsigned long jacobian_updates;
  unsigned long decompositions;
  unsigned long solves;
} chem_statistics;

const chem_mechanism *chem_generated_mechanism(void);

/* Evaluate one reaction's rate coefficient for one cell. Parameters are in
   descriptor order; pass NULL when the mechanism has no custom parameters.
   Returns 0 on success, -1 for invalid input, and -2 when a surface-rate
   species lacks its required molecular weight or diffusion coefficient. */
int chem_evaluate_rate(
    const chem_mechanism *mechanism,
    size_t reaction_index,
    const double *parameters,
    double temperature,
    double pressure,
    double *rate_constant);

/* Integrate each row of y in place for one chemistry interval.  State layout
   is cell-major: y[cell * species_count + species].  Parameters use the same
   cell-major convention.  A separate call owns a separate adaptive history. */
int chem_integrate_group(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y,
    double *restrict parameters,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    chem_statistics *statistics);

int chem_integrate_group_steps(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y,
    double *restrict parameters,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    size_t interval_count,
    chem_statistics *statistics);

#ifdef __cplusplus
}
#endif

#endif
