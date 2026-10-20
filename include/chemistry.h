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
  CHEM_BRANCHED_NITRATE,
  CHEM_TERNARY_ACTIVATION_JPL19,
  CHEM_TAYLOR_SERIES,
  CHEM_REVERSIBLE,
  CHEM_LAMBDA
} chem_rate_kind;

typedef double (*chem_rate_function)(double temperature, double pressure,
                                     double air_density,
                                     const double *rate_parameters);

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
  double p[24];
  chem_rate_function rate_function;
} chem_reaction;

typedef struct {
  const char *name;
  int third_body;
  double molecular_weight;
  double diffusion;
  double initial_concentration;
} chem_species;

typedef struct {
  const chem_species *species;
  size_t species_count;
  const chem_reaction *reactions;
  size_t reaction_count;
  const char *const *custom_labels;
  const double *default_custom_parameters;
  size_t custom_parameter_count;
  double initial_temperature;
  double initial_pressure;
} chem_mechanism;

double chem_arrhenius(double a, double b, double c, double d, double e,
                      double temperature, double air_density);
double chem_reaction_rate(const chem_reaction *reaction, double temperature,
                          double air_density, const double *rate_parameters);
int chem_forcing(const chem_mechanism *mechanism, const double *concentration,
                 double temperature, double air_density,
                 const double *rate_parameters, double *forcing);
int chem_jacobian(const chem_mechanism *mechanism, const double *concentration,
                  double temperature, double air_density,
                  const double *rate_parameters, double *jacobian);
int chem_rate_constants(const chem_mechanism *mechanism, double temperature,
                        double pressure, double air_density,
                        const double *rate_parameters,
                        double *rate_constants);
int chem_forcing_with_constants(const chem_mechanism *mechanism,
                                const double *concentration,
                                const double *rate_constants, double air_density,
                                double *forcing);
int chem_jacobian_with_constants(const chem_mechanism *mechanism,
                                 const double *concentration,
                                 const double *rate_constants,
                                 double air_density, double *jacobian);
int chem_negative_jacobian_with_constants(
    const chem_mechanism *mechanism, const double *concentration,
    const double *rate_constants, double *negative_jacobian);

typedef struct {
  size_t internal_steps;
  size_t accepted;
  size_t rejected;
  size_t function_calls;
  size_t jacobian_updates;
  size_t decompositions;
  size_t solves;
} chem_solver_stats;

typedef struct {
  double relative_tolerance;
  double absolute_tolerance;
  size_t max_steps;
} chem_solver_options;

typedef struct {
  void *implementation;
} chem_solver_workspace;

chem_solver_options chem_solver_default_options(void);
int chem_solver_workspace_init(chem_solver_workspace *workspace,
                              const chem_mechanism *mechanism);
int chem_solver_workspace_init_group(chem_solver_workspace *workspace,
                                    const chem_mechanism *mechanism,
                                    size_t cell_capacity);
void chem_solver_workspace_destroy(chem_solver_workspace *workspace);
int chem_rosenbrock3_solve_with_workspace(
    const chem_mechanism *mechanism, chem_solver_workspace *workspace,
    double *state, double temperature, double pressure, double air_density,
    const double *rate_parameters, double duration,
    const chem_solver_options *options, chem_solver_stats *statistics);
int chem_rosenbrock3_solve_group_with_workspace(
    const chem_mechanism *mechanism, chem_solver_workspace *workspace,
    double *states, const double *temperatures, const double *pressures,
    const double *air_densities, const double *rate_parameters,
    size_t cell_count, double duration, const chem_solver_options *options,
    chem_solver_stats *statistics);
int chem_rosenbrock3_solve(const chem_mechanism *mechanism, double *state,
                           double temperature, double air_density,
                           const double *rate_parameters, double duration,
                           const chem_solver_options *options,
                           chem_solver_stats *statistics);

#endif
