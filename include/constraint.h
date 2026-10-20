#ifndef CHEM_CONSTRAINT_H
#define CHEM_CONSTRAINT_H

#include <stddef.h>

typedef enum {
  CHEM_CONSTRAINT_LINEAR = 0,
  CHEM_CONSTRAINT_EQUILIBRIUM = 1,
  CHEM_CONSTRAINT_CALLBACK = 2
} chem_constraint_kind;

typedef int (*chem_constraint_residual_callback)(
    const double *concentrations, size_t species_count,
    double equilibrium_constant, double temperature,
    double *residual, void *context);
typedef int (*chem_constraint_jacobian_callback)(
    const double *concentrations, size_t species_count,
    double equilibrium_constant, double temperature,
    double *gradient, void *context);

typedef struct {
  size_t species_index;
  double coefficient;
} chem_constraint_term;

typedef struct {
  const char *name;
  chem_constraint_kind kind;
  size_t algebraic_species;
  const chem_constraint_term *terms;
  size_t term_count;
  double constant;
  const chem_constraint_term *reactants;
  size_t reactant_count;
  const chem_constraint_term *products;
  size_t product_count;
  double equilibrium_constant_reference;
  double enthalpy_change;
  double gas_constant;
  double reference_temperature;
  const size_t *dependencies;
  size_t dependency_count;
  chem_constraint_residual_callback residual_callback;
  chem_constraint_jacobian_callback jacobian_callback;
  void *callback_context;
} chem_constraint;

typedef struct {
  const chem_constraint *constraints;
  size_t constraint_count;
  size_t species_count;
} chem_constraint_set;

typedef struct {
  size_t row;
  size_t column;
} chem_jacobian_element;

int chem_constraint_linear_init(chem_constraint *constraint, const char *name,
                                size_t algebraic_species,
                                const chem_constraint_term *terms,
                                size_t term_count, double constant);
int chem_constraint_equilibrium_init(
    chem_constraint *constraint, const char *name, size_t algebraic_species,
    const chem_constraint_term *reactants, size_t reactant_count,
    const chem_constraint_term *products, size_t product_count,
    double equilibrium_constant_reference, double enthalpy_change,
    double gas_constant, double reference_temperature);
int chem_constraint_callback_init(
    chem_constraint *constraint, const char *name, size_t algebraic_species,
    const size_t *dependencies, size_t dependency_count,
    chem_constraint_residual_callback residual_callback,
    chem_constraint_jacobian_callback jacobian_callback, void *context);
int chem_constraint_set_init(chem_constraint_set *set,
                             const chem_constraint *constraints,
                             size_t constraint_count, size_t species_count);
size_t chem_constraint_set_size(const chem_constraint_set *set);
int chem_constraint_set_algebraic_species(const chem_constraint_set *set,
                                          size_t *indices, size_t capacity);
int chem_constraint_set_nonzero_jacobian_elements(
    const chem_constraint_set *set, chem_jacobian_element *elements,
    size_t capacity, size_t *element_count);
int chem_constraint_residual(const chem_constraint *constraint,
                             const double *concentrations, size_t species_count,
                             double equilibrium_constant, double temperature,
                             double *residual);
int chem_constraint_jacobian(const chem_constraint *constraint,
                             const double *concentrations, size_t species_count,
                             double equilibrium_constant, double temperature,
                             double *gradient);
int chem_constraint_set_add_forcing(const chem_constraint_set *set,
                                    size_t cell_count,
                                    const double *concentrations,
                                    const double *equilibrium_constants,
                                    const double *temperatures,
                                    double *forcing);
int chem_constraint_set_subtract_jacobian(
    const chem_constraint_set *set, size_t cell_count,
    const double *concentrations, const double *equilibrium_constants,
    const double *temperatures, double *dense_jacobians);
int chem_constraint_set_project(
    const chem_constraint_set *set, size_t cell_count, double *concentrations,
    const double *equilibrium_constants, const double *temperatures,
    const double *absolute_tolerances, double relative_tolerance,
    double correction_tolerance, size_t maximum_iterations,
    size_t maximum_backtracks, double backtrack_factor,
    double sufficient_decrease, size_t *iterations);

#endif
