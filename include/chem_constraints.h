#ifndef CHEM_CONSTRAINTS_H
#define CHEM_CONSTRAINTS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  CHEM_CONSTRAINT_LINEAR = 0,
  CHEM_CONSTRAINT_EQUILIBRIUM = 1
} chem_constraint_kind;

typedef struct {
  unsigned species;
  double coefficient;
} chem_constraint_term;

typedef struct {
  const char *name;
  chem_constraint_kind kind;
  unsigned algebraic_species;
  const chem_constraint_term *terms;
  size_t term_count;
  double constant;
  const chem_constraint_term *reactants;
  size_t reactant_count;
  const chem_constraint_term *products;
  size_t product_count;
  double equilibrium_constant_ref;
  double delta_enthalpy;
  double gas_constant;
  double reference_temperature;
  size_t parameter_index;
} chem_constraint;

typedef struct {
  size_t species_count;
  size_t constraint_count;
  size_t state_parameter_count;
  const chem_constraint *constraints;
} chem_constraint_set;

typedef enum {
  CHEM_CONSTRAINT_INIT_CONVERGED = 0,
  CHEM_CONSTRAINT_INIT_FAILED = 1,
  CHEM_CONSTRAINT_INIT_NAN = 2,
  CHEM_CONSTRAINT_INIT_INFINITY = 3,
  CHEM_CONSTRAINT_INIT_INVALID_ARGUMENT = 4
} chem_constraint_init_status;

/* Resolve a species name against the mechanism's integrated variable order.
   Return 0 on success and -1 for an unknown name or invalid arguments. */
int chem_constraint_find_species(
    const char *name,
    const char *const *species_names,
    size_t species_count,
    unsigned *species_index);

/* Validate one borrowed constraint descriptor or a borrowed constraint set.
   Return 0 on success and -1 for invalid terms, indices, parameters, or shape. */
int chem_constraint_validate(const chem_constraint *constraint, size_t species_count);
int chem_constraint_set_validate(const chem_constraint_set *set);
size_t chem_constraint_set_size(const chem_constraint_set *set);
size_t chem_constraint_set_parameter_count(const chem_constraint_set *set);

/* Fill a byte mask (one means differential) and dense row-major sparsity
   pattern for the algebraic rows. The pattern always includes each algebraic
   diagonal and every species dependency. */
int chem_constraint_set_species_mask(const chem_constraint_set *set, unsigned char *differential_species);
int chem_constraint_set_nonzero_pattern(const chem_constraint_set *set, unsigned char *pattern);

/* For equilibrium constraints, compute K(T) using the Van't Hoff equation.
   Linear constraints have no temperature dependent parameter. */
int chem_constraint_equilibrium_parameter(const chem_constraint *constraint, double temperature, double *value);

/* Overwrite algebraic rows with constraint residuals. State and forcing use
   [cell][species], state parameters use [cell][state_parameter]. */
int chem_constraint_set_apply_forcing(
    const chem_constraint_set *set,
    size_t cell_count,
    const double *state,
    const double *state_parameters,
    double *forcing);

/* Add -dG/dy to a dense [cell][row][column] matrix. */
int chem_constraint_set_accumulate_negative_jacobian(
    const chem_constraint_set *set,
    size_t cell_count,
    const double *state,
    const double *state_parameters,
    double *negative_jacobian);

/* Project algebraic species onto the constraint manifold with damped Newton
   updates. state is [cell][species]; absolute_tolerances is [species]. On
   failure, state is restored. iterations counts residual evaluations. */
chem_constraint_init_status chem_constraint_set_initialize(
    const chem_constraint_set *set,
    size_t cell_count,
    const double *state_parameters,
    double *state,
    const double *absolute_tolerances,
    double relative_tolerance,
    size_t max_updates,
    double correction_tolerance,
    size_t max_backtracks,
    double backtrack_factor,
    double sufficient_decrease,
    size_t *iterations);

const char *chem_constraint_init_status_name(chem_constraint_init_status status);

#ifdef __cplusplus
}
#endif

#endif
