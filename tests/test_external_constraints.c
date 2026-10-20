#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
  return 1; \
} } while (0)
#define NEAR(actual, expected, tolerance) CHECK(fabs((actual) - (expected)) <= (tolerance))
#define GAS_CONSTANT 8.31446261815324
#define T_REF 298.15

static chem_constraint linear_constraint(const char *name, unsigned algebraic,
    const chem_constraint_term *terms, size_t count, double constant) {
  return (chem_constraint){ .name = name, .kind = CHEM_CONSTRAINT_LINEAR,
    .algebraic_species = algebraic, .terms = terms, .term_count = count, .constant = constant };
}

static chem_constraint equilibrium_constraint(const char *name, unsigned algebraic,
    unsigned reactant, unsigned product, double keq, double delta_h) {
  static chem_constraint_term reactants[1], products[1];
  reactants[0] = (chem_constraint_term){ reactant, 1.0 };
  products[0] = (chem_constraint_term){ product, 1.0 };
  return (chem_constraint){ .name = name, .kind = CHEM_CONSTRAINT_EQUILIBRIUM,
    .algebraic_species = algebraic, .reactants = reactants, .reactant_count = 1u,
    .products = products, .product_count = 1u, .equilibrium_constant_ref = keq,
    .delta_enthalpy = delta_h, .gas_constant = GAS_CONSTANT,
    .reference_temperature = T_REF, .parameter_index = 0u };
}

static chem_constraint_set make_set(size_t species_count, const chem_constraint *items,
    size_t count, size_t parameters) {
  return (chem_constraint_set){ species_count, count, parameters, items };
}

static chem_rosenbrock_parameters dae_parameters(void) {
  chem_rosenbrock_parameters p = chem_rosenbrock_four_stage_da_parameters();
  p.max_internal_steps = 100000u;
  return p;
}

static chem_mechanism one_step_mechanism(const chem_species *species, size_t species_count,
    unsigned reactant, unsigned product, double rate) {
  /* The caller keeps the returned mechanism's one-entry arrays alive by using
     static storage for descriptors in each test. */
  static chem_term r[1], p[1];
  static chem_reaction reactions[1];
  r[0] = (chem_term){ reactant, 1.0 };
  p[0] = (chem_term){ product, 1.0 };
  reactions[0] = (chem_reaction){ 0u, 1u, 0u, 1u,
    { .kind = CHEM_RATE_ARRHENIUS, .label = "first_order", .p = { [0] = rate, [3] = T_REF } } };
  return (chem_mechanism){ .name = "external constraint equivalent", .species_count = species_count,
    .reaction_count = 1u, .reactant_term_count = 1u, .product_term_count = 1u,
    .species = species, .reactions = reactions, .reactant_terms = r, .product_terms = p,
    .temperature = T_REF, .pressure = 101325.0, .relative_tolerance = 1e-4 };
}

static int integrate_dae(const chem_mechanism *mechanism, const chem_constraint_set *set,
    size_t cells, double *state, const double *constraint_parameters,
    const double *temperatures, double interval) {
  double pressure[3] = { 101325.0, 101325.0, 101325.0 };
  chem_rosenbrock_parameters solver = dae_parameters();
  chem_statistics statistics = { 0 };
  int result = chem_integrate_group_with_constraints(mechanism, set, cells, state, NULL,
      constraint_parameters, temperatures, pressure, interval, &solver, &statistics);
  if (result != 0) fprintf(stderr, "DAE returned %d (accepted=%lu rejected=%lu steps=%lu species=%zu constraints=%zu)\n",
      result, statistics.accepted, statistics.rejected, statistics.internal_steps,
      mechanism->species_count, set->constraint_count);
  return result;
}

static int add_external_model_with_constraints(void) {
  const chem_species species[] = {
    { "A_GAS", 0.0, 0.0, 0.0, 1e-8, 0u }, { "AEROSOL.A_AQ", 0.0, 0.0, 0.0, 1e-8, 0u }
  };
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint item = linear_constraint("total_mass", 1u, terms, 2u, 1.0);
  const chem_constraint_set set = make_set(2u, &item, 1u, 0u);
  unsigned char mask[2], pattern[4];
  CHECK(chem_constraint_set_validate(&set) == 0);
  CHECK(chem_constraint_set_size(&set) == 1u);
  CHECK(strcmp(species[0].name, "A_GAS") == 0 && strcmp(species[1].name, "AEROSOL.A_AQ") == 0);
  CHECK(chem_constraint_set_species_mask(&set, mask) == 0);
  CHECK(mask[0] == 1u && mask[1] == 0u);
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  CHECK(pattern[2] && pattern[3] && !pattern[0] && !pattern[1]);
  return 0;
}

static int add_external_model_process_only(void) {
  const chem_species species[] = {
    { "A_GAS", 0.0, 0.0, 0.0, 1e-10, 0u }, { "AEROSOL.A_AQ", 0.0, 0.0, 0.0, 1e-10, 0u }
  };
  chem_mechanism mechanism = one_step_mechanism(species, 2u, 0u, 1u, 0.1);
  const chem_constraint_set empty = make_set(2u, NULL, 0u, 0u);
  unsigned char mask[2];
  CHECK(chem_constraint_set_species_mask(&empty, mask) == 0 && mask[0] && mask[1]);
  double state[] = { 1.0, 0.0 }, temperature[] = { T_REF }, pressure[] = { 101325.0 };
  CHECK(chem_integrate_group(&mechanism, 1u, state, NULL, temperature, pressure, 1.0, NULL) == 0);
  NEAR(state[0], exp(-0.1), 2e-6);
  NEAR(state[1], 1.0 - exp(-0.1), 2e-6);
  return 0;
}

static int dae_solve_enforces_conservation(void) {
  const chem_species species[] = {
    { "A_GAS", 0.0, 0.0, 0.0, 1e-4, 0u }, { "AEROSOL.A_AQ", 0.0, 0.0, 0.0, 1e-4, 0u }
  };
  chem_mechanism mechanism = one_step_mechanism(species, 2u, 0u, 1u, 0.1);
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint item = linear_constraint("total_mass", 1u, terms, 2u, 1.0);
  const chem_constraint_set set = make_set(2u, &item, 1u, 0u);
  const double temperature[] = { T_REF };
  double state[] = { 0.9, 0.1 };
  for (size_t step = 0; step < 20u; ++step) {
    CHECK(integrate_dae(&mechanism, &set, 1u, state, NULL, temperature, 10.0) == 0);
    NEAR(state[0] + state[1], 1.0, 1e-7);
  }
  CHECK(state[0] < 0.13 && state[1] > 0.87);
  return 0;
}

static int algebraic_row_holds_only_constraint_terms(void) {
  const chem_species species[] = {
    { "A_GAS", 0.0, 0.0, 0.0, 1e-8, 0u }, { "AEROSOL.A_AQ", 0.0, 0.0, 0.0, 1e-8, 0u }
  };
  chem_mechanism mechanism = one_step_mechanism(species, 2u, 0u, 1u, 10.0);
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint item = linear_constraint("total_mass", 1u, terms, 2u, 1.0);
  const chem_constraint_set set = make_set(2u, &item, 1u, 0u);
  const double state[] = { 0.7, 0.2 }, rates[] = { 10.0 };
  double forcing[2], negative_jacobian[4] = { 0.0 };
  unsigned char mask[2];
  CHECK(chem_constraint_set_species_mask(&set, mask) == 0);
  CHECK(chem_evaluate_forcing(&mechanism, 1u, state, rates, forcing) == 0);
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, state, NULL, forcing) == 0);
  NEAR(forcing[0], -7.0, 1e-14);
  NEAR(forcing[1], -0.1, 1e-14);
  CHECK(chem_accumulate_negative_jacobian(&mechanism, 1u, state, rates, mask, negative_jacobian) == 0);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 1u, state, NULL, negative_jacobian) == 0);
  NEAR(negative_jacobian[0], 10.0, 1e-14);
  NEAR(negative_jacobian[2], -1.0, 1e-14);
  NEAR(negative_jacobian[3], -1.0, 1e-14);
  return 0;
}

static int combined_builtin_and_external_constraints(void) {
  const chem_species species[] = {
    { "A_GAS", 0.0, 0.0, 0.0, 1e-8, 0u }, { "B", 0.0, 0.0, 0.0, 1e-8, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-8, 0u }, { "AEROSOL.A_AQ", 0.0, 0.0, 0.0, 1e-8, 0u }
  };
  chem_mechanism mechanism = one_step_mechanism(species, 4u, 0u, 1u, 0.05);
  const chem_constraint_term mass[] = { { 0u, 1.0 }, { 3u, 1.0 } };
  const chem_constraint_term r[] = { { 1u, 1.0 } }, p[] = { { 2u, 1.0 } };
  const chem_constraint items[] = {
    { .name = "B_C_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 2u,
      .reactants = r, .reactant_count = 1u, .products = p, .product_count = 1u,
      .equilibrium_constant_ref = 5.0, .gas_constant = GAS_CONSTANT, .reference_temperature = T_REF },
    { .name = "gas_aq_mass", .kind = CHEM_CONSTRAINT_LINEAR, .algebraic_species = 3u,
      .terms = mass, .term_count = 2u, .constant = 1.0 }
  };
  const chem_constraint_set set = make_set(4u, items, 2u, 1u);
  unsigned char mask[4];
  CHECK(chem_constraint_set_validate(&set) == 0 && chem_constraint_set_size(&set) == 2u);
  CHECK(chem_constraint_set_species_mask(&set, mask) == 0);
  CHECK(mask[0] && mask[1] && !mask[2] && !mask[3]);
  const double temperature[] = { T_REF }, parameters[] = { 5.0 };
  double state[] = { 0.8, 0.04, 0.2, 0.2 };
  CHECK(integrate_dae(&mechanism, &set, 1u, state, parameters, temperature, 0.5) == 0);
  NEAR(state[0] + state[3], 1.0, 1e-8);
  NEAR(state[2], 5.0 * state[1], 1e-9);
  return 0;
}

static int add_external_model_only_standard_rosenbrock(void) {
  const chem_species species[] = {
    { "A_GAS", 0.0, 0.0, 0.0, 1e-10, 0u }, { "AEROSOL.A_AQ", 0.0, 0.0, 0.0, 1e-10, 0u }
  };
  chem_mechanism mechanism = one_step_mechanism(species, 2u, 0u, 1u, 0.01);
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint item = linear_constraint("inactive_for_standard_rosenbrock", 1u, terms, 2u, 1.0);
  const chem_constraint_set descriptor = make_set(2u, &item, 1u, 0u);
  const double temperature[] = { T_REF }, pressure[] = { 101325.0 };
  unsigned char mask[2];
  double state[] = { 1.0, 0.0 };
  CHECK(chem_constraint_set_species_mask(&descriptor, mask) == 0 && mask[0] && !mask[1]);
  /* The standard ROS API takes only the process mechanism; constraint activation
     is explicit in C and cannot be silently enabled by attaching a descriptor. */
  CHECK(chem_integrate_group(&mechanism, 1u, state, NULL, temperature, pressure, 1.0, NULL) == 0);
  NEAR(state[0], exp(-0.01), 2e-7);
  return 0;
}

static int add_external_model_constraints_only(void) {
  const chem_species species[] = {
    { "A", 0.0, 0.0, 0.0, 1e-8, 0u }, { "B", 0.0, 0.0, 0.0, 1e-8, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-8, 0u }
  };
  chem_mechanism mechanism = one_step_mechanism(species, 3u, 0u, 1u, 0.1);
  const chem_constraint item = equilibrium_constraint("B_C_eq", 2u, 1u, 2u, 2.0, 0.0);
  const chem_constraint_set set = make_set(3u, &item, 1u, 1u);
  const double parameters[] = { 2.0 }, temperature[] = { T_REF };
  double state[] = { 1.0, 0.0, 0.0 };
  CHECK(integrate_dae(&mechanism, &set, 1u, state, parameters, temperature, 20.0) == 0);
  CHECK(state[1] > 0.7);
  NEAR(state[2], 2.0 * state[1], 1e-9);
  return 0;
}

static int multi_grid_cell(void) {
  const chem_species species[] = {
    { "A_GAS", 0.0, 0.0, 0.0, 1e-8, 0u }, { "AEROSOL.A_AQ", 0.0, 0.0, 0.0, 1e-8, 0u }
  };
  chem_mechanism mechanism = one_step_mechanism(species, 2u, 0u, 1u, 0.1);
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint item = linear_constraint("total_mass", 1u, terms, 2u, 1.0);
  const chem_constraint_set set = make_set(2u, &item, 1u, 0u);
  const double temperature[] = { T_REF, T_REF, T_REF };
  double state[] = { 0.9, 0.1, 0.7, 0.3, 0.5, 0.5 };
  CHECK(integrate_dae(&mechanism, &set, 3u, state, NULL, temperature, 50.0) == 0);
  for (size_t cell = 0; cell < 3u; ++cell) NEAR(state[cell * 2u] + state[cell * 2u + 1u], 1.0, 1e-7);
  CHECK(state[0] > state[2] && state[2] > state[4]);
  return 0;
}

static int solve_kinetic_three(double *state) {
  const chem_species species[] = {
    { "A", 0.0, 0.0, 0.0, 1e-4, 0u }, { "B", 0.0, 0.0, 0.0, 1e-4, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-4, 0u }
  };
  const chem_term r[] = { { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 } };
  const chem_term p[] = { { 1u, 1.0 }, { 2u, 1.0 }, { 1u, 1.0 } };
  const chem_reaction reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 0.1, [3] = T_REF } } },
    { 1u, 1u, 1u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 10.0, [3] = T_REF } } },
    { 2u, 1u, 2u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 2.0, [3] = T_REF } } }
  };
  const chem_mechanism mechanism = { .name = "kinetic equilibrium reference", .species_count = 3u,
    .reaction_count = 3u, .reactant_term_count = 3u, .product_term_count = 3u, .species = species,
    .reactions = reactions, .reactant_terms = r, .product_terms = p,
    .temperature = T_REF, .pressure = 101325.0, .relative_tolerance = 1e-8 };
  const double temperature[] = { T_REF }, pressure[] = { 101325.0 };
  state[0] = 1.0; state[1] = 0.0; state[2] = 0.0;
  for (size_t step = 0; step < 200u; ++step)
    if (chem_integrate_group(&mechanism, 1u, state, NULL, temperature, pressure, 1.0, NULL) != 0) return -1;
  return 0;
}

static int solve_conservative_three(double *state, int conserve) {
  const chem_species species[] = {
    { "A", 0.0, 0.0, 0.0, 1e-4, 0u }, { "B", 0.0, 0.0, 0.0, 1e-4, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-4, 0u }
  };
  chem_mechanism mechanism = one_step_mechanism(species, 3u, 0u, 1u, 0.1);
  const chem_constraint_term ratio_r[] = { { 1u, 1.0 } }, ratio_p[] = { { 2u, 1.0 } };
  const chem_constraint_term mass[] = { { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 } };
  const chem_constraint items[] = {
    { .name = "B_C_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 2u,
      .reactants = ratio_r, .reactant_count = 1u, .products = ratio_p, .product_count = 1u,
      .equilibrium_constant_ref = 5.0, .gas_constant = GAS_CONSTANT, .reference_temperature = T_REF,
      .parameter_index = 0u },
    { .name = "mass", .kind = CHEM_CONSTRAINT_LINEAR, .algebraic_species = 1u,
      .terms = mass, .term_count = 3u, .constant = 1.0 }
  };
  const chem_constraint_set set = make_set(3u, items, conserve ? 2u : 1u, 1u);
  const double parameters[] = { 5.0 }, temperature[] = { T_REF };
  state[0] = 1.0; state[1] = 0.0; state[2] = 0.0;
  for (size_t step = 0; step < 200u; ++step)
    if (integrate_dae(&mechanism, &set, 1u, state, parameters, temperature, 1.0) != 0) return -1;
  return 0;
}

static int kinetic_vs_conservative_constraint_convergence(void) {
  double kinetic[3], constrained[3];
  CHECK(solve_kinetic_three(kinetic) == 0);
  CHECK(solve_conservative_three(constrained, 1) == 0);
  NEAR(kinetic[0], 0.0, 1e-6); NEAR(constrained[0], 0.0, 1e-6);
  NEAR(kinetic[1], 1.0 / 6.0, 1e-3); NEAR(constrained[1], 1.0 / 6.0, 1e-3);
  NEAR(kinetic[2], 5.0 / 6.0, 1e-3); NEAR(constrained[2], 5.0 / 6.0, 1e-3);
  NEAR(kinetic[0] + kinetic[1] + kinetic[2], 1.0, 1e-7);
  NEAR(constrained[0] + constrained[1] + constrained[2], 1.0, 1e-7);
  return 0;
}

static int simple_constraint_preserves_ratio_not_mass(void) {
  double kinetic[3], constrained[3];
  CHECK(solve_kinetic_three(kinetic) == 0);
  CHECK(solve_conservative_three(constrained, 0) == 0);
  NEAR(kinetic[2] / kinetic[1], 5.0, 1e-3);
  NEAR(constrained[2] / constrained[1], 5.0, 1e-8);
  NEAR(constrained[0] + constrained[1] + constrained[2], 6.0, 2e-3);
  return 0;
}

static int built_in_vs_external_model_constraint_step_by_step(void) {
  const chem_species species[] = {
    { "A", 0.0, 0.0, 0.0, 1e-4, 0u }, { "B", 0.0, 0.0, 0.0, 1e-4, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-4, 0u }
  };
  chem_mechanism mechanism = one_step_mechanism(species, 3u, 0u, 1u, 0.1);
  const chem_constraint builtin = equilibrium_constraint("built_in_B_C", 2u, 1u, 2u, 5.0, 0.0);
  const chem_constraint external = equilibrium_constraint("external_B_C", 2u, 1u, 2u, 5.0, 0.0);
  const chem_constraint_set built_set = make_set(3u, &builtin, 1u, 1u);
  const chem_constraint_set external_set = make_set(3u, &external, 1u, 1u);
  const double parameters[] = { 5.0 }, temperature[] = { T_REF };
  double built_state[] = { 1.0, 0.0, 0.0 }, external_state[] = { 1.0, 0.0, 0.0 };
  for (size_t step = 0; step < 100u; ++step) {
    CHECK(integrate_dae(&mechanism, &built_set, 1u, built_state, parameters, temperature, 0.5) == 0);
    CHECK(integrate_dae(&mechanism, &external_set, 1u, external_state, parameters, temperature, 0.5) == 0);
    for (size_t s = 0; s < 3u; ++s) CHECK(built_state[s] == external_state[s]);
    NEAR(5.0 * external_state[1] - external_state[2], 0.0, 1e-9);
  }
  return 0;
}

static int multi_equilibrium_kinetic_vs_composed_constraints(void) {
  const chem_species species[] = {
    { "A", 0.0, 0.0, 0.0, 1e-4, 0u }, { "B", 0.0, 0.0, 0.0, 1e-4, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-4, 0u }, { "D", 0.0, 0.0, 0.0, 1e-4, 0u }
  };
  const chem_term kr[] = { { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 }, { 1u, 1.0 }, { 3u, 1.0 } };
  const chem_term kp[] = { { 1u, 1.0 }, { 2u, 1.0 }, { 1u, 1.0 }, { 3u, 1.0 }, { 1u, 1.0 } };
  const chem_reaction kinetic_reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 0.1, [3] = T_REF } } },
    { 1u, 1u, 1u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 10.0, [3] = T_REF } } },
    { 2u, 1u, 2u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 2.0, [3] = T_REF } } },
    { 3u, 1u, 3u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 8.0, [3] = T_REF } } },
    { 4u, 1u, 4u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 8.0 / 3.0, [3] = T_REF } } }
  };
  const chem_mechanism kinetic_mechanism = { .name = "two-equilibrium kinetic", .species_count = 4u,
    .reaction_count = 5u, .reactant_term_count = 5u, .product_term_count = 5u, .species = species,
    .reactions = kinetic_reactions, .reactant_terms = kr, .product_terms = kp,
    .temperature = T_REF, .pressure = 101325.0, .relative_tolerance = 1e-8 };
  chem_mechanism dae_mechanism = one_step_mechanism(species, 4u, 0u, 1u, 0.1);
  const chem_constraint_term br[] = { { 1u, 1.0 } }, cp[] = { { 2u, 1.0 } }, dp[] = { { 3u, 1.0 } };
  const chem_constraint_term mass[] = { { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 }, { 3u, 1.0 } };
  const chem_constraint items[] = {
    { .name = "B_C_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 2u,
      .reactants = br, .reactant_count = 1u, .products = cp, .product_count = 1u,
      .equilibrium_constant_ref = 5.0, .gas_constant = GAS_CONSTANT, .reference_temperature = T_REF,
      .parameter_index = 0u },
    { .name = "B_D_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 3u,
      .reactants = br, .reactant_count = 1u, .products = dp, .product_count = 1u,
      .equilibrium_constant_ref = 3.0, .gas_constant = GAS_CONSTANT, .reference_temperature = T_REF,
      .parameter_index = 1u },
    { .name = "mass", .kind = CHEM_CONSTRAINT_LINEAR, .algebraic_species = 1u,
      .terms = mass, .term_count = 4u, .constant = 1.0 }
  };
  const chem_constraint_set set = make_set(4u, items, 3u, 2u);
  const double parameters[] = { 5.0, 3.0 }, temperature[] = { T_REF }, pressure[] = { 101325.0 };
  double kinetic[] = { 1.0, 0.0, 0.0, 0.0 }, dae[] = { 1.0, 0.0, 0.0, 0.0 };
  for (size_t step = 0; step < 200u; ++step) {
    CHECK(chem_integrate_group(&kinetic_mechanism, 1u, kinetic, NULL, temperature, pressure, 1.0, NULL) == 0);
    CHECK(integrate_dae(&dae_mechanism, &set, 1u, dae, parameters, temperature, 1.0) == 0);
  }
  const double expected_b = 1.0 / 9.0;
  NEAR(kinetic[0], 0.0, 1e-6); NEAR(dae[0], 0.0, 1e-6);
  NEAR(kinetic[1], expected_b, 1e-3); NEAR(dae[1], expected_b, 1e-3);
  NEAR(kinetic[2], 5.0 * expected_b, 1e-3); NEAR(dae[2], 5.0 * expected_b, 1e-3);
  NEAR(kinetic[3], 3.0 * expected_b, 1e-3); NEAR(dae[3], 3.0 * expected_b, 1e-3);
  return 0;
}

static int temperature_dependent_constraint_parameter(void) {
  const chem_species species[] = {
    { "A", 0.0, 0.0, 0.0, 1e-4, 0u }, { "B", 0.0, 0.0, 0.0, 1e-4, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-4, 0u }
  };
  chem_mechanism mechanism = one_step_mechanism(species, 3u, 0u, 1u, 0.1);
  const chem_constraint item = equilibrium_constraint("B_C_eq_T", 2u, 1u, 2u, 2.0, 3000.0 * GAS_CONSTANT);
  const chem_constraint_set set = make_set(3u, &item, 1u, 1u);
  const double temperatures[] = { T_REF, 350.0 };
  for (size_t case_index = 0; case_index < 2u; ++case_index) {
    double keq = 0.0;
    CHECK(chem_constraint_equilibrium_parameter(&item, temperatures[case_index], &keq) == 0);
    double parameters[] = { keq }, state[] = { 1.0, 0.0, 0.0 };
    for (size_t step = 0; step < 200u; ++step)
      CHECK(integrate_dae(&mechanism, &set, 1u, state, parameters, &temperatures[case_index], 1.0) == 0);
    CHECK(state[1] > 0.0);
    NEAR(state[2] / state[1], keq, 2e-8);
  }
  return 0;
}

typedef struct { const char *name; int (*run)(void); } test_case;

int main(int argc, char **argv) {
  static const test_case tests[] = {
    { "AddExternalModelWithConstraints", add_external_model_with_constraints },
    { "AddExternalModelProcessOnly", add_external_model_process_only },
    { "DAESolveEnforcesConservation", dae_solve_enforces_conservation },
    { "AlgebraicRowHoldsOnlyConstraintTerms", algebraic_row_holds_only_constraint_terms },
    { "CombinedBuiltInAndExternalConstraints", combined_builtin_and_external_constraints },
    { "AddExternalModelOnlyStandardRosenbrock", add_external_model_only_standard_rosenbrock },
    { "AddExternalModelConstraintsOnly", add_external_model_constraints_only },
    { "MultiGridCell", multi_grid_cell },
    { "KineticVsConservativeConstraintConvergence", kinetic_vs_conservative_constraint_convergence },
    { "SimpleConstraintPreservesRatioNotMass", simple_constraint_preserves_ratio_not_mass },
    { "BuiltInVsExternalModelConstraintStepByStep", built_in_vs_external_model_constraint_step_by_step },
    { "MultiEquilibriumKineticVsComposedConstraints", multi_equilibrium_kinetic_vs_composed_constraints },
    { "TemperatureDependentConstraintParameter", temperature_dependent_constraint_parameter },
  };
  if (argc != 2) return 2;
  for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    if (strcmp(argv[1], tests[i].name) == 0) return tests[i].run();
  fprintf(stderr, "unknown external-constraint test: %s\n", argv[1]);
  return 2;
}
