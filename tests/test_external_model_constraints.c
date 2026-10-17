/* MICM test/integration/test_external_model_constraints.cpp: external models
 * with processes, constraints, or both (stub_aerosol.c and the models below),
 * alone and with built-in constraints, and finite-difference checks of their
 * Jacobians (a port of micm/util/jacobian_verification.hpp).
 *
 * micm::CpuSolverBuilder<RosenbrockSolverParameters>(options) is a builder
 * with MICM's Matrix and SparseMatrixStandardOrdering (L = 1) and its default
 * LuDecomposition (Doolittle); see cpu_solver_builder(). The tests that leave
 * SetIgnoreUnusedSpecies at MICM's default (true) leave the library's at its
 * default (false), so the library's unused-species check also runs, a
 * stricter build than MICM's. */
#include "harness.h"

#include "stub_aerosol.h"
#include "tmc_external_ext.h"

/* ---- builders, systems, reactions, constraints -------------------------------- */

static tmc_builder *cpu_solver_builder(tmc_rosenbrock_parameters options)
{
  tmc_builder *b = tmc_builder_create_rosenbrock(options);
  tmc_builder_set_width(b, 1);
  tmc_builder_set_lu_algorithm(b, TMC_LU_DOOLITTLE);
  return b;
}

static tmc_solver *build_at(const tmc_builder *b, int line)
{
  tmc_error err;
  tmc_solver *s = tmc_builder_build(b, &err);
  if (!s)
    tmc_test_fail(__FILE__, line, "Build() failed: %s (category %d, code %d)", err.message, (int)err.category, err.code);
  return s;
}
#define BUILD(b) build_at((b), __LINE__)

/* micm::Phase gas_phase{ "gas", { species... } } */
static tmc_phase gas_phase_of(const char *const *names, int n)
{
  tmc_phase p = tmc_phase_create("gas");
  for (int i = 0; i < n; ++i)
  {
    tmc_species s = tmc_species_create(names[i]);
    tmc_phase_add_species(&p, &s);
  }
  return p;
}

/* reactant -> product, ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 } */
static tmc_reaction first_order_reaction(const char *reactant, const char *product, double k, const tmc_phase *phase)
{
  tmc_species r = tmc_species_create(reactant), p = tmc_species_create(product);
  tmc_reaction rx = tmc_reaction_create();
  tmc_reaction_add_reactant(&rx, &r);
  tmc_reaction_add_product(&rx, &p, 1.0);
  tmc_arrhenius_parameters a = tmc_arrhenius_defaults();
  a.A = k;
  a.B = 0.0;
  a.C = 0.0;
  tmc_reaction_set_rate_constant(&rx, tmc_rate_arrhenius(a));
  tmc_reaction_set_phase(&rx, phase);
  return rx;
}

/* EquilibriumConstraint(name, algebraic, { { reactant, 1.0 } }, { { product, 1.0 } }, { K_eq, 0.0 }) */
static tmc_constraint equilibrium_constraint(const char *name, const char *algebraic, const char *reactant, const char *product, double K_eq)
{
  tmc_species alg = tmc_species_create(algebraic);
  tmc_stoich_species r = { tmc_species_create(reactant), 1.0 };
  tmc_stoich_species p = { tmc_species_create(product), 1.0 };
  tmc_vant_hoff vant_hoff = { K_eq, 0.0, 0.0, 0.0 };
  tmc_constraint c;
  tmc_error err;
  if (tmc_equilibrium_constraint(&c, name, &alg, &r, 1, &p, 1, vant_hoff, &err) != TMC_OK)
    tmc_test_fail(__FILE__, __LINE__, "EquilibriumConstraint %s: %s", name, err.message);
  return c;
}

/* state.variables_[cell][state.variable_map_.at(name)] (NaN if unknown) */
static double var(const tmc_state *state, int cell, const char *name)
{
  int i = tmc_state_variable_index(state, name);
  return i < 0 ? NAN : tmc_state_get(state, cell, i);
}

static void set_var(tmc_state *state, int cell, const char *name, double value)
{
  int i = tmc_state_variable_index(state, name);
  if (i < 0)
    tmc_test_fail(__FILE__, __LINE__, "unknown state variable %s", name);
  else
    tmc_state_set(state, cell, i, value);
}

/* ---- name and element lists ----------------------------------------------------- */

static int name_list(char (*names)[TMC_NAME_MAX], int max, const char *const *list, int n)
{
  stub_name_set set;
  set.n = 0;
  for (int i = 0; i < n; ++i)
    stub_name_set_insert(&set, list[i]);
  return stub_name_set_finish(&set, names, max);
}

static int no_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  (void)context;
  (void)names;
  (void)max;
  return 0;
}

static void no_parameter_update(void *context, const tmc_conditions *conditions, int n_cells, tmc_dense *state_parameters)
{
  (void)context;
  (void)conditions;
  (void)n_cells;
  (void)state_parameters;
}

/* ---- EquilibriumConstraintModel ---------------------------------------------------
 * A constraint-only model: K_eq [reactant] - [product] = 0 replaces the
 * product's row. No state variables, no processes. */
typedef struct equilibrium_constraint_model
{
  char reactant[TMC_NAME_MAX], product[TMC_NAME_MAX];
  double K_eq;
  int i_r, i_p, pr_flat, pp_flat;
} equilibrium_constraint_model;

static void equilibrium_constraint_model_init(equilibrium_constraint_model *m, const char *reactant, const char *product, double K_eq)
{
  memset(m, 0, sizeof *m);
  snprintf(m->reactant, sizeof m->reactant, "%s", reactant);
  snprintf(m->product, sizeof m->product, "%s", product);
  m->K_eq = K_eq;
  m->i_r = m->i_p = m->pr_flat = m->pp_flat = -1;
}

static int ecm_algebraic_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  const equilibrium_constraint_model *m = context;
  const char *list[] = { m->product };
  return name_list(names, max, list, 1);
}

static int ecm_dependencies(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  const equilibrium_constraint_model *m = context;
  const char *list[] = { m->reactant, m->product };
  return name_list(names, max, list, 2);
}

static int ecm_elements(void *context, const tmc_name_map *variables, tmc_index_pair *out, int max)
{
  const equilibrium_constraint_model *m = context;
  stub_pair_list list;
  list.n = 0;
  int i_r = tmc_name_map_lookup(variables, m->reactant), i_p = tmc_name_map_lookup(variables, m->product);
  if (i_r < 0 || i_p < 0) /* MICM's at() throws */
    return 0;
  stub_pair_list_add(&list, i_p, i_r);
  stub_pair_list_add(&list, i_p, i_p);
  return stub_pair_list_finish(&list, out, max);
}

static void ecm_finalize(void *context, const tmc_name_map *parameters, const tmc_name_map *variables, const tmc_sparse *jacobian)
{
  equilibrium_constraint_model *m = context;
  (void)parameters;
  m->i_r = tmc_name_map_lookup(variables, m->reactant);
  m->i_p = tmc_name_map_lookup(variables, m->product);
  m->pr_flat = tmc_sparse_element(jacobian, m->i_p, m->i_r);
  m->pp_flat = tmc_sparse_element(jacobian, m->i_p, m->i_p);
}

/* G = K_eq [reactant] - [product] */
static void ecm_residual(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_dense *forcing)
{
  const equilibrium_constraint_model *m = context;
  (void)state_parameters;
  const int rows = stub_dense_rows(forcing);
  for (int row = 0; row < rows; ++row)
    *stub_dense_at(forcing, row, m->i_p) = m->K_eq * stub_dense_value(state_variables, row, m->i_r) - stub_dense_value(state_variables, row, m->i_p);
}

/* jacobian -= dG/dy */
static void ecm_jacobian(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_sparse *jacobian)
{
  const equilibrium_constraint_model *m = context;
  (void)state_parameters;
  (void)state_variables;
  if (m->pr_flat < 0 || m->pp_flat < 0)
    return;
  const int blocks = stub_sparse_blocks(jacobian);
  for (int block = 0; block < blocks; ++block)
  {
    *stub_sparse_at(jacobian, block, m->pr_flat) -= m->K_eq;
    *stub_sparse_at(jacobian, block, m->pp_flat) -= -1.0;
  }
}

static tmc_external_model equilibrium_constraint_model_table(equilibrium_constraint_model *m)
{
  tmc_external_model t;
  memset(&t, 0, sizeof t);
  t.context = m;
  t.constraint_algebraic_variable_names = ecm_algebraic_names;
  t.constraint_species_dependencies = ecm_dependencies;
  t.nonzero_constraint_jacobian_elements = ecm_elements;
  t.constraint_state_parameter_names = no_names;
  t.finalize_constraint_setup = ecm_finalize;
  t.update_constraint_state_parameters = no_parameter_update;
  t.add_constraint_residual = ecm_residual;
  t.subtract_constraint_jacobian = ecm_jacobian;
  return t;
}

/* ---- ConservativeEquilibriumConstraintModel ----------------------------------------
 * For A -> B with a fast B <-> C equilibrium: conservation in the B row,
 * G_B = [A] + [B] + [C] - total, and equilibrium in the C row,
 * G_C = K_eq [B] - [C]. B and C are algebraic; only A evolves. */
typedef struct conservative_equilibrium_model
{
  char a[TMC_NAME_MAX], b[TMC_NAME_MAX], c[TMC_NAME_MAX];
  double K_eq, total;
  int i_a, i_b, i_c;
  int ba_flat, bb_flat, bc_flat, cb_flat, cc_flat;
} conservative_equilibrium_model;

static void conservative_equilibrium_model_init(conservative_equilibrium_model *m, const char *a, const char *b, const char *c, double K_eq, double total)
{
  memset(m, 0, sizeof *m);
  snprintf(m->a, sizeof m->a, "%s", a);
  snprintf(m->b, sizeof m->b, "%s", b);
  snprintf(m->c, sizeof m->c, "%s", c);
  m->K_eq = K_eq;
  m->total = total;
  m->i_a = m->i_b = m->i_c = -1;
  m->ba_flat = m->bb_flat = m->bc_flat = m->cb_flat = m->cc_flat = -1;
}

static int cem_algebraic_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  const conservative_equilibrium_model *m = context;
  const char *list[] = { m->b, m->c };
  return name_list(names, max, list, 2);
}

static int cem_dependencies(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  const conservative_equilibrium_model *m = context;
  const char *list[] = { m->a, m->b, m->c };
  return name_list(names, max, list, 3);
}

static int cem_elements(void *context, const tmc_name_map *variables, tmc_index_pair *out, int max)
{
  const conservative_equilibrium_model *m = context;
  stub_pair_list list;
  list.n = 0;
  int i_a = tmc_name_map_lookup(variables, m->a), i_b = tmc_name_map_lookup(variables, m->b), i_c = tmc_name_map_lookup(variables, m->c);
  if (i_a < 0 || i_b < 0 || i_c < 0)
    return 0;
  /* B row (conservation) */
  stub_pair_list_add(&list, i_b, i_a);
  stub_pair_list_add(&list, i_b, i_b);
  stub_pair_list_add(&list, i_b, i_c);
  /* C row (equilibrium) */
  stub_pair_list_add(&list, i_c, i_b);
  stub_pair_list_add(&list, i_c, i_c);
  return stub_pair_list_finish(&list, out, max);
}

static void cem_finalize(void *context, const tmc_name_map *parameters, const tmc_name_map *variables, const tmc_sparse *jacobian)
{
  conservative_equilibrium_model *m = context;
  (void)parameters;
  m->i_a = tmc_name_map_lookup(variables, m->a);
  m->i_b = tmc_name_map_lookup(variables, m->b);
  m->i_c = tmc_name_map_lookup(variables, m->c);
  m->ba_flat = tmc_sparse_element(jacobian, m->i_b, m->i_a);
  m->bb_flat = tmc_sparse_element(jacobian, m->i_b, m->i_b);
  m->bc_flat = tmc_sparse_element(jacobian, m->i_b, m->i_c);
  m->cb_flat = tmc_sparse_element(jacobian, m->i_c, m->i_b);
  m->cc_flat = tmc_sparse_element(jacobian, m->i_c, m->i_c);
}

static void cem_residual(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_dense *forcing)
{
  const conservative_equilibrium_model *m = context;
  (void)state_parameters;
  const int rows = stub_dense_rows(forcing);
  for (int row = 0; row < rows; ++row)
  {
    const double a = stub_dense_value(state_variables, row, m->i_a);
    const double b = stub_dense_value(state_variables, row, m->i_b);
    const double c = stub_dense_value(state_variables, row, m->i_c);
    *stub_dense_at(forcing, row, m->i_b) = a + b + c - m->total;
    *stub_dense_at(forcing, row, m->i_c) = m->K_eq * b - c;
  }
}

static void cem_jacobian(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_sparse *jacobian)
{
  const conservative_equilibrium_model *m = context;
  (void)state_parameters;
  (void)state_variables;
  if (m->ba_flat < 0 || m->bb_flat < 0 || m->bc_flat < 0 || m->cb_flat < 0 || m->cc_flat < 0)
    return;
  const int blocks = stub_sparse_blocks(jacobian);
  for (int block = 0; block < blocks; ++block)
  {
    *stub_sparse_at(jacobian, block, m->ba_flat) -= 1.0;
    *stub_sparse_at(jacobian, block, m->bb_flat) -= 1.0;
    *stub_sparse_at(jacobian, block, m->bc_flat) -= 1.0;
    *stub_sparse_at(jacobian, block, m->cb_flat) -= m->K_eq;
    *stub_sparse_at(jacobian, block, m->cc_flat) -= -1.0;
  }
}

static tmc_external_model conservative_equilibrium_model_table(conservative_equilibrium_model *m)
{
  tmc_external_model t;
  memset(&t, 0, sizeof t);
  t.context = m;
  t.constraint_algebraic_variable_names = cem_algebraic_names;
  t.constraint_species_dependencies = cem_dependencies;
  t.nonzero_constraint_jacobian_elements = cem_elements;
  t.constraint_state_parameter_names = no_names;
  t.finalize_constraint_setup = cem_finalize;
  t.update_constraint_state_parameters = no_parameter_update;
  t.add_constraint_residual = cem_residual;
  t.subtract_constraint_jacobian = cem_jacobian;
  return t;
}

/* ---- MassConservationModel ----------------------------------------------------------
 * Replaces the row of the controlled species with G = sum(all species) - total. */
#define MASS_MAX_SPECIES 8
typedef struct mass_conservation_model
{
  char controlled[TMC_NAME_MAX];
  int n_species;
  char species[MASS_MAX_SPECIES][TMC_NAME_MAX];
  double total;
  int i_ctrl;
  int indices[MASS_MAX_SPECIES];
  int flat_ids[MASS_MAX_SPECIES];
} mass_conservation_model;

static void mass_conservation_model_init(mass_conservation_model *m, const char *controlled, const char *const *all_species, int n, double total)
{
  memset(m, 0, sizeof *m);
  snprintf(m->controlled, sizeof m->controlled, "%s", controlled);
  m->n_species = n < MASS_MAX_SPECIES ? n : MASS_MAX_SPECIES;
  for (int i = 0; i < m->n_species; ++i)
    snprintf(m->species[i], sizeof m->species[i], "%s", all_species[i]);
  m->total = total;
  m->i_ctrl = -1;
}

static int mcm_algebraic_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  const mass_conservation_model *m = context;
  const char *list[] = { m->controlled };
  return name_list(names, max, list, 1);
}

static int mcm_dependencies(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  const mass_conservation_model *m = context;
  const char *list[MASS_MAX_SPECIES];
  for (int i = 0; i < m->n_species; ++i)
    list[i] = m->species[i];
  return name_list(names, max, list, m->n_species);
}

static int mcm_elements(void *context, const tmc_name_map *variables, tmc_index_pair *out, int max)
{
  const mass_conservation_model *m = context;
  stub_pair_list list;
  list.n = 0;
  int i_ctrl = tmc_name_map_lookup(variables, m->controlled);
  if (i_ctrl < 0)
    return 0;
  for (int i = 0; i < m->n_species; ++i)
  {
    int j = tmc_name_map_lookup(variables, m->species[i]);
    if (j < 0)
      return 0;
    stub_pair_list_add(&list, i_ctrl, j);
  }
  return stub_pair_list_finish(&list, out, max);
}

static void mcm_finalize(void *context, const tmc_name_map *parameters, const tmc_name_map *variables, const tmc_sparse *jacobian)
{
  mass_conservation_model *m = context;
  (void)parameters;
  m->i_ctrl = tmc_name_map_lookup(variables, m->controlled);
  for (int i = 0; i < m->n_species; ++i)
  {
    m->indices[i] = tmc_name_map_lookup(variables, m->species[i]);
    m->flat_ids[i] = tmc_sparse_element(jacobian, m->i_ctrl, m->indices[i]);
  }
}

/* sum = -total, then += each species in order; the controlled row takes the sum */
static void mcm_residual(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_dense *forcing)
{
  const mass_conservation_model *m = context;
  (void)state_parameters;
  const int rows = stub_dense_rows(forcing);
  for (int row = 0; row < rows; ++row)
  {
    double sum = -m->total;
    for (int i = 0; i < m->n_species; ++i)
      sum += stub_dense_value(state_variables, row, m->indices[i]);
    *stub_dense_at(forcing, row, m->i_ctrl) = sum;
  }
}

static void mcm_jacobian(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_sparse *jacobian)
{
  const mass_conservation_model *m = context;
  (void)state_parameters;
  (void)state_variables;
  const int blocks = stub_sparse_blocks(jacobian);
  for (int i = 0; i < m->n_species; ++i)
  {
    if (m->flat_ids[i] < 0)
      continue;
    for (int block = 0; block < blocks; ++block)
      *stub_sparse_at(jacobian, block, m->flat_ids[i]) -= 1.0;
  }
}

static tmc_external_model mass_conservation_model_table(mass_conservation_model *m)
{
  tmc_external_model t;
  memset(&t, 0, sizeof t);
  t.context = m;
  t.constraint_algebraic_variable_names = mcm_algebraic_names;
  t.constraint_species_dependencies = mcm_dependencies;
  t.nonzero_constraint_jacobian_elements = mcm_elements;
  t.constraint_state_parameter_names = no_names;
  t.finalize_constraint_setup = mcm_finalize;
  t.update_constraint_state_parameters = no_parameter_update;
  t.add_constraint_residual = mcm_residual;
  t.subtract_constraint_jacobian = mcm_jacobian;
  return t;
}

/* ---- TemperatureDependentEquilibriumModel -------------------------------------------
 * K_eq(T) [reactant] - [product] = 0, with K_eq(T) = K_eq_ref
 * exp(delta_H / R (1 / T_ref - 1 / T)) kept in the state parameter
 * "<product>_K_eq" and refreshed from the conditions. */
typedef struct temperature_equilibrium_model
{
  char reactant[TMC_NAME_MAX], product[TMC_NAME_MAX], param_name[TMC_NAME_MAX];
  double K_eq_ref, delta_H_over_R, T_ref;
  int i_r, i_p, i_K, pr_flat, pp_flat;
} temperature_equilibrium_model;

static void temperature_equilibrium_model_init(temperature_equilibrium_model *m, const char *reactant, const char *product, double K_eq_ref, double delta_H_over_R, double T_ref)
{
  memset(m, 0, sizeof *m);
  snprintf(m->reactant, sizeof m->reactant, "%s", reactant);
  snprintf(m->product, sizeof m->product, "%s", product);
  char label[2 * TMC_NAME_MAX];
  snprintf(label, sizeof label, "%s_K_eq", product);
  snprintf(m->param_name, sizeof m->param_name, "%.*s", TMC_NAME_MAX - 1, label);
  m->K_eq_ref = K_eq_ref;
  m->delta_H_over_R = delta_H_over_R;
  m->T_ref = T_ref;
  m->i_r = m->i_p = m->i_K = m->pr_flat = m->pp_flat = -1;
}

static int tem_algebraic_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  const temperature_equilibrium_model *m = context;
  const char *list[] = { m->product };
  return name_list(names, max, list, 1);
}

static int tem_dependencies(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  const temperature_equilibrium_model *m = context;
  const char *list[] = { m->reactant, m->product };
  return name_list(names, max, list, 2);
}

static int tem_elements(void *context, const tmc_name_map *variables, tmc_index_pair *out, int max)
{
  const temperature_equilibrium_model *m = context;
  stub_pair_list list;
  list.n = 0;
  int i_r = tmc_name_map_lookup(variables, m->reactant), i_p = tmc_name_map_lookup(variables, m->product);
  if (i_r < 0 || i_p < 0)
    return 0;
  stub_pair_list_add(&list, i_p, i_r);
  stub_pair_list_add(&list, i_p, i_p);
  return stub_pair_list_finish(&list, out, max);
}

static int tem_parameter_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  const temperature_equilibrium_model *m = context;
  const char *list[] = { m->param_name };
  return name_list(names, max, list, 1);
}

static void tem_finalize(void *context, const tmc_name_map *parameters, const tmc_name_map *variables, const tmc_sparse *jacobian)
{
  temperature_equilibrium_model *m = context;
  m->i_r = tmc_name_map_lookup(variables, m->reactant);
  m->i_p = tmc_name_map_lookup(variables, m->product);
  m->i_K = tmc_name_map_lookup(parameters, m->param_name);
  m->pr_flat = tmc_sparse_element(jacobian, m->i_p, m->i_r);
  m->pp_flat = tmc_sparse_element(jacobian, m->i_p, m->i_p);
}

static void tem_update(void *context, const tmc_conditions *conditions, int n_cells, tmc_dense *state_parameters)
{
  const temperature_equilibrium_model *m = context;
  if (m->i_K < 0)
    return;
  for (int cell = 0; cell < n_cells; ++cell)
    *stub_dense_at(state_parameters, cell, m->i_K) = m->K_eq_ref * exp(m->delta_H_over_R * (1.0 / m->T_ref - 1.0 / conditions[cell].temperature));
}

static void tem_residual(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_dense *forcing)
{
  const temperature_equilibrium_model *m = context;
  const int rows = stub_dense_rows(forcing);
  for (int row = 0; row < rows; ++row)
    *stub_dense_at(forcing, row, m->i_p) = stub_dense_value(state_parameters, row, m->i_K) * stub_dense_value(state_variables, row, m->i_r) - stub_dense_value(state_variables, row, m->i_p);
}

static void tem_jacobian(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_sparse *jacobian)
{
  const temperature_equilibrium_model *m = context;
  (void)state_variables;
  if (m->pr_flat < 0 || m->pp_flat < 0)
    return;
  const int blocks = stub_sparse_blocks(jacobian);
  for (int block = 0; block < blocks; ++block)
  {
    *stub_sparse_at(jacobian, block, m->pr_flat) -= stub_dense_value(state_parameters, block, m->i_K);
    *stub_sparse_at(jacobian, block, m->pp_flat) -= -1.0;
  }
}

static tmc_external_model temperature_equilibrium_model_table(temperature_equilibrium_model *m)
{
  tmc_external_model t;
  memset(&t, 0, sizeof t);
  t.context = m;
  t.constraint_algebraic_variable_names = tem_algebraic_names;
  t.constraint_species_dependencies = tem_dependencies;
  t.nonzero_constraint_jacobian_elements = tem_elements;
  t.constraint_state_parameter_names = tem_parameter_names;
  t.finalize_constraint_setup = tem_finalize;
  t.update_constraint_state_parameters = tem_update;
  t.add_constraint_residual = tem_residual;
  t.subtract_constraint_jacobian = tem_jacobian;
  return t;
}

/* ---- NoOpProcessModel -----------------------------------------------------------------
 * A process-only model with no state, no Jacobian elements, and no forcing. */
static int noop_elements(void *context, const tmc_name_map *variables, tmc_index_pair *out, int max)
{
  (void)context;
  (void)variables;
  (void)out;
  (void)max;
  return 0;
}

static void noop_finalize(void *context, const tmc_name_map *parameters, const tmc_name_map *variables, const tmc_sparse *jacobian)
{
  (void)context;
  (void)parameters;
  (void)variables;
  (void)jacobian;
}

static void noop_forcing(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_dense *forcing)
{
  (void)context;
  (void)state_parameters;
  (void)state_variables;
  (void)forcing;
}

static void noop_jacobian(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_sparse *jacobian)
{
  (void)context;
  (void)state_parameters;
  (void)state_variables;
  (void)jacobian;
}

static tmc_external_model noop_process_model_table(void)
{
  tmc_external_model t;
  memset(&t, 0, sizeof t);
  t.species_used = no_names;
  t.nonzero_jacobian_elements = noop_elements;
  t.finalize_process_setup = noop_finalize;
  t.update_state_parameters = no_parameter_update;
  t.add_forcing_terms = noop_forcing;
  t.subtract_jacobian_terms = noop_jacobian;
  return t;
}

/* ---- micm/util/jacobian_verification.hpp -----------------------------------------------
 * A dense finite-difference Jacobian (central differences, forward ones where
 * x - h would be negative) compared with an analytical one, which holds
 * -df/dx as MICM's Jacobians do. */
typedef struct jacobian_comparison
{
  bool passed;
  double max_abs_error, max_rel_error;
  int worst_block, worst_row, worst_col;
  double worst_analytical, worst_fd;
} jacobian_comparison;

typedef void (*forcing_function)(void *context, const tmc_dense *variables, tmc_dense *forcing);

/* result[block][row * num_species + col] = d f_row / d x_col */
static tmc_dense *finite_difference_jacobian(forcing_function forcing_func, void *context, const tmc_dense *base_variables, int num_species, double perturbation)
{
  const int num_blocks = base_variables->n_rows, L = base_variables->L;
  tmc_dense *result = tmc_dense_create(num_blocks, num_species * num_species, L, 0.0);
  tmc_dense *vars_plus = tmc_dense_create(num_blocks, num_species, L, 0.0);
  tmc_dense *vars_minus = tmc_dense_create(num_blocks, num_species, L, 0.0);
  tmc_dense *forcing_plus = tmc_dense_create(num_blocks, num_species, L, 0.0);
  tmc_dense *forcing_minus = tmc_dense_create(num_blocks, num_species, L, 0.0);
  tmc_dense *forcing_base = tmc_dense_create(num_blocks, num_species, L, 0.0);
  for (int col = 0; col < num_species; ++col)
  {
    tmc_dense_copy_into(vars_plus, base_variables);
    tmc_dense_copy_into(vars_minus, base_variables);
    for (int block = 0; block < num_blocks; ++block)
    {
      double x_j = tmc_dense_get(base_variables, block, col);
      double h = perturbation * fmax(1.0, fabs(x_j));
      tmc_dense_set(vars_plus, block, col, x_j + h);
      tmc_dense_set(vars_minus, block, col, x_j - h);
    }
    tmc_dense_fill(forcing_plus, 0.0);
    tmc_dense_fill(forcing_minus, 0.0);
    forcing_func(context, vars_plus, forcing_plus);
    forcing_func(context, vars_minus, forcing_minus);
    for (int block = 0; block < num_blocks; ++block)
    {
      double x_j = tmc_dense_get(base_variables, block, col);
      double h = perturbation * fmax(1.0, fabs(x_j));
      if (x_j - h >= 0.0)
      {
        double two_h = 2.0 * h;
        for (int row = 0; row < num_species; ++row)
          tmc_dense_set(result, block, row * num_species + col, (tmc_dense_get(forcing_plus, block, row) - tmc_dense_get(forcing_minus, block, row)) / two_h);
      }
      else
      {
        tmc_dense_fill(forcing_base, 0.0);
        forcing_func(context, base_variables, forcing_base);
        for (int row = 0; row < num_species; ++row)
          tmc_dense_set(result, block, row * num_species + col, (tmc_dense_get(forcing_plus, block, row) - tmc_dense_get(forcing_base, block, row)) / h);
      }
    }
  }
  tmc_dense_free(vars_plus);
  tmc_dense_free(vars_minus);
  tmc_dense_free(forcing_plus);
  tmc_dense_free(forcing_minus);
  tmc_dense_free(forcing_base);
  return result;
}

/* |analytical - fd| < atol + rtol max(|analytical|, |fd|) over the analytical pattern */
static jacobian_comparison compare_jacobian_to_finite_difference(const tmc_sparse *analytical, const tmc_dense *fd, int num_species, double atol, double rtol)
{
  jacobian_comparison r;
  memset(&r, 0, sizeof r);
  r.passed = true;
  for (int block = 0; block < analytical->n_blocks; ++block)
    for (int row = 0; row < num_species; ++row)
      for (int col = 0; col < num_species; ++col)
      {
        if (tmc_sparse_is_zero(analytical, row, col))
          continue;
        double analytical_val = -tmc_sparse_get(analytical, block, row, col);
        double fd_val = tmc_dense_get(fd, block, row * num_species + col);
        double abs_error = fabs(analytical_val - fd_val);
        double scale = fmax(fabs(analytical_val), fabs(fd_val));
        double rel_error = scale > 0.0 ? abs_error / scale : 0.0;
        if (abs_error > r.max_abs_error)
        {
          r.max_abs_error = abs_error;
          r.worst_block = block;
          r.worst_row = row;
          r.worst_col = col;
          r.worst_analytical = analytical_val;
          r.worst_fd = fd_val;
        }
        if (rel_error > r.max_rel_error)
          r.max_rel_error = rel_error;
        if (abs_error > atol + rtol * scale)
          r.passed = false;
      }
  return r;
}

/* no finite-difference entry above the threshold outside the analytical pattern */
static jacobian_comparison check_jacobian_sparsity_completeness(const tmc_sparse *analytical, const tmc_dense *fd, int num_species, double threshold)
{
  jacobian_comparison r;
  memset(&r, 0, sizeof r);
  r.passed = true;
  for (int block = 0; block < analytical->n_blocks; ++block)
    for (int row = 0; row < num_species; ++row)
      for (int col = 0; col < num_species; ++col)
      {
        if (!tmc_sparse_is_zero(analytical, row, col))
          continue;
        double fd_val = fabs(tmc_dense_get(fd, block, row * num_species + col));
        if (fd_val > threshold)
        {
          r.passed = false;
          if (fd_val > r.max_abs_error)
          {
            r.max_abs_error = fd_val;
            r.worst_block = block;
            r.worst_row = row;
            r.worst_col = col;
            r.worst_analytical = 0.0;
            r.worst_fd = fd_val;
          }
        }
      }
  return r;
}

/* MICM's defaults for double precision */
#define FD_PERTURBATION 1.0e-8
#define FD_ATOL 1.0e-7
#define FD_RTOL 1.0e-7
#define FD_SPARSITY_THRESHOLD 1.0e-6

/* A sparse matrix of n_blocks blocks (standard ordering) with the elements a
 * model declares. */
static tmc_sparse *jacobian_of(const tmc_index_pair *elements, int n_elements, int num_species, int n_blocks)
{
  tmc_sparse_builder *sb = tmc_sparse_builder_create(num_species);
  tmc_sparse_builder_set_number_of_blocks(sb, n_blocks);
  tmc_sparse_builder_set_initial_value(sb, 0.0);
  for (int i = 0; i < n_elements; ++i)
    tmc_sparse_builder_with_element(sb, elements[i].row, elements[i].col, NULL);
  tmc_sparse *m = tmc_sparse_build(sb, TMC_SPARSE_CSR, 1, false);
  tmc_sparse_builder_free(sb);
  return m;
}

typedef struct model_forcing
{
  const tmc_external_model *model;
  const tmc_dense *parameters;
} model_forcing;

static void process_forcing(void *context, const tmc_dense *variables, tmc_dense *forcing)
{
  const model_forcing *f = context;
  f->model->add_forcing_terms(f->model->context, f->parameters, variables, forcing);
}

static void constraint_residual(void *context, const tmc_dense *variables, tmc_dense *forcing)
{
  const model_forcing *f = context;
  f->model->add_constraint_residual(f->model->context, f->parameters, variables, forcing);
}

#define EXPECT_COMPARISON(c, what)                                                                       \
  do                                                                                                     \
  {                                                                                                      \
    if (!(c).passed)                                                                                     \
      tmc_test_fail(__FILE__, __LINE__, "%s: block=%d row=%d col=%d analytical=%.17g fd=%.17g", what, (c).worst_block, \
                    (c).worst_row, (c).worst_col, (c).worst_analytical, (c).worst_fd);                    \
  } while (0)

/* ---- AddExternalModelWithConstraints: both processes and constraints ---------------- */

TMC_TEST(add_external_model_with_constraints)
{
  const char *species[] = { "A_GAS" };
  tmc_phase gas = gas_phase_of(species, 1);
  tmc_system system = tmc_system_create(&gas);
  stub_aerosol_constraints aerosol;
  stub_aerosol_constraints_init_with_total(&aerosol, 0.01, 1.0);
  tmc_external_model model = stub_aerosol_constraints_model(&aerosol);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, NULL, 0);
  tmc_builder_add_external_model(b, &model);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  if (state)
  {
    TMC_EXPECT_EQ(tmc_state_num_variables(state), 2); /* A_GAS, AEROSOL.A_AQ */
    TMC_EXPECT_EQ(tmc_solver_num_constraints(solver), 1);
    int i_gas = tmc_state_variable_index(state, "A_GAS"), i_aq = tmc_state_variable_index(state, "AEROSOL.A_AQ");
    TMC_EXPECT_TRUE(i_gas >= 0 && i_aq >= 0);
    if (i_gas >= 0 && i_aq >= 0)
    {
      /* the A_AQ row is algebraic */
      TMC_EXPECT_DOUBLE_EQ(tmc_solver_mass_matrix_diagonal(solver)[i_gas], 1.0);
      TMC_EXPECT_DOUBLE_EQ(tmc_solver_mass_matrix_diagonal(solver)[i_aq], 0.0);
    }
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

/* ---- AddExternalModelProcessOnly: no total mass, no constraints --------------------- */

TMC_TEST(add_external_model_process_only)
{
  const char *species[] = { "A_GAS" };
  tmc_phase gas = gas_phase_of(species, 1);
  tmc_system system = tmc_system_create(&gas);
  stub_aerosol_constraints aerosol;
  stub_aerosol_constraints_init(&aerosol, 0.01);
  tmc_external_model model = stub_aerosol_constraints_model(&aerosol);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, NULL, 0);
  tmc_builder_add_external_model(b, &model);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  if (state)
  {
    /* pure ODE mode */
    TMC_EXPECT_EQ(tmc_solver_num_constraints(solver), 0);
    int i_gas = tmc_state_variable_index(state, "A_GAS"), i_aq = tmc_state_variable_index(state, "AEROSOL.A_AQ");
    TMC_EXPECT_TRUE(i_gas >= 0 && i_aq >= 0);
    if (i_gas >= 0 && i_aq >= 0)
    {
      TMC_EXPECT_DOUBLE_EQ(tmc_solver_mass_matrix_diagonal(solver)[i_gas], 1.0);
      TMC_EXPECT_DOUBLE_EQ(tmc_solver_mass_matrix_diagonal(solver)[i_aq], 1.0);
    }
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

/* ---- DAESolveEnforcesConservation -------------------------------------------------------- */

TMC_TEST(dae_solve_enforces_conservation)
{
  const char *species[] = { "A_GAS" };
  tmc_phase gas = gas_phase_of(species, 1);
  tmc_system system = tmc_system_create(&gas);
  const double total = 1.0, k = 0.1;
  stub_aerosol_constraints aerosol;
  stub_aerosol_constraints_init_with_total(&aerosol, k, total);
  tmc_external_model model = stub_aerosol_constraints_model(&aerosol);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, NULL, 0);
  tmc_builder_add_external_model(b, &model);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  if (state)
  {
    /* most of the mass in the gas phase */
    set_var(state, 0, "A_GAS", 0.9);
    set_var(state, 0, "AEROSOL.A_AQ", 0.1);
    tmc_state_conditions(state)[0].temperature = 298.0;
    tmc_state_conditions(state)[0].pressure = 101325.0;
    const double dt = 10.0;
    for (int step = 0; step < 20; ++step)
    {
      tmc_result result = tmc_solver_solve(solver, dt, state);
      if (result.state != TMC_CONVERGED)
        tmc_test_fail(__FILE__, __LINE__, "step %d: %s", step, tmc_solver_state_name(result.state));
      double sum = var(state, 0, "A_GAS") + var(state, 0, "AEROSOL.A_AQ");
      if (!(fabs(sum - total) <= 1e-4))
        tmc_test_fail(__FILE__, __LINE__, "conservation violated at step %d: sum %.17g", step, sum);
    }
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

/* ---- AlgebraicRowHoldsOnlyConstraintTerms -------------------------------------------------
 * The process A_GAS -> A_AQ also writes into the algebraic A_AQ row; after the
 * rates and the constraints run in the Rosenbrock sequence, the A_AQ row holds
 * the constraint only and the A_GAS row keeps the process terms. */
TMC_TEST(algebraic_row_holds_only_constraint_terms)
{
  const char *species[] = { "A_GAS" };
  tmc_phase gas = gas_phase_of(species, 1);
  tmc_system system = tmc_system_create(&gas);
  const double k = 10.0, total = 1.0;
  stub_aerosol_constraints aerosol;
  stub_aerosol_constraints_init_with_total(&aerosol, k, total);
  tmc_external_model model = stub_aerosol_constraints_model(&aerosol);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, NULL, 0);
  tmc_builder_set_reorder_state(b, false);
  tmc_builder_add_external_model(b, &model);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  tmc_dense *forcing = NULL;
  if (state)
  {
    const int i_gas = tmc_state_variable_index(state, "A_GAS"), i_aq = tmc_state_variable_index(state, "AEROSOL.A_AQ");
    TMC_EXPECT_TRUE(i_gas >= 0 && i_aq >= 0);
    if (i_gas >= 0 && i_aq >= 0)
    {
      const double gas_value = 0.7, aq_value = 0.2, tol = 1.0e-12;
      tmc_state_set(state, 0, i_gas, gas_value);
      tmc_state_set(state, 0, i_aq, aq_value);
      tmc_dense *y = tmc_state_variables(state);

      /* forcing */
      forcing = tmc_dense_create(1, tmc_state_num_variables(state), tmc_state_width(state), 0.0);
      tmc_dense_fill(forcing, 0.0);
      tmc_solver_add_forcing_terms(solver, state, y, forcing);
      tmc_solver_add_constraint_terms(solver, state, y, forcing);
      TMC_EXPECT_NEAR(tmc_dense_get(forcing, 0, i_gas), -k * gas_value, tol);
      TMC_EXPECT_NEAR(tmc_dense_get(forcing, 0, i_aq), gas_value + aq_value - total, tol);

      /* Jacobian (the solver stores -J) */
      tmc_sparse *jacobian = tmc_state_jacobian(state);
      tmc_sparse_fill(jacobian, 0.0);
      tmc_solver_subtract_jacobian_terms(solver, state, y, jacobian);
      tmc_solver_subtract_constraint_jacobian_terms(solver, state, y, jacobian);
      TMC_EXPECT_NEAR(tmc_sparse_get(jacobian, 0, i_gas, i_gas), k, tol);
      /* the external process term is not left in the algebraic row */
      TMC_EXPECT_NEAR(tmc_sparse_get(jacobian, 0, i_aq, i_gas), -1.0, tol);
      TMC_EXPECT_NEAR(tmc_sparse_get(jacobian, 0, i_aq, i_aq), -1.0, tol);
    }
  }
  tmc_dense_free(forcing);
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

/* ---- CombinedBuiltInAndExternalConstraints -------------------------------------------- */

TMC_TEST(combined_built_in_and_external_constraints)
{
  const char *species[] = { "A_GAS", "B", "C" };
  tmc_phase gas = gas_phase_of(species, 3);
  tmc_system system = tmc_system_create(&gas);
  stub_aerosol_constraints aerosol;
  stub_aerosol_constraints_init_with_total(&aerosol, 0.01, 1.0);
  tmc_external_model model = stub_aerosol_constraints_model(&aerosol);

  /* built-in constraint: B <-> C equilibrium; process A_GAS -> B */
  tmc_constraint constraint = equilibrium_constraint("B_C_eq", "C", "B", "C", 5.0);
  tmc_reaction rxn = first_order_reaction("A_GAS", "B", 0.05, &gas);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, &rxn, 1);
  tmc_builder_add_constraint(b, &constraint);
  tmc_builder_add_external_model(b, &model);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  if (state)
  {
    /* two algebraic variables: C (built-in) and AEROSOL.A_AQ (external) */
    TMC_EXPECT_EQ(tmc_solver_num_constraints(solver), 2);
    int i_aq = tmc_state_variable_index(state, "AEROSOL.A_AQ"), i_c = tmc_state_variable_index(state, "C");
    TMC_EXPECT_TRUE(i_aq >= 0 && i_c >= 0);
    if (i_aq >= 0 && i_c >= 0)
    {
      TMC_EXPECT_DOUBLE_EQ(tmc_solver_mass_matrix_diagonal(solver)[i_aq], 0.0);
      TMC_EXPECT_DOUBLE_EQ(tmc_solver_mass_matrix_diagonal(solver)[i_c], 0.0);
    }
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_reaction_free(&rxn);
  tmc_constraint_free(&constraint);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

/* ---- AddExternalModelOnlyStandardRosenbrock: processes only ---------------------------- */

TMC_TEST(add_external_model_only_standard_rosenbrock)
{
  const char *species[] = { "A_GAS" };
  tmc_phase gas = gas_phase_of(species, 1);
  tmc_system system = tmc_system_create(&gas);
  stub_aerosol_constraints aerosol;
  stub_aerosol_constraints_init(&aerosol, 0.01);
  tmc_external_model model = stub_aerosol_constraints_model(&aerosol);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, NULL, 0);
  tmc_builder_add_external_model(b, &model);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  if (state)
    TMC_EXPECT_EQ(tmc_solver_num_constraints(solver), 0);
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

/* ---- AddExternalModelConstraintsOnly ----------------------------------------------------- */

TMC_TEST(add_external_model_constraints_only)
{
  const char *species[] = { "A_GAS" };
  tmc_phase gas = gas_phase_of(species, 1);
  tmc_system system = tmc_system_create(&gas);
  const double total = 1.0;
  stub_aerosol_constraints aerosol;
  stub_aerosol_constraints_init_with_total(&aerosol, 0.01, total);
  tmc_external_model model = stub_aerosol_constraints_model(&aerosol);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, NULL, 0);
  tmc_builder_add_external_model(b, &model);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  if (state)
  {
    /* the constraint is active */
    TMC_EXPECT_EQ(tmc_solver_num_constraints(solver), 1);
    int i_aq = tmc_state_variable_index(state, "AEROSOL.A_AQ");
    TMC_EXPECT_TRUE(i_aq >= 0);
    if (i_aq >= 0)
    {
      TMC_EXPECT_DOUBLE_EQ(tmc_solver_mass_matrix_diagonal(solver)[i_aq], 0.0);
      /* solve and check conservation */
      set_var(state, 0, "A_GAS", 0.8);
      tmc_state_set(state, 0, i_aq, 0.2);
      tmc_state_conditions(state)[0].temperature = 298.0;
      tmc_state_conditions(state)[0].pressure = 101325.0;
      tmc_result result = tmc_solver_solve(solver, 50.0, state);
      TMC_EXPECT_EQ(result.state, TMC_CONVERGED);
      double sum = var(state, 0, "A_GAS") + tmc_state_get(state, 0, i_aq);
      TMC_EXPECT_NEAR(sum, total, 1e-4);
    }
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

/* ---- MultiGridCell ------------------------------------------------------------------------ */

TMC_TEST(multi_grid_cell)
{
  const char *species[] = { "A_GAS" };
  tmc_phase gas = gas_phase_of(species, 1);
  tmc_system system = tmc_system_create(&gas);
  const double total = 1.0;
  stub_aerosol_constraints aerosol;
  stub_aerosol_constraints_init_with_total(&aerosol, 0.1, total);
  tmc_external_model model = stub_aerosol_constraints_model(&aerosol);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, NULL, 0);
  tmc_builder_add_external_model(b, &model);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  const int num_cells = 3;
  tmc_state *state = solver ? tmc_solver_create_state(solver, num_cells) : NULL;
  if (state)
  {
    /* different initial conditions per cell */
    for (int c = 0; c < num_cells; ++c)
    {
      double gas_frac = 0.9 - 0.2 * c;
      set_var(state, c, "A_GAS", gas_frac);
      set_var(state, c, "AEROSOL.A_AQ", total - gas_frac);
      tmc_state_conditions(state)[c].temperature = 298.0;
      tmc_state_conditions(state)[c].pressure = 101325.0;
    }
    tmc_result result = tmc_solver_solve(solver, 50.0, state);
    TMC_EXPECT_EQ(result.state, TMC_CONVERGED);
    for (int c = 0; c < num_cells; ++c)
    {
      double sum = var(state, c, "A_GAS") + var(state, c, "AEROSOL.A_AQ");
      if (!(fabs(sum - total) <= 1e-4))
        tmc_test_fail(__FILE__, __LINE__, "conservation violated in cell %d: sum %.17g", c, sum);
    }
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

/* ---- Convergence: kinetic (process-only) against constraint systems -----------------------
 * A -> B (slow driver, k = 0.1), B <-> C (fast, K_eq = 5). A single ratio
 * constraint K_eq [B] - [C] = 0 keeps the ratio but not the mass; with a
 * conservation constraint as well, B and C are algebraic and the DAE matches
 * the kinetic system. Steady state: [A] = 0, [B] = 1 / (1 + K_eq),
 * [C] = K_eq / (1 + K_eq). */

#define K_EQ 5.0                                  /* [C] / [B] */
#define K_DRIVE 0.1                               /* A -> B */
#define K_FWD 10.0                                /* B -> C */
#define K_BWD (K_FWD / K_EQ)                      /* C -> B */
#define EXPECTED_A 0.0
#define EXPECTED_B (1.0 / (1.0 + K_EQ))
#define EXPECTED_C (K_EQ / (1.0 + K_EQ))
#define KINETIC_MASS_TOL 1.0e-6                   /* double precision */
#define DAE_FORMULATION_AGREEMENT_TOL 1.0e-8      /* double precision */

/* state A = 1, B = C = 0 at 298 K and 101325 Pa */
static void set_abc_initial(tmc_state *state)
{
  set_var(state, 0, "A", 1.0);
  set_var(state, 0, "B", 0.0);
  set_var(state, 0, "C", 0.0);
  tmc_state_conditions(state)[0].temperature = 298.0;
  tmc_state_conditions(state)[0].pressure = 101325.0;
}

/* steps of UpdateStateParameters and Solve; failures name the system */
static void integrate(tmc_solver *solver, tmc_state *state, int steps, double dt, const char *what)
{
  for (int step = 0; step < steps; ++step)
  {
    tmc_solver_update_state_parameters(solver, state);
    tmc_result result = tmc_solver_solve(solver, dt, state);
    if (result.state != TMC_CONVERGED)
      tmc_test_fail(__FILE__, __LINE__, "%s solve failed at step %d: %s", what, step, tmc_solver_state_name(result.state));
  }
}

/* SolveKineticSystem: A -> B, B -> C, C -> B, three-stage Rosenbrock */
static bool solve_kinetic_system(double out[3])
{
  const char *species[] = { "A", "B", "C" };
  tmc_phase gas = gas_phase_of(species, 3);
  tmc_system system = tmc_system_create(&gas);
  tmc_reaction reactions[3] = {
    first_order_reaction("A", "B", K_DRIVE, &gas),
    first_order_reaction("B", "C", K_FWD, &gas),
    first_order_reaction("C", "B", K_BWD, &gas),
  };
  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, reactions, 3);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  const bool ok = state != NULL;
  if (state)
  {
    set_abc_initial(state);
    /* 200 s, 20 times the slowest time scale (1 / k = 10 s) */
    integrate(solver, state, 200, 1.0, "kinetic");
    out[0] = var(state, 0, "A");
    out[1] = var(state, 0, "B");
    out[2] = var(state, 0, "C");
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  for (int i = 0; i < 3; ++i)
    tmc_reaction_free(&reactions[i]);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return ok;
}

/* SolveConservativeConstraintSystem: A -> B, ConservativeEquilibriumConstraintModel */
static bool solve_conservative_constraint_system(double out[3])
{
  const char *species[] = { "A", "B", "C" };
  tmc_phase gas = gas_phase_of(species, 3);
  tmc_system system = tmc_system_create(&gas);
  tmc_reaction rxn_ab = first_order_reaction("A", "B", K_DRIVE, &gas);
  conservative_equilibrium_model eq_model;
  conservative_equilibrium_model_init(&eq_model, "A", "B", "C", K_EQ, 1.0);
  tmc_external_model model = conservative_equilibrium_model_table(&eq_model);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, &rxn_ab, 1);
  tmc_builder_add_external_model(b, &model);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  const bool ok = state != NULL;
  if (state)
  {
    set_abc_initial(state);
    integrate(solver, state, 200, 1.0, "conservative constraint");
    out[0] = var(state, 0, "A");
    out[1] = var(state, 0, "B");
    out[2] = var(state, 0, "C");
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_reaction_free(&rxn_ab);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return ok;
}

/* SolveSimpleConstraintSystem: A -> B, EquilibriumConstraintModel(B, C) */
static bool solve_simple_constraint_system(double out[3])
{
  const char *species[] = { "A", "B", "C" };
  tmc_phase gas = gas_phase_of(species, 3);
  tmc_system system = tmc_system_create(&gas);
  tmc_reaction rxn_ab = first_order_reaction("A", "B", K_DRIVE, &gas);
  equilibrium_constraint_model eq_model;
  equilibrium_constraint_model_init(&eq_model, "B", "C", K_EQ);
  tmc_external_model model = equilibrium_constraint_model_table(&eq_model);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, &rxn_ab, 1);
  tmc_builder_add_external_model(b, &model);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  const bool ok = state != NULL;
  if (state)
  {
    set_abc_initial(state);
    integrate(solver, state, 200, 1.0, "simple constraint");
    out[0] = var(state, 0, "A");
    out[1] = var(state, 0, "B");
    out[2] = var(state, 0, "C");
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_reaction_free(&rxn_ab);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return ok;
}

/* KineticVsConservativeConstraintConvergence */
TMC_TEST(kinetic_vs_conservative_constraint_convergence)
{
  double kin[3], con[3];
  TMC_ASSERT_TRUE(solve_kinetic_system(kin));
  TMC_ASSERT_TRUE(solve_conservative_constraint_system(con));

  /* both near the analytical steady state */
  TMC_EXPECT_NEAR(kin[0], EXPECTED_A, 1e-6);
  TMC_EXPECT_NEAR(con[0], EXPECTED_A, 1e-6);
  TMC_EXPECT_NEAR(kin[1], EXPECTED_B, 1e-3);
  TMC_EXPECT_NEAR(con[1], EXPECTED_B, 1e-3);
  TMC_EXPECT_NEAR(kin[2], EXPECTED_C, 1e-3);
  TMC_EXPECT_NEAR(con[2], EXPECTED_C, 1e-3);

  /* the two systems agree */
  TMC_EXPECT_NEAR(kin[0], con[0], 1e-6);
  TMC_EXPECT_NEAR(kin[1], con[1], 1e-3);
  TMC_EXPECT_NEAR(kin[2], con[2], 1e-3);

  /* mass conservation in both */
  TMC_EXPECT_NEAR(kin[0] + kin[1] + kin[2], 1.0, KINETIC_MASS_TOL);
  TMC_EXPECT_NEAR(con[0] + con[1] + con[2], 1.0, 1e-3);
  return TMC_TEST_PASS;
}

/* SimpleConstraintPreservesRatioNotMass: the ratio constraint creates C from B */
TMC_TEST(simple_constraint_preserves_ratio_not_mass)
{
  double kin[3], sim[3];
  TMC_ASSERT_TRUE(solve_kinetic_system(kin));
  TMC_ASSERT_TRUE(solve_simple_constraint_system(sim));

  /* both keep the equilibrium ratio */
  TMC_EXPECT_NEAR(kin[2] / kin[1], K_EQ, 1e-3);
  TMC_EXPECT_NEAR(sim[2] / sim[1], K_EQ, 1e-3);
  /* the kinetic system conserves mass */
  TMC_EXPECT_NEAR(kin[0] + kin[1] + kin[2], 1.0, KINETIC_MASS_TOL);
  /* the ratio constraint does not: A -> 0, B -> 1, C -> K_eq, total 1 + K_eq */
  TMC_EXPECT_NEAR(sim[0] + sim[1] + sim[2], 1.0 + K_EQ, 0.1);
  return TMC_TEST_PASS;
}

/* ---- BuiltInVsExternalModelConstraintStepByStep ----------------------------------------- */

TMC_TEST(built_in_vs_external_model_constraint_step_by_step)
{
  const char *species[] = { "A", "B", "C" };
  tmc_phase gas = gas_phase_of(species, 3);
  tmc_system system = tmc_system_create(&gas);
  tmc_reaction rxn_ab = first_order_reaction("A", "B", K_DRIVE, &gas);
  int status = TMC_TEST_PASS;

  /* the built-in constraint */
  tmc_constraint constraint = equilibrium_constraint("B_C_eq", "C", "B", "C", K_EQ);
  tmc_builder *builtin_builder = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(builtin_builder, &system);
  tmc_builder_set_reactions(builtin_builder, &rxn_ab, 1);
  tmc_builder_add_constraint(builtin_builder, &constraint);
  tmc_builder_set_reorder_state(builtin_builder, false);
  tmc_solver *builtin_solver = BUILD(builtin_builder);

  /* the external model's constraint */
  equilibrium_constraint_model eq_model;
  equilibrium_constraint_model_init(&eq_model, "B", "C", K_EQ);
  tmc_external_model model = equilibrium_constraint_model_table(&eq_model);
  tmc_builder *ext_builder = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(ext_builder, &system);
  tmc_builder_set_reactions(ext_builder, &rxn_ab, 1);
  tmc_builder_add_external_model(ext_builder, &model);
  tmc_builder_set_reorder_state(ext_builder, false);
  tmc_solver *ext_solver = BUILD(ext_builder);

  tmc_state *state_bi = builtin_solver ? tmc_solver_create_state(builtin_solver, 1) : NULL;
  tmc_state *state_ext = ext_solver ? tmc_solver_create_state(ext_solver, 1) : NULL;
  if (state_bi && state_ext)
  {
    set_abc_initial(state_bi);
    set_abc_initial(state_ext);
    const double dt = 0.5;
    for (int step = 0; step < 100; ++step)
    {
      tmc_solver_update_state_parameters(builtin_solver, state_bi);
      tmc_solver_update_state_parameters(ext_solver, state_ext);
      tmc_result res_bi = tmc_solver_solve(builtin_solver, dt, state_bi);
      tmc_result res_ext = tmc_solver_solve(ext_solver, dt, state_ext);
      if (res_bi.state != TMC_CONVERGED || res_ext.state != TMC_CONVERGED)
      {
        tmc_test_fail(__FILE__, __LINE__, "step %d: built-in %s, external %s", step, tmc_solver_state_name(res_bi.state), tmc_solver_state_name(res_ext.state));
        status = TMC_TEST_FAIL; /* ASSERT_EQ */
        break;
      }
      const double bi_A = var(state_bi, 0, "A"), bi_B = var(state_bi, 0, "B"), bi_C = var(state_bi, 0, "C");
      const double ext_A = var(state_ext, 0, "A"), ext_B = var(state_ext, 0, "B"), ext_C = var(state_ext, 0, "C");
      if (!(fabs(bi_A - ext_A) <= DAE_FORMULATION_AGREEMENT_TOL))
        tmc_test_fail(__FILE__, __LINE__, "A diverged at step %d: %.17g vs %.17g", step, bi_A, ext_A);
      if (!(fabs(bi_B - ext_B) <= DAE_FORMULATION_AGREEMENT_TOL))
        tmc_test_fail(__FILE__, __LINE__, "B diverged at step %d: %.17g vs %.17g", step, bi_B, ext_B);
      if (!(fabs(bi_C - ext_C) <= DAE_FORMULATION_AGREEMENT_TOL))
        tmc_test_fail(__FILE__, __LINE__, "C diverged at step %d: %.17g vs %.17g", step, bi_C, ext_C);
      if (!(fabs(K_EQ * bi_B - bi_C) <= 1e-6))
        tmc_test_fail(__FILE__, __LINE__, "built-in constraint violated at step %d", step);
      if (!(fabs(K_EQ * ext_B - ext_C) <= 1e-6))
        tmc_test_fail(__FILE__, __LINE__, "external constraint violated at step %d", step);
    }
  }
  else
    status = TMC_TEST_FAIL;
  tmc_state_free(state_bi);
  tmc_state_free(state_ext);
  tmc_solver_free(builtin_solver);
  tmc_solver_free(ext_solver);
  tmc_builder_free(builtin_builder);
  tmc_builder_free(ext_builder);
  tmc_constraint_free(&constraint);
  tmc_reaction_free(&rxn_ab);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return status;
}

/* ---- MultiEquilibriumKineticVsComposedConstraints ----------------------------------------
 * A -> B (driver), B <-> C (K1 = 5), B <-> D (K2 = 3). The constraint system
 * composes two EquilibriumConstraintModels (C and D rows) and a
 * MassConservationModel (B row); A is the only ODE variable. */
TMC_TEST(multi_equilibrium_kinetic_vs_composed_constraints)
{
  const double K1 = 5.0, K2 = 3.0, k_drive = 0.1, k_f1 = 10.0, k_b1 = k_f1 / K1, k_f2 = 8.0, k_b2 = k_f2 / K2;
  const double expected_B = 1.0 / (1.0 + K1 + K2), expected_C = K1 * expected_B, expected_D = K2 * expected_B;

  const char *species[] = { "A", "B", "C", "D" };
  tmc_phase gas = gas_phase_of(species, 4);
  tmc_system system = tmc_system_create(&gas);

  /* the kinetic system */
  tmc_reaction reactions[5] = {
    first_order_reaction("A", "B", k_drive, &gas),
    first_order_reaction("B", "C", k_f1, &gas),
    first_order_reaction("C", "B", k_b1, &gas),
    first_order_reaction("B", "D", k_f2, &gas),
    first_order_reaction("D", "B", k_b2, &gas),
  };
  tmc_builder *kin_builder = cpu_solver_builder(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(kin_builder, &system);
  tmc_builder_set_reactions(kin_builder, reactions, 5);
  tmc_builder_set_reorder_state(kin_builder, false);
  tmc_solver *kin_solver = BUILD(kin_builder);

  /* the constraint system: three composed external models */
  equilibrium_constraint_model eq_bc, eq_bd;
  equilibrium_constraint_model_init(&eq_bc, "B", "C", K1); /* C row: K1 [B] - [C] = 0 */
  equilibrium_constraint_model_init(&eq_bd, "B", "D", K2); /* D row: K2 [B] - [D] = 0 */
  mass_conservation_model conservation;
  mass_conservation_model_init(&conservation, "B", species, 4, 1.0); /* B row: A + B + C + D - 1 = 0 */
  tmc_external_model m_bc = equilibrium_constraint_model_table(&eq_bc);
  tmc_external_model m_bd = equilibrium_constraint_model_table(&eq_bd);
  tmc_external_model m_cons = mass_conservation_model_table(&conservation);
  tmc_builder *ext_builder = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(ext_builder, &system);
  tmc_builder_set_reactions(ext_builder, reactions, 1); /* A -> B only */
  tmc_builder_add_external_model(ext_builder, &m_bc);
  tmc_builder_add_external_model(ext_builder, &m_bd);
  tmc_builder_add_external_model(ext_builder, &m_cons);
  tmc_builder_set_reorder_state(ext_builder, false);
  tmc_solver *ext_solver = BUILD(ext_builder);

  tmc_state *state_kin = kin_solver ? tmc_solver_create_state(kin_solver, 1) : NULL;
  tmc_state *state_ext = ext_solver ? tmc_solver_create_state(ext_solver, 1) : NULL;
  if (state_kin && state_ext)
  {
    tmc_state *states[2] = { state_kin, state_ext };
    for (int s = 0; s < 2; ++s)
    {
      set_abc_initial(states[s]);
      set_var(states[s], 0, "D", 0.0);
    }
    for (int step = 0; step < 200; ++step)
    {
      tmc_solver_update_state_parameters(kin_solver, state_kin);
      tmc_solver_update_state_parameters(ext_solver, state_ext);
      tmc_result res_kin = tmc_solver_solve(kin_solver, 1.0, state_kin);
      tmc_result res_ext = tmc_solver_solve(ext_solver, 1.0, state_ext);
      if (res_kin.state != TMC_CONVERGED)
        tmc_test_fail(__FILE__, __LINE__, "kinetic failed step %d: %s", step, tmc_solver_state_name(res_kin.state));
      if (res_ext.state != TMC_CONVERGED)
        tmc_test_fail(__FILE__, __LINE__, "constraint failed step %d: %s", step, tmc_solver_state_name(res_ext.state));
    }
    const double kin_A = var(state_kin, 0, "A"), kin_B = var(state_kin, 0, "B"), kin_C = var(state_kin, 0, "C"), kin_D = var(state_kin, 0, "D");
    const double ext_A = var(state_ext, 0, "A"), ext_B = var(state_ext, 0, "B"), ext_C = var(state_ext, 0, "C"), ext_D = var(state_ext, 0, "D");

    /* both near the analytical steady state */
    TMC_EXPECT_NEAR(kin_A, 0.0, 1e-6);
    TMC_EXPECT_NEAR(ext_A, 0.0, 1e-6);
    TMC_EXPECT_NEAR(kin_B, expected_B, 1e-3);
    TMC_EXPECT_NEAR(ext_B, expected_B, 1e-3);
    TMC_EXPECT_NEAR(kin_C, expected_C, 1e-3);
    TMC_EXPECT_NEAR(ext_C, expected_C, 1e-3);
    TMC_EXPECT_NEAR(kin_D, expected_D, 1e-3);
    TMC_EXPECT_NEAR(ext_D, expected_D, 1e-3);
    /* the kinetic and constraint systems agree */
    TMC_EXPECT_NEAR(kin_B, ext_B, 1e-3);
    TMC_EXPECT_NEAR(kin_C, ext_C, 1e-3);
    TMC_EXPECT_NEAR(kin_D, ext_D, 1e-3);
    /* mass conservation */
    TMC_EXPECT_NEAR(kin_A + kin_B + kin_C + kin_D, 1.0, KINETIC_MASS_TOL);
    TMC_EXPECT_NEAR(ext_A + ext_B + ext_C + ext_D, 1.0, 1e-3);
    /* the equilibrium constraints hold */
    TMC_EXPECT_NEAR(K1 * ext_B - ext_C, 0.0, 1e-6);
    TMC_EXPECT_NEAR(K2 * ext_B - ext_D, 0.0, 1e-6);
  }
  tmc_state_free(state_kin);
  tmc_state_free(state_ext);
  tmc_solver_free(kin_solver);
  tmc_solver_free(ext_solver);
  tmc_builder_free(kin_builder);
  tmc_builder_free(ext_builder);
  for (int i = 0; i < 5; ++i)
    tmc_reaction_free(&reactions[i]);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

/* ---- ProcessJacobianElementInAlgebraicRowSurvivesFiltering --------------------------------
 * StubAerosolWithSolvent's process declares (A_AQ, S) in the algebraic A_AQ
 * row, which its constraint does not declare; the element stays in the
 * pattern, so the model's index binding and the solve work. */
TMC_TEST(process_jacobian_element_in_algebraic_row_survives_filtering)
{
  const char *species[] = { "A_GAS", "S" };
  tmc_phase gas = gas_phase_of(species, 2);
  tmc_system system = tmc_system_create(&gas);
  const double total = 1.0, k = 0.1;
  stub_aerosol_solvent aerosol;
  stub_aerosol_solvent_init(&aerosol, k, total);
  tmc_external_model model = stub_aerosol_solvent_model(&aerosol);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, NULL, 0);
  tmc_builder_add_external_model(b, &model);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  if (state)
  {
    int i_aq = tmc_state_variable_index(state, "AEROSOL.A_AQ");
    int i_gas = tmc_state_variable_index(state, "A_GAS"), i_s = tmc_state_variable_index(state, "S");
    TMC_EXPECT_TRUE(i_aq >= 0 && i_gas >= 0 && i_s >= 0);
    if (i_aq >= 0 && i_gas >= 0 && i_s >= 0)
    {
      /* A_AQ is algebraic, and the process element (A_AQ, S) is in the pattern */
      TMC_EXPECT_DOUBLE_EQ(tmc_solver_mass_matrix_diagonal(solver)[i_aq], 0.0);
      TMC_EXPECT_FALSE(tmc_sparse_is_zero(tmc_state_jacobian(state), i_aq, i_s));
      TMC_EXPECT_TRUE(aerosol.aq_s_flat >= 0);
      set_var(state, 0, "A_GAS", 0.9);
      set_var(state, 0, "AEROSOL.A_AQ", 0.1);
      set_var(state, 0, "S", 1.0);
      tmc_state_conditions(state)[0].temperature = 298.0;
      tmc_state_conditions(state)[0].pressure = 101325.0;
      const double dt = 10.0;
      for (int step = 0; step < 10; ++step)
      {
        tmc_result result = tmc_solver_solve(solver, dt, state);
        if (result.state != TMC_CONVERGED)
          tmc_test_fail(__FILE__, __LINE__, "step %d: %s", step, tmc_solver_state_name(result.state));
        double sum = var(state, 0, "A_GAS") + var(state, 0, "AEROSOL.A_AQ");
        if (!(fabs(sum - total) <= 1e-4))
          tmc_test_fail(__FILE__, __LINE__, "conservation violated at step %d: sum %.17g", step, sum);
      }
    }
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

/* ---- ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering -------------------------
 * A built-in equilibrium constraint makes the A_AQ row algebraic; the builder
 * clears built-in algebraic rows, and the external process element
 * (A_AQ, A_GAS) must come back for the model. */
TMC_TEST(process_jacobian_element_in_built_in_algebraic_row_survives_filtering)
{
  const char *species[] = { "A_GAS", "B" };
  tmc_phase gas = gas_phase_of(species, 2);
  tmc_system system = tmc_system_create(&gas);
  stub_aerosol_constraints aerosol;
  stub_aerosol_constraints_init(&aerosol, 0.1); /* its own constraint inactive */
  tmc_external_model model = stub_aerosol_constraints_model(&aerosol);
  /* built-in: K_eq [B] - [A_AQ] = 0, A_AQ algebraic */
  tmc_constraint constraint = equilibrium_constraint("B_AQ_eq", "AEROSOL.A_AQ", "B", "AEROSOL.A_AQ", 5.0);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, NULL, 0);
  tmc_builder_add_constraint(b, &constraint);
  tmc_builder_add_external_model(b, &model);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  if (state)
  {
    int i_gas = tmc_state_variable_index(state, "A_GAS"), i_aq = tmc_state_variable_index(state, "AEROSOL.A_AQ");
    TMC_EXPECT_TRUE(i_gas >= 0 && i_aq >= 0);
    if (i_gas >= 0 && i_aq >= 0)
    {
      /* A_AQ is algebraic through the built-in constraint */
      TMC_EXPECT_EQ(tmc_solver_num_constraints(solver), 1);
      TMC_EXPECT_DOUBLE_EQ(tmc_solver_mass_matrix_diagonal(solver)[i_aq], 0.0);
      /* the external process element in the algebraic row is in the pattern */
      TMC_EXPECT_FALSE(tmc_sparse_is_zero(tmc_state_jacobian(state), i_aq, i_gas));
    }
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_constraint_free(&constraint);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

/* ---- ExternalModelFiniteDifferenceJacobian.ProcessForcingJacobian ------------------------- */

TMC_TEST(fd_process_forcing_jacobian)
{
  stub_aerosol_constraints aerosol;
  stub_aerosol_constraints_init(&aerosol, 3.5);
  tmc_external_model model = stub_aerosol_constraints_model(&aerosol);
  const char *variable_names[] = { "A_GAS", "AEROSOL.A_AQ" };
  tmc_name_map *var_map = tmc_name_map_create(variable_names, 2);
  tmc_name_map *param_map = tmc_name_map_create(NULL, 0);
  const int num_species = 2;

  tmc_index_pair elements[8];
  int n_elements = model.nonzero_jacobian_elements(model.context, var_map, NULL, 0);
  TMC_ASSERT_TRUE(n_elements > 0 && n_elements <= 8);
  n_elements = model.nonzero_jacobian_elements(model.context, var_map, elements, n_elements);
  tmc_sparse *analytical_jac = jacobian_of(elements, n_elements, num_species, 2);
  model.finalize_process_setup(model.context, param_map, var_map, analytical_jac);

  tmc_dense *variables = tmc_dense_create(2, num_species, 1, 0.0);
  tmc_dense_set(variables, 0, 0, 0.8);
  tmc_dense_set(variables, 0, 1, 0.2);
  tmc_dense_set(variables, 1, 0, 0.3);
  tmc_dense_set(variables, 1, 1, 0.7);
  tmc_dense *params = tmc_dense_create(2, 0, 1, 0.0);
  model.subtract_jacobian_terms(model.context, params, variables, analytical_jac);

  model_forcing wrapper = { &model, params };
  tmc_dense *fd_jac = finite_difference_jacobian(process_forcing, &wrapper, variables, num_species, FD_PERTURBATION);
  jacobian_comparison comparison = compare_jacobian_to_finite_difference(analytical_jac, fd_jac, num_species, FD_ATOL, FD_RTOL);
  EXPECT_COMPARISON(comparison, "process Jacobian mismatch");
  jacobian_comparison sparsity = check_jacobian_sparsity_completeness(analytical_jac, fd_jac, num_species, FD_SPARSITY_THRESHOLD);
  EXPECT_COMPARISON(sparsity, "missing sparsity");

  tmc_dense_free(fd_jac);
  tmc_dense_free(params);
  tmc_dense_free(variables);
  tmc_sparse_free(analytical_jac);
  tmc_name_map_free(param_map);
  tmc_name_map_free(var_map);
  return TMC_TEST_PASS;
}

/* ---- ExternalModelFiniteDifferenceJacobian.ConstraintResidualJacobian --------------------- */

TMC_TEST(fd_constraint_residual_jacobian)
{
  stub_aerosol_constraints aerosol;
  stub_aerosol_constraints_init_with_total(&aerosol, 0.1, 1.0);
  tmc_external_model model = stub_aerosol_constraints_model(&aerosol);
  const char *variable_names[] = { "A_GAS", "AEROSOL.A_AQ" };
  tmc_name_map *var_map = tmc_name_map_create(variable_names, 2);
  tmc_name_map *param_map = tmc_name_map_create(NULL, 0);
  const int num_species = 2;

  tmc_index_pair elements[8];
  int n_elements = model.nonzero_constraint_jacobian_elements(model.context, var_map, NULL, 0);
  TMC_ASSERT_TRUE(n_elements > 0 && n_elements <= 8);
  n_elements = model.nonzero_constraint_jacobian_elements(model.context, var_map, elements, n_elements);
  tmc_sparse *analytical_jac = jacobian_of(elements, n_elements, num_species, 2);
  model.finalize_constraint_setup(model.context, param_map, var_map, analytical_jac);

  tmc_dense *variables = tmc_dense_create(2, num_species, 1, 0.0);
  tmc_dense_set(variables, 0, 0, 0.6);
  tmc_dense_set(variables, 0, 1, 0.4);
  tmc_dense_set(variables, 1, 0, 0.2);
  tmc_dense_set(variables, 1, 1, 0.8);
  tmc_dense *dummy_params = tmc_dense_create(2, 1, 1, 0.0);
  model.subtract_constraint_jacobian(model.context, dummy_params, variables, analytical_jac);

  model_forcing wrapper = { &model, dummy_params };
  tmc_dense *fd_jac = finite_difference_jacobian(constraint_residual, &wrapper, variables, num_species, FD_PERTURBATION);
  jacobian_comparison comparison = compare_jacobian_to_finite_difference(analytical_jac, fd_jac, num_species, FD_ATOL, FD_RTOL);
  EXPECT_COMPARISON(comparison, "constraint Jacobian mismatch");

  tmc_dense_free(fd_jac);
  tmc_dense_free(dummy_params);
  tmc_dense_free(variables);
  tmc_sparse_free(analytical_jac);
  tmc_name_map_free(param_map);
  tmc_name_map_free(var_map);
  return TMC_TEST_PASS;
}

/* ---- ExternalModelFiniteDifferenceJacobian.EquilibriumConstraintModelJacobian ------------- */

TMC_TEST(fd_equilibrium_constraint_model_jacobian)
{
  equilibrium_constraint_model eq_model;
  equilibrium_constraint_model_init(&eq_model, "A", "B", 2.5);
  tmc_external_model model = equilibrium_constraint_model_table(&eq_model);
  const char *variable_names[] = { "A", "B" };
  tmc_name_map *var_map = tmc_name_map_create(variable_names, 2);
  tmc_name_map *param_map = tmc_name_map_create(NULL, 0);
  const int num_species = 2;

  tmc_index_pair elements[8];
  int n_elements = model.nonzero_constraint_jacobian_elements(model.context, var_map, NULL, 0);
  TMC_ASSERT_TRUE(n_elements > 0 && n_elements <= 8);
  n_elements = model.nonzero_constraint_jacobian_elements(model.context, var_map, elements, n_elements);
  tmc_sparse *analytical_jac = jacobian_of(elements, n_elements, num_species, 1);
  model.finalize_constraint_setup(model.context, param_map, var_map, analytical_jac);

  tmc_dense *variables = tmc_dense_create(1, num_species, 1, 0.0);
  tmc_dense_set(variables, 0, 0, 3.0);
  tmc_dense_set(variables, 0, 1, 5.0);
  tmc_dense *dummy_params = tmc_dense_create(1, 1, 1, 0.0);
  model.subtract_constraint_jacobian(model.context, dummy_params, variables, analytical_jac);

  model_forcing wrapper = { &model, dummy_params };
  tmc_dense *fd_jac = finite_difference_jacobian(constraint_residual, &wrapper, variables, num_species, FD_PERTURBATION);
  jacobian_comparison comparison = compare_jacobian_to_finite_difference(analytical_jac, fd_jac, num_species, FD_ATOL, FD_RTOL);
  EXPECT_COMPARISON(comparison, "EquilibriumConstraintModel Jacobian mismatch");
  jacobian_comparison sparsity = check_jacobian_sparsity_completeness(analytical_jac, fd_jac, num_species, FD_SPARSITY_THRESHOLD);
  EXPECT_COMPARISON(sparsity, "missing sparsity");

  tmc_dense_free(fd_jac);
  tmc_dense_free(dummy_params);
  tmc_dense_free(variables);
  tmc_sparse_free(analytical_jac);
  tmc_name_map_free(param_map);
  tmc_name_map_free(var_map);
  return TMC_TEST_PASS;
}

/* ---- TemperatureDependentConstraintParameter ----------------------------------------------
 * A -> B (kinetic) and K_eq(T) [B] - [C] = 0 (algebraic), solved at 298.15 K
 * and at 350 K: [C] / [B] = K_eq(T) at each. */
TMC_TEST(temperature_dependent_constraint_parameter)
{
  const double k_drive = 0.1, K_EQ_REF = 2.0, DELTA_H_OVER_R = 3000.0, T_REF = 298.15;
  const char *species[] = { "A", "B", "C" };
  tmc_phase gas = gas_phase_of(species, 3);
  tmc_system system = tmc_system_create(&gas);
  tmc_reaction rxn_ab = first_order_reaction("A", "B", k_drive, &gas);
  temperature_equilibrium_model eq_model;
  temperature_equilibrium_model_init(&eq_model, "B", "C", K_EQ_REF, DELTA_H_OVER_R, T_REF);
  tmc_external_model model = temperature_equilibrium_model_table(&eq_model);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, &rxn_ab, 1);
  tmc_builder_add_external_model(b, &model);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *solver = BUILD(b);
  if (solver)
  {
    TMC_EXPECT_TRUE(tmc_name_map_lookup(tmc_solver_parameter_map(solver), "C_K_eq") >= 0);
    const double temperatures[2] = { T_REF, 350.0 };
    for (int t = 0; t < 2; ++t)
    {
      tmc_state *state = tmc_solver_create_state(solver, 1);
      set_var(state, 0, "A", 1.0);
      set_var(state, 0, "B", 0.0);
      set_var(state, 0, "C", 0.0);
      tmc_state_conditions(state)[0].temperature = temperatures[t];
      tmc_state_conditions(state)[0].pressure = 101325.0;
      integrate(solver, state, 200, 1.0, t == 0 ? "T=298" : "T=350");
      const double B_val = var(state, 0, "B"), C_val = var(state, 0, "C");
      const double K_eq_expected = t == 0 ? K_EQ_REF : K_EQ_REF * exp(DELTA_H_OVER_R * (1.0 / T_REF - 1.0 / 350.0));
      TMC_EXPECT_TRUE(B_val > 0.0);
      if (t == 1)
        TMC_EXPECT_TRUE(K_eq_expected > K_EQ_REF); /* K_eq increases with T for a positive delta_H */
      TMC_EXPECT_NEAR(C_val / B_val, K_eq_expected, 1e-4);
      tmc_state_free(state);
    }
  }
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_reaction_free(&rxn_ab);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

/* ---- ConstraintActiveWhenProcessOnlyModelAddedFirst ----------------------------------------
 * A -> B with the external constraint K_EQ [B] - [C] = 0, with a no-op
 * process model before or after the constraint model. The constraint must be
 * active either way, and the order must not change a bit of the result. */
static bool solve_with_model_order(bool process_only_first, double *B_out, double *C_out)
{
  const char *species[] = { "A", "B", "C" };
  tmc_phase gas = gas_phase_of(species, 3);
  tmc_system system = tmc_system_create(&gas);
  tmc_reaction rxn_ab = first_order_reaction("A", "B", K_DRIVE, &gas);
  equilibrium_constraint_model eq_model;
  equilibrium_constraint_model_init(&eq_model, "B", "C", K_EQ);
  tmc_external_model eq = equilibrium_constraint_model_table(&eq_model);
  tmc_external_model no_op = noop_process_model_table();

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, &rxn_ab, 1);
  tmc_builder_set_reorder_state(b, false);
  tmc_builder_add_external_model(b, process_only_first ? &no_op : &eq);
  tmc_builder_add_external_model(b, process_only_first ? &eq : &no_op);
  tmc_solver *solver = BUILD(b);
  tmc_state *state = solver ? tmc_solver_create_state(solver, 1) : NULL;
  const bool ok = state != NULL;
  if (state)
  {
    TMC_EXPECT_EQ(tmc_solver_num_constraints(solver), 1);
    int i_c = tmc_state_variable_index(state, "C");
    TMC_EXPECT_TRUE(i_c >= 0);
    if (i_c >= 0)
      TMC_EXPECT_DOUBLE_EQ(tmc_solver_mass_matrix_diagonal(solver)[i_c], 0.0);
    set_abc_initial(state);
    integrate(solver, state, 20, 1.0, process_only_first ? "process-only model first" : "constraint model first");
    *B_out = var(state, 0, "B");
    *C_out = var(state, 0, "C");
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_builder_free(b);
  tmc_reaction_free(&rxn_ab);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return ok;
}

TMC_TEST(constraint_active_when_process_only_model_added_first)
{
  const double residual_tol = 1.0e-6;
  double B_ref = 0.0, C_ref = 0.0, B_val = 0.0, C_val = 0.0;
  TMC_ASSERT_TRUE(solve_with_model_order(false, &B_ref, &C_ref));
  TMC_ASSERT_TRUE(solve_with_model_order(true, &B_val, &C_val));

  /* control: the constraint model first */
  TMC_EXPECT_TRUE(B_ref > 0.0);
  TMC_EXPECT_NEAR(K_EQ * B_ref - C_ref, 0.0, residual_tol);
  /* the no-op model adds no work: the same result, bit for bit */
  TMC_EXPECT_TRUE(B_val > 0.0);
  TMC_EXPECT_NEAR(K_EQ * B_val - C_val, 0.0, residual_tol);
  TMC_EXPECT_EQ(B_val, B_ref);
  TMC_EXPECT_EQ(C_val, C_ref);
  return TMC_TEST_PASS;
}

/* ---- BuilderCanBuildMultipleSolvers ----------------------------------------------------------
 * Two builds from one builder give two solvers with the same behavior: the
 * builder keeps its external models and built-in constraints. */
static void solve_abcd(tmc_solver *solver, double values[4])
{
  tmc_state *state = tmc_solver_create_state(solver, 1);
  TMC_EXPECT_EQ(tmc_solver_num_constraints(solver), 2);
  set_abc_initial(state);
  set_var(state, 0, "D", 2.0);
  integrate(solver, state, 20, 1.0, "two-constraint");
  const char *names[4] = { "A", "B", "C", "D" };
  for (int i = 0; i < 4; ++i)
    values[i] = var(state, 0, names[i]);
  tmc_state_free(state);
}

TMC_TEST(builder_can_build_multiple_solvers)
{
  const char *species[] = { "A", "B", "C", "D" };
  tmc_phase gas = gas_phase_of(species, 4);
  tmc_system system = tmc_system_create(&gas);
  tmc_reaction rxn_ab = first_order_reaction("A", "B", K_DRIVE, &gas);
  /* built-in: 2 [A] - [D] = 0, D algebraic */
  tmc_constraint constraint = equilibrium_constraint("A_D_eq", "D", "A", "D", 2.0);
  /* external: K_EQ [B] - [C] = 0, C algebraic */
  equilibrium_constraint_model eq_model;
  equilibrium_constraint_model_init(&eq_model, "B", "C", K_EQ);
  tmc_external_model model = equilibrium_constraint_model_table(&eq_model);

  tmc_builder *b = cpu_solver_builder(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, &rxn_ab, 1);
  tmc_builder_add_constraint(b, &constraint);
  tmc_builder_set_reorder_state(b, false);
  tmc_builder_add_external_model(b, &model);
  tmc_solver *first_solver = BUILD(b);
  tmc_solver *second_solver = BUILD(b);
  if (first_solver && second_solver)
  {
    double first[4], second[4];
    solve_abcd(first_solver, first);
    solve_abcd(second_solver, second);
    /* both constraints hold in the first solver */
    const double residual_tol = 1.0e-6;
    TMC_EXPECT_NEAR(2.0 * first[0] - first[3], 0.0, residual_tol);
    TMC_EXPECT_NEAR(K_EQ * first[1] - first[2], 0.0, residual_tol);
    /* one configuration, so the results are bit-for-bit equal */
    for (int i = 0; i < 4; ++i)
      TMC_EXPECT_EQ(first[i], second[i]);
  }
  tmc_solver_free(first_solver);
  tmc_solver_free(second_solver);
  tmc_builder_free(b);
  tmc_constraint_free(&constraint);
  tmc_reaction_free(&rxn_ab);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("external_model_constraints:ExternalModelConstraints.AddExternalModelWithConstraints", add_external_model_with_constraints),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.AddExternalModelProcessOnly", add_external_model_process_only),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.DAESolveEnforcesConservation", dae_solve_enforces_conservation),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.AlgebraicRowHoldsOnlyConstraintTerms", algebraic_row_holds_only_constraint_terms),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.CombinedBuiltInAndExternalConstraints", combined_built_in_and_external_constraints),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.AddExternalModelOnlyStandardRosenbrock", add_external_model_only_standard_rosenbrock),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.AddExternalModelConstraintsOnly", add_external_model_constraints_only),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.MultiGridCell", multi_grid_cell),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.KineticVsConservativeConstraintConvergence", kinetic_vs_conservative_constraint_convergence),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.SimpleConstraintPreservesRatioNotMass", simple_constraint_preserves_ratio_not_mass),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.BuiltInVsExternalModelConstraintStepByStep", built_in_vs_external_model_constraint_step_by_step),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.MultiEquilibriumKineticVsComposedConstraints", multi_equilibrium_kinetic_vs_composed_constraints),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.ProcessJacobianElementInAlgebraicRowSurvivesFiltering", process_jacobian_element_in_algebraic_row_survives_filtering),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering", process_jacobian_element_in_built_in_algebraic_row_survives_filtering),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.TemperatureDependentConstraintParameter", temperature_dependent_constraint_parameter),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.ConstraintActiveWhenProcessOnlyModelAddedFirst", constraint_active_when_process_only_model_added_first),
  TMC_CASE("external_model_constraints:ExternalModelConstraints.BuilderCanBuildMultipleSolvers", builder_can_build_multiple_solvers),
  TMC_CASE("external_model_constraints:ExternalModelFiniteDifferenceJacobian.ProcessForcingJacobian", fd_process_forcing_jacobian),
  TMC_CASE("external_model_constraints:ExternalModelFiniteDifferenceJacobian.ConstraintResidualJacobian", fd_constraint_residual_jacobian),
  TMC_CASE("external_model_constraints:ExternalModelFiniteDifferenceJacobian.EquilibriumConstraintModelJacobian", fd_equilibrium_constraint_model_jacobian),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_external_model_constraints(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
