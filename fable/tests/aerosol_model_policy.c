/* MICM test/integration/aerosol_model_policy.hpp; see aerosol_model_policy.h. */
#include "aerosol_model_policy.h"

#include "harness.h"
#include "stub_aerosol.h"

/* GoogleTest's EXPECT macros, with the builder's name in the message */
#define P_EXPECT_TRUE(cond)                                                          \
  do                                                                                 \
  {                                                                                  \
    if (!(cond))                                                                     \
      tmc_test_fail(__FILE__, __LINE__, "[%s] expected true: %s", p->name, #cond);   \
  } while (0)
#define P_EXPECT_INT_EQ(a, b)                                                        \
  do                                                                                 \
  {                                                                                  \
    long p_a_ = (long)(a), p_b_ = (long)(b);                                         \
    if (p_a_ != p_b_)                                                                \
      tmc_test_fail(__FILE__, __LINE__, "[%s] expected %s == %s: %ld vs %ld", p->name, #a, #b, p_a_, p_b_); \
  } while (0)
#define P_EXPECT_EQ(a, b)                                                            \
  do                                                                                 \
  {                                                                                  \
    double p_a_ = (a), p_b_ = (b);                                                   \
    if (!(p_a_ == p_b_))                                                             \
      tmc_test_fail(__FILE__, __LINE__, "[%s] expected %s == %s: %.17g vs %.17g", p->name, #a, #b, p_a_, p_b_); \
  } while (0)
#define P_EXPECT_DOUBLE_EQ(a, b)                                                     \
  do                                                                                 \
  {                                                                                  \
    double p_a_ = (a), p_b_ = (b);                                                   \
    if (!tmc_test_double_eq(p_a_, p_b_))                                             \
      tmc_test_fail(__FILE__, __LINE__, "[%s] expected %s (%.17g) == %s (%.17g)", p->name, #a, p_a_, #b, p_b_); \
  } while (0)
#define P_EXPECT_NEAR(a, b, tol)                                                     \
  do                                                                                 \
  {                                                                                  \
    double p_a_ = (a), p_b_ = (b), p_t_ = (tol);                                     \
    if (!(fabs(p_a_ - p_b_) <= p_t_))                                                \
      tmc_test_fail(__FILE__, __LINE__, "[%s] expected |%s - %s| <= %g: %.17g vs %.17g", p->name, #a, #b, p_t_, p_a_, p_b_); \
  } while (0)

/* ---- CreateSystemWithStubAerosolModels ---------------------------------------- */

typedef struct fixture
{
  tmc_species foo, bar, baz, qux;
  tmc_phase gas, quux, corge; /* the phase map { GAS, QUUX, CORGE } */
  tmc_system system;
  stub_aerosol_1 aerosol_1;
  stub_aerosol_2 aerosol_2;
} fixture;

static void fixture_init(fixture *f)
{
  f->foo = tmc_species_create("FO2"); /* species that can partition to condensed phase */
  f->bar = tmc_species_create("BAR"); /* gas-phase only species */
  f->baz = tmc_species_create("BAZ"); /* condensed-phase only species */
  f->qux = tmc_species_create("QUX"); /* condensed-phase only species */

  f->gas = tmc_phase_create("GAS");
  tmc_phase_add_species(&f->gas, &f->foo);
  tmc_phase_add_species(&f->gas, &f->bar);
  f->quux = tmc_phase_create("QUUX"); /* condensed aerosol or cloud phase */
  tmc_phase_add_species(&f->quux, &f->baz);
  tmc_phase_add_species(&f->quux, &f->qux);
  f->corge = tmc_phase_create("CORGE"); /* another condensed aerosol or cloud phase */
  tmc_phase_add_species(&f->corge, &f->foo);
  tmc_phase_add_species(&f->corge, &f->baz);
  tmc_phase_add_species(&f->corge, &f->qux);

  stub1_rate_constants rate_constants = { STUB1_RATE_CONSTANT_FO2_CORGE, STUB1_RATE_CONSTANT_BAZ_QUUX };
  tmc_phase model_phases[2] = { f->quux, f->corge };
  stub_aerosol_1_init(&f->aerosol_1, "STUB1", model_phases, 2, rate_constants);
  stub_aerosol_2_init(&f->aerosol_2, "STUB2", model_phases, 2);

  /* the system holds the gas phase; the models add their state through the builder */
  f->system = tmc_system_create(&f->gas);
}

static void fixture_free(fixture *f)
{
  tmc_system_free(&f->system);
  tmc_phase_free(&f->gas);
  tmc_phase_free(&f->quux);
  tmc_phase_free(&f->corge);
}

static tmc_builder *policy_builder(const aerosol_policy *p)
{
  tmc_builder *b = NULL;
  switch (p->solver)
  {
    case AEROSOL_ROSENBROCK_2STAGE:
      b = tmc_builder_create_rosenbrock(tmc_rosenbrock_two_stage());
      break;
    case AEROSOL_ROSENBROCK_3STAGE:
      b = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
      break;
    case AEROSOL_ROSENBROCK_4STAGE:
      b = tmc_builder_create_rosenbrock(tmc_rosenbrock_four_stage());
      break;
    case AEROSOL_ROSENBROCK_4STAGE_DAE:
      b = tmc_builder_create_rosenbrock(tmc_rosenbrock_four_stage_dae());
      break;
    case AEROSOL_ROSENBROCK_6STAGE_DAE:
      b = tmc_builder_create_rosenbrock(tmc_rosenbrock_six_stage_dae());
      break;
    case AEROSOL_BACKWARD_EULER:
      b = tmc_builder_create_backward_euler(tmc_backward_euler_defaults());
      break;
  }
  tmc_builder_set_width(b, p->width);
  tmc_builder_set_lu_algorithm(b, p->lu);
  return b;
}

/* builder.SetSystem(system).AddExternalModel(aerosol_1).AddExternalModel(aerosol_2)
 *        .SetIgnoreUnusedSpecies(true).Build() */
static tmc_solver *build_solver(const aerosol_policy *p, fixture *f)
{
  tmc_builder *b = policy_builder(p);
  tmc_builder_set_system(b, &f->system);
  tmc_external_model aerosol_1 = stub_aerosol_1_model(&f->aerosol_1);
  tmc_external_model aerosol_2 = stub_aerosol_2_model(&f->aerosol_2);
  tmc_builder_add_external_model(b, &aerosol_1);
  tmc_builder_add_external_model(b, &aerosol_2);
  tmc_builder_set_ignore_unused_species(b, true);
  tmc_error err;
  tmc_solver *s = tmc_builder_build(b, &err);
  tmc_builder_free(b);
  if (!s)
    tmc_test_fail(__FILE__, __LINE__, "[%s] build failed: %s", p->name, err.message);
  return s;
}

/* state[name] for a cell (NaN when the name is unknown, which fails any check) */
static double value(const tmc_state *st, int cell, const char *name)
{
  int i = tmc_state_variable_index(st, name);
  return i < 0 ? NAN : tmc_state_get(st, cell, i);
}

/* state[index] for a cell */
static double value_at(const tmc_state *st, int cell, int index)
{
  return index < 0 ? NAN : tmc_state_get(st, cell, index);
}

/* state[name] = value (every cell) */
static void set_value(const aerosol_policy *p, tmc_state *st, const char *name, double v)
{
  tmc_error err;
  if (tmc_state_set_concentration(st, name, v, &err) != TMC_OK)
    tmc_test_fail(__FILE__, __LINE__, "[%s] setting %s: %s", p->name, name, err.message);
}

/* state[name] = std::vector{ values } */
static void set_values(const aerosol_policy *p, tmc_state *st, const char *name, const double *v, int n)
{
  tmc_error err;
  if (tmc_state_set_concentrations(st, name, v, n, &err) != TMC_OK)
    tmc_test_fail(__FILE__, __LINE__, "[%s] setting %s: %s", p->name, name, err.message);
}

/* A copy of the state's Jacobian in the policy's element numbering. */
static tmc_sparse *jacobian_in_policy_ordering(const aerosol_policy *p, const tmc_sparse *jacobian)
{
  if (jacobian->ordering == p->ordering)
    return tmc_sparse_copy(jacobian);
  tmc_sparse_builder *sb = tmc_sparse_builder_create(jacobian->n);
  for (int e = 0; e < jacobian->nnz; ++e)
    tmc_sparse_builder_with_element(sb, jacobian->row_of[e], jacobian->col_of[e], NULL);
  tmc_sparse_builder_set_number_of_blocks(sb, jacobian->n_blocks);
  tmc_sparse *m = tmc_sparse_build(sb, p->ordering, jacobian->L, false);
  tmc_sparse_builder_free(sb);
  for (int block = 0; block < jacobian->n_blocks; ++block)
    for (int e = 0; e < jacobian->nnz; ++e)
      tmc_sparse_set(m, block, jacobian->row_of[e], jacobian->col_of[e], tmc_sparse_get(jacobian, block, jacobian->row_of[e], jacobian->col_of[e]));
  return m;
}

/* ---- TestStateIncludesStubAerosolModel ---------------------------------------- */

void test_state_includes_stub_aerosol_model(const aerosol_policy *p)
{
  fixture f;
  fixture_init(&f);
  tmc_solver *solver = build_solver(p, &f);
  if (!solver)
  {
    fixture_free(&f);
    return;
  }
  tmc_state *state = tmc_solver_create_state(solver, 1);

  /* the size and labels match expectations */
  int n1_variables, n1_parameters, n2_variables, n2_parameters;
  P_EXPECT_INT_EQ(f.aerosol_1.n_phases, 2);
  P_EXPECT_INT_EQ(f.aerosol_2.n_phases, 2);
  stub_aerosol_1_state_size(&f.aerosol_1, &n1_variables, &n1_parameters);
  stub_aerosol_2_state_size(&f.aerosol_2, &n2_variables, &n2_parameters);
  const int n1_names = stub_aerosol_1_state_variable_names(&f.aerosol_1, NULL, 0);
  const int n2_names = stub_aerosol_2_state_variable_names(&f.aerosol_2, NULL, 0);
  const int gas_unique_names = tmc_phase_species_names(&f.system.gas_phase, true, NULL, 0);
  P_EXPECT_INT_EQ(tmc_name_map_size(tmc_solver_variable_map(solver)), gas_unique_names + n1_names + n2_names);
  P_EXPECT_INT_EQ(tmc_state_num_variables(state), gas_unique_names + n1_names + n2_names);
  P_EXPECT_INT_EQ(n1_variables, n1_names);
  P_EXPECT_INT_EQ(n2_variables, n2_names);
  P_EXPECT_INT_EQ(n1_parameters, stub_aerosol_1_state_parameter_names(&f.aerosol_1, NULL, 0));
  P_EXPECT_INT_EQ(n2_parameters, stub_aerosol_2_state_parameter_names(&f.aerosol_2, NULL, 0));

  /* the full list of expected variable names: the gas phase's species names,
   * then each model's state variable names */
  char expected[64][TMC_NAME_MAX];
  int n_expected = tmc_phase_species_names(&f.system.gas_phase, false, expected, 64);
  n_expected += stub_aerosol_1_state_variable_names(&f.aerosol_1, expected + n_expected, 64 - n_expected);
  n_expected += stub_aerosol_2_state_variable_names(&f.aerosol_2, expected + n_expected, 64 - n_expected);
  P_EXPECT_INT_EQ(n_expected, 2 + 7 + 13);
  for (int i = 0; i < n_expected && i < 64; ++i)
    if (tmc_state_variable_index(state, expected[i]) < 0)
      tmc_test_fail(__FILE__, __LINE__, "[%s] %s is not in the variable map", p->name, expected[i]);

  /* a few specific names, for the naming convention */
  P_EXPECT_TRUE(tmc_state_variable_index(state, "FO2") >= 0);
  P_EXPECT_TRUE(tmc_state_variable_index(state, "STUB1.MODE1.QUUX.QUX") >= 0);
  P_EXPECT_TRUE(tmc_state_variable_index(state, "STUB1.MODE2.CORGE.FO2") >= 0);
  P_EXPECT_TRUE(tmc_state_variable_index(state, "STUB2.MODE1.NUMBER") >= 0);
  P_EXPECT_TRUE(tmc_state_variable_index(state, "STUB2.MODE2.CORGE.BAZ") >= 0);

  tmc_state_free(state);
  tmc_solver_free(solver);
  fixture_free(&f);
}

/* ---- TestUpdateStateWithStubAerosolModel --------------------------------------- */

void test_update_state_with_stub_aerosol_model(const aerosol_policy *p)
{
  fixture f;
  fixture_init(&f);
  tmc_solver *solver = build_solver(p, &f);
  if (!solver)
  {
    fixture_free(&f);
    return;
  }
  tmc_state *state = tmc_solver_create_state(solver, 1);
  const tmc_name_map *variables = tmc_solver_variable_map(solver);
  char name[TMC_NAME_MAX];

  /* gas-phase species by name, and by species object */
  set_value(p, state, "FO2", 1.34);
  tmc_species bar = tmc_species_create("BAR");
  set_value(p, state, bar.name, 2.53);
  /* condensed-phase species of the first model by unique name */
  set_value(p, state, "STUB1.MODE1.QUUX.BAZ", 0.75);
  set_value(p, state, "STUB1.MODE2.CORGE.FO2", 1.23);
  /* of the second model by the model's name for a species object, and by index */
  tmc_species qux = tmc_species_create("QUX");
  stub_aerosol_2_species(&f.aerosol_2, 2, &f.corge, &qux, name);
  set_value(p, state, name, 0.42);
  int corge_baz_index = tmc_state_variable_index(state, "STUB2.MODE3.CORGE.BAZ");
  if (corge_baz_index < 0)
  {
    tmc_test_fail(__FILE__, __LINE__, "[%s] STUB2.MODE3.CORGE.BAZ is not in the variable map", p->name);
    goto done;
  }
  tmc_state_set(state, 0, corge_baz_index, 0.33);
  /* model-specific variables by unique name, and through the model */
  set_value(p, state, "STUB2.MODE1.NUMBER", 1000.0);
  stub_aerosol_2_number(&f.aerosol_2, 1, name);
  set_value(p, state, name, 500.0);

  /* every value by name */
  P_EXPECT_DOUBLE_EQ(value(state, 0, "FO2"), 1.34);
  P_EXPECT_DOUBLE_EQ(value(state, 0, "BAR"), 2.53);
  P_EXPECT_DOUBLE_EQ(value(state, 0, "STUB1.MODE1.QUUX.BAZ"), 0.75);
  P_EXPECT_DOUBLE_EQ(value(state, 0, "STUB1.MODE2.CORGE.FO2"), 1.23);
  P_EXPECT_DOUBLE_EQ(value(state, 0, "STUB2.MODE3.CORGE.QUX"), 0.42);
  P_EXPECT_DOUBLE_EQ(value(state, 0, "STUB2.MODE3.CORGE.BAZ"), 0.33);
  P_EXPECT_DOUBLE_EQ(value(state, 0, "STUB2.MODE1.NUMBER"), 1000.0);
  P_EXPECT_DOUBLE_EQ(value(state, 0, "STUB2.MODE2.NUMBER"), 500.0);

  /* again by the indices of the variable map */
  P_EXPECT_DOUBLE_EQ(value_at(state, 0, tmc_name_map_lookup(variables, "FO2")), 1.34);
  P_EXPECT_DOUBLE_EQ(value_at(state, 0, tmc_name_map_lookup(variables, "BAR")), 2.53);
  P_EXPECT_DOUBLE_EQ(value_at(state, 0, tmc_name_map_lookup(variables, "STUB1.MODE1.QUUX.BAZ")), 0.75);
  P_EXPECT_DOUBLE_EQ(value_at(state, 0, tmc_name_map_lookup(variables, "STUB1.MODE2.CORGE.FO2")), 1.23);
  P_EXPECT_DOUBLE_EQ(value_at(state, 0, tmc_name_map_lookup(variables, "STUB2.MODE3.CORGE.QUX")), 0.42);
  P_EXPECT_DOUBLE_EQ(value_at(state, 0, tmc_name_map_lookup(variables, "STUB2.MODE3.CORGE.BAZ")), 0.33);
  P_EXPECT_DOUBLE_EQ(value_at(state, 0, tmc_name_map_lookup(variables, "STUB2.MODE1.NUMBER")), 1000.0);
  P_EXPECT_DOUBLE_EQ(value_at(state, 0, tmc_name_map_lookup(variables, "STUB2.MODE2.NUMBER")), 500.0);

  /* and by species objects and the models' names */
  tmc_species fo2 = tmc_species_create("FO2"), baz = tmc_species_create("BAZ");
  P_EXPECT_DOUBLE_EQ(value(state, 0, fo2.name), 1.34);
  P_EXPECT_DOUBLE_EQ(value(state, 0, bar.name), 2.53);
  stub_aerosol_1_species(&f.aerosol_1, 0, &f.quux, &baz, name);
  P_EXPECT_DOUBLE_EQ(value(state, 0, name), 0.75);
  stub_aerosol_1_species(&f.aerosol_1, 1, &f.corge, &fo2, name);
  P_EXPECT_DOUBLE_EQ(value(state, 0, name), 1.23);
  stub_aerosol_2_species(&f.aerosol_2, 2, &f.corge, &qux, name);
  P_EXPECT_DOUBLE_EQ(value(state, 0, name), 0.42);
  stub_aerosol_2_species(&f.aerosol_2, 2, &f.corge, &baz, name);
  P_EXPECT_DOUBLE_EQ(value(state, 0, name), 0.33);
  stub_aerosol_2_number(&f.aerosol_2, 0, name);
  P_EXPECT_DOUBLE_EQ(value(state, 0, name), 1000.0);
  stub_aerosol_2_number(&f.aerosol_2, 1, name);
  P_EXPECT_DOUBLE_EQ(value(state, 0, name), 500.0);

done:
  tmc_state_free(state);
  tmc_solver_free(solver);
  fixture_free(&f);
}

/* ---- TestUpdateMultiCellStateWithStubAerosolModel ------------------------------ */

/* EXPECT_EQ(state[...], std::vector{ ... }): every cell, exactly */
#define P_EXPECT_CELLS(get, expected)                                                \
  do                                                                                 \
  {                                                                                  \
    for (int cell_ = 0; cell_ < 3; ++cell_)                                          \
    {                                                                                \
      double got_ = get;                                                             \
      if (!(got_ == (expected)[cell_]))                                              \
        tmc_test_fail(__FILE__, __LINE__, "[%s] cell %d: expected %s == %.17g, got %.17g", p->name, cell_, #get, (expected)[cell_], got_); \
    }                                                                                \
  } while (0)

void test_update_multi_cell_state_with_stub_aerosol_model(const aerosol_policy *p)
{
  fixture f;
  fixture_init(&f);
  tmc_solver *solver = build_solver(p, &f);
  if (!solver)
  {
    fixture_free(&f);
    return;
  }
  const int num_cells = 3;
  tmc_state *state = tmc_solver_create_state(solver, num_cells);
  const tmc_name_map *variables = tmc_solver_variable_map(solver);
  char name[TMC_NAME_MAX];

  const double fo2[] = { 1.34, 1.35, 1.36 };
  const double bar_values[] = { 2.53, 2.54, 2.55 };
  const double baz_mode1[] = { 0.75, 0.76, 0.77 };
  const double fo2_mode2[] = { 1.23, 1.24, 1.25 };
  const double qux_mode3[] = { 0.42, 0.43, 0.44 };
  const double baz_mode3[] = { 0.33, 0.34, 0.35 };
  const double number1[] = { 1000.0, 1001.0, 1002.0 };
  const double number2[] = { 500.0, 501.0, 502.0 };

  set_values(p, state, "FO2", fo2, num_cells);
  tmc_species bar = tmc_species_create("BAR");
  set_values(p, state, bar.name, bar_values, num_cells);
  set_values(p, state, "STUB1.MODE1.QUUX.BAZ", baz_mode1, num_cells);
  set_values(p, state, "STUB1.MODE2.CORGE.FO2", fo2_mode2, num_cells);
  tmc_species qux = tmc_species_create("QUX");
  stub_aerosol_2_species(&f.aerosol_2, 2, &f.corge, &qux, name);
  set_values(p, state, name, qux_mode3, num_cells);
  int corge_baz_index = tmc_state_variable_index(state, "STUB2.MODE3.CORGE.BAZ");
  if (corge_baz_index < 0)
  {
    tmc_test_fail(__FILE__, __LINE__, "[%s] STUB2.MODE3.CORGE.BAZ is not in the variable map", p->name);
    goto done;
  }
  for (int cell = 0; cell < num_cells; ++cell)
    tmc_state_set(state, cell, corge_baz_index, baz_mode3[cell]);
  set_values(p, state, "STUB2.MODE1.NUMBER", number1, num_cells);
  stub_aerosol_2_number(&f.aerosol_2, 1, name);
  set_values(p, state, name, number2, num_cells);

  /* by name */
  P_EXPECT_CELLS(value(state, cell_, "FO2"), fo2);
  P_EXPECT_CELLS(value(state, cell_, "BAR"), bar_values);
  P_EXPECT_CELLS(value(state, cell_, "STUB1.MODE1.QUUX.BAZ"), baz_mode1);
  P_EXPECT_CELLS(value(state, cell_, "STUB1.MODE2.CORGE.FO2"), fo2_mode2);
  P_EXPECT_CELLS(value(state, cell_, "STUB2.MODE3.CORGE.QUX"), qux_mode3);
  P_EXPECT_CELLS(value(state, cell_, "STUB2.MODE3.CORGE.BAZ"), baz_mode3);
  P_EXPECT_CELLS(value(state, cell_, "STUB2.MODE1.NUMBER"), number1);
  P_EXPECT_CELLS(value(state, cell_, "STUB2.MODE2.NUMBER"), number2);

  /* by the indices of the variable map */
  P_EXPECT_CELLS(value_at(state, cell_, tmc_name_map_lookup(variables, "FO2")), fo2);
  P_EXPECT_CELLS(value_at(state, cell_, tmc_name_map_lookup(variables, "BAR")), bar_values);
  P_EXPECT_CELLS(value_at(state, cell_, tmc_name_map_lookup(variables, "STUB1.MODE1.QUUX.BAZ")), baz_mode1);
  P_EXPECT_CELLS(value_at(state, cell_, tmc_name_map_lookup(variables, "STUB1.MODE2.CORGE.FO2")), fo2_mode2);
  P_EXPECT_CELLS(value_at(state, cell_, tmc_name_map_lookup(variables, "STUB2.MODE3.CORGE.QUX")), qux_mode3);
  P_EXPECT_CELLS(value_at(state, cell_, tmc_name_map_lookup(variables, "STUB2.MODE3.CORGE.BAZ")), baz_mode3);
  P_EXPECT_CELLS(value_at(state, cell_, tmc_name_map_lookup(variables, "STUB2.MODE1.NUMBER")), number1);
  P_EXPECT_CELLS(value_at(state, cell_, tmc_name_map_lookup(variables, "STUB2.MODE2.NUMBER")), number2);

  /* by species objects and the models' names */
  tmc_species fo2_species = tmc_species_create("FO2"), baz = tmc_species_create("BAZ");
  char baz1[TMC_NAME_MAX], fo2_2[TMC_NAME_MAX], qux3[TMC_NAME_MAX], baz3[TMC_NAME_MAX], n1[TMC_NAME_MAX], n2[TMC_NAME_MAX];
  stub_aerosol_1_species(&f.aerosol_1, 0, &f.quux, &baz, baz1);
  stub_aerosol_1_species(&f.aerosol_1, 1, &f.corge, &fo2_species, fo2_2);
  stub_aerosol_2_species(&f.aerosol_2, 2, &f.corge, &qux, qux3);
  stub_aerosol_2_species(&f.aerosol_2, 2, &f.corge, &baz, baz3);
  stub_aerosol_2_number(&f.aerosol_2, 0, n1);
  stub_aerosol_2_number(&f.aerosol_2, 1, n2);
  P_EXPECT_CELLS(value(state, cell_, fo2_species.name), fo2);
  P_EXPECT_CELLS(value(state, cell_, bar.name), bar_values);
  P_EXPECT_CELLS(value(state, cell_, baz1), baz_mode1);
  P_EXPECT_CELLS(value(state, cell_, fo2_2), fo2_mode2);
  P_EXPECT_CELLS(value(state, cell_, qux3), qux_mode3);
  P_EXPECT_CELLS(value(state, cell_, baz3), baz_mode3);
  P_EXPECT_CELLS(value(state, cell_, n1), number1);
  P_EXPECT_CELLS(value(state, cell_, n2), number2);

done:
  tmc_state_free(state);
  tmc_solver_free(solver);
  fixture_free(&f);
}

/* ---- the state of the single-cell forcing and Jacobian checks ------------------ */

static void set_forcing_state(const aerosol_policy *p, tmc_state *state)
{
  set_value(p, state, "FO2", 1.0);
  set_value(p, state, "BAR", 2.0);
  set_value(p, state, "STUB1.MODE1.QUUX.BAZ", 0.5);
  set_value(p, state, "STUB1.MODE2.CORGE.FO2", 0.8);
  set_value(p, state, "STUB2.MODE3.CORGE.QUX", 0.3);
  set_value(p, state, "STUB2.MODE3.CORGE.BAZ", 0.2);
  set_value(p, state, "STUB2.MODE1.NUMBER", 1000.0);
  set_value(p, state, "STUB2.MODE2.NUMBER", 500.0);
}

/* ---- TestSingleCellForcingWithStubAerosolModel --------------------------------- */

void test_single_cell_forcing_with_stub_aerosol_model(const aerosol_policy *p)
{
  fixture f;
  fixture_init(&f);
  tmc_solver *solver = build_solver(p, &f);
  if (!solver)
  {
    fixture_free(&f);
    return;
  }
  tmc_state *state = tmc_solver_create_state(solver, 1);
  set_forcing_state(p, state);
  tmc_dense *forcing_1 = NULL;

  /* the forcing of the first model, called directly on the test's own copy of
   * the model (the solver holds the other) */
  stub_aerosol_1 aerosol_1 = f.aerosol_1;
  tmc_external_model model = stub_aerosol_1_model(&aerosol_1);
  model.finalize_process_setup(model.context, tmc_solver_parameter_map(solver), tmc_solver_variable_map(solver), tmc_state_jacobian(state));
  forcing_1 = tmc_dense_copy(tmc_state_variables(state));
  tmc_dense_fill(forcing_1, 0.0);
  model.add_forcing_terms(model.context, tmc_state_custom_rate_parameters(state), tmc_state_variables(state), forcing_1);

  /* FO2 gas -> mode 2 CORGE: a loss of FO2 in the gas phase, an equal gain in mode 2 CORGE */
  int fo2_gas_index = tmc_state_variable_index(state, "FO2");
  int fo2_mode2_index = tmc_state_variable_index(state, "STUB1.MODE2.CORGE.FO2");
  if (fo2_gas_index < 0 || fo2_mode2_index < 0)
  {
    tmc_test_fail(__FILE__, __LINE__, "[%s] FO2 or STUB1.MODE2.CORGE.FO2 is not in the variable map", p->name);
    goto done;
  }
  double expected_fo2_loss = -STUB1_RATE_CONSTANT_FO2_CORGE * value(state, 0, "FO2");
  double expected_fo2_gain = -expected_fo2_loss;
  P_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing_1, 0, fo2_gas_index), expected_fo2_loss);
  P_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing_1, 0, fo2_mode2_index), expected_fo2_gain);

  /* BAZ mode 1 -> mode 2 QUUX */
  int baz_mode1_index = tmc_state_variable_index(state, "STUB1.MODE1.QUUX.BAZ");
  int baz_mode2_index = tmc_state_variable_index(state, "STUB1.MODE2.QUUX.BAZ");
  if (baz_mode1_index < 0 || baz_mode2_index < 0)
  {
    tmc_test_fail(__FILE__, __LINE__, "[%s] STUB1.MODE1.QUUX.BAZ or STUB1.MODE2.QUUX.BAZ is not in the variable map", p->name);
    goto done;
  }
  double expected_baz_loss = -STUB1_RATE_CONSTANT_BAZ_QUUX * value(state, 0, "STUB1.MODE1.QUUX.BAZ");
  double expected_baz_gain = -expected_baz_loss;
  P_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing_1, 0, baz_mode1_index), expected_baz_loss);
  P_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing_1, 0, baz_mode2_index), expected_baz_gain);

done:
  tmc_dense_free(forcing_1);
  tmc_state_free(state);
  tmc_solver_free(solver);
  fixture_free(&f);
}

/* ---- TestSingleCellJacobianWithStubAerosolModel -------------------------------- */

void test_single_cell_jacobian_with_stub_aerosol_model(const aerosol_policy *p)
{
  fixture f;
  fixture_init(&f);
  tmc_solver *solver = build_solver(p, &f);
  if (!solver)
  {
    fixture_free(&f);
    return;
  }
  tmc_state *state = tmc_solver_create_state(solver, 1);
  set_forcing_state(p, state);

  /* jacobian_1 = state.jacobian_ in the policy's numbering, zeroed; the model
   * binds its elements to it and subtracts its terms */
  tmc_sparse *jacobian_1 = jacobian_in_policy_ordering(p, tmc_state_jacobian(state));
  P_EXPECT_TRUE(jacobian_1->ordering == p->ordering);
  tmc_sparse_fill(jacobian_1, 0.0);
  stub_aerosol_1 aerosol_1 = f.aerosol_1;
  tmc_external_model model = stub_aerosol_1_model(&aerosol_1);
  model.finalize_process_setup(model.context, tmc_solver_parameter_map(solver), tmc_solver_variable_map(solver), jacobian_1);
  model.subtract_jacobian_terms(model.context, tmc_state_custom_rate_parameters(state), tmc_state_variables(state), jacobian_1);

  /* FO2: d(loss)/d(FO2) = -k in (FO2, FO2) and d(gain)/d(FO2) = k in (mode 2 FO2, FO2);
   * the Jacobian holds the negative partial derivatives */
  int fo2_gas_index = tmc_state_variable_index(state, "FO2");
  int fo2_mode2_index = tmc_state_variable_index(state, "STUB1.MODE2.CORGE.FO2");
  if (fo2_gas_index < 0 || fo2_mode2_index < 0)
  {
    tmc_test_fail(__FILE__, __LINE__, "[%s] FO2 or STUB1.MODE2.CORGE.FO2 is not in the variable map", p->name);
    goto done;
  }
  double expected_fo2_gas_partial = -STUB1_RATE_CONSTANT_FO2_CORGE;
  double expected_fo2_mode2_partial = STUB1_RATE_CONSTANT_FO2_CORGE;
  P_EXPECT_TRUE(!tmc_sparse_is_zero(jacobian_1, fo2_gas_index, fo2_gas_index));
  P_EXPECT_TRUE(!tmc_sparse_is_zero(jacobian_1, fo2_mode2_index, fo2_gas_index));
  P_EXPECT_DOUBLE_EQ(tmc_sparse_get(jacobian_1, 0, fo2_gas_index, fo2_gas_index), -expected_fo2_gas_partial);
  P_EXPECT_DOUBLE_EQ(tmc_sparse_get(jacobian_1, 0, fo2_mode2_index, fo2_gas_index), -expected_fo2_mode2_partial);

  /* BAZ mode 1 -> mode 2 QUUX */
  int baz_mode1_index = tmc_state_variable_index(state, "STUB1.MODE1.QUUX.BAZ");
  int baz_mode2_index = tmc_state_variable_index(state, "STUB1.MODE2.QUUX.BAZ");
  if (baz_mode1_index < 0 || baz_mode2_index < 0)
  {
    tmc_test_fail(__FILE__, __LINE__, "[%s] STUB1.MODE1.QUUX.BAZ or STUB1.MODE2.QUUX.BAZ is not in the variable map", p->name);
    goto done;
  }
  double expected_baz_mode1_partial = -STUB1_RATE_CONSTANT_BAZ_QUUX;
  double expected_baz_mode2_partial = STUB1_RATE_CONSTANT_BAZ_QUUX;
  P_EXPECT_TRUE(!tmc_sparse_is_zero(jacobian_1, baz_mode1_index, baz_mode1_index));
  P_EXPECT_TRUE(!tmc_sparse_is_zero(jacobian_1, baz_mode2_index, baz_mode1_index));
  P_EXPECT_DOUBLE_EQ(tmc_sparse_get(jacobian_1, 0, baz_mode1_index, baz_mode1_index), -expected_baz_mode1_partial);
  P_EXPECT_DOUBLE_EQ(tmc_sparse_get(jacobian_1, 0, baz_mode2_index, baz_mode1_index), -expected_baz_mode2_partial);

done:
  tmc_sparse_free(jacobian_1);
  tmc_state_free(state);
  tmc_solver_free(solver);
  fixture_free(&f);
}

/* ---- TestSolveWithStubAerosolModel1 -------------------------------------------- */

void test_solve_with_stub_aerosol_model_1(const aerosol_policy *p, double base_relative_tolerance)
{
  fixture f;
  fixture_init(&f);
  tmc_solver *solver = build_solver(p, &f);
  if (!solver)
  {
    fixture_free(&f);
    return;
  }
  tmc_state *state = tmc_solver_create_state(solver, 1);

  const double fo2_initial = 1.0;
  const double baz_mode1_initial = 0.5;
  const double fo2_mode2_initial = 0.8;
  set_value(p, state, "FO2", fo2_initial);
  set_value(p, state, "BAR", 2.0);
  set_value(p, state, "STUB1.MODE1.QUUX.BAZ", baz_mode1_initial);
  set_value(p, state, "STUB1.MODE2.CORGE.FO2", fo2_mode2_initial);
  set_value(p, state, "STUB2.MODE3.CORGE.QUX", 0.3);
  set_value(p, state, "STUB2.MODE3.CORGE.BAZ", 0.2);
  set_value(p, state, "STUB2.MODE1.NUMBER", 1000.0);
  set_value(p, state, "STUB2.MODE2.NUMBER", 500.0);

  /* the analytical solution */
  const double time_step = 10.0; /* s */
  const double stub1_rxn1_delta = fo2_initial * (1.0 - exp(-STUB1_RATE_CONSTANT_FO2_CORGE * time_step));
  const double stub1_rxn2_delta = baz_mode1_initial * (1.0 - exp(-STUB1_RATE_CONSTANT_BAZ_QUUX * time_step));

  tmc_solver_update_state_parameters(solver, state);
  tmc_result results = tmc_solver_solve(solver, time_step, state);

  P_EXPECT_INT_EQ(results.state, TMC_CONVERGED);
  P_EXPECT_NEAR(value(state, 0, "FO2"), fo2_initial - stub1_rxn1_delta, base_relative_tolerance * fo2_initial);
  P_EXPECT_EQ(value(state, 0, "BAR"), 2.0);
  P_EXPECT_NEAR(value(state, 0, "STUB1.MODE1.QUUX.BAZ"), baz_mode1_initial - stub1_rxn2_delta, base_relative_tolerance * baz_mode1_initial);
  P_EXPECT_NEAR(value(state, 0, "STUB1.MODE2.QUUX.BAZ"), stub1_rxn2_delta, base_relative_tolerance * baz_mode1_initial);
  P_EXPECT_NEAR(value(state, 0, "STUB1.MODE2.CORGE.FO2"), fo2_mode2_initial + stub1_rxn1_delta, base_relative_tolerance * fo2_initial);
  P_EXPECT_EQ(value(state, 0, "STUB2.MODE3.CORGE.QUX"), 0.3);
  P_EXPECT_EQ(value(state, 0, "STUB2.MODE3.CORGE.BAZ"), 0.2);
  P_EXPECT_EQ(value(state, 0, "STUB2.MODE1.NUMBER"), 1000.0);
  P_EXPECT_EQ(value(state, 0, "STUB2.MODE2.NUMBER"), 500.0);

  tmc_state_free(state);
  tmc_solver_free(solver);
  fixture_free(&f);
}

/* ---- TestSolveWithTwoStubAerosolModels ----------------------------------------- */

#define FO2_TO_BAZ "STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT"

void test_solve_with_two_stub_aerosol_models(const aerosol_policy *p, double base_relative_tolerance)
{
  fixture f;
  fixture_init(&f);
  tmc_solver *solver = build_solver(p, &f);
  if (!solver)
  {
    fixture_free(&f);
    return;
  }
  tmc_state *state = tmc_solver_create_state(solver, 1);

  const double fo2_initial = 1.0;
  const double baz_mode1_initial = 0.5;
  const double fo2_mode2_initial = 0.8;
  const double stub2_mode2_fo2_initial = 0.4;
  const double stub2_mode2_baz_initial = 0.1;
  const double stub2_mode3_baz_initial = 0.2;
  const double stub2_mode3_qux_initial = 0.3;
  const double temperature = 275.0;
  const double fo2_to_baz_rate_constant = 0.01;
  set_value(p, state, "FO2", fo2_initial);
  set_value(p, state, "BAR", 2.0);
  set_value(p, state, "STUB1.MODE1.QUUX.BAZ", baz_mode1_initial);
  set_value(p, state, "STUB1.MODE2.CORGE.FO2", fo2_mode2_initial);
  set_value(p, state, "STUB2.MODE2.CORGE.FO2", stub2_mode2_fo2_initial);
  set_value(p, state, "STUB2.MODE2.CORGE.BAZ", stub2_mode2_baz_initial);
  set_value(p, state, "STUB2.MODE3.QUUX.BAZ", stub2_mode3_baz_initial);
  set_value(p, state, "STUB2.MODE3.QUUX.QUX", stub2_mode3_qux_initial);
  set_value(p, state, "STUB2.MODE3.CORGE.BAZ", 0.2);
  set_value(p, state, "STUB2.MODE1.NUMBER", 1000.0);
  set_value(p, state, "STUB2.MODE2.NUMBER", 500.0);

  int param = tmc_state_custom_parameter_index(state, FO2_TO_BAZ);
  if (param < 0)
  {
    tmc_test_fail(__FILE__, __LINE__, "[%s] %s is not a state parameter", p->name, FO2_TO_BAZ);
    goto done;
  }
  tmc_dense_set(tmc_state_custom_rate_parameters(state), 0, param, fo2_to_baz_rate_constant);
  tmc_state_conditions(state)[0].temperature = temperature;

  /* the analytical solution */
  const double time_step = 10.0; /* s */
  const double stub1_rxn1_delta = fo2_initial * (1.0 - exp(-STUB1_RATE_CONSTANT_FO2_CORGE * time_step));
  const double stub1_rxn2_delta = baz_mode1_initial * (1.0 - exp(-STUB1_RATE_CONSTANT_BAZ_QUUX * time_step));
  const double stub2_rxn1_delta = stub2_mode2_fo2_initial * (1.0 - exp(-fo2_to_baz_rate_constant * time_step));
  const double stub2_rxn2_delta = stub2_mode3_baz_initial * (1.0 - exp(-temperature * 0.005 * time_step));

  tmc_solver_update_state_parameters(solver, state);
  tmc_result results = tmc_solver_solve(solver, time_step, state);

  P_EXPECT_INT_EQ(results.state, TMC_CONVERGED);
  P_EXPECT_NEAR(value(state, 0, "FO2"), fo2_initial - stub1_rxn1_delta, base_relative_tolerance * fo2_initial);
  P_EXPECT_EQ(value(state, 0, "BAR"), 2.0);
  P_EXPECT_NEAR(value(state, 0, "STUB1.MODE1.QUUX.BAZ"), baz_mode1_initial - stub1_rxn2_delta, base_relative_tolerance * baz_mode1_initial);
  P_EXPECT_NEAR(value(state, 0, "STUB1.MODE2.QUUX.BAZ"), stub1_rxn2_delta, base_relative_tolerance * baz_mode1_initial);
  P_EXPECT_NEAR(value(state, 0, "STUB1.MODE2.CORGE.FO2"), fo2_mode2_initial + stub1_rxn1_delta, base_relative_tolerance * fo2_initial);
  P_EXPECT_NEAR(value(state, 0, "STUB2.MODE2.CORGE.FO2"), stub2_mode2_fo2_initial - stub2_rxn1_delta, base_relative_tolerance * stub2_mode2_fo2_initial);
  P_EXPECT_NEAR(value(state, 0, "STUB2.MODE2.CORGE.BAZ"), stub2_mode2_baz_initial + stub2_rxn1_delta, base_relative_tolerance * stub2_mode2_fo2_initial);
  P_EXPECT_NEAR(value(state, 0, "STUB2.MODE3.QUUX.BAZ"), stub2_mode3_baz_initial - stub2_rxn2_delta, base_relative_tolerance * stub2_mode3_baz_initial);
  P_EXPECT_NEAR(value(state, 0, "STUB2.MODE3.QUUX.QUX"), stub2_mode3_qux_initial + stub2_rxn2_delta, base_relative_tolerance * stub2_mode3_baz_initial);
  P_EXPECT_EQ(value(state, 0, "STUB2.MODE1.NUMBER"), 1000.0);
  P_EXPECT_EQ(value(state, 0, "STUB2.MODE2.NUMBER"), 500.0);

done:
  tmc_state_free(state);
  tmc_solver_free(solver);
  fixture_free(&f);
}

/* ---- TestSolveWithStubAerosolModel1MultiCell ----------------------------------- */

void test_solve_with_stub_aerosol_model_1_multi_cell(const aerosol_policy *p, double base_relative_tolerance)
{
  fixture f;
  fixture_init(&f);
  tmc_solver *solver = build_solver(p, &f);
  if (!solver)
  {
    fixture_free(&f);
    return;
  }
  enum
  {
    num_cells = 3
  };
  tmc_state *state = tmc_solver_create_state(solver, num_cells);

  const double fo2_initial[num_cells] = { 1.0, 1.5, 2.0 };
  const double baz_mode1_initial[num_cells] = { 0.5, 0.7, 0.9 };
  const double fo2_mode2_initial[num_cells] = { 0.8, 1.0, 1.2 };
  const double bar[num_cells] = { 2.0, 2.5, 3.0 };
  const double stub2_qux[num_cells] = { 0.3, 0.4, 0.5 };
  const double stub2_baz[num_cells] = { 0.2, 0.3, 0.4 };
  const double stub2_num1[num_cells] = { 1000.0, 1100.0, 1200.0 };
  const double stub2_num2[num_cells] = { 500.0, 550.0, 600.0 };
  set_values(p, state, "FO2", fo2_initial, num_cells);
  set_values(p, state, "BAR", bar, num_cells);
  set_values(p, state, "STUB1.MODE1.QUUX.BAZ", baz_mode1_initial, num_cells);
  set_values(p, state, "STUB1.MODE2.CORGE.FO2", fo2_mode2_initial, num_cells);
  set_values(p, state, "STUB2.MODE3.CORGE.QUX", stub2_qux, num_cells);
  set_values(p, state, "STUB2.MODE3.CORGE.BAZ", stub2_baz, num_cells);
  set_values(p, state, "STUB2.MODE1.NUMBER", stub2_num1, num_cells);
  set_values(p, state, "STUB2.MODE2.NUMBER", stub2_num2, num_cells);

  /* the analytical solution in each cell */
  const double time_step = 10.0; /* s */
  double stub1_rxn1_delta[num_cells], stub1_rxn2_delta[num_cells];
  for (int i = 0; i < num_cells; ++i)
  {
    stub1_rxn1_delta[i] = fo2_initial[i] * (1.0 - exp(-STUB1_RATE_CONSTANT_FO2_CORGE * time_step));
    stub1_rxn2_delta[i] = baz_mode1_initial[i] * (1.0 - exp(-STUB1_RATE_CONSTANT_BAZ_QUUX * time_step));
  }

  tmc_solver_update_state_parameters(solver, state);
  tmc_result results = tmc_solver_solve(solver, time_step, state);

  P_EXPECT_INT_EQ(results.state, TMC_CONVERGED);
  for (int i = 0; i < num_cells; ++i)
  {
    P_EXPECT_NEAR(value(state, i, "FO2"), fo2_initial[i] - stub1_rxn1_delta[i], base_relative_tolerance * fo2_initial[i]);
    P_EXPECT_EQ(value(state, i, "BAR"), bar[i]);
    P_EXPECT_NEAR(value(state, i, "STUB1.MODE1.QUUX.BAZ"), baz_mode1_initial[i] - stub1_rxn2_delta[i], base_relative_tolerance * baz_mode1_initial[i]);
    P_EXPECT_NEAR(value(state, i, "STUB1.MODE2.QUUX.BAZ"), stub1_rxn2_delta[i], base_relative_tolerance * baz_mode1_initial[i]);
    P_EXPECT_NEAR(value(state, i, "STUB1.MODE2.CORGE.FO2"), fo2_mode2_initial[i] + stub1_rxn1_delta[i], base_relative_tolerance * fo2_initial[i]);
    P_EXPECT_EQ(value(state, i, "STUB2.MODE3.CORGE.QUX"), stub2_qux[i]);
    P_EXPECT_EQ(value(state, i, "STUB2.MODE3.CORGE.BAZ"), stub2_baz[i]);
    P_EXPECT_EQ(value(state, i, "STUB2.MODE1.NUMBER"), stub2_num1[i]);
    P_EXPECT_EQ(value(state, i, "STUB2.MODE2.NUMBER"), stub2_num2[i]);
  }

  tmc_state_free(state);
  tmc_solver_free(solver);
  fixture_free(&f);
}

/* ---- TestSolveWithTwoStubAerosolModelsMultiCell -------------------------------- */

void test_solve_with_two_stub_aerosol_models_multi_cell(const aerosol_policy *p, double base_relative_tolerance)
{
  fixture f;
  fixture_init(&f);
  tmc_solver *solver = build_solver(p, &f);
  if (!solver)
  {
    fixture_free(&f);
    return;
  }
  enum
  {
    num_cells = 3
  };
  tmc_state *state = tmc_solver_create_state(solver, num_cells);

  const double fo2_initial[num_cells] = { 1.0, 1.5, 2.0 };
  const double baz_mode1_initial[num_cells] = { 0.5, 0.7, 0.9 };
  const double fo2_mode2_initial[num_cells] = { 0.8, 1.0, 1.2 };
  const double stub2_mode2_fo2_initial[num_cells] = { 0.4, 0.6, 0.8 };
  const double stub2_mode2_baz_initial[num_cells] = { 0.1, 0.15, 0.2 };
  const double stub2_mode3_baz_initial[num_cells] = { 0.2, 0.3, 0.4 };
  const double stub2_mode3_qux_initial[num_cells] = { 0.3, 0.4, 0.5 };
  const double temperature[num_cells] = { 275.0, 285.0, 295.0 };
  const double fo2_to_baz_rate_constant[num_cells] = { 0.01, 0.015, 0.02 };
  const double bar[num_cells] = { 2.0, 2.5, 3.0 };
  const double stub2_corge_baz[num_cells] = { 0.2, 0.25, 0.3 };
  const double stub2_num1[num_cells] = { 1000.0, 1100.0, 1200.0 };
  const double stub2_num2[num_cells] = { 500.0, 550.0, 600.0 };

  set_values(p, state, "FO2", fo2_initial, num_cells);
  set_values(p, state, "BAR", bar, num_cells);
  set_values(p, state, "STUB1.MODE1.QUUX.BAZ", baz_mode1_initial, num_cells);
  set_values(p, state, "STUB1.MODE2.CORGE.FO2", fo2_mode2_initial, num_cells);
  set_values(p, state, "STUB2.MODE2.CORGE.FO2", stub2_mode2_fo2_initial, num_cells);
  set_values(p, state, "STUB2.MODE2.CORGE.BAZ", stub2_mode2_baz_initial, num_cells);
  set_values(p, state, "STUB2.MODE3.QUUX.BAZ", stub2_mode3_baz_initial, num_cells);
  set_values(p, state, "STUB2.MODE3.QUUX.QUX", stub2_mode3_qux_initial, num_cells);
  set_values(p, state, "STUB2.MODE3.CORGE.BAZ", stub2_corge_baz, num_cells);
  set_values(p, state, "STUB2.MODE1.NUMBER", stub2_num1, num_cells);
  set_values(p, state, "STUB2.MODE2.NUMBER", stub2_num2, num_cells);

  int param = tmc_state_custom_parameter_index(state, FO2_TO_BAZ);
  if (param < 0)
  {
    tmc_test_fail(__FILE__, __LINE__, "[%s] %s is not a state parameter", p->name, FO2_TO_BAZ);
    goto done;
  }
  for (int i = 0; i < num_cells; ++i)
  {
    tmc_dense_set(tmc_state_custom_rate_parameters(state), i, param, fo2_to_baz_rate_constant[i]);
    tmc_state_conditions(state)[i].temperature = temperature[i];
  }

  /* the analytical solution in each cell */
  const double time_step = 10.0; /* s */
  double stub1_rxn1_delta[num_cells], stub1_rxn2_delta[num_cells], stub2_rxn1_delta[num_cells], stub2_rxn2_delta[num_cells];
  for (int i = 0; i < num_cells; ++i)
  {
    stub1_rxn1_delta[i] = fo2_initial[i] * (1.0 - exp(-STUB1_RATE_CONSTANT_FO2_CORGE * time_step));
    stub1_rxn2_delta[i] = baz_mode1_initial[i] * (1.0 - exp(-STUB1_RATE_CONSTANT_BAZ_QUUX * time_step));
    stub2_rxn1_delta[i] = stub2_mode2_fo2_initial[i] * (1.0 - exp(-fo2_to_baz_rate_constant[i] * time_step));
    stub2_rxn2_delta[i] = stub2_mode3_baz_initial[i] * (1.0 - exp(-temperature[i] * 0.005 * time_step));
  }

  tmc_solver_update_state_parameters(solver, state);
  tmc_result results = tmc_solver_solve(solver, time_step, state);

  P_EXPECT_INT_EQ(results.state, TMC_CONVERGED);
  for (int i = 0; i < num_cells; ++i)
  {
    P_EXPECT_NEAR(value(state, i, "FO2"), fo2_initial[i] - stub1_rxn1_delta[i], base_relative_tolerance * fo2_initial[i]);
    P_EXPECT_EQ(value(state, i, "BAR"), bar[i]);
    P_EXPECT_NEAR(value(state, i, "STUB1.MODE1.QUUX.BAZ"), baz_mode1_initial[i] - stub1_rxn2_delta[i], base_relative_tolerance * baz_mode1_initial[i]);
    P_EXPECT_NEAR(value(state, i, "STUB1.MODE2.QUUX.BAZ"), stub1_rxn2_delta[i], base_relative_tolerance * baz_mode1_initial[i]);
    P_EXPECT_NEAR(value(state, i, "STUB1.MODE2.CORGE.FO2"), fo2_mode2_initial[i] + stub1_rxn1_delta[i], base_relative_tolerance * fo2_initial[i]);
    P_EXPECT_NEAR(value(state, i, "STUB2.MODE2.CORGE.FO2"), stub2_mode2_fo2_initial[i] - stub2_rxn1_delta[i], base_relative_tolerance * stub2_mode2_fo2_initial[i]);
    P_EXPECT_NEAR(value(state, i, "STUB2.MODE2.CORGE.BAZ"), stub2_mode2_baz_initial[i] + stub2_rxn1_delta[i], base_relative_tolerance * stub2_mode2_fo2_initial[i]);
    P_EXPECT_NEAR(value(state, i, "STUB2.MODE3.QUUX.BAZ"), stub2_mode3_baz_initial[i] - stub2_rxn2_delta[i], base_relative_tolerance * stub2_mode3_baz_initial[i]);
    P_EXPECT_NEAR(value(state, i, "STUB2.MODE3.QUUX.QUX"), stub2_mode3_qux_initial[i] + stub2_rxn2_delta[i], base_relative_tolerance * stub2_mode3_baz_initial[i]);
    P_EXPECT_EQ(value(state, i, "STUB2.MODE1.NUMBER"), stub2_num1[i]);
    P_EXPECT_EQ(value(state, i, "STUB2.MODE2.NUMBER"), stub2_num2[i]);
  }

done:
  tmc_state_free(state);
  tmc_solver_free(solver);
  fixture_free(&f);
}
