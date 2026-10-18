/* MICM test/unit/solver/test_state.cpp.
 *
 * MICM builds a State directly from StateParameters (variable names, number
 * of rate constants, custom parameter labels); in C a state comes from a
 * solver, so each test builds a mechanism with exactly those: the variables
 * foo, bar, baz, quz in order, ten reactions, and one user-defined rate
 * constant per custom label. MICM's state is standard ordering (L = 1); every
 * check also runs at L = 4, where the cells share padded blocks.
 *
 * C equivalents: copying is tmc_state_copy, a move hands over the pointer,
 * and operator[] is element access through tmc_state_variable_index. Where
 * MICM throws, the C call returns an error category and code. */
#include "harness.h"
#include "tmc_solver_ext.h"

#include <float.h>

static const int widths[2] = { 1, 4 };

static const char *const names[4] = { "foo", "bar", "baz", "quz" };
static const char *const two_labels[2] = { "quux", "corge" };
static const char *const five_labels[5] = { "O1", "O2", "O3", "AAA", "BBB" };

/* A solver whose states have the variables foo, bar, baz, quz, ten rate
 * constants, and the given custom parameter labels, in order. */
static tmc_solver *state_solver(const char *const *labels, int n_labels, int L)
{
  tmc_species species[4];
  tmc_phase gas = tmc_phase_create("gas");
  for (int i = 0; i < 4; ++i)
  {
    species[i] = tmc_species_create(names[i]);
    tmc_phase_add_species(&gas, &species[i]);
  }
  tmc_system system = tmc_system_create(&gas);
  tmc_reaction reactions[10];
  for (int r = 0; r < 10; ++r)
  {
    reactions[r] = tmc_reaction_create();
    tmc_reaction_add_reactant(&reactions[r], &species[r % 4]);
    tmc_reaction_add_product(&reactions[r], &species[(r + 1) % 4], 1.0);
    tmc_reaction_set_phase(&reactions[r], &gas);
    if (r < n_labels)
      tmc_reaction_set_rate_constant(&reactions[r], tmc_rate_user_defined(labels[r], 1.0));
    else
    {
      tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
      p.A = 1.0e-3 * (r + 1);
      tmc_reaction_set_rate_constant(&reactions[r], tmc_rate_arrhenius(p));
    }
  }
  tmc_builder *builder = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(builder, &system);
  tmc_builder_set_reactions(builder, reactions, 10);
  tmc_builder_set_reorder_state(builder, false);
  tmc_builder_set_width(builder, L);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_solver *solver = tmc_builder_build(builder, &err);
  if (!solver)
    fprintf(stderr, "  build failed: %s\n", err.message);
  tmc_builder_free(builder);
  for (int r = 0; r < 10; ++r)
    tmc_reaction_free(&reactions[r]);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return solver;
}

static int variable(const tmc_state *state, const char *name)
{
  return tmc_state_variable_index(state, name);
}

static double custom_parameter(tmc_state *state, int cell, int index)
{
  return tmc_dense_get(tmc_state_custom_rate_parameters(state), cell, index);
}

static bool is_state_error(int category, const tmc_error *err, int code)
{
  return category == TMC_ERROR_STATE && err->category == TMC_ERROR_STATE && err->code == code;
}

/* Distinct values in every element, so copies can be compared meaningfully
 * (MICM's states hold zeros). */
static void fill_distinct(tmc_state *state)
{
  for (int c = 0; c < tmc_state_num_cells(state); ++c)
  {
    for (int v = 0; v < tmc_state_num_variables(state); ++v)
      tmc_state_set(state, c, v, 10.0 * c + v + 0.5);
    tmc_dense *custom = tmc_state_custom_rate_parameters(state);
    for (int p = 0; p < custom->n_cols; ++p)
      tmc_dense_set(custom, c, p, 100.0 + 10.0 * c + p);
    tmc_dense *rates = tmc_state_rate_constants(state);
    for (int r = 0; r < rates->n_cols; ++r)
      tmc_dense_set(rates, c, r, 1000.0 + 10.0 * c + r);
    tmc_state_conditions(state)[c].temperature = 250.0 + c;
    tmc_state_conditions(state)[c].pressure = 90000.0 + c;
    tmc_state_conditions(state)[c].air_density = 30.0 + c;
  }
}

static void set_tolerances(tmc_state *state)
{
  tmc_error err;
  tmc_error_clear(&err);
  const double atol[4] = { 1e-10, 1e-10, 1e-10, 1e-10 };
  tmc_state_set_relative_tolerance(state, 1e-05);
  TMC_EXPECT_EQ(tmc_state_set_absolute_tolerances(state, atol, 4, &err), TMC_OK);
}

/* Every property MICM's copy and move tests compare, and the conditions. */
static void expect_same_state(tmc_state *a, tmc_state *b)
{
  TMC_EXPECT_EQ(tmc_state_num_cells(a), tmc_state_num_cells(b));
  TMC_EXPECT_EQ(tmc_state_num_variables(a), tmc_state_num_variables(b));
  for (int v = 0; v < tmc_state_num_variables(a) && v < tmc_state_num_variables(b); ++v)
    TMC_EXPECT_STREQ(tmc_state_variable_name(a, v), tmc_state_variable_name(b, v));
  for (int p = 0; p < 2; ++p)
    TMC_EXPECT_EQ(tmc_state_custom_parameter_index(a, two_labels[p]), tmc_state_custom_parameter_index(b, two_labels[p]));
  const tmc_dense *da[3] = { tmc_state_variables(a), tmc_state_custom_rate_parameters(a), tmc_state_rate_constants(a) };
  const tmc_dense *db[3] = { tmc_state_variables(b), tmc_state_custom_rate_parameters(b), tmc_state_rate_constants(b) };
  for (int m = 0; m < 3; ++m)
  {
    TMC_EXPECT_EQ(da[m]->n_rows, db[m]->n_rows);
    TMC_EXPECT_EQ(da[m]->n_cols, db[m]->n_cols);
    if (da[m]->n_rows != db[m]->n_rows || da[m]->n_cols != db[m]->n_cols)
      continue;
    for (int i = 0; i < da[m]->n_rows; ++i)
      for (int j = 0; j < da[m]->n_cols; ++j)
        TMC_EXPECT_EQ(tmc_dense_get(da[m], i, j), tmc_dense_get(db[m], i, j));
  }
  TMC_EXPECT_EQ(tmc_state_relative_tolerance(a), tmc_state_relative_tolerance(b));
  for (int v = 0; v < tmc_state_num_variables(a) && v < tmc_state_num_variables(b); ++v)
    TMC_EXPECT_EQ(tmc_state_absolute_tolerances(a)[v], tmc_state_absolute_tolerances(b)[v]);
  for (int c = 0; c < tmc_state_num_cells(a) && c < tmc_state_num_cells(b); ++c)
  {
    TMC_EXPECT_EQ(tmc_state_conditions(a)[c].temperature, tmc_state_conditions(b)[c].temperature);
    TMC_EXPECT_EQ(tmc_state_conditions(a)[c].pressure, tmc_state_conditions(b)[c].pressure);
    TMC_EXPECT_EQ(tmc_state_conditions(a)[c].air_density, tmc_state_conditions(b)[c].air_density);
  }
}

/* The moved-to state of MICM's move tests: 3 cells, 4 variables, 2 custom
 * parameters, 10 rate constants, rtol 1e-5, atol 1e-10, and the values held
 * before the move. */
static void expect_moved_state(tmc_state *moved, tmc_state *expected)
{
  TMC_EXPECT_EQ(tmc_state_num_cells(moved), 3);
  TMC_EXPECT_EQ(tmc_state_variables(moved)->n_rows, 3);
  TMC_EXPECT_EQ(tmc_state_variables(moved)->n_cols, 4);
  TMC_EXPECT_EQ(tmc_state_custom_rate_parameters(moved)->n_rows, 3);
  TMC_EXPECT_EQ(tmc_state_custom_rate_parameters(moved)->n_cols, 2);
  TMC_EXPECT_EQ(tmc_state_rate_constants(moved)->n_rows, 3);
  TMC_EXPECT_EQ(tmc_state_rate_constants(moved)->n_cols, 10);
  TMC_EXPECT_NEAR(tmc_state_relative_tolerance(moved), 1e-05, 1e-05 * 1e1 * DBL_EPSILON);
  for (int v = 0; v < 4; ++v)
    TMC_EXPECT_EQ(tmc_state_absolute_tolerances(moved)[v], 1e-10);
  expect_same_state(moved, expected);
}

/* ---- construction ---------------------------------------------------------- */

TMC_TEST(default_constructor)
{
  /* MICM default-constructs an empty State and a vector of ten; in C an
   * unset state pointer is NULL, freeing it is a no-op, and a state can hold
   * zero cells. */
  tmc_state_free(NULL);
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *states[10];
    for (int i = 0; i < 10; ++i)
    {
      states[i] = tmc_solver_create_state(solver, 0);
      TMC_EXPECT_TRUE(states[i] != NULL);
      TMC_EXPECT_EQ(tmc_state_num_cells(states[i]), 0);
    }
    for (int i = 0; i < 10; ++i)
      tmc_state_free(states[i]);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(constructor)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 3);
    TMC_EXPECT_EQ(tmc_state_num_cells(state), 3);
    TMC_EXPECT_EQ(variable(state, "foo"), 0);
    TMC_EXPECT_EQ(variable(state, "bar"), 1);
    TMC_EXPECT_EQ(variable(state, "baz"), 2);
    TMC_EXPECT_EQ(variable(state, "quz"), 3);
    TMC_EXPECT_EQ(tmc_state_variables(state)->n_rows, 3);
    TMC_EXPECT_EQ(tmc_state_variables(state)->n_cols, 4);
    TMC_EXPECT_EQ(tmc_state_custom_parameter_index(state, "quux"), 0);
    TMC_EXPECT_EQ(tmc_state_custom_parameter_index(state, "corge"), 1);
    TMC_EXPECT_EQ(tmc_state_custom_rate_parameters(state)->n_rows, 3);
    TMC_EXPECT_EQ(tmc_state_custom_rate_parameters(state)->n_cols, 2);
    TMC_EXPECT_EQ(tmc_state_rate_constants(state)->n_rows, 3);
    TMC_EXPECT_EQ(tmc_state_rate_constants(state)->n_cols, 10);
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

/* ---- copies and moves ------------------------------------------------------ */

TMC_TEST(copy_constructor)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *original = tmc_solver_create_state(solver, 3);
    set_tolerances(original);
    fill_distinct(original);
    tmc_state *copy = tmc_state_copy(original);
    expect_same_state(copy, original);
    /* the copy owns its data */
    tmc_state_set(copy, 1, 2, -1.0);
    TMC_EXPECT_EQ(tmc_state_get(original, 1, 2), 12.5);
    tmc_state_free(copy);
    tmc_state_free(original);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(copy_assignment_operator)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *original = tmc_solver_create_state(solver, 3);
    set_tolerances(original);
    fill_distinct(original);
    /* assignment replaces an existing state */
    tmc_state *copy = tmc_solver_create_state(solver, 1);
    tmc_state_free(copy);
    copy = tmc_state_copy(original);
    expect_same_state(copy, original);
    tmc_state_set(copy, 2, 3, -1.0);
    TMC_EXPECT_EQ(tmc_state_get(original, 2, 3), 23.5);
    tmc_state_free(copy);
    tmc_state_free(original);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(move_constructor)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *original = tmc_solver_create_state(solver, 3);
    set_tolerances(original);
    fill_distinct(original);
    tmc_state *expected = tmc_state_copy(original);
    tmc_state *moved = original;
    original = NULL;
    expect_moved_state(moved, expected);
    tmc_state_free(moved);
    tmc_state_free(expected);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(move_assignment_operator)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *original = tmc_solver_create_state(solver, 3);
    set_tolerances(original);
    fill_distinct(original);
    tmc_state *expected = tmc_state_copy(original);
    tmc_state *moved = tmc_solver_create_state(solver, 1);
    tmc_state_free(moved);
    moved = original;
    original = NULL;
    expect_moved_state(moved, expected);
    tmc_state_free(moved);
    tmc_state_free(expected);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

/* ---- concentrations ---------------------------------------------------------- */

TMC_TEST(setting_single_concentration_with_invalid_arguments_throws_exception)
{
  int failures = 0;
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 3);
    tmc_error err;
    const double two[2] = { 1.0, 2.0 }, four[4] = { 1.0, 2.0, 3.0, 4.0 };

    /* MICM rejects one value for a state of three cells
     * (MICM_STATE_ERROR_CODE_INVALID_CONCENTRATION_COUNT_MULTIGRID) */
    tmc_error_clear(&err);
    int category = tmc_state_set_concentration(state, "foo", 1.0, &err);
    if (!is_state_error(category, &err, TMC_STATE_INVALID_CONCENTRATION_COUNT_MULTIGRID))
      ++failures;

    tmc_error_clear(&err);
    TMC_EXPECT_TRUE(is_state_error(tmc_state_set_concentrations(state, "foo", two, 2, &err), &err, TMC_STATE_INVALID_CONCENTRATION_COUNT_MULTIGRID));
    tmc_error_clear(&err);
    TMC_EXPECT_TRUE(is_state_error(tmc_state_set_concentration(state, "not foo", 1.0, &err), &err, TMC_STATE_UNKNOWN_SPECIES));
    tmc_error_clear(&err);
    TMC_EXPECT_TRUE(is_state_error(tmc_state_set_concentrations(state, "not foo", four, 4, &err), &err, TMC_STATE_UNKNOWN_SPECIES));
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  if (failures)
  {
    tmc_test_set_reason("tmc_state_set_concentration sets one value in every cell of a multi-cell state; MICM's SetConcentration(species, value) throws INVALID_CONCENTRATION_COUNT_MULTIGRID");
    return TMC_TEST_FAIL;
  }
  return TMC_TEST_PASS;
}

TMC_TEST(set_single_concentration)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_error err;
    {
      tmc_state *state = tmc_solver_create_state(solver, 3);
      const double concentrations[3] = { 12.0, 42.0, 35.2 };
      tmc_error_clear(&err);
      TMC_EXPECT_EQ(tmc_state_set_concentrations(state, "bar", concentrations, 3, &err), TMC_OK);
      for (int i = 0; i < 3; ++i)
        TMC_EXPECT_EQ(tmc_state_get(state, i, variable(state, "bar")), concentrations[i]);
      tmc_state_free(state);
    }
    {
      tmc_state *state = tmc_solver_create_state(solver, 1);
      tmc_error_clear(&err);
      TMC_EXPECT_EQ(tmc_state_set_concentration(state, "bar", 324.2, &err), TMC_OK);
      TMC_EXPECT_NEAR(tmc_state_get(state, 0, variable(state, "bar")), 324.2, fabs(324.2) * 1e1 * DBL_EPSILON);
      tmc_state_free(state);
    }
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(set_concentration_by_element_single_value)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 1);
    tmc_error err;
    tmc_error_clear(&err);
    TMC_EXPECT_EQ(tmc_state_set_concentration(state, "baz", 42.0, &err), TMC_OK);
    TMC_EXPECT_EQ(tmc_state_set_concentration(state, "quz", 12.0, &err), TMC_OK);
    TMC_EXPECT_EQ(tmc_state_get(state, 0, variable(state, "baz")), 42.0);
    TMC_EXPECT_EQ(tmc_state_get(state, 0, variable(state, "quz")), 12.0);
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(set_concentration_by_element_vector)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 3);
    const double concentrations[3] = { 12.0, 42.0, 35.2 };
    tmc_error err;
    tmc_error_clear(&err);
    TMC_EXPECT_EQ(tmc_state_set_concentrations(state, "foo", concentrations, 3, &err), TMC_OK);
    for (int i = 0; i < 3; ++i)
      TMC_EXPECT_EQ(tmc_state_get(state, i, variable(state, "foo")), concentrations[i]);
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

/* MICM's SetConcentrations and SetCustomRateParameters (maps): every entry in
 * turn, stopping at the first error. */
typedef struct named_values
{
  const char *name;
  int n;
  double values[3];
} named_values;

static int set_concentration_map(tmc_state *state, const named_values *entries, int n_entries, tmc_error *err)
{
  for (int e = 0; e < n_entries; ++e)
  {
    int category = tmc_state_set_concentrations(state, entries[e].name, entries[e].values, entries[e].n, err);
    if (category != TMC_OK)
      return category;
  }
  return TMC_OK;
}

static int set_custom_rate_parameter_map(tmc_state *state, const named_values *entries, int n_entries, tmc_error *err)
{
  for (int e = 0; e < n_entries; ++e)
  {
    int category = tmc_state_set_custom_rate_parameters(state, entries[e].name, entries[e].values, entries[e].n, err);
    if (category != TMC_OK)
      return category;
  }
  return TMC_OK;
}

TMC_TEST(setting_concentrations_with_invalid_arguements_throws_exception)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 3);
    const named_values concentrations[4] = { { "FUU", 1, { 0.1 } }, { "bar", 1, { 0.2 } }, { "baz", 1, { 0.3 } }, { "quz", 1, { 0.4 } } };
    tmc_error err;
    tmc_error_clear(&err);
    TMC_EXPECT_TRUE(set_concentration_map(state, concentrations, 4, &err) != TMC_OK);
    /* each entry is invalid: an unknown species, or one value for three cells */
    tmc_error_clear(&err);
    TMC_EXPECT_TRUE(is_state_error(set_concentration_map(state, &concentrations[0], 1, &err), &err, TMC_STATE_UNKNOWN_SPECIES));
    for (int e = 1; e < 4; ++e)
    {
      tmc_error_clear(&err);
      TMC_EXPECT_TRUE(is_state_error(set_concentration_map(state, &concentrations[e], 1, &err), &err, TMC_STATE_INVALID_CONCENTRATION_COUNT_MULTIGRID));
    }
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(set_concentrations)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    const int num_grid_cells = 3, num_species = 4;
    tmc_state *state = tmc_solver_create_state(solver, num_grid_cells);
    const named_values concentrations[4] = { { "bar", 3, { 0.2, 0.22, 0.222 } },
                                             { "baz", 3, { 0.3, 0.33, 0.333 } },
                                             { "foo", 3, { 0.9, 0.99, 0.999 } },
                                             { "quz", 3, { 0.4, 0.44, 0.444 } } };
    tmc_error err;
    tmc_error_clear(&err);
    TMC_EXPECT_EQ(set_concentration_map(state, concentrations, 4, &err), TMC_OK);
    const double concentrations_in_order[12] = { 0.9, 0.2, 0.3, 0.4, 0.99, 0.22, 0.33, 0.44, 0.999, 0.222, 0.333, 0.444 };
    int idx = 0;
    for (int i = 0; i < num_grid_cells; ++i)
      for (int j = 0; j < num_species; ++j, ++idx)
        TMC_EXPECT_EQ(tmc_state_get(state, i, j), concentrations_in_order[idx]);
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

/* ---- element access (MICM's operator[]) ----------------------------------- */

TMC_TEST(set_state_with_square_bracket_operator)
{
  const double tol_scale = 1e1 * DBL_EPSILON;
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 1);
    const int foo = variable(state, "foo"), bar = variable(state, "bar"), baz = variable(state, "baz"), quz = variable(state, "quz");
    TMC_ASSERT_TRUE(foo >= 0 && bar >= 0 && baz >= 0 && quz >= 0);

    tmc_state_set(state, 0, foo, 42.0);
    tmc_state_set(state, 0, bar, 12.0);
    TMC_EXPECT_EQ(tmc_state_variables(state)->data[tmc_dense_index(tmc_state_variables(state), 0, foo)], 42.0);
    TMC_EXPECT_EQ(tmc_state_variables(state)->data[tmc_dense_index(tmc_state_variables(state), 0, bar)], 12.0);

    /* by species: the species' name */
    tmc_species s_baz = tmc_species_create("baz"), s_quz = tmc_species_create("quz");
    tmc_state_set(state, 0, variable(state, s_baz.name), 35.2);
    tmc_state_set(state, 0, variable(state, s_quz.name), 24.2);
    TMC_EXPECT_NEAR(tmc_state_get(state, 0, baz), 35.2, fabs(35.2) * tol_scale);
    TMC_EXPECT_NEAR(tmc_state_get(state, 0, quz), 24.2, fabs(24.2) * tol_scale);

    /* read access through a const state */
    const tmc_state *const_state = state;
    TMC_EXPECT_EQ(tmc_state_get(const_state, 0, tmc_state_variable_index(const_state, "foo")), 42.0);
    TMC_EXPECT_EQ(tmc_state_get(const_state, 0, tmc_state_variable_index(const_state, "bar")), 12.0);
    TMC_EXPECT_NEAR(tmc_state_get(const_state, 0, tmc_state_variable_index(const_state, s_baz.name)), 35.2, fabs(35.2) * tol_scale);
    TMC_EXPECT_NEAR(tmc_state_get(const_state, 0, tmc_state_variable_index(const_state, s_quz.name)), 24.2, fabs(24.2) * tol_scale);

    /* arithmetic in place */
    tmc_state_set(state, 0, foo, tmc_state_get(state, 0, foo) + 8.0);
    TMC_EXPECT_EQ(tmc_state_get(state, 0, foo), 50.0);
    tmc_state_set(state, 0, foo, tmc_state_get(state, 0, foo) - 20.0);
    TMC_EXPECT_EQ(tmc_state_get(state, 0, foo), 30.0);
    tmc_state_set(state, 0, foo, tmc_state_get(state, 0, foo) * 2.0);
    TMC_EXPECT_EQ(tmc_state_get(state, 0, foo), 60.0);
    tmc_state_set(state, 0, foo, tmc_state_get(state, 0, foo) / 3.0);
    TMC_EXPECT_EQ(tmc_state_get(state, 0, foo), 20.0);

    /* an unknown name has no index (MICM throws MICM_STATE_ERROR_CODE_UNKNOWN_SPECIES) */
    TMC_EXPECT_EQ(variable(state, "not foo"), -1);
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(set_multi_cell_state_with_square_brackets)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(two_labels, 2, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 3);
    const double foo_values[3] = { 1.0, 2.0, 3.0 }, bar_values[3] = { 4.0, 5.0, 6.0 };
    tmc_error err;
    tmc_error_clear(&err);
    TMC_EXPECT_EQ(tmc_state_set_concentrations(state, "foo", foo_values, 3, &err), TMC_OK);
    TMC_EXPECT_EQ(tmc_state_set_concentrations(state, "bar", bar_values, 3, &err), TMC_OK);
    const int foo = variable(state, "foo"), bar = variable(state, "bar");
    const tmc_dense *v = tmc_state_variables(state);
    for (int i = 0; i < 3; ++i)
    {
      TMC_EXPECT_EQ(v->data[tmc_dense_index(v, i, foo)], foo_values[i]);
      TMC_EXPECT_EQ(v->data[tmc_dense_index(v, i, bar)], bar_values[i]);
      TMC_EXPECT_EQ(tmc_state_get(state, i, foo), foo_values[i]);
      TMC_EXPECT_EQ(tmc_state_get(state, i, bar), bar_values[i]);
    }
    const tmc_state *const_state = state;
    for (int i = 0; i < 3; ++i)
    {
      TMC_EXPECT_EQ(tmc_state_get(const_state, i, foo), foo_values[i]);
      TMC_EXPECT_EQ(tmc_state_get(const_state, i, bar), bar_values[i]);
    }
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

/* ---- custom rate parameters ------------------------------------------------- */

TMC_TEST(setting_custom_rate_parameter_with_invalid_vector_size_throws_exception)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(five_labels, 5, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 3);
    const named_values custom_params[5] = { { "O3", 1, { 0.3 } }, { "O1", 1, { 0.1 } }, { "O2", 1, { 0.5 } }, { "BBB", 1, { 0.7 } }, { "AAA", 1, { 0.5 } } };
    tmc_error err;
    tmc_error_clear(&err);
    TMC_EXPECT_TRUE(set_custom_rate_parameter_map(state, custom_params, 5, &err) != TMC_OK);
    for (int e = 0; e < 5; ++e)
    {
      tmc_error_clear(&err);
      TMC_EXPECT_TRUE(is_state_error(set_custom_rate_parameter_map(state, &custom_params[e], 1, &err), &err, TMC_STATE_INVALID_CUSTOM_RATE_PARAM_COUNT_MULTIGRID));
    }
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(setting_custom_rate_parameter_with_invalid_label_throws_exception)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(five_labels, 5, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 1);
    const named_values custom_params[5] = { { "O3", 1, { 0.3 } }, { "O1", 1, { 0.1 } }, { "O2", 1, { 0.5 } }, { "CCC", 1, { 0.7 } }, { "AAA", 1, { 0.5 } } };
    tmc_error err;
    tmc_error_clear(&err);
    TMC_EXPECT_TRUE(set_custom_rate_parameter_map(state, custom_params, 5, &err) != TMC_OK);
    tmc_error_clear(&err);
    TMC_EXPECT_TRUE(is_state_error(set_custom_rate_parameter_map(state, &custom_params[3], 1, &err), &err, TMC_STATE_UNKNOWN_RATE_CONSTANT_PARAMETER));
    tmc_error_clear(&err);
    TMC_EXPECT_TRUE(is_state_error(tmc_state_set_custom_rate_parameter(state, "CCC", 0.7, &err), &err, TMC_STATE_UNKNOWN_RATE_CONSTANT_PARAMETER));
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(set_custom_rate_parameter)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(five_labels, 5, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 1);
    tmc_error err;
    tmc_error_clear(&err);
    TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "O2", 42.3, &err), TMC_OK);
    TMC_EXPECT_NEAR(custom_parameter(state, 0, 1), 42.3, fabs(42.3) * 1e1 * DBL_EPSILON);
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(set_custom_rate_parameters)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(five_labels, 5, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    const int num_grid_cells = 3;
    tmc_state *state = tmc_solver_create_state(solver, num_grid_cells);
    const named_values custom_params[5] = { { "O3", 3, { 0.3, 0.33, 0.333 } },
                                            { "O1", 3, { 0.1, 0.11, 0.111 } },
                                            { "O2", 3, { 0.5, 0.55, 0.555 } },
                                            { "BBB", 3, { 0.7, 0.77, 0.777 } },
                                            { "AAA", 3, { 0.5, 0.55, 0.555 } } };
    const double custom_params_in_order[15] = { 0.1, 0.5, 0.3, 0.5, 0.7, 0.11, 0.55, 0.33, 0.55, 0.77, 0.111, 0.555, 0.333, 0.555, 0.777 };
    tmc_error err;
    tmc_error_clear(&err);
    TMC_EXPECT_EQ(set_custom_rate_parameter_map(state, custom_params, 5, &err), TMC_OK);
    int idx = 0;
    for (int i = 0; i < num_grid_cells; ++i)
      for (int j = 0; j < 5; ++j, ++idx)
        TMC_EXPECT_EQ(custom_parameter(state, i, j), custom_params_in_order[idx]);
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

static void expect_unsafe_parameters(tmc_state *state, int n_cells)
{
  const double tol_scale = 1e1 * DBL_EPSILON;
  for (int i = 0; i < n_cells; ++i)
  {
    TMC_EXPECT_NEAR(custom_parameter(state, i, 0), 0.1, 0.1 * tol_scale);
    TMC_EXPECT_NEAR(custom_parameter(state, i, 1), 0.2, 0.2 * tol_scale);
    TMC_EXPECT_NEAR(custom_parameter(state, i, 2), 0.3, 0.3 * tol_scale);
    TMC_EXPECT_NEAR(custom_parameter(state, i, 3), 0.4, 0.4 * tol_scale);
    TMC_EXPECT_EQ(custom_parameter(state, i, 4), 0.5);
  }
}

TMC_TEST(unsafely_set_custom_rate_parameter_one_cell)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(five_labels, 5, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 1);
    const double parameters[1][5] = { { 0.1, 0.2, 0.3, 0.4, 0.5 } };
    tmc_error err;
    tmc_error_clear(&err);
    TMC_EXPECT_EQ(tmc_state_unsafely_set_custom_rate_parameters(state, &parameters[0][0], 1, 5, &err), TMC_OK);
    expect_unsafe_parameters(state, 1);
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(unsafely_set_custom_rate_parameter_multi_cell)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(five_labels, 5, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    const int num_grid_cells = 3;
    tmc_state *state = tmc_solver_create_state(solver, num_grid_cells);
    const double parameters[3][5] = { { 0.1, 0.2, 0.3, 0.4, 0.5 }, { 0.1, 0.2, 0.3, 0.4, 0.5 }, { 0.1, 0.2, 0.3, 0.4, 0.5 } };
    tmc_error err;
    tmc_error_clear(&err);
    TMC_EXPECT_EQ(tmc_state_unsafely_set_custom_rate_parameters(state, &parameters[0][0], num_grid_cells, 5, &err), TMC_OK);
    expect_unsafe_parameters(state, num_grid_cells);
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(unsafely_set_custom_rate_parameter_catches_too_few_grid_cells)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(five_labels, 5, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 2);
    const double parameters[1][5] = { { 0.1, 0.2, 0.3, 0.4, 0.5 } };
    tmc_error err;
    tmc_error_clear(&err);
    TMC_EXPECT_TRUE(is_state_error(tmc_state_unsafely_set_custom_rate_parameters(state, &parameters[0][0], 1, 5, &err), &err, TMC_STATE_INVALID_CUSTOM_RATE_PARAM_COUNT_MULTIGRID));
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(unsafely_set_custom_rate_parameter_catches_too_parameters)
{
  for (int w = 0; w < 2; ++w)
  {
    tmc_solver *solver = state_solver(five_labels, 5, widths[w]);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, 2);
    tmc_error err;
    /* MICM's case: one row of three values for two cells (the row count is checked first) */
    const double one_row[1][3] = { { 0.1, 0.2, 0.3 } };
    tmc_error_clear(&err);
    TMC_EXPECT_TRUE(tmc_state_unsafely_set_custom_rate_parameters(state, &one_row[0][0], 1, 3, &err) != TMC_OK);
    /* and the parameter count alone: two rows of three values for five parameters */
    const double two_rows[2][3] = { { 0.1, 0.2, 0.3 }, { 0.1, 0.2, 0.3 } };
    tmc_error_clear(&err);
    TMC_EXPECT_TRUE(is_state_error(tmc_state_unsafely_set_custom_rate_parameters(state, &two_rows[0][0], 2, 3, &err), &err, TMC_STATE_INVALID_CUSTOM_RATE_PARAM_COUNT));
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("state:State.DefaultConstructor", default_constructor),
  TMC_CASE("state:State.Constructor", constructor),
  TMC_CASE("state:State.CopyConstructor", copy_constructor),
  TMC_CASE("state:State.CopyAssignmentOperator", copy_assignment_operator),
  TMC_CASE("state:State.MoveConstructor", move_constructor),
  TMC_CASE("state:State.MoveAssignmentOperator", move_assignment_operator),
  TMC_CASE("state:State.SettingSingleConcentrationWithInvalidArgumentsThowsException", setting_single_concentration_with_invalid_arguments_throws_exception),
  TMC_CASE("state:State.SetSingleConcentration", set_single_concentration),
  TMC_CASE("state:State.SetConcentrationByElementSingleValue", set_concentration_by_element_single_value),
  TMC_CASE("state:State.SetConcentrationByElementVector", set_concentration_by_element_vector),
  TMC_CASE("state:State.SettingConcentrationsWithInvalidArguementsThrowsException", setting_concentrations_with_invalid_arguements_throws_exception),
  TMC_CASE("state:State.SetConcentrations", set_concentrations),
  TMC_CASE("state:State.SetStateWithSquareBracketOperator", set_state_with_square_bracket_operator),
  TMC_CASE("state:State.SetMultiCellStateWithSquareBrackets", set_multi_cell_state_with_square_brackets),
  TMC_CASE("state:State.SettingCustomRateParameterWithInvalidVectorSizeThrowsException", setting_custom_rate_parameter_with_invalid_vector_size_throws_exception),
  TMC_CASE("state:State.SettingCustomRateParameterWithInvalidLabelThrowsException", setting_custom_rate_parameter_with_invalid_label_throws_exception),
  TMC_CASE("state:State.SetCustomRateParameter", set_custom_rate_parameter),
  TMC_CASE("state:State.SetCustomRateParameters", set_custom_rate_parameters),
  TMC_CASE("state:State.UnsafelySetCustomRateParameterOneCell", unsafely_set_custom_rate_parameter_one_cell),
  TMC_CASE("state:State.UnsafelySetCustomRateParameterMultiCell", unsafely_set_custom_rate_parameter_multi_cell),
  TMC_CASE("state:State.UnsafelySetCustomRateParameterCatchesTooFewGridCells", unsafely_set_custom_rate_parameter_catches_too_few_grid_cells),
  TMC_CASE("state:State.UnsafelySetCustomRateParameterCatchesTooParameters", unsafely_set_custom_rate_parameter_catches_too_parameters),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_state(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
