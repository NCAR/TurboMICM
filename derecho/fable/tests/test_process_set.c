/* MICM test/unit/process/test_process_set.cpp (with test_process_set_policy.hpp)
 *
 * MICM's ProcessSet is the solver's forcing and Jacobian of the reactions: a
 * built solver's tmc_solver_add_forcing_terms and
 * tmc_solver_subtract_jacobian_terms, with the Jacobian in the solver's LU
 * pattern (tmc_lu_matrix). MICM's State is a tmc_state of that solver; the
 * cases that set rate constants by hand write the state's rate-constant
 * matrix, one column per reaction in the solver's order. MICM instantiates
 * the cases for Matrix and VectorMatrix<1..4>; here those are solver widths. */
#include "harness.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- helpers ---- */

typedef struct
{
  unsigned long long state;
} rng;

static double uniform(rng *g)
{
  g->state = g->state * 6364136223846793005ULL + 1442695040888963407ULL;
  return ((double)(g->state >> 11) + 0.5) / 9007199254740992.0;
}

static int uniform_int(rng *g, int lo, int hi) /* inclusive */
{
  return lo + (int)(uniform(g) * (double)(hi - lo + 1));
}

static double lognormal(rng *g, double mean, double sigma)
{
  double u1 = uniform(g), u2 = uniform(g);
  return exp(mean + sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2));
}

static tmc_arrhenius_parameters arrhenius_12_2(void)
{
  tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
  p.A = 12.2;
  p.C = 300.0;
  return p;
}

static tmc_reaction make_reaction(const tmc_species *const *reactants, int n_reactants, const tmc_species *const *products,
                                  const double *yields, int n_products, tmc_rate_constant rc, const tmc_phase *phase)
{
  tmc_reaction r = tmc_reaction_create();
  for (int i = 0; i < n_reactants; ++i)
    tmc_reaction_add_reactant(&r, reactants[i]);
  for (int i = 0; i < n_products; ++i)
    tmc_reaction_add_product(&r, products[i], yields[i]);
  tmc_reaction_set_rate_constant(&r, rc);
  tmc_reaction_set_phase(&r, phase);
  return r;
}

static int var(const tmc_state *st, const char *name)
{
  int i = tmc_state_variable_index(st, name);
  if (i < 0)
    tmc_test_fail(__FILE__, __LINE__, "unknown variable %s", name);
  return i;
}

/* the Jacobian element (row, col) of a block, which must be in the pattern */
static double element(const tmc_sparse *m, int block, int row, int col)
{
  if (row < 0 || col < 0 || tmc_sparse_is_zero(m, row, col))
  {
    tmc_test_fail(__FILE__, __LINE__, "element (%d, %d) is not in the Jacobian pattern", row, col);
    return NAN;
  }
  return tmc_sparse_get(m, block, row, col);
}

/* ---- ProcessSet.Matrix and ProcessSet.VectorMatrix ---- */

typedef struct
{
  const char *row, *col;
  int cell;
  double value;
} jacobian_value;

static void process_set_case(int L)
{
  tmc_species foo = tmc_species_create("foo"), bar = tmc_species_create("bar"), baz = tmc_species_create("baz");
  tmc_species quz = tmc_species_create("quz"), quuz = tmc_species_create("quuz"), qux = tmc_species_create("qux");
  tmc_species corge = tmc_species_create("corge");
  tmc_species_set_parameterization(&qux, 0.0, 0.0, 0.0, 0.72); /* MICM: { .c_rho_ = 0.72 } */

  tmc_phase gas_phase = tmc_phase_create("gas");
  const tmc_species *in_phase[] = { &foo, &bar, &baz, &quz, &quuz, &corge };
  for (int i = 0; i < 6; ++i)
    tmc_phase_add_species(&gas_phase, in_phase[i]);

  tmc_rate_constant arrhenius = tmc_rate_arrhenius(arrhenius_12_2());
  tmc_reaction r[4];
  {
    const tmc_species *re[] = { &foo, &baz }, *pr[] = { &bar, &quuz };
    const double y[] = { 1.0, 2.4 };
    r[0] = make_reaction(re, 2, pr, y, 2, arrhenius, &gas_phase);
  }
  {
    const tmc_species *re[] = { &bar, &qux }, *pr[] = { &foo, &quz };
    const double y[] = { 1.0, 1.4 };
    r[1] = make_reaction(re, 2, pr, y, 2, arrhenius, &gas_phase);
  }
  {
    const tmc_species *re[] = { &quz };
    r[2] = make_reaction(re, 1, NULL, NULL, 0, arrhenius, &gas_phase);
  }
  {
    const tmc_species *re[] = { &baz, &qux }, *pr[] = { &bar, &quz };
    const double y[] = { 1.0, 2.5 };
    r[3] = make_reaction(re, 2, pr, y, 2, arrhenius, &gas_phase);
  }

  tmc_system system = tmc_system_create(&gas_phase);
  tmc_builder *b = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, r, 4);
  tmc_builder_set_width(b, L);
  tmc_builder_set_reorder_state(b, false); /* MICM's state lists the species in declaration order */

  /* MICM's SpeciesUsed: every reactant and product, so all species but corge.
   * The builder's unused-species check names exactly the species that are not. */
  tmc_builder_set_ignore_unused_species(b, false);
  tmc_error err;
  tmc_solver *s = tmc_builder_build(b, &err);
  TMC_EXPECT_TRUE(s == NULL);
  TMC_EXPECT_EQ(err.category, TMC_ERROR_SOLVER);
  TMC_EXPECT_EQ(err.code, TMC_SOLVER_UNUSED_SPECIES);
  TMC_EXPECT_TRUE(strstr(err.message, "'corge'") != NULL);
  const char *used[] = { "'foo'", "'bar'", "'baz'", "'quz'", "'quuz'", "'qux'" };
  for (int i = 0; i < 6; ++i)
    TMC_EXPECT_TRUE(strstr(err.message, used[i]) == NULL);
  tmc_solver_free(s);

  tmc_builder_set_ignore_unused_species(b, true);
  s = tmc_builder_build(b, &err);
  TMC_EXPECT_TRUE(s != NULL);
  if (!s)
    goto done;
  {
    int n_rxn;
    const tmc_reaction *sorted = tmc_solver_reactions(s, &n_rxn);
    TMC_EXPECT_EQ(n_rxn, 4);
    const char *first_reactant[] = { "foo", "bar", "quz", "baz" };
    for (int i = 0; i < 4 && i < n_rxn; ++i)
      TMC_EXPECT_STREQ(sorted[i].reactants[0].name, first_reactant[i]);

    tmc_state *state = tmc_solver_create_state(s, 2);
    TMC_EXPECT_EQ(tmc_state_num_cells(state), 2);
    TMC_EXPECT_EQ(tmc_state_num_variables(state), 6);
    tmc_dense *variables = tmc_state_variables(state);
    TMC_EXPECT_EQ(variables->n_rows, 2);
    TMC_EXPECT_EQ(variables->n_cols, 6);
    const int i_foo = var(state, "foo"), i_bar = var(state, "bar"), i_baz = var(state, "baz");
    const int i_quz = var(state, "quz"), i_quuz = var(state, "quuz"), i_corge = var(state, "corge");
    if (i_foo < 0 || i_bar < 0 || i_baz < 0 || i_quz < 0 || i_quuz < 0 || i_corge < 0)
    {
      tmc_state_free(state);
      goto done;
    }
    const double values[2][6] = { { 0.1, 0.2, 0.3, 0.4, 0.5, 0.0 }, { 1.1, 1.2, 1.3, 1.4, 1.5, 0.0 } };
    const int ids[6] = { i_foo, i_bar, i_baz, i_quz, i_quuz, i_corge };
    for (int c = 0; c < 2; ++c)
      for (int i = 0; i < 6; ++i)
        tmc_state_set(state, c, ids[i], values[c][i]);
    tmc_conditions *conditions = tmc_state_conditions(state);
    conditions[0].air_density = 70.0;
    conditions[1].air_density = 80.0;

    /* "the rate constants will have been calculated and combined with
     * parameterized species before calculating forcing terms": the solver
     * multiplies reactions 2 and 4 by qux = 0.72 * air density */
    conditions[0].temperature = conditions[1].temperature = 298.15;
    tmc_solver_update_state_parameters(s, state);
    tmc_dense *rate_constants = tmc_state_rate_constants(state);
    for (int c = 0; c < 2; ++c)
    {
      double k = 12.2 * exp(300.0 / 298.15);
      double qux_value = 0.72 * conditions[c].air_density;
      TMC_EXPECT_REL(tmc_dense_get(rate_constants, c, 0), k, 1e-14);
      TMC_EXPECT_REL(tmc_dense_get(rate_constants, c, 1), k * qux_value, 1e-14);
      TMC_EXPECT_REL(tmc_dense_get(rate_constants, c, 2), k, 1e-14);
      TMC_EXPECT_REL(tmc_dense_get(rate_constants, c, 3), k * qux_value, 1e-14);
    }
    /* MICM's test then sets them by hand */
    const double k_hand[2][4] = { { 10.0, 20.0 * 70.0 * 0.72, 30.0, 40.0 * 70.0 * 0.72 },
                                  { 110.0, 120.0 * 80.0 * 0.72, 130.0, 140.0 * 80.0 * 0.72 } };
    for (int c = 0; c < 2; ++c)
      for (int i = 0; i < 4; ++i)
        tmc_dense_set(rate_constants, c, i, k_hand[c][i]);

    tmc_dense *forcing = tmc_dense_create(2, 6, L, 1000.0);
    tmc_solver_add_forcing_terms(s, state, variables, forcing);
#define F(c, v) tmc_dense_get(forcing, c, v)
    TMC_EXPECT_DOUBLE_EQ(F(0, i_foo), 1000.0 - 10.0 * 0.1 * 0.3 + 20.0 * 70.0 * 0.72 * 0.2);
    TMC_EXPECT_DOUBLE_EQ(F(1, i_foo), 1000.0 - 110.0 * 1.1 * 1.3 + 120.0 * 80.0 * 0.72 * 1.2);
    TMC_EXPECT_DOUBLE_EQ(F(0, i_bar), 1000.0 + 10.0 * 0.1 * 0.3 - 20.0 * 0.2 * 70.0 * 0.72 + 40.0 * 70.0 * 0.72 * 0.3);
    const double forcing_1_1 = 1000.0 + 110.0 * 1.1 * 1.3 - 120.0 * 1.2 * 80.0 * 0.72 + 140.0 * 80.0 * 0.72 * 1.3;
    TMC_EXPECT_NEAR(F(1, i_bar), forcing_1_1, forcing_1_1 * 1.0e-12);
    TMC_EXPECT_DOUBLE_EQ(F(0, i_baz), 1000.0 - 10.0 * 0.1 * 0.3 - 40.0 * 70.0 * 0.72 * 0.3);
    TMC_EXPECT_DOUBLE_EQ(F(1, i_baz), 1000.0 - 110.0 * 1.1 * 1.3 - 140.0 * 80.0 * 0.72 * 1.3);
    TMC_EXPECT_DOUBLE_EQ(F(0, i_quz), 1000.0 + 20.0 * 70.0 * 0.72 * 0.2 * 1.4 - 30.0 * 0.4 + 40.0 * 70.0 * 0.72 * 2.5 * 0.3);
    TMC_EXPECT_DOUBLE_EQ(F(1, i_quz), 1000.0 + 120.0 * 80.0 * 0.72 * 1.2 * 1.4 - 130.0 * 1.4 + 140.0 * 80.0 * 0.72 * 2.5 * 1.3);
    TMC_EXPECT_DOUBLE_EQ(F(0, i_quuz), 1000.0 + 10.0 * 0.1 * 0.3 * 2.4);
    TMC_EXPECT_DOUBLE_EQ(F(1, i_quuz), 1000.0 + 110.0 * 1.1 * 1.3 * 2.4);
    TMC_EXPECT_EQ(F(0, i_corge), 1000.0);
    TMC_EXPECT_EQ(F(1, i_corge), 1000.0);
#undef F
    tmc_dense_free(forcing);

    /* MICM's NonZeroJacobianElements: the elements the Jacobian terms write.
     *   ---- foo  bar  baz  quz  quuz
     *   foo   0    1    2    -    -
     *   bar   3    4    5    -    -
     *   baz   6    -    7    -    -
     *   quz   -    8    9   10    -
     *   quuz 11    -   12    -    -
     * The solver keeps them in its LU pattern, which adds the diagonal and
     * fill-in; running the Jacobian terms with NaN rate constants turns
     * exactly the written elements into NaN. */
    const char *nonzero[13][2] = { { "foo", "foo" }, { "foo", "bar" }, { "foo", "baz" }, { "bar", "foo" }, { "bar", "bar" },
                                   { "bar", "baz" }, { "baz", "foo" }, { "baz", "baz" }, { "quz", "bar" }, { "quz", "baz" },
                                   { "quz", "quz" }, { "quuz", "foo" }, { "quuz", "baz" } };
    {
      tmc_state *probe = tmc_solver_create_state(s, 2);
      tmc_dense_fill(tmc_state_variables(probe), 1.0);
      tmc_dense_fill(tmc_state_rate_constants(probe), NAN);
      tmc_sparse *written = tmc_lu_matrix(tmc_solver_lu(s), 2, 0.0);
      tmc_solver_subtract_jacobian_terms(s, probe, tmc_state_variables(probe), written);
      for (int row = 0; row < 6; ++row)
        for (int col = 0; col < 6; ++col)
        {
          bool expected = false;
          for (int e = 0; e < 13; ++e)
            expected = expected || (row == var(state, nonzero[e][0]) && col == var(state, nonzero[e][1]));
          for (int c = 0; c < 2; ++c)
          {
            bool is_written = !tmc_sparse_is_zero(written, row, col) && isnan(tmc_sparse_get(written, c, row, col));
            if (is_written != expected)
              tmc_test_fail(__FILE__, __LINE__, "L = %d: element (%s, %s) %s written", L, tmc_state_variable_name(state, row),
                            tmc_state_variable_name(state, col), is_written ? "is" : "is not");
          }
        }
      tmc_sparse_free(written);
      tmc_state_free(probe);
    }

    tmc_sparse *jacobian = tmc_lu_matrix(tmc_solver_lu(s), 2, 100.0);
    tmc_solver_subtract_jacobian_terms(s, state, variables, jacobian);
    const jacobian_value expected_values[] = {
      { "foo", "foo", 0, 100.0 + 10.0 * 0.3 },
      { "foo", "foo", 1, 100.0 + 110.0 * 1.3 },
      { "foo", "bar", 0, 100.0 - 20.0 * 70.0 * 0.72 },
      { "foo", "bar", 1, 100.0 - 120.0 * 80.0 * 0.72 },
      { "foo", "baz", 0, 100.0 + 10.0 * 0.1 },
      { "foo", "baz", 1, 100.0 + 110.0 * 1.1 },
      { "bar", "foo", 0, 100.0 - 10.0 * 0.3 },
      { "bar", "foo", 1, 100.0 - 110.0 * 1.3 },
      { "bar", "bar", 0, 100.0 + 20.0 * 70.0 * 0.72 },
      { "bar", "bar", 1, 100.0 + 120.0 * 80.0 * 0.72 },
      { "bar", "baz", 0, 100.0 - 10.0 * 0.1 - 40.0 * 70.0 * 0.72 },
      { "bar", "baz", 1, 100.0 - 110.0 * 1.1 - 140.0 * 80.0 * 0.72 },
      { "baz", "foo", 0, 100.0 + 10.0 * 0.3 },
      { "baz", "foo", 1, 100.0 + 110.0 * 1.3 },
      { "baz", "baz", 0, 100.0 + 10.0 * 0.1 + 40.0 * 70.0 * 0.72 },
      { "baz", "baz", 1, 100.0 + 110.0 * 1.1 + 140.0 * 80.0 * 0.72 },
      { "quz", "bar", 0, 100.0 - 1.4 * 20.0 * 70.0 * 0.72 },
      { "quz", "bar", 1, 100.0 - 1.4 * 120.0 * 80.0 * 0.72 },
      { "quz", "baz", 0, 100.0 - 2.5 * 40.0 * 70.0 * 0.72 },
      { "quz", "baz", 1, 100.0 - 2.5 * 140.0 * 80.0 * 0.72 },
      { "quz", "quz", 0, 100.0 + 30.0 },
      { "quz", "quz", 1, 100.0 + 130.0 },
      { "quuz", "foo", 0, 100.0 - 2.4 * 10.0 * 0.3 },
      { "quuz", "foo", 1, 100.0 - 2.4 * 110.0 * 1.3 },
      { "quuz", "baz", 0, 100.0 - 2.4 * 10.0 * 0.1 },
      { "quuz", "baz", 1, 100.0 - 2.4 * 110.0 * 1.1 },
    };
    for (size_t e = 0; e < sizeof expected_values / sizeof expected_values[0]; ++e)
    {
      const jacobian_value *v = &expected_values[e];
      double got = element(jacobian, v->cell, var(state, v->row), var(state, v->col));
      if (!tmc_test_double_eq(got, v->value))
        tmc_test_fail(__FILE__, __LINE__, "L = %d: J[%d][%s][%s] = %.17g, expected %.17g", L, v->cell, v->row, v->col, got, v->value);
    }
    /* the other elements of the pattern (diagonal, fill-in) are untouched */
    const tmc_sparse *pattern = tmc_solver_lu_pattern(s);
    for (int e = 0; e < pattern->nnz; ++e)
    {
      int row = pattern->row_of[e], col = pattern->col_of[e];
      bool listed = false;
      for (int k = 0; k < 13; ++k)
        listed = listed || (row == var(state, nonzero[k][0]) && col == var(state, nonzero[k][1]));
      if (!listed)
        for (int c = 0; c < 2; ++c)
          TMC_EXPECT_EQ(tmc_sparse_get(jacobian, c, row, col), 100.0);
    }
    tmc_sparse_free(jacobian);
    tmc_state_free(state);
  }
done:
  tmc_solver_free(s);
  tmc_builder_free(b);
  tmc_system_free(&system);
  for (int i = 0; i < 4; ++i)
    tmc_reaction_free(&r[i]);
  tmc_phase_free(&gas_phase);
}

TMC_TEST(process_set_matrix)
{
  process_set_case(1);
  return TMC_TEST_PASS;
}

TMC_TEST(process_set_vector_matrix)
{
  for (int L = 1; L <= 4; ++L)
    process_set_case(L);
  return TMC_TEST_PASS;
}

/* ---- RandomProcessSet.Matrix ---- */

/* MICM builds random reactions (0-3 reactants, 0-10 products of yield 1.2)
 * and evaluates the forcing with lognormal(-2, 4) variables and rate
 * constants; it checks only that the evaluation runs. Here the forcing is
 * also compared with a sum computed in the test: each species gains
 * yield * rate per appearance as a product and loses the rate per appearance
 * as a reactant. A second pass with values near 1 makes every term count. */
static void random_system(int n_cells, int n_reactions, int n_species, unsigned long long seed)
{
  rng g = { seed };
  tmc_species *species = malloc(sizeof *species * (size_t)n_species);
  tmc_phase gas_phase = tmc_phase_create("gas");
  for (int i = 0; i < n_species; ++i)
  {
    char name[16];
    snprintf(name, sizeof name, "%d", i);
    species[i] = tmc_species_create(name);
    tmc_phase_add_species(&gas_phase, &species[i]);
  }
  tmc_reaction *reactions = malloc(sizeof *reactions * (size_t)n_reactions);
  for (int r = 0; r < n_reactions; ++r)
  {
    reactions[r] = tmc_reaction_create();
    int n_react = uniform_int(&g, 0, 3);
    for (int i = 0; i < n_react; ++i)
      tmc_reaction_add_reactant(&reactions[r], &species[uniform_int(&g, 0, n_species - 1)]);
    int n_product = uniform_int(&g, 0, 10);
    for (int i = 0; i < n_product; ++i)
      tmc_reaction_add_product(&reactions[r], &species[uniform_int(&g, 0, n_species - 1)], 1.2);
    tmc_reaction_set_rate_constant(&reactions[r], tmc_rate_arrhenius(arrhenius_12_2()));
    tmc_reaction_set_phase(&reactions[r], &gas_phase);
  }
  tmc_system system = tmc_system_create(&gas_phase);
  tmc_builder *b = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, reactions, n_reactions);
  tmc_builder_set_width(b, 1);
  tmc_builder_set_ignore_unused_species(b, true); /* random reactions may leave species out */
  tmc_error err;
  tmc_solver *s = tmc_builder_build(b, &err);
  TMC_EXPECT_TRUE(s != NULL);
  if (s)
  {
    tmc_state *state = tmc_solver_create_state(s, n_cells);
    tmc_dense *variables = tmc_state_variables(state);
    tmc_dense *rate_constants = tmc_state_rate_constants(state);
    TMC_EXPECT_EQ(rate_constants->n_cols, n_reactions);
    int n_rxn;
    const tmc_reaction *sorted = tmc_solver_reactions(s, &n_rxn);
    double *expected = malloc(sizeof(double) * (size_t)n_species);
    double *scale = malloc(sizeof(double) * (size_t)n_species);
    for (int pass = 0; pass < 2; ++pass)
    {
      for (size_t i = 0; i < variables->size; ++i)
        variables->data[i] = pass == 0 ? lognormal(&g, -2.0, 4.0) : 0.5 + uniform(&g);
      for (size_t i = 0; i < rate_constants->size; ++i)
        rate_constants->data[i] = pass == 0 ? lognormal(&g, -2.0, 4.0) : 0.5 + uniform(&g);
      tmc_dense *forcing = tmc_dense_create(n_cells, n_species, 1, 1000.0);
      tmc_solver_add_forcing_terms(s, state, variables, forcing);
      for (int c = 0; c < n_cells; ++c)
      {
        for (int v = 0; v < n_species; ++v)
          expected[v] = scale[v] = 1000.0;
        for (int r = 0; r < n_rxn; ++r)
        {
          double rate = tmc_dense_get(rate_constants, c, r);
          for (int i = 0; i < sorted[r].n_reactants; ++i)
            rate *= tmc_state_get(state, c, tmc_state_variable_index(state, sorted[r].reactants[i].name));
          for (int i = 0; i < sorted[r].n_reactants; ++i)
          {
            int v = tmc_state_variable_index(state, sorted[r].reactants[i].name);
            expected[v] -= rate;
            scale[v] += fabs(rate);
          }
          for (int i = 0; i < sorted[r].n_products; ++i)
          {
            int v = tmc_state_variable_index(state, sorted[r].products[i].species.name);
            expected[v] += sorted[r].products[i].coefficient * rate;
            scale[v] += fabs(sorted[r].products[i].coefficient * rate);
          }
        }
        for (int v = 0; v < n_species; ++v)
        {
          double got = tmc_dense_get(forcing, c, v);
          if (!(fabs(got - expected[v]) <= 1e-12 * scale[v]))
          {
            tmc_test_fail(__FILE__, __LINE__, "%d cells, %d reactions, %d species, pass %d: forcing[%d][%d] = %.17g, expected %.17g", n_cells,
                          n_reactions, n_species, pass, c, v, got, expected[v]);
            c = n_cells; /* one report per system */
            break;
          }
        }
      }
      tmc_dense_free(forcing);
    }
    free(expected);
    free(scale);
    tmc_state_free(state);
  }
  tmc_solver_free(s);
  tmc_builder_free(b);
  tmc_system_free(&system);
  tmc_reactions_free(reactions, n_reactions);
  tmc_phase_free(&gas_phase);
  free(species);
}

TMC_TEST(random_process_set_matrix)
{
  random_system(200, 50, 40, 1);
  random_system(300, 30, 20, 2);
  random_system(400, 100, 80, 3);
  return TMC_TEST_PASS;
}

/* ---- ProcessSetAlgebraicVariables.CudaMatrix ---- */

/* MICM marks B and D algebraic on the process set directly
 * (SetAlgebraicVariableIds); the solver takes algebraic variables from its
 * constraints, so two linear constraints make B and D algebraic. Their
 * dependencies keep (B, A), (B, B), (D, A), and (D, B) in the Jacobian
 * pattern, which the kinetic Jacobian terms must leave untouched; the forcing
 * and Jacobian terms of the reactions do not include the constraints. */
static void algebraic_masking_case(int L)
{
  tmc_species A = tmc_species_create("A"), B = tmc_species_create("B"), C = tmc_species_create("C"), D = tmc_species_create("D");
  tmc_phase gas_phase = tmc_phase_create("gas");
  tmc_phase_add_species(&gas_phase, &A);
  tmc_phase_add_species(&gas_phase, &B);
  tmc_phase_add_species(&gas_phase, &C);
  tmc_phase_add_species(&gas_phase, &D);
  const tmc_species *re[] = { &A, &B }, *pr[] = { &C, &D };
  const double y[] = { 1.0, 1.0 };
  tmc_reaction r1 = make_reaction(re, 2, pr, y, 2, tmc_rate_arrhenius(arrhenius_12_2()), &gas_phase);

  tmc_error err;
  tmc_constraint b_constraint, d_constraint;
  const tmc_stoich_species b_terms[] = { { A, 1.0 }, { B, 1.0 } };
  const tmc_stoich_species d_terms[] = { { A, 1.0 }, { B, 1.0 }, { D, 1.0 } };
  TMC_EXPECT_EQ(tmc_linear_constraint(&b_constraint, "B balance", &B, b_terms, 2, 3.0, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_linear_constraint(&d_constraint, "D balance", &D, d_terms, 3, 10.0, &err), TMC_OK);

  tmc_system system = tmc_system_create(&gas_phase);
  tmc_builder *b = tmc_builder_create_rosenbrock(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, &r1, 1);
  tmc_builder_add_constraint(b, &b_constraint);
  tmc_builder_add_constraint(b, &d_constraint);
  tmc_builder_set_width(b, L);
  tmc_builder_set_reorder_state(b, false);
  tmc_solver *s = tmc_builder_build(b, &err);
  TMC_EXPECT_TRUE(s != NULL);
  if (s)
  {
    tmc_state *state = tmc_solver_create_state(s, 2);
    const int iA = var(state, "A"), iB = var(state, "B"), iC = var(state, "C"), iD = var(state, "D");
    /* the algebraic variables are B (index 1) and D (index 3) */
    const int *algebraic;
    int n_algebraic = tmc_solver_algebraic_variables(s, &algebraic);
    TMC_EXPECT_EQ(n_algebraic, 2);
    if (n_algebraic == 2)
    {
      TMC_EXPECT_EQ(algebraic[0], iB);
      TMC_EXPECT_EQ(algebraic[1], iD);
    }
    TMC_EXPECT_EQ(iB, 1);
    TMC_EXPECT_EQ(iD, 3);

    const double values[2][4] = { { 1.0, 2.0, 3.0, 4.0 }, { 5.0, 6.0, 7.0, 8.0 } };
    const int ids[4] = { iA, iB, iC, iD };
    for (int c = 0; c < 2; ++c)
      for (int i = 0; i < 4; ++i)
        tmc_state_set(state, c, ids[i], values[c][i]);
    tmc_dense *rate_constants = tmc_state_rate_constants(state);
    TMC_EXPECT_EQ(rate_constants->n_cols, 1);
    tmc_dense_set(rate_constants, 0, 0, 10.0);
    tmc_dense_set(rate_constants, 1, 0, 20.0);

    /* rate = k [A] [B]: 20 in cell 0, 600 in cell 1 */
    tmc_dense *forcing = tmc_dense_create(2, 4, L, 1000.0);
    tmc_solver_add_forcing_terms(s, state, tmc_state_variables(state), forcing);
    TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing, 0, iA), 1000.0 - 20.0);
    TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing, 1, iA), 1000.0 - 600.0);
    TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing, 0, iB), 1000.0);
    TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing, 1, iB), 1000.0);
    TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing, 0, iC), 1000.0 + 20.0);
    TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing, 1, iC), 1000.0 + 600.0);
    TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing, 0, iD), 1000.0);
    TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing, 1, iD), 1000.0);
    tmc_dense_free(forcing);

    /* d rate / d A = k [B]: 20 and 120; d rate / d B = k [A]: 10 and 100 */
    tmc_sparse *jacobian = tmc_lu_matrix(tmc_solver_lu(s), 2, 500.0);
    tmc_solver_subtract_jacobian_terms(s, state, tmc_state_variables(state), jacobian);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 0, iA, iA), 500.0 + 20.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 1, iA, iA), 500.0 + 120.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 0, iA, iB), 500.0 + 10.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 1, iA, iB), 500.0 + 100.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 0, iB, iA), 500.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 1, iB, iA), 500.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 0, iB, iB), 500.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 1, iB, iB), 500.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 0, iC, iA), 500.0 - 20.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 1, iC, iA), 500.0 - 120.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 0, iC, iB), 500.0 - 10.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 1, iC, iB), 500.0 - 100.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 0, iD, iA), 500.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 1, iD, iA), 500.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 0, iD, iB), 500.0);
    TMC_EXPECT_DOUBLE_EQ(element(jacobian, 1, iD, iB), 500.0);
    tmc_sparse_free(jacobian);
    tmc_state_free(state);
  }
  tmc_solver_free(s);
  tmc_builder_free(b);
  tmc_constraint_free(&b_constraint);
  tmc_constraint_free(&d_constraint);
  tmc_system_free(&system);
  tmc_reaction_free(&r1);
  tmc_phase_free(&gas_phase);
}

TMC_TEST(algebraic_masking)
{
  for (int L = 1; L <= 4; ++L)
    algebraic_masking_case(L);
  return TMC_TEST_PASS;
}

/* ---- ProcessSetFiniteDifferenceJacobian.Matrix ---- */

/* MICM's FiniteDifferenceJacobian (micm/util/jacobian_verification.hpp):
 * central differences with the step 1e-8 max(1, |x|), or forward differences
 * where x - h < 0. fd[(block * n + row) * n + col] = d forcing_row / d x_col. */
static double *finite_difference_jacobian(const tmc_solver *s, const tmc_state *state, const tmc_dense *base, int n)
{
  const int n_blocks = base->n_rows;
  double *fd = calloc((size_t)n_blocks * (size_t)n * (size_t)n, sizeof(double));
  tmc_dense *plus = tmc_dense_copy(base), *minus = tmc_dense_copy(base);
  tmc_dense *f_plus = tmc_dense_create(n_blocks, n, base->L, 0.0);
  tmc_dense *f_minus = tmc_dense_create(n_blocks, n, base->L, 0.0);
  tmc_dense *f_base = tmc_dense_create(n_blocks, n, base->L, 0.0);
  for (int col = 0; col < n; ++col)
  {
    tmc_dense_copy_into(plus, base);
    tmc_dense_copy_into(minus, base);
    for (int block = 0; block < n_blocks; ++block)
    {
      double x = tmc_dense_get(base, block, col);
      double h = 1.0e-8 * fmax(1.0, fabs(x));
      tmc_dense_set(plus, block, col, x + h);
      tmc_dense_set(minus, block, col, x - h);
    }
    tmc_dense_fill(f_plus, 0.0);
    tmc_dense_fill(f_minus, 0.0);
    tmc_solver_add_forcing_terms(s, state, plus, f_plus);
    tmc_solver_add_forcing_terms(s, state, minus, f_minus);
    for (int block = 0; block < n_blocks; ++block)
    {
      double x = tmc_dense_get(base, block, col);
      double h = 1.0e-8 * fmax(1.0, fabs(x));
      if (x - h >= 0.0)
        for (int row = 0; row < n; ++row)
          fd[((size_t)block * n + row) * n + col] = (tmc_dense_get(f_plus, block, row) - tmc_dense_get(f_minus, block, row)) / (2.0 * h);
      else
      {
        tmc_dense_fill(f_base, 0.0);
        tmc_solver_add_forcing_terms(s, state, base, f_base);
        for (int row = 0; row < n; ++row)
          fd[((size_t)block * n + row) * n + col] = (tmc_dense_get(f_plus, block, row) - tmc_dense_get(f_base, block, row)) / h;
      }
    }
  }
  tmc_dense_free(plus);
  tmc_dense_free(minus);
  tmc_dense_free(f_plus);
  tmc_dense_free(f_minus);
  tmc_dense_free(f_base);
  return fd;
}

TMC_TEST(finite_difference_jacobian_matrix)
{
  /* r1: A + B -> 2 C (k1); r2: C -> A (k2) */
  tmc_species A = tmc_species_create("A"), B = tmc_species_create("B"), C = tmc_species_create("C");
  tmc_phase gas_phase = tmc_phase_create("gas");
  tmc_phase_add_species(&gas_phase, &A);
  tmc_phase_add_species(&gas_phase, &B);
  tmc_phase_add_species(&gas_phase, &C);
  const int num_species = 3;
  tmc_arrhenius_parameters rate = tmc_arrhenius_defaults();
  rate.A = 1.0;
  tmc_reaction r[2];
  {
    const tmc_species *re[] = { &A, &B }, *pr[] = { &C };
    const double y[] = { 2.0 };
    r[0] = make_reaction(re, 2, pr, y, 1, tmc_rate_arrhenius(rate), &gas_phase);
  }
  {
    const tmc_species *re[] = { &C }, *pr[] = { &A };
    const double y[] = { 1.0 };
    r[1] = make_reaction(re, 1, pr, y, 1, tmc_rate_arrhenius(rate), &gas_phase);
  }
  tmc_system system = tmc_system_create(&gas_phase);
  tmc_builder *b = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, r, 2);
  tmc_builder_set_width(b, 1);
  tmc_builder_set_reorder_state(b, false);
  tmc_error err;
  tmc_solver *s = tmc_builder_build(b, &err);
  TMC_ASSERT_TRUE(s != NULL);
  tmc_state *state = tmc_solver_create_state(s, 2);
  const int ids[3] = { var(state, "A"), var(state, "B"), var(state, "C") };
  const double values[2][3] = { { 1.0, 2.0, 3.0 }, { 0.5, 4.0, 1.5 } };
  for (int c = 0; c < 2; ++c)
    for (int i = 0; i < 3; ++i)
      tmc_state_set(state, c, ids[i], values[c][i]);
  tmc_dense *rate_constants = tmc_state_rate_constants(state);
  tmc_dense_set(rate_constants, 0, 0, 5.0);
  tmc_dense_set(rate_constants, 0, 1, 10.0);
  tmc_dense_set(rate_constants, 1, 0, 8.0);
  tmc_dense_set(rate_constants, 1, 1, 12.0);

  tmc_sparse *analytical = tmc_lu_matrix(tmc_solver_lu(s), 2, 0.0);
  tmc_solver_subtract_jacobian_terms(s, state, tmc_state_variables(state), analytical);
  double *fd = finite_difference_jacobian(s, state, tmc_state_variables(state), num_species);

  /* MICM's CompareJacobianToFiniteDifference: |-J - fd| <= 1e-7 + 1e-7 max(|J|, |fd|)
   * over the pattern, and CheckJacobianSparsityCompleteness: |fd| <= 1e-6 off it */
  for (int block = 0; block < 2; ++block)
    for (int row = 0; row < num_species; ++row)
      for (int col = 0; col < num_species; ++col)
      {
        double fd_val = fd[((size_t)block * num_species + row) * num_species + col];
        if (tmc_sparse_is_zero(analytical, row, col))
        {
          if (fabs(fd_val) > 1.0e-6)
            tmc_test_fail(__FILE__, __LINE__, "undeclared non-zero at block %d row %d col %d: fd %.17g", block, row, col, fd_val);
          continue;
        }
        double analytical_val = -tmc_sparse_get(analytical, block, row, col);
        double abs_error = fabs(analytical_val - fd_val);
        double scale = fmax(fabs(analytical_val), fabs(fd_val));
        if (abs_error > 1.0e-7 + 1.0e-7 * scale)
          tmc_test_fail(__FILE__, __LINE__, "block %d row %d col %d: analytical %.17g, fd %.17g", block, row, col, analytical_val, fd_val);
      }
  free(fd);
  tmc_sparse_free(analytical);
  tmc_state_free(state);
  tmc_solver_free(s);
  tmc_builder_free(b);
  tmc_system_free(&system);
  tmc_reaction_free(&r[0]);
  tmc_reaction_free(&r[1]);
  tmc_phase_free(&gas_phase);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("process_set:ProcessSet.Matrix", process_set_matrix),
  TMC_CASE("process_set:ProcessSet.VectorMatrix", process_set_vector_matrix),
  TMC_CASE("process_set:RandomProcessSet.Matrix", random_process_set_matrix),
  TMC_CASE("process_set:ProcessSetAlgebraicVariables.CudaMatrix", algebraic_masking),
  TMC_CASE("process_set:ProcessSetFiniteDifferenceJacobian.Matrix", finite_difference_jacobian_matrix),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_process_set(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
