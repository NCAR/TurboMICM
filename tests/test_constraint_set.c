/* MICM test/unit/constraint/test_constraint_set.cpp with
 * test_constraint_set_policy.hpp.
 *
 * MICM builds a ConstraintSet from constraints and a variable map. The C
 * library keeps the constraint set inside a solver, so each test builds one
 * (cst_constraint_set: the species in variable-map order, no kinetics) and
 * evaluates its constraints with the set's evaluators. The tests' sparse
 * Jacobian, the diagonal plus the set's non-zero elements, is the solver's
 * pattern, which holds both (and the LU fill-in, which no test reads). The
 * policy templates become functions of the vector width L: 1 for Matrix with
 * the standard sparse ordering, L for VectorMatrix<L> with
 * SparseMatrixVectorOrdering<L>. Tolerances are MICM's for double. */
#include "constraint_support.h"

/* Every constraint here: unit stoichiometry, K_HLC_ref = 3.3e-2, delta_H = -24000 */
static tmc_constraint equilibrium(const char *name, const char *algebraic, const char *const *reactants, int n_reactants, const char *const *products, int n_products)
{
  tmc_stoich_species r[4], p[4];
  for (int i = 0; i < n_reactants; ++i)
    r[i] = cst_stoich(reactants[i], 1.0);
  for (int i = 0; i < n_products; ++i)
    p[i] = cst_stoich(products[i], 1.0);
  tmc_species a = tmc_species_create(algebraic);
  return cst_equilibrium(name, &a, r, n_reactants, p, n_products, 3.3e-2, -24000.0);
}

/* A + B <-> AB, AB algebraic, over the variables A, B, AB */
static tmc_constraint a_b_eq(void)
{
  return equilibrium("A_B_eq", "AB", (const char *[]){ "A", "B" }, 2, (const char *[]){ "AB" }, 1);
}

static const char *const a_b_ab[] = { "A", "B", "AB" };

static bool contains(const int *ids, int n, int id)
{
  for (int i = 0; i < n; ++i)
    if (ids[i] == id)
      return true;
  return false;
}

/* ---- the policy templates ------------------------------------------------ */

static void construction(int L)
{
  tmc_constraint c = a_b_eq();
  tmc_error err;
  tmc_solver *s = cst_constraint_set(a_b_ab, 3, &c, 1, L, &err);
  CST_REQUIRE_BUILT(s, err);
  TMC_EXPECT_EQ(tmc_solver_num_constraints(s), 1);
done:
  tmc_solver_free(s);
  tmc_constraint_free(&c);
}

static void replace_state_rows_maps_to_algebraic_species(int L)
{
  const char *const species[] = { "A", "B", "C" };
  tmc_constraint c = equilibrium("B_C_eq", "C", (const char *[]){ "B" }, 1, (const char *[]){ "C" }, 1);
  tmc_index_pair *pairs = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set(species, 3, &c, 1, L, &err);
  CST_REQUIRE_BUILT(s, err);
  TMC_EXPECT_EQ(tmc_solver_num_constraints(s), 1);
  const int *ids;
  int n_ids = tmc_solver_algebraic_variables(s, &ids);
  TMC_EXPECT_EQ(n_ids, 1);
  TMC_EXPECT_TRUE(contains(ids, n_ids, 2)); /* C row is algebraic */
  int n;
  pairs = cst_nonzero_elements(s, &n);
  TMC_EXPECT_EQ(n, 2);
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 2, 1)); /* row C, col B */
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 2, 2)); /* row C, col C */
done:
  free(pairs);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
}

static void nonzero_jacobian_elements(int L)
{
  tmc_constraint c = a_b_eq();
  tmc_index_pair *pairs = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set(a_b_ab, 3, &c, 1, L, &err);
  CST_REQUIRE_BUILT(s, err);
  int n;
  pairs = cst_nonzero_elements(s, &n);
  /* AB (index 2) is algebraic: dependencies A, B, AB, plus the diagonal */
  TMC_EXPECT_EQ(n, 3);
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 2, 0));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 2, 1));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 2, 2));
done:
  free(pairs);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
}

static void multiple_constraints(int L)
{
  const char *const species[] = { "A", "B", "AB", "C", "D" };
  tmc_constraint c[2] = { a_b_eq(), equilibrium("C_D_eq", "D", (const char *[]){ "C" }, 1, (const char *[]){ "D" }, 1) };
  tmc_index_pair *pairs = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set(species, 5, c, 2, L, &err);
  CST_REQUIRE_BUILT(s, err);
  TMC_EXPECT_EQ(tmc_solver_num_constraints(s), 2);
  int n;
  pairs = cst_nonzero_elements(s, &n);
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 2, 0));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 2, 1));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 2, 2));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 4, 3));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 4, 4));
  TMC_EXPECT_EQ(n, 5);
done:
  free(pairs);
  tmc_solver_free(s);
  tmc_constraint_free(&c[0]);
  tmc_constraint_free(&c[1]);
}

static void add_forcing_terms(int L)
{
  tmc_constraint c = a_b_eq();
  tmc_dense *state = NULL, *forcing = NULL, *parameters = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set(a_b_ab, 3, &c, 1, L, &err);
  CST_REQUIRE_BUILT(s, err);
  /* SetConstraintFunctions: the solver maps the constraint parameter to a column */
  CST_REQUIRE(cst_parameter_index(s, "A_B_eq") >= 0);
  state = tmc_dense_create(2, 3, L, 0.0);
  tmc_dense_set_row(state, 0, (const double[]){ 0.2, 0.4, 0.6 }); /* away from equilibrium */
  tmc_dense_set_row(state, 1, (const double[]){ 0.1, 0.3, 0.7 });
  forcing = tmc_dense_create(2, 3, L, 0.0);
  parameters = cst_parameters(s, 2, 3.3e-2);
  cst_add_forcing_terms(s, state, parameters, forcing);
  /* G = K_eq [A][B] - [AB] replaces the AB row */
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 0, 2), -0.59736, 1e-5);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 1, 2), -0.69901, 1e-5);
done:
  tmc_dense_free(state);
  tmc_dense_free(forcing);
  tmc_dense_free(parameters);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
}

static void subtract_jacobian_terms(int L)
{
  tmc_constraint c = a_b_eq();
  tmc_dense *state = NULL, *parameters = NULL;
  tmc_sparse *jacobian = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set(a_b_ab, 3, &c, 1, L, &err);
  CST_REQUIRE_BUILT(s, err);
  CST_REQUIRE(cst_parameter_index(s, "A_B_eq") >= 0);
  jacobian = cst_jacobian(s, 1);
  state = tmc_dense_create(1, 3, L, 0.0);
  tmc_dense_set_row(state, 0, (const double[]){ 0.01, 0.02, 0.05 });
  parameters = cst_parameters(s, 1, 3.3e-2);
  cst_subtract_jacobian_terms(s, state, parameters, jacobian);
  /* dG/dA = K_eq [B], dG/dB = K_eq [A], dG/dAB = -1, subtracted in row 2 */
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 2, 0), -0.00066, 1e-10);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 2, 1), -0.00033, 1e-10);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 2, 2), 1.0, 1e-10);
done:
  tmc_dense_free(state);
  tmc_dense_free(parameters);
  tmc_sparse_free(jacobian);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
}

static void empty_constraint_set(int L)
{
  /* the default-constructed set: a solver without constraints */
  tmc_dense *state = NULL, *forcing = NULL, *parameters = NULL;
  tmc_index_pair *pairs = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set((const char *[]){ "A", "B", "C" }, 3, NULL, 0, L, &err);
  CST_REQUIRE_BUILT(s, err);
  TMC_EXPECT_EQ(tmc_solver_num_constraints(s), 0);
  int n;
  pairs = cst_nonzero_elements(s, &n);
  TMC_EXPECT_EQ(n, 0);
  state = tmc_dense_create(1, 3, L, 0.0);
  tmc_dense_set_row(state, 0, (const double[]){ 0.1, 0.2, 0.3 });
  forcing = tmc_dense_create(1, 4, L, 1.0);
  parameters = tmc_dense_create(1, 0, L, 0.0);
  cst_add_forcing_terms(s, state, parameters, forcing);
  for (int j = 0; j < 4; ++j)
    TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(forcing, 0, j), 1.0);
done:
  free(pairs);
  tmc_dense_free(state);
  tmc_dense_free(forcing);
  tmc_dense_free(parameters);
  tmc_solver_free(s);
}

static void unknown_species_throws(int L)
{
  tmc_constraint c = equilibrium("invalid", "XY", (const char *[]){ "X", "Y" }, 2, (const char *[]){ "XY" }, 1);
  tmc_error err;
  tmc_solver *s = cst_constraint_set((const char *[]){ "A", "B" }, 2, &c, 1, L, &err);
  TMC_EXPECT_TRUE(s == NULL);
  CST_EXPECT_ERROR(err.category, err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_UNKNOWN_SPECIES);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
}

static void three_d_state_one_constraint(int L)
{
  const double K_eq = 3.3e-2;
  const char *const species[] = { "X", "Y", "Z" };
  tmc_constraint c = equilibrium("X_Y_eq", "Y", (const char *[]){ "X" }, 1, (const char *[]){ "Y" }, 1);
  tmc_index_pair *pairs = NULL;
  tmc_dense *state = NULL, *forcing = NULL, *parameters = NULL;
  tmc_sparse *jacobian = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set(species, 3, &c, 1, L, &err);
  CST_REQUIRE_BUILT(s, err);
  TMC_EXPECT_EQ(tmc_solver_num_constraints(s), 1);
  /* Y (index 1) is algebraic */
  int n;
  pairs = cst_nonzero_elements(s, &n);
  TMC_EXPECT_EQ(n, 2);
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 1, 0));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 1, 1));
  TMC_EXPECT_FALSE(cst_has_pair(pairs, n, 1, 2)); /* Z not involved */
  jacobian = cst_jacobian(s, 2);
  CST_REQUIRE(cst_parameter_index(s, "X_Y_eq") >= 0);
  state = tmc_dense_create(2, 3, L, 0.0);
  tmc_dense_set_row(state, 0, (const double[]){ 10.0, 0.2, 0.5 });  /* away from equilibrium */
  tmc_dense_set_row(state, 1, (const double[]){ 10.0, 0.33, 0.3 }); /* at equilibrium, Y = K_eq X */
  forcing = tmc_dense_create(2, 3, L, 0.0);
  parameters = cst_parameters(s, 2, 3.3e-2);
  cst_add_forcing_terms(s, state, parameters, forcing);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 0, 1), 0.13, 1e-10);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 1, 1), 0.0, 1e-10);
  cst_subtract_jacobian_terms(s, state, parameters, jacobian);
  /* dG/dX = K_eq, dG/dY = -1, subtracted in row 1 */
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 1, 0), -K_eq, 1e-10);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 1, 1), 1.0, 1e-10);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 1, 1, 0), -K_eq, 1e-10);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 1, 1, 1), 1.0, 1e-10);
  /* the Z row is untouched */
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 2, 2), 0.0, 1e-10);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 1, 2, 2), 0.0, 1e-10);
done:
  free(pairs);
  tmc_dense_free(state);
  tmc_dense_free(forcing);
  tmc_dense_free(parameters);
  tmc_sparse_free(jacobian);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
}

static void four_d_state_two_constraints(int L)
{
  const double K_eq1 = 3.3e-2, K_eq2 = 3.3e-2;
  const char *const species[] = { "A", "B", "C", "D" };
  /* A <-> B with B algebraic (row 1); C + D <-> A with A algebraic (row 0) */
  tmc_constraint c[2] = {
    equilibrium("A_B_eq", "B", (const char *[]){ "A" }, 1, (const char *[]){ "B" }, 1),
    equilibrium("CD_A_eq", "A", (const char *[]){ "C", "D" }, 2, (const char *[]){ "A" }, 1),
  };
  tmc_index_pair *pairs = NULL;
  tmc_dense *state = NULL, *forcing = NULL, *parameters = NULL;
  tmc_sparse *jacobian = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set(species, 4, c, 2, L, &err);
  CST_REQUIRE_BUILT(s, err);
  TMC_EXPECT_EQ(tmc_solver_num_constraints(s), 2);
  int n;
  pairs = cst_nonzero_elements(s, &n);
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 1, 0));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 1, 1));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 0, 2));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 0, 3));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 0, 0));
  TMC_EXPECT_EQ(n, 5);
  jacobian = cst_jacobian(s, 3);
  /* SetConstraintFunctions: the solver maps each constraint parameter to its own column */
  const int k1 = cst_parameter_index(s, "A_B_eq"), k2 = cst_parameter_index(s, "CD_A_eq");
  CST_REQUIRE(k1 >= 0 && k2 >= 0 && k1 != k2);
  state = tmc_dense_create(3, 4, L, 0.0);
  tmc_dense_set_row(state, 0, (const double[]){ 0.33, 0.01089, 10.0, 1.0 }); /* both satisfied */
  tmc_dense_set_row(state, 1, (const double[]){ 0.33, 0.01089, 5.0, 1.0 });  /* the first satisfied */
  tmc_dense_set_row(state, 2, (const double[]){ 0.33, 0.02, 5.0, 1.0 });     /* neither */
  forcing = tmc_dense_create(3, 4, L, 0.0);
  parameters = cst_parameters(s, 3, 0.0);
  for (int i = 0; i < 3; ++i)
  {
    tmc_dense_set(parameters, i, k1, 3.3e-2);
    tmc_dense_set(parameters, i, k2, 3.3e-2);
  }
  cst_add_forcing_terms(s, state, parameters, forcing);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 0, 1), 0.0, 1e-5);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 0, 0), 0.0, 1e-10);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 1, 1), 0.0, 1e-5);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 1, 0), -0.165, 1e-10);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 2, 1), -0.00911, 1e-5);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 2, 0), -0.165, 1e-10);
  cst_subtract_jacobian_terms(s, state, parameters, jacobian);
  /* row 1: dG1/dA = K_eq1, dG1/dB = -1; row 0: dG2/dC = K_eq2 [D], dG2/dD = K_eq2 [C], dG2/dA = -1 */
  for (int i = 0; i < 3; ++i)
  {
    TMC_EXPECT_NEAR(CST_JAC(jacobian, i, 1, 0), -K_eq1, 1e-10);
    TMC_EXPECT_NEAR(CST_JAC(jacobian, i, 1, 1), 1.0, 1e-10);
    TMC_EXPECT_NEAR(CST_JAC(jacobian, i, 0, 2), -K_eq2 * tmc_dense_get(state, i, 3), 1e-10);
    TMC_EXPECT_NEAR(CST_JAC(jacobian, i, 0, 3), -K_eq2 * tmc_dense_get(state, i, 2), 1e-10);
    TMC_EXPECT_NEAR(CST_JAC(jacobian, i, 0, 0), 1.0, 1e-10);
  }
done:
  free(pairs);
  tmc_dense_free(state);
  tmc_dense_free(forcing);
  tmc_dense_free(parameters);
  tmc_sparse_free(jacobian);
  tmc_solver_free(s);
  tmc_constraint_free(&c[0]);
  tmc_constraint_free(&c[1]);
}

static void coupled_constraints_shared_species(int L)
{
  const double K_eq1 = 3.3e-2, K_eq2 = 3.3e-2;
  const char *const species[] = { "A", "B", "C" };
  /* both depend on A: A <-> B (B algebraic), A <-> C (C algebraic) */
  tmc_constraint c[2] = {
    equilibrium("A_B_eq", "B", (const char *[]){ "A" }, 1, (const char *[]){ "B" }, 1),
    equilibrium("A_C_eq", "C", (const char *[]){ "A" }, 1, (const char *[]){ "C" }, 1),
  };
  tmc_index_pair *pairs = NULL;
  tmc_dense *state = NULL, *forcing = NULL, *parameters = NULL;
  tmc_sparse *jacobian = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set(species, 3, c, 2, L, &err);
  CST_REQUIRE_BUILT(s, err);
  TMC_EXPECT_EQ(tmc_solver_num_constraints(s), 2);
  int n;
  pairs = cst_nonzero_elements(s, &n);
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 1, 0));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 1, 1));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 2, 0));
  TMC_EXPECT_TRUE(cst_has_pair(pairs, n, 2, 2));
  TMC_EXPECT_EQ(n, 4);
  jacobian = cst_jacobian(s, 1);
  const int k1 = cst_parameter_index(s, "A_B_eq"), k2 = cst_parameter_index(s, "A_C_eq");
  CST_REQUIRE(k1 >= 0 && k2 >= 0 && k1 != k2);
  state = tmc_dense_create(1, 3, L, 0.0);
  tmc_dense_set_row(state, 0, (const double[]){ 0.1, 0.0033, 0.0033 }); /* both at equilibrium */
  forcing = tmc_dense_create(1, 3, L, 0.0);
  parameters = cst_parameters(s, 1, 0.0);
  tmc_dense_set(parameters, 0, k1, 3.3e-2);
  tmc_dense_set(parameters, 0, k2, 3.3e-2);
  cst_add_forcing_terms(s, state, parameters, forcing);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 0, 1), 0.0, 1e-10);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 0, 2), 0.0, 1e-10);
  cst_subtract_jacobian_terms(s, state, parameters, jacobian);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 1, 0), -K_eq1, 1e-10);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 1, 1), 1.0, 1e-10);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 2, 0), -K_eq2, 1e-10);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 2, 2), 1.0, 1e-10);
done:
  free(pairs);
  tmc_dense_free(state);
  tmc_dense_free(forcing);
  tmc_dense_free(parameters);
  tmc_sparse_free(jacobian);
  tmc_solver_free(s);
  tmc_constraint_free(&c[0]);
  tmc_constraint_free(&c[1]);
}

static void vectorized_matrices_respect_grid_cell_indexing(int L)
{
  tmc_constraint c = a_b_eq();
  tmc_dense *state = NULL, *forcing = NULL, *parameters = NULL;
  tmc_sparse *jacobian = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set(a_b_ab, 3, &c, 1, L, &err);
  CST_REQUIRE_BUILT(s, err);
  jacobian = cst_jacobian(s, 3);
  CST_REQUIRE(cst_parameter_index(s, "A_B_eq") >= 0);
  state = tmc_dense_create(3, 3, L, 0.0);
  tmc_dense_set_row(state, 0, (const double[]){ 0.01, 0.02, 0.05 });
  tmc_dense_set_row(state, 1, (const double[]){ 0.03, 0.01, 0.2 });
  tmc_dense_set_row(state, 2, (const double[]){ 0.001, 0.002, 0.004 });
  forcing = tmc_dense_create(3, 3, L, 0.0);
  parameters = cst_parameters(s, 3, 3.3e-2);
  cst_add_forcing_terms(s, state, parameters, forcing);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 0, 2), 3.3e-2 * 0.01 * 0.02 - 0.05, 1e-9);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 1, 2), 3.3e-2 * 0.03 * 0.01 - 0.2, 1e-9);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 2, 2), 3.3e-2 * 0.001 * 0.002 - 0.004, 1e-9);
  cst_subtract_jacobian_terms(s, state, parameters, jacobian);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 2, 0), -(3.3e-2 * 0.02), 1e-12);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 2, 1), -(3.3e-2 * 0.01), 1e-12);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 2, 2), 1.0, 1e-12);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 1, 2, 0), -(3.3e-2 * 0.01), 1e-12);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 1, 2, 1), -(3.3e-2 * 0.03), 1e-12);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 1, 2, 2), 1.0, 1e-12);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 2, 2, 0), -(3.3e-2 * 0.002), 1e-12);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 2, 2, 1), -(3.3e-2 * 0.001), 1e-12);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 2, 2, 2), 1.0, 1e-12);
done:
  tmc_dense_free(state);
  tmc_dense_free(forcing);
  tmc_dense_free(parameters);
  tmc_sparse_free(jacobian);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
}

/* ---- the cases (Matrix and the standard ordering are L = 1) -------------- */

TMC_TEST(test_construction)
{
  construction(1);
  return TMC_TEST_PASS;
}

TMC_TEST(test_replace_state_rows)
{
  replace_state_rows_maps_to_algebraic_species(1);
  return TMC_TEST_PASS;
}

TMC_TEST(test_nonzero_jacobian_elements)
{
  nonzero_jacobian_elements(1);
  return TMC_TEST_PASS;
}

TMC_TEST(test_multiple_constraints)
{
  multiple_constraints(1);
  return TMC_TEST_PASS;
}

TMC_TEST(test_add_forcing_terms)
{
  add_forcing_terms(1);
  return TMC_TEST_PASS;
}

TMC_TEST(test_subtract_jacobian_terms)
{
  subtract_jacobian_terms(1);
  return TMC_TEST_PASS;
}

TMC_TEST(test_empty_constraint_set)
{
  empty_constraint_set(1);
  return TMC_TEST_PASS;
}

TMC_TEST(test_unknown_species_throws)
{
  unknown_species_throws(1);
  return TMC_TEST_PASS;
}

TMC_TEST(test_three_d_state_one_constraint)
{
  three_d_state_one_constraint(1);
  return TMC_TEST_PASS;
}

TMC_TEST(test_four_d_state_two_constraints)
{
  four_d_state_two_constraints(1);
  return TMC_TEST_PASS;
}

TMC_TEST(test_coupled_constraints_shared_species)
{
  coupled_constraints_shared_species(1);
  return TMC_TEST_PASS;
}

TMC_TEST(test_vectorized_matrices)
{
  vectorized_matrices_respect_grid_cell_indexing(4);
  return TMC_TEST_PASS;
}

static void vector_matrix(int L)
{
  construction(L);
  nonzero_jacobian_elements(L);
  add_forcing_terms(L);
  subtract_jacobian_terms(L);
}

TMC_TEST(test_vector_matrix_1)
{
  vector_matrix(1);
  return TMC_TEST_PASS;
}

TMC_TEST(test_vector_matrix_2)
{
  vector_matrix(2);
  return TMC_TEST_PASS;
}

TMC_TEST(test_vector_matrix_3)
{
  vector_matrix(3);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("constraint_set:ConstraintSet.Construction", test_construction),
  TMC_CASE("constraint_set:ConstraintSet.ReplaceStateRowsMapsToAlgebraicSpecies", test_replace_state_rows),
  TMC_CASE("constraint_set:ConstraintSet.NonZeroJacobianElements", test_nonzero_jacobian_elements),
  TMC_CASE("constraint_set:ConstraintSet.MultipleConstraints", test_multiple_constraints),
  TMC_CASE("constraint_set:ConstraintSet.AddForcingTerms", test_add_forcing_terms),
  TMC_CASE("constraint_set:ConstraintSet.SubtractJacobianTerms", test_subtract_jacobian_terms),
  TMC_CASE("constraint_set:ConstraintSet.EmptyConstraintSet", test_empty_constraint_set),
  TMC_CASE("constraint_set:ConstraintSet.UnknownSpeciesThrows", test_unknown_species_throws),
  TMC_CASE("constraint_set:ConstraintSet.ThreeDStateOneConstraint", test_three_d_state_one_constraint),
  TMC_CASE("constraint_set:ConstraintSet.FourDStateTwoConstraints", test_four_d_state_two_constraints),
  TMC_CASE("constraint_set:ConstraintSet.CoupledConstraintsSharedSpecies", test_coupled_constraints_shared_species),
  TMC_CASE("constraint_set:ConstraintSet.VectorizedMatricesRespectGridCellIndexing", test_vectorized_matrices),
  TMC_CASE("constraint_set:ConstraintSet.VectorMatrix1", test_vector_matrix_1),
  TMC_CASE("constraint_set:ConstraintSet.VectorMatrix2", test_vector_matrix_2),
  TMC_CASE("constraint_set:ConstraintSet.VectorMatrix3", test_vector_matrix_3),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_constraint_set(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
