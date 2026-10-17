/* MICM test/unit/solver/test_rosenbrock.cpp with
 * test_rosenbrock_solver_policy.hpp. MICM's StandardBuilder is width L = 1;
 * its VectorBuilder<L> is tmc_builder_set_width(L). The solver internals the
 * tests call (NormalizedError, AlphaMinusJacobian) are in tmc_solver_ext.h. */
#include "harness.h"
#include "tmc_solver_ext.h"

static tmc_rate_constant arrhenius(double A, double B, double C)
{
  tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
  p.A = A;
  p.B = B;
  p.C = C;
  return tmc_rate_arrhenius(p);
}

static tmc_solver *build(tmc_builder *builder, tmc_system *system, tmc_reaction *reactions, int n_reactions, tmc_phase *gas)
{
  tmc_builder_set_system(builder, system);
  tmc_builder_set_reactions(builder, reactions, n_reactions);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_solver *solver = tmc_builder_build(builder, &err);
  if (!solver)
    fprintf(stderr, "  build failed: %s\n", err.message);
  tmc_builder_free(builder);
  for (int i = 0; i < n_reactions; ++i)
    tmc_reaction_free(&reactions[i]);
  tmc_system_free(system);
  tmc_phase_free(gas);
  return solver;
}

/* MICM's GetSolver: five species, three reactions, no reordering.
 *   ---- foo  bar  baz  quz  quuz
 *   foo   0    1    2    -    -
 *   bar   3    4    5    -    -
 *   baz   6    -    7    -    -
 *   quz   -    8    -    9    -
 *   quuz 10    -   11    -    12 */
static tmc_solver *get_solver(tmc_rosenbrock_parameters parameters, int L, tmc_lu_algorithm algorithm)
{
  tmc_species foo = tmc_species_create("foo"), bar = tmc_species_create("bar"), baz = tmc_species_create("baz");
  tmc_species quz = tmc_species_create("quz"), quuz = tmc_species_create("quuz");
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &foo);
  tmc_phase_add_species(&gas, &bar);
  tmc_phase_add_species(&gas, &baz);
  tmc_phase_add_species(&gas, &quz);
  tmc_phase_add_species(&gas, &quuz);
  tmc_system system = tmc_system_create(&gas);

  tmc_reaction r[3];
  r[0] = tmc_reaction_create();
  tmc_reaction_add_reactant(&r[0], &foo);
  tmc_reaction_add_reactant(&r[0], &baz);
  tmc_reaction_add_product(&r[0], &bar, 1);
  tmc_reaction_add_product(&r[0], &quuz, 2.4);
  tmc_reaction_set_phase(&r[0], &gas);
  tmc_reaction_set_rate_constant(&r[0], arrhenius(2.0e-11, 0, 110));
  r[1] = tmc_reaction_create();
  tmc_reaction_add_reactant(&r[1], &bar);
  tmc_reaction_add_product(&r[1], &foo, 1);
  tmc_reaction_add_product(&r[1], &quz, 1.4);
  tmc_reaction_set_phase(&r[1], &gas);
  tmc_reaction_set_rate_constant(&r[1], arrhenius(1.0e-6, 0, 0));
  r[2] = tmc_reaction_create();
  tmc_reaction_add_reactant(&r[2], &quz);
  tmc_reaction_set_phase(&r[2], &gas);
  tmc_reaction_set_rate_constant(&r[2], arrhenius(3.5e-6, 0, 0));

  tmc_builder *builder = tmc_builder_create_rosenbrock(parameters);
  tmc_builder_set_reorder_state(builder, false);
  tmc_builder_set_width(builder, L);
  tmc_builder_set_lu_algorithm(builder, algorithm);
  return build(builder, &system, r, 3, &gas);
}

/* ---- AlphaMinusJacobian -------------------------------------------------- */

static void test_alpha_minus_jacobian(int L, int n_cells)
{
  tmc_solver *solver = get_solver(tmc_rosenbrock_three_stage(), L, TMC_LU_MOZART);
  TMC_EXPECT_TRUE(solver != NULL);
  if (!solver)
    return;
  tmc_state *state = tmc_solver_create_state(solver, n_cells);
  tmc_sparse *J = tmc_state_jacobian(state);
  TMC_EXPECT_EQ(J->n_blocks, n_cells);
  TMC_EXPECT_EQ(J->n, 5); /* rows; the blocks are square */
  TMC_EXPECT_EQ(J->L, L);
  TMC_EXPECT_TRUE(J->size >= (size_t)13 * (size_t)n_cells);

  static const int elements[13][2] = { { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 0 }, { 1, 1 }, { 1, 2 }, { 2, 0 },
                                       { 2, 2 }, { 3, 1 }, { 3, 3 }, { 4, 0 }, { 4, 2 }, { 4, 4 } };
  for (int e = 0; e < 13; ++e)
    TMC_EXPECT_FALSE(tmc_sparse_is_zero(J, elements[e][0], elements[e][1]));

  /* a negative Jacobian */
  for (int c = 0; c < n_cells; ++c)
  {
    tmc_sparse_set(J, c, 0, 0, -12.2);
    tmc_sparse_set(J, c, 0, 1, -24.3 * (c + 2));
    tmc_sparse_set(J, c, 0, 2, -42.3);
    tmc_sparse_set(J, c, 1, 0, -0.43);
    tmc_sparse_set(J, c, 1, 1, -23.4);
    tmc_sparse_set(J, c, 1, 2, -83.4 / (c + 3));
    tmc_sparse_set(J, c, 2, 0, -4.74);
    tmc_sparse_set(J, c, 2, 2, -6.91);
    tmc_sparse_set(J, c, 3, 1, -59.1);
    tmc_sparse_set(J, c, 3, 3, -83.4);
    tmc_sparse_set(J, c, 4, 0, -78.5);
    tmc_sparse_set(J, c, 4, 2, -53.6);
    tmc_sparse_set(J, c, 4, 4, -1.0);
  }

  tmc_solver_alpha_minus_jacobian(solver, state, 42.042);

  for (int c = 0; c < n_cells; ++c)
  {
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 0, 0), 42.042 - 12.2);
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 0, 1), -24.3 * (c + 2));
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 0, 2), -42.3);
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 1, 0), -0.43);
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 1, 1), 42.042 - 23.4);
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 1, 2), -83.4 / (c + 3));
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 2, 0), -4.74);
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 2, 2), 42.042 - 6.91);
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 3, 1), -59.1);
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 3, 3), 42.042 - 83.4);
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 4, 0), -78.5);
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 4, 2), -53.6);
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(J, c, 4, 4), 42.042 - 1.0);
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
}

TMC_TEST(standard_alpha_minus_jacobian)
{
  for (int n_cells = 1; n_cells <= 4; ++n_cells)
    test_alpha_minus_jacobian(1, n_cells);
  return TMC_TEST_PASS;
}

TMC_TEST(vector_alpha_minus_jacobian)
{
  test_alpha_minus_jacobian(1, 1);
  test_alpha_minus_jacobian(2, 4);
  test_alpha_minus_jacobian(3, 3);
  test_alpha_minus_jacobian(4, 2);
  return TMC_TEST_PASS;
}

/* ---- tolerances ---------------------------------------------------------- */

TMC_TEST(can_set_tolerances)
{
  for (int n_cells = 1; n_cells <= 10; ++n_cells)
  {
    tmc_species foo = tmc_species_create_with_property("foo", TMC_PROPERTY_ABSOLUTE_TOLERANCE, 1.0e-07);
    tmc_species bar = tmc_species_create_with_property("bar", TMC_PROPERTY_ABSOLUTE_TOLERANCE, 1.0e-08);
    tmc_phase gas = tmc_phase_create("gas");
    tmc_phase_add_species(&gas, &foo);
    tmc_phase_add_species(&gas, &bar);
    tmc_system system = tmc_system_create(&gas);
    tmc_reaction r1 = tmc_reaction_create();
    tmc_reaction_add_reactant(&r1, &foo);
    tmc_reaction_add_product(&r1, &bar, 1);
    tmc_reaction_set_rate_constant(&r1, arrhenius(2.0e-11, 0, 110));
    tmc_reaction_set_phase(&r1, &gas);
    tmc_solver *solver = build(tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage()), &system, &r1, 1, &gas);
    TMC_ASSERT_TRUE(solver != NULL);
    tmc_state *state = tmc_solver_create_state(solver, n_cells);
    const double *absolute_tolerances = tmc_state_absolute_tolerances(state);
    TMC_EXPECT_EQ(tmc_state_num_variables(state), 2);
    TMC_EXPECT_DOUBLE_EQ(absolute_tolerances[0], 1.0e-07);
    TMC_EXPECT_DOUBLE_EQ(absolute_tolerances[1], 1.0e-08);
    tmc_state_free(state);
    tmc_solver_free(solver);
  }
  return TMC_TEST_PASS;
}

/* ---- NormalizedError ----------------------------------------------------- */

/* MICM's TestNormalizedErrorDiff: elements differ across the matrices, so the
 * root mean square changes with their size; the padded lanes hold non-zero
 * values that must not enter the norm. */
static void test_normalized_error_diff(int L, int n_cells)
{
  tmc_solver *solver = get_solver(tmc_rosenbrock_three_stage(), L, TMC_LU_MOZART);
  TMC_EXPECT_TRUE(solver != NULL);
  if (!solver)
    return;
  tmc_state *state = tmc_solver_create_state(solver, n_cells);
  const double *atol = tmc_state_absolute_tolerances(state);
  const double rtol = tmc_state_relative_tolerance(state);
  const int n_vars = tmc_state_num_variables(state);

  tmc_dense *y_old = tmc_dense_create(n_cells, n_vars, L, 7.7);
  tmc_dense *y_new = tmc_dense_create(n_cells, n_vars, L, -13.9);
  tmc_dense *errors = tmc_dense_create(n_cells, n_vars, L, 81.57);
  double expected_error = 0.0;
  for (int i = 0; i < n_cells; ++i)
    for (int j = 0; j < n_vars; ++j)
    {
      tmc_dense_set(y_old, i, j, tmc_dense_get(y_old, i, j) * (double)i + (double)j);
      tmc_dense_set(y_new, i, j, tmc_dense_get(y_new, i, j) / (double)(j + 1) - (double)i);
      tmc_dense_set(errors, i, j, tmc_dense_get(errors, i, j) / (double)(i + 7) / (double)(j + 3));
      double ymax = fmax(fabs(tmc_dense_get(y_old, i, j)), fabs(tmc_dense_get(y_new, i, j)));
      double scale = atol[j] + rtol * ymax;
      double e = tmc_dense_get(errors, i, j);
      expected_error += e * e / (scale * scale);
    }
  expected_error = fmax(sqrt(expected_error / (double)(n_cells * n_vars)), 1.0e-10);

  double computed_error = tmc_solver_normalized_error(solver, state, y_old, y_new, errors);
  double relative_error = fabs(computed_error - expected_error) / fmax(fabs(computed_error), fabs(expected_error));
  if (!(relative_error <= 1.0e-11))
    tmc_test_fail(__FILE__, __LINE__, "L = %d, %d cells: computed %.12g, expected %.12g, relative error %.3g", L, n_cells, computed_error, expected_error, relative_error);

  tmc_dense_free(y_old);
  tmc_dense_free(y_new);
  tmc_dense_free(errors);
  tmc_state_free(state);
  tmc_solver_free(solver);
}

TMC_TEST(standard_normalized_error)
{
  for (int n_cells = 1; n_cells <= 4; ++n_cells)
    test_normalized_error_diff(1, n_cells);
  return TMC_TEST_PASS;
}

TMC_TEST(vector_normalized_error)
{
  /* exact fits */
  test_normalized_error_diff(1, 1);
  test_normalized_error_diff(2, 2);
  test_normalized_error_diff(3, 3);
  test_normalized_error_diff(4, 4);
  /* inexact fits */
  test_normalized_error_diff(2, 1);
  test_normalized_error_diff(3, 2);
  test_normalized_error_diff(4, 3);
  test_normalized_error_diff(8, 5);
  test_normalized_error_diff(10, 3);
  return TMC_TEST_PASS;
}

/* MICM's TestNormalizedErrorIncludesAllVariables builds A -> B plus an
 * equilibrium constraint B <-> C, so C is an algebraic variable, checks that
 * the state has one constraint, and checks that the error norm includes every
 * variable. The library has no algebraic constraints (DAE), so only the norm
 * over the three variables can be checked; the case reports the gap. */
static void test_normalized_error_includes_all_variables(tmc_rosenbrock_parameters parameters, int L, int n_cells)
{
  tmc_species a = tmc_species_create("A"), b = tmc_species_create("B"), c = tmc_species_create("C");
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &a);
  tmc_phase_add_species(&gas, &b);
  tmc_phase_add_species(&gas, &c);
  tmc_system system = tmc_system_create(&gas);
  tmc_reaction r = tmc_reaction_create();
  tmc_reaction_add_reactant(&r, &a);
  tmc_reaction_add_product(&r, &b, 1);
  tmc_reaction_set_rate_constant(&r, arrhenius(0.5, 0.0, 0.0));
  tmc_reaction_set_phase(&r, &gas);
  tmc_builder *builder = tmc_builder_create_rosenbrock(parameters);
  tmc_builder_set_reorder_state(builder, false);
  tmc_builder_set_width(builder, L);
  /* C would be the constraint's algebraic variable; without the constraint it is unused */
  tmc_builder_set_ignore_unused_species(builder, true);
  tmc_solver *solver = build(builder, &system, &r, 1, &gas);
  TMC_EXPECT_TRUE(solver != NULL);
  if (!solver)
    return;
  tmc_state *state = tmc_solver_create_state(solver, n_cells);
  const int n_vars = tmc_state_num_variables(state);
  TMC_EXPECT_EQ(n_vars, 3);

  tmc_dense *y_old = tmc_dense_create(n_cells, n_vars, L, 0.0);
  tmc_dense *y_new = tmc_dense_create(n_cells, n_vars, L, 0.0);
  tmc_dense *errors = tmc_dense_create(n_cells, n_vars, L, 0.0);
  const double *atol = tmc_state_absolute_tolerances(state);
  const double rtol = tmc_state_relative_tolerance(state);
  double expected_error = 0.0;
  for (int i = 0; i < n_cells; ++i)
    for (int j = 0; j < n_vars; ++j)
    {
      tmc_dense_set(y_old, i, j, 1.0 + i + 0.1 * j);
      tmc_dense_set(y_new, i, j, 0.8 + 0.5 * i + 0.2 * j);
      tmc_dense_set(errors, i, j, 0.01 * (1 + i + j));
      const double ymax = fmax(fabs(tmc_dense_get(y_old, i, j)), fabs(tmc_dense_get(y_new, i, j)));
      const double scale = atol[j] + rtol * ymax;
      const double e = tmc_dense_get(errors, i, j);
      expected_error += e * e / (scale * scale);
    }
  expected_error = sqrt(expected_error / (double)(n_cells * n_vars));
  expected_error = fmax(expected_error, 1.0e-10);

  double computed_error = tmc_solver_normalized_error(solver, state, y_old, y_new, errors);
  TMC_EXPECT_NEAR(computed_error, expected_error, 1e-12);

  tmc_dense_free(y_old);
  tmc_dense_free(y_new);
  tmc_dense_free(errors);
  tmc_state_free(state);
  tmc_solver_free(solver);
}

static const char *const no_constraints = "needs algebraic constraints (an EquilibriumConstraint and the state's constraint count), which the library does not implement; the error norm over all variables passes";

TMC_TEST(standard_normalized_error_with_constraints)
{
  test_normalized_error_includes_all_variables(tmc_rosenbrock_two_stage(), 1, 2);
  test_normalized_error_includes_all_variables(tmc_rosenbrock_three_stage(), 1, 2);
  test_normalized_error_includes_all_variables(tmc_rosenbrock_four_stage(), 1, 2);
  test_normalized_error_includes_all_variables(tmc_rosenbrock_four_stage_dae(), 1, 2);
  test_normalized_error_includes_all_variables(tmc_rosenbrock_six_stage_dae(), 1, 2);
  test_normalized_error_includes_all_variables(tmc_rosenbrock_three_stage(), 1, 5);
  tmc_test_set_reason(no_constraints);
  return TMC_TEST_FAIL;
}

TMC_TEST(vector_normalized_error_with_constraints)
{
  test_normalized_error_includes_all_variables(tmc_rosenbrock_two_stage(), 4, 3);
  test_normalized_error_includes_all_variables(tmc_rosenbrock_three_stage(), 4, 3);
  test_normalized_error_includes_all_variables(tmc_rosenbrock_four_stage(), 4, 3);
  test_normalized_error_includes_all_variables(tmc_rosenbrock_four_stage_dae(), 4, 3);
  test_normalized_error_includes_all_variables(tmc_rosenbrock_six_stage_dae(), 4, 3);
  test_normalized_error_includes_all_variables(tmc_rosenbrock_three_stage(), 4, 5);
  tmc_test_set_reason(no_constraints);
  return TMC_TEST_FAIL;
}

/* ---- rejected steps -------------------------------------------------------- */

typedef struct step_counts
{
  tmc_solver_state state;
  long steps, accepted, rejected;
} step_counts;

/* MICM's Rosenbrock loop with its default, non-in-place linear solver
 * (CpuSolverBuilder): the Jacobian is formed once per step, a rejected attempt
 * shifts its diagonal by the change in alpha only, and the shifted matrix is
 * factored into separate L and U matrices (Doolittle). The library's solver is
 * MICM's in-place variant: it rebuilds alpha I - J for every attempt. Forcing,
 * Jacobian, error norm, and step-size rule are the library's. */
static step_counts reference_standard_solve(const tmc_solver *solver, tmc_state *state, double time_step)
{
  const tmc_rosenbrock_parameters *p = tmc_solver_rosenbrock_parameters(solver);
  const int n_cells = tmc_state_num_cells(state), n_vars = tmc_state_num_variables(state), L = tmc_state_width(state);
  step_counts counts = { TMC_RUNNING, 0, 0, 0 };
  tmc_error err;
  tmc_error_clear(&err);
  tmc_lu *lu = tmc_lu_create(tmc_solver_lu_pattern(solver), TMC_LU_DOOLITTLE, &err);
  if (!lu)
  {
    counts.state = TMC_NOT_YET_CALLED;
    return counts;
  }
  tmc_sparse *jacobian = tmc_lu_matrix(lu, n_cells, 0.0);
  tmc_sparse *lower = tmc_lu_lower_matrix(lu, n_cells, 0.0);
  tmc_sparse *upper = tmc_lu_upper_matrix(lu, n_cells, 0.0);
  tmc_dense *Y = tmc_dense_copy(tmc_state_variables(state));
  tmc_dense *Ynew = tmc_dense_create(n_cells, n_vars, L, 0.0);
  tmc_dense *Yerror = tmc_dense_create(n_cells, n_vars, L, 0.0);
  tmc_dense *initial_forcing = tmc_dense_create(n_cells, n_vars, L, 0.0);
  tmc_dense *K[TMC_MAX_STAGES];
  for (int i = 0; i < p->stages; ++i)
    K[i] = tmc_dense_create(n_cells, n_vars, L, 0.0);

  const double h_min = p->h_min == 0.0 ? 1.0e-15 * time_step : p->h_min;
  const double h_max = p->h_max == 0.0 ? time_step : fmin(time_step, p->h_max);
  const double h_start = p->h_start == 0.0 ? 1.0e-6 * time_step : fmin(h_max, p->h_start);
  double H = fmin(fmax(h_min, fabs(h_start)), fabs(h_max));
  double present_time = 0.0;
  bool reject_last_h = false, reject_more_h = false;

  while ((present_time - time_step + p->round_off) <= 0 && counts.state == TMC_RUNNING)
  {
    if (counts.steps > p->max_number_of_steps)
    {
      counts.state = TMC_CONVERGENCE_EXCEEDED_MAX_STEPS;
      break;
    }
    if (((present_time + 0.1 * H) == present_time) || (H <= fmin(p->round_off, h_min)))
    {
      counts.state = TMC_STEP_SIZE_TOO_SMALL;
      break;
    }
    H = fmin(H, fabs(time_step - present_time));
    tmc_dense_fill(initial_forcing, 0.0);
    tmc_solver_add_forcing_terms(solver, state, Y, initial_forcing);
    tmc_sparse_fill(jacobian, 0.0);
    tmc_solver_subtract_jacobian_terms(solver, state, Y, jacobian);

    bool accepted = false;
    double last_alpha = 0.0;
    while (!accepted)
    {
      double alpha = 1.0 / (H * p->gamma[0]);
      const double cumulative_alpha = alpha;
      alpha -= last_alpha;
      last_alpha = cumulative_alpha;
      tmc_sparse_add_to_diagonal(jacobian, alpha);
      tmc_lu_decompose(lu, jacobian, lower, upper);

      for (int stage = 0; stage < p->stages; ++stage)
      {
        const int stage_combinations = stage * (stage - 1) / 2;
        if (stage == 0)
          tmc_dense_copy_into(K[0], initial_forcing);
        else if (p->new_function_evaluation[stage])
        {
          tmc_dense_copy_into(Ynew, Y);
          for (int j = 0; j < stage; ++j)
            tmc_dense_axpy(Ynew, p->a[stage_combinations + j], K[j]);
          tmc_dense_fill(K[stage], 0.0);
          tmc_solver_add_forcing_terms(solver, state, Ynew, K[stage]);
        }
        if (stage + 1 < p->stages && !p->new_function_evaluation[stage + 1])
          tmc_dense_copy_into(K[stage + 1], K[stage]);
        for (int j = 0; j < stage; ++j)
          tmc_dense_axpy(K[stage], p->c[stage_combinations + j] / H, K[j]);
        tmc_lu_solve(lu, lower, upper, K[stage]);
      }

      tmc_dense_copy_into(Ynew, Y);
      for (int stage = 0; stage < p->stages; ++stage)
        tmc_dense_axpy(Ynew, p->m[stage], K[stage]);
      tmc_dense_fill(Yerror, 0.0);
      for (int stage = 0; stage < p->stages; ++stage)
        tmc_dense_axpy(Yerror, p->e[stage], K[stage]);

      const double error = tmc_solver_normalized_error(solver, state, Y, Ynew, Yerror);
      double Hnew = H * tmc_default_step_controller(error, p, NULL);
      counts.steps += 1;
      if (isnan(error))
      {
        counts.state = TMC_NAN_DETECTED;
        break;
      }
      if (isinf(error))
      {
        counts.state = TMC_INF_DETECTED;
        break;
      }
      if ((error < 1) || (H < h_min))
      {
        counts.accepted += 1;
        present_time = present_time + H;
        tmc_dense_copy_into(Y, Ynew);
        Hnew = fmax(h_min, fmin(Hnew, h_max));
        if (reject_last_h)
          Hnew = fmin(Hnew, H);
        reject_last_h = false;
        reject_more_h = false;
        H = Hnew;
        accepted = true;
      }
      else
      {
        if (reject_more_h)
          Hnew = H * p->rejection_factor_decrease;
        reject_more_h = reject_last_h;
        reject_last_h = true;
        H = Hnew;
        if (counts.accepted >= 1)
          counts.rejected += 1;
      }
    }
  }
  if (counts.state == TMC_RUNNING)
    counts.state = TMC_CONVERGED;

  for (int i = 0; i < p->stages; ++i)
    tmc_dense_free(K[i]);
  tmc_dense_free(Y);
  tmc_dense_free(Ynew);
  tmc_dense_free(Yerror);
  tmc_dense_free(initial_forcing);
  tmc_sparse_free(jacobian);
  tmc_sparse_free(lower);
  tmc_sparse_free(upper);
  tmc_lu_free(lu);
  return counts;
}

static tmc_state *rejected_step_state(const tmc_solver *solver)
{
  tmc_state *state = tmc_solver_create_state(solver, 1);
  for (int v = 0; v < 5; ++v)
    tmc_state_set(state, 0, v, 1.0e12);
  tmc_state_conditions(state)[0].temperature = 298.15;
  tmc_solver_update_state_parameters(solver, state);
  return state;
}

TMC_TEST(rejected_step_alpha_matches_in_place_solver)
{
  tmc_rosenbrock_parameters options = tmc_rosenbrock_three_stage();
  options.h_start = 1.0e3;
  /* MICM's two builders are standard ordering (L = 1) */
  tmc_solver *solver = get_solver(options, 1, TMC_LU_MOZART);
  TMC_ASSERT_TRUE(solver != NULL);

  /* the standard (non-in-place) solver */
  tmc_state *state = rejected_step_state(solver);
  step_counts standard = reference_standard_solve(solver, state, options.h_start);
  TMC_EXPECT_EQ(standard.state, TMC_CONVERGED);
  TMC_EXPECT_TRUE(standard.steps - standard.accepted >= 2);
  tmc_state_free(state);

  /* the library's (in-place) solver */
  state = rejected_step_state(solver);
  tmc_result result = tmc_solver_solve(solver, options.h_start, state);
  TMC_EXPECT_EQ(result.state, TMC_CONVERGED);
  TMC_EXPECT_TRUE(result.stats.number_of_steps - result.stats.accepted >= 2);
  tmc_state_free(state);

  if (standard.steps != result.stats.number_of_steps || standard.accepted != result.stats.accepted || standard.rejected != result.stats.rejected)
    tmc_test_fail(__FILE__, __LINE__, "steps/accepted/rejected: standard %ld/%ld/%ld, in place %ld/%ld/%ld", standard.steps, standard.accepted, standard.rejected,
                  result.stats.number_of_steps, result.stats.accepted, result.stats.rejected);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("rosenbrock:RosenbrockSolver.StandardAlphaMinusJacobian", standard_alpha_minus_jacobian),
  TMC_CASE("rosenbrock:RosenbrockSolver.VectorAlphaMinusJacobian", vector_alpha_minus_jacobian),
  TMC_CASE("rosenbrock:RosenbrockSolver.CanSetTolerances", can_set_tolerances),
  TMC_CASE("rosenbrock:RosenbrockSolver.StandardNormalizedError", standard_normalized_error),
  TMC_CASE("rosenbrock:RosenbrockSolver.VectorNormalizedError", vector_normalized_error),
  TMC_CASE("rosenbrock:RosenbrockSolver.StandardNormalizedErrorWithConstraints", standard_normalized_error_with_constraints),
  TMC_CASE("rosenbrock:RosenbrockSolver.VectorNormalizedErrorWithConstraints", vector_normalized_error_with_constraints),
  TMC_CASE("rosenbrock:RosenbrockSolver.RejectedStepAlphaMatchesInPlaceSolver", rejected_step_alpha_matches_in_place_solver),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_rosenbrock(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
