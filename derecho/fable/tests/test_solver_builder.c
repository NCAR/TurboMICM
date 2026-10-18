/* MICM test/unit/solver/test_solver_builder.cpp
 *
 * MICM's CpuSolverBuilder with the default policies is a builder of width 1
 * (Matrix and the standard sparse ordering); with VectorMatrix<4> and
 * SparseMatrixVectorOrdering<4> it is a builder of width 4. Solve(dt, state,
 * options) is tmc_solver_set_*_parameters followed by tmc_solver_solve. */
#include "harness.h"

typedef struct
{
  tmc_species a, b, c;
  tmc_phase gas_phase;
  tmc_system system;
  tmc_reaction reactions[2];
} mechanism;

/* A -> B and B -> C, as in MICM's test */
static void mechanism_create(mechanism *m)
{
  m->a = tmc_species_create("A");
  m->b = tmc_species_create("B");
  m->c = tmc_species_create("C");
  m->gas_phase = tmc_phase_create("gas");
  tmc_phase_add_species(&m->gas_phase, &m->a);
  tmc_phase_add_species(&m->gas_phase, &m->b);
  tmc_phase_add_species(&m->gas_phase, &m->c);

  tmc_arrhenius_parameters p1 = tmc_arrhenius_defaults();
  p1.A = 2.15e-11;
  p1.B = 0;
  p1.C = 110;
  m->reactions[0] = tmc_reaction_create();
  tmc_reaction_add_reactant(&m->reactions[0], &m->a);
  tmc_reaction_add_product(&m->reactions[0], &m->b, 1);
  tmc_reaction_set_rate_constant(&m->reactions[0], tmc_rate_arrhenius(p1));
  tmc_reaction_set_phase(&m->reactions[0], &m->gas_phase);

  tmc_arrhenius_parameters p2 = tmc_arrhenius_defaults();
  p2.A = 3.3e-11;
  p2.B = 0;
  p2.C = 55;
  m->reactions[1] = tmc_reaction_create();
  tmc_reaction_add_reactant(&m->reactions[1], &m->b);
  tmc_reaction_add_product(&m->reactions[1], &m->c, 1);
  tmc_reaction_set_rate_constant(&m->reactions[1], tmc_rate_arrhenius(p2));
  tmc_reaction_set_phase(&m->reactions[1], &m->gas_phase);

  m->system = tmc_system_create(&m->gas_phase);
}

static void mechanism_free(mechanism *m)
{
  tmc_reaction_free(&m->reactions[0]);
  tmc_reaction_free(&m->reactions[1]);
  tmc_system_free(&m->system);
  tmc_phase_free(&m->gas_phase);
}

static tmc_solver *build(const mechanism *m, tmc_builder *b, int L)
{
  tmc_builder_set_system(b, &m->system);
  tmc_builder_set_reactions(b, m->reactions, 2);
  tmc_builder_set_width(b, L);
  tmc_error err;
  tmc_solver *s = tmc_builder_build(b, &err);
  if (!s)
    tmc_test_fail(__FILE__, __LINE__, "build failed: %s", err.message);
  tmc_builder_free(b);
  return s;
}

TMC_TEST(throws_missing_system)
{
  tmc_builder *b = tmc_builder_create_backward_euler(tmc_backward_euler_defaults());
  tmc_error err;
  tmc_solver *s = tmc_builder_build(b, &err);
  TMC_EXPECT_TRUE(s == NULL);
  TMC_EXPECT_EQ(err.category, TMC_ERROR_SOLVER);
  TMC_EXPECT_EQ(err.code, TMC_SOLVER_MISSING_CHEMICAL_SYSTEM);
  tmc_solver_free(s);
  tmc_builder_free(b);
  return TMC_TEST_PASS;
}

TMC_TEST(throws_unused_species)
{
  mechanism m;
  mechanism_create(&m);
  tmc_species d = tmc_species_create("D");
  tmc_phase phase_with_unused = tmc_phase_copy(&m.gas_phase);
  tmc_phase_add_species(&phase_with_unused, &d);
  tmc_system system_with_unused = tmc_system_create(&phase_with_unused);

  tmc_builder *b = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(b, &system_with_unused);
  tmc_builder_set_reactions(b, m.reactions, 2);
  tmc_builder_set_ignore_unused_species(b, false);
  tmc_error err;
  tmc_solver *s = tmc_builder_build(b, &err);
  TMC_EXPECT_TRUE(s == NULL);
  TMC_EXPECT_EQ(err.category, TMC_ERROR_SOLVER);
  TMC_EXPECT_STREQ(tmc_error_category_name(err.category), "Solver");
  TMC_EXPECT_EQ(err.code, TMC_SOLVER_UNUSED_SPECIES);
  TMC_EXPECT_TRUE(strstr(err.message, "'D'") != NULL);
  tmc_solver_free(s);

  /* and the builder accepts the system when told to ignore unused species */
  tmc_builder_set_ignore_unused_species(b, true);
  s = tmc_builder_build(b, &err);
  TMC_EXPECT_TRUE(s != NULL);
  tmc_solver_free(s);

  tmc_builder_free(b);
  tmc_system_free(&system_with_unused);
  tmc_phase_free(&phase_with_unused);
  mechanism_free(&m);
  return TMC_TEST_PASS;
}

/* MICM's MaximumNumberOfGridCells is the number of cells in a group of a
 * vectorized solver (4 for VectorMatrix<4>) and unlimited for the standard
 * solver. In C a group is a block of tmc_solver_width cells; width 1 is
 * MICM's standard ordering, where a state holds any number of cells. */
static void check_built(tmc_solver *standard, tmc_solver *vector, const mechanism *m, tmc_method method)
{
  TMC_EXPECT_EQ(tmc_solver_method(standard), method);
  TMC_EXPECT_EQ(tmc_solver_method(vector), method);
  TMC_EXPECT_STREQ(tmc_solver_system(vector)->gas_phase.name, m->system.gas_phase.name);
  TMC_EXPECT_EQ(tmc_solver_system(vector)->gas_phase.n_species, m->system.gas_phase.n_species);

  TMC_EXPECT_EQ(tmc_solver_width(standard), 1);
  TMC_EXPECT_EQ(tmc_solver_width(vector), 4);
  tmc_state *big = tmc_solver_create_state(standard, 1000);
  TMC_EXPECT_EQ(tmc_state_num_cells(big), 1000);
  TMC_EXPECT_EQ(tmc_state_width(big), 1);
  tmc_state_free(big);
  tmc_state *four = tmc_solver_create_state(vector, 4), *five = tmc_solver_create_state(vector, 5);
  TMC_EXPECT_EQ(tmc_state_width(four), 4);
  TMC_EXPECT_EQ(tmc_state_num_blocks(four), 1);
  TMC_EXPECT_EQ(tmc_state_num_blocks(five), 2);
  tmc_state_free(four);
  tmc_state_free(five);
}

TMC_TEST(can_build_backward_euler)
{
  mechanism m;
  mechanism_create(&m);
  tmc_solver *backward_euler = build(&m, tmc_builder_create_backward_euler(tmc_backward_euler_defaults()), 1);
  tmc_solver *backward_euler_vector = build(&m, tmc_builder_create_backward_euler(tmc_backward_euler_defaults()), 4);
  if (backward_euler && backward_euler_vector)
    check_built(backward_euler, backward_euler_vector, &m, TMC_METHOD_BACKWARD_EULER);
  tmc_solver_free(backward_euler);
  tmc_solver_free(backward_euler_vector);
  mechanism_free(&m);
  return TMC_TEST_PASS;
}

TMC_TEST(can_build_rosenbrock)
{
  mechanism m;
  mechanism_create(&m);
  tmc_solver *rosenbrock = build(&m, tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage()), 1);
  tmc_solver *rosenbrock_vector = build(&m, tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage()), 4);
  if (rosenbrock && rosenbrock_vector)
    check_built(rosenbrock, rosenbrock_vector, &m, TMC_METHOD_ROSENBROCK);
  tmc_solver_free(rosenbrock);
  tmc_solver_free(rosenbrock_vector);
  mechanism_free(&m);
  return TMC_TEST_PASS;
}

TMC_TEST(can_build_backward_euler_overloaded_solve)
{
  mechanism m;
  mechanism_create(&m);
  tmc_solver *solver = build(&m, tmc_builder_create_backward_euler(tmc_backward_euler_defaults()), 1);
  TMC_ASSERT_TRUE(solver != NULL);
  tmc_state *state = tmc_solver_create_state(solver, 1);
  tmc_backward_euler_parameters options = tmc_backward_euler_defaults();
  tmc_solver_set_backward_euler_parameters(solver, &options);
  tmc_result solve = tmc_solver_solve(solver, 5, state);

  TMC_EXPECT_EQ(solve.stats.final_time, 5);
  TMC_EXPECT_EQ(solve.stats.function_calls, 2);
  TMC_EXPECT_EQ(solve.stats.jacobian_updates, 2);
  TMC_EXPECT_EQ(solve.stats.number_of_steps, 2);
  TMC_EXPECT_EQ(solve.stats.solves, 2);

  options.small = 1.0;
  options.max_number_of_steps = 1;
  tmc_solver_set_backward_euler_parameters(solver, &options);
  solve = tmc_solver_solve(solver, 5, state);

  TMC_EXPECT_EQ(solve.stats.final_time, 0.03125);
  TMC_EXPECT_EQ(solve.stats.function_calls, 6);
  TMC_EXPECT_EQ(solve.stats.jacobian_updates, 6);
  TMC_EXPECT_EQ(solve.stats.number_of_steps, 6);
  TMC_EXPECT_EQ(solve.stats.solves, 6);

  tmc_state_free(state);
  tmc_solver_free(solver);
  mechanism_free(&m);
  return TMC_TEST_PASS;
}

/* MICM sets state.variables_[0] = { 1.0, 0.0, 0.0 }: by position, in the
 * solver's variable order. The rate constants are never computed (0), so
 * the solution, and the step statistics checked, do not depend on which
 * species holds the 1. */
static void set_first_variable(tmc_state *state)
{
  for (int i = 0; i < tmc_state_num_variables(state); ++i)
    tmc_state_set(state, 0, i, i == 0 ? 1.0 : 0.0);
}

TMC_TEST(can_build_rosenbrock_overloaded_solve)
{
  mechanism m;
  mechanism_create(&m);
  tmc_rosenbrock_parameters options = tmc_rosenbrock_three_stage();
  tmc_solver *solver = build(&m, tmc_builder_create_rosenbrock(options), 1);
  TMC_ASSERT_TRUE(solver != NULL);
  tmc_state *state = tmc_solver_create_state(solver, 1);
  set_first_variable(state);

  tmc_result solve = tmc_solver_solve(solver, 5, state);
  TMC_EXPECT_EQ(solve.stats.final_time, 5);
  TMC_EXPECT_EQ(solve.stats.function_calls, 18);
  TMC_EXPECT_EQ(solve.stats.jacobian_updates, 9);
  TMC_EXPECT_EQ(solve.stats.number_of_steps, 9);
  TMC_EXPECT_EQ(solve.stats.solves, 27);

  options.h_min = 15.0;
  options.max_number_of_steps = 6;
  set_first_variable(state);
  tmc_solver_set_rosenbrock_parameters(solver, &options);
  solve = tmc_solver_solve(solver, 5, state);
  TMC_EXPECT_EQ(solve.stats.final_time, 5);
  TMC_EXPECT_EQ(solve.stats.function_calls, 2);
  TMC_EXPECT_EQ(solve.stats.jacobian_updates, 1);
  TMC_EXPECT_EQ(solve.stats.number_of_steps, 1);
  TMC_EXPECT_EQ(solve.stats.solves, 3);

  /* the solver keeps the parameters it was last given */
  TMC_EXPECT_EQ(tmc_solver_rosenbrock_parameters(solver)->h_min, 15.0);
  TMC_EXPECT_EQ(tmc_solver_rosenbrock_parameters(solver)->max_number_of_steps, 6);

  tmc_state_free(state);
  tmc_solver_free(solver);
  mechanism_free(&m);
  return TMC_TEST_PASS;
}

TMC_TEST(finalize_receives_dense_matrix_policy)
{
  TMC_NA("C++ template dispatch: checks by typeid that FinalizeProcessSetup receives the DenseMatrixPolicy template argument; "
         "C has no matrix-policy types, and finalize_process_setup has one signature");
}

TMC_TEST(finalize_supports_sparse_only_signature)
{
  TMC_NA("C++ overload selection: checks that the builder calls the FinalizeProcessSetup overload without a DenseMatrixPolicy "
         "parameter; the C callback has one signature");
}

static const tmc_test_case cases[] = {
  TMC_CASE("solver_builder:SolverBuilder.ThrowsMissingSystem", throws_missing_system),
  TMC_CASE("solver_builder:SolverBuilder.ThrowsUnusedSpecies", throws_unused_species),
  TMC_CASE("solver_builder:SolverBuilder.CanBuildBackwardEuler", can_build_backward_euler),
  TMC_CASE("solver_builder:SolverBuilder.CanBuildRosenbrock", can_build_rosenbrock),
  TMC_CASE("solver_builder:SolverBuilder.CanBuildBackwardEulerOverloadedSolverMethod", can_build_backward_euler_overloaded_solve),
  TMC_CASE("solver_builder:SolverBuilder.CanBuildRosenbrockOverloadedSolveMethod", can_build_rosenbrock_overloaded_solve),
  TMC_CASE("solver_builder:SolverBuilder.FinalizeReceivesDenseMatrixPolicy", finalize_receives_dense_matrix_policy),
  TMC_CASE("solver_builder:SolverBuilder.FinalizeSupportsSparseOnlySignature", finalize_supports_sparse_only_signature),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_solver_builder(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
