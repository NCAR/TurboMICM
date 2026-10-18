/* MICM test/unit/openmp/test_openmp.cpp: one solver shared by eight OpenMP
 * threads, each with its own state, must give every thread the same result
 * (run_solver.hpp's ten 200 s steps of the Robertson mechanism with
 * user-defined rate constants). */
#include "harness.h"

#include <omp.h>

#define N_THREADS 8
#define N_VARS 3

typedef struct abc_mechanism
{
  tmc_phase gas;
  tmc_system system;
  tmc_reaction reactions[3];
} abc_mechanism;

static abc_mechanism abc_create(void)
{
  abc_mechanism m;
  tmc_species a = tmc_species_create("A"), b = tmc_species_create("B"), c = tmc_species_create("C");
  m.gas = tmc_phase_create("gas");
  tmc_phase_add_species(&m.gas, &a);
  tmc_phase_add_species(&m.gas, &b);
  tmc_phase_add_species(&m.gas, &c);
  m.reactions[0] = tmc_reaction_create();
  tmc_reaction_add_reactant(&m.reactions[0], &a);
  tmc_reaction_add_product(&m.reactions[0], &b, 1);
  tmc_reaction_set_rate_constant(&m.reactions[0], tmc_rate_user_defined("r1", 1.0));
  tmc_reaction_set_phase(&m.reactions[0], &m.gas);
  m.reactions[1] = tmc_reaction_create();
  tmc_reaction_add_reactant(&m.reactions[1], &b);
  tmc_reaction_add_reactant(&m.reactions[1], &b);
  tmc_reaction_add_product(&m.reactions[1], &b, 1);
  tmc_reaction_add_product(&m.reactions[1], &c, 1);
  tmc_reaction_set_rate_constant(&m.reactions[1], tmc_rate_user_defined("r2", 1.0));
  tmc_reaction_set_phase(&m.reactions[1], &m.gas);
  m.reactions[2] = tmc_reaction_create();
  tmc_reaction_add_reactant(&m.reactions[2], &b);
  tmc_reaction_add_reactant(&m.reactions[2], &c);
  tmc_reaction_add_product(&m.reactions[2], &a, 1);
  tmc_reaction_add_product(&m.reactions[2], &c, 1);
  tmc_reaction_set_rate_constant(&m.reactions[2], tmc_rate_user_defined("r3", 1.0));
  tmc_reaction_set_phase(&m.reactions[2], &m.gas);
  m.system = tmc_system_create(&m.gas);
  return m;
}

static void abc_free(abc_mechanism *m)
{
  for (int i = 0; i < 3; ++i)
    tmc_reaction_free(&m->reactions[i]);
  tmc_system_free(&m->system);
  tmc_phase_free(&m->gas);
}

/* run_solver.hpp's run_solver_on_thread_with_own_state: the final
 * concentrations of a one-cell state after ten 200 s steps; returns false
 * when a solve makes no progress. */
static bool run_solver_on_thread_with_own_state(const tmc_solver *solver, tmc_state *state, double result[N_VARS])
{
  const double initial[N_VARS] = { 1, 0, 0 }; /* mol m-3 */
  for (int v = 0; v < N_VARS; ++v)
    tmc_state_set(state, 0, v, initial[v]);
  if (tmc_state_set_custom_rate_parameter(state, "r1", 0.04, NULL) != TMC_OK ||
      tmc_state_set_custom_rate_parameter(state, "r2", 3e7, NULL) != TMC_OK ||
      tmc_state_set_custom_rate_parameter(state, "r3", 1e4, NULL) != TMC_OK)
    return false;
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = 272.5;  /* K */
  conditions[0].pressure = 101253.3;  /* Pa */
  conditions[0].air_density = 1e6;    /* mol m-3 */

  const double time_step = 200; /* s */
  bool progressed = true;
  for (int i = 0; i < 10 && progressed; ++i)
  {
    double elapsed_solve_time = 0;
    tmc_solver_update_state_parameters(solver, state);
    while (elapsed_solve_time < time_step)
    {
      tmc_result r = tmc_solver_solve(solver, time_step - elapsed_solve_time, state);
      if (r.stats.final_time <= 0)
      {
        progressed = false;
        break;
      }
      elapsed_solve_time = r.stats.final_time;
    }
  }
  for (int v = 0; v < N_VARS; ++v)
    result[v] = tmc_state_get(state, 0, v);
  return progressed;
}

TMC_TEST(one_solver_many_states)
{
  abc_mechanism m = abc_create();
  tmc_builder *builder = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(builder, &m.system);
  tmc_builder_set_reactions(builder, m.reactions, 3);
  tmc_builder_set_width(builder, 1);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_solver *solver = tmc_builder_build(builder, &err);
  tmc_builder_free(builder);
  abc_free(&m);
  TMC_ASSERT_TRUE(solver != NULL);

  double results[N_THREADS][N_VARS];
  bool progressed[N_THREADS];
  int threads_seen = 0;
#pragma omp parallel num_threads(N_THREADS)
  {
    const int t = omp_get_thread_num();
    tmc_state *state = tmc_solver_create_state(solver, 1);
    progressed[t] = run_solver_on_thread_with_own_state(solver, state, results[t]);
    tmc_state_free(state);
#pragma omp atomic
    ++threads_seen;
  }
  TMC_EXPECT_EQ(threads_seen, N_THREADS);
  for (int t = 0; t < N_THREADS; ++t)
    TMC_EXPECT_TRUE(progressed[t]);

  /* compare each thread to thread 0, bit for bit */
  for (int t = 1; t < N_THREADS; ++t)
    for (int v = 0; v < N_VARS; ++v)
      TMC_EXPECT_EQ(results[0][v], results[t][v]);
  /* and the run did something: A was consumed */
  TMC_EXPECT_TRUE(results[0][0] < 1.0 && results[0][0] > 0.0);

  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("openmp:OpenMP.OneSolverManyStates", one_solver_many_states),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_openmp(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
