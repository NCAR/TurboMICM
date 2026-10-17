/* Checks of this port's own, beyond MICM's cases: the generated TS1 kernels
 * against the generic kernels, and the TS1 gate on the reference input. They
 * run with `test_cpu --extra` (the Makefile's test-cpu target does). */
#include "harness.h"
#include "ts1.h"

#include <stdio.h>

typedef struct csv
{
  int n_columns, n_rows;
  char (*names)[64];
  double *values;
} csv;

static csv read_csv(const char *path)
{
  csv t = { 0, 0, NULL, NULL };
  FILE *f = fopen(path, "r");
  if (!f)
    return t;
  char line[65536];
  if (!fgets(line, sizeof line, f))
  {
    fclose(f);
    return t;
  }
  int cap = 512;
  t.names = malloc((size_t)cap * sizeof *t.names);
  for (char *tok = strtok(line, ",\r\n"); tok; tok = strtok(NULL, ",\r\n"))
    snprintf(t.names[t.n_columns++], sizeof t.names[0], "%s", tok);
  int rows_cap = 128;
  t.values = malloc((size_t)rows_cap * t.n_columns * sizeof(double));
  while (fgets(line, sizeof line, f))
  {
    if (t.n_rows == rows_cap)
    {
      rows_cap *= 2;
      t.values = realloc(t.values, (size_t)rows_cap * t.n_columns * sizeof(double));
    }
    int k = 0;
    for (char *tok = strtok(line, ",\r\n"); tok && k < t.n_columns; tok = strtok(NULL, ",\r\n"))
      t.values[(size_t)t.n_rows * t.n_columns + k++] = strtod(tok, NULL);
    if (k == t.n_columns)
      ++t.n_rows;
  }
  fclose(f);
  return t;
}

static void free_csv(csv *t)
{
  free(t->names);
  free(t->values);
}

static tmc_solver *build_ts1(bool kernels, int width, tmc_error *err)
{
  tmc_system system = ts1_system();
  tmc_reaction *reactions;
  int n = ts1_reactions(&reactions);
  tmc_builder *b = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, reactions, n);
  tmc_builder_set_ignore_unused_species(b, true);
  tmc_builder_set_width(b, width);
  if (kernels)
    tmc_builder_set_kernels(b, ts1_kernels());
  tmc_solver *s = tmc_builder_build(b, err);
  tmc_builder_free(b);
  tmc_reactions_free(reactions, n);
  tmc_system_free(&system);
  return s;
}

/* Cells first..first+n-1 of the table into a state, as the benchmark does. */
static void load_cells(const tmc_solver *s, tmc_state *st, const csv *t, int first, int n)
{
  (void)s;
  ts1_initial_state(st);
  tmc_conditions *cond = tmc_state_conditions(st);
  for (int c = 0; c < n; ++c)
  {
    const double *row = t->values + (size_t)(first + c) * t->n_columns;
    for (int k = 0; k < t->n_columns; ++k)
    {
      if (strcmp(t->names[k], "temperature") == 0)
        cond[c].temperature = row[k];
      else if (strcmp(t->names[k], "pressure") == 0)
        cond[c].pressure = row[k];
      else
      {
        int v = tmc_state_variable_index(st, t->names[k]);
        if (v >= 0)
          tmc_state_set(st, c, v, row[k]);
        else
        {
          int p = tmc_state_custom_parameter_index(st, t->names[k]);
          if (p >= 0)
            tmc_dense_set(tmc_state_custom_rate_parameters(st), c, p, row[k]);
        }
      }
    }
    tmc_conditions_calculate_ideal_air_density(&cond[c]);
  }
}

/* The generated kernels and the generic kernels on 8 cells of the reference input, 3 steps. */
TMC_TEST(generated_matches_generic)
{
  csv t = read_csv(tmc_test_reference_path("ts1_cells_100_seed1.csv"));
  TMC_ASSERT_TRUE(t.n_rows >= 8);
  tmc_error err;
  tmc_solver *gen = build_ts1(true, TMC_WIDTH, &err);
  TMC_ASSERT_TRUE(gen != NULL);
  tmc_solver *ref = build_ts1(false, TMC_WIDTH, &err);
  TMC_ASSERT_TRUE(ref != NULL);
  TMC_EXPECT_TRUE(tmc_solver_uses_kernels(gen));
  TMC_EXPECT_FALSE(tmc_solver_uses_kernels(ref));
  tmc_state *a = tmc_solver_create_state(gen, 8), *b = tmc_solver_create_state(ref, 8);
  load_cells(gen, a, &t, 0, 8);
  load_cells(ref, b, &t, 0, 8);
  tmc_solver_update_state_parameters(gen, a);
  tmc_solver_update_state_parameters(ref, b);
  double worst = 0.0;
  for (int step = 0; step < 3; ++step)
  {
    tmc_result ra = tmc_solver_solve(gen, 30.0, a), rb = tmc_solver_solve(ref, 30.0, b);
    TMC_EXPECT_EQ(ra.state, TMC_CONVERGED);
    TMC_EXPECT_EQ(rb.state, TMC_CONVERGED);
    TMC_EXPECT_EQ(ra.stats.number_of_steps, rb.stats.number_of_steps);
  }
  for (int c = 0; c < 8; ++c)
    for (int v = 0; v < tmc_state_num_variables(a); ++v)
    {
      double x = tmc_state_get(a, c, v), y = tmc_state_get(b, c, v);
      TMC_EXPECT_TRUE(fabs(x - y) <= 1e-30 + 1e-10 * fabs(y));
      if (y != 0.0 && fabs(x - y) / fabs(y) > worst)
        worst = fabs(x - y) / fabs(y);
    }
  char reason[64];
  snprintf(reason, sizeof reason, "largest relative difference %.2e", worst);
  tmc_test_set_reason(reason);
  tmc_state_free(a);
  tmc_state_free(b);
  tmc_solver_free(gen);
  tmc_solver_free(ref);
  free_csv(&t);
  return TMC_TEST_PASS;
}

/* The experiment's gate on the 100-cell reference input: one warm-up step and
 * 100 steps of 30 s against MICM's final state, in groups of 8 cells. */
TMC_TEST(reference_gate)
{
  csv t = read_csv(tmc_test_reference_path("ts1_cells_100_seed1.csv"));
  csv micm = read_csv(tmc_test_reference_path("micm_ts1_final_100_seed1.csv"));
  TMC_ASSERT_TRUE(t.n_rows == 100 && micm.n_rows == 100);
  tmc_error err;
  tmc_solver *s = build_ts1(true, TMC_WIDTH, &err);
  TMC_ASSERT_TRUE(s != NULL);
  double worst = 0.0;
  long failures = 0;
  for (int first = 0; first < 100; first += 8)
  {
    int n = first + 8 <= 100 ? 8 : 100 - first;
    tmc_state *st = tmc_solver_create_state(s, n);
    load_cells(s, st, &t, first, n);
    tmc_solver_update_state_parameters(s, st);
    for (int step = 0; step < 101; ++step)
    {
      tmc_result r = tmc_solver_solve(s, 30.0, st);
      TMC_EXPECT_EQ(r.state, TMC_CONVERGED);
    }
    for (int c = 0; c < n; ++c)
      for (int k = 0; k < micm.n_columns; ++k)
      {
        int v = tmc_state_variable_index(st, micm.names[k]);
        if (v < 0)
          continue; /* temperature, pressure */
        double x = tmc_state_get(st, c, v), y = micm.values[(size_t)(first + c) * micm.n_columns + k];
        if (!(fabs(x - y) <= 1e-30 + 1e-10 * fabs(y)))
          ++failures;
        if (y != 0.0 && fabs(x - y) / fabs(y) > worst)
          worst = fabs(x - y) / fabs(y);
      }
    tmc_state_free(st);
  }
  TMC_EXPECT_EQ(failures, 0);
  char reason[64];
  snprintf(reason, sizeof reason, "largest relative difference %.2e", worst);
  tmc_test_set_reason(reason);
  tmc_solver_free(s);
  free_csv(&t);
  free_csv(&micm);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("ts1:GeneratedKernels.MatchGeneric", generated_matches_generic),
  TMC_CASE("ts1:GeneratedKernels.ReferenceGate", reference_gate),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_ts1_kernels(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
