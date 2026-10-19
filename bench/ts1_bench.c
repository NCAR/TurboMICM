/* TS1 benchmark: MICM's three-stage Rosenbrock solver on CELLS cells read from
 * a CSV file, one untimed warm-up step, then STEPS timed steps of DT seconds.
 *
 *   ts1_bench_cpu [--threads N] [--group-size G] [--schedule static|dynamic|guided]
 *                 [--chunk K] [--generic] CELLS STEPS DT CELLS_CSV
 *
 * Prints elapsed_ms and the solver statistics summed over the timed steps,
 * and writes the final state to ts1_final.csv. Cells are solved in groups of
 * G (default TMC_DEFAULT_GROUP_SIZE), each group a state with its own step
 * size; an OpenMP loop hands groups to threads with the chosen schedule
 * (default dynamic, one group at a time). --generic runs the library's
 * generic kernels instead of the generated TS1 kernels. */
#include "csv_table.h"
#include "ts1.h"
#include "tmc.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TMC_DEFAULT_GROUP_SIZE
#define TMC_DEFAULT_GROUP_SIZE 32
#endif


typedef struct
{
  tmc_state *state;
  int first_cell, n_cells;
} group;

int main(int argc, char **argv)
{
  tmc_program = "ts1_bench_cpu";
  int threads = 0, group_size = TMC_DEFAULT_GROUP_SIZE, chunk = 1;
  bool generic = false;
  omp_sched_t schedule = omp_sched_dynamic;
  int i = 1;
  for (; i < argc && argv[i][0] == '-'; ++i)
  {
    if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc)
      threads = atoi(argv[++i]);
    else if (strcmp(argv[i], "--group-size") == 0 && i + 1 < argc)
      group_size = atoi(argv[++i]);
    else if (strcmp(argv[i], "--chunk") == 0 && i + 1 < argc)
      chunk = atoi(argv[++i]);
    else if (strcmp(argv[i], "--generic") == 0)
      generic = true;
    else if (strcmp(argv[i], "--schedule") == 0 && i + 1 < argc)
    {
      ++i;
      if (strcmp(argv[i], "static") == 0)
        schedule = omp_sched_static;
      else if (strcmp(argv[i], "dynamic") == 0)
        schedule = omp_sched_dynamic;
      else if (strcmp(argv[i], "guided") == 0)
        schedule = omp_sched_guided;
      else
        fail("--schedule takes static, dynamic, or guided");
    }
    else
      fail("unknown option");
  }
  if (argc - i != 4)
  {
    fprintf(stderr, "usage: ts1_bench_cpu [--threads N] [--group-size G] [--schedule S] [--chunk K] [--generic] CELLS STEPS DT CELLS_CSV\n");
    return 1;
  }
  const int n_cells = atoi(argv[i]);
  const int n_steps = atoi(argv[i + 1]);
  const double dt = atof(argv[i + 2]);
  const char *csv_path = argv[i + 3];
  if (n_cells <= 0 || n_steps < 0 || dt <= 0 || group_size <= 0 || chunk <= 0)
    fail("CELLS, STEPS, DT, group size, and chunk must be positive");
  if (threads > 0)
    omp_set_num_threads(threads);
  omp_set_schedule(schedule, chunk);

  /* The solver: TS1 from the generated configuration tables */
  tmc_system system = ts1_system();
  tmc_reaction *reactions;
  int n_reactions = ts1_reactions(&reactions);
  tmc_builder *builder = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(builder, &system);
  tmc_builder_set_reactions(builder, reactions, n_reactions);
  tmc_builder_set_ignore_unused_species(builder, true);
  if (!generic)
    tmc_builder_set_kernels(builder, ts1_kernels());
  tmc_error error;
  tmc_solver *solver = tmc_builder_build(builder, &error);
  if (!solver)
  {
    fprintf(stderr, "ts1_bench_cpu: %s\n", error.message);
    return 1;
  }
  tmc_builder_free(builder);
  tmc_reactions_free(reactions, n_reactions);

  /* The input: every column must name temperature, pressure, a species, or a rate parameter */
  table t = read_csv(csv_path);
  if (t.n_rows < n_cells)
    fail("CELLS_CSV has fewer rows than CELLS");
  int *column_kind = malloc((size_t)t.n_columns * sizeof(int));  /* 0 T, 1 P, 2 species, 3 parameter */
  int *column_index = malloc((size_t)t.n_columns * sizeof(int));
  {
    tmc_state *probe = tmc_solver_create_state(solver, 1);
    for (int k = 0; k < t.n_columns; ++k)
    {
      if (strcmp(t.names[k], "temperature") == 0)
        column_kind[k] = 0;
      else if (strcmp(t.names[k], "pressure") == 0)
        column_kind[k] = 1;
      else if ((column_index[k] = tmc_state_variable_index(probe, t.names[k])) >= 0)
        column_kind[k] = 2;
      else if ((column_index[k] = tmc_state_custom_parameter_index(probe, t.names[k])) >= 0)
        column_kind[k] = 3;
      else
      {
        fprintf(stderr, "ts1_bench_cpu: CELLS_CSV column '%s' is not a TS1 species or rate parameter\n", t.names[k]);
        return 1;
      }
    }
    tmc_state_free(probe);
  }

  const int n_groups = (n_cells + group_size - 1) / group_size;
  group *groups = malloc((size_t)n_groups * sizeof *groups);
  for (int g = 0; g < n_groups; ++g)
  {
    groups[g].first_cell = g * group_size;
    groups[g].n_cells = (g == n_groups - 1) ? n_cells - g * group_size : group_size;
    groups[g].state = NULL;
  }
  int n_threads = 1;
  volatile int failed = 0;
#pragma omp parallel
  {
#pragma omp single
    n_threads = omp_get_num_threads();
  }
  tmc_stats *thread_stats = calloc((size_t)n_threads, sizeof *thread_stats);

  /* Initial state, rate constants, and the untimed warm-up step, per group */
#pragma omp parallel for schedule(runtime)
  for (int g = 0; g < n_groups; ++g)
  {
    if (failed)
      continue;
    tmc_state *state = tmc_solver_create_state(solver, groups[g].n_cells);
    ts1_initial_state(state);
    tmc_conditions *cond = tmc_state_conditions(state);
    for (int c = 0; c < groups[g].n_cells; ++c)
    {
      const double *row = t.values + (size_t)(groups[g].first_cell + c) * (size_t)t.n_columns;
      for (int k = 0; k < t.n_columns; ++k)
      {
        switch (column_kind[k])
        {
          case 0: cond[c].temperature = row[k]; break;
          case 1: cond[c].pressure = row[k]; break;
          case 2: tmc_state_set(state, c, column_index[k], row[k]); break;
          default: tmc_dense_set(tmc_state_custom_rate_parameters(state), c, column_index[k], row[k]); break;
        }
      }
      tmc_conditions_calculate_ideal_air_density(&cond[c]);
    }
    tmc_solver_update_state_parameters(solver, state);
    tmc_result r = tmc_solver_solve(solver, dt, state);
    if (r.state != TMC_CONVERGED)
    {
      fprintf(stderr, "ts1_bench_cpu: warm-up solve of cells %d-%d: %s\n", groups[g].first_cell, groups[g].first_cell + groups[g].n_cells - 1, tmc_solver_state_name(r.state));
      failed = 1;
    }
    groups[g].state = state;
  }
  if (failed)
    return 1;

  const double t0 = omp_get_wtime();
#pragma omp parallel for schedule(runtime)
  for (int g = 0; g < n_groups; ++g)
  {
    if (failed)
      continue;
    tmc_stats *stats = &thread_stats[omp_get_thread_num()];
    tmc_state *state = groups[g].state;
    for (int step = 0; step < n_steps; ++step)
    {
      tmc_result r = tmc_solver_solve(solver, dt, state);
      stats->number_of_steps += r.stats.number_of_steps;
      stats->accepted += r.stats.accepted;
      stats->rejected += r.stats.rejected;
      stats->function_calls += r.stats.function_calls;
      stats->jacobian_updates += r.stats.jacobian_updates;
      stats->decompositions += r.stats.decompositions;
      stats->solves += r.stats.solves;
      if (r.state != TMC_CONVERGED)
      {
        fprintf(stderr, "ts1_bench_cpu: timed step %d of cells %d-%d: %s\n", step, groups[g].first_cell, groups[g].first_cell + groups[g].n_cells - 1, tmc_solver_state_name(r.state));
        failed = 1;
        break;
      }
    }
  }
  const double t1 = omp_get_wtime();
  if (failed)
    return 1;

  tmc_stats total;
  memset(&total, 0, sizeof total);
  for (int k = 0; k < n_threads; ++k)
  {
    total.number_of_steps += thread_stats[k].number_of_steps;
    total.accepted += thread_stats[k].accepted;
    total.rejected += thread_stats[k].rejected;
    total.function_calls += thread_stats[k].function_calls;
    total.jacobian_updates += thread_stats[k].jacobian_updates;
    total.decompositions += thread_stats[k].decompositions;
    total.solves += thread_stats[k].solves;
  }

  FILE *out = fopen("ts1_final.csv", "w");
  if (!out)
    fail("cannot write ts1_final.csv");
  const int n_vars = tmc_solver_num_variables(solver);
  fprintf(out, "temperature,pressure");
  for (int v = 0; v < n_vars; ++v)
    fprintf(out, ",%s", tmc_solver_variable_name(solver, v));
  fprintf(out, "\n");
  for (int g = 0; g < n_groups; ++g)
  {
    const tmc_conditions *cond = tmc_state_conditions(groups[g].state);
    for (int c = 0; c < groups[g].n_cells; ++c)
    {
      fprintf(out, "%.17g,%.17g", cond[c].temperature, cond[c].pressure);
      for (int v = 0; v < n_vars; ++v)
        fprintf(out, ",%.17g", tmc_state_get(groups[g].state, c, v));
      fprintf(out, "\n");
    }
  }
  fclose(out);

  printf("elapsed_ms=%.3f\n", (t1 - t0) * 1000.0);
  printf("internal_steps=%ld accepted=%ld rejected=%ld function_calls=%ld jacobian_updates=%ld decompositions=%ld solves=%ld\n",
         total.number_of_steps, total.accepted, total.rejected, total.function_calls, total.jacobian_updates, total.decompositions, total.solves);

  for (int g = 0; g < n_groups; ++g)
    tmc_state_free(groups[g].state);
  free(groups);
  free(thread_stats);
  free(column_kind);
  free(column_index);
  free(t.names);
  free(t.values);
  tmc_solver_free(solver);
  tmc_system_free(&system);
  return 0;
}
