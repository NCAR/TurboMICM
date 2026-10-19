/* TS1 scenario benchmark (Tier 3). A scenario file names the cells (a
 * CELLS_CSV), the step length and count, the photolysis factor f(t), and the
 * snapshot interval (see ../scenarios/README.md).
 *
 *   ts1_scenario [--threads N] [--method NAME] [--option KEY=VALUE ...] SCENARIO
 *
 * Every cell starts from TS1's InitState, takes its CELLS_CSV row, and has its
 * air density recomputed. At the start of step k, every photolysis rate
 * parameter (a label beginning with "j") is its CELLS_CSV value times f(k dt),
 * and only the rate constants that read those parameters are recomputed. The
 * timer covers all the steps, including those updates, and nothing else.
 * Prints elapsed_ms, the solver statistics summed over the steps and groups,
 * and this program's own counters; writes ts1_final.csv and, when the scenario
 * asks for snapshots, ts1_snapshots.csv, in the working directory.
 *
 * Methods and options are listed by --help and documented in DESIGN.md. The
 * default, micm3, is MICM's three-stage Rosenbrock in MICM's loop with its
 * parameters and tolerances (relative 1e-6, absolute 1e-3 mol m-3). */
#include "csv_table.h"
#include "ts1.h"
#include "tmc.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TMC_DEFAULT_GROUP_SIZE
#define TMC_DEFAULT_GROUP_SIZE 32
#endif

static const char *USAGE =
    "usage: ts1_scenario [--threads N] [--method NAME] [--option KEY=VALUE ...] SCENARIO\n"
    "methods: micm3 (default) ros2 ros3 ros4 rodas3 rodas4 ros34pw2 extrap\n"
    "options: rtol=R atol=A group_size=G schedule=static|dynamic|guided chunk=K generic=0|1\n"
    "         warm_start=0|1 h_start=SECONDS h_start_factor=F h_max=SECONDS max_steps=M\n"
    "         norm=rms|cellmax|max controller=micm|predictive|h211b clamp=0|1 quantize=0|1\n"
    "         reuse_age=STEPS reuse_change=FRACTION freeze_band=FACTOR\n"
    "         safety=DELTA fac_min=Q fac_max=Q h211_b=B h211_k=K filter_error=0|1\n"
    "         order=K sequence=harmonic|shifted|even (extrap)\n";

typedef struct
{
  double rtol, atol, h_start, h_start_factor, h_max, reuse_change, freeze_band, safety, fac_min, fac_max, h211_b, h211_k;
  int group_size, chunk, generic, warm_start, norm, controller, clamp, quantize, reuse_age, order, sequence, max_steps, filter_error;
  omp_sched_t schedule;
} options;

typedef struct
{
  char cells[4096], photolysis[4096];
  double dt;
  int steps, snapshot_every;
} scenario;

typedef struct
{
  int n;
  double *time, *factor;
} photolysis_table;

typedef struct
{
  tmc_state *state;
  int first_cell, n_cells;
  double *base;      /* [n_cells * n_photolysis] the CELLS_CSV photolysis values */
  double *snapshots; /* [n_snapshots * n_cells * n_vars] */
  long rate_updates;
  char error[256];
} group;

static int parse_bool(const char *v)
{
  if (!strcmp(v, "1") || !strcmp(v, "yes") || !strcmp(v, "true"))
    return 1;
  if (!strcmp(v, "0") || !strcmp(v, "no") || !strcmp(v, "false"))
    return 0;
  return -1;
}

static void option(options *o, const char *key, const char *value)
{
  char message[512];
  int b;
  if (!strcmp(key, "rtol"))
    o->rtol = atof(value);
  else if (!strcmp(key, "atol"))
    o->atol = atof(value);
  else if (!strcmp(key, "group_size"))
    o->group_size = atoi(value);
  else if (!strcmp(key, "chunk"))
    o->chunk = atoi(value);
  else if (!strcmp(key, "schedule"))
  {
    if (!strcmp(value, "static"))
      o->schedule = omp_sched_static;
    else if (!strcmp(value, "dynamic"))
      o->schedule = omp_sched_dynamic;
    else if (!strcmp(value, "guided"))
      o->schedule = omp_sched_guided;
    else
      fail("schedule takes static, dynamic, or guided");
  }
  else if (!strcmp(key, "generic"))
    o->generic = parse_bool(value);
  else if (!strcmp(key, "warm_start"))
    o->warm_start = parse_bool(value);
  else if (!strcmp(key, "h_start"))
    o->h_start = atof(value);
  else if (!strcmp(key, "h_start_factor"))
    o->h_start_factor = atof(value);
  else if (!strcmp(key, "h_max"))
    o->h_max = atof(value);
  else if (!strcmp(key, "max_steps"))
    o->max_steps = atoi(value);
  else if (!strcmp(key, "norm"))
  {
    if (!strcmp(value, "rms"))
      o->norm = TMC_NORM_RMS;
    else if (!strcmp(value, "cellmax"))
      o->norm = TMC_NORM_CELL_MAX;
    else if (!strcmp(value, "max"))
      o->norm = TMC_NORM_MAX;
    else
      fail("norm takes rms, cellmax, or max");
  }
  else if (!strcmp(key, "controller"))
  {
    if (!strcmp(value, "micm"))
      o->controller = TMC_CONTROLLER_MICM;
    else if (!strcmp(value, "predictive"))
      o->controller = TMC_CONTROLLER_PREDICTIVE;
    else if (!strcmp(value, "h211b"))
      o->controller = TMC_CONTROLLER_H211B;
    else
      fail("controller takes micm, predictive, or h211b");
  }
  else if (!strcmp(key, "clamp"))
    o->clamp = parse_bool(value);
  else if (!strcmp(key, "quantize"))
    o->quantize = parse_bool(value);
  else if (!strcmp(key, "reuse_age"))
    o->reuse_age = atoi(value);
  else if (!strcmp(key, "reuse_change"))
    o->reuse_change = atof(value);
  else if (!strcmp(key, "freeze_band"))
    o->freeze_band = atof(value);
  else if (!strcmp(key, "order"))
    o->order = atoi(value);
  else if (!strcmp(key, "safety"))
    o->safety = atof(value);
  else if (!strcmp(key, "fac_min"))
    o->fac_min = atof(value);
  else if (!strcmp(key, "fac_max"))
    o->fac_max = atof(value);
  else if (!strcmp(key, "h211_b"))
    o->h211_b = atof(value);
  else if (!strcmp(key, "h211_k"))
    o->h211_k = atof(value);
  else if (!strcmp(key, "filter_error"))
    o->filter_error = parse_bool(value);
  else if (!strcmp(key, "sequence"))
  {
    if (!strcmp(value, "harmonic"))
      o->sequence = TMC_SEQUENCE_HARMONIC;
    else if (!strcmp(value, "shifted"))
      o->sequence = TMC_SEQUENCE_SHIFTED;
    else if (!strcmp(value, "even"))
      o->sequence = TMC_SEQUENCE_EVEN;
    else
      fail("sequence takes harmonic, shifted, or even");
  }
  else
  {
    snprintf(message, sizeof message, "unknown option '%s'\n%s", key, USAGE);
    fail(message);
  }
  b = 0;
  (void)b;
  if (o->generic < 0 || o->warm_start < -1 || o->clamp < 0 || o->quantize < 0 || o->filter_error < 0)
    fail("a 0/1 option was given another value");
}

/* The method's parameter set, in MICM's loop (micm3) or the ROW or extrapolation loop. */
static bool method_parameters(const char *name, tmc_rosenbrock_parameters *p, int *warm_start_default)
{
  *warm_start_default = 1;
  if (!strcmp(name, "micm3"))
  {
    *p = tmc_rosenbrock_three_stage();
    *warm_start_default = 0;
    return true;
  }
  if (!strcmp(name, "ros2"))
    *p = tmc_rosenbrock_two_stage();
  else if (!strcmp(name, "ros3"))
    *p = tmc_rosenbrock_three_stage();
  else if (!strcmp(name, "ros4"))
    *p = tmc_rosenbrock_four_stage();
  else if (!strcmp(name, "rodas3"))
    *p = tmc_rosenbrock_four_stage_dae();
  else if (!strcmp(name, "rodas4"))
    *p = tmc_rosenbrock_six_stage_dae();
  else if (!strcmp(name, "ros34pw2"))
    *p = tmc_rosenbrock_ros34pw2();
  else if (!strcmp(name, "extrap"))
  {
    *p = tmc_rosenbrock_three_stage();
    p->stages = 1; /* one work vector */
    p->gamma[0] = 1.0;
    p->loop = TMC_LOOP_EXTRAPOLATION;
    return true;
  }
  else
    return false;
  p->loop = TMC_LOOP_ROW;
  return true;
}

static scenario read_scenario(const char *path)
{
  FILE *f = fopen(path, "r");
  if (!f)
    fail("cannot open SCENARIO");
  scenario sc;
  memset(&sc, 0, sizeof sc);
  char directory[4096] = "";
  const char *slash = strrchr(path, '/');
  if (slash)
    snprintf(directory, sizeof directory, "%.*s", (int)(slash - path + 1), path);
  char line[8192];
  bool have_cells = false, have_dt = false, have_steps = false, have_photolysis = false;
  while (fgets(line, sizeof line, f))
  {
    char *hash = strchr(line, '#');
    if (hash)
      *hash = '\0';
    char *equals = strchr(line, '=');
    if (!equals)
      continue;
    *equals = '\0';
    char *key = line, *value = equals + 1;
    while (*key == ' ' || *key == '\t')
      ++key;
    char *end = key + strlen(key);
    while (end > key && (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\r' || end[-1] == '\n'))
      *--end = '\0';
    while (*value == ' ' || *value == '\t')
      ++value;
    end = value + strlen(value);
    while (end > value && (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\r' || end[-1] == '\n'))
      *--end = '\0';
    if (!strcmp(key, "cells"))
    {
      snprintf(sc.cells, sizeof sc.cells, "%s%s", value[0] == '/' ? "" : directory, value);
      have_cells = true;
    }
    else if (!strcmp(key, "photolysis"))
    {
      snprintf(sc.photolysis, sizeof sc.photolysis, "%s%s", value[0] == '/' ? "" : directory, value);
      have_photolysis = true;
    }
    else if (!strcmp(key, "dt"))
    {
      sc.dt = atof(value);
      have_dt = true;
    }
    else if (!strcmp(key, "steps"))
    {
      sc.steps = atoi(value);
      have_steps = true;
    }
    else if (!strcmp(key, "snapshot_every"))
      sc.snapshot_every = atoi(value);
    else
      fail("unknown key in SCENARIO");
  }
  fclose(f);
  if (!have_cells || !have_dt || !have_steps || !have_photolysis)
    fail("SCENARIO needs cells, dt, steps, and photolysis");
  if (sc.dt <= 0 || sc.steps < 1 || sc.snapshot_every < 0)
    fail("SCENARIO: dt must be positive, steps at least 1, snapshot_every non-negative");
  return sc;
}

static photolysis_table read_photolysis(const char *path)
{
  table t = read_csv(path);
  if (t.n_columns != 2 || t.n_rows < 1 || strcmp(t.names[0], "time_s") != 0 || strcmp(t.names[1], "factor") != 0)
    fail("the photolysis file needs a time_s,factor header and rows");
  photolysis_table ph;
  ph.n = t.n_rows;
  ph.time = malloc((size_t)ph.n * sizeof(double));
  ph.factor = malloc((size_t)ph.n * sizeof(double));
  for (int i = 0; i < ph.n; ++i)
  {
    ph.time[i] = t.values[(size_t)i * 2];
    ph.factor[i] = t.values[(size_t)i * 2 + 1];
    if (i > 0 && ph.time[i] <= ph.time[i - 1])
      fail("photolysis times must increase");
  }
  free(t.names);
  free(t.values);
  return ph;
}

/* f(t): linear between the rows, constant before the first and after the last. */
static double photolysis_at(const photolysis_table *ph, double time)
{
  if (time <= ph->time[0])
    return ph->factor[0];
  if (time >= ph->time[ph->n - 1])
    return ph->factor[ph->n - 1];
  int lo = 0, hi = ph->n - 1; /* time[lo] <= time < time[hi] */
  while (hi - lo > 1)
  {
    int mid = (lo + hi) / 2;
    if (ph->time[mid] <= time)
      lo = mid;
    else
      hi = mid;
  }
  const double w = (time - ph->time[lo]) / (ph->time[hi] - ph->time[lo]);
  return (1.0 - w) * ph->factor[lo] + w * ph->factor[hi];
}

static void add_stats(tmc_stats *total, const tmc_stats *part)
{
  total->number_of_steps += part->number_of_steps;
  total->accepted += part->accepted;
  total->rejected += part->rejected;
  total->function_calls += part->function_calls;
  total->jacobian_updates += part->jacobian_updates;
  total->decompositions += part->decompositions;
  total->solves += part->solves;
  total->matrix_reuses += part->matrix_reuses;
}

int main(int argc, char **argv)
{
  tmc_program = "ts1_scenario";
  options o = { -1.0, -1.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, TMC_DEFAULT_GROUP_SIZE, 1, 0, -1, TMC_NORM_RMS, TMC_CONTROLLER_MICM, 0, 0, 0, 4, TMC_SEQUENCE_HARMONIC, 0, 0, omp_sched_dynamic };
  const char *method = "micm3";
  int threads = 0;
  int i = 1;
  for (; i < argc && argv[i][0] == '-'; ++i)
  {
    if (!strcmp(argv[i], "--threads") && i + 1 < argc)
      threads = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--method") && i + 1 < argc)
      method = argv[++i];
    else if (!strcmp(argv[i], "--option") && i + 1 < argc)
    {
      char *spec = argv[++i];
      char *equals = strchr(spec, '=');
      if (!equals)
        fail("--option takes KEY=VALUE");
      *equals = '\0';
      option(&o, spec, equals + 1);
    }
    else if (!strcmp(argv[i], "--help"))
    {
      fputs(USAGE, stdout);
      return 0;
    }
    else
    {
      fprintf(stderr, "%s", USAGE);
      return 2;
    }
  }
  if (argc - i != 1)
  {
    fprintf(stderr, "%s", USAGE);
    return 2;
  }
  if (o.group_size <= 0 || o.chunk <= 0)
    fail("group_size and chunk must be positive");
  if (threads > 0)
    omp_set_num_threads(threads);
  omp_set_schedule(o.schedule, o.chunk);

  /* The method's parameters with the options applied */
  tmc_rosenbrock_parameters p;
  int warm_default;
  if (!method_parameters(method, &p, &warm_default))
  {
    fprintf(stderr, "ts1_scenario: unknown method '%s'\n%s", method, USAGE);
    return 2;
  }
  p.warm_start = o.warm_start < 0 ? warm_default : o.warm_start;
  p.h_start = o.h_start;
  p.h_start_factor = o.h_start_factor;
  p.h_max = o.h_max;
  p.norm = o.norm;
  p.controller_type = o.controller;
  p.clamp_internal = o.clamp;
  p.quantize = o.quantize;
  p.max_matrix_age = o.reuse_age;
  p.matrix_rate_change = o.reuse_change;
  p.freeze_band = o.freeze_band;
  if (o.max_steps > 0)
    p.max_number_of_steps = o.max_steps;
  if (o.safety > 0.0)
    p.safety_factor = o.safety;
  if (o.fac_min > 0.0)
    p.factor_min = o.fac_min;
  if (o.fac_max > 0.0)
    p.factor_max = o.fac_max;
  p.h211_b = o.h211_b;
  p.h211_k = o.h211_k;
  p.filter_error = o.filter_error;
  if (p.loop == TMC_LOOP_EXTRAPOLATION)
  {
    if (o.order < 2 || o.order > TMC_MAX_EXTRAPOLATION)
      fail("order must be between 2 and 8");
    p.extrapolation_order = o.order;
    p.extrapolation_sequence = o.sequence;
    p.estimator_of_local_order = o.order;
  }

  /* The solver: TS1 from the generated configuration tables */
  tmc_system system = ts1_system();
  tmc_reaction *reactions;
  int n_reactions = ts1_reactions(&reactions);
  tmc_builder *builder = tmc_builder_create_rosenbrock(p);
  tmc_builder_set_system(builder, &system);
  tmc_builder_set_reactions(builder, reactions, n_reactions);
  tmc_builder_set_ignore_unused_species(builder, true);
  if (!o.generic)
    tmc_builder_set_kernels(builder, ts1_kernels());
  tmc_error error;
  tmc_solver *solver = tmc_builder_build(builder, &error);
  if (!solver)
  {
    fprintf(stderr, "ts1_scenario: %s\n", error.message);
    return 1;
  }
  tmc_builder_free(builder);
  tmc_reactions_free(reactions, n_reactions);
  const int n_vars = tmc_solver_num_variables(solver);
  const int n_custom = tmc_solver_num_custom_parameters(solver);

  /* The scenario and its inputs */
  const scenario sc = read_scenario(argv[i]);
  const table t = read_csv(sc.cells);
  const photolysis_table ph = read_photolysis(sc.photolysis);
  const int n_cells = t.n_rows;
  if (n_cells < 1)
    fail("the cells file has no cells");
  int *column_kind = malloc((size_t)t.n_columns * sizeof(int)); /* 0 T, 1 P, 2 species, 3 parameter */
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
        fprintf(stderr, "ts1_scenario: CELLS_CSV column '%s' is not a TS1 species or rate parameter\n", t.names[k]);
        return 1;
      }
    }
    tmc_state_free(probe);
  }
  /* the photolysis parameters: labels beginning with "j" */
  int n_photolysis = 0;
  int *photolysis_index = malloc((size_t)(n_custom > 0 ? n_custom : 1) * sizeof(int));
  for (int k = 0; k < n_custom; ++k)
    if (tmc_solver_custom_parameter_label(solver, k)[0] == 'j')
      photolysis_index[n_photolysis++] = k;

  /* the snapshot steps: every snapshot_every steps, and always the last */
  int n_snapshots = 0;
  int *snapshot_steps = malloc((size_t)(sc.steps + 1) * sizeof(int));
  if (sc.snapshot_every > 0)
    for (int step = sc.snapshot_every; step <= sc.steps; step += sc.snapshot_every)
      snapshot_steps[n_snapshots++] = step;
  if (n_snapshots == 0 || snapshot_steps[n_snapshots - 1] != sc.steps)
    snapshot_steps[n_snapshots++] = sc.steps;

  /* The groups: a state each, initialized and with their rate constants, before the timer */
  const int n_groups = (n_cells + o.group_size - 1) / o.group_size;
  group *groups = calloc((size_t)n_groups, sizeof *groups);
  double *atol = NULL;
  if (o.atol > 0.0)
  {
    atol = malloc((size_t)n_vars * sizeof(double));
    for (int v = 0; v < n_vars; ++v)
      atol[v] = o.atol;
  }
  int n_threads = 1;
#pragma omp parallel
  {
#pragma omp single
    n_threads = omp_get_num_threads();
  }
  volatile int failed = 0;
#pragma omp parallel for schedule(runtime)
  for (int g = 0; g < n_groups; ++g)
  {
    group *grp = &groups[g];
    grp->first_cell = g * o.group_size;
    grp->n_cells = g == n_groups - 1 ? n_cells - g * o.group_size : o.group_size;
    tmc_state *state = tmc_solver_create_state(solver, grp->n_cells);
    ts1_initial_state(state);
    tmc_conditions *cond = tmc_state_conditions(state);
    for (int c = 0; c < grp->n_cells; ++c)
    {
      const double *row = t.values + (size_t)(grp->first_cell + c) * (size_t)t.n_columns;
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
    if (o.rtol > 0.0)
      tmc_state_set_relative_tolerance(state, o.rtol);
    if (atol)
      tmc_state_set_absolute_tolerances(state, atol, n_vars, NULL);
    tmc_solver_update_state_parameters(solver, state);
    grp->base = malloc((size_t)grp->n_cells * (size_t)(n_photolysis > 0 ? n_photolysis : 1) * sizeof(double));
    tmc_dense *params = tmc_state_custom_rate_parameters(state);
    for (int c = 0; c < grp->n_cells; ++c)
      for (int j = 0; j < n_photolysis; ++j)
        grp->base[(size_t)c * n_photolysis + j] = tmc_dense_get(params, c, photolysis_index[j]);
    grp->snapshots = malloc((size_t)n_snapshots * (size_t)grp->n_cells * (size_t)n_vars * sizeof(double));
    grp->state = state;
  }
  tmc_stats *thread_stats = calloc((size_t)n_threads, sizeof *thread_stats);

  /* The timed steps */
  const double t0 = omp_get_wtime();
#pragma omp parallel for schedule(runtime)
  for (int g = 0; g < n_groups; ++g)
  {
    if (failed)
      continue;
    group *grp = &groups[g];
    tmc_state *state = grp->state;
    tmc_stats *stats = &thread_stats[omp_get_thread_num()];
    tmc_dense *params = tmc_state_custom_rate_parameters(state);
    double last_factor = -1.0;
    int next_snapshot = 0;
    for (int step = 0; step < sc.steps; ++step)
    {
      const double factor = photolysis_at(&ph, (double)step * sc.dt);
      if (factor != last_factor)
      {
        for (int c = 0; c < grp->n_cells; ++c)
          for (int j = 0; j < n_photolysis; ++j)
            tmc_dense_set(params, c, photolysis_index[j], grp->base[(size_t)c * n_photolysis + j] * factor);
        tmc_solver_update_custom_rate_constants(solver, state);
        grp->rate_updates += 1;
        last_factor = factor;
      }
      tmc_result r = tmc_solver_solve(solver, sc.dt, state);
      add_stats(stats, &r.stats);
      if (r.state != TMC_CONVERGED)
      {
        snprintf(grp->error, sizeof grp->error, "step %d of cells %d-%d: %s", step, grp->first_cell, grp->first_cell + grp->n_cells - 1, tmc_solver_state_name(r.state));
        failed = 1;
        break;
      }
      if (next_snapshot < n_snapshots && step + 1 == snapshot_steps[next_snapshot])
      {
        double *snap = grp->snapshots + (size_t)next_snapshot * (size_t)grp->n_cells * (size_t)n_vars;
        for (int c = 0; c < grp->n_cells; ++c)
          for (int v = 0; v < n_vars; ++v)
            snap[(size_t)c * n_vars + v] = tmc_state_get(state, c, v);
        ++next_snapshot;
      }
    }
  }
  const double t1 = omp_get_wtime();
  if (failed)
  {
    for (int g = 0; g < n_groups; ++g)
      if (groups[g].error[0])
        fprintf(stderr, "ts1_scenario: %s\n", groups[g].error);
    return 1;
  }

  tmc_stats total;
  memset(&total, 0, sizeof total);
  long rate_updates = 0;
  for (int k = 0; k < n_threads; ++k)
    add_stats(&total, &thread_stats[k]);
  for (int g = 0; g < n_groups; ++g)
    rate_updates += groups[g].rate_updates;

  /* The outputs */
  FILE *out = fopen("ts1_final.csv", "w");
  if (!out)
    fail("cannot write ts1_final.csv");
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
  if (sc.snapshot_every > 0)
  {
    out = fopen("ts1_snapshots.csv", "w");
    if (!out)
      fail("cannot write ts1_snapshots.csv");
    fprintf(out, "step,cell,temperature,pressure");
    for (int v = 0; v < n_vars; ++v)
      fprintf(out, ",%s", tmc_solver_variable_name(solver, v));
    fprintf(out, "\n");
    for (int sn = 0; sn < n_snapshots; ++sn)
      for (int g = 0; g < n_groups; ++g)
      {
        const tmc_conditions *cond = tmc_state_conditions(groups[g].state);
        const double *snap = groups[g].snapshots + (size_t)sn * (size_t)groups[g].n_cells * (size_t)n_vars;
        for (int c = 0; c < groups[g].n_cells; ++c)
        {
          fprintf(out, "%d,%d,%.17g,%.17g", snapshot_steps[sn], groups[g].first_cell + c, cond[c].temperature, cond[c].pressure);
          for (int v = 0; v < n_vars; ++v)
            fprintf(out, ",%.17g", snap[(size_t)c * n_vars + v]);
          fprintf(out, "\n");
        }
      }
    fclose(out);
  }

  printf("elapsed_ms=%.3f\n", (t1 - t0) * 1000.0);
  printf("internal_steps=%ld accepted=%ld rejected=%ld function_calls=%ld jacobian_updates=%ld decompositions=%ld solves=%ld\n",
         total.number_of_steps, total.accepted, total.rejected, total.function_calls, total.jacobian_updates, total.decompositions, total.solves);
  printf("matrix_reuses=%ld rate_updates=%ld groups=%d threads=%d\n", total.matrix_reuses, rate_updates, n_groups, n_threads);

  for (int g = 0; g < n_groups; ++g)
  {
    tmc_state_free(groups[g].state);
    free(groups[g].base);
    free(groups[g].snapshots);
  }
  free(groups);
  free(thread_stats);
  free(atol);
  free(column_kind);
  free(column_index);
  free(photolysis_index);
  free(snapshot_steps);
  free(ph.time);
  free(ph.factor);
  free(t.names);
  free(t.values);
  tmc_solver_free(solver);
  tmc_system_free(&system);
  return 0;
}
