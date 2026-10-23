#include "chemistry.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

extern const chem_mechanism chem_generated_mechanism;

typedef struct {
  char **names;
  double **rows;
  size_t columns;
  size_t row_count;
} csv_table;

typedef enum { COLUMN_UNKNOWN, COLUMN_TEMPERATURE, COLUMN_PRESSURE,
               COLUMN_SPECIES, COLUMN_PARAMETER } column_kind;

typedef struct {
  column_kind kind;
  size_t index;
} column_map;

static void free_table(csv_table *table) {
  if (!table) return;
  for (size_t i = 0; i < table->columns; ++i) free(table->names[i]);
  for (size_t i = 0; i < table->row_count; ++i) free(table->rows[i]);
  free(table->names);
  free(table->rows);
  memset(table, 0, sizeof(*table));
}

static int split_header(char *line, csv_table *table) {
  size_t count = 1;
  for (const char *p = line; *p; ++p) if (*p == ',') ++count;
  table->names = calloc(count, sizeof(*table->names));
  if (!table->names) return 1;
  table->columns = count;
  size_t index = 0;
  char *start = line;
  for (char *p = line;; ++p) {
    if (*p == ',' || *p == '\0' || *p == '\n' || *p == '\r') {
      char end = *p;
      *p = '\0';
      table->names[index] = malloc(strlen(start) + 1);
      if (!table->names[index]) return 1;
      strcpy(table->names[index], start);
      ++index;
      if (end != ',') break;
      start = p + 1;
    }
  }
  return index == count ? 0 : 1;
}

static int read_csv(const char *path, csv_table *table) {
  FILE *file = fopen(path, "r");
  if (!file) {
    fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
    return 1;
  }
  char *line = malloc(262144);
  if (!line) { fclose(file); return 1; }
  if (!fgets(line, 262144, file) || split_header(line, table)) {
    fprintf(stderr, "%s: missing or invalid header\n", path);
    free(line); fclose(file); return 1;
  }
  size_t capacity = 0;
  while (fgets(line, 262144, file)) {
    if (line[0] == '\n' || line[0] == '\r' || line[0] == '\0') continue;
    if (table->row_count == capacity) {
      size_t next = capacity ? capacity * 2 : 64;
      double **rows = realloc(table->rows, next * sizeof(*rows));
      if (!rows) { free(line); fclose(file); return 1; }
      table->rows = rows;
      capacity = next;
    }
    double *row = calloc(table->columns, sizeof(*row));
    if (!row) { free(line); fclose(file); return 1; }
    size_t column = 0;
    char *start = line;
    for (char *p = line;; ++p) {
      if (*p == ',' || *p == '\0' || *p == '\n' || *p == '\r') {
        char end = *p;
        *p = '\0';
        if (column >= table->columns) {
          fprintf(stderr, "%s: row %zu has too many fields\n", path, table->row_count + 2);
          free(row); free(line); fclose(file); return 1;
        }
        char *tail = NULL;
        errno = 0;
        row[column] = strtod(start, &tail);
        if (tail == start || errno == ERANGE || *tail != '\0' || !isfinite(row[column])) {
          fprintf(stderr, "%s: invalid number in row %zu column %zu\n", path,
                  table->row_count + 2, column + 1);
          free(row); free(line); fclose(file); return 1;
        }
        ++column;
        if (end != ',') break;
        start = p + 1;
      }
    }
    if (column != table->columns) {
      fprintf(stderr, "%s: row %zu has %zu fields; expected %zu\n", path,
              table->row_count + 2, column, table->columns);
      free(row); free(line); fclose(file); return 1;
    }
    table->rows[table->row_count++] = row;
  }
  free(line);
  fclose(file);
  return 0;
}

static size_t species_index(const chem_mechanism *m, const char *name) {
  for (size_t i = 0; i < m->species_count; ++i)
    if (strcmp(m->species[i].name, name) == 0) return i;
  return (size_t)-1;
}

static size_t parameter_index(const chem_mechanism *m, const char *name) {
  for (size_t i = 0; i < m->custom_parameter_count; ++i)
    if (strcmp(m->custom_labels[i], name) == 0) return i;
  return (size_t)-1;
}

static int map_columns(const chem_mechanism *m, const csv_table *table,
                       column_map **result) {
  column_map *map = calloc(table->columns, sizeof(*map));
  if (!map) return 1;
  for (size_t i = 0; i < table->columns; ++i) {
    if (strcmp(table->names[i], "temperature") == 0) {
      map[i].kind = COLUMN_TEMPERATURE;
    } else if (strcmp(table->names[i], "pressure") == 0) {
      map[i].kind = COLUMN_PRESSURE;
    } else {
      size_t species = species_index(m, table->names[i]);
      size_t parameter = parameter_index(m, table->names[i]);
      if (species != (size_t)-1 && !m->species[species].third_body) {
        map[i].kind = COLUMN_SPECIES;
        map[i].index = species;
      } else if (parameter != (size_t)-1) {
        map[i].kind = COLUMN_PARAMETER;
        map[i].index = parameter;
      } else {
        fprintf(stderr, "unknown TS1 CSV column '%s'\n", table->names[i]);
        free(map);
        return 1;
      }
    }
    for (size_t j = 0; j < i; ++j) {
      if (strcmp(table->names[i], table->names[j]) == 0) {
        fprintf(stderr, "duplicate TS1 CSV column '%s'\n", table->names[i]);
        free(map);
        return 1;
      }
    }
  }
  *result = map;
  return 0;
}

static int parse_size(const char *text, size_t *value) {
  char *tail = NULL;
  errno = 0;
  unsigned long long parsed = strtoull(text, &tail, 10);
  if (tail == text || *tail != '\0' || errno || parsed == 0 ||
      parsed > (unsigned long long)(size_t)-1) return 1;
  *value = (size_t)parsed;
  return 0;
}

static int parse_double(const char *text, double *value) {
  char *tail = NULL;
  errno = 0;
  *value = strtod(text, &tail);
  return tail == text || *tail != '\0' || errno || !isfinite(*value);
}

static double now_seconds(void) {
  struct timespec ts;
  timespec_get(&ts, TIME_UTC);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1.0e-9;
}

static void destroy_workspaces(chem_solver_workspace *workspaces, size_t count) {
  if (!workspaces) return;
  for (size_t i = 0; i < count; ++i) chem_solver_workspace_destroy(&workspaces[i]);
  free(workspaces);
}

typedef struct {
  char *cells_path;
  char *photolysis_path;
  double dt;
  size_t steps;
  size_t snapshot_every;
} scenario_config;

static char *trim(char *text) {
  while (*text == ' ' || *text == '\t') ++text;
  char *end = text + strlen(text);
  while (end > text && (end[-1] == ' ' || end[-1] == '\t' ||
                        end[-1] == '\r' || end[-1] == '\n')) --end;
  *end = '\0';
  return text;
}

static char *scenario_relative_path(const char *scenario, const char *leaf) {
  const char *slash = strrchr(scenario, '/');
  size_t directory_length = slash ? (size_t)(slash - scenario + 1) : 0;
  size_t leaf_length = strlen(leaf);
  char *path = malloc(directory_length + leaf_length + 1);
  if (!path) return NULL;
  if (directory_length) memcpy(path, scenario, directory_length);
  memcpy(path + directory_length, leaf, leaf_length + 1);
  return path;
}

static void free_scenario(scenario_config *scenario) {
  free(scenario->cells_path);
  free(scenario->photolysis_path);
  memset(scenario, 0, sizeof(*scenario));
}

static int read_scenario(const char *path, scenario_config *scenario) {
  FILE *file = fopen(path, "r");
  if (!file) {
    fprintf(stderr, "cannot open scenario %s: %s\n", path, strerror(errno));
    return 1;
  }
  char *line = malloc(4096);
  if (!line) { fclose(file); return 1; }
  int failed = 0;
  while (fgets(line, 4096, file)) {
    char *comment = strchr(line, '#');
    if (comment) *comment = '\0';
    char *entry = trim(line);
    if (!*entry) continue;
    char *equals = strchr(entry, '=');
    if (!equals) { failed = 1; break; }
    *equals = '\0';
    char *key = trim(entry), *value = trim(equals + 1);
    if (strcmp(key, "cells") == 0)
      scenario->cells_path = scenario_relative_path(path, value);
    else if (strcmp(key, "photolysis") == 0)
      scenario->photolysis_path = scenario_relative_path(path, value);
    else if (strcmp(key, "dt") == 0)
      failed = parse_double(value, &scenario->dt) || scenario->dt <= 0.0;
    else if (strcmp(key, "steps") == 0)
      failed = parse_size(value, &scenario->steps);
    else if (strcmp(key, "snapshot_every") == 0)
      failed = parse_size(value, &scenario->snapshot_every);
    if (failed) break;
  }
  free(line);
  fclose(file);
  if (!scenario->cells_path || !scenario->photolysis_path ||
      scenario->dt <= 0.0 || scenario->steps == 0 || failed) {
    fprintf(stderr, "%s: scenario needs valid cells, photolysis, dt, and steps\n", path);
    free_scenario(scenario);
    return 1;
  }
  return 0;
}

static double scenario_factor(const csv_table *table, double time) {
  if (table->row_count == 0 || table->columns != 2) return NAN;
  if (time <= table->rows[0][0]) return table->rows[0][1];
  for (size_t i = 1; i < table->row_count; ++i) {
    if (time <= table->rows[i][0]) {
      const double fraction = (time - table->rows[i-1][0]) /
          (table->rows[i][0] - table->rows[i-1][0]);
      return table->rows[i-1][1] +
          fraction * (table->rows[i][1] - table->rows[i-1][1]);
    }
  }
  return table->rows[table->row_count-1][1];
}

static void write_snapshot_row(FILE *file, size_t step, size_t cell,
                               const chem_mechanism *m, const double *state,
                               double temperature, double pressure) {
  fprintf(file, "%zu,%zu,%.17g,%.17g", step, cell, temperature, pressure);
  for (size_t s = 0; s < m->species_count; ++s)
    if (!m->species[s].third_body) fprintf(file, ",%.17g", state[s]);
  fputc('\n', file);
}

static int run_scenario(const char *scenario_path, size_t threads,
                        size_t requested_group_size, double rtol,
                        double atol) {
  scenario_config scenario = {0};
  csv_table input = {0}, factors = {0};
  column_map *columns = NULL;
  chem_solver_workspace *workspaces = NULL;
  size_t workspace_count = 0;
  double *states = NULL, *base_rates = NULL, *rates = NULL;
  double *temperature = NULL, *pressure = NULL, *density = NULL;
  chem_solver_stats *group_stats = NULL;
  int *status = NULL;
  int result = 1;
  const chem_mechanism *m = &chem_generated_mechanism;
  const size_t n = m->species_count, nr = m->custom_parameter_count;
  if (read_scenario(scenario_path, &scenario) ||
      read_csv(scenario.cells_path, &input) ||
      read_csv(scenario.photolysis_path, &factors)) goto cleanup;
  const size_t cells = input.row_count;
  if (cells == 0 || factors.columns != 2 || factors.row_count == 0 ||
      map_columns(m, &input, &columns)) goto cleanup;
  if (cells > (size_t)-1 / n ||
      cells > (size_t)-1 / (nr ? nr : 1) ||
      cells*n > (size_t)-1 / sizeof(*states) ||
      cells*(nr ? nr : 1) > (size_t)-1 / sizeof(*rates)) {
    fprintf(stderr, "scenario cell count is too large\n");
    goto cleanup;
  }
  for (size_t i = 1; i < factors.row_count; ++i)
    if (!(factors.rows[i][0] > factors.rows[i-1][0])) {
      fprintf(stderr, "photolysis times must be strictly increasing\n");
      goto cleanup;
    }
  states = calloc(cells*n, sizeof(*states));
  base_rates = calloc(cells*(nr ? nr : 1), sizeof(*base_rates));
  rates = calloc(cells*(nr ? nr : 1), sizeof(*rates));
  temperature = malloc(cells*sizeof(*temperature));
  pressure = malloc(cells*sizeof(*pressure));
  density = malloc(cells*sizeof(*density));
  if (!states || !base_rates || !rates || !temperature || !pressure || !density)
    goto cleanup;
  for (size_t c = 0; c < cells; ++c) {
    temperature[c] = m->initial_temperature;
    pressure[c] = m->initial_pressure;
    for (size_t s = 0; s < n; ++s)
      states[c*n+s] = m->species[s].initial_concentration;
    if (nr) memcpy(base_rates+c*nr, m->default_custom_parameters, nr*sizeof(double));
    for (size_t j = 0; j < input.columns; ++j) {
      const double value = input.rows[c][j];
      switch (columns[j].kind) {
      case COLUMN_TEMPERATURE: temperature[c] = value; break;
      case COLUMN_PRESSURE: pressure[c] = value; break;
      case COLUMN_SPECIES: states[c*n+columns[j].index] = value; break;
      case COLUMN_PARAMETER: base_rates[c*nr+columns[j].index] = value; break;
      default: break;
      }
    }
    if (!(temperature[c] > 0.0) || pressure[c] < 0.0) goto cleanup;
    density[c] = pressure[c] / (8.31446261815324 * temperature[c]);
  }
  if (threads > cells) threads = cells;
  const size_t group_size = requested_group_size ? requested_group_size :
      cells / threads + (cells % threads != 0);
  const size_t groups = cells / group_size + (cells % group_size != 0);
  const size_t capacity = group_size < cells ? group_size : cells;
  workspaces = calloc(groups, sizeof(*workspaces));
  group_stats = calloc(groups, sizeof(*group_stats));
  status = calloc(groups, sizeof(*status));
  if (!workspaces || !group_stats || !status) goto cleanup;
  workspace_count = groups;
  for (size_t g = 0; g < groups; ++g)
    if (chem_solver_workspace_init_group(&workspaces[g], m, capacity)) goto cleanup;
  chem_solver_options options = chem_solver_default_options();
  options.relative_tolerance = rtol;
  options.absolute_tolerance = atol;
  FILE *snapshots = fopen("ts1_snapshots.csv", "w");
  if (!snapshots) goto cleanup;
  fputs("step,cell,temperature,pressure", snapshots);
  for (size_t s = 0; s < n; ++s)
    if (!m->species[s].third_body) fprintf(snapshots, ",%s", m->species[s].name);
  fputc('\n', snapshots);
  double excluded_io = 0.0;
  double start = now_seconds();
  for (size_t step = 0; step < scenario.steps; ++step) {
    const double factor = scenario_factor(&factors, (double)step*scenario.dt);
    if (!isfinite(factor)) { fclose(snapshots); goto cleanup; }
    for (size_t c = 0; c < cells; ++c) {
      if (nr) memcpy(rates+c*nr, base_rates+c*nr, nr*sizeof(double));
      for (size_t p = 0; p < nr; ++p)
        if (m->custom_labels[p][0] == 'j') rates[c*nr+p] *= factor;
    }
    #pragma omp parallel for schedule(static) num_threads(threads)
    for (size_t g = 0; g < groups; ++g) {
      const size_t first = g*group_size;
      const size_t count = cells-first < group_size ? cells-first : group_size;
      status[g] = chem_rosenbrock3_solve_group_with_workspace(
          m, &workspaces[g], states+first*n, temperature+first, pressure+first,
          density+first, nr ? rates+first*nr : NULL, count, scenario.dt,
          &options, &group_stats[g]);
    }
    for (size_t g = 0; g < groups; ++g) if (status[g]) {
      fprintf(stderr, "scenario group %zu failed with status %d\n", g, status[g]);
      fclose(snapshots); goto cleanup;
    }
    const size_t completed = step + 1;
    if ((scenario.snapshot_every &&
         completed % scenario.snapshot_every == 0) || completed == scenario.steps) {
      const double before_io = now_seconds();
      for (size_t c = 0; c < cells; ++c)
        write_snapshot_row(snapshots, completed, c, m, states+c*n,
                           temperature[c], pressure[c]);
      excluded_io += now_seconds() - before_io;
    }
  }
  const double elapsed_ms = (now_seconds() - start - excluded_io)*1000.0;
  fclose(snapshots);
  chem_solver_stats total = {0};
  for (size_t g = 0; g < groups; ++g) {
    total.internal_steps += group_stats[g].internal_steps;
    total.accepted += group_stats[g].accepted;
    total.rejected += group_stats[g].rejected;
    total.function_calls += group_stats[g].function_calls;
    total.jacobian_updates += group_stats[g].jacobian_updates;
    total.decompositions += group_stats[g].decompositions;
    total.solves += group_stats[g].solves;
  }
  FILE *final = fopen("ts1_final.csv", "w");
  if (!final) goto cleanup;
  fputs("temperature,pressure", final);
  for (size_t s = 0; s < n; ++s)
    if (!m->species[s].third_body) fprintf(final, ",%s", m->species[s].name);
  fputc('\n', final);
  for (size_t c = 0; c < cells; ++c) {
    fprintf(final, "%.17g,%.17g", temperature[c], pressure[c]);
    for (size_t s = 0; s < n; ++s)
      if (!m->species[s].third_body) fprintf(final, ",%.17g", states[c*n+s]);
    fputc('\n', final);
  }
  fclose(final);
  printf("elapsed_ms=%.3f\ninternal_steps=%zu accepted=%zu rejected=%zu function_calls=%zu jacobian_updates=%zu decompositions=%zu solves=%zu\n",
         elapsed_ms, total.internal_steps, total.accepted, total.rejected,
         total.function_calls, total.jacobian_updates, total.decompositions,
         total.solves);
  result = 0;
cleanup:
  destroy_workspaces(workspaces, workspace_count);
  free(group_stats); free(status); free(states); free(base_rates); free(rates);
  free(temperature); free(pressure); free(density); free(columns);
  free_table(&input); free_table(&factors); free_scenario(&scenario);
  return result;
}

int main(int argc, char **argv) {
  size_t threads = 1;
  size_t group_size = 0;
  double rtol = 1.0e-6, atol = 1.0e-3;
#ifdef _OPENMP
  threads = (size_t)omp_get_max_threads();
#endif
  const char *positional[4];
  size_t npos = 0;
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
      if (parse_size(argv[++i], &threads)) {
        fprintf(stderr, "--threads expects a positive integer\n");
        return 2;
      }
    } else if (strcmp(argv[i], "--group-size") == 0 && i + 1 < argc) {
      if (parse_size(argv[++i], &group_size)) {
        fprintf(stderr, "--group-size expects a positive integer\n");
        return 2;
      }
    } else if (strcmp(argv[i], "--rtol") == 0 && i + 1 < argc) {
      if (parse_double(argv[++i], &rtol) || rtol <= 0.0) {
        fprintf(stderr, "--rtol expects a positive number\n");
        return 2;
      }
    } else if (strcmp(argv[i], "--atol") == 0 && i + 1 < argc) {
      if (parse_double(argv[++i], &atol) || atol <= 0.0) {
        fprintf(stderr, "--atol expects a positive number\n");
        return 2;
      }
    } else if (npos < 4) {
      positional[npos++] = argv[i];
    } else {
      fprintf(stderr, "usage: %s [--threads N] CELLS STEPS DT CELLS_CSV | [--threads N] [--group-size N] [--rtol R] [--atol A] SCENARIO\n", argv[0]);
      return 2;
    }
  }
  if (npos == 1)
    return run_scenario(positional[0], threads, group_size, rtol, atol);
  if (npos != 4) {
    fprintf(stderr, "usage: %s [--threads N] CELLS STEPS DT CELLS_CSV | [--threads N] [--group-size N] [--rtol R] [--atol A] SCENARIO\n", argv[0]);
    return 2;
  }
  size_t cells, steps;
  double dt;
  if (parse_size(positional[0], &cells) || parse_size(positional[1], &steps) ||
      parse_double(positional[2], &dt) || dt <= 0.0) {
    fprintf(stderr, "CELLS and STEPS must be positive integers; DT must be positive\n");
    return 2;
  }
  csv_table input = {0};
  if (read_csv(positional[3], &input)) { free_table(&input); return 1; }
  if (input.row_count < cells) {
    fprintf(stderr, "CELLS is %zu but CELLS_CSV has %zu data rows\n", cells, input.row_count);
    free_table(&input); return 1;
  }
  column_map *columns = NULL;
  if (map_columns(&chem_generated_mechanism, &input, &columns)) {
    free_table(&input); return 1;
  }
  const chem_mechanism *m = &chem_generated_mechanism;
  const size_t n = m->species_count, nr = m->custom_parameter_count;
  if (cells > (size_t)-1 / n || cells > (size_t)-1 / nr) {
    fprintf(stderr, "requested cell count is too large\n");
    free(columns); free_table(&input); return 1;
  }
  double *states = calloc(cells*n, sizeof(*states));
  double *rates = calloc(cells*nr, sizeof(*rates));
  double *temperature = malloc(cells*sizeof(*temperature));
  double *pressure = malloc(cells*sizeof(*pressure));
  double *density = malloc(cells*sizeof(*density));
  chem_solver_stats *stats = calloc(cells, sizeof(*stats));
  int *status = calloc(cells, sizeof(*status));
  if (!states || !rates || !temperature || !pressure || !density || !stats || !status) {
    fprintf(stderr, "allocation failed for %zu cells\n", cells);
    free(states); free(rates); free(temperature); free(pressure); free(density); free(stats); free(status);
    free(columns); free_table(&input); return 1;
  }
  for (size_t c = 0; c < cells; ++c) {
    temperature[c] = m->initial_temperature;
    pressure[c] = m->initial_pressure;
    for (size_t s = 0; s < n; ++s)
      states[c*n+s] = m->species[s].initial_concentration;
    memcpy(rates + c*nr, m->default_custom_parameters, nr*sizeof(double));
    for (size_t j = 0; j < input.columns; ++j) {
      const double value = input.rows[c][j];
      switch (columns[j].kind) {
      case COLUMN_TEMPERATURE: temperature[c] = value; break;
      case COLUMN_PRESSURE: pressure[c] = value; break;
      case COLUMN_SPECIES: states[c*n + columns[j].index] = value; break;
      case COLUMN_PARAMETER: rates[c*nr + columns[j].index] = value; break;
      default: break;
      }
    }
    if (temperature[c] <= 0.0 || pressure[c] < 0.0) {
      fprintf(stderr, "cell %zu has invalid temperature or pressure\n", c);
      status[c] = 1;
    } else {
      density[c] = pressure[c] / (8.31446261815324 * temperature[c]);
    }
  }
  for (size_t c = 0; c < cells; ++c) {
    if (status[c]) {
      free(columns); free_table(&input);
      free(states); free(rates); free(temperature); free(pressure);
      free(density); free(stats); free(status);
      return 1;
    }
  }
  free(columns);
  free_table(&input);
  if (threads > cells) threads = cells;
  const size_t group_capacity = (cells + threads - 1) / threads;
  chem_solver_options options = chem_solver_default_options();
  chem_solver_workspace *workspaces = calloc(threads, sizeof(*workspaces));
  if (!workspaces) {
    fprintf(stderr, "allocation failed for solver workspaces\n");
    free(states); free(rates); free(temperature); free(pressure); free(density); free(stats); free(status);
    return 1;
  }
  for (size_t i = 0; i < threads; ++i) {
    if (chem_solver_workspace_init_group(&workspaces[i], m, group_capacity)) {
      fprintf(stderr, "allocation failed for solver workspace %zu\n", i);
      destroy_workspaces(workspaces, threads);
      free(states); free(rates); free(temperature); free(pressure); free(density); free(stats); free(status);
      return 1;
    }
  }
  #pragma omp parallel for schedule(static) num_threads(threads)
  for (size_t worker = 0; worker < threads; ++worker) {
    const size_t first = cells*worker/threads;
    const size_t last = cells*(worker+1)/threads;
    for (size_t c = first; c < last; ++c) if (status[c]) status[worker] = 1;
    if (!status[worker]) {
      chem_solver_stats warmup_stats;
      status[worker] = chem_rosenbrock3_solve_group_with_workspace(
          m, &workspaces[worker], states+first*n, temperature+first,
          pressure+first, density+first, nr ? rates+first*nr : NULL,
          last-first, dt, &options, &warmup_stats);
    }
  }
  for (size_t worker = 0; worker < threads; ++worker) {
    if (status[worker]) {
      fprintf(stderr, "group %zu failed warm-up solve with status %d\n", worker, status[worker]);
      destroy_workspaces(workspaces, threads);
      free(states); free(rates); free(temperature); free(pressure); free(density); free(stats); free(status);
      return 1;
    }
  }
  const double start = now_seconds();
  #pragma omp parallel for schedule(static) num_threads(threads)
  for (size_t worker = 0; worker < threads; ++worker) {
    const size_t first = cells*worker/threads;
    const size_t last = cells*(worker+1)/threads;
    for (size_t step = 0; step < steps && status[worker] == 0; ++step) {
        chem_solver_stats one = {0};
        status[worker] = chem_rosenbrock3_solve_group_with_workspace(
            m, &workspaces[worker], states+first*n, temperature+first,
            pressure+first, density+first, nr ? rates+first*nr : NULL,
            last-first, dt, &options, &one);
        stats[worker].internal_steps += one.internal_steps;
        stats[worker].accepted += one.accepted;
        stats[worker].rejected += one.rejected;
        stats[worker].function_calls += one.function_calls;
        stats[worker].jacobian_updates += one.jacobian_updates;
        stats[worker].decompositions += one.decompositions;
        stats[worker].solves += one.solves;
    }
  }
  const double elapsed_ms = (now_seconds() - start)*1000.0;
  chem_solver_stats total = {0};
  for (size_t worker = 0; worker < threads; ++worker) {
    if (status[worker]) {
      fprintf(stderr, "group %zu failed integration with status %d\n", worker, status[worker]);
      destroy_workspaces(workspaces, threads);
      free(states); free(rates); free(temperature); free(pressure); free(density); free(stats); free(status);
      return 1;
    }
    total.internal_steps += stats[worker].internal_steps;
    total.accepted += stats[worker].accepted;
    total.rejected += stats[worker].rejected;
    total.function_calls += stats[worker].function_calls;
    total.jacobian_updates += stats[worker].jacobian_updates;
    total.decompositions += stats[worker].decompositions;
    total.solves += stats[worker].solves;
  }
  FILE *output = fopen("ts1_final.csv", "w");
  if (!output) {
    fprintf(stderr, "cannot write ts1_final.csv: %s\n", strerror(errno));
    destroy_workspaces(workspaces, threads);
    free(states); free(rates); free(temperature); free(pressure); free(density); free(stats); free(status);
    return 1;
  }
  fputs("temperature,pressure", output);
  for (size_t s = 0; s < n; ++s)
    if (!m->species[s].third_body) fprintf(output, ",%s", m->species[s].name);
  fputc('\n', output);
  for (size_t c = 0; c < cells; ++c) {
    fprintf(output, "%.17g,%.17g", temperature[c], pressure[c]);
    for (size_t s = 0; s < n; ++s)
      if (!m->species[s].third_body) fprintf(output, ",%.17g", states[c*n+s]);
    fputc('\n', output);
  }
  fclose(output);
  printf("elapsed_ms=%.3f\ninternal_steps=%zu accepted=%zu rejected=%zu function_calls=%zu jacobian_updates=%zu decompositions=%zu solves=%zu\n",
         elapsed_ms, total.internal_steps, total.accepted, total.rejected,
         total.function_calls, total.jacobian_updates, total.decompositions,
         total.solves);
  free(states); free(rates); free(temperature); free(pressure); free(density); free(stats); free(status);
  destroy_workspaces(workspaces, threads);
  return 0;
}
