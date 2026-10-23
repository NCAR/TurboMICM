#include "chem.h"
#ifdef CHEM_USE_SPECIALIZED
#include "chem_specialized.h"
#endif

#include <ctype.h>
#include <errno.h>
#include <omp.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
  char **names;
  double **rows;
  size_t column_count;
  size_t row_count;
} csv_table;

static void csv_destroy(csv_table *table) {
  if (!table)
    return;
  for (size_t i = 0; i < table->column_count; ++i)
    free(table->names[i]);
  for (size_t i = 0; i < table->row_count; ++i)
    free(table->rows[i]);
  free(table->names);
  free(table->rows);
  memset(table, 0, sizeof(*table));
}

static char *copy_string(const char *start, size_t length) {
  char *value = malloc(length + 1);
  if (!value)
    return NULL;
  memcpy(value, start, length);
  value[length] = '\0';
  return value;
}

static char *read_line(FILE *file, size_t *length) {
  size_t capacity = 128, used = 0;
  char *line = malloc(capacity);
  if (!line)
    return NULL;
  int ch;
  while ((ch = fgetc(file)) != EOF) {
    if (used + 1 >= capacity) {
      if (capacity > (size_t)-1 / 2) {
        free(line);
        return NULL;
      }
      capacity *= 2;
      char *grown = realloc(line, capacity);
      if (!grown) {
        free(line);
        return NULL;
      }
      line = grown;
    }
    line[used++] = (char)ch;
    if (ch == '\n')
      break;
  }
  if (used == 0 && ch == EOF) {
    free(line);
    return NULL;
  }
  line[used] = '\0';
  *length = used;
  return line;
}

static size_t split_csv(char *line, char **fields, size_t capacity) {
  size_t count = 0;
  char *start = line;
  for (char *p = line;; ++p) {
    if (*p == ',' || *p == '\0' || *p == '\n' || *p == '\r') {
      if (count >= capacity)
        return capacity + 1;
      fields[count++] = copy_string(start, (size_t)(p - start));
      if (!fields[count - 1])
        return capacity + 1;
      if (*p == '\0' || *p == '\n' || *p == '\r')
        break;
      start = p + 1;
    }
  }
  return count;
}

static int read_csv(const char *path, csv_table *table) {
  FILE *file = fopen(path, "r");
  if (!file) {
    fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
    return -1;
  }

  size_t length = 0;
  char *line = read_line(file, &length);
  if (!line) {
    fprintf(stderr, "%s is empty or could not be read\n", path);
    fclose(file);
    return -1;
  }
  char **fields = malloc((length + 1u) * sizeof(char *));
  if (!fields) {
    fclose(file);
    free(line);
    return -1;
  }
  size_t count = split_csv(line, fields, length + 1u);
  if (!count || count > length + 1u) {
    fprintf(stderr, "invalid CSV header in %s\n", path);
    free(fields);
    fclose(file);
    free(line);
    return -1;
  }
  table->names = calloc(count, sizeof(char *));
  if (!table->names) {
    for (size_t i = 0; i < count; ++i)
      free(fields[i]);
    free(fields);
    fclose(file);
    free(line);
    return -1;
  }
  for (size_t i = 0; i < count; ++i)
    table->names[i] = fields[i];
  table->column_count = count;
  table->rows = NULL;
  table->row_count = 0;
  free(fields);
  free(line);

  size_t row_capacity = 0;
  while ((line = read_line(file, &length)) != NULL) {
    if (length == 0 || line[0] == '\n' || line[0] == '\r') {
      free(line);
      continue;
    }
    fields = malloc(count * sizeof(char *));
    if (!fields) {
      free(line);
      goto error;
    }
    const size_t actual = split_csv(line, fields, count);
    if (actual != count) {
      for (size_t i = 0; i < count; ++i)
        free(fields[i]);
      free(fields);
      free(line);
      fprintf(stderr, "%s row %zu has %zu fields; expected %zu\n",
              path, table->row_count + 2, actual, count);
      goto error;
    }
    double *row = malloc(count * sizeof(double));
    if (!row) {
      for (size_t i = 0; i < count; ++i)
        free(fields[i]);
      free(fields);
      free(line);
      goto error;
    }
    int valid = 1;
    for (size_t i = 0; i < count; ++i) {
      char *end = NULL;
      row[i] = strtod(fields[i], &end);
      if (end == fields[i] || *end != '\0') {
        fprintf(stderr, "%s row %zu column %s is not numeric\n",
                path, table->row_count + 2, table->names[i]);
        valid = 0;
      }
    }
    for (size_t i = 0; i < count; ++i)
      free(fields[i]);
    free(fields);
    free(line);
    if (!valid) {
      free(row);
      goto error;
    }
    if (table->row_count == row_capacity) {
      row_capacity = row_capacity ? row_capacity * 2 : 16;
      double **grown = realloc(table->rows, row_capacity * sizeof(double *));
      if (!grown) {
        free(row);
        goto error;
      }
      table->rows = grown;
    }
    table->rows[table->row_count++] = row;
  }
  fclose(file);
  return 0;

error:
  fclose(file);
  csv_destroy(table);
  return -1;
}

static int species_index(const chem_mechanism *mechanism, const char *name) {
  for (size_t i = 0; i < mechanism->species_count; ++i)
    if (strcmp(mechanism->species[i].name, name) == 0)
      return (int)i;
  return -1;
}

static int parameter_index(const chem_mechanism *mechanism, const char *name) {
  for (size_t i = 0; i < mechanism->parameter_count; ++i)
    if (strcmp(mechanism->parameters[i].name, name) == 0)
      return (int)i;
  return -1;
}

static double parse_double(const char *text, const char *what) {
  char *end = NULL;
  errno = 0;
  const double value = strtod(text, &end);
  if (errno || end == text || *end != '\0') {
    fprintf(stderr, "invalid %s: %s\n", what, text);
    exit(2);
  }
  return value;
}

static int write_final(
    const chem_mechanism *mechanism,
    size_t cells,
    const double *state,
    const double *temperature,
    const double *pressure) {
  FILE *file = fopen("ts1_final.csv", "w");
  if (!file) {
    fprintf(stderr, "cannot write ts1_final.csv: %s\n", strerror(errno));
    return -1;
  }
  fputs("temperature,pressure", file);
  for (size_t s = 0; s < mechanism->species_count; ++s)
    fprintf(file, ",%s", mechanism->species[s].name);
  fputc('\n', file);
  for (size_t c = 0; c < cells; ++c) {
    fprintf(file, "%.17g,%.17g", temperature[c], pressure[c]);
    for (size_t s = 0; s < mechanism->species_count; ++s)
      fprintf(file, ",%.17g", state[c * mechanism->species_count + s]);
    fputc('\n', file);
  }
  return fclose(file);
}

static double wall_time(void) {
  struct timespec value;
  timespec_get(&value, TIME_UTC);
  return (double)value.tv_sec + (double)value.tv_nsec * 1.0e-9;
}

#ifndef CHEM_SCENARIO
static int integrate_group(
    const chem_mechanism *mechanism, size_t cells, double *state, double *parameters,
    const double *temperature, const double *pressure, double interval, size_t intervals,
    chem_statistics *statistics) {
#ifdef CHEM_USE_SPECIALIZED
  (void)mechanism;
  return chem_generated_integrate_group_steps(
      cells, state, parameters, temperature, pressure, interval, intervals, -1.0, -1.0, statistics);
#else
  if (intervals == 1u)
    return chem_integrate_group(mechanism, cells, state, parameters, temperature, pressure, interval, statistics);
  return chem_integrate_group_steps(
      mechanism, cells, state, parameters, temperature, pressure, interval, intervals, statistics);
#endif
}
#endif

#ifndef CHEM_SCENARIO
int main(int argc, char **argv) {
  int threads = omp_get_max_threads();
  int arg = 1;
  if (argc > 1 && strcmp(argv[1], "--threads") == 0) {
    if (argc < 3) {
      fprintf(stderr, "usage: %s [--threads N] CELLS STEPS DT CELLS_CSV\n", argv[0]);
      return 2;
    }
    threads = (int)parse_double(argv[2], "thread count");
    arg = 3;
  }
  if (argc - arg != 4 || threads < 1) {
    fprintf(stderr, "usage: %s [--threads N] CELLS STEPS DT CELLS_CSV\n", argv[0]);
    return 2;
  }
  const size_t cells = (size_t)parse_double(argv[arg], "cell count");
  const size_t steps = (size_t)parse_double(argv[arg + 1], "step count");
  const double interval = parse_double(argv[arg + 2], "time step");
  const char *csv_path = argv[arg + 3];
  if (!cells || interval <= 0.0) {
    fprintf(stderr, "cells and DT must be positive\n");
    return 2;
  }

  const chem_mechanism *mechanism = chem_generated_mechanism();
  csv_table table = { 0 };
  if (read_csv(csv_path, &table) != 0)
    return 2;
  if (table.row_count < cells) {
    fprintf(stderr, "%s has %zu rows; requested %zu cells\n", csv_path, table.row_count, cells);
    csv_destroy(&table);
    return 2;
  }

  double *state = calloc(cells * mechanism->species_count, sizeof(double));
  double *parameters = calloc(cells * (mechanism->parameter_count ? mechanism->parameter_count : 1u), sizeof(double));
  double *temperature = malloc(cells * sizeof(double));
  double *pressure = malloc(cells * sizeof(double));
  if (!state || !parameters || !temperature || !pressure) {
    fprintf(stderr, "not enough memory for %zu cells\n", cells);
    free(state); free(parameters); free(temperature); free(pressure); csv_destroy(&table);
    return 2;
  }
  for (size_t c = 0; c < cells; ++c) {
    temperature[c] = mechanism->temperature;
    pressure[c] = mechanism->pressure;
    for (size_t s = 0; s < mechanism->species_count; ++s)
      state[c * mechanism->species_count + s] = mechanism->species[s].initial_value;
    for (size_t p = 0; p < mechanism->parameter_count; ++p)
      parameters[c * mechanism->parameter_count + p] = mechanism->parameters[p].initial_value;
  }
  for (size_t col = 0; col < table.column_count; ++col) {
    const char *name = table.names[col];
    const int s = species_index(mechanism, name);
    const int p = parameter_index(mechanism, name);
    const int is_temperature = strcmp(name, "temperature") == 0;
    const int is_pressure = strcmp(name, "pressure") == 0;
    if (s < 0 && p < 0 && !is_temperature && !is_pressure) {
      fprintf(stderr, "unknown TS1 CSV column: %s\n", name);
      free(state); free(parameters); free(temperature); free(pressure); csv_destroy(&table);
      return 2;
    }
    for (size_t c = 0; c < cells; ++c) {
      const double value = table.rows[c][col];
      if (is_temperature) temperature[c] = value;
      else if (is_pressure) pressure[c] = value;
      else if (s >= 0) state[c * mechanism->species_count + (size_t)s] = value;
      else parameters[c * mechanism->parameter_count + (size_t)p] = value;
    }
  }
  csv_destroy(&table);

  const int group_count = threads < (int)cells ? threads : (int)cells;
  chem_statistics *group_stats = calloc((size_t)group_count, sizeof(*group_stats));
  int *group_status = calloc((size_t)group_count, sizeof(int));
  if (!group_stats || !group_status) {
    free(state); free(parameters); free(temperature); free(pressure); free(group_stats); free(group_status);
    return 2;
  }
  omp_set_dynamic(0);
  #pragma omp parallel for schedule(static) num_threads(threads)
  for (int g = 0; g < group_count; ++g) {
    const size_t first = cells * (size_t)g / (size_t)group_count;
    const size_t last = cells * (size_t)(g + 1) / (size_t)group_count;
    const size_t count = last - first;
    double *group_y = state + first * mechanism->species_count;
    double *group_parameters = mechanism->parameter_count ? parameters + first * mechanism->parameter_count : NULL;
    if (integrate_group(mechanism, count, group_y, group_parameters, temperature + first, pressure + first,
                        interval, 1u, NULL) != 0)
      group_status[g] = -1;
  }
  for (int g = 0; g < group_count; ++g) {
    if (group_status[g]) {
      fprintf(stderr, "warm-up integration failed for group %d\n", g);
      free(state); free(parameters); free(temperature); free(pressure); free(group_stats); free(group_status);
      return 1;
    }
  }

  const double start = wall_time();
  if (steps) {
    #pragma omp parallel for schedule(static) num_threads(threads)
    for (int g = 0; g < group_count; ++g) {
      const size_t first = cells * (size_t)g / (size_t)group_count;
      const size_t last = cells * (size_t)(g + 1) / (size_t)group_count;
      const size_t count = last - first;
      double *group_y = state + first * mechanism->species_count;
      double *group_parameters = mechanism->parameter_count ? parameters + first * mechanism->parameter_count : NULL;
      if (integrate_group(mechanism, count, group_y, group_parameters, temperature + first,
                          pressure + first, interval, steps, &group_stats[g]) != 0)
        group_status[g] = -1;
    }
  }
  const double elapsed_ms = (wall_time() - start) * 1000.0;
  for (int g = 0; g < group_count; ++g) {
    if (group_status[g]) {
      fprintf(stderr, "timed integration failed for group %d\n", g);
      free(state); free(parameters); free(temperature); free(pressure); free(group_stats); free(group_status);
      return 1;
    }
  }
  chem_statistics total = { 0 };
  for (int g = 0; g < group_count; ++g) {
    total.internal_steps += group_stats[g].internal_steps;
    total.accepted += group_stats[g].accepted;
    total.rejected += group_stats[g].rejected;
    total.function_calls += group_stats[g].function_calls;
    total.jacobian_updates += group_stats[g].jacobian_updates;
    total.decompositions += group_stats[g].decompositions;
    total.solves += group_stats[g].solves;
  }
  printf("elapsed_ms=%.3f\n", elapsed_ms);
  printf("internal_steps=%lu accepted=%lu rejected=%lu function_calls=%lu jacobian_updates=%lu decompositions=%lu solves=%lu\n",
         total.internal_steps, total.accepted, total.rejected, total.function_calls, total.jacobian_updates,
         total.decompositions, total.solves);
  const int output_status = write_final(mechanism, cells, state, temperature, pressure);
  free(state); free(parameters); free(temperature); free(pressure); free(group_stats); free(group_status);
  return output_status == 0 ? 0 : 1;
}
#else

static char *trim(char *value) {
  while (isspace((unsigned char)*value)) ++value;
  char *end = value + strlen(value);
  while (end > value && isspace((unsigned char)end[-1])) --end;
  *end = '\0';
  return value;
}

static char *join_path(const char *scenario, const char *relative) {
  if (relative[0] == '/') return copy_string(relative, strlen(relative));
  const char *slash = strrchr(scenario, '/');
  if (!slash) return copy_string(relative, strlen(relative));
  const size_t directory_length = (size_t)(slash - scenario) + 1u;
  const size_t relative_length = strlen(relative);
  char *joined = malloc(directory_length + relative_length + 1u);
  if (!joined) return NULL;
  memcpy(joined, scenario, directory_length);
  memcpy(joined + directory_length, relative, relative_length + 1u);
  return joined;
}

static int read_scenario(
    const char *path, char *cells_path, size_t cells_capacity, char *photolysis_path,
    size_t photolysis_capacity, double *dt, size_t *steps, size_t *snapshot_every) {
  FILE *file = fopen(path, "r");
  if (!file) {
    fprintf(stderr, "cannot open scenario %s: %s\n", path, strerror(errno));
    return -1;
  }
  char line[4096];
  int have_cells = 0, have_photolysis = 0, have_dt = 0, have_steps = 0;
  *snapshot_every = 0;
  while (fgets(line, sizeof(line), file)) {
    char *comment = strchr(line, '#');
    if (comment) *comment = '\0';
    char *entry = trim(line);
    if (!*entry) continue;
    char *equals = strchr(entry, '=');
    if (!equals) {
      fprintf(stderr, "invalid scenario line: %s\n", entry);
      fclose(file);
      return -1;
    }
    *equals = '\0';
    char *key = trim(entry), *value = trim(equals + 1);
    if (strcmp(key, "cells") == 0) {
      if (strlen(value) >= cells_capacity) { fclose(file); return -1; }
      strcpy(cells_path, value); have_cells = 1;
    } else if (strcmp(key, "photolysis") == 0) {
      if (strlen(value) >= photolysis_capacity) { fclose(file); return -1; }
      strcpy(photolysis_path, value); have_photolysis = 1;
    } else if (strcmp(key, "dt") == 0) {
      *dt = parse_double(value, "scenario time step"); have_dt = 1;
    } else if (strcmp(key, "steps") == 0) {
      *steps = (size_t)parse_double(value, "scenario step count"); have_steps = 1;
    } else if (strcmp(key, "snapshot_every") == 0) {
      *snapshot_every = (size_t)parse_double(value, "snapshot interval");
    }
  }
  fclose(file);
  if (!have_cells || !have_photolysis || !have_dt || !have_steps || *dt <= 0.0 || !*steps) {
    fputs("scenario needs cells, photolysis, positive dt, and positive steps\n", stderr);
    return -1;
  }
  return 0;
}

static double photolysis_factor(const csv_table *table, double time) {
  if (table->row_count == 1 || time <= table->rows[0][0]) return table->rows[0][1];
  if (time >= table->rows[table->row_count - 1u][0]) return table->rows[table->row_count - 1u][1];
  for (size_t row = 0; row + 1u < table->row_count; ++row) {
    const double left_time = table->rows[row][0], right_time = table->rows[row + 1u][0];
    if (time <= right_time) {
      const double fraction = (time - left_time) / (right_time - left_time);
      return table->rows[row][1] + fraction * (table->rows[row + 1u][1] - table->rows[row][1]);
    }
  }
  return table->rows[table->row_count - 1u][1];
}

static int write_snapshots(
    const chem_mechanism *mechanism, size_t cells, size_t snapshot_count, const size_t *snapshot_steps,
    const double *states, const double *temperature, const double *pressure) {
  FILE *file = fopen("ts1_snapshots.csv", "w");
  if (!file) { fprintf(stderr, "cannot write ts1_snapshots.csv: %s\n", strerror(errno)); return -1; }
  fputs("step,cell,temperature,pressure", file);
  for (size_t species = 0; species < mechanism->species_count; ++species)
    fprintf(file, ",%s", mechanism->species[species].name);
  fputc('\n', file);
  for (size_t snapshot = 0; snapshot < snapshot_count; ++snapshot) {
    const double *state = states + snapshot * cells * mechanism->species_count;
    for (size_t cell = 0; cell < cells; ++cell) {
      fprintf(file, "%zu,%zu,%.17g,%.17g", snapshot_steps[snapshot], cell, temperature[cell], pressure[cell]);
      for (size_t species = 0; species < mechanism->species_count; ++species)
        fprintf(file, ",%.17g", state[cell * mechanism->species_count + species]);
      fputc('\n', file);
    }
  }
  return fclose(file);
}

int main(int argc, char **argv) {
  int threads = omp_get_max_threads();
  size_t group_size = 32u;
  double relative_tolerance = -1.0, absolute_tolerance = -1.0;
  int arg = 1;
  while (arg < argc && argv[arg][0] == '-') {
    if (arg + 1 >= argc) break;
    const char *option = argv[arg++];
    const char *value = argv[arg++];
    if (strcmp(option, "--threads") == 0) threads = (int)parse_double(value, "thread count");
    else if (strcmp(option, "--group-size") == 0) group_size = (size_t)parse_double(value, "group size");
    else if (strcmp(option, "--rtol") == 0) relative_tolerance = parse_double(value, "relative tolerance");
    else if (strcmp(option, "--atol") == 0) absolute_tolerance = parse_double(value, "absolute tolerance");
    else { fprintf(stderr, "unknown option: %s\n", option); return 2; }
  }
  if (argc - arg != 1 || threads < 1 || !group_size || relative_tolerance < -1.0 || absolute_tolerance < -1.0) {
    fprintf(stderr, "usage: %s [--threads N] [--rtol R] [--atol A] [--group-size G] SCENARIO\n", argv[0]);
    return 2;
  }
  const char *scenario_path = argv[arg];
  char cells_relative[4096], photolysis_relative[4096];
  double interval = 0.0;
  size_t steps = 0, snapshot_every = 0;
  if (read_scenario(scenario_path, cells_relative, sizeof(cells_relative), photolysis_relative,
                    sizeof(photolysis_relative), &interval, &steps, &snapshot_every) != 0) return 2;
  char *cells_path = join_path(scenario_path, cells_relative);
  char *photolysis_path = join_path(scenario_path, photolysis_relative);
  if (!cells_path || !photolysis_path) { free(cells_path); free(photolysis_path); return 2; }
  const chem_mechanism *mechanism = chem_generated_mechanism();
  csv_table cells_table = { 0 }, photo_table = { 0 };
  if (read_csv(cells_path, &cells_table) || read_csv(photolysis_path, &photo_table)) {
    free(cells_path); free(photolysis_path); csv_destroy(&cells_table); csv_destroy(&photo_table); return 2;
  }
  free(cells_path); free(photolysis_path);
  if (cells_table.column_count < 2u || photo_table.column_count != 2u || !photo_table.row_count ||
      strcmp(photo_table.names[0], "time_s") != 0 || strcmp(photo_table.names[1], "factor") != 0) {
    fputs("scenario CSV headers are invalid\n", stderr);
    csv_destroy(&cells_table); csv_destroy(&photo_table); return 2;
  }
  for (size_t row = 1; row < photo_table.row_count; ++row)
    if (!(photo_table.rows[row][0] > photo_table.rows[row - 1u][0])) {
      fputs("photolysis times must increase\n", stderr);
      csv_destroy(&cells_table); csv_destroy(&photo_table); return 2;
    }
  const size_t cells = cells_table.row_count, n = mechanism->species_count, pcount = mechanism->parameter_count;
  if (!cells) {
    fputs("scenario cells CSV is empty\n", stderr);
    csv_destroy(&cells_table); csv_destroy(&photo_table); return 2;
  }
  if (group_size > cells) group_size = cells;
  double *state = calloc(cells * n, sizeof(double));
  double *parameters = calloc(cells * (pcount ? pcount : 1u), sizeof(double));
  double *photo_base = calloc(cells * (pcount ? pcount : 1u), sizeof(double));
  double *temperature = malloc(cells * sizeof(double)), *pressure = malloc(cells * sizeof(double));
  size_t photo_count = 0;
  for (size_t p = 0; p < pcount; ++p) if (mechanism->parameters[p].name[0] == 'j') ++photo_count;
  size_t *photo_indices = malloc((photo_count ? photo_count : 1u) * sizeof(size_t));
  size_t photo_pos = 0;
  for (size_t p = 0; p < pcount; ++p)
    if (mechanism->parameters[p].name[0] == 'j') photo_indices[photo_pos++] = p;
  if (!state || !parameters || !photo_base || !temperature || !pressure || !photo_indices) {
    fputs("not enough memory for scenario cells\n", stderr);
    free(state); free(parameters); free(photo_base); free(temperature); free(pressure); free(photo_indices);
    csv_destroy(&cells_table); csv_destroy(&photo_table); return 2;
  }
  for (size_t cell = 0; cell < cells; ++cell) {
    temperature[cell] = mechanism->temperature; pressure[cell] = mechanism->pressure;
    for (size_t species = 0; species < n; ++species)
      state[cell * n + species] = mechanism->species[species].initial_value;
    for (size_t parameter = 0; parameter < pcount; ++parameter)
      parameters[cell * pcount + parameter] = mechanism->parameters[parameter].initial_value;
  }
  for (size_t column = 0; column < cells_table.column_count; ++column) {
    const char *name = cells_table.names[column];
    const int species = species_index(mechanism, name), parameter = parameter_index(mechanism, name);
    const int is_temperature = strcmp(name, "temperature") == 0, is_pressure = strcmp(name, "pressure") == 0;
    if (species < 0 && parameter < 0 && !is_temperature && !is_pressure) {
      fprintf(stderr, "unknown TS1 CSV column: %s\n", name);
      free(state); free(parameters); free(photo_base); free(temperature); free(pressure); free(photo_indices);
      csv_destroy(&cells_table); csv_destroy(&photo_table); return 2;
    }
    for (size_t cell = 0; cell < cells; ++cell) {
      const double value = cells_table.rows[cell][column];
      if (is_temperature) temperature[cell] = value;
      else if (is_pressure) pressure[cell] = value;
      else if (species >= 0) state[cell * n + (size_t)species] = value;
      else parameters[cell * pcount + (size_t)parameter] = value;
    }
  }
  memcpy(photo_base, parameters, cells * pcount * sizeof(double));
  csv_destroy(&cells_table);

  const size_t groups = (cells + group_size - 1u) / group_size;
  chem_statistics *group_stats = calloc(groups, sizeof(*group_stats));
  int *group_status = calloc(groups, sizeof(*group_status));
  const size_t snapshots_capacity = snapshot_every ? (steps / snapshot_every + (steps % snapshot_every != 0u)) : 1u;
  if (snapshots_capacity > SIZE_MAX / cells / n / sizeof(double)) {
    fputs("snapshot set is too large\n", stderr);
    free(state); free(parameters); free(photo_base); free(temperature); free(pressure); free(photo_indices);
    free(group_stats); free(group_status); csv_destroy(&photo_table); return 2;
  }
  double *snapshot_states = malloc(snapshots_capacity * cells * n * sizeof(double));
  size_t *snapshot_steps = malloc(snapshots_capacity * sizeof(size_t));
  if (!group_stats || !group_status || !snapshot_states || !snapshot_steps) {
    fputs("not enough memory for scenario groups and snapshots\n", stderr);
    free(state); free(parameters); free(photo_base); free(temperature); free(pressure); free(photo_indices);
    free(group_stats); free(group_status); free(snapshot_states); free(snapshot_steps); csv_destroy(&photo_table); return 2;
  }
  double elapsed = 0.0;
  size_t snapshot_count = 0;
  omp_set_dynamic(0);
  for (size_t step = 0; step < steps; ++step) {
    const double start = wall_time();
    const double factor = photolysis_factor(&photo_table, (double)step * interval);
    for (size_t cell = 0; cell < cells; ++cell)
      for (size_t p = 0; p < photo_count; ++p) {
        const size_t index = cell * pcount + photo_indices[p];
        parameters[index] = photo_base[index] * factor;
      }
    #pragma omp parallel for schedule(static) num_threads(threads)
    for (size_t group = 0; group < groups; ++group) {
      const size_t first = group * group_size;
      const size_t count = cells - first < group_size ? cells - first : group_size;
      double *group_parameters = pcount ? parameters + first * pcount : NULL;
      chem_statistics step_statistics = { 0 };
      if (chem_generated_integrate_group_steps(
              count, state + first * n, group_parameters, temperature + first, pressure + first,
              interval, 1u, relative_tolerance, absolute_tolerance, &step_statistics) != 0) {
        group_status[group] = -1;
      } else {
        group_stats[group].internal_steps += step_statistics.internal_steps;
        group_stats[group].accepted += step_statistics.accepted;
        group_stats[group].rejected += step_statistics.rejected;
        group_stats[group].function_calls += step_statistics.function_calls;
        group_stats[group].jacobian_updates += step_statistics.jacobian_updates;
        group_stats[group].decompositions += step_statistics.decompositions;
        group_stats[group].solves += step_statistics.solves;
      }
    }
    elapsed += wall_time() - start;
    for (size_t group = 0; group < groups; ++group)
      if (group_status[group]) {
        fprintf(stderr, "scenario integration failed at step %zu group %zu\n", step + 1u, group);
        free(state); free(parameters); free(photo_base); free(temperature); free(pressure); free(photo_indices);
        free(group_stats); free(group_status); free(snapshot_states); free(snapshot_steps); csv_destroy(&photo_table); return 1;
      }
    const size_t completed = step + 1u;
    if ((snapshot_every && completed % snapshot_every == 0u) || completed == steps) {
      snapshot_steps[snapshot_count] = completed;
      memcpy(snapshot_states + snapshot_count * cells * n, state, cells * n * sizeof(double));
      ++snapshot_count;
    }
  }
  chem_statistics total = { 0 };
  for (size_t group = 0; group < groups; ++group) {
    total.internal_steps += group_stats[group].internal_steps;
    total.accepted += group_stats[group].accepted;
    total.rejected += group_stats[group].rejected;
    total.function_calls += group_stats[group].function_calls;
    total.jacobian_updates += group_stats[group].jacobian_updates;
    total.decompositions += group_stats[group].decompositions;
    total.solves += group_stats[group].solves;
  }
  printf("elapsed_ms=%.3f\n", elapsed * 1000.0);
  printf("internal_steps=%lu accepted=%lu rejected=%lu function_calls=%lu jacobian_updates=%lu decompositions=%lu solves=%lu\n",
         total.internal_steps, total.accepted, total.rejected, total.function_calls, total.jacobian_updates,
         total.decompositions, total.solves);
  const int final_status = write_final(mechanism, cells, state, temperature, pressure);
  const int snapshot_status = write_snapshots(
      mechanism, cells, snapshot_count, snapshot_steps, snapshot_states, temperature, pressure);
  free(state); free(parameters); free(photo_base); free(temperature); free(pressure); free(photo_indices);
  free(group_stats); free(group_status); free(snapshot_states); free(snapshot_steps); csv_destroy(&photo_table);
  return final_status == 0 && snapshot_status == 0 ? 0 : 1;
}
#endif
