#include "kinetics.h"

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define LINE_CAP 262144
#define PATH_CAP 8192

typedef struct {
  char **names;
  size_t count;
} csv_header;

typedef struct {
  double *times;
  double *factors;
  size_t count;
} photo_series;

typedef struct {
  char cells[PATH_CAP];
  char photolysis[PATH_CAP];
  double dt;
  size_t steps;
  size_t snapshot_every;
} scenario_config;

static char *trim(char *value)
{
  char *end;
  while (isspace((unsigned char)*value))
    ++value;
  end = value + strlen(value);
  while (end > value && isspace((unsigned char)end[-1]))
    --end;
  *end = '\0';
  return value;
}

static char *next_field(char **cursor)
{
  char *field, *end;
  if (!cursor || !*cursor)
    return NULL;
  field = *cursor;
  end = field;
  while (*end && *end != ',' && *end != '\r' && *end != '\n')
    ++end;
  if (*end == ',') {
    *end = '\0';
    *cursor = end + 1;
  } else {
    *end = '\0';
    *cursor = NULL;
  }
  return field;
}

static char *copy_string(const char *source)
{
  size_t length = strlen(source) + 1;
  char *copy = (char *)malloc(length);
  if (copy)
    memcpy(copy, source, length);
  return copy;
}

static void free_header(csv_header *header)
{
  size_t i;
  for (i = 0; i < header->count; ++i)
    free(header->names[i]);
  free(header->names);
  header->names = NULL;
  header->count = 0;
}

static int read_header(FILE *file, csv_header *header)
{
  char line[LINE_CAP];
  char *cursor, *field;
  size_t capacity = 0;
  if (!fgets(line, sizeof(line), file))
    return -1;
  header->names = NULL;
  header->count = 0;
  cursor = line;
  while ((field = next_field(&cursor)) != NULL) {
    if (header->count == capacity) {
      size_t next = capacity ? capacity * 2 : 16;
      char **names = (char **)realloc(header->names, next * sizeof(*names));
      if (!names) {
        free_header(header);
        return -1;
      }
      header->names = names;
      capacity = next;
    }
    header->names[header->count] = copy_string(trim(field));
    if (!header->names[header->count]) {
      free_header(header);
      return -1;
    }
    ++header->count;
  }
  return header->count ? 0 : -1;
}

static int find_species(const chem_mechanism *mechanism, const char *name)
{
  size_t i;
  for (i = 0; i < mechanism->species_count; ++i)
    if (strcmp(mechanism->species[i].name, name) == 0)
      return (int)i;
  return -1;
}

static int find_custom(const chem_mechanism *mechanism, const char *name)
{
  size_t i;
  for (i = 0; i < mechanism->custom_count; ++i)
    if (strcmp(mechanism->custom_parameters[i].name, name) == 0)
      return (int)i;
  return -1;
}

static int column_index(const chem_mechanism *mechanism, const char *name)
{
  int index;
  if (strcmp(name, "temperature") == 0)
    return -2;
  if (strcmp(name, "pressure") == 0)
    return -3;
  index = find_species(mechanism, name);
  if (index >= 0)
    return index;
  index = find_custom(mechanism, name);
  return index >= 0 ? -1000 - index : -1;
}

static int count_cells(const char *path, size_t *count)
{
  FILE *file = fopen(path, "r");
  char line[LINE_CAP];
  size_t rows = 0;
  if (!file) {
    fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
    return -1;
  }
  if (!fgets(line, sizeof(line), file)) {
    fclose(file);
    return -1;
  }
  while (fgets(line, sizeof(line), file))
    if (*trim(line))
      ++rows;
  fclose(file);
  if (!rows)
    return -1;
  *count = rows;
  return 0;
}

static int read_cells(const chem_mechanism *mechanism, const char *path,
                      size_t cells, double *temperature, double *pressure,
                      double *custom, double *concentrations)
{
  FILE *file = fopen(path, "r");
  csv_header header = {0};
  int *columns = NULL;
  char line[LINE_CAP];
  size_t cell, i;
  int status = -1;
  if (!file) {
    fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
    return -1;
  }
  if (read_header(file, &header) != 0)
    goto done;
  columns = (int *)malloc(header.count * sizeof(*columns));
  if (!columns)
    goto done;
  for (i = 0; i < header.count; ++i) {
    columns[i] = column_index(mechanism, header.names[i]);
    if (columns[i] == -1) {
      fprintf(stderr, "%s: unknown species or rate parameter '%s'\n",
              path, header.names[i]);
      goto done;
    }
  }
  for (cell = 0; cell < cells; ++cell) {
    char *cursor, *field;
    size_t column = 0;
    if (!fgets(line, sizeof(line), file)) {
      fprintf(stderr, "%s: only %zu rows for %zu cells\n", path, cell, cells);
      goto done;
    }
    cursor = line;
    while ((field = next_field(&cursor)) != NULL) {
      char *start = trim(field), *end;
      double value;
      if (column >= header.count) {
        fprintf(stderr, "%s: row %zu has extra fields\n", path, cell + 1);
        goto done;
      }
      errno = 0;
      value = strtod(start, &end);
      if (errno || end == start || *end || !isfinite(value)) {
        fprintf(stderr, "%s: invalid value at row %zu, column %zu\n",
                path, cell + 1, column + 1);
        goto done;
      }
      if (columns[column] == -2)
        temperature[cell] = value;
      else if (columns[column] == -3)
        pressure[cell] = value;
      else if (columns[column] <= -1000)
        custom[cell * mechanism->custom_count + (size_t)(-1000 - columns[column])] = value;
      else
        concentrations[cell * mechanism->species_count + (size_t)columns[column]] = value;
      ++column;
    }
    if (column != header.count) {
      fprintf(stderr, "%s: row %zu has %zu fields; expected %zu\n",
              path, cell + 1, column, header.count);
      goto done;
    }
  }
  status = 0;
done:
  free(columns);
  free_header(&header);
  fclose(file);
  return status;
}

static int parse_positive_size(const char *text, size_t *value)
{
  char *end;
  unsigned long long parsed;
  errno = 0;
  parsed = strtoull(text, &end, 10);
  if (errno || end == text || *end || parsed == 0 || parsed > SIZE_MAX)
    return -1;
  *value = (size_t)parsed;
  return 0;
}

static int parse_positive_double(const char *text, double *value)
{
  char *end;
  errno = 0;
  *value = strtod(text, &end);
  return errno || end == text || *end || !isfinite(*value) || !(*value > 0.0)
             ? -1 : 0;
}

static int load_scenario(const char *path, scenario_config *config)
{
  FILE *file = fopen(path, "r");
  char line[LINE_CAP];
  memset(config, 0, sizeof(*config));
  if (!file) {
    fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
    return -1;
  }
  while (fgets(line, sizeof(line), file)) {
    char *comment = strchr(line, '#');
    char *equals, *key, *value;
    if (comment)
      *comment = '\0';
    key = trim(line);
    if (!*key)
      continue;
    equals = strchr(key, '=');
    if (!equals) {
      fclose(file);
      return -1;
    }
    *equals = '\0';
    key = trim(key);
    value = trim(equals + 1);
    if (strcmp(key, "cells") == 0)
      snprintf(config->cells, sizeof(config->cells), "%s", value);
    else if (strcmp(key, "photolysis") == 0)
      snprintf(config->photolysis, sizeof(config->photolysis), "%s", value);
    else if (strcmp(key, "dt") == 0) {
      if (parse_positive_double(value, &config->dt) != 0) {
        fclose(file);
        return -1;
      }
    } else if (strcmp(key, "steps") == 0) {
      if (parse_positive_size(value, &config->steps) != 0) {
        fclose(file);
        return -1;
      }
    } else if (strcmp(key, "snapshot_every") == 0) {
      if (parse_positive_size(value, &config->snapshot_every) != 0) {
        fclose(file);
        return -1;
      }
    }
  }
  fclose(file);
  if (!config->cells[0] || !config->photolysis[0] ||
      !(config->dt > 0.0) || !config->steps)
    return -1;
  if (!config->snapshot_every)
    config->snapshot_every = config->steps;
  return 0;
}

static int resolve_path(const char *scenario, const char *relative,
                        char *result, size_t capacity)
{
  const char *slash;
  size_t directory_length;
  if (relative[0] == '/') {
    if (snprintf(result, capacity, "%s", relative) >= (int)capacity)
      return -1;
    return 0;
  }
  slash = strrchr(scenario, '/');
  directory_length = slash ? (size_t)(slash - scenario + 1) : 0;
  if (directory_length + strlen(relative) + 1 > capacity)
    return -1;
  if (directory_length)
    memcpy(result, scenario, directory_length);
  memcpy(result + directory_length, relative, strlen(relative) + 1);
  return 0;
}

static int read_photo_series(const char *path, photo_series *series)
{
  FILE *file = fopen(path, "r");
  char line[LINE_CAP];
  size_t capacity = 0;
  memset(series, 0, sizeof(*series));
  if (!file) {
    fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
    return -1;
  }
  if (!fgets(line, sizeof(line), file)) {
    fclose(file);
    return -1;
  }
  while (fgets(line, sizeof(line), file)) {
    char *cursor = line, *first = next_field(&cursor), *second = next_field(&cursor);
    char *end1, *end2;
    double time, factor;
    if (!first || !second)
      continue;
    errno = 0;
    time = strtod(trim(first), &end1);
    factor = strtod(trim(second), &end2);
    if (errno || *end1 || *end2 || !isfinite(time) || !isfinite(factor) ||
        (series->count && time <= series->times[series->count - 1])) {
      fclose(file);
      free(series->times);
      free(series->factors);
      memset(series, 0, sizeof(*series));
      return -1;
    }
    if (series->count == capacity) {
      size_t next = capacity ? capacity * 2 : 64;
      double *times = (double *)realloc(series->times, next * sizeof(*times));
      double *factors;
      if (!times) {
        fclose(file);
        free(series->times);
        free(series->factors);
        memset(series, 0, sizeof(*series));
        return -1;
      }
      series->times = times;
      factors = (double *)realloc(series->factors, next * sizeof(*factors));
      if (!factors) {
        fclose(file);
        free(series->times);
        free(series->factors);
        memset(series, 0, sizeof(*series));
        return -1;
      }
      series->factors = factors;
      capacity = next;
    }
    series->times[series->count] = time;
    series->factors[series->count] = factor;
    ++series->count;
  }
  fclose(file);
  return series->count ? 0 : -1;
}

static double photo_factor(const photo_series *series, double time)
{
  size_t low, high;
  if (time <= series->times[0])
    return series->factors[0];
  if (time >= series->times[series->count - 1])
    return series->factors[series->count - 1];
  low = 0;
  high = series->count - 1;
  while (high - low > 1) {
    size_t middle = low + (high - low) / 2;
    if (series->times[middle] <= time)
      low = middle;
    else
      high = middle;
  }
  return series->factors[low] +
         (time - series->times[low]) /
             (series->times[high] - series->times[low]) *
             (series->factors[high] - series->factors[low]);
}

static double wall_time(void)
{
#ifdef _OPENMP
  return omp_get_wtime();
#else
  struct timespec now;
  timespec_get(&now, TIME_UTC);
  return (double)now.tv_sec + (double)now.tv_nsec * 1.0e-9;
#endif
}

static int write_final(size_t cells, const double *temperature,
                       const double *pressure, const double *state)
{
  FILE *file = fopen("ts1_final.csv", "w");
  size_t cell, species;
  if (!file)
    return -1;
  fputs("temperature,pressure", file);
  for (species = 0; species < chem_model.species_count; ++species)
    fprintf(file, ",%s", chem_model.species[species].name);
  fputc('\n', file);
  for (cell = 0; cell < cells; ++cell) {
    fprintf(file, "%.17g,%.17g", temperature[cell], pressure[cell]);
    for (species = 0; species < chem_model.species_count; ++species)
      fprintf(file, ",%.17g", state[cell * chem_model.species_count + species]);
    fputc('\n', file);
  }
  return fclose(file) == 0 ? 0 : -1;
}

static int write_snapshot_header(FILE *file)
{
  size_t species;
  fputs("step,cell,temperature,pressure", file);
  for (species = 0; species < chem_model.species_count; ++species)
    fprintf(file, ",%s", chem_model.species[species].name);
  return fputc('\n', file) == EOF ? -1 : 0;
}

static int write_snapshot(FILE *file, size_t step, size_t cells,
                          const double *temperature, const double *pressure,
                          const double *state)
{
  size_t cell, species;
  for (cell = 0; cell < cells; ++cell) {
    fprintf(file, "%zu,%zu,%.17g,%.17g", step, cell,
            temperature[cell], pressure[cell]);
    for (species = 0; species < chem_model.species_count; ++species)
      fprintf(file, ",%.17g", state[cell * chem_model.species_count + species]);
    if (fputc('\n', file) == EOF)
      return -1;
  }
  return ferror(file) ? -1 : 0;
}

int main(int argc, char **argv)
{
  size_t threads = 0, group_size = 0, cells = 0, arg = 1, cell, i, step;
  double relative_tolerance = 1.0e-6, absolute_tolerance = 1.0e-3;
  char scenario_path[PATH_CAP], cells_path[PATH_CAP], photo_path[PATH_CAP];
  const char *scenario_argument = NULL;
  scenario_config config;
  photo_series series = {0};
  double *temperature = NULL, *pressure = NULL, *custom = NULL, *custom_base = NULL;
  double *state = NULL;
  unsigned char *photolysis_parameter = NULL;
  FILE *snapshots = NULL;
  chem_solver *solver = NULL;
  chem_rosenbrock_parameters parameters;
  chem_solver_stats total = {0};
  double elapsed = 0.0;
  int status = 1;

  while (arg < (size_t)argc && argv[arg][0] == '-') {
    if (strcmp(argv[arg], "--threads") == 0 ||
        strcmp(argv[arg], "--group-size") == 0) {
      size_t value;
      int is_threads = strcmp(argv[arg], "--threads") == 0;
      if (arg + 1 >= (size_t)argc || parse_positive_size(argv[arg + 1], &value) != 0)
        goto usage;
      if (is_threads)
        threads = value;
      else
        group_size = value;
      arg += 2;
    } else if (strcmp(argv[arg], "--rtol") == 0 || strcmp(argv[arg], "--atol") == 0) {
      double value;
      int is_rtol = strcmp(argv[arg], "--rtol") == 0;
      if (arg + 1 >= (size_t)argc || parse_positive_double(argv[arg + 1], &value) != 0)
        goto usage;
      if (is_rtol)
        relative_tolerance = value;
      else
        absolute_tolerance = value;
      arg += 2;
    } else {
      goto usage;
    }
  }
  if (arg + 1 != (size_t)argc)
    goto usage;
  scenario_argument = argv[arg];
  if (snprintf(scenario_path, sizeof(scenario_path), "%s", scenario_argument) >=
          (int)sizeof(scenario_path) ||
      load_scenario(scenario_path, &config) != 0 ||
      resolve_path(scenario_argument, config.cells, cells_path, sizeof(cells_path)) != 0 ||
      resolve_path(scenario_argument, config.photolysis, photo_path, sizeof(photo_path)) != 0 ||
      count_cells(cells_path, &cells) != 0 || read_photo_series(photo_path, &series) != 0) {
    fprintf(stderr, "invalid scenario or missing input\n");
    goto done;
  }
#ifdef _OPENMP
  if (threads)
    omp_set_num_threads((int)threads);
#else
  if (threads > 1) {
    fprintf(stderr, "this build has no OpenMP support\n");
    goto done;
  }
#endif
  temperature = (double *)malloc(cells * sizeof(*temperature));
  pressure = (double *)malloc(cells * sizeof(*pressure));
  custom = (double *)malloc((cells && chem_model.custom_count
                                  ? cells * chem_model.custom_count : 1) * sizeof(*custom));
  custom_base = (double *)malloc((cells && chem_model.custom_count
                                       ? cells * chem_model.custom_count : 1) * sizeof(*custom_base));
  photolysis_parameter = (unsigned char *)calloc(chem_model.custom_count
                                                      ? chem_model.custom_count : 1, 1);
  state = (double *)malloc(cells * chem_model.species_count * sizeof(*state));
  if (!temperature || !pressure || !custom || !custom_base ||
      !photolysis_parameter || !state) {
    fprintf(stderr, "cannot allocate scenario state for %zu cells\n", cells);
    goto done;
  }
  for (cell = 0; cell < cells; ++cell) {
    temperature[cell] = chem_model.default_temperature;
    pressure[cell] = chem_model.default_pressure;
    for (i = 0; i < chem_model.custom_count; ++i)
      custom[cell * chem_model.custom_count + i] =
          chem_model.custom_parameters[i].default_value;
    for (i = 0; i < chem_model.species_count; ++i)
      state[cell * chem_model.species_count + i] = chem_model.species[i].initial_value;
  }
  if (read_cells(&chem_model, cells_path, cells, temperature, pressure,
                 custom, state) != 0)
    goto done;
  memcpy(custom_base, custom, cells * chem_model.custom_count * sizeof(*custom));
  for (i = 0; i < chem_model.custom_count; ++i)
    photolysis_parameter[i] = chem_model.custom_parameters[i].name[0] == 'j';
  chem_rosenbrock_default(&parameters);
  parameters.relative_tolerance = relative_tolerance;
  parameters.default_absolute_tolerance = absolute_tolerance;
  solver = chem_solver_create(&chem_model, cells, temperature, pressure, custom,
                              &parameters, group_size);
  if (!solver) {
    fprintf(stderr, "could not create specialized scenario solver\n");
    goto done;
  }
  snapshots = fopen("ts1_snapshots.csv", "w");
  if (!snapshots || write_snapshot_header(snapshots) != 0) {
    fprintf(stderr, "cannot write ts1_snapshots.csv\n");
    goto done;
  }
  for (step = 0; step < config.steps; ++step) {
    const double start = wall_time();
    const double factor = photo_factor(&series, (double)step * config.dt);
    chem_solver_stats current;
    for (cell = 0; cell < cells; ++cell)
      for (i = 0; i < chem_model.custom_count; ++i)
        if (photolysis_parameter[i])
          custom[cell * chem_model.custom_count + i] =
              custom_base[cell * chem_model.custom_count + i] * factor;
    if (chem_solver_update_custom_rates(solver, temperature, pressure, custom) != 0 ||
        chem_solver_step(solver, state, config.dt, &current) != 0) {
      fprintf(stderr, "scenario integration failed at step %zu\n", step + 1);
      goto done;
    }
    elapsed += wall_time() - start;
    total.internal_steps += current.internal_steps;
    total.accepted += current.accepted;
    total.rejected += current.rejected;
    total.function_calls += current.function_calls;
    total.jacobian_updates += current.jacobian_updates;
    total.decompositions += current.decompositions;
    total.solves += current.solves;
    if ((step + 1) % config.snapshot_every == 0 || step + 1 == config.steps)
      if (write_snapshot(snapshots, step + 1, cells, temperature, pressure, state) != 0) {
        fprintf(stderr, "cannot write scenario snapshot at step %zu\n", step + 1);
        goto done;
      }
  }
  if (fclose(snapshots) != 0) {
    snapshots = NULL;
    fprintf(stderr, "cannot finish ts1_snapshots.csv\n");
    goto done;
  }
  snapshots = NULL;
  if (write_final(cells, temperature, pressure, state) != 0) {
    fprintf(stderr, "cannot write ts1_final.csv\n");
    goto done;
  }
  printf("elapsed_ms=%.3f\n", elapsed * 1000.0);
  printf("internal_steps=%" PRIu64 " accepted=%" PRIu64 " rejected=%" PRIu64
         " function_calls=%" PRIu64 " jacobian_updates=%" PRIu64
         " decompositions=%" PRIu64 " solves=%" PRIu64 "\n",
         total.internal_steps, total.accepted, total.rejected, total.function_calls,
         total.jacobian_updates, total.decompositions, total.solves);
  status = 0;
  goto done;

usage:
  fprintf(stderr, "usage: %s [--threads N] [--group-size G] [--rtol R] [--atol A] SCENARIO\n",
          argv[0]);
done:
  if (snapshots)
    fclose(snapshots);
  chem_solver_destroy(solver);
  free(series.times);
  free(series.factors);
  free(temperature);
  free(pressure);
  free(custom);
  free(custom_base);
  free(photolysis_parameter);
  free(state);
  return status;
}
