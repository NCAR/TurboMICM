#include "chem.h"

#include <errno.h>
#include <omp.h>
#include <stdio.h>
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
    if (chem_integrate_group(mechanism, count, group_y, group_parameters, temperature + first, pressure + first,
                             interval, NULL) != 0)
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
      if (chem_integrate_group_steps(mechanism, count, group_y, group_parameters, temperature + first,
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
