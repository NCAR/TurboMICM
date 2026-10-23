#include "kinetics.h"

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define LINE_CAP 262144

typedef struct {
  char **names;
  size_t count;
} csv_header;

static char *trim(char *s)
{
  char *end;
  while (isspace((unsigned char)*s))
    ++s;
  end = s + strlen(s);
  while (end > s && isspace((unsigned char)end[-1]))
    --end;
  *end = '\0';
  return s;
}

static char *next_csv_field(char **cursor)
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
  const size_t length = strlen(source) + 1;
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
  while ((field = next_csv_field(&cursor)) != NULL) {
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

static int find_species(const chem_mechanism *m, const char *name)
{
  size_t i;
  for (i = 0; i < m->species_count; ++i)
    if (strcmp(m->species[i].name, name) == 0)
      return (int)i;
  return -1;
}

static int find_custom(const chem_mechanism *m, const char *name)
{
  size_t i;
  for (i = 0; i < m->custom_count; ++i)
    if (strcmp(m->custom_parameters[i].name, name) == 0)
      return (int)i;
  return -1;
}

static int parse_column(const chem_mechanism *m, const char *name)
{
  int index;
  if (strcmp(name, "temperature") == 0)
    return -2;
  if (strcmp(name, "pressure") == 0)
    return -3;
  index = find_species(m, name);
  if (index >= 0)
    return index;
  index = find_custom(m, name);
  if (index >= 0)
    return -1000 - index;
  return -1;
}

static int read_cells(const chem_mechanism *m, const char *path, size_t count,
                      double *temperature, double *pressure, double *custom,
                      double *y)
{
  FILE *file = fopen(path, "r");
  csv_header header = {0};
  int *columns = NULL;
  char line[LINE_CAP];
  size_t i, cell;
  int status = -1;
  if (!file) {
    fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
    return -1;
  }
  if (read_header(file, &header) != 0) {
    fprintf(stderr, "%s: missing CSV header\n", path);
    goto done;
  }
  columns = (int *)malloc(header.count * sizeof(*columns));
  if (!columns)
    goto done;
  for (i = 0; i < header.count; ++i) {
    columns[i] = parse_column(m, header.names[i]);
    if (columns[i] == -1) {
      fprintf(stderr, "%s: unknown TS1 species or rate parameter '%s'\n", path, header.names[i]);
      goto done;
    }
  }
  for (cell = 0; cell < count; ++cell) {
    char *cursor, *field;
    size_t column = 0;
    if (!fgets(line, sizeof(line), file)) {
      fprintf(stderr, "%s: only %zu data rows for %zu cells\n", path, cell, count);
      goto done;
    }
    cursor = line;
    while ((field = next_csv_field(&cursor)) != NULL) {
      char *end;
      char *start;
      double value;
      if (column >= header.count) {
        fprintf(stderr, "%s: row %zu has more fields than its header\n", path, cell + 1);
        goto done;
      }
      errno = 0;
      start = trim(field);
      value = strtod(start, &end);
      if (errno || end == start || *end != '\0' || !isfinite(value)) {
        fprintf(stderr, "%s: invalid number at row %zu, column %zu\n", path, cell + 1, column + 1);
        goto done;
      }
      if (columns[column] == -2)
        temperature[cell] = value;
      else if (columns[column] == -3)
        pressure[cell] = value;
      else if (columns[column] <= -1000)
        custom[cell * m->custom_count + (size_t)(-1000 - columns[column])] = value;
      else
        y[cell * m->species_count + (size_t)columns[column]] = value;
      ++column;
    }
    if (column != header.count) {
      fprintf(stderr, "%s: row %zu has %zu fields; expected %zu\n", path, cell + 1, column, header.count);
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

static int write_final(const chem_mechanism *m, size_t cells,
                       const double *temperature, const double *pressure,
                       const double *y)
{
  FILE *file = fopen("ts1_final.csv", "w");
  size_t cell, variable;
  if (!file) {
    fprintf(stderr, "cannot write ts1_final.csv: %s\n", strerror(errno));
    return -1;
  }
  fputs("temperature,pressure", file);
  for (variable = 0; variable < m->species_count; ++variable)
    fprintf(file, ",%s", m->species[variable].name);
  fputc('\n', file);
  for (cell = 0; cell < cells; ++cell) {
    fprintf(file, "%.17g,%.17g", temperature[cell], pressure[cell]);
    for (variable = 0; variable < m->species_count; ++variable)
      fprintf(file, ",%.17g", y[cell * m->species_count + variable]);
    fputc('\n', file);
  }
  return fclose(file) == 0 ? 0 : -1;
}

static int positive_size(const char *text, size_t *value)
{
  char *end;
  unsigned long long parsed;
  errno = 0;
  parsed = strtoull(text, &end, 10);
  if (errno || end == text || *end || parsed == 0)
    return -1;
  *value = (size_t)parsed;
  return 0;
}

int main(int argc, char **argv)
{
  size_t threads = 0, group_size = 0, cells, steps, cell;
  double dt, elapsed_start, elapsed_end;
  char *end;
  const char *csv_path;
  double *temperature = NULL, *pressure = NULL, *custom = NULL, *y = NULL;
  chem_solver *solver = NULL;
  chem_rosenbrock_parameters parameters;
  chem_solver_stats warmup, total = {0};
  int arg = 1;
  if (argc > 1 && strcmp(argv[1], "--threads") == 0) {
    if (argc <= 2 || positive_size(argv[2], &threads) != 0) {
      fprintf(stderr, "usage: %s [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV\n", argv[0]);
      return 2;
    }
    arg = 3;
  }
  if (arg < argc && strcmp(argv[arg], "--group-size") == 0) {
    if (arg + 1 >= argc || positive_size(argv[arg + 1], &group_size) != 0) {
      fprintf(stderr, "invalid --group-size\n");
      return 2;
    }
    arg += 2;
  }
  if (argc - arg != 4 || positive_size(argv[arg], &cells) != 0 ||
      positive_size(argv[arg + 1], &steps) != 0) {
    fprintf(stderr, "usage: %s [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV\n", argv[0]);
    return 2;
  }
  errno = 0;
  dt = strtod(argv[arg + 2], &end);
  if (errno || end == argv[arg + 2] || *end || !isfinite(dt) || dt <= 0.0) {
    fprintf(stderr, "DT must be a positive number\n");
    return 2;
  }
  csv_path = argv[arg + 3];
#ifdef _OPENMP
  if (threads)
    omp_set_num_threads((int)threads);
#else
  if (threads > 1) {
    fprintf(stderr, "this build has no OpenMP support\n");
    return 2;
  }
#endif
  temperature = (double *)malloc(cells * sizeof(*temperature));
  pressure = (double *)malloc(cells * sizeof(*pressure));
  custom = (double *)malloc((cells != 0 && chem_model.custom_count != 0
                                 ? cells * chem_model.custom_count : 1) * sizeof(*custom));
  y = (double *)malloc(cells * chem_model.species_count * sizeof(*y));
  if (!temperature || !pressure || !custom || !y) {
    fprintf(stderr, "cannot allocate benchmark state for %zu cells\n", cells);
    goto error;
  }
  for (cell = 0; cell < cells; ++cell) {
    size_t i;
    temperature[cell] = chem_model.default_temperature;
    pressure[cell] = chem_model.default_pressure;
    for (i = 0; i < chem_model.custom_count; ++i)
      custom[cell * chem_model.custom_count + i] = chem_model.custom_parameters[i].default_value;
    for (i = 0; i < chem_model.species_count; ++i)
      y[cell * chem_model.species_count + i] = chem_model.species[i].initial_value;
  }
  if (read_cells(&chem_model, csv_path, cells, temperature, pressure, custom, y) != 0)
    goto error;
  chem_rosenbrock_default(&parameters);
  solver = chem_solver_create(&chem_model, cells, temperature, pressure, custom,
                              &parameters, group_size);
  if (!solver) {
    fprintf(stderr, "could not build solver state\n");
    goto error;
  }
  if (chem_solver_step(solver, y, dt, &warmup) != 0) {
    fprintf(stderr, "warm-up step failed\n");
    goto error;
  }
  elapsed_start = wall_time();
  for (cell = 0; cell < steps; ++cell) {
    chem_solver_stats current;
    if (chem_solver_step(solver, y, dt, &current) != 0) {
      fprintf(stderr, "timed step %zu failed\n", cell + 1);
      goto error;
    }
    total.internal_steps += current.internal_steps;
    total.accepted += current.accepted;
    total.rejected += current.rejected;
    total.function_calls += current.function_calls;
    total.jacobian_updates += current.jacobian_updates;
    total.decompositions += current.decompositions;
    total.solves += current.solves;
  }
  elapsed_end = wall_time();
  if (write_final(&chem_model, cells, temperature, pressure, y) != 0)
    goto error;
  printf("elapsed_ms=%.3f\n", (elapsed_end - elapsed_start) * 1000.0);
  printf("internal_steps=%" PRIu64 " accepted=%" PRIu64 " rejected=%" PRIu64
         " function_calls=%" PRIu64 " jacobian_updates=%" PRIu64
         " decompositions=%" PRIu64 " solves=%" PRIu64 "\n",
         total.internal_steps, total.accepted, total.rejected, total.function_calls,
         total.jacobian_updates, total.decompositions, total.solves);
  chem_solver_destroy(solver);
  free(temperature);
  free(pressure);
  free(custom);
  free(y);
  return 0;
error:
  chem_solver_destroy(solver);
  free(temperature);
  free(pressure);
  free(custom);
  free(y);
  return 1;
}
