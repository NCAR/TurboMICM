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

int main(int argc, char **argv) {
  size_t threads = 1;
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
    } else if (npos < 4) {
      positional[npos++] = argv[i];
    } else {
      fprintf(stderr, "usage: %s [--threads N] CELLS STEPS DT CELLS_CSV\n", argv[0]);
      return 2;
    }
  }
  if (npos != 4) {
    fprintf(stderr, "usage: %s [--threads N] CELLS STEPS DT CELLS_CSV\n", argv[0]);
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
  chem_solver_stats *stats = calloc(cells, sizeof(*stats));
  int *status = calloc(cells, sizeof(*status));
  if (!states || !rates || !temperature || !pressure || !stats || !status) {
    fprintf(stderr, "allocation failed for %zu cells\n", cells);
    free(states); free(rates); free(temperature); free(pressure); free(stats); free(status);
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
    }
  }
  free(columns);
  free_table(&input);
  if (threads > cells) threads = cells;
  chem_solver_options options = chem_solver_default_options();
  chem_solver_workspace *workspaces = calloc(threads, sizeof(*workspaces));
  if (!workspaces) {
    fprintf(stderr, "allocation failed for solver workspaces\n");
    free(states); free(rates); free(temperature); free(pressure); free(stats); free(status);
    return 1;
  }
  for (size_t i = 0; i < threads; ++i) {
    if (chem_solver_workspace_init(&workspaces[i], m)) {
      fprintf(stderr, "allocation failed for solver workspace %zu\n", i);
      destroy_workspaces(workspaces, threads);
      free(states); free(rates); free(temperature); free(pressure); free(stats); free(status);
      return 1;
    }
  }
  #pragma omp parallel num_threads(threads)
  {
#ifdef _OPENMP
    const size_t worker = (size_t)omp_get_thread_num();
#else
    const size_t worker = 0;
#endif
    #pragma omp for schedule(static)
    for (size_t c = 0; c < cells; ++c) {
      if (status[c]) continue;
      const double density = pressure[c] /
                             (8.31446261815324 * temperature[c]);
      chem_solver_stats warmup_stats;
      status[c] = chem_rosenbrock3_solve_with_workspace(
          m, &workspaces[worker], states+c*n, temperature[c], density,
          rates+c*nr, dt, &options, &warmup_stats);
    }
  }
  for (size_t c = 0; c < cells; ++c) {
    if (status[c]) {
      fprintf(stderr, "cell %zu failed warm-up solve with status %d\n", c, status[c]);
      destroy_workspaces(workspaces, threads);
      free(states); free(rates); free(temperature); free(pressure); free(stats); free(status);
      return 1;
    }
  }
  const double start = now_seconds();
  #pragma omp parallel num_threads(threads)
  {
#ifdef _OPENMP
    const size_t worker = (size_t)omp_get_thread_num();
#else
    const size_t worker = 0;
#endif
    #pragma omp for schedule(static)
    for (size_t c = 0; c < cells; ++c) {
      const double density = pressure[c] /
                             (8.31446261815324 * temperature[c]);
      for (size_t step = 0; step < steps && status[c] == 0; ++step) {
        chem_solver_stats one = {0};
        status[c] = chem_rosenbrock3_solve_with_workspace(
            m, &workspaces[worker], states+c*n, temperature[c], density,
            rates+c*nr, dt, &options, &one);
        stats[c].internal_steps += one.internal_steps;
        stats[c].accepted += one.accepted;
        stats[c].rejected += one.rejected;
        stats[c].function_calls += one.function_calls;
        stats[c].jacobian_updates += one.jacobian_updates;
        stats[c].decompositions += one.decompositions;
        stats[c].solves += one.solves;
      }
    }
  }
  const double elapsed_ms = (now_seconds() - start)*1000.0;
  chem_solver_stats total = {0};
  for (size_t c = 0; c < cells; ++c) {
    if (status[c]) {
      fprintf(stderr, "cell %zu failed integration with status %d\n", c, status[c]);
      destroy_workspaces(workspaces, threads);
      free(states); free(rates); free(temperature); free(pressure); free(stats); free(status);
      return 1;
    }
    total.internal_steps += stats[c].internal_steps;
    total.accepted += stats[c].accepted;
    total.rejected += stats[c].rejected;
    total.function_calls += stats[c].function_calls;
    total.jacobian_updates += stats[c].jacobian_updates;
    total.decompositions += stats[c].decompositions;
    total.solves += stats[c].solves;
  }
  FILE *output = fopen("ts1_final.csv", "w");
  if (!output) {
    fprintf(stderr, "cannot write ts1_final.csv: %s\n", strerror(errno));
    destroy_workspaces(workspaces, threads);
    free(states); free(rates); free(temperature); free(pressure); free(stats); free(status);
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
  free(states); free(rates); free(temperature); free(pressure); free(stats); free(status);
  destroy_workspaces(workspaces, threads);
  return 0;
}
