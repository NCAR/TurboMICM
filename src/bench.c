#include "chem.h"
#include <errno.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

ChemModel *chem_load_ts1(int specialized);
#ifndef CHEM_GROUP_SIZE
#define CHEM_GROUP_SIZE 32
#endif

typedef struct { char **names; double *values; size_t columns, rows; } Table;
static char *line_read(FILE *file) {
  size_t capacity = 1024, length = 0; char *line = malloc(capacity);
  if (!line) return NULL;
  int ch;
  while ((ch = fgetc(file)) != EOF && ch != '\n') {
    if (length+1 >= capacity) {
      if (capacity > SIZE_MAX/2) { free(line); return NULL; }
      capacity *= 2; char *p = realloc(line, capacity);
      if (!p) { free(line); return NULL; } line = p;
    }
    line[length++] = (char)ch;
  }
  if (!length && ch == EOF) { free(line); return NULL; }
  if (length && line[length-1] == '\r') --length;
  line[length] = '\0'; return line;
}
static void table_free(Table *table) {
  for (size_t c = 0; c < table->columns; ++c) free(table->names[c]);
  free(table->names); free(table->values);
}
static int table_read(const char *path, Table *t, size_t needed) {
  FILE *file = fopen(path, "r");
  if (!file) { perror(path); return 0; }
  char *header = line_read(file);
  if (!header || !*header) { free(header); fclose(file); return 0; }
  t->columns = 1;
  for (char *p = header; *p; ++p) if (*p == ',') ++t->columns;
  t->names = calloc(t->columns, sizeof(char *));
  if (!t->names || needed > SIZE_MAX/t->columns/sizeof(double)) { free(header); fclose(file); return 0; }
  t->values = malloc(needed*t->columns*sizeof(double));
  if (!t->values) { free(header); fclose(file); return 0; }
  char *field = header;
  for (size_t c = 0; c < t->columns; ++c) {
    char *end = strchr(field, ','); if (end) *end = '\0';
    t->names[c] = malloc(strlen(field)+1);
    if (!t->names[c] || !*field) { free(header); fclose(file); return 0; }
    strcpy(t->names[c], field); field = end ? end+1 : field+strlen(field);
  }
  free(header);
  while (t->rows < needed) {
    char *line = line_read(file);
    if (!line) break;
    if (!*line) { free(line); continue; }
    field = line;
    for (size_t c = 0; c < t->columns; ++c) {
      char *end; errno = 0; double val = strtod(field, &end);
      if (end == field || errno || !isfinite(val) || (c+1 == t->columns ? *end != '\0' : *end != ',')) {
        fprintf(stderr, "invalid CSV row %zu, column %zu\n", t->rows+1, c+1); free(line); fclose(file); return 0;
      }
      t->values[t->rows*t->columns+c] = val; field = *end ? end+1 : end;
    }
    ++t->rows; free(line);
  }
  fclose(file);
  if (t->rows != needed) { fprintf(stderr, "CSV has %zu cells, need %zu\n", t->rows, needed); return 0; }
  return 1;
}
static int positive_size(const char *text, size_t *value, int allow_zero) {
  char *end; errno = 0;
  if (!*text || *text == '-') return 0;
  unsigned long long v = strtoull(text, &end, 10);
  if (*end || errno || v > SIZE_MAX || (!allow_zero && !v)) return 0;
  *value = (size_t)v; return 1;
}
static void add_stats(ChemStats *sum, const ChemStats *v) {
  sum->internal_steps += v->internal_steps; sum->accepted += v->accepted; sum->rejected += v->rejected;
  sum->function_calls += v->function_calls; sum->jacobian_updates += v->jacobian_updates;
  sum->decompositions += v->decompositions; sum->solves += v->solves;
}
int main(int argc, char **argv) {
  size_t threads = (size_t)omp_get_max_threads(), group_size = CHEM_GROUP_SIZE, cells, steps;
  int arg = 1;
  while (arg+1 < argc && !strncmp(argv[arg], "--", 2)) {
    size_t v; if (!positive_size(argv[arg+1], &v, 0)) goto usage;
    if (!strcmp(argv[arg], "--threads")) threads = v;
    else if (!strcmp(argv[arg], "--group")) group_size = v;
    else goto usage;
    arg += 2;
  }
  if (argc-arg != 4 || !positive_size(argv[arg], &cells, 0) || !positive_size(argv[arg+1], &steps, 1) ||
      threads > 2147483647 || !group_size) goto usage;
  char *end; errno = 0; double dt = strtod(argv[arg+2], &end);
  if (errno || *end || !isfinite(dt) || dt <= 0.0) goto usage;
  omp_set_num_threads((int)threads);
  Table table = {0};
  if (!table_read(argv[arg+3], &table, cells)) { table_free(&table); return 1; }
  ChemModel *m = chem_load_ts1(1);
  if (!m) { fprintf(stderr, "cannot build mechanism\n"); table_free(&table); return 1; }
  for (size_t c = 0; c < table.columns; ++c) {
    const char *label = table.names[c];
    if (strcmp(label, "temperature") && strcmp(label, "pressure") && chem_variable_id(m, label) == SIZE_MAX && chem_custom_id(m, label) == SIZE_MAX) {
      fprintf(stderr, "unknown CSV column: %s\n", label); table_free(&table); chem_model_free(m); return 1;
    }
  }
  size_t groups = 1+(cells-1)/group_size;
  ChemState **states = calloc(groups, sizeof *states);
  ChemStats *statistics = calloc(groups, sizeof *statistics);
  ChemStatus *statuses = calloc(groups, sizeof *statuses);
  if (!states || !statistics || !statuses) { free(states); free(statistics); free(statuses); table_free(&table); chem_model_free(m); return 1; }
  int failed = 0;
  #pragma omp parallel for schedule(static)
  for (size_t g = 0; g < groups; ++g) {
    size_t first = g*group_size, count = cells-first < group_size ? cells-first : group_size;
    ChemState *s = states[g] = chem_state_create(m, count);
    if (!s) { statuses[g] = CHEM_MEMORY; continue; }
    for (size_t c = 0; c < table.columns; ++c) {
      const char *label = table.names[c];
      size_t v = chem_variable_id(m, label), u = chem_custom_id(m, label);
      for (size_t i = 0; i < count; ++i) {
        double value = table.values[(first+i)*table.columns+c];
        if (!strcmp(label, "temperature")) s->conditions[i].temperature = value;
        else if (!strcmp(label, "pressure")) s->conditions[i].pressure = value;
        else if (v != SIZE_MAX) s->variables[v*s->stride+i] = value;
        else s->custom[u*s->stride+i] = value;
      }
    }
    for (size_t c = 0; c < count; ++c) chem_ideal_air_density(&s->conditions[c]);
    chem_calculate_rates(s);
    ChemSolver solver = chem_solver(3); ChemStats warmup;
    statuses[g] = chem_solve(&solver, dt, s, &warmup);
  }
  for (size_t g = 0; g < groups; ++g) if (statuses[g]) { fprintf(stderr, "warm-up group %zu: %s\n", g, chem_status_string(statuses[g])); failed = 1; }
  double start = omp_get_wtime();
  if (!failed) {
    #pragma omp parallel for schedule(static)
    for (size_t g = 0; g < groups; ++g) {
      ChemSolver solver = chem_solver(3);
      for (size_t step = 0; step < steps; ++step) {
        ChemStats stats;
        statuses[g] = chem_solve(&solver, dt, states[g], &stats);
        add_stats(&statistics[g], &stats);
        if (statuses[g]) break;
      }
    }
  }
  double elapsed = (omp_get_wtime()-start)*1000.0;
  ChemStats sum = {0};
  for (size_t g = 0; g < groups; ++g) {
    add_stats(&sum, &statistics[g]);
    if (statuses[g]) { fprintf(stderr, "group %zu: %s\n", g, chem_status_string(statuses[g])); failed = 1; }
  }
  if (!failed) {
    printf("elapsed_ms=%.3f\n", elapsed);
    printf("internal_steps=%llu accepted=%llu rejected=%llu function_calls=%llu jacobian_updates=%llu decompositions=%llu solves=%llu\n",
      (unsigned long long)sum.internal_steps, (unsigned long long)sum.accepted, (unsigned long long)sum.rejected,
      (unsigned long long)sum.function_calls, (unsigned long long)sum.jacobian_updates,
      (unsigned long long)sum.decompositions, (unsigned long long)sum.solves);
    FILE *out = fopen("ts1_final.csv", "w");
    if (!out) { perror("ts1_final.csv"); failed = 1; }
    else {
      fprintf(out, "temperature,pressure");
      for (size_t i = 0; i < m->nvariables; ++i) fprintf(out, ",%s", m->species[m->variable_species[i]].name);
      fputc('\n', out);
      for (size_t c = 0; c < cells; ++c) {
        ChemState *s = states[c/group_size]; size_t i = c%group_size;
        fprintf(out, "%.17g,%.17g", s->conditions[i].temperature, s->conditions[i].pressure);
        for (size_t v = 0; v < m->nvariables; ++v) fprintf(out, ",%.17g", s->variables[v*s->stride+i]);
        fputc('\n', out);
      }
      if (ferror(out)) failed = 1;
      if (fclose(out)) failed = 1;
    }
  }
  for (size_t g = 0; g < groups; ++g) chem_state_free(states[g]);
  free(states); free(statistics); free(statuses); table_free(&table); chem_model_free(m);
  return failed ? 1 : 0;
usage:
  fprintf(stderr, "usage: %s [--threads N] [--group N] CELLS STEPS DT CELLS_CSV\n", argv[0]);
  return 2;
}
