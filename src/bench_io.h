/* Shared C11 CSV, state initialization and output for the TS1 drivers. */
#ifndef CHEM_BENCH_IO_H
#define CHEM_BENCH_IO_H
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
  size_t capacity = needed ? needed : 64;
  if (capacity > SIZE_MAX/t->columns/sizeof(double)) { free(header); fclose(file); return 0; }
  t->values = malloc(capacity*t->columns*sizeof(double));
  if (!t->values) { free(header); fclose(file); return 0; }
  char *field = header;
  for (size_t c = 0; c < t->columns; ++c) {
    char *end = strchr(field, ','); if (end) *end = '\0';
    t->names[c] = malloc(strlen(field)+1);
    if (!t->names[c] || !*field) { free(header); fclose(file); return 0; }
    strcpy(t->names[c], field); field = end ? end+1 : field+strlen(field);
  }
  free(header);
  while (!needed || t->rows < needed) {
    char *line = line_read(file);
    if (!line) break;
    if (!*line) { free(line); continue; }
    if (t->rows == capacity) {
      if (capacity > SIZE_MAX/2/t->columns/sizeof(double)) { free(line); fclose(file); return 0; }
      capacity *= 2;
      double *p = realloc(t->values,capacity*t->columns*sizeof(double));
      if (!p) { free(line); fclose(file); return 0; }
      t->values = p;
    }
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
  if (!t->rows || (needed && t->rows != needed)) { fprintf(stderr, "CSV has %zu cells, need %zu\n", t->rows, needed); return 0; }
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
static int state_columns(const ChemModel *m, const Table *t) {
  for (size_t c = 0; c < t->columns; ++c) {
    const char *label = t->names[c];
    if (strcmp(label,"temperature") && strcmp(label,"pressure") &&
        chem_variable_id(m,label) == SIZE_MAX && chem_custom_id(m,label) == SIZE_MAX) {
      fprintf(stderr,"unknown CSV column: %s\n",label); return 0;
    }
  }
  return 1;
}
static void state_initialize(ChemState *s, const Table *t, size_t first) {
  for (size_t col = 0; col < t->columns; ++col) {
    const char *label = t->names[col];
    size_t v = chem_variable_id(s->model,label), u = chem_custom_id(s->model,label);
    for (size_t c = 0; c < s->cells; ++c) {
      double value = t->values[(first+c)*t->columns+col];
      if (!strcmp(label,"temperature")) s->conditions[c].temperature = value;
      else if (!strcmp(label,"pressure")) s->conditions[c].pressure = value;
      else if (v != SIZE_MAX) s->variables[v*s->stride+c] = value;
      else s->custom[u*s->stride+c] = value;
    }
  }
  for (size_t c = 0; c < s->cells; ++c) chem_ideal_air_density(s->conditions+c);
}
static void state_header(FILE *file, const ChemModel *m, int snapshot) {
  fprintf(file,"%stemperature,pressure",snapshot ? "step,cell," : "");
  for (size_t v = 0; v < m->nvariables; ++v) fprintf(file,",%s",m->species[m->variable_species[v]].name);
  fputc('\n',file);
}
static int state_write(FILE *file, ChemState **states, size_t groups, size_t step, int snapshot) {
  size_t cell = 0;
  for (size_t g = 0; g < groups; ++g) {
    const ChemState *s = states[g];
    for (size_t c = 0; c < s->cells; ++c,++cell) {
      if (snapshot) fprintf(file,"%zu,%zu,",step,cell);
      fprintf(file,"%.17g,%.17g",s->conditions[c].temperature,s->conditions[c].pressure);
      for (size_t v = 0; v < s->model->nvariables; ++v) fprintf(file,",%.17g",s->variables[v*s->stride+c]);
      fputc('\n',file);
    }
  }
  return !ferror(file);
}
static void print_stats(double elapsed, const ChemStats *sum) {
  printf("elapsed_ms=%.3f\n",elapsed);
  printf("internal_steps=%llu accepted=%llu rejected=%llu function_calls=%llu jacobian_updates=%llu decompositions=%llu solves=%llu\n",
    (unsigned long long)sum->internal_steps,(unsigned long long)sum->accepted,(unsigned long long)sum->rejected,
    (unsigned long long)sum->function_calls,(unsigned long long)sum->jacobian_updates,
    (unsigned long long)sum->decompositions,(unsigned long long)sum->solves);
}
#endif
