/* Shared standard-C CSV input helpers for the two native drivers. */
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
