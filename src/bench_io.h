/* Shared ISO C CSV input/output for the two benchmark drivers. */
typedef struct
{
  ck_state *state;
  ck_stats stats;
  ck_status status;
  size_t first;
} group;

static int size_argument(const char *text, size_t *value)
{
  if (!text || !isdigit((unsigned char)*text))
    return 0;
  char *end;
  errno = 0;
  unsigned long long parsed = strtoull(text, &end, 10);
  if (errno || *end || parsed > SIZE_MAX)
    return 0;
  *value = (size_t)parsed;
  return 1;
}
static char *read_line(FILE *file)
{
  size_t cap = 1024, size = 0;
  char *line = malloc(cap);
  if (!line)
    return NULL;
  int c;
  while ((c = fgetc(file)) != EOF && c != '\n')
  {
    if (size + 1 >= cap)
    {
      if (cap > SIZE_MAX / 2)
      {
        free(line);
        return NULL;
      }
      char *next = realloc(line, cap * 2);
      if (!next)
      {
        free(line);
        return NULL;
      }
      cap *= 2;
      line = next;
    }
    line[size++] = (char)c;
  }
  if (c == EOF && size == 0)
  {
    free(line);
    return NULL;
  }
  if (size && line[size - 1] == '\r')
    --size;
  line[size] = '\0';
  return line;
}
static int read_cells(const char *path, const ck_model *model, group *groups, size_t ngroups, size_t ncells,
                      size_t width)
{
  FILE *file = fopen(path, "r");
  if (!file)
  {
    fprintf(stderr, "Cannot open %s\n", path);
    return 0;
  }
  char *header = read_line(file);
  if (!header)
  {
    fclose(file);
    return 0;
  }
  size_t ncols = 1;
  for (char *p = header; *p; ++p)
    if (*p == ',')
      ++ncols;
  size_t *index = calloc(ncols, sizeof *index);
  int *kind = calloc(ncols, sizeof *kind);
  if (!index || !kind)
  {
    free(index);
    free(kind);
    free(header);
    fclose(file);
    return 0;
  }
  int ok = 0;
  char *field = header;
  for (size_t col = 0; col < ncols; ++col)
  {
    char *end = strchr(field, ',');
    if (end)
      *end = '\0';
    if (!strcmp(field, "temperature"))
      kind[col] = 0;
    else if (!strcmp(field, "pressure"))
      kind[col] = 1;
    else if (ck_species_index(model, field, &index[col]) == CK_OK)
      kind[col] = 2;
    else if (ck_custom_index(model, field, &index[col]) == CK_OK)
      kind[col] = 3;
    else
    {
      fprintf(stderr, "Unknown CSV column: %s\n", field);
      goto done;
    }
    /* Repeated columns are ambiguous and usually indicate a malformed input. */
    for (size_t prev = 0; prev < col; ++prev)
      if (kind[prev] == kind[col] && index[prev] == index[col])
      {
        fprintf(stderr, "Duplicate CSV column\n");
        goto done;
      }
    field = end ? end + 1 : field + strlen(field);
  }
  for (size_t cell = 0; cell < ncells; ++cell)
  {
    char *line;
    do
    {
      line = read_line(file);
      if (!line)
      {
        fprintf(stderr, "Too few CSV rows\n");
        goto done;
      }
    } while (*line == '\0' && (free(line), 1));
    ck_state *s = groups[cell / width].state;
    size_t local = cell - groups[cell / width].first;
    field = line;
    for (size_t col = 0; col < ncols; ++col)
    {
      char *end = strchr(field, ',');
      if (end)
        *end = '\0';
      char *tail;
      errno = 0;
      double v = strtod(field, &tail);
      if (errno || tail == field || *tail || !isfinite(v) || ((col + 1 < ncols) != (end != NULL)))
      {
        fprintf(stderr, "Invalid CSV value or field count in row %zu\n", cell + 1);
        free(line);
        goto done;
      }
      if (kind[col] == 0)
        s->conditions[local].temperature = v;
      else if (kind[col] == 1)
        s->conditions[local].pressure = v;
      else if (kind[col] == 2)
        s->y[index[col] * s->cells + local] = v;
      else
        s->custom[index[col] * s->cells + local] = v;
      field = end ? end + 1 : tail;
    }
    if (!(s->conditions[local].temperature > 0) || !(s->conditions[local].pressure > 0))
    {
      fprintf(stderr, "Temperature and pressure must be positive\n");
      free(line);
      goto done;
    }
    ck_ideal_air_density(&s->conditions[local]);
    free(line);
  }
  for (size_t g = 0; g < ngroups; ++g)
    if (ck_update_rates(groups[g].state))
      goto done;
  ok = 1;
done:
  free(index);
  free(kind);
  free(header);
  fclose(file);
  return ok;
}
static int write_final(const ck_model *model, const group *groups, size_t ngroups)
{
  FILE *file = fopen("ts1_final.csv", "w");
  if (!file)
    return 0;
  fputs("temperature,pressure", file);
  size_t n = ck_model_nspecies(model);
  for (size_t i = 0; i < n; ++i)
    fprintf(file, ",%s", ck_species_name(model, i));
  fputc('\n', file);
  for (size_t g = 0; g < ngroups; ++g)
  {
    const ck_state *s = groups[g].state;
    for (size_t c = 0; c < s->cells; ++c)
    {
      fprintf(file, "%.17g,%.17g", s->conditions[c].temperature, s->conditions[c].pressure);
      for (size_t i = 0; i < n; ++i)
        fprintf(file, ",%.17g", s->y[i * s->cells + c]);
      fputc('\n', file);
    }
  }
  int ok = !ferror(file);
  if (fclose(file))
    ok = 0;
  return ok;
}
