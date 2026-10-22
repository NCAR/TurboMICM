/* Shared CSV and group utilities for the TS1 drivers. */
typedef struct { ChemState *state; size_t cells; ChemStats stats; ChemStatus status; } Group;

static char *read_line(FILE *f)
{
    size_t capacity = 1024, used = 0;
    char *line = malloc(capacity);
    if (!line) return NULL;
    int c;
    while ((c = fgetc(f)) != EOF && c != '\n') {
        if (used + 1 >= capacity) {
            if (capacity > SIZE_MAX / 2) { free(line); return NULL; }
            capacity *= 2;
            char *grown = realloc(line, capacity);
            if (!grown) { free(line); return NULL; }
            line = grown;
        }
        line[used++] = (char)c;
    }
    if (!used && c == EOF) { free(line); return NULL; }
    if (used && line[used - 1] == '\r') --used;
    line[used] = 0;
    return line;
}
static size_t split(char *line, char ***out)
{
    size_t count = 1;
    for (char *p = line; *p; ++p) if (*p == ',') ++count;
    char **fields = malloc(count * sizeof(*fields));
    if (!fields) return 0;
    size_t i = 0;
    fields[i++] = line;
    for (char *p = line; *p; ++p) if (*p == ',') { *p = 0; fields[i++] = p + 1; }
    *out = fields;
    return count;
}
static int parse_size(const char *text, size_t *value)
{
    char *end;
    if (!*text || *text == '-') return 0;
    errno = 0;
    unsigned long long n = strtoull(text, &end, 10);
    if (errno || *end || n > SIZE_MAX) return 0;
    *value = (size_t)n;
    return 1;
}
static int parse_double(const char *text, double *value)
{
    char *end;
    errno = 0;
    *value = strtod(text, &end);
    return *text && !*end && !errno && isfinite(*value);
}
static void add_stats(ChemStats *sum, const ChemStats *a)
{
    sum->internal_steps += a->internal_steps; sum->accepted += a->accepted; sum->rejected += a->rejected;
    sum->function_calls += a->function_calls; sum->jacobian_updates += a->jacobian_updates;
    sum->decompositions += a->decompositions; sum->solves += a->solves;
}
/* Maps the header to a compact index: T=0, P=1, then species, then parameters. */
static int column_index(const char *name)
{
    if (!strcmp(name, "temperature")) return 0;
    if (!strcmp(name, "pressure")) return 1;
    for (size_t i = 0; i < ts1_data.species_count; ++i)
        if (!strcmp(name, ts1_data.species_names[i])) return (int)i + 2;
    for (size_t i = 0; i < ts1_data.parameter_count; ++i)
        if (!strcmp(name, ts1_data.parameter_names[i])) return (int)(i + ts1_data.species_count) + 2;
    return -1;
}
static int read_cells(const char *path, Group *groups, size_t cells, size_t group_size)
{
    FILE *file = fopen(path, "r");
    if (!file) { fprintf(stderr, "cannot open %s\n", path); return 0; }
    char *line = read_line(file), **header = NULL;
    if (!line) { fclose(file); return 0; }
    size_t columns = split(line, &header);
    int *indices = malloc(columns * sizeof(*indices));
    int ok = 0;
    if (!columns || !indices) goto done;
    for (size_t j = 0; j < columns; ++j) {
        indices[j] = column_index(header[j]);
        if (indices[j] < 0) { fprintf(stderr, "unknown CSV column: %s\n", header[j]); goto done; }
        for (size_t k = 0; k < j; ++k)
            if (indices[k] == indices[j]) { fprintf(stderr, "duplicate CSV column: %s\n", header[j]); goto done; }
    }
    free(header); header = NULL; free(line); line = NULL;
    for (size_t i = 0; i < cells;) {
        line = read_line(file);
        if (!line) { fprintf(stderr, "CSV has fewer than %zu rows\n", cells); goto done; }
        if (!*line) { free(line); line = NULL; continue; }
        char **fields = NULL;
        size_t count = split(line, &fields);
        if (count != columns) { free(fields); fprintf(stderr, "CSV row %zu has wrong field count\n", i + 1); goto done; }
        Group *g = &groups[i / group_size];
        size_t c = i % group_size, nc = g->cells;
        for (size_t j = 0; j < columns; ++j) {
            double value;
            if (!parse_double(fields[j], &value)) { free(fields); fprintf(stderr, "invalid number at row %zu\n", i + 1); goto done; }
            int index = indices[j];
            if (index == 0) chem_temperatures(g->state)[c] = value;
            else if (index == 1) chem_pressures(g->state)[c] = value;
            else if ((size_t)index < ts1_data.species_count + 2)
                chem_concentrations(g->state)[((size_t)index - 2) * nc + c] = value;
            else chem_parameters(g->state)[((size_t)index - 2 - ts1_data.species_count) * nc + c] = value;
        }
        free(fields); free(line); line = NULL;
        ++i;
    }
    ok = 1;
done:
    free(indices); free(header); free(line); fclose(file);
    return ok;
}
static int write_final(const Group *groups, size_t ng)
{
    FILE *file = fopen("ts1_final.csv", "w");
    if (!file) return 0;
    fprintf(file, "temperature,pressure");
    for (size_t i = 0; i < ts1_data.species_count; ++i) fprintf(file, ",%s", ts1_data.species_names[i]);
    fputc('\n', file);
    for (size_t g = 0; g < ng; ++g) {
        size_t nc = groups[g].cells;
        ChemState *s = groups[g].state;
        for (size_t c = 0; c < nc; ++c) {
            fprintf(file, "%.17g,%.17g", chem_temperatures(s)[c], chem_pressures(s)[c]);
            for (size_t i = 0; i < ts1_data.species_count; ++i)
                fprintf(file, ",%.17g", chem_concentrations(s)[i * nc + c]);
            fputc('\n', file);
        }
    }
    int ok = !ferror(file);
    return fclose(file) == 0 && ok;
}
