/* SPDX-License-Identifier: Apache-2.0 */
#include "chem.h"
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef CH_SPECIALIZED
#include "ts1_specialized.h"
#define ch_build_ts1 ch_build_ts1_specialized
#define ch_init_ts1 ch_init_ts1_specialized
#define ch_update_rates ch_rates_ts1_specialized
#define ch_solve ch_solve_ts1_specialized
#endif

#ifndef CH_DEFAULT_GROUP_SIZE
#define CH_DEFAULT_GROUP_SIZE 32
#endif

ch_model *ch_build_ts1(void);
int ch_init_ts1(ch_state *);

typedef struct { ch_state *state; size_t first; int status; ch_stats stats; } group;
typedef struct { char *name; int kind, index; } column;

static char *line_read(FILE *file) {
  size_t n = 0, cap = 1024;
  char *line = malloc(cap);
  if (!line) return NULL;
  int ch;
  while ((ch = fgetc(file)) != EOF && ch != '\n') {
    if (n+1 == cap) {
      if (cap > SIZE_MAX/2) { free(line); return NULL; }
      char *next = realloc(line,cap*2);
      if (!next) { free(line); return NULL; }
      line = next; cap *= 2;
    }
    line[n++] = (char)ch;
  }
  if (ch == EOF && n == 0) { free(line); return NULL; }
  if (n && line[n-1] == '\r') --n;
  line[n] = 0;
  return line;
}

/* The reference format has unquoted names and decimal values, no embedded commas. */
static char *field_next(char **cursor) {
  if (!*cursor) return NULL;
  char *start = *cursor, *comma = strchr(start,',');
  if (comma) { *comma = 0; *cursor = comma+1; }
  else *cursor = NULL;
  return start;
}
static int size_parse(const char *text, size_t *value) {
  char *end; errno = 0;
  if (!text[0] || text[0] == '-') return 0;
  unsigned long long n = strtoull(text,&end,10);
  if (errno || *end || n != (size_t)n) return 0;
  *value = (size_t)n; return 1;
}
static int double_parse(const char *text, double *value) {
  char *end; errno = 0; *value = strtod(text,&end);
  return end != text && !*end && isfinite(*value) && errno != ERANGE;
}
static int initialize_csv(const char *path, const ch_model *m, group *groups, size_t cells, size_t group_size) {
  FILE *file = fopen(path,"r");
  if (!file) { fprintf(stderr,"Cannot open input: %s\n",path); return CH_INVALID; }
  char *header = line_read(file), *line = NULL;
  column *cols = NULL;
  size_t ncols = 0;
  int status = CH_INVALID;
  if (!header) goto done;
  char *cursor = header, *field;
  while ((field = field_next(&cursor))) {
    column col = {field,0,0};
    if (!strcmp(field,"temperature")) col.kind = 0;
    else if (!strcmp(field,"pressure")) col.kind = 1;
    else if ((col.index = ch_model_species_index(m,field)) >= 0) col.kind = 2;
    else if ((col.index = ch_model_parameter_index(m,field)) >= 0) col.kind = 3;
    else { fprintf(stderr,"Unknown input column: %s\n",field); goto done; }
    for (size_t i = 0; i < ncols; ++i) if (!strcmp(cols[i].name,field)) goto done;
    column *next = realloc(cols,(ncols+1)*sizeof(*cols));
    if (!next) { status = CH_MEMORY; goto done; }
    cols = next; cols[ncols++] = col;
  }
  for (size_t c = 0; c < cells; ++c) {
    line = line_read(file);
    if (!line) { fprintf(stderr,"Input has fewer than %zu cells\n",cells); goto done; }
    cursor = line;
    ch_state *s = groups[c/group_size].state;
    size_t local = c-groups[c/group_size].first, count = ch_state_cells(s);
    for (size_t i = 0; i < ncols; ++i) {
      field = field_next(&cursor);
      double value;
      if (!field || !double_parse(field,&value)) {
        fprintf(stderr,"Invalid value in cell %zu, column %s\n",c,cols[i].name); goto done;
      }
      switch (cols[i].kind) {
        case 0: ch_state_conditions(s)[local].temperature = value; break;
        case 1: ch_state_conditions(s)[local].pressure = value; break;
        case 2: ch_state_values(s)[(size_t)cols[i].index*count+local] = value; break;
        case 3: ch_state_parameters(s)[(size_t)cols[i].index*count+local] = value; break;
      }
    }
    if (cursor) { fprintf(stderr,"Extra value in cell %zu\n",c); goto done; }
    ch_conditions *condition = &ch_state_conditions(s)[local];
    if (condition->temperature <= 0 || condition->pressure <= 0) goto done;
    ch_ideal_air_density(condition);
    free(line); line = NULL;
  }
  status = CH_OK;
done:
  free(line); free(header); free(cols); fclose(file); return status;
}

