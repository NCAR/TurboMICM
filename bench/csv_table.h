/* A CSV table with a header row and numeric rows, as the benchmark drivers
 * read CELLS_CSV and photolysis files. `tmc_program` names the program in
 * error messages. */
#ifndef TMC_BENCH_CSV_TABLE_H
#define TMC_BENCH_CSV_TABLE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *tmc_program = "bench";

typedef struct
{
  int n_columns;
  char (*names)[128];
  int n_rows;
  double *values; /* n_rows x n_columns */
} table;

static void fail(const char *message)
{
  fprintf(stderr, "%s: %s\n", tmc_program, message);
  exit(1);
}

static table read_csv(const char *path)
{
  FILE *f = fopen(path, "r");
  if (!f)
  {
    fprintf(stderr, "%s: cannot open %s\n", tmc_program, path);
    exit(1);
  }
  table t = { 0, NULL, 0, NULL };
  char *line = NULL;
  size_t cap = 0;
  ssize_t len;
  int capacity_rows = 0;
  while ((len = getline(&line, &cap, f)) > 0)
  {
    while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
      line[--len] = '\0';
    if (len == 0)
      continue;
    if (t.n_columns == 0)
    {
      int cap_cols = 16;
      t.names = malloc((size_t)cap_cols * sizeof *t.names);
      for (char *tok = strtok(line, ","); tok; tok = strtok(NULL, ","))
      {
        if (t.n_columns == cap_cols)
        {
          cap_cols *= 2;
          t.names = realloc(t.names, (size_t)cap_cols * sizeof *t.names);
        }
        snprintf(t.names[t.n_columns++], sizeof t.names[0], "%s", tok);
      }
      continue;
    }
    if (t.n_rows == capacity_rows)
    {
      capacity_rows = capacity_rows ? 2 * capacity_rows : 1024;
      t.values = realloc(t.values, (size_t)capacity_rows * (size_t)t.n_columns * sizeof(double));
    }
    double *row = t.values + (size_t)t.n_rows * (size_t)t.n_columns;
    int k = 0;
    for (char *tok = strtok(line, ","); tok; tok = strtok(NULL, ","))
    {
      if (k == t.n_columns)
        fail("a CSV row has more fields than the header");
      row[k++] = strtod(tok, NULL);
    }
    if (k != t.n_columns)
      fail("a CSV row has fewer fields than the header");
    ++t.n_rows;
  }
  free(line);
  fclose(f);
  return t;
}

#endif
