/* Diagnostic only: does the solver select the generated kernels?
 *
 *   kernel_check CELLS STEPS DT CELLS_CSV
 *
 * Uses only the entry's public API (kinetics.h). It copies the generated
 * chem_model, keeps its kernel table's sizes, and replaces each kernel by a
 * wrapper that counts the call and forwards to the generated function. The
 * solver uses a kernel only when its own selection test passes (src/solver.c),
 * so nonzero counts mean the generated kernels ran; zero means generic code.
 * Same problem as the benchmark: one thread, default grouping (group size 0),
 * MICM's default Rosenbrock settings, one warm-up step, then STEPS steps of
 * DT. Writes ts1_final.csv in the benchmark's format. */
#include "kinetics.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

extern const chem_mechanism chem_model;
static const chem_specialized_kernels *generated;
static unsigned long long rates_calls, custom_calls, evaluate_calls, factor_calls, solve_calls;

static void count_rates(double t, double p, const double *c, double *k)
{ ++rates_calls; generated->rates(t, p, c, k); }
static void count_custom(double t, double p, const double *c, double *k)
{ ++custom_calls; generated->custom_rates(t, p, c, k); }
static int count_evaluate(const double *y, const double *k, double *f, double *j,
                          const chem_sparse_pattern *pattern)
{ ++evaluate_calls; return generated->evaluate(y, k, f, j, pattern); }
static int count_factor(const chem_sparse_pattern *pattern, double alpha,
                        const double *j, double *lu, void *context)
{ ++factor_calls; return generated->factor(pattern, alpha, j, lu, context); }
static void count_solve(const chem_sparse_pattern *pattern, const double *lu,
                        double *b, void *context)
{ ++solve_calls; generated->solve(pattern, lu, b, context); }

static void add(chem_solver_stats *total, const chem_solver_stats *s)
{
  total->internal_steps += s->internal_steps; total->function_calls += s->function_calls;
  total->decompositions += s->decompositions; total->solves += s->solves;
}

int main(int argc, char **argv)
{
  static char line[262144];
  size_t cells, steps, n, nc, c, i, columns = 0;
  int column_of[4096];
  double dt, *t, *p, *custom, *y;
  chem_specialized_kernels wrapped;
  chem_mechanism m;
  chem_rosenbrock_parameters parameters;
  chem_solver_stats stats, total = {0};
  chem_solver *solver;
  FILE *file;
  char *field, *cursor;
  if (argc != 5) { fprintf(stderr, "usage: %s CELLS STEPS DT CELLS_CSV\n", argv[0]); return 2; }
  cells = strtoul(argv[1], NULL, 10); steps = strtoul(argv[2], NULL, 10); dt = strtod(argv[3], NULL);
  if (!chem_model.kernels) { fprintf(stderr, "no generated kernels in this model\n"); return 1; }
  generated = chem_model.kernels;
  wrapped = *generated; /* same species_count and pattern_nonzeros */
  wrapped.rates = generated->rates ? count_rates : NULL;
  wrapped.custom_rates = generated->custom_rates ? count_custom : NULL;
  wrapped.evaluate = generated->evaluate ? count_evaluate : NULL;
  wrapped.factor = generated->factor ? count_factor : NULL;
  wrapped.solve = generated->solve ? count_solve : NULL;
  m = chem_model;
  m.kernels = &wrapped;
  n = m.species_count; nc = m.custom_count;
  t = malloc(cells * sizeof *t); p = malloc(cells * sizeof *p);
  custom = malloc((cells * nc ? cells * nc : 1) * sizeof *custom); y = malloc(cells * n * sizeof *y);
  for (c = 0; c < cells; ++c) {
    t[c] = m.default_temperature; p[c] = m.default_pressure;
    for (i = 0; i < nc; ++i) custom[c * nc + i] = m.custom_parameters[i].default_value;
    for (i = 0; i < n; ++i) y[c * n + i] = m.species[i].initial_value;
  }
  if (!(file = fopen(argv[4], "r")) || !fgets(line, sizeof line, file)) { fprintf(stderr, "cannot read %s\n", argv[4]); return 1; }
  line[strcspn(line, "\r\n")] = '\0';
  for (cursor = line; (field = strsep(&cursor, ",")) != NULL; ++columns) {
    column_of[columns] = -1;
    if (!strcmp(field, "temperature")) column_of[columns] = -2;
    else if (!strcmp(field, "pressure")) column_of[columns] = -3;
    for (i = 0; i < n; ++i) if (!strcmp(field, m.species[i].name)) column_of[columns] = (int)i;
    for (i = 0; i < nc; ++i) if (!strcmp(field, m.custom_parameters[i].name)) column_of[columns] = -1000 - (int)i;
    if (column_of[columns] == -1) { fprintf(stderr, "unknown column %s\n", field); return 1; }
  }
  for (c = 0; c < cells; ++c) {
    size_t k = 0;
    if (!fgets(line, sizeof line, file)) { fprintf(stderr, "too few rows\n"); return 1; }
    line[strcspn(line, "\r\n")] = '\0';
    for (cursor = line; (field = strsep(&cursor, ",")) != NULL && k < columns; ++k) {
      const double v = strtod(field, NULL);
      const int col = column_of[k];
      if (col == -2) t[c] = v; else if (col == -3) p[c] = v;
      else if (col <= -1000) custom[c * nc + (size_t)(-1000 - col)] = v;
      else y[c * n + (size_t)col] = v;
    }
  }
  fclose(file);
#ifdef _OPENMP
  omp_set_num_threads(1);
#endif
  chem_rosenbrock_default(&parameters);
  if (!(solver = chem_solver_create(&m, cells, t, p, custom, &parameters, 0))) return 1;
  for (i = 0; i <= steps; ++i) { /* i == 0 is the warm-up step */
    if (chem_solver_step(solver, y, dt, &stats) != 0) { fprintf(stderr, "step %zu failed\n", i); return 1; }
    add(&total, &stats);
  }
  file = fopen("ts1_final.csv", "w");
  fputs("temperature,pressure", file);
  for (i = 0; i < n; ++i) fprintf(file, ",%s", m.species[i].name);
  fputc('\n', file);
  for (c = 0; c < cells; ++c) {
    fprintf(file, "%.17g,%.17g", t[c], p[c]);
    for (i = 0; i < n; ++i) fprintf(file, ",%.17g", y[c * n + i]);
    fputc('\n', file);
  }
  fclose(file);
  printf("kernel table: species_count=%zu pattern_nonzeros=%zu\n", generated->species_count, generated->pattern_nonzeros);
  printf("solver totals incl. warm-up: internal_steps=%llu function_calls=%llu decompositions=%llu solves=%llu\n",
         (unsigned long long)total.internal_steps, (unsigned long long)total.function_calls,
         (unsigned long long)total.decompositions, (unsigned long long)total.solves);
  printf("generated kernel calls: rates=%llu custom_rates=%llu evaluate=%llu factor=%llu solve=%llu\n",
         rates_calls, custom_calls, evaluate_calls, factor_calls, solve_calls);
  printf("expected if every call is generated: rates=%zu evaluate=%llu factor=%llu solve=%llu\n", cells,
         (unsigned long long)(total.function_calls * cells), (unsigned long long)(total.decompositions * cells),
         (unsigned long long)(total.solves * cells));
  chem_solver_destroy(solver);
  return 0;
}
