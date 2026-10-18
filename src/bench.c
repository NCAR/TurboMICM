#include "chem.h"
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef CK_DEFAULT_GROUP_SIZE
#define CK_DEFAULT_GROUP_SIZE 32
#endif
extern const ck_mechanism ts1;
#include "bench_io.h"

int main(int argc, char **argv)
{
  size_t threads = (size_t)omp_get_max_threads(), width = CK_DEFAULT_GROUP_SIZE;
  const char *arguments[4];
  size_t nargs = 0;
  const char *env_width = getenv("CK_GROUP_SIZE");
  int explicit_width = env_width != NULL;
  if (env_width && (!size_argument(env_width, &width) || !width))
    goto usage;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "--threads"))
    {
      if (++a >= argc || !size_argument(argv[a], &threads) || !threads || threads > INT_MAX)
        goto usage;
    }
    else if (!strcmp(argv[a], "--group-size"))
    {
      if (++a >= argc || !size_argument(argv[a], &width) || !width)
        goto usage;
      explicit_width = 1;
    }
    else if (nargs < 4)
      arguments[nargs++] = argv[a];
    else
      goto usage;
  }
  size_t cells, steps;
  if (nargs != 4 || !size_argument(arguments[0], &cells) || !cells || !size_argument(arguments[1], &steps))
    goto usage;
  char *end;
  errno = 0;
  double dt = strtod(arguments[2], &end);
  if (errno || *end || !(dt > 0) || !isfinite(dt))
    goto usage;
  if (!width)
    goto usage;
  if (!getenv("OMP_SCHEDULE"))
    omp_set_schedule(omp_sched_static, 0);
  if (!explicit_width)
  {
    size_t balanced = (cells - 1) / threads + 1;
    if (width > balanced)
      width = balanced;
  }
  if (width > cells)
    width = cells;
  size_t ngroups = (cells - 1) / width + 1;
  if (threads > ngroups)
    threads = ngroups;
  omp_set_num_threads((int)threads);
  ck_model *model = NULL;
  group *groups = NULL;
  ck_status status = ck_model_create(&ts1, &model);
  if (status)
  {
    fprintf(stderr, "%s\n", ck_status_string(status));
    return 1;
  }
  groups = calloc(ngroups, sizeof *groups);
  int exit_status = 1;
  if (!groups)
    goto done;
  for (size_t g = 0; g < ngroups; ++g)
  {
    groups[g].first = g * width;
    size_t count = cells - g * width;
    if (count > width)
      count = width;
    status = ck_state_create(model, count, &groups[g].state);
    if (status)
    {
      fprintf(stderr, "%s\n", ck_status_string(status));
      goto done;
    }
  }
  if (!read_cells(arguments[3], model, groups, ngroups, cells, width))
    goto done;
  ck_rosenbrock_parameters parameters;
  ck_rosenbrock_parameters_init(3, &parameters);
#pragma omp parallel for schedule(runtime)
  for (size_t g = 0; g < ngroups; ++g)
  {
    ck_stats ignored;
    groups[g].status = ck_solve(groups[g].state, dt, &parameters, &ignored);
  }
  for (size_t g = 0; g < ngroups; ++g)
    if (groups[g].status)
    {
      fprintf(stderr, "Warm-up failed: %s\n", ck_status_string(groups[g].status));
      goto done;
    }
  double t0 = omp_get_wtime();
#pragma omp parallel for schedule(runtime)
  for (size_t g = 0; g < ngroups; ++g)
  {
    for (size_t step = 0; step < steps; ++step)
    {
      ck_stats stats;
      groups[g].status = ck_solve(groups[g].state, dt, &parameters, &stats);
      ck_stats_add(&groups[g].stats, &stats);
      if (groups[g].status)
        break;
    }
  }
  double elapsed_ms = (omp_get_wtime() - t0) * 1000;
  ck_stats sum = {0};
  for (size_t g = 0; g < ngroups; ++g)
  {
    if (groups[g].status)
    {
      fprintf(stderr, "Solve failed: %s\n", ck_status_string(groups[g].status));
      goto done;
    }
    ck_stats_add(&sum, &groups[g].stats);
  }
  if (!write_final(model, groups, ngroups))
  {
    fprintf(stderr, "Cannot write ts1_final.csv\n");
    goto done;
  }
  printf("elapsed_ms=%.3f\n", elapsed_ms);
  printf("internal_steps=%llu accepted=%llu rejected=%llu function_calls=%llu jacobian_updates=%llu "
         "decompositions=%llu solves=%llu\n",
         (unsigned long long)sum.internal_steps, (unsigned long long)sum.accepted,
         (unsigned long long)sum.rejected, (unsigned long long)sum.function_calls,
         (unsigned long long)sum.jacobian_updates, (unsigned long long)sum.decompositions,
         (unsigned long long)sum.solves);
  exit_status = 0;
done:
  if (groups)
    for (size_t g = 0; g < ngroups; ++g)
      ck_state_free(groups[g].state);
  free(groups);
  ck_model_free(model);
  return exit_status;
usage:
  fprintf(stderr, "usage: %s [--threads N] [--group-size N] CELLS STEPS DT CELLS_CSV\n", argv[0]);
  return 2;
}
