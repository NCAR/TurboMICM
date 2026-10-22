#include "bench_common.h"

typedef struct
{
  char *cells_path, *photo_path;
  size_t steps, every;
  double dt;
} scenario;
typedef struct
{
  double time, factor;
} photo_point;

static char *trim(char *text)
{
  while (isspace((unsigned char)*text))
    ++text;
  size_t n = strlen(text);
  while (n && isspace((unsigned char)text[n - 1]))
    text[--n] = '\0';
  return text;
}
static char *relative_path(const char *base, const char *name)
{
  const char *slash = strrchr(base, '/');
  size_t prefix = slash && *name != '/' ? (size_t)(slash - base + 1) : 0;
  char *path = malloc(prefix + strlen(name) + 1);
  if (path)
  {
    memcpy(path, base, prefix);
    strcpy(path + prefix, name);
  }
  return path;
}
static int real_argument(const char *text, double *value)
{
  char *end;
  errno = 0;
  *value = strtod(text, &end);
  return !errno && end != text && !*end && isfinite(*value);
}
static int read_scenario(const char *path, scenario *sc)
{
  FILE *file = fopen(path, "r");
  if (!file)
    return 0;
  int ok = 1;
  unsigned seen = 0;
  char *line;
  while ((line = read_line(file)))
  {
    char *comment = strchr(line, '#');
    if (comment)
      *comment = '\0';
    char *key = trim(line);
    if (*key)
    {
      char *equals = strchr(key, '=');
      if (!equals)
        ok = 0;
      else
      {
        *equals = '\0';
        key = trim(key);
        char *value = trim(equals + 1);
        unsigned bit = 0;
        if (!strcmp(key, "cells"))
        {
          bit = 1;
          if (!(seen & bit))
            sc->cells_path = relative_path(path, value);
          ok = ok && *value && sc->cells_path;
        }
        else if (!strcmp(key, "photolysis"))
        {
          bit = 2;
          if (!(seen & bit))
            sc->photo_path = relative_path(path, value);
          ok = ok && *value && sc->photo_path;
        }
        else if (!strcmp(key, "dt"))
        {
          bit = 4;
          ok = ok && real_argument(value, &sc->dt) && sc->dt > 0;
        }
        else if (!strcmp(key, "steps"))
        {
          bit = 8;
          ok = ok && size_argument(value, &sc->steps) && sc->steps > 0;
        }
        else if (!strcmp(key, "snapshot_every"))
        {
          bit = 16;
          ok = ok && size_argument(value, &sc->every);
        }
        else
          ok = 0;
        if (seen & bit)
          ok = 0;
        seen |= bit;
      }
    }
    free(line);
    if (!ok)
      break;
  }
  if (ferror(file))
    ok = 0;
  fclose(file);
  return ok && (seen & 15) == 15;
}
static size_t count_cells(const char *path)
{
  FILE *file = fopen(path, "r");
  if (!file)
    return 0;
  char *line = read_line(file);
  free(line);
  size_t count = 0;
  while ((line = read_line(file)))
  {
    if (*line)
      ++count;
    free(line);
  }
  if (ferror(file))
    count = 0;
  fclose(file);
  return count;
}
static int read_photolysis(const char *path, photo_point **points, size_t *count)
{
  FILE *file = fopen(path, "r");
  if (!file)
    return 0;
  char *line = read_line(file);
  int ok = line && !strcmp(line, "time_s,factor");
  free(line);
  size_t capacity = 0;
  while (ok && (line = read_line(file)))
  {
    if (!*line)
    {
      free(line);
      continue;
    }
    char *comma = strchr(line, ',');
    photo_point p;
    if (!comma)
      ok = 0;
    else
    {
      *comma = '\0';
      ok = real_argument(line, &p.time) && real_argument(comma + 1, &p.factor) && p.factor >= 0;
      if (ok && *count && p.time <= (*points)[*count - 1].time)
        ok = 0;
      if (ok && *count == capacity)
      {
        if (capacity > SIZE_MAX / 2 / sizeof **points)
          ok = 0;
        else
        {
          capacity = capacity ? 2 * capacity : 128;
          photo_point *next = realloc(*points, capacity * sizeof **points);
          if (next)
            *points = next;
          else
            ok = 0;
        }
      }
      if (ok)
        (*points)[(*count)++] = p;
    }
    free(line);
  }
  if (ferror(file))
    ok = 0;
  fclose(file);
  return ok && *count;
}
static double photo_factor(const photo_point *p, size_t count, double time)
{
  if (time <= p[0].time)
    return p[0].factor;
  if (time >= p[count - 1].time)
    return p[count - 1].factor;
  size_t lo = 1, hi = count - 1;
  while (lo < hi)
  {
    size_t mid = lo + (hi - lo) / 2;
    if (time < p[mid].time)
      hi = mid;
    else
      lo = mid + 1;
  }
  double weight = (time - p[lo - 1].time) / (p[lo].time - p[lo - 1].time);
  return (1.0 - weight) * p[lo - 1].factor + weight * p[lo].factor;
}
static int write_snapshots(const ck_model *model, const group *groups, size_t ngroups, double **snapshots,
                           const scenario *sc, size_t nsnap)
{
  FILE *file = fopen("ts1_snapshots.csv", "w");
  if (!file)
    return 0;
  fputs("step,cell,temperature,pressure", file);
  size_t n = ck_model_nspecies(model);
  for (size_t v = 0; v < n; ++v)
    fprintf(file, ",%s", ck_species_name(model, v));
  fputc('\n', file);
  for (size_t snap = 0; snap < nsnap; ++snap)
  {
    size_t step = snap + 1 == nsnap ? sc->steps : (snap + 1) * sc->every;
    for (size_t g = 0; g < ngroups; ++g)
    {
      const ck_state *s = groups[g].state;
      const double *values = snapshots[g] + snap * n * s->cells;
      for (size_t c = 0; c < s->cells; ++c)
      {
        fprintf(file, "%zu,%zu,%.17g,%.17g", step, groups[g].first + c, s->conditions[c].temperature,
                s->conditions[c].pressure);
        for (size_t v = 0; v < n; ++v)
          fprintf(file, ",%.17g", values[v * s->cells + c]);
        fputc('\n', file);
      }
    }
  }
  int ok = !ferror(file);
  if (fclose(file))
    ok = 0;
  return ok;
}
int main(int argc, char **argv)
{
  const char *path = NULL;
  size_t threads = (size_t)omp_get_max_threads(), width = CK_DEFAULT_GROUP_SIZE;
  double rtol = 1e-6, atol = 0;
  const char *env_width = getenv("CK_GROUP_SIZE");
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
    }
    else if (!strcmp(argv[a], "--rtol"))
    {
      if (++a >= argc || !real_argument(argv[a], &rtol) || rtol < 0)
        goto usage;
    }
    else if (!strcmp(argv[a], "--atol"))
    {
      if (++a >= argc || !real_argument(argv[a], &atol) || atol <= 0)
        goto usage;
    }
    else if (!path && argv[a][0] != '-')
      path = argv[a];
    else
      goto usage;
  }
  if (!path || !width)
    goto usage;
  int result = 1;
  scenario sc = {0};
  photo_point *points = NULL;
  size_t npoints = 0, ngroups = 0;
  ck_model *model = NULL;
  group *groups = NULL;
  double **snapshots = NULL, **base = NULL;
  size_t *photo_indices = NULL, nphoto = 0;
  if (!read_scenario(path, &sc) || !read_photolysis(sc.photo_path, &points, &npoints))
  {
    fprintf(stderr, "Cannot read scenario or photolysis table\n");
    goto done;
  }
  size_t cells = count_cells(sc.cells_path);
  if (!cells || ck_model_create(&BENCH_MECHANISM, &model))
    goto done;
  ngroups = (cells - 1) / width + 1;
  if (threads > ngroups)
    threads = ngroups;
  omp_set_num_threads((int)threads);
  if (!getenv("OMP_SCHEDULE"))
    omp_set_schedule(omp_sched_static, 0);
  groups = calloc(ngroups, sizeof *groups);
  snapshots = calloc(ngroups, sizeof *snapshots);
  base = calloc(ngroups, sizeof *base);
  if (!groups || !snapshots || !base)
    goto done;
  size_t n = ck_model_nspecies(model), nc = ck_model_ncustom(model);
  photo_indices = calloc(nc + !nc, sizeof *photo_indices);
  if (!photo_indices)
    goto done;
  for (size_t j = 0; j < nc; ++j)
    if (ck_custom_name(model, j)[0] == 'j')
      photo_indices[nphoto++] = j;
  size_t nsnap = sc.every ? (sc.steps - 1) / sc.every + 1 : 0;
  for (size_t g = 0; g < ngroups; ++g)
  {
    groups[g].first = g * width;
    size_t count = cells - groups[g].first;
    if (count > width)
      count = width;
    if (ck_state_create(model, count, &groups[g].state))
      goto done;
    ck_state *s = groups[g].state;
    s->relative_tolerance = rtol;
    if (atol)
      for (size_t v = 0; v < n; ++v)
        s->absolute_tolerance[v] = atol;
    if (nc > SIZE_MAX / count / sizeof(double) || (nsnap && nsnap > SIZE_MAX / n / count / sizeof(double)))
      goto done;
    base[g] = calloc(nc * count + !nc, sizeof(double));
    snapshots[g] = calloc(nsnap ? nsnap * n * count : 1, sizeof(double));
    if (!base[g] || !snapshots[g])
      goto done;
  }
  if (!read_cells(sc.cells_path, model, groups, ngroups, cells, width))
    goto done;
  for (size_t g = 0; g < ngroups; ++g)
    memcpy(base[g], groups[g].state->custom, nc * groups[g].state->cells * sizeof(double));
  ck_rosenbrock_parameters parameters;
  ck_rosenbrock_parameters_init(3, &parameters);
  /* Start the OpenMP workers without advancing the chemistry. */
#pragma omp parallel
  {
#pragma omp barrier
  }
  double elapsed = 0;
  /* Time each group's steps; buffer snapshots outside each timed solve. The
     shared clock measures the parallel wall time, with output copying excluded
     by pausing at each common snapshot boundary. No chemistry warm-up. */
  size_t first = 0, snap = 0;
  while (first < sc.steps)
  {
    size_t last = sc.steps;
    if (sc.every && sc.every < sc.steps - first)
      last = first + sc.every;
    double t0 = omp_get_wtime();
#pragma omp parallel for schedule(runtime)
    for (size_t g = 0; g < ngroups; ++g)
    {
      ck_state *s = groups[g].state;
      for (size_t step = first; step < last; ++step)
      {
        double factor = photo_factor(points, npoints, (double)step * sc.dt);
        for (size_t p = 0; p < nphoto; ++p)
        {
          size_t j = photo_indices[p];
#pragma omp simd
          for (size_t c = 0; c < s->cells; ++c)
            s->custom[j * s->cells + c] = base[g][j * s->cells + c] * factor;
        }
        groups[g].status = bench_update(s);
        if (groups[g].status)
          break;
        ck_stats stats;
        groups[g].status = bench_solve(s, sc.dt, &parameters, &stats);
        ck_stats_add(&groups[g].stats, &stats);
        if (groups[g].status)
          break;
      }
    }
    elapsed += omp_get_wtime() - t0;
    for (size_t g = 0; g < ngroups; ++g)
    {
      if (groups[g].status)
      {
        fprintf(stderr, "Solve failed in group %zu: %s\n", g, ck_status_string(groups[g].status));
        goto done;
      }
      if (sc.every)
      {
        ck_state *s = groups[g].state;
        memcpy(snapshots[g] + snap * n * s->cells, s->y, n * s->cells * sizeof(double));
      }
    }
    ++snap;
    first = last;
  }
  if (!write_final(model, groups, ngroups) ||
      (sc.every && !write_snapshots(model, groups, ngroups, snapshots, &sc, nsnap)))
    goto done;
  ck_stats sum = {0};
  for (size_t g = 0; g < ngroups; ++g)
    ck_stats_add(&sum, &groups[g].stats);
  printf("elapsed_ms=%.3f\n", elapsed * 1000);
  printf("internal_steps=%llu accepted=%llu rejected=%llu function_calls=%llu jacobian_updates=%llu "
         "decompositions=%llu solves=%llu\n",
         (unsigned long long)sum.internal_steps, (unsigned long long)sum.accepted,
         (unsigned long long)sum.rejected, (unsigned long long)sum.function_calls,
         (unsigned long long)sum.jacobian_updates, (unsigned long long)sum.decompositions,
         (unsigned long long)sum.solves);
  result = 0;
done:
  for (size_t g = 0; g < ngroups; ++g)
  {
    if (groups)
      ck_state_free(groups[g].state);
    if (base)
      free(base[g]);
    if (snapshots)
      free(snapshots[g]);
  }
  free(groups);
  free(base);
  free(snapshots);
  ck_model_free(model);
  free(points);
  free(photo_indices);
  free(sc.cells_path);
  free(sc.photo_path);
  return result;
usage:
  fprintf(stderr, "usage: %s [--threads N] [--group-size N] [--rtol R] [--atol A] SCENARIO\n", argv[0]);
  return 2;
}
