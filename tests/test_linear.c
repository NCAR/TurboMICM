/* Port of the LU reconstruction and A*x=b tests, using a portable seeded PRNG.
   Numeric shapes, distributions and tolerances follow the reference policies. */
#include "test.h"
#include <stdint.h>

static uint32_t random_state;
static double uniform(void)
{
  random_state = random_state * UINT32_C(1664525) + UINT32_C(1013904223);
  return (random_state + 0.5) / 4294967296.0;
}
static double lognormal(double deviation)
{
  return exp(-2 + deviation * sqrt(-2 * log(uniform())) * cos(2 * CK_PI * uniform()));
}
static double entry(const ck_model *m, const double *a, size_t cells, size_t row, size_t col, size_t c)
{
  size_t off = ck_jacobian_index(m, row, col);
  return off == SIZE_MAX ? 0 : a[off * cells + c];
}
static int numeric_case(size_t n, size_t cells, int kind, int solve, int separate, ck_lu_algorithm algorithm,
                        double initial)
{
  unsigned char *pattern = calloc(n * n, 1);
  double *dense = calloc(n * n * cells, sizeof(double));
  CHECK(pattern && dense);
  random_state = 12345;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      pattern[i * n + j] = i == j || kind == 0 || (kind != 2 && uniform() < 0.5);
  ck_model *m = NULL;
  CHECK(ck_linear_model_create(n, pattern, &m) == CK_OK);
  size_t nnz = ck_jacobian_size(m);
  double *a = calloc(nnz * cells, sizeof(double)), *lu = malloc(nnz * cells * sizeof(double));
  double *l = malloc(nnz * cells * sizeof(double)), *u = malloc(nnz * cells * sizeof(double));
  double *b = calloc(n * cells, sizeof(double)), *x = calloc(n * cells, sizeof(double));
  CHECK(a && lu && l && u && b && x);
  for (size_t i = 0; i < nnz * cells; ++i)
    lu[i] = l[i] = u[i] = initial;
  const double example[] = {2, -1, -2, -4, 6, 3, -4, -2, 8};
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      if (pattern[i * n + j])
      {
        size_t off = ck_jacobian_index(m, i, j);
        for (size_t c = 0; c < cells; ++c)
        {
          double value = kind == 0 ? example[i * n + j] : lognormal(kind == 2 ? 4 : 2);
          dense[(i * n + j) * cells + c] = a[off * cells + c] = value;
        }
      }
  for (size_t i = 0; i < n; ++i)
    for (size_t c = 0; c < cells; ++c)
      b[i * cells + c] = kind == 0 ? (i == 0 ? 23 : i == 1 ? 42 : 9) : lognormal(2);
  memcpy(x, b, n * cells * sizeof(double));
  if (separate)
    CHECK(ck_factor_separate(m, cells, a, l, u, algorithm) == CK_OK);
  else
  {
    memcpy(lu, a, nnz * cells * sizeof(double));
    ck_factor(m, cells, lu, algorithm);
  }
  if (solve)
  {
    if (separate)
      ck_linear_solve_separate(m, cells, l, u, x);
    else
      ck_linear_solve(m, cells, lu, x);
    for (size_t c = 0; c < cells; ++c)
      for (size_t i = 0; i < n; ++i)
      {
        double result = 0;
        for (size_t j = 0; j < n; ++j)
          if (pattern[i * n + j])
            result += dense[(i * n + j) * cells + c] * x[j * cells + c];
        double tolerance = kind == 0 || kind == 2 ? 1e-5 : n == 30 ? 2e-6 : 1e-6;
        NEAR(result, b[i * cells + c], tolerance);
      }
  }
  else
  {
    for (size_t c = 0; c < cells; ++c)
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
        {
          double result = 0;
          for (size_t k = 0; k < n; ++k)
          {
            double left = i < k ? 0 : i == k ? 1 : entry(m, separate ? l : lu, cells, i, k, c);
            double right = k > j ? 0 : entry(m, separate ? u : lu, cells, k, j, c);
            result += left * right;
          }
          NEAR(result, dense[(i * n + j) * cells + c], kind == 0 || kind == 2 ? 1e-10 : 1e-9);
          if (separate)
          {
            if (i < j)
              NEAR(entry(m, l, cells, i, j, c), 0, 0);
            if (j < i)
              NEAR(entry(m, u, cells, i, j, c), 0, 0);
          }
        }
  }
  free(pattern);
  free(dense);
  free(a);
  free(lu);
  free(l);
  free(u);
  free(b);
  free(x);
  ck_model_free(m);
  return 0;
}
static int markowitz(void)
{
  size_t n = 50;
  unsigned char pattern[2500], ordered[2500];
  size_t perm[50];
  random_state = 1;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      pattern[i * n + j] = i == j || uniform() < 0.5;
  CHECK(ck_markowitz_order(n, pattern, perm) == CK_OK);
  size_t sum = 0;
  for (size_t i = 0; i < n; ++i)
  {
    CHECK(perm[i] < n);
    sum += perm[i];
    for (size_t j = 0; j < i; ++j)
      CHECK(perm[i] != perm[j]);
  }
  CHECK(sum == n * (n - 1) / 2);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      ordered[i * n + j] = pattern[perm[i] * n + perm[j]];
  ck_model *a, *b;
  CHECK(ck_linear_model_create(n, pattern, &a) == CK_OK);
  CHECK(ck_linear_model_create(n, ordered, &b) == CK_OK);
  CHECK(ck_jacobian_size(a) > ck_jacobian_size(b));
  ck_model_free(a);
  ck_model_free(b);
  return 0;
}
int test_linear(const char *suite, const char *test)
{
  int solve = !strcmp(suite, "LinearSolver") || !strcmp(suite, "LinearSolverInPlace");
  int lu = !strncmp(suite, "LuDecomposition", 15);
  if (!solve && !lu)
    return 2;
  int separate = strstr(suite, "InPlace") == NULL;
  ck_lu_algorithm algorithm = strstr(suite, "Doolittle") ? CK_LU_DOOLITTLE : CK_LU_MOZART;
  if (IS("DiagonalMarkowitzReorder") || IS("VectorDiagonalMarkowitzReordering"))
    return markowitz();
  int kind, n;
  if (!strncmp(test, "DenseMatrix", 11))
  {
    kind = 0;
    n = 3;
  }
  else if (!strncmp(test, "DiagonalMatrix", 14))
  {
    kind = 2;
    n = 6;
  }
  else if (!strncmp(test, "RandomMatrix", 12))
  {
    kind = 1;
    n = 10;
  }
  else if (strstr(test, "AgnosticToInitialValue"))
  {
    kind = 3;
    n = solve ? 30 : 10;
  }
  else
    return 2;
  int vector = strstr(test, "Vector") != NULL;
  size_t widths[] = {1, 2, 3, 4};
  double initial[] = {-INFINITY, -1, 0, 1, INFINITY};
  for (size_t v = 0; v < (vector ? 4u : 1u); ++v)
  {
    size_t width = vector ? widths[v] : 5;
    /* A vector ordering is represented by independent state groups, including a
       final partial group, instead of separate C++ matrix specializations. */
    size_t cells = kind == 0 ? 1 : 5;
    for (size_t init = 0; init < (kind == 3 ? 5u : 1u); ++init)
    {
      for (size_t first = 0; first < cells; first += width)
      {
        size_t count = cells - first;
        if (count > width)
          count = width;
        CHECK(numeric_case((size_t)n, count, kind == 3 ? 1 : kind, solve, separate, algorithm,
                           initial[init]) == 0);
      }
    }
  }
  if (lu && kind == 1 && !vector)
    CHECK(numeric_case(10, 1, 1, 0, separate, algorithm, 0) == 0);
  return 0;
}
