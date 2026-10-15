#include "test.h"
#include "chem_matrix.h"
#define D(m, r, c) (*(double *)ck_matrix_at((m), (r), (c)))
#define S(m, b, r, c) (*(double *)ck_sparse_at((m), (b), (r), (c)))
static ck_status arithmetic(ck_view *v, size_t n, size_t row, void *context)
{
  int op = *(int *)context;
  double *out = ck_view_at(v, row);
  if (!out)
    return CK_INVALID;
  double a = n > 1 ? *(const double *)ck_view_const_at(v + 1, row) : 0;
  double b = n > 2 ? *(const double *)ck_view_const_at(v + 2, row) : 0;
  *out = op == 0   ? a * 2
         : op == 1 ? a * 3 + 1
         : op == 2 ? a - b
         : op == 3 ? a + 5
         : op == 4 ? *out + a
         : op == 5 ? a + b
                   : 55.5;
  return CK_OK;
}
static ck_status sum_value(const void *v, void *total, void *context)
{
  (void)context;
  *(double *)total += *(const double *)v;
  return CK_OK;
}
static int dense_grouped(const char *test, size_t width)
{
  int roundtrip = IS("MatrixColumnViewRoundTrip") || IS("VectorMatrixColumnViewContiguousBlock");
  int plain = IS("PlainMatrixGroupedRoundTrip");
  size_t rows = roundtrip ? (width ? 8 : 3)
                : plain   ? 4
                          : 3 * width,
         cols = roundtrip ? (width ? 3 : 4)
                : plain   ? 3
                          : 4;
  double zero = 0;
  ck_matrix m;
  CHECK(ck_matrix_create(rows, cols, width, sizeof(double), &zero, &m) == CK_OK);
  for (size_t r = 0; r < rows; ++r)
    for (size_t c = 0; c < cols; ++c)
      D(&m, r, c) = plain ? r + 5 * c : (roundtrip && width ? 100 : 0) + 10 * r + c;
  size_t groups = roundtrip ? 1 : (rows / (width ? width : 1));
  int op = roundtrip ? 0 : plain ? 2 : 1;
  for (size_t group = 0; group < groups; ++group)
  {
    ck_view raw[3], v[3];
    size_t dest = roundtrip ? (width ? 1 : 2) : plain ? 2 : 3, src = roundtrip ? (width ? 2 : 3) : 0;
    CHECK(ck_matrix_column(&m, dest, 1, raw) == CK_OK);
    CHECK(ck_matrix_column(&m, src, 0, raw + 1) == CK_OK);
    if (plain)
      CHECK(ck_matrix_column(&m, 1, 0, raw + 2) == CK_OK);
    for (size_t j = 0; j < (plain ? 3u : 2u); ++j)
      CHECK(ck_view_group(raw[j], roundtrip ? 1 : group, v + j) == CK_OK);
    if (roundtrip)
      for (size_t i = 0; i < v[0].rows; ++i)
      {
        size_t r = (width ? width : 1) + i;
        CHECK(*(const double *)ck_view_const_at(v, i) == (width ? 100 : 0) + 10 * r + dest);
        CHECK(*(const double *)ck_view_const_at(v + 1, i) == (width ? 100 : 0) + 10 * r + src);
        CHECK(ck_view_at(v + 1, i) == NULL);
      }
    CHECK(ck_view_apply(v, plain ? 3 : 2, 1, arithmetic, &op) == CK_OK);
  }
  for (size_t r = 0; r < rows; ++r)
    for (size_t c = 0; c < cols; ++c)
    {
      double expected = plain ? r + 5 * c : (roundtrip && width ? 100 : 0) + 10 * r + c;
      if (roundtrip)
      {
        if (r >= (width ? width : 1) && r < 2 * (width ? width : 1) && c == (width ? 1u : 2u))
          expected = 2 * ((width ? 100 : 0) + 10 * r + (width ? 2 : 3));
      }
      else if (c == (plain ? 2u : 3u))
        expected = plain ? -5.0 : 30 * r + 1;
      CHECK(D(&m, r, c) == expected);
    }
  ck_matrix_free(&m);
  return 0;
}
static int sparse_grouped(const char *test, size_t width)
{
  ck_sparse_order order = strstr(test, "CSC") ? CK_SPARSE_COLUMN : CK_SPARSE_ROW;
  size_t L = width ? width : 1, blocks = 2 * L;
  double zero = 0;
  int mixed = strstr(test, "Mixed") != NULL, fill = strstr(test, "Fill") != NULL,
      constant = IS("ConstGroupViewProducesGroupedViews"),
      sum = IS("SparseVectorRowColOverloadStillFunctional");
  if (fill)
    blocks = L + 1;
  if (constant)
    blocks = 4;
  unsigned char pattern[] = {0, 1, 0, 0, 0, 1, 1, 0, 1};
  ck_sparse m;
  CHECK(ck_sparse_create(3, blocks, width, order, pattern, sizeof(double), &zero, &m) == CK_OK);
  for (size_t b = 0; b < blocks; ++b)
  {
    S(&m, b, 0, 1) = fill || mixed || constant || sum ? b + 1 : width ? 100 + b : b ? 10 : 1;
    S(&m, b, 1, 2) = width ? 200 + b : b ? 20 : 2;
    S(&m, b, 2, 0) = 100 + b;
    S(&m, b, 2, 2) = b ? 30 : 3;
  }
  ck_view raw[3], v[3];
  CHECK(ck_sparse_block_view(&m, m.map[1], 1, raw) == CK_OK);
  CHECK(ck_sparse_block_view(&m, m.map[5], 0, raw + 1) == CK_OK);
  if (mixed)
  {
    ck_matrix dense;
    CHECK(ck_matrix_create(blocks, 3, width, sizeof(double), &zero, &dense) == CK_OK);
    for (size_t b = 0; b < blocks; ++b)
      for (size_t c = 0; c < 3; ++c)
        D(&dense, b, c) = 10 + b + c;
    ck_view dv[2];
    CHECK(ck_matrix_column(&dense, 0, 1, dv) == CK_OK);
    CHECK(ck_matrix_column(&dense, 1, 0, dv + 1) == CK_OK);
    int op = 5;
    for (size_t g = 0; g < 2; ++g)
    {
      CHECK(ck_view_group(dv[0], g, v) == CK_OK);
      CHECK(ck_view_group(dv[1], g, v + 1) == CK_OK);
      CHECK(ck_view_group(raw[0], g, v + 2) == CK_OK);
      CHECK(ck_view_apply(v, 3, 1, arithmetic, &op) == CK_OK);
    }
    for (size_t b = 0; b < blocks; ++b)
      for (size_t c = 0; c < 3; ++c)
        CHECK(D(&dense, b, c) == (c ? 10 + b + c : 12 + 2 * b));
    ck_matrix_free(&dense);
  }
  else if (constant)
  {
    raw[0].writable = 0;
    CHECK(ck_view_group(raw[0], 0, v) == CK_OK);
    for (size_t i = 0; i < 4; ++i)
    {
      CHECK(*(const double *)ck_view_const_at(v, i) == i + 1);
      CHECK(!ck_view_at(v, i));
    }
    ck_matrix dense;
    CHECK(ck_matrix_create(4, 2, 4, sizeof(double), &zero, &dense) == CK_OK);
    for (size_t b = 0; b < 4; ++b)
      D(&dense, b, 1) = 100 + b;
    ck_view d, g;
    CHECK(ck_matrix_column(&dense, 1, 0, &d) == CK_OK);
    CHECK(ck_view_group(d, 0, &g) == CK_OK);
    for (size_t i = 0; i < 4; ++i)
    {
      CHECK(*(const double *)ck_view_const_at(&g, i) == 100 + i);
      CHECK(!ck_view_at(&g, i));
    }
    ck_matrix_free(&dense);
  }
  else if (sum)
  {
    CHECK(ck_view_group(raw[0], 1, v) == CK_OK);
    double result = 0;
    CHECK(ck_view_reduce(v[0], 1, sum_value, &result, NULL) == CK_OK);
    CHECK(result == 26);
  }
  else if (fill)
  {
    size_t group = !width && order == CK_SPARSE_COLUMN ? 1 : 0;
    double value = group ? 9 : 42;
    CHECK(ck_view_group(raw[0], group, v) == CK_OK);
    CHECK(ck_view_fill(v[0], &value, 0) == CK_OK);
    for (size_t b = 0; b < blocks; ++b)
    {
      CHECK(S(&m, b, 0, 1) == (b / L == group ? value : b + 1));
      CHECK(S(&m, b, 2, 0) == 100 + b);
    }
    ck_view dest;
    CHECK(ck_sparse_block_view(&m, m.map[6], 1, &dest) == CK_OK);
    CHECK(ck_view_group(dest, group, v + 1) == CK_OK);
    CHECK(ck_view_copy(v[1], v[0], 1) == CK_OK);
    for (size_t b = 0; b < blocks; ++b)
      CHECK(S(&m, b, 2, 0) == (b / L == group ? value : 100 + b));
    if (width)
    {
      CHECK(ck_view_group(raw[0], 1, v) == CK_OK);
      value = -7;
      CHECK(ck_view_fill(v[0], &value, 0) == CK_OK);
      CHECK(S(&m, L, 0, 1) == -7);
    }
  }
  else
  {
    CHECK(ck_view_group(raw[0], 1, v) == CK_OK);
    CHECK(ck_view_group(raw[1], 1, v + 1) == CK_OK);
    for (size_t i = 0; i < L; ++i)
    {
      CHECK(*(const double *)ck_view_const_at(v, i) == (width ? 100 + L + i : 10));
      CHECK(*(const double *)ck_view_const_at(v + 1, i) == (width ? 200 + L + i : 20));
    }
    int op = width ? 4 : 3;
    CHECK(ck_view_apply(v, 2, 1, arithmetic, &op) == CK_OK);
    for (size_t b = 0; b < blocks; ++b)
      CHECK(S(&m, b, 0, 1) == (b >= L ? (width ? 300 + 2 * b : 25) : (width ? 100 + b : 1)));
  }
  ck_sparse_free(&m);
  return 0;
}
static int cross_sources(void)
{
  ck_sparse src, dst, ref;
  unsigned char pattern[] = {0, 1, 0, 0};
  double zero = 0;
  CHECK(ck_sparse_create(2, 4, 4, CK_SPARSE_ROW, pattern, sizeof(double), &zero, &src) == CK_OK);
  CHECK(ck_sparse_create(2, 4, 4, CK_SPARSE_ROW, pattern, sizeof(double), &zero, &dst) == CK_OK);
  CHECK(ck_sparse_create(2, 4, 4, CK_SPARSE_ROW, pattern, sizeof(double), &zero, &ref) == CK_OK);
  for (size_t b = 0; b < 4; ++b)
  {
    S(&src, b, 0, 1) = b + 1;
    S(&ref, b, 0, 1) = 999;
  }
  ck_view sv, dv, rv, g[3];
  CHECK(ck_sparse_block_view(&src, 0, 0, &sv) == CK_OK);
  CHECK(ck_sparse_block_view(&dst, 0, 1, &dv) == CK_OK);
  CHECK(ck_sparse_block_view(&ref, 0, 1, &rv) == CK_OK);
  CHECK(ck_view_group(sv, 0, g) == CK_OK);
  CHECK(ck_view_group(dv, 0, g + 1) == CK_OK);
  CHECK(ck_view_group(rv, 0, g + 2) == CK_OK);
  CHECK(ck_view_copy(g[1], g[0], 1) == CK_OK);
  for (size_t b = 0; b < 4; ++b)
    CHECK(S(&dst, b, 0, 1) == b + 1);
  double value = 55.5;
  int op = 6;
  CHECK(ck_view_fill(g[1], &value, 0) == CK_OK);
  CHECK(ck_view_apply(g + 2, 1, 1, arithmetic, &op) == CK_OK);
  for (size_t b = 0; b < 4; ++b)
    CHECK(S(&dst, b, 0, 1) == S(&ref, b, 0, 1));
  ck_sparse_free(&src);
  ck_sparse_free(&dst);
  ck_sparse_free(&ref);
  return 0;
}
int test_grouped(const char *suite, const char *test)
{
  if (strcmp(suite, "GroupedView"))
    return 2;
  if (IS("CategoryTagsAreDistinct"))
    return 3;
  size_t width = 0;
  const char *l = strrchr(test, 'L');
  if (l)
    width = (size_t)strtoul(l + 1, NULL, 10);
  if (strstr(test, "Vector") || IS("ConstGroupViewProducesGroupedViews"))
    if (!width)
      width = 4;
  if (IS("FillCopyCrossGroupViewSources"))
    return cross_sources();
  if (strstr(test, "Sparse") || IS("ConstGroupViewProducesGroupedViews"))
    return sparse_grouped(test, width);
  return dense_grouped(test, width);
}
