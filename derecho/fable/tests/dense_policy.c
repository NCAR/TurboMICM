/* MICM test/unit/util/test_matrix_policy.hpp: the checks MICM runs on each
 * dense matrix policy, as functions of the width L (MICM's Matrix is L = 1,
 * its VectorMatrix<L> is L).
 *
 * MICM's Function lambdas become kernels: a GroupView's column views, row
 * variables, and vector views become the lane pointers of
 * tmc_array_group_column, _variable, and _vector, and each ForEachRow lambda a
 * lane function whose x[k] is the current lane of operand k. Where MICM
 * expects a throw, the C call returns TMC_ERROR_MATRIX with MICM's code. MICM's
 * int-valued matrices and vectors hold the same integers in doubles here
 * (tmc_dense stores doubles; every value involved is exact). */
#include "dense_policy.h"

#include "tmc_dense_ext.h"

#include <float.h>

#define N_ARGS(a) ((int)(sizeof(a) / sizeof((a)[0])))
#define AT(m, i, j) tmc_dense_get((m), (i), (j))

#define EXPECT_MATRIX_ERROR(rc, err, expected_code)                                           \
  do                                                                                          \
  {                                                                                           \
    int rc_ = (rc);                                                                           \
    TMC_EXPECT_EQ(rc_, TMC_ERROR_MATRIX);                                                     \
    TMC_EXPECT_EQ((err).category, TMC_ERROR_MATRIX);                                          \
    TMC_EXPECT_EQ((err).code, (expected_code));                                               \
  } while (0)

/* f applied to args succeeds (MICM: no throw). */
#define EXPECT_APPLY_OK(f, args)                                                              \
  do                                                                                          \
  {                                                                                           \
    tmc_error e_;                                                                             \
    tmc_error_clear(&e_);                                                                     \
    if (tmc_array_function_apply((f), N_ARGS(args), (args), &e_) != TMC_OK)                   \
      tmc_test_fail(__FILE__, __LINE__, "L = %d: apply failed: %s", L, e_.message);           \
  } while (0)

/* f applied to args fails with the matrix error code (MICM: EXPECT_ANY_THROW). */
#define EXPECT_APPLY_ERROR(f, args, expected_code)                                            \
  do                                                                                          \
  {                                                                                           \
    tmc_error e_;                                                                             \
    tmc_error_clear(&e_);                                                                     \
    EXPECT_MATRIX_ERROR(tmc_array_function_apply((f), N_ARGS(args), (args), &e_), e_, (expected_code)); \
  } while (0)

static tmc_array_function *make_function(tmc_array_kernel kernel, void *ctx, int n_args, const tmc_array_arg *args, int L)
{
  tmc_error err;
  tmc_error_clear(&err);
  tmc_array_function *f = tmc_array_function_create(kernel, ctx, n_args, args, &err);
  if (!f)
    tmc_test_fail(__FILE__, __LINE__, "L = %d: function creation failed: %s", L, err.message);
  return f;
}

/* Lane functions: x[k] is the current lane of operand k. */
static void sum3_into_last(double *const *x, void *ctx) /* t = a + b + c */
{
  (void)ctx;
  *x[3] = *x[0] + *x[1] + *x[2];
}

static void sum2_into_last(double *const *x, void *ctx) /* c = a + b */
{
  (void)ctx;
  *x[2] = *x[0] + *x[1];
}

static void copy_first_to_second(double *const *x, void *ctx) /* b = a */
{
  (void)ctx;
  *x[1] = *x[0];
}

static void copy_second_to_first(double *const *x, void *ctx) /* a = b */
{
  (void)ctx;
  *x[0] = *x[1];
}

static void twice_first_into_second(double *const *x, void *ctx) /* b = a * 2 */
{
  (void)ctx;
  *x[1] = *x[0] * 2.0;
}

/* ---- construction and element access ------------------------------------ */

tmc_dense *dense_policy_small_matrix(int L)
{
  tmc_dense *m = tmc_dense_create(3, 5, L, 0.0);
  tmc_dense_set(m, 1, 3, 64.7);
  tmc_dense_set(m, 0, 0, 41.2);
  tmc_dense_set(m, 2, 4, 102.3);
  TMC_EXPECT_EQ(AT(m, 1, 3), 64.7);
  TMC_EXPECT_EQ(AT(m, 0, 0), 41.2);
  TMC_EXPECT_EQ(AT(m, 2, 4), 102.3);
  TMC_EXPECT_TRUE(m->size >= 3);
  return m;
}

tmc_dense *dense_policy_small_const_matrix(int L)
{
  tmc_dense *m = tmc_dense_create(3, 5, L, 0.0);
  tmc_dense_set(m, 1, 3, 64.7);
  tmc_dense_set(m, 0, 0, 41.2);
  tmc_dense_set(m, 2, 4, 102.3);
  tmc_dense *copy = tmc_dense_copy(m);
  tmc_dense_free(m);
  const tmc_dense *const_matrix = copy;
  TMC_EXPECT_EQ(AT(const_matrix, 1, 3), 64.7);
  TMC_EXPECT_EQ(AT(const_matrix, 0, 0), 41.2);
  TMC_EXPECT_EQ(AT(const_matrix, 2, 4), 102.3);
  TMC_EXPECT_TRUE(const_matrix->size >= 3);
  return copy;
}

int dense_policy_initialize_matrix(int L)
{
  tmc_dense *m = tmc_dense_create(2, 3, L, 12.4);
  TMC_EXPECT_EQ(AT(m, 0, 0), 12.4);
  TMC_EXPECT_EQ(AT(m, 1, 0), 12.4);
  TMC_EXPECT_EQ(AT(m, 1, 2), 12.4);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

int dense_policy_initialize_const_matrix(int L)
{
  tmc_dense *storage = tmc_dense_create(2, 3, L, 12.4);
  const tmc_dense *m = storage;
  TMC_EXPECT_EQ(AT(m, 0, 0), 12.4);
  TMC_EXPECT_EQ(AT(m, 1, 0), 12.4);
  TMC_EXPECT_EQ(AT(m, 1, 2), 12.4);
  tmc_dense_free(storage);
  return TMC_TEST_PASS;
}

static tmc_dense *looped_matrix(int L)
{
  tmc_dense *m = tmc_dense_create(3, 4, L, 0.0);
  for (int i = 0; i < m->n_rows; ++i)
    for (int j = 0; j < m->n_cols; ++j)
      tmc_dense_set(m, i, j, i * 100 + j);
  return m;
}

int dense_policy_loop_over_matrix(int L)
{
  tmc_dense *m = looped_matrix(L);
  TMC_EXPECT_EQ(AT(m, 0, 0), 0);
  TMC_EXPECT_EQ(AT(m, 1, 2), 102);
  TMC_EXPECT_EQ(AT(m, 2, 3), 203);
  TMC_EXPECT_EQ(AT(m, 0, 3), 3);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

int dense_policy_loop_over_const_matrix(int L)
{
  tmc_dense *m = looped_matrix(L);
  tmc_dense *copy = tmc_dense_copy(m);
  const tmc_dense *const_matrix = copy;
  TMC_EXPECT_EQ(AT(const_matrix, 0, 0), 0);
  TMC_EXPECT_EQ(AT(const_matrix, 1, 2), 102);
  TMC_EXPECT_EQ(AT(const_matrix, 2, 3), 203);
  TMC_EXPECT_EQ(AT(const_matrix, 0, 3), 3);
  tmc_dense_free(m);
  tmc_dense_free(copy);
  return TMC_TEST_PASS;
}

tmc_dense *dense_policy_strides(int L)
{
  tmc_dense *m = tmc_dense_create(3, 4, L, 0.0);
  size_t row_stride = tmc_dense_row_stride(m), column_stride = tmc_dense_column_stride(m);
  for (int i = 0; i < m->n_rows; ++i)
    for (int j = 0; j < m->n_cols; ++j)
      m->data[(size_t)i * row_stride + (size_t)j * column_stride] = i * 100 + j;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 4; ++j)
      TMC_EXPECT_EQ(AT(m, i, j), i * 100 + j);
  return m;
}

/* ---- rows as vectors ------------------------------------------------------ */

static tmc_dense *row_one_matrix(int L)
{
  tmc_dense *m = tmc_dense_create(2, 3, L, 0.0);
  tmc_dense_set(m, 1, 0, 13.2);
  tmc_dense_set(m, 1, 1, 31.2);
  tmc_dense_set(m, 1, 2, 314.2);
  return m;
}

int dense_policy_conversion_to_vector(int L)
{
  tmc_dense *m = row_one_matrix(L);
  double slice[3];
  tmc_dense_get_row(m, 1, slice);
  TMC_EXPECT_EQ(slice[0], 13.2);
  TMC_EXPECT_EQ(slice[1], 31.2);
  TMC_EXPECT_EQ(slice[2], 314.2);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

int dense_policy_const_conversion_to_vector(int L)
{
  tmc_dense *m = row_one_matrix(L);
  tmc_dense *copy = tmc_dense_copy(m);
  const tmc_dense *const_matrix = copy;
  double slice[3];
  tmc_dense_get_row(const_matrix, 1, slice);
  TMC_EXPECT_EQ(slice[0], 13.2);
  TMC_EXPECT_EQ(slice[1], 31.2);
  TMC_EXPECT_EQ(slice[2], 314.2);
  tmc_dense_free(m);
  tmc_dense_free(copy);
  return TMC_TEST_PASS;
}

int dense_policy_conversion_from_vector(int L)
{
  tmc_error err;
  tmc_error_clear(&err);
  tmc_dense *zero_matrix = tmc_dense_create_from_rows(0, NULL, NULL, L, &err);
  TMC_ASSERT_TRUE(zero_matrix != NULL);
  TMC_EXPECT_EQ(zero_matrix->n_rows, 0);

  const double row0[] = { 412.3, 32.4, 41.3 }, row1[] = { 5.33, -0.3, 31.2 };
  const double *rows[] = { row0, row1 };
  const int lengths[] = { 3, 3 };
  tmc_dense *m = tmc_dense_create_from_rows(2, rows, lengths, L, &err);
  TMC_ASSERT_TRUE(m != NULL);
  TMC_EXPECT_EQ(m->n_rows, 2);
  TMC_EXPECT_EQ(m->n_cols, 3); /* each row's Size() */
  TMC_EXPECT_EQ(AT(m, 0, 0), 412.3);
  TMC_EXPECT_EQ(AT(m, 0, 1), 32.4);
  TMC_EXPECT_EQ(AT(m, 0, 2), 41.3);
  TMC_EXPECT_EQ(AT(m, 1, 0), 5.33);
  TMC_EXPECT_EQ(AT(m, 1, 1), -0.3);
  TMC_EXPECT_EQ(AT(m, 1, 2), 31.2);

  /* rows of different lengths: { { 3 }, { 4, 5 }, { 5 } } */
  const double bad0[] = { 3 }, bad1[] = { 4, 5 }, bad2[] = { 5 };
  const double *bad_rows[] = { bad0, bad1, bad2 };
  const int bad_lengths[] = { 1, 2, 1 };
  tmc_dense *bad_matrix = tmc_dense_create_from_rows(3, bad_rows, bad_lengths, L, &err);
  TMC_EXPECT_TRUE(bad_matrix == NULL);
  TMC_EXPECT_EQ(err.category, TMC_ERROR_MATRIX);
  TMC_EXPECT_EQ(err.code, TMC_MATRIX_INVALID_VECTOR);

  tmc_dense_free(zero_matrix);
  tmc_dense_free(m);
  tmc_dense_free(bad_matrix);
  return TMC_TEST_PASS;
}

int dense_policy_assignment_from_vector(int L)
{
  tmc_dense *m = tmc_dense_create(4, 3, L, 0.0);
  const double other[] = { 12.3, 15.1, 24.3 };
  const double big_other[] = { 14.3, 52.3, 65.7, 16.34 };
  const double small_other[] = { 13.2, 52.8 };
  tmc_error err;
  tmc_error_clear(&err);

  TMC_EXPECT_EQ(tmc_dense_assign_row(m, 2, other, 3, &err), TMC_OK);
  TMC_EXPECT_EQ(AT(m, 0, 0), 0.0);
  TMC_EXPECT_EQ(AT(m, 2, 0), 12.3);
  TMC_EXPECT_EQ(AT(m, 2, 1), 15.1);
  TMC_EXPECT_EQ(AT(m, 2, 2), 24.3);
  TMC_EXPECT_EQ(AT(m, 3, 0), 0.0);

  TMC_EXPECT_EQ(tmc_dense_assign_row(m, 2, big_other, 4, &err), TMC_OK);
  TMC_EXPECT_EQ(AT(m, 0, 0), 0.0);
  TMC_EXPECT_EQ(AT(m, 2, 0), 14.3);
  TMC_EXPECT_EQ(AT(m, 2, 1), 52.3);
  TMC_EXPECT_EQ(AT(m, 2, 2), 65.7);
  TMC_EXPECT_EQ(AT(m, 3, 0), 0.0);

  EXPECT_MATRIX_ERROR(tmc_dense_assign_row(m, 2, small_other, 2, &err), err, TMC_MATRIX_ROW_SIZE_MISMATCH);
  /* the failed assignment leaves the row as it was */
  TMC_EXPECT_EQ(AT(m, 2, 0), 14.3);
  TMC_EXPECT_EQ(AT(m, 2, 1), 52.3);
  TMC_EXPECT_EQ(AT(m, 2, 2), 65.7);

  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

/* ---- whole-matrix operations ---------------------------------------------- */

int dense_policy_axpy(int L)
{
  const int num_rows = 4, num_columns = 3;
  tmc_dense *m = tmc_dense_create(num_rows, num_columns, L, 100.0);
  tmc_dense *other = tmc_dense_create(num_rows, num_columns, L, 200.0);
  const double alpha = 1.39;
  double sum = 0.0, result = 0.0;
  for (int i = 0; i < num_rows; ++i)
    for (int j = 0; j < num_columns; ++j)
    {
      double y = i * 10.3 + j * 100.5;
      double x = i * 1.7 + j * 10.2;
      tmc_dense_set(m, i, j, y);
      tmc_dense_set(other, i, j, x);
      sum += y + alpha * x;
    }
  tmc_dense_axpy(m, alpha, other);
  for (int i = 0; i < num_rows; ++i)
    for (int j = 0; j < num_columns; ++j)
      result += AT(m, i, j);
  TMC_EXPECT_NEAR(sum, result, 1.0e-5);
  tmc_dense_free(m);
  tmc_dense_free(other);
  return TMC_TEST_PASS;
}

typedef struct for_each_count
{
  int calls;
} for_each_count;

static void add_into(double *y, const double *a, void *ctx) /* a += b */
{
  ++((for_each_count *)ctx)->calls;
  *y += *a;
}

static void add_subtract_into(double *y, const double *a, const double *b, void *ctx) /* a = a + b - c */
{
  ++((for_each_count *)ctx)->calls;
  *y = *y + *a - *b;
}

int dense_policy_for_each(int L)
{
  tmc_dense *m = tmc_dense_create(4, 3, L, 0.0);
  tmc_dense *other = tmc_dense_create(4, 3, L, 0.0);
  tmc_dense *other2 = tmc_dense_create(4, 3, L, 0.0);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 3; ++j)
    {
      tmc_dense_set(m, i, j, i * 10.3 + j * 100.5);
      tmc_dense_set(other, i, j, i * 1.7 + j * 10.2);
      tmc_dense_set(other2, i, j, i * 19.5 + j * 32.2);
    }

  for_each_count count = { 0 };
  tmc_dense_for_each(m, add_into, &count, other);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 3; ++j)
      TMC_EXPECT_NEAR(AT(m, i, j), (i * 10.3 + j * 100.5) + (i * 1.7 + j * 10.2), 1.0e-5);
  TMC_EXPECT_EQ(count.calls, 4 * 3); /* once per element, padding skipped */

  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 3; ++j)
      tmc_dense_set(m, i, j, i * 10.3 + j * 100.5);

  count.calls = 0;
  tmc_dense_for_each2(m, add_subtract_into, &count, other, other2);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 3; ++j)
      TMC_EXPECT_NEAR(AT(m, i, j), (i * 10.3 + j * 100.5) + (i * 1.7 + j * 10.2) - (i * 19.5 + j * 32.2), 1.0e-5);
  TMC_EXPECT_EQ(count.calls, 4 * 3);

  tmc_dense_free(m);
  tmc_dense_free(other);
  tmc_dense_free(other2);
  return TMC_TEST_PASS;
}

int dense_policy_set_scalar(int L)
{
  tmc_dense *m = tmc_dense_create(2, 3, L, 0.0);
  tmc_dense_fill(m, 2.0);
  for (size_t k = 0; k < m->size; ++k)
    TMC_EXPECT_EQ(m->data[k], 2.0);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

int dense_policy_max(int L)
{
  tmc_dense *m = tmc_dense_create(2, 3, L, 0.0);
  tmc_dense_max(m, 2.0);
  for (size_t k = 0; k < m->size; ++k)
    TMC_EXPECT_EQ(m->data[k], 2.0);
  tmc_dense_fill(m, 1.0);
  tmc_dense_set(m, 1, 1, 3.0);
  tmc_dense_max(m, 2.0);
  TMC_EXPECT_EQ(AT(m, 0, 0), 2.0);
  TMC_EXPECT_EQ(AT(m, 0, 1), 2.0);
  TMC_EXPECT_EQ(AT(m, 0, 2), 2.0);
  TMC_EXPECT_EQ(AT(m, 1, 0), 2.0);
  TMC_EXPECT_EQ(AT(m, 1, 1), 3.0);
  TMC_EXPECT_EQ(AT(m, 1, 2), 2.0);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

int dense_policy_min(int L)
{
  tmc_dense *m = tmc_dense_create(2, 3, L, 0.0);
  tmc_dense_min(m, 2.0);
  for (size_t k = 0; k < m->size; ++k)
    TMC_EXPECT_EQ(m->data[k], 0.0);
  tmc_dense_fill(m, 1.0);
  tmc_dense_set(m, 1, 1, 3.0);
  tmc_dense_min(m, 2.0);
  TMC_EXPECT_EQ(AT(m, 0, 0), 1.0);
  TMC_EXPECT_EQ(AT(m, 0, 1), 1.0);
  TMC_EXPECT_EQ(AT(m, 0, 2), 1.0);
  TMC_EXPECT_EQ(AT(m, 1, 0), 1.0);
  TMC_EXPECT_EQ(AT(m, 1, 1), 2.0);
  TMC_EXPECT_EQ(AT(m, 1, 2), 1.0);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

int dense_policy_print(int L)
{
  tmc_dense *m = tmc_dense_create(2, 3, L, 0.0);
  tmc_dense_set(m, 1, 1, 3.0);
  const char *expected = "0,0,0\n0,3,0\n";

  char text[64];
  size_t length = tmc_dense_format(m, text, sizeof text);
  TMC_EXPECT_STREQ(text, expected);
  TMC_EXPECT_EQ(length, strlen(expected));

  /* a short buffer gets the start of the text; the full length is returned */
  char short_text[5];
  TMC_EXPECT_EQ(tmc_dense_format(m, short_text, sizeof short_text), strlen(expected));
  TMC_EXPECT_STREQ(short_text, "0,0,");

  FILE *stream = tmpfile();
  TMC_ASSERT_TRUE(stream != NULL);
  TMC_EXPECT_EQ(tmc_dense_print(m, stream), (int)strlen(expected));
  rewind(stream);
  char printed[64] = { 0 };
  size_t n = fread(printed, 1, sizeof printed - 1, stream);
  fclose(stream);
  TMC_EXPECT_EQ(n, strlen(expected));
  TMC_EXPECT_STREQ(printed, expected);

  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

/* ---- Function ------------------------------------------------------------- */

static void four_times_second_into_first(double *const *x, void *ctx) /* c = 4 t */
{
  (void)ctx;
  *x[0] = 4.0 * *x[1];
}

/* tmp = col0 + col1 + col2; col2 = 4 tmp */
static void array_function_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, 1), tmc_array_group_column(g, 0, 2), tmp };
  tmc_array_group_for_each(g, sum3_into_last, NULL, 4, sum);
  double *scale[] = { tmc_array_group_column(g, 0, 2), tmp };
  tmc_array_group_for_each(g, four_times_second_into_first, NULL, 2, scale);
}

int dense_policy_array_function(int L)
{
  tmc_dense *m = tmc_dense_create(5, 3, L, -1.0);
  for (int i = 0; i < m->n_rows; ++i)
    for (int j = 0; j < m->n_cols; ++j)
      tmc_dense_set(m, i, j, i - 2 + 10 * j);
  tmc_array_arg args[] = { tmc_array_dense(m) };
  tmc_array_function *f = make_function(array_function_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 2), 4.0 * (-2 + 8 + 18));
  TMC_EXPECT_EQ(AT(m, 1, 2), 4.0 * (-1 + 9 + 19));
  TMC_EXPECT_EQ(AT(m, 2, 2), 4.0 * (0 + 10 + 20));
  TMC_EXPECT_EQ(AT(m, 3, 2), 4.0 * (1 + 11 + 21));
  TMC_EXPECT_EQ(AT(m, 4, 2), 4.0 * (2 + 12 + 22));
  TMC_EXPECT_EQ(AT(m, 0, 0), -2.0);
  TMC_EXPECT_EQ(AT(m, 1, 0), -1.0);
  TMC_EXPECT_EQ(AT(m, 2, 0), 0.0);
  TMC_EXPECT_EQ(AT(m, 3, 0), 1.0);
  TMC_EXPECT_EQ(AT(m, 4, 0), 2.0);
  TMC_EXPECT_EQ(AT(m, 0, 1), 8.0);
  TMC_EXPECT_EQ(AT(m, 1, 1), 9.0);
  TMC_EXPECT_EQ(AT(m, 2, 1), 10.0);
  TMC_EXPECT_EQ(AT(m, 3, 1), 11.0);
  TMC_EXPECT_EQ(AT(m, 4, 1), 12.0);

  /* the same function on a matrix with the same columns and fewer rows */
  tmc_dense *m2 = tmc_dense_create(3, 3, L, -1.0);
  tmc_array_arg args2[] = { tmc_array_dense(m2) };
  EXPECT_APPLY_OK(f, args2);
  for (int i = 0; i < 3; ++i)
  {
    TMC_EXPECT_EQ(AT(m2, i, 2), 4.0 * (-1 + -1 + -1));
    TMC_EXPECT_EQ(AT(m2, i, 0), -1.0);
    TMC_EXPECT_EQ(AT(m2, i, 1), -1.0);
  }

  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_dense_free(m2);
  return TMC_TEST_PASS;
}

/* tmp = A[0] + B[2]; A[1] = tmp, with A argument 0 and B argument 1 */
static void multi_matrix_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 1, 2), tmp };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
  double *copy[] = { tmp, tmc_array_group_column(g, 0, 1) };
  tmc_array_group_for_each(g, copy_first_to_second, NULL, 2, copy);
}

int dense_policy_multi_matrix_array_function(int L)
{
  tmc_dense *a = tmc_dense_create(3, 2, L, 1.0);
  tmc_dense *b = tmc_dense_create(3, 3, L, 2.0);
  for (int i = 0; i < a->n_rows; ++i)
  {
    for (int j = 0; j < a->n_cols; ++j)
      tmc_dense_set(a, i, j, i + 10 * j);
    for (int j = 0; j < b->n_cols; ++j)
      tmc_dense_set(b, i, j, i * 2 + 20 * j);
  }
  for (int i = 0; i < b->n_rows; ++i)
    tmc_dense_set(b, i, 2, i * 4);

  tmc_array_arg args[] = { tmc_array_dense(a), tmc_array_const_dense(b) };
  tmc_array_function *f = make_function(multi_matrix_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(a, 0, 1), 0 + 0);
  TMC_EXPECT_EQ(AT(a, 1, 1), 1 + 4);
  TMC_EXPECT_EQ(AT(a, 2, 1), 2 + 8);
  TMC_EXPECT_EQ(AT(a, 0, 0), 0.0);
  TMC_EXPECT_EQ(AT(a, 1, 0), 1.0);
  TMC_EXPECT_EQ(AT(a, 2, 0), 2.0);
  TMC_EXPECT_EQ(AT(b, 0, 0), 0.0);
  TMC_EXPECT_EQ(AT(b, 1, 0), 2.0);
  TMC_EXPECT_EQ(AT(b, 2, 0), 4.0);
  tmc_array_function_free(f);
  tmc_dense_free(a);
  tmc_dense_free(b);
  return TMC_TEST_PASS;
}

static void product_of_last_two_into_first(double *const *x, void *ctx) /* c = d * t */
{
  (void)ctx;
  *x[0] = *x[1] * *x[2];
}

/* tmp = col0 + col1 + col2; col2 = v * tmp, with the vector argument 1 */
static void vector_in_matrix_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, 1), tmc_array_group_column(g, 0, 2), tmp };
  tmc_array_group_for_each(g, sum3_into_last, NULL, 4, sum);
  double *scale[] = { tmc_array_group_column(g, 0, 2), tmc_array_group_vector(g, 1), tmp };
  tmc_array_group_for_each(g, product_of_last_two_into_first, NULL, 3, scale);
}

int dense_policy_vector_in_matrix_function(int L)
{
  tmc_dense *m = tmc_dense_create(5, 3, L, -1.0);
  for (int i = 0; i < m->n_rows; ++i)
    for (int j = 0; j < m->n_cols; ++j)
      tmc_dense_set(m, i, j, i - 2 + 10 * j);
  tmc_padded_vector *v = tmc_padded_vector_create(m->n_rows, L, 0.0);
  v->data[0] = 100.0;
  v->data[1] = 200.0;
  v->data[2] = 300.0;
  v->data[3] = 400.0;
  v->data[4] = 500.0;

  tmc_array_arg args[] = { tmc_array_dense(m), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(vector_in_matrix_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 2), 100.0 * (-2 + 8 + 18));
  TMC_EXPECT_EQ(AT(m, 1, 2), 200.0 * (-1 + 9 + 19));
  TMC_EXPECT_EQ(AT(m, 2, 2), 300.0 * (0 + 10 + 20));
  TMC_EXPECT_EQ(AT(m, 3, 2), 400.0 * (1 + 11 + 21));
  TMC_EXPECT_EQ(AT(m, 4, 2), 500.0 * (2 + 12 + 22));
  TMC_EXPECT_EQ(AT(m, 0, 0), -2.0);
  TMC_EXPECT_EQ(AT(m, 1, 0), -1.0);
  TMC_EXPECT_EQ(AT(m, 2, 0), 0.0);
  TMC_EXPECT_EQ(AT(m, 3, 0), 1.0);
  TMC_EXPECT_EQ(AT(m, 4, 0), 2.0);
  TMC_EXPECT_EQ(AT(m, 0, 1), 8.0);
  TMC_EXPECT_EQ(AT(m, 1, 1), 9.0);
  TMC_EXPECT_EQ(AT(m, 2, 1), 10.0);
  TMC_EXPECT_EQ(AT(m, 3, 1), 11.0);
  TMC_EXPECT_EQ(AT(m, 4, 1), 12.0);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

/* ---- rows at creation and at invocation ------------------------------------ */

/* tmp = A[0] + B[2]; A[1] = 2 tmp */
static void different_rows_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 1, 2), tmp };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
  double *twice[] = { tmp, tmc_array_group_column(g, 0, 1) };
  tmc_array_group_for_each(g, twice_first_into_second, NULL, 2, twice);
}

int dense_policy_multi_matrix_different_rows_from_creation(int L)
{
  tmc_dense *a_create = tmc_dense_create(3, 2, L, 0.0);
  tmc_dense *b_create = tmc_dense_create(3, 3, L, 0.0);
  tmc_array_arg create_args[] = { tmc_array_dense(a_create), tmc_array_const_dense(b_create) };
  tmc_array_function *f = make_function(different_rows_kernel, NULL, N_ARGS(create_args), create_args, L);
  if (!f)
    return TMC_TEST_FAIL;

  /* more rows than at creation */
  tmc_dense *a = tmc_dense_create(5, 2, L, 0.0);
  tmc_dense *b = tmc_dense_create(5, 3, L, 0.0);
  for (int i = 0; i < 5; ++i)
  {
    tmc_dense_set(a, i, 0, i + 1);
    tmc_dense_set(b, i, 2, i * 10);
  }
  tmc_array_arg args[] = { tmc_array_dense(a), tmc_array_const_dense(b) };
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(a, 0, 1), (1.0 + 0.0) * 2.0);
  TMC_EXPECT_EQ(AT(a, 1, 1), (2.0 + 10.0) * 2.0);
  TMC_EXPECT_EQ(AT(a, 2, 1), (3.0 + 20.0) * 2.0);
  TMC_EXPECT_EQ(AT(a, 3, 1), (4.0 + 30.0) * 2.0);
  TMC_EXPECT_EQ(AT(a, 4, 1), (5.0 + 40.0) * 2.0);

  /* fewer rows than at creation */
  tmc_dense *a2 = tmc_dense_create(2, 2, L, 0.0);
  tmc_dense *b2 = tmc_dense_create(2, 3, L, 0.0);
  tmc_dense_set(a2, 0, 0, 10.0);
  tmc_dense_set(a2, 1, 0, 20.0);
  tmc_dense_set(b2, 0, 2, 5.0);
  tmc_dense_set(b2, 1, 2, 15.0);
  tmc_array_arg args2[] = { tmc_array_dense(a2), tmc_array_const_dense(b2) };
  EXPECT_APPLY_OK(f, args2);
  TMC_EXPECT_EQ(AT(a2, 0, 1), (10.0 + 5.0) * 2.0);
  TMC_EXPECT_EQ(AT(a2, 1, 1), (20.0 + 15.0) * 2.0);

  tmc_array_function_free(f);
  tmc_dense_free(a_create);
  tmc_dense_free(b_create);
  tmc_dense_free(a);
  tmc_dense_free(b);
  tmc_dense_free(a2);
  tmc_dense_free(b2);
  return TMC_TEST_PASS;
}

static void three_times_first_into_second(double *const *x, void *ctx) /* c = 3 t */
{
  (void)ctx;
  *x[1] = *x[0] * 3.0;
}

/* tmp = col0 + v; col1 = 3 tmp */
static void matrix_vector_rows_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_vector(g, 1), tmp };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
  double *thrice[] = { tmp, tmc_array_group_column(g, 0, 1) };
  tmc_array_group_for_each(g, three_times_first_into_second, NULL, 2, thrice);
}

int dense_policy_matrix_vector_different_rows_from_creation(int L)
{
  tmc_dense *m_create = tmc_dense_create(3, 3, L, 0.0);
  tmc_padded_vector *v_create = tmc_padded_vector_create(3, L, 0.0);
  tmc_array_arg create_args[] = { tmc_array_dense(m_create), tmc_array_const_vector(v_create) };
  tmc_array_function *f = make_function(matrix_vector_rows_kernel, NULL, N_ARGS(create_args), create_args, L);
  if (!f)
    return TMC_TEST_FAIL;

  tmc_dense *m = tmc_dense_create(5, 3, L, 0.0);
  tmc_padded_vector *v = tmc_padded_vector_create(5, L, 0.0);
  for (int i = 0; i < 5; ++i)
  {
    tmc_dense_set(m, i, 0, i + 1);
    v->data[i] = i * 10;
  }
  tmc_array_arg args[] = { tmc_array_dense(m), tmc_array_const_vector(v) };
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 1), (1.0 + 0.0) * 3.0);
  TMC_EXPECT_EQ(AT(m, 1, 1), (2.0 + 10.0) * 3.0);
  TMC_EXPECT_EQ(AT(m, 2, 1), (3.0 + 20.0) * 3.0);
  TMC_EXPECT_EQ(AT(m, 3, 1), (4.0 + 30.0) * 3.0);
  TMC_EXPECT_EQ(AT(m, 4, 1), (5.0 + 40.0) * 3.0);

  tmc_array_function_free(f);
  tmc_dense_free(m_create);
  tmc_padded_vector_free(v_create);
  tmc_dense_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

/* col0 = 2 v */
static void vector_doubled_into_column_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *ops[] = { tmc_array_group_vector(g, 1), tmc_array_group_column(g, 0, 0) };
  tmc_array_group_for_each(g, twice_first_into_second, NULL, 2, ops);
}

int dense_policy_mismatched_rows_at_invocation(int L)
{
  tmc_dense *m_create = tmc_dense_create(3, 2, L, 0.0);
  tmc_padded_vector *v_create = tmc_padded_vector_create(3, L, 0.0);
  tmc_array_arg create_args[] = { tmc_array_dense(m_create), tmc_array_const_vector(v_create) };
  tmc_array_function *f = make_function(vector_doubled_into_column_kernel, NULL, N_ARGS(create_args), create_args, L);
  if (!f)
    return TMC_TEST_FAIL;

  /* matrix of 5 rows, vector of 3 */
  tmc_dense *m = tmc_dense_create(5, 2, L, 0.0);
  tmc_padded_vector *v = tmc_padded_vector_create(3, L, 0.0);
  tmc_array_arg args[] = { tmc_array_dense(m), tmc_array_const_vector(v) };
  EXPECT_APPLY_ERROR(f, args, TMC_MATRIX_INVALID_VECTOR);

  /* matrix of 3 rows, vector of 5 */
  tmc_dense *m2 = tmc_dense_create(3, 2, L, 0.0);
  tmc_padded_vector *v2 = tmc_padded_vector_create(5, L, 0.0);
  tmc_array_arg args2[] = { tmc_array_dense(m2), tmc_array_const_vector(v2) };
  EXPECT_APPLY_ERROR(f, args2, TMC_MATRIX_INVALID_VECTOR);

  tmc_array_function_free(f);
  tmc_dense_free(m_create);
  tmc_padded_vector_free(v_create);
  tmc_dense_free(m);
  tmc_padded_vector_free(v);
  tmc_dense_free(m2);
  tmc_padded_vector_free(v2);
  return TMC_TEST_PASS;
}

/* A[0] = 2 B[0] */
static void b_doubled_into_a_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *ops[] = { tmc_array_group_column(g, 1, 0), tmc_array_group_column(g, 0, 0) };
  tmc_array_group_for_each(g, twice_first_into_second, NULL, 2, ops);
}

int dense_policy_multiple_matrices_mismatched_rows_at_invocation(int L)
{
  tmc_dense *a_create = tmc_dense_create(3, 2, L, 0.0);
  tmc_dense *b_create = tmc_dense_create(3, 3, L, 0.0);
  tmc_array_arg create_args[] = { tmc_array_dense(a_create), tmc_array_const_dense(b_create) };
  tmc_array_function *f = make_function(b_doubled_into_a_kernel, NULL, N_ARGS(create_args), create_args, L);
  if (!f)
    return TMC_TEST_FAIL;

  tmc_dense *a = tmc_dense_create(5, 2, L, 0.0);
  tmc_dense *b = tmc_dense_create(3, 3, L, 0.0);
  tmc_array_arg args[] = { tmc_array_dense(a), tmc_array_const_dense(b) };
  EXPECT_APPLY_ERROR(f, args, TMC_MATRIX_INVALID_VECTOR);

  tmc_array_function_free(f);
  tmc_dense_free(a_create);
  tmc_dense_free(b_create);
  tmc_dense_free(a);
  tmc_dense_free(b);
  return TMC_TEST_PASS;
}

/* tmp = col0 + col1 + col2 */
static void sum_of_three_columns_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, 1), tmc_array_group_column(g, 0, 2), tmp };
  tmc_array_group_for_each(g, sum3_into_last, NULL, 4, sum);
}

int dense_policy_wrong_column_count_at_invocation(int L)
{
  tmc_dense *m_create = tmc_dense_create(4, 3, L, 0.0);
  tmc_array_arg create_args[] = { tmc_array_dense(m_create) };
  tmc_array_function *f = make_function(sum_of_three_columns_kernel, NULL, N_ARGS(create_args), create_args, L);
  if (!f)
    return TMC_TEST_FAIL;

  tmc_dense *wrong_cols = tmc_dense_create(4, 4, L, 0.0);
  tmc_array_arg args[] = { tmc_array_dense(wrong_cols) };
  EXPECT_APPLY_ERROR(f, args, TMC_MATRIX_INVALID_VECTOR);

  tmc_dense *wrong_cols2 = tmc_dense_create(4, 2, L, 0.0);
  tmc_array_arg args2[] = { tmc_array_dense(wrong_cols2) };
  EXPECT_APPLY_ERROR(f, args2, TMC_MATRIX_INVALID_VECTOR);

  /* other rows, the same columns */
  tmc_dense *ok = tmc_dense_create(7, 3, L, 0.0);
  tmc_array_arg args3[] = { tmc_array_dense(ok) };
  EXPECT_APPLY_OK(f, args3);

  tmc_array_function_free(f);
  tmc_dense_free(m_create);
  tmc_dense_free(wrong_cols);
  tmc_dense_free(wrong_cols2);
  tmc_dense_free(ok);
  return TMC_TEST_PASS;
}

/* A[1] = A[0] + B[0] */
static void column_sum_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *ops[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 1, 0), tmc_array_group_column(g, 0, 1) };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, ops);
}

int dense_policy_mismatched_row_dimensions(int L)
{
  /* matrices of different row counts at creation are accepted */
  tmc_dense *a = tmc_dense_create(3, 3, L, 1.0);
  tmc_dense *b = tmc_dense_create(4, 3, L, 2.0);
  tmc_array_arg create_args[] = { tmc_array_dense(a), tmc_array_const_dense(b) };
  tmc_array_function *f = make_function(column_sum_kernel, NULL, N_ARGS(create_args), create_args, L);
  if (!f)
    return TMC_TEST_FAIL;

  tmc_dense *c = tmc_dense_create(5, 3, L, 0.0);
  tmc_dense *d = tmc_dense_create(5, 3, L, 1.0);
  tmc_array_arg args[] = { tmc_array_dense(c), tmc_array_const_dense(d) };
  EXPECT_APPLY_OK(f, args);

  tmc_dense *e = tmc_dense_create(3, 3, L, 0.0);
  tmc_dense *f_matrix = tmc_dense_create(4, 3, L, 1.0);
  tmc_array_arg args2[] = { tmc_array_dense(e), tmc_array_const_dense(f_matrix) };
  EXPECT_APPLY_ERROR(f, args2, TMC_MATRIX_INVALID_VECTOR);

  tmc_array_function_free(f);
  tmc_dense_free(a);
  tmc_dense_free(b);
  tmc_dense_free(c);
  tmc_dense_free(d);
  tmc_dense_free(e);
  tmc_dense_free(f_matrix);
  return TMC_TEST_PASS;
}

/* col[*ctx] = 2 col0, for a column the matrix lacks */
static void missing_column_kernel(tmc_array_group *g, void *ctx)
{
  const int *column = ctx;
  double *ops[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, *column) };
  tmc_array_group_for_each(g, twice_first_into_second, NULL, 2, ops);
}

int dense_policy_mismatched_column_dimensions(int L)
{
  tmc_dense *m = tmc_dense_create(3, 4, L, 1.0);
  tmc_array_arg args[] = { tmc_array_dense(m) };
  /* MICM's column 5 and (tighter) column 4, the first that does not exist */
  int missing[] = { 5, 4 };
  for (int k = 0; k < 2; ++k)
  {
    tmc_array_function *f = make_function(missing_column_kernel, &missing[k], N_ARGS(args), args, L);
    if (!f)
      return TMC_TEST_FAIL;
    /* MICM asserts on the missing column (EXPECT_DEATH in debug builds); the
     * C accessor reports it and the application fails without writing */
    EXPECT_APPLY_ERROR(f, args, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
    tmc_array_function_free(f);
  }
  for (size_t k = 0; k < m->size; ++k)
    TMC_EXPECT_EQ(m->data[k], 1.0);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

/* col3 = 2 col0 */
static void column_three_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *ops[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, 3) };
  tmc_array_group_for_each(g, twice_first_into_second, NULL, 2, ops);
}

int dense_policy_wrong_matrix_dimensions(int L)
{
  tmc_dense *m1 = tmc_dense_create(3, 4, L, 1.0);
  tmc_dense *m2 = tmc_dense_create(3, 5, L, 2.0);
  tmc_array_arg args1[] = { tmc_array_dense(m1) };
  tmc_array_function *f = make_function(column_three_kernel, NULL, N_ARGS(args1), args1, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args1);
  EXPECT_APPLY_OK(f, args1);

  tmc_array_arg args2[] = { tmc_array_dense(m2) };
  EXPECT_APPLY_ERROR(f, args2, TMC_MATRIX_INVALID_VECTOR);

  tmc_dense *m3 = tmc_dense_create(7, 4, L, 1.0);
  tmc_array_arg args3[] = { tmc_array_dense(m3) };
  EXPECT_APPLY_OK(f, args3);

  tmc_array_function_free(f);
  tmc_dense_free(m1);
  tmc_dense_free(m2);
  tmc_dense_free(m3);
  return TMC_TEST_PASS;
}

/* ---- row variables and column views ------------------------------------- */

static void product_into_last(double *const *x, void *ctx) /* t = a * b */
{
  (void)ctx;
  *x[2] = *x[0] * *x[1];
}

static void difference_into_last(double *const *x, void *ctx) /* c = t1 - t2 */
{
  (void)ctx;
  *x[2] = *x[0] - *x[1];
}

/* tmp1 = col0 * col1; tmp2 = col0 + col1; col2 = tmp1 + tmp2;
 * col3 = tmp1 - tmp2; col4 = tmp1 * tmp2 */
static void multiple_temporaries_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp1 = tmc_array_group_variable(g);
  double *tmp2 = tmc_array_group_variable(g);
  double *product[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, 1), tmp1 };
  tmc_array_group_for_each(g, product_into_last, NULL, 3, product);
  double *sum[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, 1), tmp2 };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
  double *col2[] = { tmp1, tmp2, tmc_array_group_column(g, 0, 2) };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, col2);
  double *col3[] = { tmp1, tmp2, tmc_array_group_column(g, 0, 3) };
  tmc_array_group_for_each(g, difference_into_last, NULL, 3, col3);
  double *col4[] = { tmp1, tmp2, tmc_array_group_column(g, 0, 4) };
  tmc_array_group_for_each(g, product_into_last, NULL, 3, col4);
}

int dense_policy_multiple_temporaries(int L)
{
  tmc_dense *m = tmc_dense_create(4, 5, L, 0.0);
  for (int i = 0; i < m->n_rows; ++i)
  {
    tmc_dense_set(m, i, 0, i + 1);
    tmc_dense_set(m, i, 1, (i + 1) * 10);
  }
  tmc_array_arg args[] = { tmc_array_dense(m) };
  tmc_array_function *f = make_function(multiple_temporaries_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 2), 10.0 + 11.0);
  TMC_EXPECT_EQ(AT(m, 0, 3), 10.0 - 11.0);
  TMC_EXPECT_EQ(AT(m, 0, 4), 10.0 * 11.0);
  TMC_EXPECT_EQ(AT(m, 1, 2), 40.0 + 22.0);
  TMC_EXPECT_EQ(AT(m, 1, 3), 40.0 - 22.0);
  TMC_EXPECT_EQ(AT(m, 1, 4), 40.0 * 22.0);
  TMC_EXPECT_EQ(AT(m, 3, 2), 160.0 + 44.0);
  TMC_EXPECT_EQ(AT(m, 3, 3), 160.0 - 44.0);
  TMC_EXPECT_EQ(AT(m, 3, 4), 160.0 * 44.0);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

/* the column lanes are looked up once and reused:
 * col1 = 2 col0; col2 = col0 + col1; col3 = col2 * col1 */
static void column_view_reuse_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *col0 = tmc_array_group_column(g, 0, 0);
  double *col1 = tmc_array_group_column(g, 0, 1);
  double *col2 = tmc_array_group_column(g, 0, 2);
  double *col3 = tmc_array_group_column(g, 0, 3);
  double *twice[] = { col0, col1 };
  tmc_array_group_for_each(g, twice_first_into_second, NULL, 2, twice);
  double *sum[] = { col0, col1, col2 };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
  double *product[] = { col2, col1, col3 };
  tmc_array_group_for_each(g, product_into_last, NULL, 3, product);
}

int dense_policy_column_view_reuse(int L)
{
  tmc_dense *m = tmc_dense_create(3, 4, L, 0.0);
  for (int i = 0; i < m->n_rows; ++i)
    tmc_dense_set(m, i, 0, i + 1);
  tmc_array_arg args[] = { tmc_array_dense(m) };
  tmc_array_function *f = make_function(column_view_reuse_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 1), 2.0);
  TMC_EXPECT_EQ(AT(m, 0, 2), 3.0);
  TMC_EXPECT_EQ(AT(m, 0, 3), 6.0);
  TMC_EXPECT_EQ(AT(m, 1, 1), 4.0);
  TMC_EXPECT_EQ(AT(m, 1, 2), 6.0);
  TMC_EXPECT_EQ(AT(m, 1, 3), 24.0);
  TMC_EXPECT_EQ(AT(m, 2, 1), 6.0);
  TMC_EXPECT_EQ(AT(m, 2, 2), 9.0);
  TMC_EXPECT_EQ(AT(m, 2, 3), 54.0);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

static void twice_second_into_first(double *const *x, void *ctx) /* c = 2 t */
{
  (void)ctx;
  *x[0] = 2.0 * *x[1];
}

/* tmp = col0 + col1 + col2; col2 = 2 tmp */
static void function_reusability_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, 1), tmc_array_group_column(g, 0, 2), tmp };
  tmc_array_group_for_each(g, sum3_into_last, NULL, 4, sum);
  double *scale[] = { tmc_array_group_column(g, 0, 2), tmp };
  tmc_array_group_for_each(g, twice_second_into_first, NULL, 2, scale);
}

int dense_policy_function_reusability(int L)
{
  tmc_dense *m1 = tmc_dense_create(2, 3, L, 1.0);
  tmc_array_arg args1[] = { tmc_array_dense(m1) };
  tmc_array_function *f = make_function(function_reusability_kernel, NULL, N_ARGS(args1), args1, L);
  if (!f)
    return TMC_TEST_FAIL;

  for (int i = 0; i < m1->n_rows; ++i)
    for (int j = 0; j < m1->n_cols; ++j)
      tmc_dense_set(m1, i, j, i + j);
  EXPECT_APPLY_OK(f, args1);
  TMC_EXPECT_EQ(AT(m1, 0, 2), 2.0 * (0 + 1 + 2));
  TMC_EXPECT_EQ(AT(m1, 1, 2), 2.0 * (1 + 2 + 3));

  tmc_dense *m2 = tmc_dense_create(2, 3, L, 5.0);
  tmc_array_arg args2[] = { tmc_array_dense(m2) };
  EXPECT_APPLY_OK(f, args2);
  TMC_EXPECT_EQ(AT(m2, 0, 2), 2.0 * (5 + 5 + 5));
  TMC_EXPECT_EQ(AT(m2, 1, 2), 2.0 * (5 + 5 + 5));

  tmc_dense *m3 = tmc_dense_create(2, 3, L, 0.0);
  for (int i = 0; i < m3->n_rows; ++i)
    tmc_dense_set(m3, i, 0, i * 10);
  tmc_array_arg args3[] = { tmc_array_dense(m3) };
  EXPECT_APPLY_OK(f, args3);
  TMC_EXPECT_EQ(AT(m3, 0, 2), 2.0 * (0 + 0 + 0));
  TMC_EXPECT_EQ(AT(m3, 1, 2), 2.0 * (10 + 0 + 0));

  tmc_array_function_free(f);
  tmc_dense_free(m1);
  tmc_dense_free(m2);
  tmc_dense_free(m3);
  return TMC_TEST_PASS;
}

static void accumulate_first(double *const *x, void *ctx) /* sum += val */
{
  *(double *)ctx += *x[0];
}

/* reads only: tmp = col0 + col1; a local sum over col2 */
static void const_matrix_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, 1), tmp };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
  double local_sum = 0.0;
  double *col2[] = { tmc_array_group_column(g, 0, 2) };
  tmc_array_group_for_each(g, accumulate_first, &local_sum, 1, col2);
}

int dense_policy_const_matrix_function(int L)
{
  tmc_dense *m = tmc_dense_create(3, 4, L, 0.0);
  for (int i = 0; i < m->n_rows; ++i)
    for (int j = 0; j < m->n_cols; ++j)
      tmc_dense_set(m, i, j, i * 10 + j);
  const tmc_dense *const_matrix = m;
  tmc_array_arg args[] = { tmc_array_const_dense(const_matrix) };
  tmc_array_function *f = make_function(const_matrix_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 0), 0.0);
  TMC_EXPECT_EQ(AT(m, 1, 2), 12.0);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

static void set_ninety_nine(double *const *x, void *ctx)
{
  (void)ctx;
  *x[0] = 99.0;
}

/* would write 99 to col0; counts its calls */
static void empty_rows_kernel(tmc_array_group *g, void *ctx)
{
  ++*(int *)ctx;
  double *col0[] = { tmc_array_group_column(g, 0, 0) };
  tmc_array_group_for_each(g, set_ninety_nine, NULL, 1, col0);
}

static void no_op_kernel(tmc_array_group *g, void *ctx)
{
  (void)g;
  (void)ctx;
}

int dense_policy_empty_matrix_function(int L)
{
  tmc_dense *empty_rows = tmc_dense_create(0, 3, L, 1.0);
  int calls = 0;
  tmc_array_arg args[] = { tmc_array_dense(empty_rows) };
  tmc_array_function *f = make_function(empty_rows_kernel, &calls, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(calls, 0); /* no rows, no kernel calls */

  tmc_dense *empty_cols = tmc_dense_create(3, 0, L, 1.0);
  tmc_array_arg args2[] = { tmc_array_dense(empty_cols) };
  tmc_array_function *f2 = make_function(no_op_kernel, NULL, N_ARGS(args2), args2, L);
  if (!f2)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f2, args2);

  tmc_array_function_free(f);
  tmc_array_function_free(f2);
  tmc_dense_free(empty_rows);
  tmc_dense_free(empty_cols);
  return TMC_TEST_PASS;
}

/* ---- vectors -------------------------------------------------------------- */

/* tmp = col0 + v, read only */
static void column_plus_vector_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_vector(g, 1), tmp };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
}

/* A function of a 5 x 3 matrix and a vector of `size` values: creation
 * succeeds, application fails. */
static int mismatched_vector_rows(int L, int size)
{
  tmc_dense *m = tmc_dense_create(5, 3, L, 1.0);
  tmc_padded_vector *v = tmc_padded_vector_create(size, L, 0.0);
  tmc_array_arg args[] = { tmc_array_const_dense(m), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(column_plus_vector_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_ERROR(f, args, TMC_MATRIX_INVALID_VECTOR);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

int dense_policy_vector_too_small(int L)
{
  return mismatched_vector_rows(L, 3);
}

int dense_policy_vector_too_large(int L)
{
  return mismatched_vector_rows(L, 10);
}

/* col0 = v */
static void vector_into_column_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *ops[] = { tmc_array_group_vector(g, 1), tmc_array_group_column(g, 0, 0) };
  tmc_array_group_for_each(g, copy_first_to_second, NULL, 2, ops);
}

/* A function of a matrix of `rows` rows and a vector of `size` values;
 * returns the application's result. */
static int apply_vector_into_column(int L, int rows, int size, tmc_error *err)
{
  tmc_dense *m = tmc_dense_create(rows, 3, L, 1.0);
  tmc_padded_vector *v = tmc_padded_vector_create(size, L, 0.0);
  tmc_array_arg args[] = { tmc_array_dense(m), tmc_array_const_vector(v) };
  int rc = TMC_ERROR_INTERNAL;
  tmc_array_function *f = make_function(vector_into_column_kernel, NULL, N_ARGS(args), args, L);
  if (f)
    rc = tmc_array_function_apply(f, N_ARGS(args), args, err);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_padded_vector_free(v);
  return rc;
}

int dense_policy_empty_vector_non_empty_matrix(int L)
{
  tmc_error err;
  tmc_error_clear(&err);
  EXPECT_MATRIX_ERROR(apply_vector_into_column(L, 5, 0, &err), err, TMC_MATRIX_INVALID_VECTOR);
  return TMC_TEST_PASS;
}

int dense_policy_non_empty_vector_empty_matrix(int L)
{
  tmc_error err;
  tmc_error_clear(&err);
  EXPECT_MATRIX_ERROR(apply_vector_into_column(L, 0, 5, &err), err, TMC_MATRIX_INVALID_VECTOR);
  return TMC_TEST_PASS;
}

int dense_policy_empty_vector_empty_matrix(int L)
{
  tmc_error err;
  tmc_error_clear(&err);
  TMC_EXPECT_EQ(apply_vector_into_column(L, 0, 0, &err), TMC_OK);
  return TMC_TEST_PASS;
}

/* col0 = v1 + v2, the vectors arguments 1 and 2 */
static void two_vectors_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *ops[] = { tmc_array_group_vector(g, 1), tmc_array_group_vector(g, 2), tmc_array_group_column(g, 0, 0) };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, ops);
}

int dense_policy_multiple_vectors_different_sizes(int L)
{
  tmc_dense *m = tmc_dense_create(5, 3, L, 1.0);
  tmc_padded_vector *v1 = tmc_padded_vector_create(5, L, 0.0);
  tmc_padded_vector *v2 = tmc_padded_vector_create(3, L, 0.0);
  tmc_array_arg args[] = { tmc_array_dense(m), tmc_array_const_vector(v1), tmc_array_const_vector(v2) };
  tmc_array_function *f = make_function(two_vectors_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_ERROR(f, args, TMC_MATRIX_INVALID_VECTOR);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_padded_vector_free(v1);
  tmc_padded_vector_free(v2);
  return TMC_TEST_PASS;
}

int dense_policy_multiple_vectors_same_size(int L)
{
  tmc_dense *m = tmc_dense_create(5, 3, L, 0.0);
  tmc_padded_vector *v1 = tmc_padded_vector_create(5, L, 0.0);
  tmc_padded_vector *v2 = tmc_padded_vector_create(5, L, 0.0);
  for (int i = 0; i < 5; ++i)
  {
    v1->data[i] = i + 1;
    v2->data[i] = (i + 1) * 10;
  }
  tmc_array_arg args[] = { tmc_array_dense(m), tmc_array_const_vector(v1), tmc_array_const_vector(v2) };
  tmc_array_function *f = make_function(two_vectors_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 0), 1.0 + 10.0);
  TMC_EXPECT_EQ(AT(m, 1, 0), 2.0 + 20.0);
  TMC_EXPECT_EQ(AT(m, 2, 0), 3.0 + 30.0);
  TMC_EXPECT_EQ(AT(m, 3, 0), 4.0 + 40.0);
  TMC_EXPECT_EQ(AT(m, 4, 0), 5.0 + 50.0);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_padded_vector_free(v1);
  tmc_padded_vector_free(v2);
  return TMC_TEST_PASS;
}

/* A[1] = A[0] + v; B[1] = B[0] + v, with A, B, v arguments 0, 1, 2 */
static void two_matrices_one_vector_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *a_ops[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_vector(g, 2), tmc_array_group_column(g, 0, 1) };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, a_ops);
  double *b_ops[] = { tmc_array_group_column(g, 1, 0), tmc_array_group_vector(g, 2), tmc_array_group_column(g, 1, 1) };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, b_ops);
}

int dense_policy_multiple_matrices_one_vector(int L)
{
  tmc_dense *a = tmc_dense_create(4, 2, L, 0.0);
  tmc_dense *b = tmc_dense_create(4, 3, L, 0.0);
  tmc_padded_vector *v = tmc_padded_vector_create(4, L, 0.0);
  for (int i = 0; i < 4; ++i)
  {
    v->data[i] = i * 2;
    tmc_dense_set(a, i, 0, i + 1);
    tmc_dense_set(b, i, 0, i * 10);
  }
  tmc_array_arg args[] = { tmc_array_dense(a), tmc_array_dense(b), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(two_matrices_one_vector_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(a, 0, 1), 1.0 + 0.0);
  TMC_EXPECT_EQ(AT(a, 1, 1), 2.0 + 2.0);
  TMC_EXPECT_EQ(AT(a, 2, 1), 3.0 + 4.0);
  TMC_EXPECT_EQ(AT(a, 3, 1), 4.0 + 6.0);
  TMC_EXPECT_EQ(AT(b, 0, 1), 0.0 + 0.0);
  TMC_EXPECT_EQ(AT(b, 1, 1), 10.0 + 2.0);
  TMC_EXPECT_EQ(AT(b, 2, 1), 20.0 + 4.0);
  TMC_EXPECT_EQ(AT(b, 3, 1), 30.0 + 6.0);
  tmc_array_function_free(f);
  tmc_dense_free(a);
  tmc_dense_free(b);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

/* A[0] = v, with A, B, v arguments 0, 1, 2 */
static void vector_into_first_matrix_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *ops[] = { tmc_array_group_vector(g, 2), tmc_array_group_column(g, 0, 0) };
  tmc_array_group_for_each(g, copy_first_to_second, NULL, 2, ops);
}

/* Matrices of a_rows and b_rows rows and a vector of `size` values: creation
 * succeeds, application fails. */
static int two_matrices_vector_mismatch(int L, int a_rows, int b_rows, int size)
{
  tmc_dense *a = tmc_dense_create(a_rows, 2, L, 0.0);
  tmc_dense *b = tmc_dense_create(b_rows, 3, L, 0.0);
  tmc_padded_vector *v = tmc_padded_vector_create(size, L, 0.0);
  tmc_array_arg args[] = { tmc_array_dense(a), tmc_array_const_dense(b), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(vector_into_first_matrix_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_ERROR(f, args, TMC_MATRIX_INVALID_VECTOR);
  tmc_array_function_free(f);
  tmc_dense_free(a);
  tmc_dense_free(b);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

int dense_policy_multiple_matrices_different_rows_vector(int L)
{
  return two_matrices_vector_mismatch(L, 4, 5, 4);
}

int dense_policy_vector_length_matches_one_matrix_only(int L)
{
  return two_matrices_vector_mismatch(L, 5, 5, 4);
}

int dense_policy_const_vector(int L)
{
  tmc_dense *m = tmc_dense_create(3, 2, L, 0.0);
  const double values[] = { 10.0, 20.0, 30.0 };
  tmc_padded_vector *storage = tmc_padded_vector_from_values(values, 3, L);
  const tmc_padded_vector *const_vec = storage;
  tmc_array_arg args[] = { tmc_array_dense(m), tmc_array_const_vector(const_vec) };
  tmc_array_function *f = make_function(vector_doubled_into_column_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 0), 20.0);
  TMC_EXPECT_EQ(AT(m, 1, 0), 40.0);
  TMC_EXPECT_EQ(AT(m, 2, 0), 60.0);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_padded_vector_free(storage);
  return TMC_TEST_PASS;
}

/* v = 3 col0: the matrix is read, the vector written */
static void column_tripled_into_vector_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *ops[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_vector(g, 1) };
  tmc_array_group_for_each(g, three_times_first_into_second, NULL, 2, ops);
}

int dense_policy_mutable_vector(int L)
{
  tmc_dense *m = tmc_dense_create(3, 2, L, 0.0);
  const double values[] = { 5.0, 10.0, 15.0 };
  tmc_padded_vector *v = tmc_padded_vector_from_values(values, 3, L);
  for (int i = 0; i < 3; ++i)
    tmc_dense_set(m, i, 0, i + 1);
  tmc_array_arg args[] = { tmc_array_const_dense(m), tmc_array_vector(v) };
  tmc_array_function *f = make_function(column_tripled_into_vector_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(v->data[0], 3.0);
  TMC_EXPECT_EQ(v->data[1], 6.0);
  TMC_EXPECT_EQ(v->data[2], 9.0);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

static void ten_times_first_into_second(double *const *x, void *ctx) /* b = 10 a */
{
  (void)ctx;
  *x[1] = *x[0] * 10.0;
}

/* col0 = 10 v */
static void vector_times_ten_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *ops[] = { tmc_array_group_vector(g, 1), tmc_array_group_column(g, 0, 0) };
  tmc_array_group_for_each(g, ten_times_first_into_second, NULL, 2, ops);
}

int dense_policy_function_reusability_with_vectors(int L)
{
  tmc_dense *m = tmc_dense_create(3, 2, L, 0.0);
  const double values1[] = { 1.0, 2.0, 3.0 }, values2[] = { 5.0, 6.0, 7.0 };
  tmc_padded_vector *v1 = tmc_padded_vector_from_values(values1, 3, L);
  tmc_array_arg args1[] = { tmc_array_dense(m), tmc_array_const_vector(v1) };
  tmc_array_function *f = make_function(vector_times_ten_kernel, NULL, N_ARGS(args1), args1, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args1);
  TMC_EXPECT_EQ(AT(m, 0, 0), 10.0);
  TMC_EXPECT_EQ(AT(m, 1, 0), 20.0);
  TMC_EXPECT_EQ(AT(m, 2, 0), 30.0);

  tmc_padded_vector *v2 = tmc_padded_vector_from_values(values2, 3, L);
  tmc_array_arg args2[] = { tmc_array_dense(m), tmc_array_const_vector(v2) };
  EXPECT_APPLY_OK(f, args2);
  TMC_EXPECT_EQ(AT(m, 0, 0), 50.0);
  TMC_EXPECT_EQ(AT(m, 1, 0), 60.0);
  TMC_EXPECT_EQ(AT(m, 2, 0), 70.0);

  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_padded_vector_free(v1);
  tmc_padded_vector_free(v2);
  return TMC_TEST_PASS;
}

int dense_policy_function_invocation_with_wrong_sized_vector(int L)
{
  tmc_dense *m = tmc_dense_create(3, 2, L, 0.0);
  tmc_padded_vector *v_correct = tmc_padded_vector_create(3, L, 0.0);
  tmc_padded_vector *v_wrong = tmc_padded_vector_create(5, L, 0.0);
  tmc_array_arg args[] = { tmc_array_dense(m), tmc_array_const_vector(v_correct) };
  tmc_array_function *f = make_function(vector_into_column_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  /* MICM leaves the wrong-sized call commented out (it warns at compile
   * time); the C function rejects it */
  tmc_array_arg wrong[] = { tmc_array_dense(m), tmc_array_const_vector(v_wrong) };
  EXPECT_APPLY_ERROR(f, wrong, TMC_MATRIX_INVALID_VECTOR);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_padded_vector_free(v_correct);
  tmc_padded_vector_free(v_wrong);
  return TMC_TEST_PASS;
}

/* tmp = col0 + col1; col2 = tmp + v */
static void mixed_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, 1), tmp };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
  double *result[] = { tmp, tmc_array_group_vector(g, 1), tmc_array_group_column(g, 0, 2) };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, result);
}

int dense_policy_mixed_vector_column_view_row_variable(int L)
{
  tmc_dense *m = tmc_dense_create(4, 3, L, 0.0);
  tmc_padded_vector *v = tmc_padded_vector_create(4, L, 0.0);
  for (int i = 0; i < 4; ++i)
  {
    tmc_dense_set(m, i, 0, i + 1);
    tmc_dense_set(m, i, 1, (i + 1) * 10);
    v->data[i] = (i + 1) * 100;
  }
  tmc_array_arg args[] = { tmc_array_dense(m), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(mixed_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 2), 111.0);
  TMC_EXPECT_EQ(AT(m, 1, 2), 222.0);
  TMC_EXPECT_EQ(AT(m, 2, 2), 333.0);
  TMC_EXPECT_EQ(AT(m, 3, 2), 444.0);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

int dense_policy_integer_vector(int L)
{
  /* MICM's Matrix<int> and VectorType<int>: the same integers held in doubles */
  tmc_dense *m = tmc_dense_create(3, 2, L, 0);
  const double values[] = { 10, 20, 30 };
  tmc_padded_vector *v = tmc_padded_vector_from_values(values, 3, L);
  tmc_array_arg args[] = { tmc_array_dense(m), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(vector_doubled_into_column_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 0), 20);
  TMC_EXPECT_EQ(AT(m, 1, 0), 40);
  TMC_EXPECT_EQ(AT(m, 2, 0), 60);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

int dense_policy_function_with_const_signature(int L)
{
  tmc_dense *m = tmc_dense_create(3, 2, L, 0.0);
  tmc_padded_vector *v = tmc_padded_vector_create(3, L, 0.0);
  tmc_array_arg create_args[] = { tmc_array_dense(m), tmc_array_vector(v) };
  tmc_array_function *f = make_function(vector_doubled_into_column_kernel, NULL, N_ARGS(create_args), create_args, L);
  if (!f)
    return TMC_TEST_FAIL;
  /* MICM wraps the function in a std::function taking (Matrix&, const
   * Vector&); here, a pointer to the apply function, called with the vector
   * read-only */
  int (*func_std)(const tmc_array_function *, int, const tmc_array_arg *, tmc_error *) = tmc_array_function_apply;
  tmc_array_arg args[] = { tmc_array_dense(m), tmc_array_const_vector(v) };
  tmc_error err;
  tmc_error_clear(&err);
  TMC_EXPECT_EQ(func_std(f, N_ARGS(args), args, &err), TMC_OK);
  for (int i = 0; i < 3; ++i)
    TMC_EXPECT_EQ(AT(m, i, 0), 0.0);
  tmc_array_function_free(f);
  tmc_dense_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

/* ---- Fill and Copy -------------------------------------------------------- */

static void fill_column_one_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  tmc_array_group_fill(g, tmc_array_group_column(g, 0, 1), 3.2);
}

static void fill_vector_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  tmc_array_group_fill(g, tmc_array_group_vector(g, 1), 3.2);
}

/* tmp = 9.9, then col0 = tmp */
static void fill_variable_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp = tmc_array_group_variable(g);
  tmc_array_group_fill(g, tmp, 9.9);
  double *ops[] = { tmc_array_group_column(g, 0, 0), tmp };
  tmc_array_group_for_each(g, copy_second_to_first, NULL, 2, ops);
}

int dense_policy_fill(int L)
{
  tmc_dense *m = tmc_dense_create(3, 2, L, 0.0);

  /* a matrix column */
  tmc_array_arg args[] = { tmc_array_dense(m) };
  tmc_array_function *f = make_function(fill_column_one_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  for (int i = 0; i < 3; ++i)
  {
    TMC_EXPECT_EQ(AT(m, i, 0), 0.0);
    TMC_EXPECT_EQ(AT(m, i, 1), 3.2);
  }
  tmc_array_function_free(f);

  /* a caller's vector, through a read-only matrix */
  tmc_padded_vector *v = tmc_padded_vector_create(3, L, 0.0);
  tmc_array_arg vector_args[] = { tmc_array_const_dense(m), tmc_array_vector(v) };
  f = make_function(fill_vector_kernel, NULL, N_ARGS(vector_args), vector_args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, vector_args);
  TMC_EXPECT_EQ(v->data[0], 3.2);
  TMC_EXPECT_EQ(v->data[1], 3.2);
  TMC_EXPECT_EQ(v->data[2], 3.2);
  tmc_array_function_free(f);

  /* a row variable, broadcast into column 0 */
  f = make_function(fill_variable_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 0), 9.9);
  TMC_EXPECT_EQ(AT(m, 1, 0), 9.9);
  TMC_EXPECT_EQ(AT(m, 2, 0), 9.9);
  tmc_array_function_free(f);

  tmc_dense_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

/* vector into column 1 (only the real rows, as MICM's dense GroupView copies a vector) */
static void copy_vector_into_column_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  tmc_array_group_copy_strict(g, tmc_array_group_column(g, 0, 1), tmc_array_group_vector(g, 1));
}

static void copy_column_into_vector_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  tmc_array_group_copy(g, tmc_array_group_vector(g, 1), tmc_array_group_column(g, 0, 1));
}

static void copy_column_one_into_zero_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  tmc_array_group_copy(g, tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, 1));
}

/* column 1 -> row variable; column 1 = 0; column 1 = row variable */
static void copy_round_trip_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *tmp = tmc_array_group_variable(g);
  tmc_array_group_copy(g, tmp, tmc_array_group_column(g, 0, 1));
  tmc_array_group_fill(g, tmc_array_group_column(g, 0, 1), 0.0);
  double *ops[] = { tmc_array_group_column(g, 0, 1), tmp };
  tmc_array_group_for_each(g, copy_second_to_first, NULL, 2, ops);
}

int dense_policy_copy(int L)
{
  tmc_dense *m = tmc_dense_create(3, 2, L, 0.0);
  const double values[] = { 3.2, 4.2, 1.3 };
  tmc_padded_vector *v = tmc_padded_vector_from_values(values, 3, L);
  /* (tighter than MICM) a sentinel in the vector's padding lanes, to see which
   * copies carry padding along */
  for (int i = v->size; i < v->padded_size; ++i)
    v->data[i] = -1.0;
  const int padded_rows = m->n_groups * L;

  /* a vector into a matrix column: the real rows only, so the column's
   * padding keeps its initial 0 */
  tmc_array_arg vector_args[] = { tmc_array_dense(m), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(copy_vector_into_column_kernel, NULL, N_ARGS(vector_args), vector_args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, vector_args);
  for (int i = 0; i < 3; ++i)
  {
    TMC_EXPECT_EQ(AT(m, i, 0), 0.0);
    TMC_EXPECT_EQ(AT(m, i, 1), values[i]);
  }
  for (int i = 3; i < padded_rows; ++i)
    TMC_EXPECT_EQ(m->data[tmc_dense_index(m, i, 1)], 0.0);
  tmc_array_function_free(f);

  /* a read-only matrix column into a vector: all lanes, so the vector's
   * padding takes the column's padding (0) */
  tmc_padded_vector *v2 = tmc_padded_vector_create(3, L, -1.0);
  tmc_array_arg vector2_args[] = { tmc_array_const_dense(m), tmc_array_vector(v2) };
  f = make_function(copy_column_into_vector_kernel, NULL, N_ARGS(vector2_args), vector2_args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, vector2_args);
  TMC_EXPECT_EQ(v2->data[0], 3.2);
  TMC_EXPECT_EQ(v2->data[1], 4.2);
  TMC_EXPECT_EQ(v2->data[2], 1.3);
  for (int i = v2->size; i < v2->padded_size; ++i)
    TMC_EXPECT_EQ(v2->data[i], 0.0);
  tmc_array_function_free(f);

  /* one matrix column into another */
  tmc_array_arg args[] = { tmc_array_dense(m) };
  f = make_function(copy_column_one_into_zero_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  for (int i = 0; i < 3; ++i)
  {
    TMC_EXPECT_EQ(AT(m, i, 0), values[i]);
    TMC_EXPECT_EQ(AT(m, i, 1), values[i]);
  }

  /* again, into a cleared column 0: MICM reads column 1 through a const view
   * this time, which in C is the same accessor */
  for (int i = 0; i < m->n_rows; ++i)
    tmc_dense_set(m, i, 0, 0.0);
  EXPECT_APPLY_OK(f, args);
  for (int i = 0; i < 3; ++i)
    TMC_EXPECT_EQ(AT(m, i, 0), values[i]);
  tmc_array_function_free(f);

  /* column -> row variable -> column */
  f = make_function(copy_round_trip_kernel, NULL, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  for (int i = 0; i < 3; ++i)
    TMC_EXPECT_EQ(AT(m, i, 1), values[i]);
  tmc_array_function_free(f);

  tmc_dense_free(m);
  tmc_padded_vector_free(v);
  tmc_padded_vector_free(v2);
  return TMC_TEST_PASS;
}

/* ---- Reduce ---------------------------------------------------------------- */

/* Runs `kernel` on a read-only matrix with ctx as its context. */
static int apply_read_only(tmc_dense *m, tmc_array_kernel kernel, void *ctx, int L)
{
  tmc_array_arg args[] = { tmc_array_const_dense(m) };
  tmc_array_function *f = make_function(kernel, ctx, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  tmc_array_function_free(f);
  return TMC_TEST_PASS;
}

/* Reduces columns 0 and 1 with `reduce` into the accumulator in ctx. */
typedef struct two_column_reduction
{
  tmc_array_reduce_fn reduce;
  void *accumulator;
} two_column_reduction;

static void reduce_two_columns_kernel(tmc_array_group *g, void *ctx)
{
  const two_column_reduction *r = ctx;
  double *col0[] = { tmc_array_group_column(g, 0, 0) };
  tmc_array_group_reduce(g, r->reduce, NULL, 1, col0, r->accumulator);
  double *col1[] = { tmc_array_group_column(g, 0, 1) };
  tmc_array_group_reduce(g, r->reduce, NULL, 1, col1, r->accumulator);
}

static void sum_of_squares(double *const *x, void *accumulator, void *ctx)
{
  (void)ctx;
  *(double *)accumulator += *x[0] * *x[0];
}

int dense_policy_reduce_sum(int L)
{
  /* the padding lanes (L > 1) hold the initial 0, the sum's identity */
  tmc_dense *m = tmc_dense_create(3, 2, L, 0.0);
  tmc_dense_set(m, 0, 0, 1.0);
  tmc_dense_set(m, 0, 1, 2.0);
  tmc_dense_set(m, 1, 0, 3.0);
  tmc_dense_set(m, 1, 1, 4.0);
  tmc_dense_set(m, 2, 0, 5.0);
  tmc_dense_set(m, 2, 1, 6.0);
  double total = 0.0;
  two_column_reduction r = { sum_of_squares, &total };
  int rc = apply_read_only(m, reduce_two_columns_kernel, &r, L);
  TMC_EXPECT_EQ(total, 91.0);
  tmc_dense_free(m);
  return rc;
}

static void maximum(double *const *x, void *accumulator, void *ctx)
{
  (void)ctx;
  double *acc = accumulator;
  if (*x[0] > *acc)
    *acc = *x[0];
}

int dense_policy_reduce_max(int L)
{
  tmc_dense *m = tmc_dense_create(3, 2, L, 0.0);
  tmc_dense_set(m, 0, 0, 1.0);
  tmc_dense_set(m, 0, 1, 8.0);
  tmc_dense_set(m, 1, 0, 5.0);
  tmc_dense_set(m, 1, 1, 2.0);
  tmc_dense_set(m, 2, 0, 3.0);
  tmc_dense_set(m, 2, 1, 4.0);
  double max_val = -DBL_MAX; /* std::numeric_limits<double>::lowest() */
  two_column_reduction r = { maximum, &max_val };
  int rc = apply_read_only(m, reduce_two_columns_kernel, &r, L);
  TMC_EXPECT_EQ(max_val, 8.0);
  tmc_dense_free(m);
  return rc;
}

static void any_above_ten(double *const *x, void *accumulator, void *ctx)
{
  (void)ctx;
  bool *acc = accumulator;
  *acc = *acc || (*x[0] > 10.0);
}

int dense_policy_reduce_lor(int L)
{
  int rc = TMC_TEST_PASS;
  /* no element above the threshold */
  {
    tmc_dense *m = tmc_dense_create(3, 2, L, 1.0);
    bool any_large = false;
    two_column_reduction r = { any_above_ten, &any_large };
    if (apply_read_only(m, reduce_two_columns_kernel, &r, L) != TMC_TEST_PASS)
      rc = TMC_TEST_FAIL;
    TMC_EXPECT_FALSE(any_large);
    tmc_dense_free(m);
  }
  /* one element above it */
  {
    tmc_dense *m = tmc_dense_create(3, 2, L, 1.0);
    tmc_dense_set(m, 2, 1, 42.0);
    bool any_large = false;
    two_column_reduction r = { any_above_ten, &any_large };
    if (apply_read_only(m, reduce_two_columns_kernel, &r, L) != TMC_TEST_PASS)
      rc = TMC_TEST_FAIL;
    TMC_EXPECT_TRUE(any_large);
    tmc_dense_free(m);
  }
  return rc;
}

static void all_finite(double *const *x, void *accumulator, void *ctx)
{
  (void)ctx;
  bool *acc = accumulator;
  *acc = *acc && isfinite(*x[0]);
}

int dense_policy_reduce_land(int L)
{
  int rc = TMC_TEST_PASS;
  /* every element finite; so are the padding lanes (the initial 1.0) */
  {
    tmc_dense *m = tmc_dense_create(3, 2, L, 1.0);
    bool finite = true;
    two_column_reduction r = { all_finite, &finite };
    if (apply_read_only(m, reduce_two_columns_kernel, &r, L) != TMC_TEST_PASS)
      rc = TMC_TEST_FAIL;
    TMC_EXPECT_TRUE(finite);
    tmc_dense_free(m);
  }
  /* one element not finite */
  {
    tmc_dense *m = tmc_dense_create(3, 2, L, 1.0);
    tmc_dense_set(m, 1, 0, NAN);
    bool finite = true;
    two_column_reduction r = { all_finite, &finite };
    if (apply_read_only(m, reduce_two_columns_kernel, &r, L) != TMC_TEST_PASS)
      rc = TMC_TEST_FAIL;
    TMC_EXPECT_FALSE(finite);
    tmc_dense_free(m);
  }
  return rc;
}

static void count_lanes(double *const *x, void *accumulator, void *ctx)
{
  (void)x;
  (void)ctx;
  *(double *)accumulator += 1.0;
}

static void count_visits(double *const *x, void *ctx)
{
  (void)x;
  *(double *)ctx += 1.0;
}

typedef struct lane_counts
{
  double strict, all;                   /* reduce_strict, reduce */
  double for_each_strict, for_each_all; /* for_each_strict, for_each */
} lane_counts;

static void count_lanes_kernel(tmc_array_group *g, void *ctx)
{
  lane_counts *counts = ctx;
  double *col0[] = { tmc_array_group_column(g, 0, 0) };
  tmc_array_group_reduce_strict(g, count_lanes, NULL, 1, col0, &counts->strict);
  tmc_array_group_reduce(g, count_lanes, NULL, 1, col0, &counts->all);
  tmc_array_group_for_each_strict(g, count_visits, &counts->for_each_strict, 1, col0);
  tmc_array_group_for_each(g, count_visits, &counts->for_each_all, 1, col0);
}

int dense_policy_reduce_strict(int L)
{
  /* 3 rows: for L > 1 the last group is part padding */
  tmc_dense *m = tmc_dense_create(3, 1, L, 0.0);
  lane_counts counts = { 0.0, 0.0, 0.0, 0.0 };
  int rc = apply_read_only(m, count_lanes_kernel, &counts, L);
  TMC_EXPECT_EQ(counts.strict, 3.0);
  /* (tighter than MICM) the non-strict variants visit the padding lanes too,
   * as MICM's Reduce and ForEachRow do; the strict ones only the 3 rows */
  TMC_EXPECT_EQ(counts.all, (double)(m->n_groups * L));
  TMC_EXPECT_EQ(counts.for_each_strict, 3.0);
  TMC_EXPECT_EQ(counts.for_each_all, (double)(m->n_groups * L));
  tmc_dense_free(m);
  return rc;
}

/* ---- captured vectors ------------------------------------------------------ */

typedef struct captured_indices
{
  const size_t *vec1; /* { 2, 0, 1 } */
  int n1;
  const size_t *vec2; /* { 10, 20, 30 } */
} captured_indices;

static void add_captured_value(double *const *x, void *ctx)
{
  *x[0] += *(const double *)ctx;
}

/* for each v1 in vec1: col0 += vec2[v1]; col1 += vec2[vec1[1]] */
static void vector_capture_kernel(tmc_array_group *g, void *ctx)
{
  const captured_indices *c = ctx;
  for (int k = 0; k < c->n1; ++k)
  {
    double value = (double)c->vec2[c->vec1[k]];
    double *col0[] = { tmc_array_group_column(g, 0, 0) };
    tmc_array_group_for_each(g, add_captured_value, &value, 1, col0);
  }
  double value = (double)c->vec2[c->vec1[1]];
  double *col1[] = { tmc_array_group_column(g, 0, 1) };
  tmc_array_group_for_each(g, add_captured_value, &value, 1, col1);
}

int dense_policy_vector_capture(int L)
{
  tmc_dense *m = tmc_dense_create(3, 3, L, 0.0);
  const size_t vec1[] = { 2, 0, 1 }, vec2[] = { 10, 20, 30 };
  captured_indices c = { vec1, 3, vec2 };
  tmc_array_arg args[] = { tmc_array_dense(m) };
  tmc_array_function *f = make_function(vector_capture_kernel, &c, N_ARGS(args), args, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  for (int i = 0; i < 3; ++i)
  {
    TMC_EXPECT_EQ(AT(m, i, 0), 60.0);
    TMC_EXPECT_EQ(AT(m, i, 1), 10.0);
    TMC_EXPECT_EQ(AT(m, i, 2), 0.0);
  }
  tmc_array_function_free(f);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}
