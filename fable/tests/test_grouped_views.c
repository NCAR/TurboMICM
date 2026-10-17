/* MICM test/unit/util/test_grouped_views.cpp: the views a GroupView of one
 * group of a matrix hands out (L rows of a dense matrix, L blocks of a sparse
 * one). A GroupView of group g becomes a kernel that tmc_array_function_apply_group
 * runs on group g; its grouped column and block views become the lane pointers
 * of tmc_array_group_column and tmc_array_group_element, which point at the
 * group's L contiguous lanes, as MICM's base_ and group_base_ + block_offset_
 * do. MICM's vector index of a sparse element, VectorIndex(0, row, col), is
 * the element id from tmc_sparse_element here. */
#include "harness.h"

#include "tmc_dense_ext.h"

#define N_ARGS(a) ((int)(sizeof(a) / sizeof((a)[0])))

/* Runs kernel on group `group` of the arguments, or on every group when group
 * is negative; returns the result. */
static int run(tmc_array_kernel kernel, void *ctx, int n_args, const tmc_array_arg *args, int group)
{
  tmc_error err;
  tmc_error_clear(&err);
  tmc_array_function *f = tmc_array_function_create(kernel, ctx, n_args, args, &err);
  if (!f)
  {
    tmc_test_fail(__FILE__, __LINE__, "function creation failed: %s", err.message);
    return err.category;
  }
  int rc = group < 0 ? tmc_array_function_apply(f, n_args, args, &err) : tmc_array_function_apply_group(f, n_args, args, group, &err);
  if (rc != TMC_OK)
    tmc_test_fail(__FILE__, __LINE__, "applying the function failed: %s", err.message);
  tmc_array_function_free(f);
  return rc;
}

enum
{
  ALL_GROUPS = -1
};

static tmc_sparse *sparse_matrix(int n, int n_elements, const int (*elements)[2], int n_blocks, tmc_sparse_ordering ordering, int L)
{
  tmc_sparse_builder *b = tmc_sparse_builder_create(n);
  for (int k = 0; k < n_elements; ++k)
    tmc_sparse_builder_with_element(b, elements[k][0], elements[k][1], NULL);
  tmc_sparse_builder_set_number_of_blocks(b, n_blocks);
  tmc_sparse *m = tmc_sparse_build(b, ordering, L, false);
  tmc_sparse_builder_free(b);
  return m;
}

static void second_doubled_into_first(double *const *x, void *ctx) /* a = b * 2 */
{
  (void)ctx;
  *x[0] = *x[1] * 2.0;
}

/* ---- category tags --------------------------------------------------------- */

TMC_TEST(category_tags_are_distinct)
{
  TMC_NA("compile-time only: static_asserts on C++ view-category tag types and concepts; C has no type-based dispatch to check");
}

/* ---- dense matrices -------------------------------------------------------- */

typedef struct matrix_round_trip
{
  double *col2, *col2_const, *col3_const; /* the lanes of the views */
  double col2_value, col2_const_value, col3_const_value;
} matrix_round_trip;

static void matrix_round_trip_kernel(tmc_array_group *g, void *ctx)
{
  matrix_round_trip *r = ctx;
  r->col2 = tmc_array_group_column(g, 0, 2);
  r->col2_const = tmc_array_group_column(g, 0, 2);
  r->col3_const = tmc_array_group_column(g, 0, 3);
  if (!r->col2 || !r->col2_const || !r->col3_const)
    return;
  r->col2_value = *r->col2;
  r->col2_const_value = *r->col2_const;
  r->col3_const_value = *r->col3_const;
  double *ops[] = { r->col2, r->col3_const };
  tmc_array_group_for_each(g, second_doubled_into_first, NULL, 2, ops);
}

TMC_TEST(matrix_column_view_round_trip)
{
  tmc_dense *m = tmc_dense_create(3, 4, 1, 0.0);
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 4; ++c)
      tmc_dense_set(m, r, c, 10 * r + c);

  /* row 1 is the group */
  matrix_round_trip r = { .col2 = NULL };
  tmc_array_arg args[] = { tmc_array_dense(m) };
  TMC_ASSERT_TRUE(run(matrix_round_trip_kernel, &r, N_ARGS(args), args, 1) == TMC_OK);
  /* the lanes start at the exact element (row 1, column 2) */
  TMC_EXPECT_TRUE(r.col2 == m->data + tmc_dense_index(m, 1, 2));
  TMC_EXPECT_TRUE(r.col3_const == m->data + tmc_dense_index(m, 1, 3));
  TMC_EXPECT_EQ(r.col2_value, 12.0);
  TMC_EXPECT_EQ(r.col2_const_value, 12.0);
  TMC_EXPECT_EQ(r.col3_const_value, 13.0);
  /* col2 = col3 * 2 in row 1 only */
  TMC_EXPECT_EQ(tmc_dense_get(m, 1, 2), 26.0);
  TMC_EXPECT_EQ(tmc_dense_get(m, 0, 2), 2.0);
  TMC_EXPECT_EQ(tmc_dense_get(m, 2, 2), 22.0);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

typedef struct contiguous_columns
{
  double *col1, *col2_const;
  double col1_values[4], col2_values[4];
} contiguous_columns;

static void contiguous_columns_kernel(tmc_array_group *g, void *ctx)
{
  contiguous_columns *c = ctx;
  c->col1 = tmc_array_group_column(g, 0, 1);
  c->col2_const = tmc_array_group_column(g, 0, 2);
  if (!c->col1 || !c->col2_const)
    return;
  for (int i = 0; i < 4; ++i)
  {
    c->col1_values[i] = c->col1[i];
    c->col2_values[i] = c->col2_const[i];
  }
  double *ops[] = { c->col1, c->col2_const };
  tmc_array_group_for_each(g, second_doubled_into_first, NULL, 2, ops);
}

TMC_TEST(vector_matrix_column_view_contiguous_block)
{
  /* two full groups of L = 4 rows, 3 columns */
  tmc_dense *m = tmc_dense_create(8, 3, 4, 0.0);
  for (int r = 0; r < 8; ++r)
    for (int c = 0; c < 3; ++c)
      tmc_dense_set(m, r, c, 100 + 10 * r + c);

  /* group 1 spans rows 4 to 7 */
  contiguous_columns c = { .col1 = NULL };
  tmc_array_arg args[] = { tmc_array_dense(m) };
  TMC_ASSERT_TRUE(run(contiguous_columns_kernel, &c, N_ARGS(args), args, 1) == TMC_OK);
  /* lane i is row 4 + i of the column */
  TMC_EXPECT_TRUE(c.col1 == m->data + tmc_dense_index(m, 4, 1));
  TMC_EXPECT_TRUE(c.col2_const == m->data + tmc_dense_index(m, 4, 2));
  for (int i = 0; i < 4; ++i)
  {
    TMC_EXPECT_EQ(c.col1_values[i], 100 + 10 * (4 + i) + 1);
    TMC_EXPECT_EQ(c.col2_values[i], 100 + 10 * (4 + i) + 2);
  }
  /* col1 = col2 * 2 over the group; the other group untouched */
  for (int r = 4; r < 8; ++r)
    TMC_EXPECT_EQ(tmc_dense_get(m, r, 1), 2.0 * (100 + 10 * r + 2));
  for (int r = 0; r < 4; ++r)
    TMC_EXPECT_EQ(tmc_dense_get(m, r, 1), 100 + 10 * r + 1);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

static void thrice_plus_one(double *const *x, void *ctx) /* out = a * 3 + 1 */
{
  (void)ctx;
  *x[1] = *x[0] * 3.0 + 1.0;
}

static void grouped_columns_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *ops[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, 3) };
  tmc_array_group_for_each(g, thrice_plus_one, NULL, 2, ops);
}

/* The grouped path (a function over the matrix) gives what the same
 * arithmetic on the matrix's elements gives. */
static int dense_grouped_matches_raw(int L)
{
  const int rows = 3 * L, cols = 4;
  tmc_dense *raw = tmc_dense_create(rows, cols, L, 0.0);
  tmc_dense *grouped = tmc_dense_create(rows, cols, L, 0.0);
  for (int r = 0; r < rows; ++r)
    for (int c = 0; c < cols; ++c)
    {
      tmc_dense_set(raw, r, c, r * 10 + c);
      tmc_dense_set(grouped, r, c, r * 10 + c);
    }
  tmc_array_arg args[] = { tmc_array_dense(grouped) };
  int rc = run(grouped_columns_kernel, NULL, N_ARGS(args), args, ALL_GROUPS) == TMC_OK ? TMC_TEST_PASS : TMC_TEST_FAIL;
  for (int r = 0; r < rows; ++r)
    tmc_dense_set(raw, r, 3, tmc_dense_get(raw, r, 0) * 3.0 + 1.0);
  for (int r = 0; r < rows; ++r)
    for (int c = 0; c < cols; ++c)
      TMC_EXPECT_EQ(tmc_dense_get(raw, r, c), tmc_dense_get(grouped, r, c));
  tmc_dense_free(raw);
  tmc_dense_free(grouped);
  return rc;
}

TMC_TEST(dense_grouped_matches_raw_l1)
{
  return dense_grouped_matches_raw(1);
}

TMC_TEST(dense_grouped_matches_raw_l2)
{
  return dense_grouped_matches_raw(2);
}

TMC_TEST(dense_grouped_matches_raw_l4)
{
  return dense_grouped_matches_raw(4);
}

TMC_TEST(dense_grouped_matches_raw_l8)
{
  return dense_grouped_matches_raw(8);
}

static void difference_into_last(double *const *x, void *ctx) /* out = a - b */
{
  (void)ctx;
  *x[2] = *x[0] - *x[1];
}

static void plain_matrix_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *ops[] = { tmc_array_group_column(g, 0, 0), tmc_array_group_column(g, 0, 1), tmc_array_group_column(g, 0, 2) };
  tmc_array_group_for_each(g, difference_into_last, NULL, 3, ops);
}

TMC_TEST(plain_matrix_grouped_round_trip)
{
  tmc_dense *raw = tmc_dense_create(4, 3, 1, 0.0);
  tmc_dense *grouped = tmc_dense_create(4, 3, 1, 0.0);
  for (int r = 0; r < 4; ++r)
    for (int c = 0; c < 3; ++c)
    {
      tmc_dense_set(raw, r, c, r + c * 5);
      tmc_dense_set(grouped, r, c, r + c * 5);
    }
  tmc_array_arg args[] = { tmc_array_dense(grouped) };
  int rc = run(plain_matrix_kernel, NULL, N_ARGS(args), args, ALL_GROUPS) == TMC_OK ? TMC_TEST_PASS : TMC_TEST_FAIL;
  for (int r = 0; r < 4; ++r)
    tmc_dense_set(raw, r, 2, tmc_dense_get(raw, r, 0) - tmc_dense_get(raw, r, 1));
  for (int r = 0; r < 4; ++r)
    for (int c = 0; c < 3; ++c)
      TMC_EXPECT_EQ(tmc_dense_get(raw, r, c), tmc_dense_get(grouped, r, c));
  tmc_dense_free(raw);
  tmc_dense_free(grouped);
  return rc;
}

/* ---- sparse matrices ------------------------------------------------------- */

typedef struct sparse_round_trip
{
  int e01, e12; /* element ids */
  double *v01, *v12_const;
  double v01_values[4], v12_values[4];
} sparse_round_trip;

static void second_plus_five_into_first(double *const *x, void *ctx) /* a = b + 5 */
{
  (void)ctx;
  *x[0] = *x[1] + 5.0;
}

static void sum_into_first(double *const *x, void *ctx) /* a = a + b */
{
  (void)ctx;
  *x[0] = *x[0] + *x[1];
}

/* Records the lanes of elements e01 and e12 of the sparse argument 0, then
 * applies the lane function in ctx's `update` to them. */
typedef struct sparse_lanes_job
{
  sparse_round_trip *r;
  tmc_array_lane_fn update;
} sparse_lanes_job;

static void sparse_lanes_kernel(tmc_array_group *g, void *ctx)
{
  sparse_lanes_job *job = ctx;
  sparse_round_trip *r = job->r;
  r->v01 = tmc_array_group_element(g, 0, r->e01);
  r->v12_const = tmc_array_group_element(g, 0, r->e12);
  if (!r->v01 || !r->v12_const)
    return;
  for (int i = 0; i < tmc_array_group_width(g); ++i)
  {
    r->v01_values[i] = r->v01[i];
    r->v12_values[i] = r->v12_const[i];
  }
  double *ops[] = { r->v01, r->v12_const };
  tmc_array_group_for_each(g, job->update, NULL, 2, ops);
}

TMC_TEST(sparse_standard_block_view_round_trip)
{
  const int elements[][2] = { { 0, 1 }, { 1, 2 }, { 2, 2 } };
  tmc_sparse *m = sparse_matrix(3, 3, elements, 2, TMC_SPARSE_CSR, 1);
  tmc_sparse_set(m, 0, 0, 1, 1.0);
  tmc_sparse_set(m, 0, 1, 2, 2.0);
  tmc_sparse_set(m, 0, 2, 2, 3.0);
  tmc_sparse_set(m, 1, 0, 1, 10.0);
  tmc_sparse_set(m, 1, 1, 2, 20.0);
  tmc_sparse_set(m, 1, 2, 2, 30.0);

  /* group 1 is block 1 (standard ordering, L = 1) */
  sparse_round_trip r = { .e01 = tmc_sparse_element(m, 0, 1), .e12 = tmc_sparse_element(m, 1, 2) };
  sparse_lanes_job job = { &r, second_plus_five_into_first };
  tmc_array_arg args[] = { tmc_array_sparse(m) };
  TMC_ASSERT_TRUE(run(sparse_lanes_kernel, &job, N_ARGS(args), args, 1) == TMC_OK);
  TMC_EXPECT_TRUE(r.v01 == m->data + tmc_sparse_data_index(m, 1, r.e01));
  TMC_EXPECT_EQ(r.v01_values[0], 10.0);
  TMC_EXPECT_EQ(r.v12_values[0], 20.0);
  TMC_EXPECT_EQ(tmc_sparse_get(m, 1, 0, 1), 25.0);
  TMC_EXPECT_EQ(tmc_sparse_get(m, 0, 0, 1), 1.0); /* untouched */
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

TMC_TEST(sparse_vector_block_view_contiguous_block)
{
  const int elements[][2] = { { 0, 1 }, { 1, 2 } };
  tmc_sparse *m = sparse_matrix(3, 2, elements, 8, TMC_SPARSE_CSR, 4);
  for (int b = 0; b < 8; ++b)
  {
    tmc_sparse_set(m, b, 0, 1, 100.0 + b);
    tmc_sparse_set(m, b, 1, 2, 200.0 + b);
  }

  /* group 1 spans blocks 4 to 7 */
  sparse_round_trip r = { .e01 = tmc_sparse_element(m, 0, 1), .e12 = tmc_sparse_element(m, 1, 2) };
  sparse_lanes_job job = { &r, sum_into_first };
  tmc_array_arg args[] = { tmc_array_sparse(m) };
  TMC_ASSERT_TRUE(run(sparse_lanes_kernel, &job, N_ARGS(args), args, 1) == TMC_OK);
  /* lane i is block 4 + i */
  TMC_EXPECT_TRUE(r.v01 == m->data + tmc_sparse_data_index(m, 4, r.e01));
  for (int i = 0; i < 4; ++i)
  {
    TMC_EXPECT_EQ(r.v01_values[i], 100.0 + (4 + i));
    TMC_EXPECT_EQ(r.v12_values[i], 200.0 + (4 + i));
  }
  /* only the group's 4 blocks change */
  for (int b = 4; b < 8; ++b)
    TMC_EXPECT_EQ(tmc_sparse_get(m, b, 0, 1), (100.0 + b) + (200.0 + b));
  for (int b = 0; b < 4; ++b)
    TMC_EXPECT_EQ(tmc_sparse_get(m, b, 0, 1), 100.0 + b);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

static void sum_last_two_into_first(double *const *x, void *ctx) /* out = a + b */
{
  (void)ctx;
  *x[0] = *x[1] + *x[2];
}

/* dense[:, 0] = dense[:, 1] + sparse[:, 0, 1], sparse argument 0 and dense argument 1 */
static void mixed_kernel(tmc_array_group *g, void *ctx)
{
  const int *e01 = ctx;
  double *ops[] = { tmc_array_group_column(g, 1, 0), tmc_array_group_column(g, 1, 1), tmc_array_group_element(g, 0, *e01) };
  tmc_array_group_for_each(g, sum_last_two_into_first, NULL, 3, ops);
}

static int mixed_sparse_and_dense(int L)
{
  const int blocks = 2 * L;
  const int elements[][2] = { { 0, 1 } };
  tmc_sparse *sparse_ref = sparse_matrix(3, 1, elements, blocks, TMC_SPARSE_CSR, L);
  tmc_sparse *sparse = sparse_matrix(3, 1, elements, blocks, TMC_SPARSE_CSR, L);
  tmc_dense *dense_ref = tmc_dense_create(blocks, 3, L, 0.0);
  tmc_dense *dense = tmc_dense_create(blocks, 3, L, 0.0);
  for (int b = 0; b < blocks; ++b)
  {
    tmc_sparse_set(sparse_ref, b, 0, 1, 1.0 + b);
    tmc_sparse_set(sparse, b, 0, 1, 1.0 + b);
    for (int c = 0; c < 3; ++c)
    {
      tmc_dense_set(dense_ref, b, c, 10.0 + b + c);
      tmc_dense_set(dense, b, c, 10.0 + b + c);
    }
  }

  int e01 = tmc_sparse_element(sparse, 0, 1);
  tmc_array_arg args[] = { tmc_array_sparse(sparse), tmc_array_dense(dense) };
  int rc = run(mixed_kernel, &e01, N_ARGS(args), args, ALL_GROUPS) == TMC_OK ? TMC_TEST_PASS : TMC_TEST_FAIL;

  /* (tighter) SparseMatrix::Function checks at creation that every matrix
   * shares the sparse matrix's L: a standard-ordered dense matrix is refused
   * when L > 1 */
  if (L > 1)
  {
    tmc_dense *standard = tmc_dense_create(blocks, 3, 1, 0.0);
    tmc_array_arg mismatched[] = { tmc_array_sparse(sparse), tmc_array_dense(standard) };
    tmc_error err;
    tmc_error_clear(&err);
    tmc_array_function *f = tmc_array_function_create(mixed_kernel, &e01, N_ARGS(mismatched), mismatched, &err);
    TMC_EXPECT_TRUE(f == NULL);
    TMC_EXPECT_EQ(err.category, TMC_ERROR_MATRIX);
    TMC_EXPECT_EQ(err.code, TMC_MATRIX_INVALID_VECTOR);
    tmc_array_function_free(f);
    tmc_dense_free(standard);
  }

  for (int b = 0; b < blocks; ++b)
    tmc_dense_set(dense_ref, b, 0, tmc_dense_get(dense_ref, b, 1) + tmc_sparse_get(sparse_ref, b, 0, 1));
  for (int b = 0; b < blocks; ++b)
    for (int c = 0; c < 3; ++c)
      TMC_EXPECT_EQ(tmc_dense_get(dense, b, c), tmc_dense_get(dense_ref, b, c));

  tmc_sparse_free(sparse_ref);
  tmc_sparse_free(sparse);
  tmc_dense_free(dense_ref);
  tmc_dense_free(dense);
  return rc;
}

TMC_TEST(mixed_sparse_and_dense_grouped_l1)
{
  return mixed_sparse_and_dense(1);
}

TMC_TEST(mixed_sparse_and_dense_grouped_l2)
{
  return mixed_sparse_and_dense(2);
}

TMC_TEST(mixed_sparse_and_dense_grouped_l4)
{
  return mixed_sparse_and_dense(4);
}

TMC_TEST(mixed_sparse_and_dense_grouped_l8)
{
  return mixed_sparse_and_dense(8);
}

typedef struct row_column_sum
{
  const tmc_sparse *m;
  double accumulator;
} row_column_sum;

static void accumulate_into_ctx(double *const *x, void *ctx)
{
  *(double *)ctx += *x[0];
}

/* the element looked up by its row and column, read only */
static void row_column_sum_kernel(tmc_array_group *g, void *ctx)
{
  row_column_sum *s = ctx;
  double *lanes[] = { tmc_array_group_element(g, 0, tmc_sparse_element(s->m, 0, 1)) };
  tmc_array_group_for_each(g, accumulate_into_ctx, &s->accumulator, 1, lanes);
}

TMC_TEST(sparse_vector_row_col_overload_still_functional)
{
  const int elements[][2] = { { 0, 1 } };
  tmc_sparse *m = sparse_matrix(3, 1, elements, 8, TMC_SPARSE_CSR, 4);
  for (int b = 0; b < 8; ++b)
    tmc_sparse_set(m, b, 0, 1, b + 1);
  row_column_sum s = { m, 0.0 };
  tmc_array_arg args[] = { tmc_array_const_sparse(m) };
  TMC_ASSERT_TRUE(run(row_column_sum_kernel, &s, N_ARGS(args), args, 1) == TMC_OK);
  /* group 1 covers blocks 4 to 7, values 5 to 8 */
  TMC_EXPECT_DOUBLE_EQ(s.accumulator, 26.0);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

/* Copies the 4 lanes of element `element` (sparse argument) or column
 * `column` (dense argument) of argument 0 out. */
typedef struct lanes_out
{
  int element, column;
  double values[4];
} lanes_out;

static void element_lanes_kernel(tmc_array_group *g, void *ctx)
{
  lanes_out *out = ctx;
  const double *lanes = tmc_array_group_element(g, 0, out->element);
  for (int i = 0; lanes && i < 4; ++i)
    out->values[i] = lanes[i];
}

static void column_lanes_kernel(tmc_array_group *g, void *ctx)
{
  lanes_out *out = ctx;
  const double *lanes = tmc_array_group_column(g, 0, out->column);
  for (int i = 0; lanes && i < 4; ++i)
    out->values[i] = lanes[i];
}

TMC_TEST(const_group_view_produces_grouped_views)
{
  const int elements[][2] = { { 0, 1 } };
  tmc_sparse *m = sparse_matrix(3, 1, elements, 4, TMC_SPARSE_CSR, 4);
  for (int b = 0; b < 4; ++b)
    tmc_sparse_set(m, b, 0, 1, b + 1);
  const tmc_sparse *cmatrix = m;
  lanes_out sparse_out = { .element = tmc_sparse_element(m, 0, 1), .column = -1 };
  tmc_array_arg sparse_args[] = { tmc_array_const_sparse(cmatrix) };
  TMC_ASSERT_TRUE(run(element_lanes_kernel, &sparse_out, N_ARGS(sparse_args), sparse_args, 0) == TMC_OK);
  for (int i = 0; i < 4; ++i)
    TMC_EXPECT_EQ(sparse_out.values[i], i + 1);

  tmc_dense *dense = tmc_dense_create(4, 2, 4, 0.0);
  for (int b = 0; b < 4; ++b)
    tmc_dense_set(dense, b, 1, 100 + b);
  const tmc_dense *cdense = dense;
  lanes_out dense_out = { .element = -1, .column = 1 };
  tmc_array_arg dense_args[] = { tmc_array_const_dense(cdense) };
  TMC_ASSERT_TRUE(run(column_lanes_kernel, &dense_out, N_ARGS(dense_args), dense_args, 0) == TMC_OK);
  for (int i = 0; i < 4; ++i)
    TMC_EXPECT_EQ(dense_out.values[i], 100 + i);

  tmc_sparse_free(m);
  tmc_dense_free(dense);
  return TMC_TEST_PASS;
}

/* ---- Fill and Copy on sparse groups ---------------------------------------- */

typedef struct fill_element
{
  int element;
  double value;
} fill_element;

static void fill_element_kernel(tmc_array_group *g, void *ctx)
{
  const fill_element *op = ctx;
  tmc_array_group_fill(g, tmc_array_group_element(g, 0, op->element), op->value);
}

typedef struct copy_element
{
  int dst_arg, dst_element, src_arg, src_element;
} copy_element;

static void copy_element_kernel(tmc_array_group *g, void *ctx)
{
  const copy_element *op = ctx;
  tmc_array_group_copy(g, tmc_array_group_element(g, op->dst_arg, op->dst_element), tmc_array_group_element(g, op->src_arg, op->src_element));
}

/* Fills element `element` of group `group` of m with `value`. */
static int fill_on_group(tmc_sparse *m, int group, int element, double value)
{
  fill_element op = { element, value };
  tmc_array_arg args[] = { tmc_array_sparse(m) };
  return run(fill_element_kernel, &op, N_ARGS(args), args, group);
}

/* Copies element src into element dst within group `group` of m. */
static int copy_on_group(tmc_sparse *m, int group, int dst, int src)
{
  copy_element op = { 0, dst, 0, src };
  tmc_array_arg args[] = { tmc_array_sparse(m) };
  return run(copy_element_kernel, &op, N_ARGS(args), args, group);
}

/* Two non-zero elements per block; a full group and one partial group. */
static int sparse_fill_copy(tmc_sparse_ordering ordering, int L)
{
  const int elements[][2] = { { 0, 1 }, { 2, 0 } };
  tmc_sparse *m = sparse_matrix(3, 2, elements, L + 1, ordering, L);
  for (int b = 0; b < L + 1; ++b)
  {
    tmc_sparse_set(m, b, 0, 1, b + 1);
    tmc_sparse_set(m, b, 2, 0, 100 + b);
  }
  int e01 = tmc_sparse_element(m, 0, 1), e20 = tmc_sparse_element(m, 2, 0);
  int rc = TMC_TEST_PASS;

  /* fill on the full group: only the (0, 1) lanes change */
  if (fill_on_group(m, 0, e01, 42.0) != TMC_OK)
    rc = TMC_TEST_FAIL;
  for (int b = 0; b < L; ++b)
  {
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, b, 0, 1), 42.0);
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, b, 2, 0), 100 + b);
  }
  /* the partial group is unaffected */
  TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, L, 0, 1), L + 1);

  /* copy on the full group: (2, 0) <- (0, 1) */
  if (copy_on_group(m, 0, e20, e01) != TMC_OK)
    rc = TMC_TEST_FAIL;
  for (int b = 0; b < L; ++b)
  {
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, b, 2, 0), 42.0);
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, b, 0, 1), 42.0);
  }

  /* fill on the partial (last) group; its padding lanes are written too */
  const int last_group = m->n_blocks / L;
  if (fill_on_group(m, last_group, e01, -7.0) != TMC_OK)
    rc = TMC_TEST_FAIL;
  TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, L, 0, 1), -7.0);
  /* (tighter than MICM, which says so in a comment) every lane of the group,
   * padding included, holds the fill */
  for (int lane = 0; lane < L; ++lane)
    TMC_EXPECT_DOUBLE_EQ(m->data[tmc_sparse_data_index(m, last_group * L + lane, e01)], -7.0);

  tmc_sparse_free(m);
  return rc;
}

TMC_TEST(fill_and_copy_sparse_standard_csr)
{
  /* standard ordering: L = 1, one group per block */
  const int elements[][2] = { { 0, 1 }, { 2, 0 } };
  tmc_sparse *m = sparse_matrix(3, 2, elements, 2, TMC_SPARSE_CSR, 1);
  for (int b = 0; b < 2; ++b)
  {
    tmc_sparse_set(m, b, 0, 1, b + 1);
    tmc_sparse_set(m, b, 2, 0, 100 + b);
  }
  int e01 = tmc_sparse_element(m, 0, 1), e20 = tmc_sparse_element(m, 2, 0);
  TMC_ASSERT_TRUE(fill_on_group(m, 0, e01, 42.0) == TMC_OK);
  TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, 0, 0, 1), 42.0);
  TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, 0, 2, 0), 100.0); /* untouched */
  TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, 1, 0, 1), 2.0);   /* other group untouched */
  TMC_ASSERT_TRUE(copy_on_group(m, 0, e20, e01) == TMC_OK);
  TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, 0, 2, 0), 42.0);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

TMC_TEST(fill_and_copy_sparse_standard_csc)
{
  const int elements[][2] = { { 0, 1 }, { 2, 0 } };
  tmc_sparse *m = sparse_matrix(3, 2, elements, 2, TMC_SPARSE_CSC, 1);
  for (int b = 0; b < 2; ++b)
  {
    tmc_sparse_set(m, b, 0, 1, b + 1);
    tmc_sparse_set(m, b, 2, 0, 100 + b);
  }
  int e01 = tmc_sparse_element(m, 0, 1), e20 = tmc_sparse_element(m, 2, 0);
  TMC_ASSERT_TRUE(fill_on_group(m, 1, e01, 9.0) == TMC_OK);
  TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, 1, 0, 1), 9.0);
  TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, 0, 0, 1), 1.0); /* other group untouched */
  TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, 1, 2, 0), 101.0);
  TMC_ASSERT_TRUE(copy_on_group(m, 1, e20, e01) == TMC_OK);
  TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(m, 1, 2, 0), 9.0);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

TMC_TEST(fill_and_copy_sparse_vector_csr_l4)
{
  return sparse_fill_copy(TMC_SPARSE_CSR, 4);
}

TMC_TEST(fill_and_copy_sparse_vector_csc_l4)
{
  return sparse_fill_copy(TMC_SPARSE_CSC, 4);
}

TMC_TEST(fill_and_copy_sparse_vector_csr_l8)
{
  return sparse_fill_copy(TMC_SPARSE_CSR, 8);
}

static void set_fifty_five_point_five(double *const *x, void *ctx)
{
  (void)ctx;
  *x[0] = 55.5;
}

static void reference_fill_kernel(tmc_array_group *g, void *ctx)
{
  const int *element = ctx;
  double *lanes[] = { tmc_array_group_element(g, 0, *element) };
  tmc_array_group_for_each(g, set_fifty_five_point_five, NULL, 1, lanes);
}

/* Copy between two matrices' groups: the destination argument 0, the source
 * (read only) argument 1. */
TMC_TEST(fill_copy_cross_group_view_sources)
{
  const int elements[][2] = { { 0, 1 } };
  tmc_sparse *src = sparse_matrix(2, 1, elements, 4, TMC_SPARSE_CSR, 4);
  tmc_sparse *dst = sparse_matrix(2, 1, elements, 4, TMC_SPARSE_CSR, 4);
  for (int b = 0; b < 4; ++b)
  {
    tmc_sparse_set(src, b, 0, 1, b + 1);
    tmc_sparse_set(dst, b, 0, 1, 0.0);
  }
  int e01 = tmc_sparse_element(dst, 0, 1);
  copy_element op = { 0, e01, 1, tmc_sparse_element(src, 0, 1) };
  tmc_array_arg args[] = { tmc_array_sparse(dst), tmc_array_const_sparse(src) };
  TMC_ASSERT_TRUE(run(copy_element_kernel, &op, N_ARGS(args), args, 0) == TMC_OK);
  for (int b = 0; b < 4; ++b)
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(dst, b, 0, 1), b + 1);

  /* the fill matches the equivalent ForEachBlock lambda */
  tmc_sparse *ref = sparse_matrix(2, 1, elements, 4, TMC_SPARSE_CSR, 4);
  for (int b = 0; b < 4; ++b)
    tmc_sparse_set(ref, b, 0, 1, 999.0);
  tmc_array_arg ref_args[] = { tmc_array_sparse(ref) };
  TMC_ASSERT_TRUE(run(reference_fill_kernel, &e01, N_ARGS(ref_args), ref_args, 0) == TMC_OK);
  TMC_ASSERT_TRUE(fill_on_group(dst, 0, e01, 55.5) == TMC_OK);
  for (int b = 0; b < 4; ++b)
    TMC_EXPECT_DOUBLE_EQ(tmc_sparse_get(dst, b, 0, 1), tmc_sparse_get(ref, b, 0, 1));

  tmc_sparse_free(src);
  tmc_sparse_free(dst);
  tmc_sparse_free(ref);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("grouped_views:GroupedView.CategoryTagsAreDistinct", category_tags_are_distinct),
  TMC_CASE("grouped_views:GroupedView.MatrixColumnViewRoundTrip", matrix_column_view_round_trip),
  TMC_CASE("grouped_views:GroupedView.VectorMatrixColumnViewContiguousBlock", vector_matrix_column_view_contiguous_block),
  TMC_CASE("grouped_views:GroupedView.DenseGroupedMatchesRawL1", dense_grouped_matches_raw_l1),
  TMC_CASE("grouped_views:GroupedView.DenseGroupedMatchesRawL2", dense_grouped_matches_raw_l2),
  TMC_CASE("grouped_views:GroupedView.DenseGroupedMatchesRawL4", dense_grouped_matches_raw_l4),
  TMC_CASE("grouped_views:GroupedView.DenseGroupedMatchesRawL8", dense_grouped_matches_raw_l8),
  TMC_CASE("grouped_views:GroupedView.PlainMatrixGroupedRoundTrip", plain_matrix_grouped_round_trip),
  TMC_CASE("grouped_views:GroupedView.SparseStandardBlockViewRoundTrip", sparse_standard_block_view_round_trip),
  TMC_CASE("grouped_views:GroupedView.SparseVectorBlockViewContiguousBlock", sparse_vector_block_view_contiguous_block),
  TMC_CASE("grouped_views:GroupedView.MixedSparseAndDenseGroupedL1", mixed_sparse_and_dense_grouped_l1),
  TMC_CASE("grouped_views:GroupedView.MixedSparseAndDenseGroupedL2", mixed_sparse_and_dense_grouped_l2),
  TMC_CASE("grouped_views:GroupedView.MixedSparseAndDenseGroupedL4", mixed_sparse_and_dense_grouped_l4),
  TMC_CASE("grouped_views:GroupedView.MixedSparseAndDenseGroupedL8", mixed_sparse_and_dense_grouped_l8),
  TMC_CASE("grouped_views:GroupedView.SparseVectorRowColOverloadStillFunctional", sparse_vector_row_col_overload_still_functional),
  TMC_CASE("grouped_views:GroupedView.ConstGroupViewProducesGroupedViews", const_group_view_produces_grouped_views),
  TMC_CASE("grouped_views:GroupedView.FillAndCopySparseStandardCSR", fill_and_copy_sparse_standard_csr),
  TMC_CASE("grouped_views:GroupedView.FillAndCopySparseStandardCSC", fill_and_copy_sparse_standard_csc),
  TMC_CASE("grouped_views:GroupedView.FillAndCopySparseVectorCSR_L4", fill_and_copy_sparse_vector_csr_l4),
  TMC_CASE("grouped_views:GroupedView.FillAndCopySparseVectorCSC_L4", fill_and_copy_sparse_vector_csc_l4),
  TMC_CASE("grouped_views:GroupedView.FillAndCopySparseVectorCSR_L8", fill_and_copy_sparse_vector_csr_l8),
  TMC_CASE("grouped_views:GroupedView.FillCopyCrossGroupViewSources", fill_copy_cross_group_view_sources),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_grouped_views(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
