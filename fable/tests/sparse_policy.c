/* MICM test/unit/util/test_sparse_matrix_policy.hpp: the checks MICM runs on
 * each sparse matrix policy, as functions of the element numbering (MICM's
 * compressed sparse row or column ordering) and the width L (1 for MICM's
 * standard orderings, L for its vector orderings).
 *
 * MICM's Function lambdas become kernels: a GroupView's block views, block
 * variables, column views, and vector views become the lane pointers of
 * tmc_array_group_element, _variable, _column, and _vector, and each
 * ForEachBlock lambda a lane function whose x[k] is the current lane of
 * operand k. MICM's GetBlockView takes VectorIndex(0, row, column), the
 * element's offset within group 0 (e at L = 1, e * L in the vector
 * orderings); tmc_array_group_element takes the element id e itself, as
 * tmc_sparse_element returns it, so the "Scalar" indices of MICM's tests are
 * element ids here. Where MICM throws, the C call
 * returns TMC_ERROR_MATRIX with MICM's code; where MICM asserts (EXPECT_DEATH,
 * compiled only into debug builds) the C accessor returns an error code, which
 * is checked. MICM's int-valued matrices and vectors hold the same integers in
 * doubles here (tmc_sparse stores doubles; every value involved is exact). */
#include "sparse_policy.h"

#include "tmc_dense_ext.h"
#include "tmc_sparse_ext.h"

#define N_ARGS(a) ((int)(sizeof(a) / sizeof((a)[0])))
#define AT(m, b, i, j) tmc_sparse_get((m), (b), (i), (j))
#define SET(m, b, i, j, v) tmc_sparse_set((m), (b), (i), (j), (v))

/* the lanes of a kernel's operands */
#define ELEM(g, arg, e) tmc_array_group_element((g), (arg), (e))
#define COL(g, arg, c) tmc_array_group_column((g), (arg), (c))
#define VEC(g, arg) tmc_array_group_vector((g), (arg))

static const char *ordering_name(tmc_sparse_ordering ordering)
{
  return ordering == TMC_SPARSE_CSR ? "CSR" : "CSC";
}

/* a failure of the check at hand, tagged with its policy */
#define POLICY_FAIL(format, ...) tmc_test_fail(__FILE__, __LINE__, "%s L = %d: " format, ordering_name(ordering), L, __VA_ARGS__)

#define EXPECT_MATRIX_ERROR(rc, err, expected_code)                                           \
  do                                                                                          \
  {                                                                                           \
    int rc_ = (rc);                                                                           \
    TMC_EXPECT_EQ(rc_, TMC_ERROR_MATRIX);                                                     \
    TMC_EXPECT_EQ((err).category, TMC_ERROR_MATRIX);                                          \
    TMC_EXPECT_EQ((err).code, (expected_code));                                               \
  } while (0)

/* an index lookup fails with the matrix error code (MICM: an assertion or a throw) */
#define EXPECT_INDEX_ERROR(call, err, expected_code)                                          \
  do                                                                                          \
  {                                                                                           \
    tmc_error_clear(&(err));                                                                  \
    size_t index_ = (call);                                                                   \
    TMC_EXPECT_EQ(index_, (size_t)-1);                                                        \
    TMC_EXPECT_EQ((err).category, TMC_ERROR_MATRIX);                                          \
    TMC_EXPECT_EQ((err).code, (expected_code));                                               \
  } while (0)

/* IsZero of an element outside the block (MICM: an assertion) */
#define EXPECT_IS_ZERO_OUT_OF_RANGE(m, i, j, err)                                             \
  do                                                                                          \
  {                                                                                           \
    tmc_error_clear(&(err));                                                                  \
    TMC_EXPECT_TRUE(tmc_sparse_is_zero_checked((m), (i), (j), &(err)));                       \
    TMC_EXPECT_EQ((err).category, TMC_ERROR_MATRIX);                                          \
    TMC_EXPECT_EQ((err).code, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);                               \
  } while (0)

/* f applied to args succeeds (MICM: no throw). */
#define EXPECT_APPLY_OK(f, args)                                                              \
  do                                                                                          \
  {                                                                                           \
    tmc_error e_;                                                                             \
    tmc_error_clear(&e_);                                                                     \
    if (tmc_array_function_apply((f), N_ARGS(args), (args), &e_) != TMC_OK)                   \
      POLICY_FAIL("apply failed: %s", e_.message);                                            \
  } while (0)

/* f applied to args fails with the matrix error code (MICM: EXPECT_ANY_THROW). */
#define EXPECT_APPLY_ERROR(f, args, expected_code)                                            \
  do                                                                                          \
  {                                                                                           \
    tmc_error e_;                                                                             \
    tmc_error_clear(&e_);                                                                     \
    EXPECT_MATRIX_ERROR(tmc_array_function_apply((f), N_ARGS(args), (args), &e_), e_, (expected_code)); \
  } while (0)

static tmc_array_function *make_function(tmc_array_kernel kernel, void *ctx, int n_args, const tmc_array_arg *args, tmc_sparse_ordering ordering, int L)
{
  tmc_error err;
  tmc_error_clear(&err);
  tmc_array_function *f = tmc_array_function_create(kernel, ctx, n_args, args, &err);
  if (!f)
    POLICY_FAIL("function creation failed: %s", err.message);
  return f;
}

/* Creating a function of args fails with the matrix error code (MICM: a throw
 * from Function). */
static void expect_create_error(tmc_array_kernel kernel, int n_args, const tmc_array_arg *args, int expected_code, tmc_sparse_ordering ordering, int L)
{
  tmc_error err;
  tmc_error_clear(&err);
  tmc_array_function *f = tmc_array_function_create(kernel, NULL, n_args, args, &err);
  if (f)
    POLICY_FAIL("function creation succeeded with %d incompatible arguments", n_args);
  TMC_EXPECT_EQ(err.category, TMC_ERROR_MATRIX);
  TMC_EXPECT_EQ(err.code, expected_code);
  tmc_array_function_free(f);
}

/* A matrix of n_blocks n x n blocks with the given non-zero elements, every
 * value initial_value (MICM's Create(n).WithElement(...)...). */
static tmc_sparse *build(tmc_sparse_ordering ordering, int L, int n, int n_elements, const int (*elements)[2], int n_blocks, double initial_value)
{
  tmc_sparse_builder *b = tmc_sparse_builder_create(n);
  for (int k = 0; k < n_elements; ++k)
  {
    tmc_error err;
    tmc_error_clear(&err);
    if (tmc_sparse_builder_with_element(b, elements[k][0], elements[k][1], &err) != TMC_OK)
      POLICY_FAIL("element (%d, %d) refused: %s", elements[k][0], elements[k][1], err.message);
  }
  tmc_sparse_builder_set_number_of_blocks(b, n_blocks);
  tmc_sparse_builder_set_initial_value(b, initial_value);
  tmc_sparse *m = tmc_sparse_build(b, ordering, L, false);
  tmc_sparse_builder_free(b);
  return m;
}

/* The element id of (row, col), which must be a non-zero element. */
static int elem(const tmc_sparse *m, int row, int col)
{
  int e = tmc_sparse_element(m, row, col);
  if (e < 0)
    tmc_test_fail(__FILE__, __LINE__, "element (%d, %d) is not in the pattern", row, col);
  return e;
}

/* A padded vector of n values of width L. */
static tmc_padded_vector *vector_of(int L, int n, const double *values)
{
  return tmc_padded_vector_from_values(values, n, L);
}

/* Lane functions: x[k] is the current lane of operand k. */
static void sum2_into_last(double *const *x, void *ctx) /* c = a + b */
{
  (void)ctx;
  *x[2] = *x[0] + *x[1];
}

static void sum3_into_last(double *const *x, void *ctx) /* t = a + b + c */
{
  (void)ctx;
  *x[3] = *x[0] + *x[1] + *x[2];
}

static void sum4_into_last(double *const *x, void *ctx) /* t = a + b + c + d */
{
  (void)ctx;
  *x[4] = *x[0] + *x[1] + *x[2] + *x[3];
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

static void twice_first_into_second(double *const *x, void *ctx) /* b = 2 a */
{
  (void)ctx;
  *x[1] = *x[0] * 2.0;
}

static void twice_second_into_first(double *const *x, void *ctx) /* a = 2 b */
{
  (void)ctx;
  *x[0] = 2.0 * *x[1];
}

static void thrice_first_into_second(double *const *x, void *ctx) /* b = 3 a */
{
  (void)ctx;
  *x[1] = *x[0] * 3.0;
}

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

static void accumulate_first(double *const *x, void *ctx) /* sum += val */
{
  *(double *)ctx += *x[0];
}

/* The element ids a kernel works on; each kernel says what its entries are. */
typedef struct element_ids
{
  int e[6];
} element_ids;

/* ---- construction and element access ------------------------------------ */

tmc_sparse *sparse_policy_zero_matrix(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse_builder *b = tmc_sparse_builder_create(3);
  TMC_EXPECT_EQ(tmc_sparse_builder_number_of_elements(b), 0);
  tmc_sparse *m = tmc_sparse_build(b, ordering, L, false);
  tmc_sparse_builder_free(b);
  TMC_EXPECT_EQ(tmc_sparse_flat_block_size(m), (size_t)0);
  TMC_EXPECT_EQ(m->size, (size_t)0);

  tmc_error err;
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 0, 0, &err), err, TMC_MATRIX_ZERO_ELEMENT_ACCESS);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 6, 0, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 1, 3, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 6, 3, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_IS_ZERO_OUT_OF_RANGE(m, 6, 0, err);
  EXPECT_IS_ZERO_OUT_OF_RANGE(m, 1, 3, err);
  EXPECT_IS_ZERO_OUT_OF_RANGE(m, 6, 3, err);
  /* the element proxies matrix[block][row][column] of MICM's writes */
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 0, 0, 4, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 1, 0, 0, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 0, 3, 0, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 0, 1, 1, &err), err, TMC_MATRIX_ZERO_ELEMENT_ACCESS);
  return m;
}

tmc_sparse *sparse_policy_const_zero_matrix(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse_builder *b = tmc_sparse_builder_create(3);
  TMC_EXPECT_EQ(tmc_sparse_builder_number_of_elements(b), 0);
  tmc_sparse *storage = tmc_sparse_build(b, ordering, L, false);
  tmc_sparse_builder_free(b);
  const tmc_sparse *m = storage;
  TMC_EXPECT_EQ(tmc_sparse_flat_block_size(m), (size_t)0);

  tmc_error err;
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 0, 0, &err), err, TMC_MATRIX_ZERO_ELEMENT_ACCESS);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 6, 0, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 1, 3, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 6, 3, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_IS_ZERO_OUT_OF_RANGE(m, 6, 0, err);
  EXPECT_IS_ZERO_OUT_OF_RANGE(m, 1, 3, err);
  EXPECT_IS_ZERO_OUT_OF_RANGE(m, 6, 3, err);
  return storage;
}

/* 0 X 0 0
 * 0 X 0 0
 * 0 X 0 X
 * 0 0 X 0 */
static const int single_block_elements[5][2] = { { 0, 1 }, { 3, 2 }, { 2, 3 }, { 1, 1 }, { 2, 1 } };

/* The diagonal of block `block` holds one non-zero element, (1, 1). */
static void expect_single_diagonal(const tmc_sparse *m, int block)
{
  tmc_error err;
  tmc_error_clear(&err);
  size_t diagonal[4];
  int count = tmc_sparse_diagonal_indices(m, block, diagonal, &err);
  TMC_EXPECT_EQ(count, 1);
  if (count == 1)
    TMC_EXPECT_EQ(diagonal[0], tmc_sparse_vector_index(m, block, 1, 1, &err));
  TMC_EXPECT_EQ(err.category, TMC_OK);
}

static void expect_single_block_pattern(const tmc_sparse *m)
{
  TMC_EXPECT_EQ(tmc_sparse_flat_block_size(m), (size_t)5);
  TMC_EXPECT_EQ(tmc_sparse_is_zero(m, 0, 1), false);
  TMC_EXPECT_EQ(tmc_sparse_is_zero(m, 3, 2), false);
  TMC_EXPECT_EQ(tmc_sparse_is_zero(m, 2, 1), false);
  TMC_EXPECT_EQ(tmc_sparse_is_zero(m, 2, 2), true);
  TMC_EXPECT_EQ(tmc_sparse_is_zero(m, 0, 0), true);
  TMC_EXPECT_EQ(tmc_sparse_is_zero(m, 3, 3), true);
}

/* MICM's debug-build assertions on a one-block matrix of the pattern above */
static void expect_single_block_index_errors(const tmc_sparse *m)
{
  tmc_error err;
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 4, 2, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 1, 5, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 1, 0, 2, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 2, 2, &err), err, TMC_MATRIX_ZERO_ELEMENT_ACCESS);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 0, 2, &err), err, TMC_MATRIX_ZERO_ELEMENT_ACCESS);
}

tmc_sparse *sparse_policy_single_block_matrix(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse_builder *b = tmc_sparse_builder_create(4);
  for (int k = 0; k < 5; ++k)
    TMC_EXPECT_EQ(tmc_sparse_builder_with_element(b, single_block_elements[k][0], single_block_elements[k][1], NULL), TMC_OK);
  TMC_EXPECT_EQ(tmc_sparse_builder_number_of_elements(b), 5);
  tmc_sparse *m = tmc_sparse_build(b, ordering, L, false);
  tmc_sparse_builder_free(b);

  expect_single_diagonal(m, 0);
  expect_single_block_pattern(m);

  TMC_EXPECT_EQ(AT(m, 0, 2, 1), 0);
  SET(m, 0, 2, 1, 45);
  TMC_EXPECT_EQ(AT(m, 0, 2, 1), 45);

  expect_single_block_index_errors(m);
  /* the element proxies of MICM's writes */
  tmc_error err;
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 0, 0, 4, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 1, 0, 0, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 0, 5, 0, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 0, 3, 3, &err), err, TMC_MATRIX_ZERO_ELEMENT_ACCESS);
  return m;
}

tmc_sparse *sparse_policy_const_single_block_matrix(tmc_sparse_ordering ordering, int L)
{
  /* MICM's SparseMatrix<int>: the same integers in doubles */
  tmc_sparse *orig_matrix = build(ordering, L, 4, 5, single_block_elements, 1, 0);
  SET(orig_matrix, 0, 2, 1, 45);
  SET(orig_matrix, 0, 3, 2, 42);
  SET(orig_matrix, 0, 2, 3, 21);
  tmc_sparse *copy = tmc_sparse_copy(orig_matrix);
  tmc_sparse_free(orig_matrix);
  const tmc_sparse *m = copy;

  expect_single_diagonal(m, 0);
  expect_single_block_pattern(m);
  TMC_EXPECT_EQ(AT(m, 0, 2, 1), 45);
  TMC_EXPECT_EQ(AT(m, 0, 3, 2), 42);
  TMC_EXPECT_EQ(AT(m, 0, 2, 3), 21);
  expect_single_block_index_errors(m);
  return copy;
}

tmc_sparse *sparse_policy_multi_block_matrix(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse_builder *b = tmc_sparse_builder_create(4);
  for (int k = 0; k < 5; ++k)
    TMC_EXPECT_EQ(tmc_sparse_builder_with_element(b, single_block_elements[k][0], single_block_elements[k][1], NULL), TMC_OK);
  tmc_sparse_builder_set_initial_value(b, 24);
  tmc_sparse_builder_set_number_of_blocks(b, 53);
  /* MICM's NumberOfElements counts the elements of every block (5 * 53); the
   * C builder counts the elements of one block */
  TMC_EXPECT_EQ(tmc_sparse_builder_number_of_elements(b), 5);
  TMC_EXPECT_EQ(tmc_sparse_builder_number_of_elements(b) * 53, 5 * 53);
  tmc_sparse *m = tmc_sparse_build(b, ordering, L, false);
  tmc_sparse_builder_free(b);

  TMC_EXPECT_EQ(tmc_sparse_flat_block_size(m), (size_t)5);
  TMC_EXPECT_EQ(m->n_blocks, 53);
  TMC_EXPECT_EQ(m->size, (size_t)((53 + L - 1) / L * L * 5));

  TMC_EXPECT_EQ(AT(m, 0, 2, 1), 24);
  SET(m, 0, 2, 1, 45);
  TMC_EXPECT_EQ(AT(m, 0, 2, 1), 45);
  SET(m, 43, 3, 2, 29);
  TMC_EXPECT_EQ(AT(m, 43, 3, 2), 29);
  SET(m, 33, 2, 3, 79);
  TMC_EXPECT_EQ(AT(m, 33, 2, 3), 79);
  SET(m, 2, 2, 3, 64);
  TMC_EXPECT_EQ(AT(m, 2, 2, 3), 64);
  /* each write landed in its own block */
  TMC_EXPECT_EQ(AT(m, 43, 2, 3), 24);
  TMC_EXPECT_EQ(AT(m, 33, 3, 2), 24);
  TMC_EXPECT_EQ(AT(m, 3, 2, 3), 24);
  TMC_EXPECT_EQ(AT(m, 1, 2, 1), 24);

  for (int block = 0; block < 3; ++block)
    expect_single_diagonal(m, block);

  tmc_error err;
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 0, 4, 2, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 2, 1, 5, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 54, 0, 2, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 1, 2, 2, &err), err, TMC_MATRIX_ZERO_ELEMENT_ACCESS);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 2, 0, 2, &err), err, TMC_MATRIX_ZERO_ELEMENT_ACCESS);
  /* the two-argument VectorIndex needs a one-block matrix (MICM throws) */
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index_one_block(m, 0, 1, &err), err, TMC_MATRIX_MISSING_BLOCK_INDEX);
  /* the element proxies of MICM's writes */
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 0, 0, 4, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 53, 0, 0, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 0, 5, 0, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 0, 3, 3, &err), err, TMC_MATRIX_ZERO_ELEMENT_ACCESS);
  return m;
}

int sparse_policy_set_scalar(tmc_sparse_ordering ordering, int L)
{
  /* MICM's matrix has no non-zero elements, so its check runs over no values */
  tmc_sparse *m = build(ordering, L, 3, 0, NULL, 1, 0.0);
  tmc_sparse_fill(m, 2.0);
  TMC_EXPECT_EQ(m->size, (size_t)0);
  for (size_t k = 0; k < m->size; ++k)
    TMC_EXPECT_EQ(m->data[k], 2.0);
  tmc_sparse_free(m);

  /* (tighter than MICM) a matrix with values: every stored value, the
   * padding of a partial last group included, takes the scalar */
  const int elements[2][2] = { { 0, 1 }, { 2, 2 } };
  tmc_sparse *m2 = build(ordering, L, 3, 2, elements, 5, 0.0);
  tmc_sparse_fill(m2, 2.0);
  TMC_EXPECT_EQ(m2->size, (size_t)((5 + L - 1) / L * L * 2));
  for (size_t k = 0; k < m2->size; ++k)
    TMC_EXPECT_EQ(m2->data[k], 2.0);
  for (int block = 0; block < 5; ++block)
  {
    TMC_EXPECT_EQ(AT(m2, block, 0, 1), 2.0);
    TMC_EXPECT_EQ(AT(m2, block, 2, 2), 2.0);
  }
  tmc_sparse_free(m2);
  return TMC_TEST_PASS;
}

/* the pattern of single_block_elements, with (0, 1) given twice */
static const int diagonal_test_elements[6][2] = { { 0, 1 }, { 3, 2 }, { 0, 1 }, { 2, 3 }, { 1, 1 }, { 2, 1 } };

int sparse_policy_add_to_diagonal(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 4, 6, diagonal_test_elements, 3, 24);
  TMC_EXPECT_EQ(m->nnz, 5); /* the repeated element counts once */

  SET(m, 0, 1, 1, 1);
  SET(m, 1, 1, 1, 3);
  SET(m, 2, 1, 1, 5);
  tmc_sparse_add_to_diagonal(m, 7);

  TMC_EXPECT_EQ(AT(m, 0, 1, 1), 8);
  TMC_EXPECT_EQ(AT(m, 1, 1, 1), 10);
  TMC_EXPECT_EQ(AT(m, 2, 1, 1), 12);
  for (int block = 0; block < 3; ++block)
  {
    TMC_EXPECT_EQ(AT(m, block, 0, 1), 24);
    TMC_EXPECT_EQ(AT(m, block, 2, 1), 24);
    TMC_EXPECT_EQ(AT(m, block, 2, 3), 24);
    TMC_EXPECT_EQ(AT(m, block, 3, 2), 24);
  }
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

/* The text of a matrix, through the buffer and the stream functions, equals
 * `expected`. */
static void expect_text(const tmc_sparse *m, const char *expected, size_t (*format)(const tmc_sparse *, char *, size_t), int (*print)(const tmc_sparse *, FILE *))
{
  char text[1024];
  size_t length = format(m, text, sizeof text);
  TMC_EXPECT_STREQ(text, expected);
  TMC_EXPECT_EQ(length, strlen(expected));

  /* a short buffer gets the start of the text; the full length is returned */
  char short_text[9];
  TMC_EXPECT_EQ(format(m, short_text, sizeof short_text), strlen(expected));
  TMC_EXPECT_STREQ(short_text, "Block 0\n");
  TMC_EXPECT_EQ(format(m, NULL, 0), strlen(expected));

  FILE *stream = tmpfile();
  if (!stream)
  {
    tmc_test_fail(__FILE__, __LINE__, "tmpfile failed");
    return;
  }
  TMC_EXPECT_EQ(print(m, stream), (int)strlen(expected));
  rewind(stream);
  char printed[1024] = { 0 };
  size_t n = fread(printed, 1, sizeof printed - 1, stream);
  fclose(stream);
  TMC_EXPECT_EQ(n, strlen(expected));
  TMC_EXPECT_STREQ(printed, expected);
}

int sparse_policy_print_non_zero(tmc_sparse_ordering ordering, int L)
{
  const int elements[5][2] = { { 0, 1 }, { 1, 1 }, { 2, 1 }, { 2, 3 }, { 3, 2 } };
  tmc_sparse *m = build(ordering, L, 4, 5, elements, 3, 32);
  SET(m, 0, 1, 1, 2);
  SET(m, 1, 1, 1, 4);
  SET(m, 2, 1, 1, 6);
  tmc_sparse_add_to_diagonal(m, 5);

  const char *expected =
      "Block 0\n"
      "0, 1, 32\n"
      "1, 1, 7\n"
      "2, 1, 32\n"
      "2, 3, 32\n"
      "3, 2, 32\n"
      "Block 1\n"
      "0, 1, 32\n"
      "1, 1, 9\n"
      "2, 1, 32\n"
      "2, 3, 32\n"
      "3, 2, 32\n"
      "Block 2\n"
      "0, 1, 32\n"
      "1, 1, 11\n"
      "2, 1, 32\n"
      "2, 3, 32\n"
      "3, 2, 32\n";
  expect_text(m, expected, tmc_sparse_format_non_zero, tmc_sparse_print_non_zero);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

int sparse_policy_print(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 4, 6, diagonal_test_elements, 3, 24);
  SET(m, 0, 1, 1, 1);
  SET(m, 1, 1, 1, 3);
  SET(m, 2, 1, 1, 5);
  tmc_sparse_add_to_diagonal(m, 7);

  const char *expected =
      "Block 0\n"
      "0,24,0,0\n"
      "0,8,0,0\n"
      "0,24,0,24\n"
      "0,0,24,0\n"
      "Block 1\n"
      "0,24,0,0\n"
      "0,10,0,0\n"
      "0,24,0,24\n"
      "0,0,24,0\n"
      "Block 2\n"
      "0,24,0,0\n"
      "0,12,0,0\n"
      "0,24,0,24\n"
      "0,0,24,0\n";
  expect_text(m, expected, tmc_sparse_format, tmc_sparse_print);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

int sparse_policy_builder_bad_configuration(void)
{
  const int bad[3][2] = { { 3, 0 }, { 2, 4 }, { 6, 7 } };
  for (int k = 0; k < 3; ++k)
  {
    tmc_sparse_builder *b = tmc_sparse_builder_create(3);
    tmc_error err;
    tmc_error_clear(&err);
    EXPECT_MATRIX_ERROR(tmc_sparse_builder_with_element(b, bad[k][0], bad[k][1], &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
    TMC_EXPECT_EQ(tmc_sparse_builder_number_of_elements(b), 0); /* nothing was added */
    tmc_sparse_builder_free(b);
  }
  return TMC_TEST_PASS;
}

/* ---- Function ------------------------------------------------------------- */

/* X 0 0 0
 * 0 X 0 0
 * 0 0 X X
 * 0 0 0 X */
static const int array_function_elements[5][2] = { { 0, 0 }, { 1, 1 }, { 2, 2 }, { 2, 3 }, { 3, 3 } };

/* tmp = (0,0) + (1,1) + (2,2) + (2,3); (2,3) = 2 tmp; e = { e00, e11, e22, e23 } */
static void array_function_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 0, ids->e[1]), ELEM(g, 0, ids->e[2]), ELEM(g, 0, ids->e[3]), tmp };
  tmc_array_group_for_each(g, sum4_into_last, NULL, 5, sum);
  double *scale[] = { ELEM(g, 0, ids->e[3]), tmp };
  tmc_array_group_for_each(g, twice_second_into_first, NULL, 2, scale);
}

int sparse_policy_array_function(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 4, 5, array_function_elements, 3, 0.0);
  tmc_sparse_fill(m, 1.0);
  SET(m, 0, 0, 0, 1.0);
  SET(m, 1, 1, 1, 2.0);
  SET(m, 2, 2, 2, 3.0);
  SET(m, 2, 2, 3, 4.0);
  SET(m, 2, 3, 3, 5.0);

  element_ids ids = { { elem(m, 0, 0), elem(m, 1, 1), elem(m, 2, 2), elem(m, 2, 3), 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m) };
  tmc_array_function *f = make_function(array_function_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);

  TMC_EXPECT_EQ(AT(m, 0, 2, 3), 2.0 * (1.0 + 1.0 + 1.0 + 1.0));
  TMC_EXPECT_EQ(AT(m, 1, 2, 3), 2.0 * (1.0 + 2.0 + 1.0 + 1.0));
  TMC_EXPECT_EQ(AT(m, 2, 2, 3), 2.0 * (1.0 + 1.0 + 3.0 + 4.0));
  TMC_EXPECT_EQ(AT(m, 0, 0, 0), 1.0);
  TMC_EXPECT_EQ(AT(m, 1, 0, 0), 1.0);
  TMC_EXPECT_EQ(AT(m, 2, 0, 0), 1.0);
  TMC_EXPECT_EQ(AT(m, 0, 1, 1), 1.0);
  TMC_EXPECT_EQ(AT(m, 1, 1, 1), 2.0);
  TMC_EXPECT_EQ(AT(m, 2, 1, 1), 1.0);
  TMC_EXPECT_EQ(AT(m, 0, 2, 2), 1.0);
  TMC_EXPECT_EQ(AT(m, 1, 2, 2), 1.0);
  TMC_EXPECT_EQ(AT(m, 2, 2, 2), 3.0);
  TMC_EXPECT_EQ(AT(m, 0, 3, 3), 1.0);
  TMC_EXPECT_EQ(AT(m, 1, 3, 3), 1.0);
  TMC_EXPECT_EQ(AT(m, 2, 3, 3), 5.0);

  /* a different matrix with the same dimensions */
  tmc_sparse *m2 = build(ordering, L, 4, 5, array_function_elements, 3, 0.0);
  tmc_sparse_fill(m2, -1.0);
  tmc_array_arg args2[] = { tmc_array_sparse(m2) };
  EXPECT_APPLY_OK(f, args2);
  for (int block = 0; block < 3; ++block)
  {
    TMC_EXPECT_EQ(AT(m2, block, 2, 3), 2.0 * (-1.0 + -1.0 + -1.0 + -1.0));
    TMC_EXPECT_EQ(AT(m2, block, 0, 0), -1.0);
    TMC_EXPECT_EQ(AT(m2, block, 1, 1), -1.0);
    TMC_EXPECT_EQ(AT(m2, block, 2, 2), -1.0);
    TMC_EXPECT_EQ(AT(m2, block, 3, 3), -1.0);
  }

  tmc_array_function_free(f);
  tmc_sparse_free(m);
  tmc_sparse_free(m2);
  return TMC_TEST_PASS;
}

/* tmp = A(0,1) + B(2,2); A(1,2) = tmp, with A argument 0 and B argument 1;
 * e = { a01, b22, a12 } */
static void multi_matrix_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 1, ids->e[1]), tmp };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
  double *copy[] = { tmp, ELEM(g, 0, ids->e[2]) };
  tmc_array_group_for_each(g, copy_first_to_second, NULL, 2, copy);
}

int sparse_policy_multi_matrix_array_function(tmc_sparse_ordering ordering, int L)
{
  /* 0 X 0        X 0 0
   * 0 0 X        0 X 0
   * 0 0 0        0 0 X */
  const int elements_a[2][2] = { { 0, 1 }, { 1, 2 } };
  const int elements_b[3][2] = { { 0, 0 }, { 1, 1 }, { 2, 2 } };
  tmc_sparse *a = build(ordering, L, 3, 2, elements_a, 3, 0.0);
  tmc_sparse *b = build(ordering, L, 3, 3, elements_b, 3, 0.0);
  for (int block = 0; block < 3; ++block)
  {
    SET(a, block, 0, 1, block + 10);
    SET(a, block, 1, 2, block * 2 + 20);
    SET(b, block, 0, 0, block * 4);
    SET(b, block, 1, 1, block * 3);
    SET(b, block, 2, 2, block * 5);
  }

  element_ids ids = { { elem(a, 0, 1), elem(b, 2, 2), elem(a, 1, 2), 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(a), tmc_array_const_sparse(b) };
  tmc_array_function *f = make_function(multi_matrix_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);

  TMC_EXPECT_EQ(AT(a, 0, 1, 2), 10 + 0);
  TMC_EXPECT_EQ(AT(a, 1, 1, 2), 11 + 5);
  TMC_EXPECT_EQ(AT(a, 2, 1, 2), 12 + 10);
  TMC_EXPECT_EQ(AT(a, 0, 0, 1), 10.0);
  TMC_EXPECT_EQ(AT(a, 1, 0, 1), 11.0);
  TMC_EXPECT_EQ(AT(a, 2, 0, 1), 12.0);
  TMC_EXPECT_EQ(AT(b, 0, 0, 0), 0.0);
  TMC_EXPECT_EQ(AT(b, 1, 0, 0), 4.0);
  TMC_EXPECT_EQ(AT(b, 2, 0, 0), 8.0);

  tmc_array_function_free(f);
  tmc_sparse_free(a);
  tmc_sparse_free(b);
  return TMC_TEST_PASS;
}

/* A(1,1) = A(0,1) + B(0,1); e = { e01, e11 } */
static void block_sum_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 1, ids->e[0]), ELEM(g, 0, ids->e[1]) };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, ops);
}

int sparse_policy_mismatched_block_dimensions(tmc_sparse_ordering ordering, int L)
{
  const int elements[2][2] = { { 0, 1 }, { 1, 1 } };
  tmc_sparse *a = build(ordering, L, 3, 2, elements, 3, 0.0);
  tmc_sparse *b = build(ordering, L, 3, 2, elements, 4, 0.0); /* a different number of blocks */

  element_ids ids = { { elem(a, 0, 1), elem(a, 1, 1), 0, 0, 0, 0 } };
  /* creation succeeds whatever the block counts */
  tmc_array_arg create_args[] = { tmc_array_sparse(a), tmc_array_const_sparse(a) };
  tmc_array_function *f = make_function(block_sum_kernel, &ids, N_ARGS(create_args), create_args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, create_args);

  /* the application fails when the block counts differ */
  tmc_array_arg args[] = { tmc_array_sparse(a), tmc_array_const_sparse(b) };
  EXPECT_APPLY_ERROR(f, args, TMC_MATRIX_INVALID_VECTOR);

  tmc_array_function_free(f);
  tmc_sparse_free(a);
  tmc_sparse_free(b);
  return TMC_TEST_PASS;
}

/* (3,3) = 2 (0,1), for an element the matrix lacks; e = { e01, e33 } */
static void missing_element_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 0, ids->e[1]) };
  tmc_array_group_for_each(g, twice_first_into_second, NULL, 2, ops);
}

int sparse_policy_mismatched_element_dimensions(tmc_sparse_ordering ordering, int L)
{
  /* 0 X 0 0
   * 0 X 0 0
   * 0 0 X 0
   * 0 0 0 0 */
  const int elements[3][2] = { { 0, 1 }, { 1, 1 }, { 2, 2 } };
  tmc_sparse *m = build(ordering, L, 4, 3, elements, 3, 0.0);
  for (int block = 0; block < 3; ++block)
    SET(m, block, 0, 1, block + 1);

  /* MICM asserts on the missing element (EXPECT_DEATH in debug builds) when
   * its index is taken */
  tmc_error err;
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 0, 3, 3, &err), err, TMC_MATRIX_ZERO_ELEMENT_ACCESS);
  element_ids ids = { { elem(m, 0, 1), tmc_sparse_element(m, 3, 3), 0, 0, 0, 0 } };
  TMC_EXPECT_EQ(ids.e[1], -1);

  /* creating the function succeeds; applying it (which MICM's release build
   * never does) fails at the missing element, without writing */
  tmc_array_arg args[] = { tmc_array_sparse(m) };
  tmc_array_function *f = make_function(missing_element_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_ERROR(f, args, TMC_MATRIX_ZERO_ELEMENT_ACCESS);
  for (int block = 0; block < 3; ++block)
  {
    TMC_EXPECT_EQ(AT(m, block, 0, 1), block + 1);
    TMC_EXPECT_EQ(AT(m, block, 1, 1), 0.0);
    TMC_EXPECT_EQ(AT(m, block, 2, 2), 0.0);
  }

  tmc_array_function_free(f);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

/* tmp1 = (0,1) * (1,2); tmp2 = (0,1) + (1,2); (2,3) = tmp1 + tmp2;
 * (3,4) = tmp1 - tmp2; (4,4) = tmp1 * tmp2; e = { e01, e12, e23, e34, e44 } */
static void multiple_temporaries_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *tmp1 = tmc_array_group_variable(g);
  double *tmp2 = tmc_array_group_variable(g);
  double *product[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 0, ids->e[1]), tmp1 };
  tmc_array_group_for_each(g, product_into_last, NULL, 3, product);
  double *sum[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 0, ids->e[1]), tmp2 };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
  double *e23[] = { tmp1, tmp2, ELEM(g, 0, ids->e[2]) };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, e23);
  double *e34[] = { tmp1, tmp2, ELEM(g, 0, ids->e[3]) };
  tmc_array_group_for_each(g, difference_into_last, NULL, 3, e34);
  double *e44[] = { tmp1, tmp2, ELEM(g, 0, ids->e[4]) };
  tmc_array_group_for_each(g, product_into_last, NULL, 3, e44);
}

int sparse_policy_multiple_temporaries(tmc_sparse_ordering ordering, int L)
{
  /* 0 X 0 0 0
   * 0 0 X 0 0
   * 0 0 0 X 0
   * 0 0 0 0 X
   * 0 0 0 0 X */
  const int elements[5][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 4 } };
  tmc_sparse *m = build(ordering, L, 5, 5, elements, 4, 0.0);
  for (int block = 0; block < 4; ++block)
  {
    SET(m, block, 0, 1, block + 1);
    SET(m, block, 1, 2, (block + 1) * 10);
  }

  element_ids ids = { { elem(m, 0, 1), elem(m, 1, 2), elem(m, 2, 3), elem(m, 3, 4), elem(m, 4, 4), 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m) };
  tmc_array_function *f = make_function(multiple_temporaries_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);

  /* block b: (0,1) = b + 1, (1,2) = 10 (b + 1); product = 10 (b + 1)^2, sum = 11 (b + 1) */
  for (int block = 0; block < 4; ++block)
  {
    double product = 10.0 * (block + 1) * (block + 1), sum = 11.0 * (block + 1);
    TMC_EXPECT_EQ(AT(m, block, 2, 3), product + sum);
    TMC_EXPECT_EQ(AT(m, block, 3, 4), product - sum);
    TMC_EXPECT_EQ(AT(m, block, 4, 4), product * sum);
  }
  TMC_EXPECT_EQ(AT(m, 0, 2, 3), 21.0);
  TMC_EXPECT_EQ(AT(m, 1, 4, 4), 880.0);
  TMC_EXPECT_EQ(AT(m, 3, 3, 4), 116.0);

  tmc_array_function_free(f);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

/* the block lanes are looked up once and reused: (1,2) = 2 (0,1);
 * (2,3) = (0,1) + (1,2); (3,3) = (2,3) * (1,2); e = { e01, e12, e23, e33 } */
static void block_view_reuse_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *elem01 = ELEM(g, 0, ids->e[0]);
  double *elem12 = ELEM(g, 0, ids->e[1]);
  double *elem23 = ELEM(g, 0, ids->e[2]);
  double *elem33 = ELEM(g, 0, ids->e[3]);
  double *twice[] = { elem01, elem12 };
  tmc_array_group_for_each(g, twice_first_into_second, NULL, 2, twice);
  double *sum[] = { elem01, elem12, elem23 };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
  double *product[] = { elem23, elem12, elem33 };
  tmc_array_group_for_each(g, product_into_last, NULL, 3, product);
}

/* 0 X 0 0
 * 0 0 X 0
 * 0 0 0 X
 * 0 0 0 X */
static const int chain_elements[4][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 3 } };

int sparse_policy_block_view_reuse(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 4, 4, chain_elements, 3, 0.0);
  for (int block = 0; block < 3; ++block)
    SET(m, block, 0, 1, block + 1);

  element_ids ids = { { elem(m, 0, 1), elem(m, 1, 2), elem(m, 2, 3), elem(m, 3, 3), 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m) };
  tmc_array_function *f = make_function(block_view_reuse_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);

  TMC_EXPECT_EQ(AT(m, 0, 1, 2), 2.0);
  TMC_EXPECT_EQ(AT(m, 0, 2, 3), 3.0);
  TMC_EXPECT_EQ(AT(m, 0, 3, 3), 6.0);
  TMC_EXPECT_EQ(AT(m, 1, 1, 2), 4.0);
  TMC_EXPECT_EQ(AT(m, 1, 2, 3), 6.0);
  TMC_EXPECT_EQ(AT(m, 1, 3, 3), 24.0);
  TMC_EXPECT_EQ(AT(m, 2, 1, 2), 6.0);
  TMC_EXPECT_EQ(AT(m, 2, 2, 3), 9.0);
  TMC_EXPECT_EQ(AT(m, 2, 3, 3), 54.0);

  tmc_array_function_free(f);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

/* tmp = (0,0) + (1,1) + (2,2); (2,2) = 2 tmp; e = { e00, e11, e22 } */
static void function_reusability_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 0, ids->e[1]), ELEM(g, 0, ids->e[2]), tmp };
  tmc_array_group_for_each(g, sum3_into_last, NULL, 4, sum);
  double *scale[] = { ELEM(g, 0, ids->e[2]), tmp };
  tmc_array_group_for_each(g, twice_second_into_first, NULL, 2, scale);
}

int sparse_policy_function_reusability(tmc_sparse_ordering ordering, int L)
{
  const int elements[3][2] = { { 0, 0 }, { 1, 1 }, { 2, 2 } };
  tmc_sparse *m1 = build(ordering, L, 3, 3, elements, 2, 0.0);
  element_ids ids = { { elem(m1, 0, 0), elem(m1, 1, 1), elem(m1, 2, 2), 0, 0, 0 } };
  tmc_array_arg args1[] = { tmc_array_sparse(m1) };
  tmc_array_function *f = make_function(function_reusability_kernel, &ids, N_ARGS(args1), args1, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;

  for (int block = 0; block < 2; ++block)
  {
    SET(m1, block, 0, 0, block);
    SET(m1, block, 1, 1, block + 1);
    SET(m1, block, 2, 2, block + 2);
  }
  EXPECT_APPLY_OK(f, args1);
  TMC_EXPECT_EQ(AT(m1, 0, 2, 2), 2.0 * (0 + 1 + 2));
  TMC_EXPECT_EQ(AT(m1, 1, 2, 2), 2.0 * (1 + 2 + 3));

  tmc_sparse *m2 = build(ordering, L, 3, 3, elements, 2, 0.0);
  tmc_sparse_fill(m2, 5.0);
  tmc_array_arg args2[] = { tmc_array_sparse(m2) };
  EXPECT_APPLY_OK(f, args2);
  TMC_EXPECT_EQ(AT(m2, 0, 2, 2), 2.0 * (5 + 5 + 5));
  TMC_EXPECT_EQ(AT(m2, 1, 2, 2), 2.0 * (5 + 5 + 5));

  tmc_sparse *m3 = build(ordering, L, 3, 3, elements, 2, 0.0);
  tmc_sparse_fill(m3, 0.0);
  for (int block = 0; block < 2; ++block)
    SET(m3, block, 0, 0, block * 10);
  tmc_array_arg args3[] = { tmc_array_sparse(m3) };
  EXPECT_APPLY_OK(f, args3);
  TMC_EXPECT_EQ(AT(m3, 0, 2, 2), 2.0 * (0 + 0 + 0));
  TMC_EXPECT_EQ(AT(m3, 1, 2, 2), 2.0 * (10 + 0 + 0));

  tmc_array_function_free(f);
  tmc_sparse_free(m1);
  tmc_sparse_free(m2);
  tmc_sparse_free(m3);
  return TMC_TEST_PASS;
}

/* tmp = A(0,1) + A(1,1) + B(0,0); A(2,2) = tmp; e = { a01, a11, b00, a22 } */
static void different_structure_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 0, ids->e[1]), ELEM(g, 1, ids->e[2]), tmp };
  tmc_array_group_for_each(g, sum3_into_last, NULL, 4, sum);
  double *copy[] = { tmp, ELEM(g, 0, ids->e[3]) };
  tmc_array_group_for_each(g, copy_first_to_second, NULL, 2, copy);
}

int sparse_policy_two_sparse_matrices_different_structure(tmc_sparse_ordering ordering, int L)
{
  /* 0 X 0        X 0 0 0 0
   * 0 X 0        0 0 X 0 0
   * 0 0 X        0 0 0 X 0
   *              0 0 0 0 X
   *              0 0 0 0 0 */
  const int elements_a[3][2] = { { 0, 1 }, { 1, 1 }, { 2, 2 } };
  const int elements_b[4][2] = { { 0, 0 }, { 1, 2 }, { 2, 3 }, { 3, 4 } };
  tmc_sparse *a = build(ordering, L, 3, 3, elements_a, 4, 0.0);
  tmc_sparse *b = build(ordering, L, 5, 4, elements_b, 4, 0.0);
  for (int block = 0; block < 4; ++block)
  {
    SET(a, block, 0, 1, block * 2 + 1);
    SET(a, block, 1, 1, block * 3 + 2);
    SET(a, block, 2, 2, block * 5 + 3);
    SET(b, block, 0, 0, block + 10);
    SET(b, block, 1, 2, block + 20);
    SET(b, block, 2, 3, block + 30);
    SET(b, block, 3, 4, block + 40);
  }

  element_ids ids = { { elem(a, 0, 1), elem(a, 1, 1), elem(b, 0, 0), elem(a, 2, 2), 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(a), tmc_array_const_sparse(b) };
  tmc_array_function *f = make_function(different_structure_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);

  TMC_EXPECT_EQ(AT(a, 0, 2, 2), 1 + 2 + 10);
  TMC_EXPECT_EQ(AT(a, 1, 2, 2), 3 + 5 + 11);
  TMC_EXPECT_EQ(AT(a, 2, 2, 2), 5 + 8 + 12);
  TMC_EXPECT_EQ(AT(a, 3, 2, 2), 7 + 11 + 13);
  TMC_EXPECT_EQ(AT(a, 0, 0, 1), 1.0);
  TMC_EXPECT_EQ(AT(a, 1, 0, 1), 3.0);
  TMC_EXPECT_EQ(AT(a, 2, 0, 1), 5.0);
  TMC_EXPECT_EQ(AT(a, 3, 0, 1), 7.0);
  TMC_EXPECT_EQ(AT(b, 0, 0, 0), 10.0);
  TMC_EXPECT_EQ(AT(b, 1, 1, 2), 21.0);
  TMC_EXPECT_EQ(AT(b, 2, 2, 3), 32.0);
  TMC_EXPECT_EQ(AT(b, 3, 3, 4), 43.0);

  tmc_array_function_free(f);
  tmc_sparse_free(a);
  tmc_sparse_free(b);
  return TMC_TEST_PASS;
}

/* tmp = s(0,1) + s(1,2) + dense[0] + dense[1]; s(3,3) = tmp, with the sparse
 * matrix argument 0 and the dense matrix argument 1; e = { s01, s12, s33 } */
static void sparse_and_dense_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 0, ids->e[1]), COL(g, 1, 0), COL(g, 1, 1), tmp };
  tmc_array_group_for_each(g, sum4_into_last, NULL, 5, sum);
  double *copy[] = { tmp, ELEM(g, 0, ids->e[2]) };
  tmc_array_group_for_each(g, copy_first_to_second, NULL, 2, copy);
}

int sparse_policy_sparse_and_dense_matrix_function(tmc_sparse_ordering ordering, int L, int dense_L)
{
  tmc_sparse *sparse = build(ordering, L, 4, 4, chain_elements, 3, 0.0);
  /* 3 rows, one per block of the sparse matrix, by 4 columns */
  tmc_dense *dense = tmc_dense_create(3, 4, dense_L, 0.0);
  for (int block = 0; block < 3; ++block)
  {
    SET(sparse, block, 0, 1, block + 1);
    SET(sparse, block, 1, 2, block + 10);
    SET(sparse, block, 2, 3, block + 20);
    SET(sparse, block, 3, 3, block + 30);
  }
  for (int row = 0; row < 3; ++row)
    for (int col = 0; col < 4; ++col)
      tmc_dense_set(dense, row, col, row * 100 + col * 10);

  element_ids ids = { { elem(sparse, 0, 1), elem(sparse, 1, 2), elem(sparse, 3, 3), 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(sparse), tmc_array_const_dense(dense) };
  tmc_array_function *f = make_function(sparse_and_dense_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);

  TMC_EXPECT_EQ(AT(sparse, 0, 3, 3), 1 + 10 + 0 + 10);
  TMC_EXPECT_EQ(AT(sparse, 1, 3, 3), 2 + 11 + 100 + 110);
  TMC_EXPECT_EQ(AT(sparse, 2, 3, 3), 3 + 12 + 200 + 210);
  TMC_EXPECT_EQ(AT(sparse, 0, 0, 1), 1.0);
  TMC_EXPECT_EQ(AT(sparse, 1, 1, 2), 11.0);
  TMC_EXPECT_EQ(AT(sparse, 2, 2, 3), 22.0);
  TMC_EXPECT_EQ(tmc_dense_get(dense, 0, 0), 0.0);
  TMC_EXPECT_EQ(tmc_dense_get(dense, 1, 2), 120.0);
  TMC_EXPECT_EQ(tmc_dense_get(dense, 2, 3), 230.0);

  tmc_array_function_free(f);
  tmc_sparse_free(sparse);
  tmc_dense_free(dense);
  return TMC_TEST_PASS;
}

static void no_op_kernel(tmc_array_group *g, void *ctx)
{
  (void)g;
  (void)ctx;
}

/* A sparse matrix of width L and a dense matrix of width dense_L cannot share
 * a function. */
static int incompatible_dense(tmc_sparse_ordering ordering, int L, int dense_L)
{
  const int elements[2][2] = { { 0, 1 }, { 1, 2 } };
  tmc_sparse *sparse = build(ordering, L, 4, 2, elements, 3, 0.0);
  tmc_dense *dense = tmc_dense_create(3, 4, dense_L, 0.0);
  tmc_array_arg args[] = { tmc_array_const_sparse(sparse), tmc_array_const_dense(dense) };
  expect_create_error(no_op_kernel, N_ARGS(args), args, TMC_MATRIX_INVALID_VECTOR, ordering, L);
  tmc_sparse_free(sparse);
  tmc_dense_free(dense);
  return TMC_TEST_PASS;
}

int sparse_policy_incompatible_ordering(tmc_sparse_ordering ordering, int L)
{
  /* MICM runs this only for L > 1: a vector-ordered sparse matrix with a
   * standard-ordered (L = 1) dense matrix */
  if (L <= 1)
    return TMC_TEST_PASS;
  return incompatible_dense(ordering, L, 1);
}

int sparse_policy_incompatible_vector_ordering(tmc_sparse_ordering ordering, int L, int dense_L)
{
  if (L == dense_L)
    return TMC_TEST_PASS;
  return incompatible_dense(ordering, L, dense_L);
}

int sparse_policy_incompatible_sparse_ordering(tmc_sparse_ordering ordering, int L, int other_L)
{
  if (L == other_L)
    return TMC_TEST_PASS;
  const int elements[2][2] = { { 0, 1 }, { 1, 2 } };
  tmc_sparse *m1 = build(ordering, L, 4, 2, elements, 3, 0.0);
  tmc_sparse *m2 = build(ordering, other_L, 4, 2, elements, 3, 0.0);
  tmc_array_arg args[] = { tmc_array_const_sparse(m1), tmc_array_const_sparse(m2) };
  expect_create_error(no_op_kernel, N_ARGS(args), args, TMC_MATRIX_INVALID_VECTOR, ordering, L);
  tmc_sparse_free(m1);
  tmc_sparse_free(m2);
  return TMC_TEST_PASS;
}

/* reads only: tmp = (0,1) + (1,2); a local sum over (2,3); e = { e01, e12, e23 } */
static void const_matrix_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 0, ids->e[1]), tmp };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
  double local_sum = 0.0;
  double *e23[] = { ELEM(g, 0, ids->e[2]) };
  tmc_array_group_for_each(g, accumulate_first, &local_sum, 1, e23);
}

int sparse_policy_const_sparse_matrix_function(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 4, 4, chain_elements, 3, 0.0);
  for (int block = 0; block < 3; ++block)
  {
    SET(m, block, 0, 1, block + 1);
    SET(m, block, 1, 2, block + 10);
  }
  const tmc_sparse *const_matrix = m;
  element_ids ids = { { elem(m, 0, 1), elem(m, 1, 2), elem(m, 2, 3), 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_const_sparse(const_matrix) };
  tmc_array_function *f = make_function(const_matrix_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 0, 1), 1.0);
  TMC_EXPECT_EQ(AT(m, 1, 1, 2), 11.0);
  tmc_array_function_free(f);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

/* takes a block variable and counts its calls; touches no element */
static void variable_only_kernel(tmc_array_group *g, void *ctx)
{
  ++*(int *)ctx;
  double *tmp = tmc_array_group_variable(g);
  (void)tmp;
}

static void counting_kernel(tmc_array_group *g, void *ctx)
{
  (void)g;
  ++*(int *)ctx;
}

int sparse_policy_empty_sparse_matrix_function(tmc_sparse_ordering ordering, int L)
{
  /* no non-zero elements: the kernel runs once per group of blocks and has
   * nothing to access */
  tmc_sparse *empty_matrix = build(ordering, L, 4, 0, NULL, 3, 0.0);
  int calls = 0;
  tmc_array_arg args[] = { tmc_array_const_sparse(empty_matrix) };
  tmc_array_function *f = make_function(variable_only_kernel, &calls, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(calls, (3 + L - 1) / L);

  /* no blocks: the kernel never runs */
  const int elements[1][2] = { { 0, 1 } };
  tmc_sparse *zero_blocks = build(ordering, L, 4, 1, elements, 0, 0.0);
  TMC_EXPECT_EQ(zero_blocks->n_blocks, 0);
  int calls2 = 0;
  tmc_array_arg args2[] = { tmc_array_const_sparse(zero_blocks) };
  tmc_array_function *f2 = make_function(counting_kernel, &calls2, N_ARGS(args2), args2, ordering, L);
  if (!f2)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f2, args2);
  TMC_EXPECT_EQ(calls2, 0);

  tmc_array_function_free(f);
  tmc_array_function_free(f2);
  tmc_sparse_free(empty_matrix);
  tmc_sparse_free(zero_blocks);
  return TMC_TEST_PASS;
}

/* 0 X 0 0
 * 0 0 X 0
 * 0 0 0 X
 * 0 0 0 0 */
static const int superdiagonal_elements[3][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 } };

/* tmp = A(0,1) + A(1,2) + B(2,3); A(2,3) = tmp; e = { e01, e12, e23 } */
static void different_blocks_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 0, ids->e[1]), ELEM(g, 1, ids->e[2]), tmp };
  tmc_array_group_for_each(g, sum3_into_last, NULL, 4, sum);
  double *copy[] = { tmp, ELEM(g, 0, ids->e[2]) };
  tmc_array_group_for_each(g, copy_first_to_second, NULL, 2, copy);
}

int sparse_policy_multiple_sparse_matrices_different_blocks_from_creation(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *a3 = build(ordering, L, 4, 3, superdiagonal_elements, 3, 0.0);
  tmc_sparse *b3 = build(ordering, L, 4, 3, superdiagonal_elements, 3, 0.0);
  for (int block = 0; block < 3; ++block)
  {
    SET(a3, block, 0, 1, block + 1);
    SET(a3, block, 1, 2, block + 10);
    SET(b3, block, 2, 3, block + 100);
  }
  element_ids ids = { { elem(a3, 0, 1), elem(a3, 1, 2), elem(a3, 2, 3), 0, 0, 0 } };
  tmc_array_arg create_args[] = { tmc_array_sparse(a3), tmc_array_const_sparse(b3) };
  tmc_array_function *f = make_function(different_blocks_kernel, &ids, N_ARGS(create_args), create_args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;

  /* four blocks, unlike the matrices the function was created with */
  tmc_sparse *a4 = build(ordering, L, 4, 3, superdiagonal_elements, 4, 0.0);
  tmc_sparse *b4 = build(ordering, L, 4, 3, superdiagonal_elements, 4, 0.0);
  for (int block = 0; block < 4; ++block)
  {
    SET(a4, block, 0, 1, block + 1);
    SET(a4, block, 1, 2, block + 10);
    SET(b4, block, 2, 3, block + 100);
  }
  tmc_array_arg args[] = { tmc_array_sparse(a4), tmc_array_const_sparse(b4) };
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(a4, 0, 2, 3), 1 + 10 + 100);
  TMC_EXPECT_EQ(AT(a4, 1, 2, 3), 2 + 11 + 101);
  TMC_EXPECT_EQ(AT(a4, 2, 2, 3), 3 + 12 + 102);
  TMC_EXPECT_EQ(AT(a4, 3, 2, 3), 4 + 13 + 103);

  tmc_array_function_free(f);
  tmc_sparse_free(a3);
  tmc_sparse_free(b3);
  tmc_sparse_free(a4);
  tmc_sparse_free(b4);
  return TMC_TEST_PASS;
}

/* tmp = (0,1) + v; (1,2) = tmp, with the vector argument 1; e = { e01, e12 } */
static void element_plus_vector_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { ELEM(g, 0, ids->e[0]), VEC(g, 1), tmp };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, sum);
  double *copy[] = { tmp, ELEM(g, 0, ids->e[1]) };
  tmc_array_group_for_each(g, copy_first_to_second, NULL, 2, copy);
}

/* 0 X 0
 * 0 0 X
 * 0 0 0 */
static const int two_element_pattern[2][2] = { { 0, 1 }, { 1, 2 } };

int sparse_policy_sparse_matrix_vector_different_blocks_from_creation(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m3 = build(ordering, L, 3, 2, two_element_pattern, 3, 0.0);
  const double values3[] = { 1.0, 2.0, 3.0 };
  tmc_padded_vector *v3 = vector_of(L, 3, values3);
  element_ids ids = { { elem(m3, 0, 1), elem(m3, 1, 2), 0, 0, 0, 0 } };
  tmc_array_arg create_args[] = { tmc_array_sparse(m3), tmc_array_const_vector(v3) };
  tmc_array_function *f = make_function(element_plus_vector_kernel, &ids, N_ARGS(create_args), create_args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;

  /* five blocks and five values, unlike at creation */
  tmc_sparse *m5 = build(ordering, L, 3, 2, two_element_pattern, 5, 0.0);
  const double values5[] = { 10.0, 20.0, 30.0, 40.0, 50.0 };
  tmc_padded_vector *v5 = vector_of(L, 5, values5);
  for (int block = 0; block < 5; ++block)
    SET(m5, block, 0, 1, block + 1);
  tmc_array_arg args[] = { tmc_array_sparse(m5), tmc_array_const_vector(v5) };
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m5, 0, 1, 2), 1 + 10);
  TMC_EXPECT_EQ(AT(m5, 1, 1, 2), 2 + 20);
  TMC_EXPECT_EQ(AT(m5, 2, 1, 2), 3 + 30);
  TMC_EXPECT_EQ(AT(m5, 3, 1, 2), 4 + 40);
  TMC_EXPECT_EQ(AT(m5, 4, 1, 2), 5 + 50);

  tmc_array_function_free(f);
  tmc_sparse_free(m3);
  tmc_sparse_free(m5);
  tmc_padded_vector_free(v3);
  tmc_padded_vector_free(v5);
  return TMC_TEST_PASS;
}

/* A(1,2) = 2 B(0,1), with A argument 0 and B argument 1; e = { e01, e12 } */
static void twice_other_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { ELEM(g, 1, ids->e[0]), ELEM(g, 0, ids->e[1]) };
  tmc_array_group_for_each(g, twice_first_into_second, NULL, 2, ops);
}

int sparse_policy_mismatched_blocks_at_invocation(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m3 = build(ordering, L, 3, 2, two_element_pattern, 3, 0.0);
  tmc_sparse *m4 = build(ordering, L, 3, 2, two_element_pattern, 4, 0.0);
  element_ids ids = { { elem(m3, 0, 1), elem(m3, 1, 2), 0, 0, 0, 0 } };
  tmc_array_arg args33[] = { tmc_array_sparse(m3), tmc_array_const_sparse(m3) };
  tmc_array_function *f = make_function(twice_other_kernel, &ids, N_ARGS(args33), args33, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args33);

  tmc_array_arg args34[] = { tmc_array_sparse(m3), tmc_array_const_sparse(m4) };
  EXPECT_APPLY_ERROR(f, args34, TMC_MATRIX_INVALID_VECTOR);
  tmc_array_arg args43[] = { tmc_array_sparse(m4), tmc_array_const_sparse(m3) };
  EXPECT_APPLY_ERROR(f, args43, TMC_MATRIX_INVALID_VECTOR);

  tmc_array_function_free(f);
  tmc_sparse_free(m3);
  tmc_sparse_free(m4);
  return TMC_TEST_PASS;
}

/* A(2,3) = A(0,1) + B(1,2) + C(2,3), with A, B, C arguments 0, 1, 2;
 * e = { e01, e12, e23 } */
static void three_matrix_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 1, ids->e[1]), ELEM(g, 2, ids->e[2]), ELEM(g, 0, ids->e[2]) };
  tmc_array_group_for_each(g, sum3_into_last, NULL, 4, ops);
}

int sparse_policy_multiple_sparse_matrices_mismatched_blocks_at_invocation(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *a3 = build(ordering, L, 4, 3, superdiagonal_elements, 3, 0.0);
  tmc_sparse *b3 = build(ordering, L, 4, 3, superdiagonal_elements, 3, 0.0);
  tmc_sparse *c4 = build(ordering, L, 4, 3, superdiagonal_elements, 4, 0.0);
  element_ids ids = { { elem(a3, 0, 1), elem(a3, 1, 2), elem(a3, 2, 3), 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(a3), tmc_array_const_sparse(b3), tmc_array_const_sparse(b3) };
  tmc_array_function *f = make_function(three_matrix_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);

  tmc_array_arg args_abc[] = { tmc_array_sparse(a3), tmc_array_const_sparse(b3), tmc_array_const_sparse(c4) };
  EXPECT_APPLY_ERROR(f, args_abc, TMC_MATRIX_INVALID_VECTOR);
  tmc_array_arg args_cbb[] = { tmc_array_sparse(c4), tmc_array_const_sparse(b3), tmc_array_const_sparse(b3) };
  EXPECT_APPLY_ERROR(f, args_cbb, TMC_MATRIX_INVALID_VECTOR);

  tmc_array_function_free(f);
  tmc_sparse_free(a3);
  tmc_sparse_free(b3);
  tmc_sparse_free(c4);
  return TMC_TEST_PASS;
}

/* ---- vectors -------------------------------------------------------------- */

/* 0 X 0
 * 0 0 0
 * 0 0 0 */
static const int one_element_pattern[1][2] = { { 0, 1 } };

/* (0,1) = 2 v, with the vector argument 1; e = { e01 } */
static void vector_doubled_into_element_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { VEC(g, 1), ELEM(g, 0, ids->e[0]) };
  tmc_array_group_for_each(g, twice_first_into_second, NULL, 2, ops);
}

int sparse_policy_vector_in_sparse_matrix_function(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 4, 3, superdiagonal_elements, 3, 0.0);
  const double values[] = { 5.0, 10.0, 15.0 };
  tmc_padded_vector *v = vector_of(L, 3, values);
  for (int block = 0; block < 3; ++block)
  {
    SET(m, block, 0, 1, block + 1);
    SET(m, block, 1, 2, block + 2);
    SET(m, block, 2, 3, block + 3);
  }
  element_ids ids = { { elem(m, 0, 1), 0, 0, 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(vector_doubled_into_element_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 0, 1), 10.0);
  TMC_EXPECT_EQ(AT(m, 1, 0, 1), 20.0);
  TMC_EXPECT_EQ(AT(m, 2, 0, 1), 30.0);
  /* the other elements are untouched */
  for (int block = 0; block < 3; ++block)
  {
    TMC_EXPECT_EQ(AT(m, block, 1, 2), block + 2);
    TMC_EXPECT_EQ(AT(m, block, 2, 3), block + 3);
  }
  tmc_array_function_free(f);
  tmc_sparse_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

/* (0,1) = v, with the vector argument 1; e = { e01 } */
static void vector_into_element_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { VEC(g, 1), ELEM(g, 0, ids->e[0]) };
  tmc_array_group_for_each(g, copy_first_to_second, NULL, 2, ops);
}

/* A function of a 3 x 3 matrix of n_blocks blocks and a vector of `size`
 * values: creation succeeds (block counts may differ then), and the
 * application returns rc. */
static int apply_vector_into_element(tmc_sparse_ordering ordering, int L, int n_blocks, int size, tmc_error *err)
{
  tmc_sparse *m = build(ordering, L, 3, 1, one_element_pattern, n_blocks, 0.0);
  tmc_padded_vector *v = tmc_padded_vector_create(size, L, 0.0);
  element_ids ids = { { elem(m, 0, 1), 0, 0, 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m), tmc_array_const_vector(v) };
  int rc = TMC_ERROR_INTERNAL;
  tmc_array_function *f = make_function(vector_into_element_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (f)
    rc = tmc_array_function_apply(f, N_ARGS(args), args, err);
  tmc_array_function_free(f);
  tmc_sparse_free(m);
  tmc_padded_vector_free(v);
  return rc;
}

static int vector_size_mismatch(tmc_sparse_ordering ordering, int L, int n_blocks, int size)
{
  tmc_error err;
  tmc_error_clear(&err);
  EXPECT_MATRIX_ERROR(apply_vector_into_element(ordering, L, n_blocks, size, &err), err, TMC_MATRIX_INVALID_VECTOR);
  return TMC_TEST_PASS;
}

int sparse_policy_vector_too_small(tmc_sparse_ordering ordering, int L)
{
  return vector_size_mismatch(ordering, L, 3, 2);
}

int sparse_policy_vector_too_large(tmc_sparse_ordering ordering, int L)
{
  return vector_size_mismatch(ordering, L, 3, 4);
}

int sparse_policy_empty_vector_non_empty_sparse_matrix(tmc_sparse_ordering ordering, int L)
{
  return vector_size_mismatch(ordering, L, 3, 0);
}

int sparse_policy_non_empty_vector_empty_sparse_matrix(tmc_sparse_ordering ordering, int L)
{
  /* MICM asserts on the index of an element in block 0 of a matrix with no
   * blocks (EXPECT_DEATH in debug builds) */
  tmc_sparse *m = build(ordering, L, 3, 1, one_element_pattern, 0, 0.0);
  tmc_error err;
  EXPECT_INDEX_ERROR(tmc_sparse_vector_index(m, 0, 0, 1, &err), err, TMC_MATRIX_ELEMENT_OUT_OF_RANGE);
  TMC_EXPECT_EQ(tmc_sparse_element(m, 0, 1), 0); /* the pattern itself has the element */
  tmc_sparse_free(m);
  /* the function can be created; applying it (which MICM's release build
   * never does) fails on the vector's size */
  return vector_size_mismatch(ordering, L, 0, 3);
}

int sparse_policy_empty_vector_empty_sparse_matrix(tmc_sparse_ordering ordering, int L)
{
  /* MICM's lambda asserts it is never called */
  tmc_sparse *m = build(ordering, L, 3, 1, one_element_pattern, 0, 0.0);
  tmc_padded_vector *v = tmc_padded_vector_create(0, L, 0.0);
  int calls = 0;
  tmc_array_arg args[] = { tmc_array_const_sparse(m), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(counting_kernel, &calls, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(calls, 0);
  tmc_array_function_free(f);
  tmc_sparse_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

/* (0,1) = v1 + v2, the vectors arguments 1 and 2; e = { e01 } */
static void two_vectors_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { VEC(g, 1), VEC(g, 2), ELEM(g, 0, ids->e[0]) };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, ops);
}

int sparse_policy_multiple_vectors_different_sizes(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  const double values1[] = { 1.0, 2.0, 3.0 }, values2[] = { 4.0, 5.0 };
  tmc_padded_vector *v1 = vector_of(L, 3, values1);
  tmc_padded_vector *v2 = vector_of(L, 2, values2);
  element_ids ids = { { elem(m, 0, 1), 0, 0, 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m), tmc_array_const_vector(v1), tmc_array_const_vector(v2) };
  tmc_array_function *f = make_function(two_vectors_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_ERROR(f, args, TMC_MATRIX_INVALID_VECTOR);
  tmc_array_function_free(f);
  tmc_sparse_free(m);
  tmc_padded_vector_free(v1);
  tmc_padded_vector_free(v2);
  return TMC_TEST_PASS;
}

int sparse_policy_multiple_vectors_same_size(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  const double values1[] = { 1.0, 2.0, 3.0 }, values2[] = { 4.0, 5.0, 6.0 };
  tmc_padded_vector *v1 = vector_of(L, 3, values1);
  tmc_padded_vector *v2 = vector_of(L, 3, values2);
  element_ids ids = { { elem(m, 0, 1), 0, 0, 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m), tmc_array_const_vector(v1), tmc_array_const_vector(v2) };
  tmc_array_function *f = make_function(two_vectors_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 0, 1), 5.0);
  TMC_EXPECT_EQ(AT(m, 1, 0, 1), 7.0);
  TMC_EXPECT_EQ(AT(m, 2, 0, 1), 9.0);
  tmc_array_function_free(f);
  tmc_sparse_free(m);
  tmc_padded_vector_free(v1);
  tmc_padded_vector_free(v2);
  return TMC_TEST_PASS;
}

/* B(0,1) = A(0,1) + v, with A, B, v arguments 0, 1, 2; e = { e01 } */
static void two_matrices_one_vector_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { ELEM(g, 0, ids->e[0]), VEC(g, 2), ELEM(g, 1, ids->e[0]) };
  tmc_array_group_for_each(g, sum2_into_last, NULL, 3, ops);
}

int sparse_policy_multiple_sparse_matrices_one_vector(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m1 = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  tmc_sparse *m2 = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  const double values[] = { 10.0, 20.0, 30.0 };
  tmc_padded_vector *v = vector_of(L, 3, values);
  for (int block = 0; block < 3; ++block)
  {
    SET(m1, block, 0, 1, block + 1);
    SET(m2, block, 0, 1, 0.0);
  }
  element_ids ids = { { elem(m1, 0, 1), 0, 0, 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_const_sparse(m1), tmc_array_sparse(m2), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(two_matrices_one_vector_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m2, 0, 0, 1), 11.0);
  TMC_EXPECT_EQ(AT(m2, 1, 0, 1), 22.0);
  TMC_EXPECT_EQ(AT(m2, 2, 0, 1), 33.0);
  tmc_array_function_free(f);
  tmc_sparse_free(m1);
  tmc_sparse_free(m2);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

/* A(0,1) = v, with A, B, v arguments 0, 1, 2; e = { e01 } */
static void vector_into_first_matrix_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { VEC(g, 2), ELEM(g, 0, ids->e[0]) };
  tmc_array_group_for_each(g, copy_first_to_second, NULL, 2, ops);
}

int sparse_policy_multiple_sparse_matrices_different_blocks_vector(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m1 = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  tmc_sparse *m2 = build(ordering, L, 3, 1, one_element_pattern, 4, 0.0);
  const double values[] = { 1.0, 2.0, 3.0 };
  tmc_padded_vector *v = vector_of(L, 3, values);
  element_ids ids = { { elem(m1, 0, 1), 0, 0, 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m1), tmc_array_const_sparse(m2), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(vector_into_first_matrix_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_ERROR(f, args, TMC_MATRIX_INVALID_VECTOR);
  tmc_array_function_free(f);
  tmc_sparse_free(m1);
  tmc_sparse_free(m2);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

int sparse_policy_vector_size_matches_one_sparse_matrix_only(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m1 = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  tmc_sparse *m2 = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  const double values1[] = { 1.0, 2.0, 3.0 }, values2[] = { 4.0, 5.0 };
  tmc_padded_vector *v1 = vector_of(L, 3, values1);
  tmc_padded_vector *v2 = vector_of(L, 2, values2);
  element_ids ids = { { elem(m1, 0, 1), 0, 0, 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m1), tmc_array_const_sparse(m2), tmc_array_const_vector(v1), tmc_array_const_vector(v2) };
  tmc_array_function *f = make_function(vector_into_first_matrix_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_ERROR(f, args, TMC_MATRIX_INVALID_VECTOR);
  tmc_array_function_free(f);
  tmc_sparse_free(m1);
  tmc_sparse_free(m2);
  tmc_padded_vector_free(v1);
  tmc_padded_vector_free(v2);
  return TMC_TEST_PASS;
}

/* (0,1) = 3 v; e = { e01 } */
static void vector_tripled_into_element_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { VEC(g, 1), ELEM(g, 0, ids->e[0]) };
  tmc_array_group_for_each(g, thrice_first_into_second, NULL, 2, ops);
}

int sparse_policy_const_vector_sparse(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  const double values[] = { 5.0, 10.0, 15.0 };
  tmc_padded_vector *storage = vector_of(L, 3, values);
  const tmc_padded_vector *v = storage;
  element_ids ids = { { elem(m, 0, 1), 0, 0, 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(vector_tripled_into_element_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 0, 1), 15.0);
  TMC_EXPECT_EQ(AT(m, 1, 0, 1), 30.0);
  TMC_EXPECT_EQ(AT(m, 2, 0, 1), 45.0);
  tmc_array_function_free(f);
  tmc_sparse_free(m);
  tmc_padded_vector_free(storage);
  return TMC_TEST_PASS;
}

/* v = 3 (0,1): the matrix is read, the vector written; e = { e01 } */
static void element_tripled_into_vector_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { ELEM(g, 0, ids->e[0]), VEC(g, 1) };
  tmc_array_group_for_each(g, thrice_first_into_second, NULL, 2, ops);
}

int sparse_policy_mutable_vector_sparse(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  const double values[] = { 5.0, 10.0, 15.0 };
  tmc_padded_vector *v = vector_of(L, 3, values);
  for (int block = 0; block < 3; ++block)
    SET(m, block, 0, 1, block + 1);
  element_ids ids = { { elem(m, 0, 1), 0, 0, 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_const_sparse(m), tmc_array_vector(v) };
  tmc_array_function *f = make_function(element_tripled_into_vector_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(v->data[0], 3.0);
  TMC_EXPECT_EQ(v->data[1], 6.0);
  TMC_EXPECT_EQ(v->data[2], 9.0);
  tmc_array_function_free(f);
  tmc_sparse_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

int sparse_policy_function_reusability_with_vectors_sparse(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m1 = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  tmc_sparse *m2 = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  const double values1[] = { 1.0, 2.0, 3.0 }, values2[] = { 10.0, 20.0, 30.0 };
  tmc_padded_vector *v1 = vector_of(L, 3, values1);
  tmc_padded_vector *v2 = vector_of(L, 3, values2);
  element_ids ids = { { elem(m1, 0, 1), 0, 0, 0, 0, 0 } };
  tmc_array_arg args1[] = { tmc_array_sparse(m1), tmc_array_const_vector(v1) };
  tmc_array_function *f = make_function(vector_doubled_into_element_kernel, &ids, N_ARGS(args1), args1, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args1);
  TMC_EXPECT_EQ(AT(m1, 0, 0, 1), 2.0);
  TMC_EXPECT_EQ(AT(m1, 1, 0, 1), 4.0);
  TMC_EXPECT_EQ(AT(m1, 2, 0, 1), 6.0);

  tmc_array_arg args2[] = { tmc_array_sparse(m2), tmc_array_const_vector(v2) };
  EXPECT_APPLY_OK(f, args2);
  TMC_EXPECT_EQ(AT(m2, 0, 0, 1), 20.0);
  TMC_EXPECT_EQ(AT(m2, 1, 0, 1), 40.0);
  TMC_EXPECT_EQ(AT(m2, 2, 0, 1), 60.0);

  tmc_array_function_free(f);
  tmc_sparse_free(m1);
  tmc_sparse_free(m2);
  tmc_padded_vector_free(v1);
  tmc_padded_vector_free(v2);
  return TMC_TEST_PASS;
}

int sparse_policy_function_invocation_with_wrong_sized_vector_sparse(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  const double values1[] = { 1.0, 2.0, 3.0 }, values2[] = { 1.0, 2.0 };
  tmc_padded_vector *v1 = vector_of(L, 3, values1);
  tmc_padded_vector *v2 = vector_of(L, 2, values2);
  tmc_array_arg args[] = { tmc_array_const_sparse(m), tmc_array_const_vector(v1) };
  tmc_array_function *f = make_function(no_op_kernel, NULL, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  tmc_array_arg wrong[] = { tmc_array_const_sparse(m), tmc_array_const_vector(v2) };
  EXPECT_APPLY_ERROR(f, wrong, TMC_MATRIX_INVALID_VECTOR);
  tmc_array_function_free(f);
  tmc_sparse_free(m);
  tmc_padded_vector_free(v1);
  tmc_padded_vector_free(v2);
  return TMC_TEST_PASS;
}

/* tmp = (0,1) + (1,2) + v; (2,3) = 2 tmp; e = { e01, e12, e23 } */
static void mixed_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *tmp = tmc_array_group_variable(g);
  double *sum[] = { ELEM(g, 0, ids->e[0]), ELEM(g, 0, ids->e[1]), VEC(g, 1), tmp };
  tmc_array_group_for_each(g, sum3_into_last, NULL, 4, sum);
  double *scale[] = { tmp, ELEM(g, 0, ids->e[2]) };
  tmc_array_group_for_each(g, twice_first_into_second, NULL, 2, scale);
}

int sparse_policy_mixed_vector_block_view_block_variable(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 4, 3, superdiagonal_elements, 3, 0.0);
  const double values[] = { 100.0, 200.0, 300.0 };
  tmc_padded_vector *v = vector_of(L, 3, values);
  for (int block = 0; block < 3; ++block)
  {
    SET(m, block, 0, 1, block + 1);
    SET(m, block, 1, 2, block + 2);
    SET(m, block, 2, 3, block + 3);
  }
  element_ids ids = { { elem(m, 0, 1), elem(m, 1, 2), elem(m, 2, 3), 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(mixed_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_EQ(AT(m, 0, 2, 3), 206.0);
  TMC_EXPECT_EQ(AT(m, 1, 2, 3), 410.0);
  TMC_EXPECT_EQ(AT(m, 2, 2, 3), 614.0);
  tmc_array_function_free(f);
  tmc_sparse_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

static void one_and_a_half_first_into_second(double *const *x, void *ctx) /* out = 1.5 i */
{
  (void)ctx;
  *x[1] = *x[0] * 1.5;
}

/* (0,1) = 1.5 i, with the integer vector argument 1; e = { e01 } */
static void integer_vector_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { VEC(g, 1), ELEM(g, 0, ids->e[0]) };
  tmc_array_group_for_each(g, one_and_a_half_first_into_second, NULL, 2, ops);
}

int sparse_policy_integer_vector_sparse(tmc_sparse_ordering ordering, int L)
{
  /* MICM's VectorType<int>: the same integers held in doubles */
  tmc_sparse *m = build(ordering, L, 3, 1, one_element_pattern, 3, 0.0);
  const double values[] = { 5, 10, 15 };
  tmc_padded_vector *v = vector_of(L, 3, values);
  element_ids ids = { { elem(m, 0, 1), 0, 0, 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(integer_vector_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  TMC_EXPECT_DOUBLE_EQ(AT(m, 0, 0, 1), 7.5);
  TMC_EXPECT_DOUBLE_EQ(AT(m, 1, 0, 1), 15.0);
  TMC_EXPECT_DOUBLE_EQ(AT(m, 2, 0, 1), 22.5);
  tmc_array_function_free(f);
  tmc_sparse_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

int sparse_policy_function_with_const_signature_sparse(tmc_sparse_ordering ordering, int L)
{
  const int elements[1][2] = { { 0, 0 } };
  tmc_sparse *m = build(ordering, L, 2, 1, elements, 3, 0.0);
  const double values[] = { 1.0, 2.0, 3.0 };
  tmc_padded_vector *v = vector_of(L, 3, values);
  element_ids ids = { { elem(m, 0, 0), 0, 0, 0, 0, 0 } };
  tmc_array_arg args[] = { tmc_array_sparse(m), tmc_array_const_vector(v) };
  tmc_array_function *f = make_function(vector_doubled_into_element_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  /* MICM wraps the function in a std::function taking (SparseMatrix&, const
   * Vector&); here, a pointer to the apply function, called with the vector
   * read-only */
  int (*func_std)(const tmc_array_function *, int, const tmc_array_arg *, tmc_error *) = tmc_array_function_apply;
  tmc_error err;
  tmc_error_clear(&err);
  TMC_EXPECT_EQ(func_std(f, N_ARGS(args), args, &err), TMC_OK);
  /* (tighter than MICM, which checks nothing) the call ran */
  TMC_EXPECT_EQ(AT(m, 0, 0, 0), 2.0);
  TMC_EXPECT_EQ(AT(m, 1, 0, 0), 4.0);
  TMC_EXPECT_EQ(AT(m, 2, 0, 0), 6.0);
  tmc_array_function_free(f);
  tmc_sparse_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

static void vector_plus_value_into_second(double *const *x, void *ctx) /* elem = val + value */
{
  *x[1] = *x[0] + *(const double *)ctx;
}

/* for each of the six elements, in the order given: element = v + 10 k */
static void set_by_vector_index_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double value = 0.0;
  for (int k = 0; k < 6; ++k)
  {
    double *ops[] = { VEC(g, 1), ELEM(g, 0, ids->e[k]) };
    tmc_array_group_for_each(g, vector_plus_value_into_second, &value, 2, ops);
    value += 10.0;
  }
}

int sparse_policy_get_block_view_by_vector_index(tmc_sparse_ordering ordering, int L)
{
  /* X 0 X
   * X X 0
   * 0 X X */
  const int elements[6][2] = { { 0, 0 }, { 0, 2 }, { 1, 0 }, { 1, 1 }, { 2, 1 }, { 2, 2 } };
  tmc_sparse *sparse = build(ordering, L, 3, 6, elements, 3, 0.0);
  /* MICM's vector of VectorIndex(0, row, column) values: the element ids */
  element_ids ids;
  for (int k = 0; k < 6; ++k)
    ids.e[k] = elem(sparse, elements[k][0], elements[k][1]);
  const double block_values[] = { 0, 1, 2 };
  tmc_padded_vector *v = vector_of(L, 3, block_values);

  tmc_array_arg args[] = { tmc_array_sparse(sparse), tmc_array_const_vector(v) };
  tmc_array_function *set_func = make_function(set_by_vector_index_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!set_func)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(set_func, args);

  for (int block = 0; block < 3; ++block)
    for (int k = 0; k < 6; ++k)
      TMC_EXPECT_EQ(AT(sparse, block, elements[k][0], elements[k][1]), block + 10.0 * k);
  TMC_EXPECT_EQ(AT(sparse, 0, 0, 0), 0.0);
  TMC_EXPECT_EQ(AT(sparse, 1, 1, 0), 21.0);
  TMC_EXPECT_EQ(AT(sparse, 2, 2, 2), 52.0);

  tmc_array_function_free(set_func);
  tmc_sparse_free(sparse);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

/* ---- Fill and Copy -------------------------------------------------------- */

/* 0 X 0
 * 0 0 0
 * X 0 0 */
static const int fill_copy_elements[2][2] = { { 0, 1 }, { 2, 0 } };

/* e = { e01, e20 } for the fill and copy kernels */
static void fill_element_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  tmc_array_group_fill(g, ELEM(g, 0, ids->e[0]), 3.2);
}

static void fill_vector_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  tmc_array_group_fill(g, VEC(g, 1), 3.2);
}

/* tmp = 9.9, then (2,0) = tmp */
static void fill_variable_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *tmp = tmc_array_group_variable(g);
  tmc_array_group_fill(g, tmp, 9.9);
  double *ops[] = { ELEM(g, 0, ids->e[1]), tmp };
  tmc_array_group_for_each(g, copy_second_to_first, NULL, 2, ops);
}

int sparse_policy_fill(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 3, 2, fill_copy_elements, 3, 0.0);
  element_ids ids = { { elem(m, 0, 1), elem(m, 2, 0), 0, 0, 0, 0 } };

  /* a block element across all blocks */
  tmc_array_arg args[] = { tmc_array_sparse(m) };
  tmc_array_function *f = make_function(fill_element_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  for (int block = 0; block < m->n_blocks; ++block)
  {
    TMC_EXPECT_EQ(AT(m, block, 0, 1), 3.2);
    TMC_EXPECT_EQ(AT(m, block, 2, 0), 0.0); /* untouched */
  }
  tmc_array_function_free(f);

  /* a caller's vector of one value per block, through a read-only matrix */
  tmc_padded_vector *v = tmc_padded_vector_create(m->n_blocks, L, 0.0);
  tmc_array_arg vector_args[] = { tmc_array_const_sparse(m), tmc_array_vector(v) };
  f = make_function(fill_vector_kernel, NULL, N_ARGS(vector_args), vector_args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, vector_args);
  for (int i = 0; i < v->size; ++i)
    TMC_EXPECT_EQ(v->data[i], 3.2);
  /* (tighter than MICM) a sparse fill writes the whole group, padding included */
  for (int i = v->size; i < v->padded_size; ++i)
    TMC_EXPECT_EQ(v->data[i], 3.2);
  tmc_array_function_free(f);

  /* a block variable, broadcast into an unrelated element */
  f = make_function(fill_variable_kernel, &ids, N_ARGS(args), args, ordering, L);
  if (!f)
    return TMC_TEST_FAIL;
  EXPECT_APPLY_OK(f, args);
  for (int block = 0; block < m->n_blocks; ++block)
    TMC_EXPECT_EQ(AT(m, block, 2, 0), 9.9);
  tmc_array_function_free(f);

  tmc_sparse_free(m);
  tmc_padded_vector_free(v);
  return TMC_TEST_PASS;
}

/* A fresh matrix with distinct values in (0, 1) of each block. */
static tmc_sparse *copy_test_matrix(tmc_sparse_ordering ordering, int L)
{
  tmc_sparse *m = build(ordering, L, 3, 2, fill_copy_elements, 3, 0.0);
  SET(m, 0, 0, 1, 3.2);
  SET(m, 1, 0, 1, 4.2);
  SET(m, 2, 0, 1, 1.3);
  return m;
}

static const double copy_test_values[3] = { 3.2, 4.2, 1.3 };

/* e = { e01, e20 } for the copy kernels; (dst, src) as MICM's Copy */
static void copy_vector_into_element_kernel(tmc_array_group *g, void *ctx) /* (0,1) <- v */
{
  const element_ids *ids = ctx;
  tmc_array_group_copy(g, ELEM(g, 0, ids->e[0]), VEC(g, 1));
}

static void copy_element_into_vector_kernel(tmc_array_group *g, void *ctx) /* v <- (0,1) */
{
  const element_ids *ids = ctx;
  tmc_array_group_copy(g, VEC(g, 1), ELEM(g, 0, ids->e[0]));
}

static void copy_element_into_element_kernel(tmc_array_group *g, void *ctx) /* (2,0) <- (0,1) */
{
  const element_ids *ids = ctx;
  tmc_array_group_copy(g, ELEM(g, 0, ids->e[1]), ELEM(g, 0, ids->e[0]));
}

/* (0,1) -> block variable; (0,1) = 0; (0,1) = block variable */
static void copy_round_trip_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *tmp = tmc_array_group_variable(g);
  tmc_array_group_copy(g, tmp, ELEM(g, 0, ids->e[0]));
  tmc_array_group_fill(g, ELEM(g, 0, ids->e[0]), 0.0);
  double *ops[] = { ELEM(g, 0, ids->e[0]), tmp };
  tmc_array_group_for_each(g, copy_second_to_first, NULL, 2, ops);
}

/* (2,0) <- dense column 1, element by element, the dense matrix argument 1 */
static void copy_dense_column_into_element_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  double *ops[] = { ELEM(g, 0, ids->e[1]), COL(g, 1, 1) };
  tmc_array_group_for_each(g, copy_second_to_first, NULL, 2, ops);
}

/* a dense function: column 0 <- v, row by row, the vector argument 1 */
static void copy_vector_into_dense_column_kernel(tmc_array_group *g, void *ctx)
{
  (void)ctx;
  double *ops[] = { COL(g, 0, 0), VEC(g, 1) };
  tmc_array_group_for_each(g, copy_second_to_first, NULL, 2, ops);
}

/* dst(1,2) <- src(0,1), with dst argument 0 and src argument 1; e = { d12, s01 } */
static void copy_across_matrices_kernel(tmc_array_group *g, void *ctx)
{
  const element_ids *ids = ctx;
  tmc_array_group_copy(g, ELEM(g, 0, ids->e[0]), ELEM(g, 1, ids->e[1]));
}

int sparse_policy_copy(tmc_sparse_ordering ordering, int L)
{
  element_ids ids;
  {
    tmc_sparse *pattern = copy_test_matrix(ordering, L);
    ids.e[0] = elem(pattern, 0, 1);
    ids.e[1] = elem(pattern, 2, 0);
    tmc_sparse_free(pattern);
  }

  /* a vector into a block element */
  {
    tmc_sparse *m = copy_test_matrix(ordering, L);
    for (int block = 0; block < m->n_blocks; ++block)
      SET(m, block, 0, 1, 0.0);
    const double values[] = { 10.0, 20.0, 30.0 };
    tmc_padded_vector *v = vector_of(L, 3, values);
    tmc_array_arg args[] = { tmc_array_sparse(m), tmc_array_const_vector(v) };
    tmc_array_function *f = make_function(copy_vector_into_element_kernel, &ids, N_ARGS(args), args, ordering, L);
    if (!f)
      return TMC_TEST_FAIL;
    EXPECT_APPLY_OK(f, args);
    TMC_EXPECT_EQ(AT(m, 0, 0, 1), 10.0);
    TMC_EXPECT_EQ(AT(m, 1, 0, 1), 20.0);
    TMC_EXPECT_EQ(AT(m, 2, 0, 1), 30.0);
    TMC_EXPECT_EQ(AT(m, 0, 2, 0), 0.0); /* untouched */
    tmc_array_function_free(f);
    tmc_sparse_free(m);
    tmc_padded_vector_free(v);
  }

  /* a read-only block element into a vector */
  {
    tmc_sparse *m = copy_test_matrix(ordering, L);
    tmc_padded_vector *v = tmc_padded_vector_create(m->n_blocks, L, -1.0);
    tmc_array_arg args[] = { tmc_array_const_sparse(m), tmc_array_vector(v) };
    tmc_array_function *f = make_function(copy_element_into_vector_kernel, &ids, N_ARGS(args), args, ordering, L);
    if (!f)
      return TMC_TEST_FAIL;
    EXPECT_APPLY_OK(f, args);
    TMC_EXPECT_EQ(v->data[0], 3.2);
    TMC_EXPECT_EQ(v->data[1], 4.2);
    TMC_EXPECT_EQ(v->data[2], 1.3);
    /* (tighter than MICM) the whole group is copied, so the vector's padding
     * takes the matrix's padding, which holds the initial 0 */
    for (int i = v->size; i < v->padded_size; ++i)
      TMC_EXPECT_EQ(v->data[i], 0.0);
    tmc_array_function_free(f);
    tmc_sparse_free(m);
    tmc_padded_vector_free(v);
  }

  /* one block element into another; MICM reads the source through a mutable
   * view and then through a const one, which in C is the same accessor */
  for (int pass = 0; pass < 2; ++pass)
  {
    tmc_sparse *m = copy_test_matrix(ordering, L);
    tmc_array_arg args[] = { tmc_array_sparse(m) };
    tmc_array_function *f = make_function(copy_element_into_element_kernel, &ids, N_ARGS(args), args, ordering, L);
    if (!f)
      return TMC_TEST_FAIL;
    EXPECT_APPLY_OK(f, args);
    TMC_EXPECT_EQ(AT(m, 0, 2, 0), 3.2);
    TMC_EXPECT_EQ(AT(m, 1, 2, 0), 4.2);
    TMC_EXPECT_EQ(AT(m, 2, 2, 0), 1.3);
    TMC_EXPECT_EQ(AT(m, 0, 0, 1), 3.2); /* the source is unchanged */
    tmc_array_function_free(f);
    tmc_sparse_free(m);
  }

  /* block element -> block variable -> block element */
  {
    tmc_sparse *m = copy_test_matrix(ordering, L);
    tmc_array_arg args[] = { tmc_array_sparse(m) };
    tmc_array_function *f = make_function(copy_round_trip_kernel, &ids, N_ARGS(args), args, ordering, L);
    if (!f)
      return TMC_TEST_FAIL;
    EXPECT_APPLY_OK(f, args);
    TMC_EXPECT_EQ(AT(m, 0, 0, 1), 3.2);
    TMC_EXPECT_EQ(AT(m, 1, 0, 1), 4.2);
    TMC_EXPECT_EQ(AT(m, 2, 0, 1), 1.3);
    tmc_array_function_free(f);
    tmc_sparse_free(m);
  }

  /* a dense column into a block element, element by element */
  {
    tmc_sparse *m = copy_test_matrix(ordering, L);
    for (int block = 0; block < m->n_blocks; ++block)
      SET(m, block, 2, 0, 0.0);
    tmc_dense *dense = tmc_dense_create(m->n_blocks, 2, L, 0.0);
    for (int row = 0; row < dense->n_rows; ++row)
      tmc_dense_set(dense, row, 1, (row + 1) * 100.0);
    tmc_array_arg args[] = { tmc_array_sparse(m), tmc_array_const_dense(dense) };
    tmc_array_function *f = make_function(copy_dense_column_into_element_kernel, &ids, N_ARGS(args), args, ordering, L);
    if (!f)
      return TMC_TEST_FAIL;
    EXPECT_APPLY_OK(f, args);
    TMC_EXPECT_EQ(AT(m, 0, 2, 0), 100.0);
    TMC_EXPECT_EQ(AT(m, 1, 2, 0), 200.0);
    TMC_EXPECT_EQ(AT(m, 2, 2, 0), 300.0);
    tmc_array_function_free(f);
    tmc_sparse_free(m);
    tmc_dense_free(dense);
  }

  /* a block element into a dense column, through a caller's vector: a sparse
   * function copies the element into the vector, a dense function the vector
   * into the column */
  {
    tmc_sparse *m = copy_test_matrix(ordering, L);
    tmc_dense *dense = tmc_dense_create(m->n_blocks, 2, L, 0.0);
    tmc_padded_vector *scratch = tmc_padded_vector_create(m->n_blocks, L, 0.0);
    tmc_array_arg sparse_args[] = { tmc_array_const_sparse(m), tmc_array_vector(scratch) };
    tmc_array_function *sparse_to_scratch = make_function(copy_element_into_vector_kernel, &ids, N_ARGS(sparse_args), sparse_args, ordering, L);
    tmc_array_arg dense_args[] = { tmc_array_dense(dense), tmc_array_const_vector(scratch) };
    tmc_array_function *scratch_to_dense = make_function(copy_vector_into_dense_column_kernel, NULL, N_ARGS(dense_args), dense_args, ordering, L);
    if (!sparse_to_scratch || !scratch_to_dense)
      return TMC_TEST_FAIL;
    EXPECT_APPLY_OK(sparse_to_scratch, sparse_args);
    EXPECT_APPLY_OK(scratch_to_dense, dense_args);
    TMC_EXPECT_EQ(tmc_dense_get(dense, 0, 0), 3.2);
    TMC_EXPECT_EQ(tmc_dense_get(dense, 1, 0), 4.2);
    TMC_EXPECT_EQ(tmc_dense_get(dense, 2, 0), 1.3);
    for (int row = 0; row < 3; ++row)
      TMC_EXPECT_EQ(tmc_dense_get(dense, row, 1), 0.0);
    tmc_array_function_free(sparse_to_scratch);
    tmc_array_function_free(scratch_to_dense);
    tmc_sparse_free(m);
    tmc_dense_free(dense);
    tmc_padded_vector_free(scratch);
  }

  /* between two sparse matrices of the same block count and different patterns */
  {
    tmc_sparse *src = copy_test_matrix(ordering, L);
    const int dst_elements[1][2] = { { 1, 2 } };
    tmc_sparse *dst = build(ordering, L, 3, 1, dst_elements, 3, 0.0);
    element_ids cross = { { elem(dst, 1, 2), elem(src, 0, 1), 0, 0, 0, 0 } };
    tmc_array_arg args[] = { tmc_array_sparse(dst), tmc_array_const_sparse(src) };
    tmc_array_function *f = make_function(copy_across_matrices_kernel, &cross, N_ARGS(args), args, ordering, L);
    if (!f)
      return TMC_TEST_FAIL;
    EXPECT_APPLY_OK(f, args);
    TMC_EXPECT_EQ(AT(dst, 0, 1, 2), 3.2);
    TMC_EXPECT_EQ(AT(dst, 1, 1, 2), 4.2);
    TMC_EXPECT_EQ(AT(dst, 2, 1, 2), 1.3);
    for (int block = 0; block < 3; ++block)
      TMC_EXPECT_EQ(AT(src, block, 0, 1), copy_test_values[block]);
    tmc_array_function_free(f);
    tmc_sparse_free(src);
    tmc_sparse_free(dst);
  }
  return TMC_TEST_PASS;
}
