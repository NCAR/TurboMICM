/* MICM's dense container interface beyond what the solver uses: construction
 * from rows, checked row assignment, ForEach, printing, padded vectors, and
 * the Function facility over dense, sparse, and vector arguments. See
 * include/tmc_dense_ext.h. */
#include "internal.h"
#include "tmc_dense_ext.h"

#include <stdarg.h>
#include <stdio.h>

/* ------------------------------------------------------------------------- */
/* Dense matrices                                                            */
/* ------------------------------------------------------------------------- */

tmc_dense *tmc_dense_create_from_rows(int n_rows, const double *const *rows, const int *row_lengths, int L, tmc_error *err)
{
  if (n_rows < 0)
  {
    tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR, "In matrix constructor from rows. Got %d rows", n_rows);
    return NULL;
  }
  int n_cols = n_rows > 0 ? row_lengths[0] : 0;
  for (int i = 0; i < n_rows; ++i)
    if (row_lengths[i] != n_cols || row_lengths[i] < 0)
    {
      tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR, "In matrix constructor from rows. Got %d columns in row %d, but expected %d", row_lengths[i], i, n_cols);
      return NULL;
    }
  tmc_dense *m = tmc_dense_create(n_rows, n_cols, L, 0.0);
  for (int i = 0; i < n_rows; ++i)
    for (int j = 0; j < n_cols; ++j)
      m->data[tmc_dense_index(m, i, j)] = rows[i][j];
  return m;
}

int tmc_dense_assign_row(tmc_dense *m, int row, const double *values, int n_values, tmc_error *err)
{
  if (row < 0 || row >= m->n_rows)
    return tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_ELEMENT_OUT_OF_RANGE, "Row %d is out of range for a matrix of %d rows", row, m->n_rows);
  if (n_values < m->n_cols)
    return tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_ROW_SIZE_MISMATCH, "In matrix row assignment. Got %d elements, but expected %d", n_values, m->n_cols);
  for (int j = 0; j < m->n_cols; ++j)
    m->data[tmc_dense_index(m, row, j)] = values[j];
  return TMC_OK;
}

int tmc_dense_number_of_groups(const tmc_dense *m)
{
  return m->n_groups;
}

int tmc_dense_group_size(const tmc_dense *m)
{
  return m->L * m->n_cols;
}

int tmc_dense_group_width(const tmc_dense *m)
{
  return m->L;
}

/* Real rows of group g (L, or fewer in the last group). */
static int dense_group_cells(const tmc_dense *m, int g)
{
  int cells = m->n_rows - g * m->L;
  return cells < m->L ? cells : m->L;
}

/* Elements in storage order, padding lanes skipped (MICM's VectorMatrix::ForEach). */
void tmc_dense_for_each(tmc_dense *m, tmc_dense_unary_fn f, void *ctx, const tmc_dense *a)
{
  for (int g = 0; g < m->n_groups; ++g)
  {
    int cells = dense_group_cells(m, g);
    for (int j = 0; j < m->n_cols; ++j)
    {
      size_t base = ((size_t)g * (size_t)m->n_cols + (size_t)j) * (size_t)m->L;
      for (int i = 0; i < cells; ++i)
        f(&m->data[base + i], &a->data[base + i], ctx);
    }
  }
}

void tmc_dense_for_each2(tmc_dense *m, tmc_dense_binary_fn f, void *ctx, const tmc_dense *a, const tmc_dense *b)
{
  for (int g = 0; g < m->n_groups; ++g)
  {
    int cells = dense_group_cells(m, g);
    for (int j = 0; j < m->n_cols; ++j)
    {
      size_t base = ((size_t)g * (size_t)m->n_cols + (size_t)j) * (size_t)m->L;
      for (int i = 0; i < cells; ++i)
        f(&m->data[base + i], &a->data[base + i], &b->data[base + i], ctx);
    }
  }
}

/* The text of element (i, j) and the character that follows it. */
static int format_element(const tmc_dense *m, int i, int j, char *text, size_t size)
{
  char separator = j + 1 < m->n_cols ? ',' : '\n';
  return snprintf(text, size, "%g%c", tmc_dense_get(m, i, j), separator);
}

size_t tmc_dense_format(const tmc_dense *m, char *buffer, size_t size)
{
  size_t length = 0;
  char text[64];
  for (int i = 0; i < m->n_rows; ++i)
  {
    if (m->n_cols == 0)
    {
      if (length + 1 < size)
        buffer[length] = '\n';
      ++length;
      continue;
    }
    for (int j = 0; j < m->n_cols; ++j)
    {
      int n = format_element(m, i, j, text, sizeof text);
      for (int k = 0; k < n; ++k, ++length)
        if (length + 1 < size)
          buffer[length] = text[k];
    }
  }
  if (size > 0)
    buffer[length < size ? length : size - 1] = '\0';
  return length;
}

int tmc_dense_print(const tmc_dense *m, FILE *stream)
{
  int total = 0;
  char text[64];
  for (int i = 0; i < m->n_rows; ++i)
  {
    if (m->n_cols == 0)
    {
      if (fputc('\n', stream) == EOF)
        return -1;
      ++total;
      continue;
    }
    for (int j = 0; j < m->n_cols; ++j)
    {
      int n = format_element(m, i, j, text, sizeof text);
      if (fputs(text, stream) == EOF)
        return -1;
      total += n;
    }
  }
  return total;
}

/* ------------------------------------------------------------------------- */
/* Padded vectors                                                            */
/* ------------------------------------------------------------------------- */

static tmc_padded_vector *padded_vector_alloc(int size, int L)
{
  if (L < 1)
    L = 1;
  if (size < 0)
    size = 0;
  tmc_padded_vector *v = tmc_calloc(1, sizeof *v);
  v->size = size;
  v->L = L;
  v->padded_size = (size + L - 1) / L * L;
  v->data = tmc_alloc_doubles((size_t)v->padded_size);
  return v;
}

tmc_padded_vector *tmc_padded_vector_create(int size, int L, double initial_value)
{
  tmc_padded_vector *v = padded_vector_alloc(size, L);
  for (int i = 0; i < v->padded_size; ++i)
    v->data[i] = initial_value;
  return v;
}

tmc_padded_vector *tmc_padded_vector_from_values(const double *values, int size, int L)
{
  tmc_padded_vector *v = padded_vector_alloc(size, L);
  for (int i = 0; i < v->padded_size; ++i)
    v->data[i] = i < v->size ? values[i] : 0.0;
  return v;
}

void tmc_padded_vector_free(tmc_padded_vector *v)
{
  if (!v)
    return;
  tmc_free_doubles(v->data);
  free(v);
}

/* ------------------------------------------------------------------------- */
/* Functions over arrays                                                     */
/* ------------------------------------------------------------------------- */

tmc_array_arg tmc_array_dense(tmc_dense *m)
{
  tmc_array_arg a = { TMC_ARRAY_DENSE, false, m, NULL, NULL };
  return a;
}

tmc_array_arg tmc_array_const_dense(const tmc_dense *m)
{
  tmc_array_arg a = { TMC_ARRAY_DENSE, true, (tmc_dense *)m, NULL, NULL };
  return a;
}

tmc_array_arg tmc_array_sparse(tmc_sparse *m)
{
  tmc_array_arg a = { TMC_ARRAY_SPARSE, false, NULL, m, NULL };
  return a;
}

tmc_array_arg tmc_array_const_sparse(const tmc_sparse *m)
{
  tmc_array_arg a = { TMC_ARRAY_SPARSE, true, NULL, (tmc_sparse *)m, NULL };
  return a;
}

tmc_array_arg tmc_array_vector(tmc_padded_vector *v)
{
  tmc_array_arg a = { TMC_ARRAY_VECTOR, false, NULL, NULL, v };
  return a;
}

tmc_array_arg tmc_array_const_vector(const tmc_padded_vector *v)
{
  tmc_array_arg a = { TMC_ARRAY_VECTOR, true, NULL, NULL, (tmc_padded_vector *)v };
  return a;
}

struct tmc_array_function
{
  tmc_array_kernel kernel;
  void *ctx;
  int n_args;
  int L;
  tmc_array_kind *kinds; /* [n_args] */
  int *columns;          /* [n_args] column counts of the dense arguments, 0 for the others */
};

/* Row variables of one application: slots of L doubles, handed out afresh to
 * each group (the same slots, so a group sees stale values, as MICM's
 * uninitialized stack variables would hold). */
typedef struct variable_pool
{
  double **slots;
  int n, cap;
  int L;
} variable_pool;

struct tmc_array_group
{
  const tmc_array_function *f;
  const tmc_array_arg *args;
  int index; /* group number */
  int cells; /* real lanes */
  variable_pool *pool;
  int next_variable;
  tmc_error error; /* the first misuse inside the kernel */
};

static const char *kind_name(tmc_array_kind kind)
{
  switch (kind)
  {
    case TMC_ARRAY_DENSE: return "dense matrix";
    case TMC_ARRAY_SPARSE: return "sparse matrix";
    case TMC_ARRAY_VECTOR: return "vector";
  }
  return "unknown array";
}

static bool arg_present(const tmc_array_arg *a)
{
  switch (a->kind)
  {
    case TMC_ARRAY_DENSE: return a->dense != NULL;
    case TMC_ARRAY_SPARSE: return a->sparse != NULL;
    case TMC_ARRAY_VECTOR: return a->vector != NULL;
  }
  return false;
}

/* The argument's group width L (the argument must be present). */
static int arg_width(const tmc_array_arg *a)
{
  switch (a->kind)
  {
    case TMC_ARRAY_DENSE: return a->dense->L;
    case TMC_ARRAY_SPARSE: return a->sparse->L;
    case TMC_ARRAY_VECTOR: return a->vector->L;
  }
  return 0;
}

/* The argument's lane count: rows, blocks, or values. */
static int arg_lanes(const tmc_array_arg *a)
{
  switch (a->kind)
  {
    case TMC_ARRAY_DENSE: return a->dense->n_rows;
    case TMC_ARRAY_SPARSE: return a->sparse->n_blocks;
    case TMC_ARRAY_VECTOR: return a->vector->size;
  }
  return 0;
}

tmc_array_function *tmc_array_function_create(tmc_array_kernel kernel, void *ctx, int n_args, const tmc_array_arg *args, tmc_error *err)
{
  if (!kernel)
  {
    tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR, "A function needs a kernel");
    return NULL;
  }
  if (n_args < 0 || (n_args > 0 && !args))
  {
    tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR, "A function needs its %d arguments", n_args);
    return NULL;
  }
  int L = 1;
  for (int k = 0; k < n_args; ++k)
  {
    if (!arg_present(&args[k]))
    {
      tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR, "Argument %d (a %s) is missing", k, kind_name(args[k].kind));
      return NULL;
    }
    int width = arg_width(&args[k]);
    if (k == 0)
      L = width;
    else if (width != L)
    {
      tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR,
                    "Incompatible arguments: argument %d (a %s) has L = %d but argument 0 has L = %d", k, kind_name(args[k].kind), width, L);
      return NULL;
    }
  }
  tmc_array_function *f = tmc_calloc(1, sizeof *f);
  f->kernel = kernel;
  f->ctx = ctx;
  f->n_args = n_args;
  f->L = L;
  f->kinds = tmc_malloc((size_t)n_args * sizeof *f->kinds);
  f->columns = tmc_malloc((size_t)n_args * sizeof *f->columns);
  for (int k = 0; k < n_args; ++k)
  {
    f->kinds[k] = args[k].kind;
    f->columns[k] = args[k].kind == TMC_ARRAY_DENSE ? args[k].dense->n_cols : 0;
  }
  return f;
}

void tmc_array_function_free(tmc_array_function *f)
{
  if (!f)
    return;
  free(f->kinds);
  free(f->columns);
  free(f);
}

/* The invocation checks of MICM's Function, in argument order; sets *n_lanes. */
static int validate(const tmc_array_function *f, int n_args, const tmc_array_arg *args, int *n_lanes, tmc_error *err)
{
  if (n_args != f->n_args || (n_args > 0 && !args))
    return tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR, "The function takes %d arguments but was given %d", f->n_args, n_args);
  int lanes = 0;
  for (int k = 0; k < n_args; ++k)
  {
    const tmc_array_arg *a = &args[k];
    if (a->kind != f->kinds[k])
      return tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR, "Argument %d is a %s but the function was created with a %s", k, kind_name(a->kind),
                           kind_name(f->kinds[k]));
    if (!arg_present(a))
      return tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR, "Argument %d (a %s) is missing", k, kind_name(a->kind));
    if (arg_width(a) != f->L)
      return tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR, "Argument %d has L = %d but the function was created with L = %d", k, arg_width(a), f->L);
    if (a->kind == TMC_ARRAY_SPARSE && !a->sparse->data)
      return tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR, "Sparse argument %d holds no values (it is indexing only)", k);
    int count = arg_lanes(a);
    if (k == 0)
      lanes = count;
    else if (count != lanes)
    {
      if (a->kind == TMC_ARRAY_VECTOR)
        return tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR, "Vector size must match matrix row count. Expected %d elements but got %d", lanes, count);
      return tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR,
                           "All matrices must have the same number of rows (blocks) when invoking function. Expected %d but got %d", lanes, count);
    }
    if (a->kind == TMC_ARRAY_DENSE && a->dense->n_cols != f->columns[k])
      return tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR, "Matrix column count does not match. Expected %d columns but got %d", f->columns[k],
                           a->dense->n_cols);
  }
  *n_lanes = lanes;
  return TMC_OK;
}

static int run_group(const tmc_array_function *f, const tmc_array_arg *args, int n_lanes, int group, variable_pool *pool, tmc_error *err)
{
  tmc_array_group g;
  g.f = f;
  g.args = args;
  g.index = group;
  int cells = n_lanes - group * f->L;
  g.cells = cells < f->L ? cells : f->L;
  g.pool = pool;
  g.next_variable = 0;
  tmc_error_clear(&g.error);
  f->kernel(&g, f->ctx);
  if (g.error.category != TMC_OK)
  {
    if (err)
      *err = g.error;
    return g.error.category;
  }
  return TMC_OK;
}

static void pool_free(variable_pool *pool)
{
  for (int i = 0; i < pool->n; ++i)
    tmc_free_doubles(pool->slots[i]);
  free(pool->slots);
}

/* Runs groups [first, last). */
static int run_groups(const tmc_array_function *f, const tmc_array_arg *args, int n_lanes, int first, int last, tmc_error *err)
{
  variable_pool pool = { NULL, 0, 0, f->L };
  int rc = TMC_OK;
  for (int group = first; group < last && rc == TMC_OK; ++group)
    rc = run_group(f, args, n_lanes, group, &pool, err);
  pool_free(&pool);
  return rc;
}

int tmc_array_function_apply(const tmc_array_function *f, int n_args, const tmc_array_arg *args, tmc_error *err)
{
  int n_lanes = 0;
  int rc = validate(f, n_args, args, &n_lanes, err);
  if (rc != TMC_OK)
    return rc;
  int n_groups = (n_lanes + f->L - 1) / f->L;
  return run_groups(f, args, n_lanes, 0, n_groups, err);
}

int tmc_array_function_apply_group(const tmc_array_function *f, int n_args, const tmc_array_arg *args, int group, tmc_error *err)
{
  int n_lanes = 0;
  int rc = validate(f, n_args, args, &n_lanes, err);
  if (rc != TMC_OK)
    return rc;
  int n_groups = (n_lanes + f->L - 1) / f->L;
  if (group < 0 || group >= n_groups)
    return tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_ELEMENT_OUT_OF_RANGE, "Group %d is out of range for %d groups", group, n_groups);
  return run_groups(f, args, n_lanes, group, group + 1, err);
}

/* ---- inside a kernel ---------------------------------------------------- */

/* Records the first misuse in a kernel; later ones are dropped. */
static void group_fail(tmc_array_group *g, int code, const char *format, ...)
{
  if (g->error.category != TMC_OK)
    return;
  g->error.category = TMC_ERROR_MATRIX;
  g->error.code = code;
  va_list args;
  va_start(args, format);
  vsnprintf(g->error.message, sizeof g->error.message, format, args);
  va_end(args);
}

int tmc_array_group_index(const tmc_array_group *g)
{
  return g->index;
}

int tmc_array_group_width(const tmc_array_group *g)
{
  return g->f->L;
}

int tmc_array_group_cells(const tmc_array_group *g)
{
  return g->cells;
}

static const tmc_array_arg *group_arg(tmc_array_group *g, int arg, tmc_array_kind kind)
{
  if (g->error.category != TMC_OK)
    return NULL;
  if (arg < 0 || arg >= g->f->n_args)
  {
    group_fail(g, TMC_MATRIX_ELEMENT_OUT_OF_RANGE, "Argument %d is out of range for a function of %d arguments", arg, g->f->n_args);
    return NULL;
  }
  const tmc_array_arg *a = &g->args[arg];
  if (a->kind != kind)
  {
    group_fail(g, TMC_MATRIX_INVALID_VECTOR, "Argument %d is a %s, not a %s", arg, kind_name(a->kind), kind_name(kind));
    return NULL;
  }
  return a;
}

double *tmc_array_group_column(tmc_array_group *g, int arg, int column)
{
  const tmc_array_arg *a = group_arg(g, arg, TMC_ARRAY_DENSE);
  if (!a)
    return NULL;
  const tmc_dense *m = a->dense;
  if (column < 0 || column >= m->n_cols)
  {
    group_fail(g, TMC_MATRIX_ELEMENT_OUT_OF_RANGE, "Column %d is out of range for a matrix of %d columns", column, m->n_cols);
    return NULL;
  }
  return m->data + ((size_t)g->index * (size_t)m->n_cols + (size_t)column) * (size_t)m->L;
}

double *tmc_array_group_element(tmc_array_group *g, int arg, int element)
{
  const tmc_array_arg *a = group_arg(g, arg, TMC_ARRAY_SPARSE);
  if (!a)
    return NULL;
  const tmc_sparse *m = a->sparse;
  if (element < 0)
  {
    group_fail(g, TMC_MATRIX_ZERO_ELEMENT_ACCESS, "Element %d is a structural zero", element);
    return NULL;
  }
  if (element >= m->nnz)
  {
    group_fail(g, TMC_MATRIX_ELEMENT_OUT_OF_RANGE, "Element %d is out of range for a matrix of %d non-zero elements", element, m->nnz);
    return NULL;
  }
  return m->data + ((size_t)g->index * (size_t)m->nnz + (size_t)element) * (size_t)m->L;
}

double *tmc_array_group_vector(tmc_array_group *g, int arg)
{
  const tmc_array_arg *a = group_arg(g, arg, TMC_ARRAY_VECTOR);
  if (!a)
    return NULL;
  return a->vector->data + (size_t)g->index * (size_t)a->vector->L;
}

double *tmc_array_group_variable(tmc_array_group *g)
{
  if (g->error.category != TMC_OK)
    return NULL;
  variable_pool *pool = g->pool;
  if (g->next_variable == pool->n)
  {
    if (pool->n == pool->cap)
    {
      pool->cap = pool->cap ? 2 * pool->cap : 4;
      pool->slots = tmc_realloc(pool->slots, (size_t)pool->cap * sizeof *pool->slots);
    }
    double *slot = tmc_alloc_doubles((size_t)pool->L);
    for (int i = 0; i < pool->L; ++i)
      slot[i] = 0.0;
    pool->slots[pool->n++] = slot;
  }
  return pool->slots[g->next_variable++];
}

/* Checks the operands of for_each and reduce; false once an error is recorded. */
static bool operands_ok(tmc_array_group *g, int n, double *const *lanes)
{
  if (g->error.category != TMC_OK)
    return false;
  if (n < 0 || (n > 0 && !lanes))
  {
    group_fail(g, TMC_MATRIX_INVALID_VECTOR, "Expected %d operands", n);
    return false;
  }
  for (int k = 0; k < n; ++k)
    if (!lanes[k])
    {
      group_fail(g, TMC_MATRIX_INVALID_VECTOR, "Operand %d has no lanes", k);
      return false;
    }
  return true;
}

enum
{
  MAX_STACK_OPERANDS = 16
};

static void for_each_lane(tmc_array_group *g, int count, tmc_array_lane_fn f, void *ctx, int n, double *const *lanes)
{
  if (!operands_ok(g, n, lanes))
    return;
  double *stack[MAX_STACK_OPERANDS];
  double **x = n <= MAX_STACK_OPERANDS ? stack : tmc_malloc((size_t)n * sizeof *x);
  for (int i = 0; i < count; ++i)
  {
    for (int k = 0; k < n; ++k)
      x[k] = lanes[k] + i;
    f(x, ctx);
  }
  if (x != stack)
    free(x);
}

static void reduce_lanes(tmc_array_group *g, int count, tmc_array_reduce_fn f, void *ctx, int n, double *const *lanes, void *accumulator)
{
  if (!operands_ok(g, n, lanes))
    return;
  double *stack[MAX_STACK_OPERANDS];
  double **x = n <= MAX_STACK_OPERANDS ? stack : tmc_malloc((size_t)n * sizeof *x);
  for (int i = 0; i < count; ++i)
  {
    for (int k = 0; k < n; ++k)
      x[k] = lanes[k] + i;
    f(x, accumulator, ctx);
  }
  if (x != stack)
    free(x);
}

void tmc_array_group_for_each(tmc_array_group *g, tmc_array_lane_fn f, void *ctx, int n, double *const *lanes)
{
  for_each_lane(g, g->f->L, f, ctx, n, lanes);
}

void tmc_array_group_for_each_strict(tmc_array_group *g, tmc_array_lane_fn f, void *ctx, int n, double *const *lanes)
{
  for_each_lane(g, g->cells, f, ctx, n, lanes);
}

void tmc_array_group_reduce(tmc_array_group *g, tmc_array_reduce_fn f, void *ctx, int n, double *const *lanes, void *accumulator)
{
  reduce_lanes(g, g->f->L, f, ctx, n, lanes, accumulator);
}

void tmc_array_group_reduce_strict(tmc_array_group *g, tmc_array_reduce_fn f, void *ctx, int n, double *const *lanes, void *accumulator)
{
  reduce_lanes(g, g->cells, f, ctx, n, lanes, accumulator);
}

void tmc_array_group_fill(tmc_array_group *g, double *lanes, double value)
{
  if (!operands_ok(g, 1, &lanes))
    return;
  for (int i = 0; i < g->f->L; ++i)
    lanes[i] = value;
}

static void copy_lanes(tmc_array_group *g, int count, double *dst, const double *src)
{
  if (g->error.category != TMC_OK)
    return;
  if (!dst || !src)
  {
    group_fail(g, TMC_MATRIX_INVALID_VECTOR, "Copy has no %s lanes", dst ? "source" : "destination");
    return;
  }
  for (int i = 0; i < count; ++i)
    dst[i] = src[i];
}

void tmc_array_group_copy(tmc_array_group *g, double *dst, const double *src)
{
  copy_lanes(g, g->f->L, dst, src);
}

void tmc_array_group_copy_strict(tmc_array_group *g, double *dst, const double *src)
{
  copy_lanes(g, g->cells, dst, src);
}
