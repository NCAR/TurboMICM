#include "test.h"
#include "chem_matrix.h"
#define S(m, b, r, c) (*(double *)ck_sparse_at(&(m), (b), (r), (c)))
#define N(m, r, c) ((m).map[(r) * (m).order + (c)])

static int sparse_new(ck_sparse *m, size_t n, size_t blocks, size_t width, ck_sparse_order ordering,
                      const size_t (*entries)[2], size_t count, double initial)
{
  unsigned char *pattern = calloc(n * n, 1);
  CHECK(pattern);
  for (size_t j = 0; j < count; ++j)
  {
    CHECK(entries[j][0] < n && entries[j][1] < n);
    pattern[entries[j][0] * n + entries[j][1]] = 1;
  }
  ck_status status = ck_sparse_create(n, blocks, width, ordering, pattern, sizeof(double), &initial, m);
  free(pattern);
  return status;
}
static int text_equal(FILE *file, const char *expected)
{
  CHECK(fflush(file) == 0);
  rewind(file);
  size_t len = strlen(expected);
  char *text = calloc(len + 2, 1);
  CHECK(text);
  size_t got = fread(text, 1, len + 1, file);
  int equal = got == len && !strcmp(text, expected);
  if (!equal)
    fprintf(stderr, "sparse print mismatch: %s\n", text);
  free(text);
  CHECK(equal);
  return 0;
}
static int sparse_basic(const char *test, size_t width, ck_sparse_order ordering)
{
  const size_t entries[][2] = {{0, 1}, {3, 2}, {2, 3}, {1, 1}, {2, 1}};
  if (IS("ZeroMatrix") || IS("ConstZeroMatrix") || IS("SetScalar"))
  {
    ck_sparse m;
    CHECK(sparse_new(&m, 3, 1, width, ordering, NULL, 0, 0) == CK_OK);
    CHECK(m.nnz == 0);
    CHECK(ck_sparse_index(&m, 0, 0, 0) == SIZE_MAX);
    CHECK(ck_sparse_at(&m, 0, 1, 1) == NULL);
    double two = 2;
    ck_matrix_fill(&m.values, &two);
    CHECK(m.values.size == 0);
    ck_sparse_free(&m);
    return 0;
  }
  if (IS("SingleBlockMatrix") || IS("ConstSingleBlockMatrix") || IS("MultiBlockMatrix"))
  {
    int multi = IS("MultiBlockMatrix");
    ck_sparse m;
    CHECK(sparse_new(&m, 4, multi ? 53 : 1, width, ordering, entries, 5, multi ? 24 : 0) == CK_OK);
    CHECK(m.nnz == 5);
    CHECK(N(m, 0, 1) != SIZE_MAX && N(m, 3, 2) != SIZE_MAX && N(m, 2, 1) != SIZE_MAX);
    CHECK(N(m, 2, 2) == SIZE_MAX && N(m, 0, 0) == SIZE_MAX && N(m, 3, 3) == SIZE_MAX);
    NEAR(S(m, 0, 2, 1), multi ? 24 : 0, 0);
    S(m, 0, 2, 1) = 45;
    NEAR(S(m, 0, 2, 1), 45, 0);
    if (multi)
    {
      S(m, 43, 3, 2) = 29;
      S(m, 33, 2, 3) = 79;
      S(m, 2, 2, 3) = 64;
      NEAR(S(m, 43, 3, 2), 29, 0);
      NEAR(S(m, 33, 2, 3), 79, 0);
      NEAR(S(m, 2, 2, 3), 64, 0);
    }
    else if (IS("ConstSingleBlockMatrix"))
    {
      S(m, 0, 3, 2) = 42;
      S(m, 0, 2, 3) = 21;
      const ck_sparse *cm = &m;
      NEAR(*(const double *)ck_sparse_const_at(cm, 0, 2, 1), 45, 0);
      NEAR(*(const double *)ck_sparse_const_at(cm, 0, 3, 2), 42, 0);
      NEAR(*(const double *)ck_sparse_const_at(cm, 0, 2, 3), 21, 0);
    }
    size_t diagonal_count = 0;
    for (size_t i = 0; i < 4; ++i)
      if (N(m, i, i) != SIZE_MAX)
        ++diagonal_count;
    CHECK(diagonal_count == 1 && N(m, 1, 1) != SIZE_MAX);
    CHECK(ck_sparse_index(&m, m.values.rows, 0, 1) == SIZE_MAX);
    CHECK(ck_sparse_at(&m, 0, 0, 4) == NULL);
    CHECK(ck_sparse_index(&m, 0, 2, 2) == SIZE_MAX);
    ck_sparse_free(&m);
    return 0;
  }
  if (IS("AddToDiagonal") || IS("Print") || IS("PrintNonZero"))
  {
    int nonzero = IS("PrintNonZero");
    ck_sparse m;
    CHECK(sparse_new(&m, 4, 3, width, ordering, entries, 5, nonzero ? 32 : 24) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
      S(m, b, 1, 1) = b * 2 + (nonzero ? 2 : 1);
    CHECK(ck_sparse_diagonal(&m, nonzero ? 5 : 7) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
    {
      NEAR(S(m, b, 1, 1), b * 2 + (nonzero ? 7 : 8), 0);
      for (size_t j = 0; j < 5; ++j)
        if (entries[j][0] != 1 || entries[j][1] != 1)
          NEAR(S(m, b, entries[j][0], entries[j][1]), nonzero ? 32 : 24, 0);
    }
    if (!IS("AddToDiagonal"))
    {
      FILE *file = tmpfile();
      CHECK(file);
      CHECK(ck_sparse_print(&m, file, nonzero) == CK_OK);
      const char *expected =
          nonzero
              ? "Block 0\n0, 1, 32\n1, 1, 7\n2, 1, 32\n2, 3, 32\n3, 2, 32\nBlock 1\n0, 1, 32\n1, 1, 9\n2, 1, "
                "32\n2, 3, 32\n3, 2, 32\nBlock 2\n0, 1, 32\n1, 1, 11\n2, 1, 32\n2, 3, 32\n3, 2, 32\n"
              : "Block 0\n0,24,0,0\n0,8,0,0\n0,24,0,24\n0,0,24,0\nBlock "
                "1\n0,24,0,0\n0,10,0,0\n0,24,0,24\n0,0,24,0\nBlock "
                "2\n0,24,0,0\n0,12,0,0\n0,24,0,24\n0,0,24,0\n";
      CHECK(text_equal(file, expected) == 0);
      fclose(file);
    }
    ck_sparse_free(&m);
    return 0;
  }
  return 2;
}
/* A small instruction representation is the C equivalent of reusable captured
   test lambdas. It drives the library's checked array functions and block views. */
typedef struct
{
  size_t argument, entry;
} slot;
typedef struct
{
  char operation;
  slot target, sources[4];
  size_t count;
  double factor;
} instruction;
typedef struct
{
  instruction steps[8];
  size_t count, temporaries;
} program;
static slot at(size_t arg, size_t entry)
{
  return (slot){arg, entry};
}
static void step(program *p, char op, slot target, size_t n, slot a, slot b, slot c, slot d, double factor)
{
  p->steps[p->count++] =
      (instruction){.operation = op, .target = target, .sources = {a, b, c, d}, .count = n, .factor = factor};
}
static ck_status slot_view(const ck_argument *args, slot s, int writable, ck_view *view)
{
  const ck_argument *a = &args[s.argument];
  if (writable && !a->writable)
    return CK_INVALID;
  if (a->kind == CK_ARGUMENT_VECTOR)
  {
    *view = a->vector;
    if (!writable)
      view->writable = 0;
    return CK_OK;
  }
  if (a->kind == CK_ARGUMENT_MATRIX)
    return ck_matrix_column(a->matrix, s.entry, writable, view);
  return ck_sparse_block_view(a->sparse, s.entry, writable, view);
}
static ck_status instruction_row(ck_view *v, size_t count, size_t row, void *context)
{
  instruction *ins = context;
  if (v[0].element_size == sizeof(int))
  {
    int *out = ck_view_at(&v[0], row);
    const int *in = ck_view_const_at(&v[1], row);
    if (!out || !in)
      return CK_INVALID;
    *out = (int)(*in * ins->factor);
    return CK_OK;
  }
  double value = ins->operation == '*' ? 1 : 0;
  for (size_t j = 1; j < count; ++j)
  {
    const double *in = ck_view_const_at(&v[j], row);
    if (!in)
      return CK_INVALID;
    if (ins->operation == '*')
      value *= *in;
    else if (ins->operation == '-' && j > 1)
      value -= *in;
    else
      value += v[j].element_size == sizeof(int) ? *(const int *)in : *in;
  }
  double *out = ck_view_at(&v[0], row);
  if (!out)
    return CK_INVALID;
  *out = ins->factor * value;
  return CK_OK;
}
static ck_status execute(const ck_argument *arguments, size_t count, void *context)
{
  program *p = context;
  ck_argument args[8];
  if (count + p->temporaries > 8)
    return CK_INVALID;
  memcpy(args, arguments, count * sizeof *args);
  ck_matrix temps[2] = {{0}};
  size_t rows =
      arguments[0].kind == CK_ARGUMENT_SPARSE ? arguments[0].sparse->values.rows : arguments[0].matrix->rows;
  ck_status status = CK_OK;
  for (size_t j = 0; j < p->temporaries; ++j)
  {
    status = ck_matrix_create(rows, 1, 0, sizeof(double), NULL, &temps[j]);
    if (status)
      goto done;
    args[count + j] = ck_matrix_argument(&temps[j], 1);
  }
  for (size_t j = 0; j < p->count; ++j)
  {
    instruction *ins = &p->steps[j];
    ck_view views[5];
    status = slot_view(args, ins->target, 1, &views[0]);
    if (status)
      goto done;
    for (size_t k = 0; k < ins->count; ++k)
    {
      status = slot_view(args, ins->sources[k], 0, &views[k + 1]);
      if (status)
        goto done;
    }
    status = ck_view_apply(views, ins->count + 1, 1, instruction_row, ins);
    if (status)
      goto done;
  }
done:
  for (size_t j = 0; j < p->temporaries; ++j)
    ck_matrix_free(&temps[j]);
  return status;
}
static int execute_once(program *p, ck_argument *args, size_t count)
{
  ck_array_function f;
  CHECK(ck_array_function_create(args, count, execute, p, &f) == CK_OK);
  ck_status status = ck_array_function_invoke(&f, args, count);
  ck_array_function_free(&f);
  return status;
}
static int sparse_functions(const char *test, size_t width, ck_sparse_order order)
{
  slot unused = {0};
  if (IS("ArrayFunction"))
  {
    const size_t entries[][2] = {{0, 0}, {1, 1}, {2, 2}, {2, 3}, {3, 3}};
    ck_sparse m, other;
    CHECK(sparse_new(&m, 4, 3, width, order, entries, 5, 1) == CK_OK);
    CHECK(sparse_new(&other, 4, 3, width, order, entries, 5, -1) == CK_OK);
    S(m, 0, 0, 0) = 1;
    S(m, 1, 1, 1) = 2;
    S(m, 2, 2, 2) = 3;
    S(m, 2, 2, 3) = 4;
    S(m, 2, 3, 3) = 5;
    program p = {0};
    step(&p, '+', at(0, N(m, 2, 3)), 4, at(0, N(m, 0, 0)), at(0, N(m, 1, 1)), at(0, N(m, 2, 2)),
         at(0, N(m, 2, 3)), 2);
    ck_argument a[] = {ck_sparse_argument(&m, 1)};
    ck_array_function f;
    CHECK(ck_array_function_create(a, 1, execute, &p, &f) == CK_OK);
    CHECK(ck_array_function_invoke(&f, a, 1) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
    {
      NEAR(S(m, b, 2, 3), b == 0 ? 8 : b == 1 ? 10 : 18, 0);
      NEAR(S(m, b, 0, 0), 1, 0);
      NEAR(S(m, b, 1, 1), b == 1 ? 2 : 1, 0);
      NEAR(S(m, b, 2, 2), b == 2 ? 3 : 1, 0);
      NEAR(S(m, b, 3, 3), b == 2 ? 5 : 1, 0);
    }
    a[0] = ck_sparse_argument(&other, 1);
    CHECK(ck_array_function_invoke(&f, a, 1) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
      for (size_t j = 0; j < 5; ++j)
        NEAR(S(other, b, entries[j][0], entries[j][1]), j == 3 ? -8 : -1, 0);
    ck_array_function_free(&f);
    ck_sparse_free(&m);
    ck_sparse_free(&other);
    return 0;
  }
  if (IS("MultiMatrixArrayFunction") || IS("TwoSparseMatricesDifferentStructure") ||
      IS("MultipleSparseMatricesDifferentBlocksFromCreation") ||
      IS("SparseMatrixVectorDifferentBlocksFromCreation"))
  {
    int different = IS("TwoSparseMatricesDifferentStructure"),
        blocks = IS("MultipleSparseMatricesDifferentBlocksFromCreation"),
        vector = IS("SparseMatrixVectorDifferentBlocksFromCreation");
    const size_t e1[][2] = {{0, 1}, {1, 2}}, e2[][2] = {{0, 0}, {1, 1}, {2, 2}},
                 e3[][2] = {{0, 1}, {1, 1}, {2, 2}}, e4[][2] = {{0, 0}, {1, 2}, {2, 3}, {3, 4}},
                 e5[][2] = {{0, 1}, {1, 2}, {2, 3}};
    size_t rows = vector ? 5 : different || blocks ? 4 : 3;
    ck_sparse m, other;
    CHECK(sparse_new(&m, blocks ? 4 : 3, rows, width, order,
                     different ? e3
                     : blocks  ? e5
                               : e1,
                     different || blocks ? 3 : 2, 0) == CK_OK);
    CHECK(sparse_new(&other,
                     different ? 5
                     : blocks  ? 4
                               : 3,
                     rows, width, order,
                     different ? e4
                     : blocks  ? e5
                               : e2,
                     different ? 4 : 3, 0) == CK_OK);
    for (size_t b = 0; b < rows; ++b)
    {
      S(m, b, 0, 1) = different ? b * 2 + 1 : blocks || vector ? b + 1 : b + 10;
      S(m, b, different ? 1 : 1, different ? 1 : 2) = different ? b * 3 + 2
                                                      : blocks  ? b + 10
                                                      : vector  ? 0
                                                                : b * 2 + 20;
      if (different)
        S(m, b, 2, 2) = b * 5 + 3;
      if (different)
      {
        S(other, b, 0, 0) = b + 10;
        S(other, b, 1, 2) = b + 20;
        S(other, b, 2, 3) = b + 30;
        S(other, b, 3, 4) = b + 40;
      }
      else if (blocks)
        S(other, b, 2, 3) = b + 100;
      else
      {
        S(other, b, 0, 0) = b * 4;
        S(other, b, 1, 1) = b * 3;
        S(other, b, 2, 2) = b * 5;
      }
    }
    double values[] = {10, 20, 30, 40, 50};
    ck_view v;
    CHECK(ck_vector_view(values, rows, sizeof(double), 0, &v) == CK_OK);
    program p = {0};
    if (different)
      step(&p, '+', at(0, N(m, 2, 2)), 3, at(0, N(m, 0, 1)), at(0, N(m, 1, 1)), at(1, N(other, 0, 0)), unused,
           1);
    else if (blocks)
      step(&p, '+', at(0, N(m, 2, 3)), 3, at(0, N(m, 0, 1)), at(0, N(m, 1, 2)), at(1, N(other, 2, 3)), unused,
           1);
    else
      step(&p, '+', at(0, N(m, 1, 2)), 2, at(0, N(m, 0, 1)), vector ? at(1, 0) : at(1, N(other, 2, 2)),
           unused, unused, 1);
    ck_argument args[] = {ck_sparse_argument(&m, 1),
                          vector ? ck_vector_argument(v) : ck_sparse_argument(&other, 0)};
    if (blocks || vector)
    {
      ck_sparse proto;
      CHECK(sparse_new(&proto, blocks ? 4 : 3, 3, width, order, blocks ? e5 : e1, blocks ? 3 : 2, 0) ==
            CK_OK);
      ck_argument prototypes[] = {ck_sparse_argument(&proto, 1),
                                  blocks ? ck_sparse_argument(&proto, 0) : ck_vector_argument(v)};
      ck_array_function f;
      CHECK(ck_array_function_create(prototypes, 2, execute, &p, &f) == CK_OK);
      CHECK(ck_array_function_invoke(&f, args, 2) == CK_OK);
      ck_array_function_free(&f);
      ck_sparse_free(&proto);
    }
    else
      CHECK(execute_once(&p, args, 2) == CK_OK);
    for (size_t b = 0; b < rows; ++b)
    {
      double expected = different ? 13 + 6 * b : blocks ? 111 + 3 * b : vector ? 11 * (b + 1) : 10 + 6 * b;
      NEAR(S(m, b, blocks ? 2 : different ? 2 : 1, blocks ? 3 : 2), expected, 0);
      NEAR(S(m, b, 0, 1), different ? 2 * b + 1 : blocks || vector ? b + 1 : b + 10, 0);
    }
    ck_sparse_free(&m);
    ck_sparse_free(&other);
    return 0;
  }
  if (IS("MultipleTemporaries") || IS("BlockViewReuse") || IS("FunctionReusability") ||
      IS("ConstSparseMatrixFunction") || IS("EmptySparseMatrixFunction"))
  {
    int temp = IS("MultipleTemporaries"), reuse = IS("FunctionReusability"),
        empty = IS("EmptySparseMatrixFunction"), constant = IS("ConstSparseMatrixFunction");
    const size_t entries[][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 4}},
                 four[][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 3}}, diag[][2] = {{0, 0}, {1, 1}, {2, 2}};
    size_t blocks = temp    ? 4
                    : reuse ? 2
                            : 3,
           ne = temp    ? 5
                : reuse ? 3
                : empty ? 0
                        : 4,
           n = temp    ? 5
               : reuse ? 3
                       : 4;
    ck_sparse m;
    CHECK(sparse_new(&m, n, blocks, width, order, temp ? entries : reuse ? diag : four, ne, 0) == CK_OK);
    for (size_t b = 0; b < blocks && !empty; ++b)
    {
      S(m, b, 0, reuse ? 0 : 1) = reuse ? b : b + 1;
      S(m, b, 1, reuse ? 1 : 2) = temp ? (b + 1) * 10.0 : reuse ? b + 1 : constant ? b + 10 : 0;
      if (reuse)
        S(m, b, 2, 2) = b + 2;
    }
    program p = {0};
    if (temp)
    {
      p.temporaries = 2;
      step(&p, '*', at(1, 0), 2, at(0, N(m, 0, 1)), at(0, N(m, 1, 2)), unused, unused, 1);
      step(&p, '+', at(2, 0), 2, at(0, N(m, 0, 1)), at(0, N(m, 1, 2)), unused, unused, 1);
      step(&p, '+', at(0, N(m, 2, 3)), 2, at(1, 0), at(2, 0), unused, unused, 1);
      step(&p, '-', at(0, N(m, 3, 4)), 2, at(1, 0), at(2, 0), unused, unused, 1);
      step(&p, '*', at(0, N(m, 4, 4)), 2, at(1, 0), at(2, 0), unused, unused, 1);
    }
    else if (reuse)
      step(&p, '+', at(0, N(m, 2, 2)), 3, at(0, N(m, 0, 0)), at(0, N(m, 1, 1)), at(0, N(m, 2, 2)), unused, 2);
    else if (!empty && !constant)
    {
      step(&p, '+', at(0, N(m, 1, 2)), 1, at(0, N(m, 0, 1)), unused, unused, unused, 2);
      step(&p, '+', at(0, N(m, 2, 3)), 2, at(0, N(m, 0, 1)), at(0, N(m, 1, 2)), unused, unused, 1);
      step(&p, '*', at(0, N(m, 3, 3)), 2, at(0, N(m, 2, 3)), at(0, N(m, 1, 2)), unused, unused, 1);
    }
    ck_argument args[] = {ck_sparse_argument(&m, !constant)};
    ck_array_function f;
    CHECK(ck_array_function_create(args, 1, execute, &p, &f) == CK_OK);
    CHECK(ck_array_function_invoke(&f, args, 1) == CK_OK);
    for (size_t b = 0; b < blocks && !empty; ++b)
    {
      double x = b + 1.0, y = x * 10;
      if (temp)
      {
        NEAR(S(m, b, 2, 3), x * y + x + y, 0);
        NEAR(S(m, b, 3, 4), x * y - (x + y), 0);
        NEAR(S(m, b, 4, 4), x * y * (x + y), 0);
      }
      else if (reuse)
        NEAR(S(m, b, 2, 2), 2 * (3 * b + 3), 0);
      else if (constant)
      {
        NEAR(S(m, b, 0, 1), b + 1, 0);
        NEAR(S(m, b, 1, 2), b + 10, 0);
      }
      else
      {
        NEAR(S(m, b, 1, 2), 2 * x, 0);
        NEAR(S(m, b, 2, 3), 3 * x, 0);
        NEAR(S(m, b, 3, 3), 6 * x * x, 0);
      }
    }
    if (reuse)
    {
      double value = 5;
      ck_matrix_fill(&m.values, &value);
      CHECK(ck_array_function_invoke(&f, args, 1) == CK_OK);
      NEAR(S(m, 0, 2, 2), 30, 0);
      NEAR(S(m, 1, 2, 2), 30, 0);
      value = 0;
      ck_matrix_fill(&m.values, &value);
      S(m, 1, 0, 0) = 10;
      CHECK(ck_array_function_invoke(&f, args, 1) == CK_OK);
      NEAR(S(m, 0, 2, 2), 0, 0);
      NEAR(S(m, 1, 2, 2), 20, 0);
    }
    if (empty)
    {
      ck_sparse zero;
      CHECK(sparse_new(&zero, 4, 0, width, order, four, 1, 0) == CK_OK);
      ck_argument a[] = {ck_sparse_argument(&zero, 0)};
      CHECK(execute_once(&p, a, 1) == CK_OK);
      ck_sparse_free(&zero);
    }
    ck_array_function_free(&f);
    ck_sparse_free(&m);
    return 0;
  }
  return 2;
}
static int sparse_validation(const char *test, size_t width, ck_sparse_order order)
{
  const size_t entries[][2] = {{0, 1}, {1, 1}, {2, 2}}, chain[][2] = {{0, 1}, {1, 2}, {2, 3}};
  slot unused = {0};
  if (IS("MismatchedElementDimensions"))
  {
    ck_sparse m;
    CHECK(sparse_new(&m, 4, 3, width, order, entries, 3, 0) == CK_OK);
    CHECK(ck_sparse_index(&m, 0, 3, 3) == SIZE_MAX);
    ck_view view;
    CHECK(ck_sparse_block_view(&m, SIZE_MAX, 1, &view) == CK_INVALID);
    ck_sparse_free(&m);
    return 0;
  }
  if (IS("IncompatibleOrdering") || IS("IncompatibleVectorOrdering") || IS("IncompatibleSparseOrdering"))
  {
    ck_sparse a, b;
    CHECK(sparse_new(&a, 4, 3, width, order, chain, 2, 0) == CK_OK);
    size_t other_width = IS("IncompatibleOrdering") ? 0 : width == 4 ? 2 : 4;
    CHECK(sparse_new(&b, 4, 3, other_width, order, chain, 2, 0) == CK_OK);
    ck_matrix dense;
    double zero = 0;
    CHECK(ck_matrix_create(3, 4, other_width, sizeof(double), &zero, &dense) == CK_OK);
    ck_argument args[] = {ck_sparse_argument(&a, 0), IS("IncompatibleSparseOrdering")
                                                         ? ck_sparse_argument(&b, 0)
                                                         : ck_matrix_argument(&dense, 0)};
    program p = {0};
    ck_array_function fn;
    ck_status status = ck_array_function_create(args, 2, execute, &p, &fn);
    if ((width ? width : 1) != (other_width ? other_width : 1))
      CHECK(status == CK_INVALID);
    else
    {
      CHECK(status == CK_OK);
      ck_array_function_free(&fn);
    }
    ck_sparse_free(&a);
    ck_sparse_free(&b);
    ck_matrix_free(&dense);
    return 0;
  }
  if (IS("MismatchedBlockDimensions") || IS("MismatchedBlocksAtInvocation") ||
      IS("MultipleSparseMatricesMismatchedBlocksAtInvocation"))
  {
    int three = IS("MultipleSparseMatricesMismatchedBlocksAtInvocation");
    ck_sparse a, b;
    CHECK(sparse_new(&a, three ? 4 : 3, 3, width, order, three ? chain : entries, three ? 3 : 2, 0) == CK_OK);
    CHECK(sparse_new(&b, three ? 4 : 3, 4, width, order, three ? chain : entries, three ? 3 : 2, 0) == CK_OK);
    program p = {0};
    step(&p, '+', at(0, 1), 1, at(1, 0), unused, unused, unused, 2);
    ck_argument args[] = {ck_sparse_argument(&a, 1), ck_sparse_argument(&a, 0), ck_sparse_argument(&a, 0)};
    size_t count = three ? 3 : 2;
    ck_array_function f;
    CHECK(ck_array_function_create(args, count, execute, &p, &f) == CK_OK);
    CHECK(ck_array_function_invoke(&f, args, count) == CK_OK);
    args[three ? 2 : 1] = ck_sparse_argument(&b, 0);
    CHECK(ck_array_function_invoke(&f, args, count) == CK_INVALID);
    args[three ? 2 : 1] = ck_sparse_argument(&a, 0);
    args[0] = ck_sparse_argument(&b, 1);
    CHECK(ck_array_function_invoke(&f, args, count) == CK_INVALID);
    ck_array_function_free(&f);
    ck_sparse_free(&a);
    ck_sparse_free(&b);
    return 0;
  }
  size_t blocks = 3, vec_rows = 2, second_blocks = 3;
  int vectors = 1, matrices = 1, fail = 1;
  if (IS("VectorTooSmall"))
  {
  }
  else if (IS("VectorTooLarge"))
    vec_rows = 4;
  else if (IS("EmptyVectorNonEmptySparseMatrix"))
    vec_rows = 0;
  else if (IS("NonEmptyVectorEmptySparseMatrix"))
  {
    blocks = 0;
    vec_rows = 3;
  }
  else if (IS("EmptyVectorEmptySparseMatrix"))
  {
    blocks = vec_rows = 0;
    fail = 0;
  }
  else if (IS("MultipleVectorsDifferentSizes"))
  {
    vectors = 2;
    vec_rows = 3;
  }
  else if (IS("MultipleSparseMatricesDifferentBlocksVector"))
  {
    matrices = 2;
    vec_rows = 3;
    second_blocks = 4;
  }
  else if (IS("VectorSizeMatchesOneSparseMatrixOnly"))
  {
    matrices = vectors = 2;
    vec_rows = 3;
  }
  else if (IS("FunctionInvocationWithWrongSizedVectorSparse"))
    vec_rows = 2;
  else
    return 2;
  ck_sparse a, b;
  CHECK(sparse_new(&a, 3, blocks, width, order, entries, 1, 0) == CK_OK);
  CHECK(sparse_new(&b, 3, second_blocks, width, order, entries, 1, 0) == CK_OK);
  double values[] = {1, 2, 3, 4};
  ck_view v1, v2;
  CHECK(ck_vector_view(values, vec_rows, sizeof(double), 0, &v1) == CK_OK);
  CHECK(ck_vector_view(values, 2, sizeof(double), 0, &v2) == CK_OK);
  ck_argument args[4];
  args[0] = ck_sparse_argument(&a, 1);
  if (matrices == 2)
    args[1] = ck_sparse_argument(&b, 0);
  args[matrices] = ck_vector_argument(v1);
  if (vectors == 2)
    args[matrices + 1] = ck_vector_argument(v2);
  program p = {0};
  ck_array_function f;
  size_t count = (size_t)(matrices + vectors);
  CHECK(ck_array_function_create(args, count, execute, &p, &f) == CK_OK);
  CHECK(ck_array_function_invoke(&f, args, count) == (fail ? CK_INVALID : CK_OK));
  ck_array_function_free(&f);
  ck_sparse_free(&a);
  ck_sparse_free(&b);
  return 0;
}
static int sparse_vectors(const char *test, size_t width, ck_sparse_order order)
{
  int plain = IS("VectorInSparseMatrixFunction"), constant = IS("ConstVectorSparse"),
      mutable = IS("MutableVectorSparse"), reuse = IS("FunctionReusabilityWithVectorsSparse"),
      multi = IS("MultipleVectorsSameSize"), two = IS("MultipleSparseMatricesOneVector"),
      mixed = IS("MixedVectorBlockViewBlockVariable"), integer = IS("IntegerVectorSparse"),
      signature = IS("FunctionWithConstSignatureSparse");
  if (!plain && !constant && !mutable && !reuse && !multi && !two && !mixed && !integer && !signature)
    return 2;
  const size_t entries[][2] = {{0, 1}, {1, 2}, {2, 3}}, zero[][2] = {{0, 0}};
  ck_sparse m, other;
  CHECK(sparse_new(&m,
                   signature        ? 2
                   : plain || mixed ? 4
                                    : 3,
                   3, width, order, signature ? zero : entries, plain || mixed ? 3 : 1, 0) == CK_OK);
  CHECK(sparse_new(&other, 3, 3, width, order, entries, 1, 0) == CK_OK);
  double vec[] = {5, 10, 15}, vec2[] = {4, 5, 6};
  int iv[] = {5, 10, 15};
  if (reuse || multi || signature)
    for (size_t j = 0; j < 3; ++j)
      vec[j] = j + 1;
  if (two)
    for (size_t j = 0; j < 3; ++j)
      vec[j] = (j + 1) * 10.0;
  if (mixed)
    for (size_t j = 0; j < 3; ++j)
      vec[j] = (j + 1) * 100.0;
  for (size_t b = 0; b < 3; ++b)
  {
    S(m, b, 0, signature ? 0 : 1) = b + 1;
    if (plain || mixed)
    {
      S(m, b, 1, 2) = b + 2;
      S(m, b, 2, 3) = b + 3;
    }
  }
  ck_view v, v2;
  CHECK(ck_vector_view(integer ? (void *)iv : (void *)vec, 3, integer ? sizeof(int) : sizeof(double), mutable,
                       &v) == CK_OK);
  CHECK(ck_vector_view(vec2, 3, sizeof(double), 0, &v2) == CK_OK);
  slot unused = {0};
  program p = {0};
  size_t target = N(m, 0, signature ? 0 : 1);
  if (mixed)
    step(&p, '+', at(0, N(m, 2, 3)), 3, at(0, N(m, 0, 1)), at(0, N(m, 1, 2)), at(1, 0), unused, 2);
  else if (two)
    step(&p, '+', at(1, N(other, 0, 1)), 2, at(0, N(m, 0, 1)), at(2, 0), unused, unused, 1);
  else if (multi)
    step(&p, '+', at(0, target), 2, at(1, 0), at(2, 0), unused, unused, 1);
  else if (mutable)
    step(&p, '+', at(1, 0), 1, at(0, target), unused, unused, unused, 3);
  else
    step(&p, '+', at(0, target), 1, at(1, 0), unused, unused, unused, constant ? 3 : integer ? 1.5 : 2);
  ck_argument args[] = {ck_sparse_argument(&m, !mutable),
                        two ? ck_sparse_argument(&other, 1) : ck_vector_argument(v),
                        two ? ck_vector_argument(v) : ck_vector_argument(v2)};
  size_t count = multi || two ? 3 : 2;
  ck_array_function f;
  CHECK(ck_array_function_create(args, count, execute, &p, &f) == CK_OK);
  CHECK(ck_array_function_invoke(&f, args, count) == CK_OK);
  for (size_t b = 0; b < 3; ++b)
  {
    if (mutable)
      NEAR(vec[b], 3 * (b + 1), 0);
    else if (mixed)
      NEAR(S(m, b, 2, 3), 206 + 204 * b, 0);
    else if (two)
      NEAR(S(other, b, 0, 1), 11 * (b + 1), 0);
    else if (multi)
      NEAR(S(m, b, 0, 1), 5 + 2 * b, 0);
    else
      NEAR(S(m, b, 0, signature ? 0 : 1), (integer ? iv[b] : vec[b]) * (constant ? 3 : integer ? 1.5 : 2), 0);
  }
  if (reuse)
  {
    double second[] = {10, 20, 30};
    CHECK(ck_vector_view(second, 3, sizeof(double), 0, &v) == CK_OK);
    args[0] = ck_sparse_argument(&other, 1);
    args[1] = ck_vector_argument(v);
    CHECK(ck_array_function_invoke(&f, args, 2) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
      NEAR(S(other, b, 0, 1), 20 * (b + 1), 0);
  }
  ck_array_function_free(&f);
  ck_sparse_free(&m);
  ck_sparse_free(&other);
  return 0;
}
static int sparse_dense(const char *test, size_t width, ck_sparse_order order)
{
  if (!IS("SparseAndDenseMatrixFunction") && !IS("SparseAndVectorMatrixFunction"))
    return 2;
  const size_t entries[][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 3}};
  ck_sparse m;
  CHECK(sparse_new(&m, 4, 3, width, order, entries, 4, 0) == CK_OK);
  ck_matrix d;
  double zero = 0;
  CHECK(ck_matrix_create(3, 4, width, sizeof(double), &zero, &d) == CK_OK);
  for (size_t b = 0; b < 3; ++b)
  {
    for (size_t j = 0; j < 4; ++j)
    {
      S(m, b, entries[j][0], entries[j][1]) = b + (j == 0 ? 1 : j * 10);
      *(double *)ck_matrix_at(&d, b, j) = b * 100 + j * 10;
    }
  }
  program p = {0};
  step(&p, '+', at(0, N(m, 3, 3)), 4, at(0, N(m, 0, 1)), at(0, N(m, 1, 2)), at(1, 0), at(1, 1), 1);
  ck_argument a[] = {ck_sparse_argument(&m, 1), ck_matrix_argument(&d, 0)};
  CHECK(execute_once(&p, a, 2) == CK_OK);
  for (size_t b = 0; b < 3; ++b)
  {
    NEAR(S(m, b, 3, 3), 21 + 202 * b, 0);
    NEAR(S(m, b, 0, 1), b + 1, 0);
    NEAR(S(m, b, 1, 2), b + 10, 0);
    NEAR(S(m, b, 2, 3), b + 20, 0);
    for (size_t j = 0; j < 4; ++j)
      NEAR(*(double *)ck_matrix_at(&d, b, j), b * 100 + j * 10, 0);
  }
  ck_sparse_free(&m);
  ck_matrix_free(&d);
  return 0;
}
static int sparse_views(const char *test, size_t width, ck_sparse_order order)
{
  if (IS("GetBlockViewByVectorIndex"))
  {
    const size_t entries[][2] = {{0, 0}, {0, 2}, {1, 0}, {1, 1}, {2, 1}, {2, 2}};
    ck_sparse m;
    CHECK(sparse_new(&m, 3, 3, width, order, entries, 6, 0) == CK_OK);
    for (size_t j = 0; j < 6; ++j)
    {
      ck_view view;
      CHECK(ck_sparse_block_view(&m, N(m, entries[j][0], entries[j][1]), 1, &view) == CK_OK);
      for (size_t b = 0; b < 3; ++b)
        *(double *)ck_view_at(&view, b) = b + j * 10;
    }
    for (size_t j = 0; j < 6; ++j)
      for (size_t b = 0; b < 3; ++b)
        NEAR(S(m, b, entries[j][0], entries[j][1]), b + j * 10, 0);
    ck_sparse_free(&m);
    return 0;
  }
  if (!IS("TestFill") && !IS("TestCopy"))
    return 2;
  const size_t entries[][2] = {{0, 1}, {2, 0}}, other_entries[][2] = {{1, 2}};
  ck_sparse m, other;
  CHECK(sparse_new(&m, 3, 3, width, order, entries, 2, 0) == CK_OK);
  CHECK(sparse_new(&other, 3, 3, width, order, other_entries, 1, 0) == CK_OK);
  ck_view v01, v20, raw, temp, output, other_view;
  double values[] = {3.2, 4.2, 1.3}, scratch[3] = {0}, out[3] = {0};
  CHECK(ck_sparse_block_view(&m, N(m, 0, 1), 1, &v01) == CK_OK);
  CHECK(ck_sparse_block_view(&m, N(m, 2, 0), 1, &v20) == CK_OK);
  CHECK(ck_sparse_block_view(&other, N(other, 1, 2), 1, &other_view) == CK_OK);
  CHECK(ck_vector_view(values, 3, sizeof(double), 1, &raw) == CK_OK);
  CHECK(ck_vector_view(scratch, 3, sizeof(double), 1, &temp) == CK_OK);
  CHECK(ck_vector_view(out, 3, sizeof(double), 1, &output) == CK_OK);
  if (IS("TestFill"))
  {
    double val = 3.2;
    CHECK(ck_view_fill(v01, &val, 1) == CK_OK);
    CHECK(ck_view_fill(raw, &val, 1) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
    {
      NEAR(S(m, b, 0, 1), 3.2, 0);
      NEAR(S(m, b, 2, 0), 0, 0);
      NEAR(values[b], 3.2, 0);
    }
    val = 9.9;
    CHECK(ck_view_fill(temp, &val, 1) == CK_OK);
    CHECK(ck_view_copy(v20, temp, 1) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
      NEAR(S(m, b, 2, 0), 9.9, 0);
  }
  else
  {
    double ten[] = {10, 20, 30};
    ck_view input;
    CHECK(ck_vector_view(ten, 3, sizeof(double), 0, &input) == CK_OK);
    CHECK(ck_view_copy(v01, input, 1) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
    {
      NEAR(S(m, b, 0, 1), ten[b], 0);
      NEAR(S(m, b, 2, 0), 0, 0);
    }
    CHECK(ck_view_copy(v01, raw, 1) == CK_OK);
    v01.writable = 0;
    CHECK(ck_view_copy(output, v01, 1) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
      NEAR(out[b], values[b], 0);
    v01.writable = 1;
    CHECK(ck_view_copy(v20, v01, 1) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
      NEAR(S(m, b, 2, 0), values[b], 0);
    double zero = 0;
    CHECK(ck_view_fill(v20, &zero, 1) == CK_OK);
    v01.writable = 0;
    CHECK(ck_view_copy(v20, v01, 1) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
      NEAR(S(m, b, 2, 0), values[b], 0);
    CHECK(ck_view_copy(temp, v01, 1) == CK_OK);
    v01.writable = 1;
    CHECK(ck_view_fill(v01, &zero, 1) == CK_OK);
    CHECK(ck_view_copy(v01, temp, 1) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
      NEAR(S(m, b, 0, 1), values[b], 0);
    ck_matrix dense;
    CHECK(ck_matrix_create(3, 2, width, sizeof(double), &zero, &dense) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
      *(double *)ck_matrix_at(&dense, b, 1) = (b + 1) * 100.0;
    ck_view d1, d0;
    CHECK(ck_matrix_column(&dense, 1, 0, &d1) == CK_OK);
    CHECK(ck_matrix_column(&dense, 0, 1, &d0) == CK_OK);
    CHECK(ck_view_copy(v20, d1, 1) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
      NEAR(S(m, b, 2, 0), (b + 1) * 100.0, 0);
    CHECK(ck_view_copy(temp, v01, 1) == CK_OK);
    CHECK(ck_view_copy(d0, temp, 1) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
      NEAR(*(double *)ck_matrix_at(&dense, b, 0), values[b], 0);
    CHECK(ck_view_copy(other_view, v01, 1) == CK_OK);
    for (size_t b = 0; b < 3; ++b)
      NEAR(S(other, b, 1, 2), values[b], 0);
    ck_matrix_free(&dense);
  }
  ck_sparse_free(&m);
  ck_sparse_free(&other);
  return 0;
}
int test_sparse(const char *suite, const char *test)
{
  if (!strcmp(suite, "SparseMatrixBuilder") && IS("BadConfiguration"))
  {
    const size_t rows[] = {3, 2, 6}, cols[] = {0, 4, 7};
    ck_sparse m;
    double zero = 0;
    for (size_t j = 0; j < 3; ++j)
      CHECK(ck_sparse_from_entries(3, 1, 0, CK_SPARSE_ROW, &rows[j], &cols[j], 1, sizeof(double), &zero,
                                   &m) == CK_INVALID);
    return 0;
  }
  if (strncmp(suite, "Sparse", 6) || !strstr(suite, "Matrix"))
    return 2;
  int vector = strstr(suite, "Vector") != NULL;
  ck_sparse_order order = strstr(suite, "Column") ? CK_SPARSE_COLUMN : CK_SPARSE_ROW;
  for (size_t width = vector ? 1 : 0; width <= (vector ? 4u : 0u); ++width)
  {
    int status = sparse_basic(test, width, order);
    if (status == 2)
      status = sparse_functions(test, width, order);
    if (status == 2)
      status = sparse_validation(test, width, order);
    if (status == 2)
      status = sparse_vectors(test, width, order);
    if (status == 2)
      status = sparse_dense(test, width, order);
    if (status == 2)
      status = sparse_views(test, width, order);
    if (status)
      return status;
  }
  return 0;
}
