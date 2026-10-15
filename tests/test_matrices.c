#include "test.h"
#include "chem_matrix.h"
#include <float.h>
#define D(m, r, c) (*(double *)ck_matrix_at(&(m), (r), (c)))
#define I(m, r, c) (*(int *)ck_matrix_at(&(m), (r), (c)))
#define NEW(m, r, c, w, v)                                                                                   \
  ck_matrix m;                                                                                               \
  double m##_initial = (v);                                                                                  \
  CHECK(ck_matrix_create((r), (c), (w), sizeof(double), &m##_initial, &m) == CK_OK)

enum action
{
  SUM3,
  MULTI,
  VECTOR_PRODUCT,
  VECTOR_SUM,
  VECTOR_SCALE,
  VECTOR_FROM_MATRIX,
  TWO_MATRIX_VECTOR,
  MIXED,
  TEMPORARIES,
  COLUMNS,
  READ_ONLY,
  NO_OP,
  BAD_COLUMN,
  INTEGER_SCALE
};
typedef struct
{
  enum action action;
  double factor;
} operation;
static ck_status arithmetic_row(ck_view *v, size_t count, size_t row, void *context)
{
  (void)count;
  operation *op = context;
  if (op->action == INTEGER_SCALE)
  {
    int *out = ck_view_at(&v[0], row);
    const int *in = ck_view_const_at(&v[1], row);
    if (!out || !in)
      return CK_INVALID;
    *out = *in * 2;
    return CK_OK;
  }
  double x[5] = {0};
  for (size_t j = 1; j < count; ++j)
  {
    const double *value = ck_view_const_at(&v[j], row);
    if (!value)
      return CK_INVALID;
    x[j] = *value;
  }
  double *out = ck_view_at(&v[0], row);
  if (!out)
    return CK_INVALID;
  switch (op->action)
  {
  case SUM3:
    *out = op->factor * (x[1] + x[2] + x[3]);
    break;
  case MULTI:
    *out = op->factor * (x[1] + x[2]);
    break;
  case VECTOR_PRODUCT:
    *out = x[1] * (x[2] + x[3] + x[4]);
    break;
  case VECTOR_SUM:
    *out = x[1] + x[2];
    break;
  case VECTOR_SCALE:
  case VECTOR_FROM_MATRIX:
    *out = op->factor * x[1];
    break;
  case MIXED:
    *out = (x[1] + x[2]) + x[3];
    break;
  default:
    return CK_INVALID;
  }
  return CK_OK;
}
static ck_status matrix_callback(const ck_argument *a, size_t count, void *context)
{
  operation *op = context;
  const ck_matrix *m = a[0].matrix;
  ck_view v[5];
  ck_status status;
  if (op->action == NO_OP)
    return CK_OK;
  if (op->action == BAD_COLUMN)
    return ck_matrix_column(m, 5, a[0].writable, &v[0]);
  if (op->action == READ_ONLY)
  {
    ck_view ro;
    if (ck_matrix_column(m, 0, 0, &ro))
      return CK_INVALID;
    if (ck_view_at(&ro, 0))
      return CK_INVALID;
    return CK_OK;
  }
  if (op->action == TEMPORARIES || op->action == COLUMNS)
  {
    size_t required = op->action == TEMPORARIES ? 5 : 4;
    for (size_t j = 0; j < required; ++j)
      if (ck_matrix_column(m, j, a[0].writable, &v[j]))
        return CK_INVALID;
    ck_matrix tmp1, tmp2;
    double zero = 0;
    if (ck_matrix_create(m->rows, 1, m->width, sizeof(double), &zero, &tmp1))
      return CK_MEMORY;
    if (ck_matrix_create(m->rows, 1, m->width, sizeof(double), &zero, &tmp2))
    {
      ck_matrix_free(&tmp1);
      return CK_MEMORY;
    }
    for (size_t r = 0; r < m->rows; ++r)
    {
      double x = *(const double *)ck_view_const_at(&v[0], r), y = *(const double *)ck_view_const_at(&v[1], r);
      if (op->action == TEMPORARIES)
      {
        D(tmp1, r, 0) = x * y;
        D(tmp2, r, 0) = x + y;
        *(double *)ck_view_at(&v[2], r) = D(tmp1, r, 0) + D(tmp2, r, 0);
        *(double *)ck_view_at(&v[3], r) = D(tmp1, r, 0) - D(tmp2, r, 0);
        *(double *)ck_view_at(&v[4], r) = D(tmp1, r, 0) * D(tmp2, r, 0);
      }
      else
      {
        *(double *)ck_view_at(&v[1], r) = x * 2;
        *(double *)ck_view_at(&v[2], r) = x + *(const double *)ck_view_const_at(&v[1], r);
        *(double *)ck_view_at(&v[3], r) =
            *(const double *)ck_view_const_at(&v[2], r) * *(const double *)ck_view_const_at(&v[1], r);
      }
    }
    ck_matrix_free(&tmp1);
    ck_matrix_free(&tmp2);
    return CK_OK;
  }
  if (op->action == TWO_MATRIX_VECTOR)
  {
    for (size_t j = 0; j < 2; ++j)
    {
      if (ck_matrix_column(a[j].matrix, 1, a[j].writable, &v[0]) ||
          ck_matrix_column(a[j].matrix, 0, 0, &v[1]))
        return CK_INVALID;
      v[2] = a[2].vector;
      operation sum = {VECTOR_SUM, 1};
      status = ck_view_apply(v, 3, 1, arithmetic_row, &sum);
      if (status)
        return status;
    }
    return CK_OK;
  }
  size_t nviews = 2;
  if (op->action == VECTOR_FROM_MATRIX)
  {
    v[0] = a[1].vector;
    if (ck_matrix_column(m, 0, 0, &v[1]))
      return CK_INVALID;
  }
  else
  {
    size_t target = (op->action == SUM3 || op->action == VECTOR_PRODUCT || op->action == MIXED) ? 2
                    : op->action == MULTI                                                       ? 1
                                                                                                : 0;
    if (ck_matrix_column(m, target, a[0].writable, &v[0]))
      return CK_INVALID;
    if (op->action == SUM3)
    {
      nviews = 4;
      for (size_t j = 0; j < 3; ++j)
        if (ck_matrix_column(m, j, 0, &v[j + 1]))
          return CK_INVALID;
    }
    else if (op->action == MULTI)
    {
      nviews = 3;
      if (ck_matrix_column(m, 0, 0, &v[1]))
        return CK_INVALID;
      if (count > 1 && a[1].kind == CK_ARGUMENT_MATRIX)
      {
        if (ck_matrix_column(a[1].matrix, 2, 0, &v[2]))
          return CK_INVALID;
      }
      else
        v[2] = a[1].vector;
    }
    else if (op->action == VECTOR_PRODUCT)
    {
      nviews = 5;
      v[1] = a[1].vector;
      for (size_t j = 0; j < 3; ++j)
        if (ck_matrix_column(m, j, 0, &v[j + 2]))
          return CK_INVALID;
    }
    else if (op->action == VECTOR_SUM)
    {
      nviews = 3;
      v[1] = a[1].vector;
      v[2] = a[2].vector;
    }
    else if (op->action == MIXED)
    {
      nviews = 4;
      if (ck_matrix_column(m, 0, 0, &v[1]) || ck_matrix_column(m, 1, 0, &v[2]))
        return CK_INVALID;
      v[3] = a[1].vector;
    }
    else
      v[1] = a[1].vector;
  }
  return ck_view_apply(v, nviews, 1, arithmetic_row, op);
}
static int invoke(operation *op, const ck_argument *args, size_t count)
{
  ck_array_function function;
  CHECK(ck_array_function_create(args, count, matrix_callback, op, &function) == CK_OK);
  ck_status status = ck_array_function_invoke(&function, args, count);
  ck_array_function_free(&function);
  return status;
}
static int dense_basic(const char *test, size_t width)
{
  if (IS("SmallMatrix") || IS("SmallConstMatrix") || IS("SmallVectorMatrix") || IS("SmallConstVectorMatrix"))
  {
    NEW(m, 3, 5, width, 0);
    D(m, 1, 3) = 64.7;
    D(m, 0, 0) = 41.2;
    D(m, 2, 4) = 102.3;
    const ck_matrix *view = &m;
    NEAR(*(const double *)ck_matrix_const_at(view, 1, 3), 64.7, 0);
    NEAR(D(m, 0, 0), 41.2, 0);
    NEAR(D(m, 2, 4), 102.3, 0);
    CHECK(m.size == (width ? ((3 + width - 1) / width) * width * 5 : 15));
    double *data = m.data;
    NEAR(data[0], 41.2, 0);
    size_t index = width ? (2 / width) * width * 5 + 4 * width + 2 % width : 14;
    NEAR(data[index], 102.3, 0);
    index = width ? (1 / width) * width * 5 + 3 * width + 1 % width : 8;
    NEAR(data[index], 64.7, 0);
    ck_matrix_free(&m);
    return 0;
  }
  if (IS("InitializeMatrix") || IS("InitializeConstMatrix") || IS("InitializeVectorMatrix") ||
      IS("InitializeConstVectorMatrix"))
  {
    NEW(m, 2, 3, width, 12.4);
    NEAR(D(m, 0, 0), 12.4, 0);
    NEAR(D(m, 1, 0), 12.4, 0);
    NEAR(D(m, 1, 2), 12.4, 0);
    ck_matrix_free(&m);
    return 0;
  }
  if (IS("LoopOverMatrix") || IS("LoopOverConstMatrix") || IS("LoopOverVectorMatrix") ||
      IS("LoopOverConstVectorMatrix") || IS("IterateOverMatrix"))
  {
    ck_matrix m;
    int zero = 0;
    CHECK(ck_matrix_create(3, 4, width, sizeof(int), &zero, &m) == CK_OK);
    if (IS("IterateOverMatrix"))
    {
      for (size_t j = 0; j < 4; ++j)
        I(m, 1, j) = 42 + (int)j;
      NEAR(I(m, 0, 3), 0, 0);
      NEAR(I(m, 1, 0), 42, 0);
      NEAR(I(m, 1, 1), 43, 0);
      NEAR(I(m, 1, 3), 45, 0);
      NEAR(I(m, 2, 0), 0, 0);
    }
    else
    {
      for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 4; ++j)
          I(m, i, j) = (int)(i * 100 + j);
      CHECK(I(m, 0, 0) == 0 && I(m, 1, 2) == 102 && I(m, 2, 3) == 203 && I(m, 0, 3) == 3);
    }
    ck_matrix_free(&m);
    return 0;
  }
  if (IS("Strides"))
  {
    if (width && width < 3)
      return 0;
    ck_matrix m;
    int zero = 0;
    CHECK(ck_matrix_create(3, 4, width, sizeof(int), &zero, &m) == CK_OK);
    size_t row_stride = width ? 1 : 4, column_stride = width ? width : 1;
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 4; ++j)
        ((int *)m.data)[i * row_stride + j * column_stride] = (int)(i * 100 + j);
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 4; ++j)
        CHECK(I(m, i, j) == (int)(i * 100 + j));
    ck_matrix_free(&m);
    return 0;
  }
  if (IS("ConversionToVector") || IS("ConstConversionToVector") || IS("AssignmentFromVector"))
  {
    NEW(m, IS("AssignmentFromVector") ? 4 : 2, 3, width, 0);
    double values[] = {IS("AssignmentFromVector") ? 12.3 : 13.2, IS("AssignmentFromVector") ? 15.1 : 31.2,
                       IS("AssignmentFromVector") ? 24.3 : 314.2},
           row[3];
    size_t r = IS("AssignmentFromVector") ? 2 : 1;
    CHECK(ck_matrix_set_row(&m, r, values, 3) == CK_OK);
    CHECK(ck_matrix_get_row(&m, r, row, 3) == CK_OK);
    for (size_t j = 0; j < 3; ++j)
      NEAR(row[j], values[j], 0);
    if (IS("AssignmentFromVector"))
    {
      double big[] = {14.3, 52.3, 65.7, 16.34}, small[] = {13.2, 52.8};
      CHECK(ck_matrix_set_row(&m, 2, big, 4) == CK_OK);
      for (size_t j = 0; j < 3; ++j)
        NEAR(D(m, 2, j), big[j], 0);
      CHECK(ck_matrix_set_row(&m, 2, small, 2) == CK_INVALID);
      NEAR(D(m, 0, 0), 0, 0);
      NEAR(D(m, 3, 0), 0, 0);
    }
    ck_matrix_free(&m);
    return 0;
  }
  if (IS("ConversionFromVector"))
  {
    ck_matrix m, empty;
    CHECK(ck_matrix_from_rows(0, NULL, NULL, width, sizeof(double), &empty) == CK_OK);
    CHECK(empty.rows == 0);
    ck_matrix_free(&empty);
    double row0[] = {412.3, 32.4, 41.3}, row1[] = {5.33, -0.3, 31.2};
    const void *rows[] = {row0, row1};
    size_t lengths[] = {3, 3};
    CHECK(ck_matrix_from_rows(2, rows, lengths, width, sizeof(double), &m) == CK_OK);
    CHECK(m.rows == 2 && m.columns == 3);
    for (size_t j = 0; j < 3; ++j)
    {
      NEAR(D(m, 0, j), row0[j], 0);
      NEAR(D(m, 1, j), row1[j], 0);
    }
    size_t bad[] = {1, 2};
    ck_matrix other;
    CHECK(ck_matrix_from_rows(2, rows, bad, width, sizeof(double), &other) == CK_INVALID);
    ck_matrix_free(&m);
    return 0;
  }
  if (IS("Axpy") || IS("ForEach"))
  {
    NEW(y, 4, 3, width, 0);
    NEW(x, 4, 3, width, 0);
    NEW(z, 4, 3, width, 0);
    double expected = 0;
    for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 3; ++j)
      {
        D(y, i, j) = i * 10.3 + j * 100.5;
        D(x, i, j) = i * 1.7 + j * 10.2;
        D(z, i, j) = i * 19.5 + j * 32.2;
        expected += D(y, i, j) + 1.39 * D(x, i, j);
      }
    if (IS("Axpy"))
    {
      CHECK(ck_matrix_axpy(&y, 1.39, &x) == CK_OK);
      double sum = 0;
      for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < 3; ++j)
          sum += D(y, i, j);
      NEAR(sum, expected, 1e-5);
    }
    else
    {
      CHECK(ck_matrix_axpy(&y, 1, &x) == CK_OK);
      for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < 3; ++j)
          NEAR(D(y, i, j), (i * 10.3 + j * 100.5) + (i * 1.7 + j * 10.2), 1e-5);
      CHECK(ck_matrix_axpy(&y, -1, &z) == CK_OK);
      for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < 3; ++j)
          NEAR(D(y, i, j), (i * 10.3 + j * 100.5) + (i * 1.7 + j * 10.2) - (i * 19.5 + j * 32.2), 1e-5);
    }
    ck_matrix_free(&y);
    ck_matrix_free(&x);
    ck_matrix_free(&z);
    return 0;
  }
  if (IS("SetScaler") || IS("Max") || IS("Min") || IS("Print"))
  {
    NEW(m, 2, 3, width, 0);
    double value = 2;
    if (IS("SetScaler"))
    {
      ck_matrix_fill(&m, &value);
      for (size_t i = 0; i < m.size; ++i)
        NEAR(((double *)m.data)[i], 2, 0);
    }
    else if (IS("Print"))
    {
      D(m, 1, 1) = 3;
      FILE *file = tmpfile();
      CHECK(file);
      CHECK(ck_matrix_print(&m, file) == CK_OK);
      rewind(file);
      char text[100] = {0};
      CHECK(fread(text, 1, 99, file) == 12);
      CHECK(!strcmp(text, "0,0,0\n0,3,0\n"));
      fclose(file);
    }
    else
    {
      int lower = IS("Max");
      ck_matrix_bound(&m, 2, lower);
      for (size_t i = 0; i < m.size; ++i)
        NEAR(((double *)m.data)[i], lower ? 2 : 0, 0);
      value = 1;
      ck_matrix_fill(&m, &value);
      D(m, 1, 1) = 3;
      ck_matrix_bound(&m, 2, lower);
      for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
          NEAR(D(m, i, j), i == 1 && j == 1 ? (lower ? 3 : 2) : (lower ? 2 : 1), 0);
    }
    ck_matrix_free(&m);
    return 0;
  }
  return 2;
}
static int dense_functions(const char *test, size_t width)
{
  if (IS("ArrayFunction") || IS("VectorInMatrixFunction") || IS("FunctionReusability"))
  {
    NEW(m, IS("FunctionReusability") ? 2 : 5, 3, width, IS("FunctionReusability") ? 1 : -1);
    for (size_t i = 0; i < m.rows; ++i)
      for (size_t j = 0; j < 3; ++j)
        D(m, i, j) = IS("FunctionReusability") ? (double)(i + j) : (double)i - 2 + 10 * j;
    double vector[] = {100, 200, 300, 400, 500};
    ck_view vv;
    CHECK(ck_vector_view(vector, m.rows, sizeof(double), 0, &vv) == CK_OK);
    operation op = {IS("VectorInMatrixFunction") ? VECTOR_PRODUCT : SUM3, IS("FunctionReusability") ? 2 : 4};
    ck_argument args[] = {ck_matrix_argument(&m, 1), ck_vector_argument(vv)};
    ck_array_function fn;
    CHECK(ck_array_function_create(args, IS("VectorInMatrixFunction") ? 2 : 1, matrix_callback, &op, &fn) ==
          CK_OK);
    CHECK(ck_array_function_invoke(&fn, args, IS("VectorInMatrixFunction") ? 2 : 1) == CK_OK);
    for (size_t i = 0; i < m.rows; ++i)
    {
      double sum = IS("FunctionReusability") ? 3 * i + 3 : 3 * ((double)i - 2) + 30;
      NEAR(D(m, i, 2), (IS("VectorInMatrixFunction") ? vector[i] : op.factor) * sum, 0);
      NEAR(D(m, i, 0), IS("FunctionReusability") ? (double)i : (double)i - 2, 0);
      NEAR(D(m, i, 1), IS("FunctionReusability") ? i + 1.0 : (double)i + 8, 0);
    }
    if (IS("ArrayFunction") || IS("FunctionReusability"))
    {
      NEW(other, IS("ArrayFunction") ? 3 : 2, 3, width, IS("ArrayFunction") ? -1 : 5);
      ck_argument second[] = {ck_matrix_argument(&other, 1)};
      CHECK(ck_array_function_invoke(&fn, second, 1) == CK_OK);
      for (size_t i = 0; i < other.rows; ++i)
        NEAR(D(other, i, 2), IS("ArrayFunction") ? -12 : 30, 0);
      if (IS("FunctionReusability"))
      {
        double zero = 0;
        ck_matrix_fill(&other, &zero);
        D(other, 1, 0) = 10;
        CHECK(ck_array_function_invoke(&fn, second, 1) == CK_OK);
        NEAR(D(other, 0, 2), 0, 0);
        NEAR(D(other, 1, 2), 20, 0);
      }
      ck_matrix_free(&other);
    }
    ck_array_function_free(&fn);
    ck_matrix_free(&m);
    return 0;
  }
  if (IS("MultiMatrixArrayFunction") || IS("MultiMatrixDifferentRowsFromCreation") ||
      IS("MatrixVectorDifferentRowsFromCreation"))
  {
    int vector = IS("MatrixVectorDifferentRowsFromCreation"), different = !IS("MultiMatrixArrayFunction");
    NEW(proto, 3, vector ? 3 : 2, width, 0);
    NEW(other_proto, 3, 3, width, 0);
    double proto_vec[3] = {0};
    ck_view pv;
    CHECK(ck_vector_view(proto_vec, 3, sizeof(double), 0, &pv) == CK_OK);
    operation op = {MULTI, vector ? 3 : different ? 2 : 1};
    ck_argument prototypes[] = {ck_matrix_argument(&proto, 1),
                                vector ? ck_vector_argument(pv) : ck_matrix_argument(&other_proto, 0)};
    ck_array_function fn;
    CHECK(ck_array_function_create(prototypes, 2, matrix_callback, &op, &fn) == CK_OK);
    size_t rows = different ? 5 : 3;
    NEW(a, rows, vector ? 3 : 2, width, 0);
    NEW(b, rows, 3, width, 0);
    double vec[5] = {0};
    for (size_t i = 0; i < rows; ++i)
    {
      D(a, i, 0) = different ? i + 1.0 : i;
      D(a, i, 1) = different ? 0 : i + 10.0;
      D(b, i, 2) = different ? i * 10.0 : i * 4.0;
      vec[i] = i * 10.0;
    }
    ck_view v;
    CHECK(ck_vector_view(vec, rows, sizeof(double), 0, &v) == CK_OK);
    ck_argument args[] = {ck_matrix_argument(&a, 1),
                          vector ? ck_vector_argument(v) : ck_matrix_argument(&b, 0)};
    CHECK(ck_array_function_invoke(&fn, args, 2) == CK_OK);
    for (size_t i = 0; i < rows; ++i)
      NEAR(D(a, i, 1), different ? (i + 1 + i * 10) * op.factor : i * 5.0, 0);
    if (IS("MultiMatrixDifferentRowsFromCreation"))
    {
      NEW(a2, 2, 2, width, 0);
      NEW(b2, 2, 3, width, 0);
      D(a2, 0, 0) = 10;
      D(a2, 1, 0) = 20;
      D(b2, 0, 2) = 5;
      D(b2, 1, 2) = 15;
      ck_argument second[] = {ck_matrix_argument(&a2, 1), ck_matrix_argument(&b2, 0)};
      CHECK(ck_array_function_invoke(&fn, second, 2) == CK_OK);
      NEAR(D(a2, 0, 1), 30, 0);
      NEAR(D(a2, 1, 1), 70, 0);
      ck_matrix_free(&a2);
      ck_matrix_free(&b2);
    }
    ck_array_function_free(&fn);
    ck_matrix_free(&proto);
    ck_matrix_free(&other_proto);
    ck_matrix_free(&a);
    ck_matrix_free(&b);
    return 0;
  }
  if (IS("MultipleTemporaries") || IS("ColumnViewReuse") || IS("ConstMatrixFunction") ||
      IS("EmptyMatrixFunction"))
  {
    size_t columns = IS("MultipleTemporaries")                            ? 5
                     : IS("ColumnViewReuse") || IS("ConstMatrixFunction") ? 4
                                                                          : 3,
           rows = IS("ColumnViewReuse") || IS("ConstMatrixFunction") ? 3
                  : IS("EmptyMatrixFunction")                        ? 0
                                                                     : 4;
    NEW(m, rows, columns, width, 0);
    for (size_t i = 0; i < rows; ++i)
    {
      D(m, i, 0) = IS("ConstMatrixFunction") ? i * 10.0 : i + 1.0;
      D(m, i, 1) = (i + 1) * 10.0;
      if (IS("ConstMatrixFunction"))
        for (size_t j = 0; j < columns; ++j)
          D(m, i, j) = i * 10.0 + j;
    }
    operation op = {IS("MultipleTemporaries")   ? TEMPORARIES
                    : IS("ColumnViewReuse")     ? COLUMNS
                    : IS("ConstMatrixFunction") ? READ_ONLY
                                                : NO_OP,
                    1};
    ck_argument args[] = {ck_matrix_argument(&m, op.action != READ_ONLY)};
    CHECK(invoke(&op, args, 1) == CK_OK);
    for (size_t i = 0; i < rows; ++i)
    {
      double x = i + 1.0, y = x * 10;
      if (op.action == TEMPORARIES)
      {
        NEAR(D(m, i, 2), x * y + x + y, 0);
        NEAR(D(m, i, 3), x * y - (x + y), 0);
        NEAR(D(m, i, 4), x * y * (x + y), 0);
      }
      else if (op.action == COLUMNS)
      {
        NEAR(D(m, i, 1), 2 * x, 0);
        NEAR(D(m, i, 2), 3 * x, 0);
        NEAR(D(m, i, 3), 6 * x * x, 0);
      }
      else
        for (size_t j = 0; j < columns; ++j)
          NEAR(D(m, i, j), i * 10.0 + j, 0);
    }
    if (op.action == NO_OP)
    {
      NEW(zero_columns, 3, 0, width, 1);
      ck_argument empty[] = {ck_matrix_argument(&zero_columns, 1)};
      CHECK(invoke(&op, empty, 1) == CK_OK);
      ck_matrix_free(&zero_columns);
    }
    ck_matrix_free(&m);
    return 0;
  }
  return 2;
}
static int dense_validation(const char *test, size_t width)
{
  if (IS("MismatchedColumnDimensions"))
  {
    NEW(m, 3, 4, width, 1);
    operation op = {BAD_COLUMN, 0};
    ck_argument a[] = {ck_matrix_argument(&m, 1)};
    CHECK(invoke(&op, a, 1) == CK_INVALID);
    ck_matrix_free(&m);
    return 0;
  }
  if (IS("WrongMatrixDimensions") || IS("WrongColumnCountAtInvocation"))
  {
    size_t columns = IS("WrongMatrixDimensions") ? 4 : 3;
    NEW(proto, 4, columns, width, 1);
    operation op = {NO_OP, 0};
    ck_argument a[] = {ck_matrix_argument(&proto, 1)};
    ck_array_function fn;
    CHECK(ck_array_function_create(a, 1, matrix_callback, &op, &fn) == CK_OK);
    CHECK(ck_array_function_invoke(&fn, a, 1) == CK_OK);
    NEW(wrong, 4, columns + 1, width, 1);
    a[0] = ck_matrix_argument(&wrong, 1);
    CHECK(ck_array_function_invoke(&fn, a, 1) == CK_INVALID);
    ck_matrix_free(&wrong);
    if (IS("WrongColumnCountAtInvocation"))
    {
      NEW(shorter, 4, columns - 1, width, 1);
      a[0] = ck_matrix_argument(&shorter, 1);
      CHECK(ck_array_function_invoke(&fn, a, 1) == CK_INVALID);
      ck_matrix_free(&shorter);
    }
    NEW(ok, 7, columns, width, 1);
    a[0] = ck_matrix_argument(&ok, 1);
    CHECK(ck_array_function_invoke(&fn, a, 1) == CK_OK);
    ck_matrix_free(&ok);
    ck_matrix_free(&proto);
    ck_array_function_free(&fn);
    return 0;
  }
  int row_mismatch = IS("MismatchedRowDimensions"), multi = IS("MultipleMatricesMismatchedRowsAtInvocation");
  if (row_mismatch || multi)
  {
    NEW(a, 3, row_mismatch ? 3 : 2, width, 1);
    NEW(b, row_mismatch ? 4 : 3, 3, width, 2);
    operation op = {NO_OP, 0};
    ck_argument prototypes[] = {ck_matrix_argument(&a, 1), ck_matrix_argument(&b, 0)};
    ck_array_function fn;
    CHECK(ck_array_function_create(prototypes, 2, matrix_callback, &op, &fn) == CK_OK);
    if (row_mismatch)
    {
      NEW(c, 5, 3, width, 0);
      NEW(d, 5, 3, width, 1);
      ck_argument ok[] = {ck_matrix_argument(&c, 1), ck_matrix_argument(&d, 0)};
      CHECK(ck_array_function_invoke(&fn, ok, 2) == CK_OK);
      CHECK(ck_array_function_invoke(&fn, prototypes, 2) == CK_INVALID);
      ck_matrix_free(&c);
      ck_matrix_free(&d);
    }
    else
    {
      NEW(c, 5, 2, width, 0);
      ck_argument bad[] = {ck_matrix_argument(&c, 1), ck_matrix_argument(&b, 0)};
      CHECK(ck_array_function_invoke(&fn, bad, 2) == CK_INVALID);
      ck_matrix_free(&c);
    }
    ck_matrix_free(&a);
    ck_matrix_free(&b);
    ck_array_function_free(&fn);
    return 0;
  }
  size_t rows = 5, vec_rows = 3, second_rows = 5;
  int vectors = 1, matrices = 1, should_fail = 1;
  if (IS("VectorTooSmall"))
  {
  }
  else if (IS("VectorTooLarge"))
    vec_rows = 10;
  else if (IS("EmptyVectorNonEmptyMatrix"))
    vec_rows = 0;
  else if (IS("NonEmptyVectorEmptyMatrix"))
  {
    rows = 0;
    vec_rows = 5;
  }
  else if (IS("EmptyVectorEmptyMatrix"))
  {
    rows = vec_rows = 0;
    should_fail = 0;
  }
  else if (IS("MultipleVectorsDifferentSizes"))
  {
    vec_rows = 5;
    vectors = 2;
  }
  else if (IS("MultipleMatricesDifferentRowsVector"))
  {
    rows = vec_rows = 4;
    second_rows = 5;
    matrices = 2;
  }
  else if (IS("VectorSizeMatchesOneMatrixOnly"))
  {
    rows = second_rows = 5;
    vec_rows = 4;
    matrices = 2;
  }
  else if (IS("MismatchedRowsAtInvocation"))
  {
    rows = 5;
    vec_rows = 3;
  }
  else if (IS("FunctionInvocationWithWrongSizedVector"))
  {
    rows = vec_rows = 3;
    should_fail = 0;
  }
  else
    return 2;
  NEW(m, rows, matrices == 2 ? 2 : 3, width, 1);
  NEW(other, second_rows, 3, width, 1);
  double vec[10] = {0};
  ck_view v1, v2;
  CHECK(ck_vector_view(vec, vec_rows, sizeof(double), 0, &v1) == CK_OK);
  CHECK(ck_vector_view(vec, 3, sizeof(double), 0, &v2) == CK_OK);
  ck_argument args[3];
  args[0] = ck_matrix_argument(&m, 1);
  args[1] = matrices == 2 ? ck_matrix_argument(&other, 0) : ck_vector_argument(v1);
  if (matrices == 2)
    args[2] = ck_vector_argument(v1);
  else
    args[2] = ck_vector_argument(v2);
  size_t count = (size_t)(matrices + vectors);
  operation op = {NO_OP, 0};
  ck_array_function fn;
  CHECK(ck_array_function_create(args, count, matrix_callback, &op, &fn) == CK_OK);
  CHECK(ck_array_function_invoke(&fn, args, count) == (should_fail ? CK_INVALID : CK_OK));
  if (IS("MismatchedRowsAtInvocation"))
  {
    ck_matrix_free(&m);
    double zero = 0;
    CHECK(ck_matrix_create(3, 3, width, sizeof(double), &zero, &m) == CK_OK);
    CHECK(ck_vector_view(vec, 5, sizeof(double), 0, &v1) == CK_OK);
    args[0] = ck_matrix_argument(&m, 1);
    args[1] = ck_vector_argument(v1);
    CHECK(ck_array_function_invoke(&fn, args, 2) == CK_INVALID);
  }
  ck_array_function_free(&fn);
  ck_matrix_free(&m);
  ck_matrix_free(&other);
  return 0;
}
static int dense_vectors(const char *test, size_t width)
{
  int const_vector = IS("ConstVector"), mutable = IS("MutableVector"),
      reuse = IS("FunctionReusabilityWithVectors"), integer = IS("IntegerVector"),
      signature = IS("FunctionWithConstSignature"), multiple = IS("MultipleVectorsSameSize"),
      two = IS("MultipleMatricesOneVector"), mixed = IS("MixedVectorColumnViewRowVariable");
  if (!const_vector && !mutable && !reuse && !integer && !signature && !multiple && !two && !mixed)
    return 2;
  size_t rows = multiple ? 5 : two || mixed ? 4 : 3, columns = multiple || mixed ? 3 : 2;
  double zero = 0;
  int izero = 0;
  ck_matrix m;
  CHECK(ck_matrix_create(rows, columns, width, integer ? sizeof(int) : sizeof(double),
                         integer ? (const void *)&izero : (const void *)&zero, &m) == CK_OK);
  NEW(other, rows, 3, width, 0);
  double vec1[5] = {0}, vec2[5] = {0};
  int int_vec[] = {10, 20, 30};
  for (size_t i = 0; i < rows; ++i)
  {
    vec1[i] = const_vector ? (i + 1) * 10.0
              : mutable    ? (i + 1) * 5.0
              : reuse      ? i + 1.0
              : multiple   ? i + 1.0
              : two        ? i * 2.0
              : mixed      ? (i + 1) * 100.0
                           : 0;
    vec2[i] = (i + 1) * 10.0;
    if (!integer)
    {
      D(m, i, 0) = mutable || two || mixed ? i + 1.0 : 0;
      if (mixed)
        D(m, i, 1) = (i + 1) * 10.0;
      D(other, i, 0) = i * 10.0;
    }
  }
  ck_view v1, v2;
  CHECK(ck_vector_view(integer ? (void *)int_vec : (void *)vec1, rows, integer ? sizeof(int) : sizeof(double),
                       mutable, &v1) == CK_OK);
  CHECK(ck_vector_view(vec2, rows, sizeof(double), 0, &v2) == CK_OK);
  operation op = {integer    ? INTEGER_SCALE
                  : mutable  ? VECTOR_FROM_MATRIX
                  : multiple ? VECTOR_SUM
                  : two      ? TWO_MATRIX_VECTOR
                  : mixed    ? MIXED
                             : VECTOR_SCALE,
                  mutable ? 3
                  : reuse ? 10
                          : 2};
  ck_argument args[] = {ck_matrix_argument(&m, !mutable),
                        two ? ck_matrix_argument(&other, 1) : ck_vector_argument(v1),
                        two ? ck_vector_argument(v1) : ck_vector_argument(v2)};
  size_t count = multiple || two ? 3 : 2;
  ck_array_function fn;
  CHECK(ck_array_function_create(args, count, matrix_callback, &op, &fn) == CK_OK);
  CHECK(ck_array_function_invoke(&fn, args, count) == CK_OK);
  for (size_t i = 0; i < rows; ++i)
  {
    if (integer)
      CHECK(I(m, i, 0) == (int)(20 * (i + 1)));
    else if (mutable)
      NEAR(vec1[i], 3 * (i + 1), 0);
    else if (multiple)
      NEAR(D(m, i, 0), 11 * (i + 1), 0);
    else if (two)
    {
      NEAR(D(m, i, 1), i * 3 + 1, 0);
      NEAR(D(other, i, 1), i * 12, 0);
    }
    else if (mixed)
      NEAR(D(m, i, 2), 111 * (i + 1), 0);
    else
      NEAR(D(m, i, 0), vec1[i] * op.factor, 0);
  }
  if (reuse)
  {
    double second[] = {5, 6, 7};
    CHECK(ck_vector_view(second, 3, sizeof(double), 0, &v1) == CK_OK);
    args[1] = ck_vector_argument(v1);
    CHECK(ck_array_function_invoke(&fn, args, 2) == CK_OK);
    for (size_t i = 0; i < 3; ++i)
      NEAR(D(m, i, 0), second[i] * 10, 0);
  }
  ck_array_function_free(&fn);
  ck_matrix_free(&m);
  ck_matrix_free(&other);
  return 0;
}
static ck_status reducer(const void *value, void *accumulator, void *context)
{
  double value_d = *(const double *)value;
  int kind = *(int *)context;
  if (kind == 0)
    *(double *)accumulator += value_d * value_d;
  else if (kind == 1)
  {
    if (value_d > *(double *)accumulator)
      *(double *)accumulator = value_d;
  }
  else if (kind == 2)
    *(int *)accumulator = *(int *)accumulator || value_d > 10;
  else if (kind == 3)
    *(int *)accumulator = *(int *)accumulator && isfinite(value_d);
  else
    *(double *)accumulator += 1;
  return CK_OK;
}
static ck_status captured_vectors(const ck_argument *a, size_t count, void *context)
{
  (void)count;
  const ck_view *vectors = context;
  ck_view col0, col1;
  if (ck_matrix_column(a[0].matrix, 0, 1, &col0) || ck_matrix_column(a[0].matrix, 1, 1, &col1))
    return CK_INVALID;
  for (size_t j = 0; j < 3; ++j)
  {
    size_t index = *(const size_t *)ck_view_const_at(&vectors[0], j);
    double val = (double)*(const size_t *)ck_view_const_at(&vectors[1], index);
    for (size_t i = 0; i < col0.rows; ++i)
      *(double *)ck_view_at(&col0, i) += val;
  }
  size_t index = *(const size_t *)ck_view_const_at(&vectors[0], 1);
  double val = (double)*(const size_t *)ck_view_const_at(&vectors[1], index);
  for (size_t i = 0; i < col1.rows; ++i)
    *(double *)ck_view_at(&col1, i) += val;
  return CK_OK;
}
static int dense_views(const char *test, size_t width)
{
  if (IS("TestFill") || IS("TestCopy"))
  {
    NEW(m, 3, 2, width, 0);
    double vec[] = {3.2, 4.2, 1.3}, out[3] = {0}, temp[3] = {0};
    ck_view c0, c1, v, ov, t;
    CHECK(ck_matrix_column(&m, 0, 1, &c0) == CK_OK);
    CHECK(ck_matrix_column(&m, 1, 1, &c1) == CK_OK);
    CHECK(ck_vector_view(vec, 3, sizeof(double), 1, &v) == CK_OK);
    CHECK(ck_vector_view(out, 3, sizeof(double), 1, &ov) == CK_OK);
    CHECK(ck_vector_view(temp, 3, sizeof(double), 1, &t) == CK_OK);
    if (IS("TestFill"))
    {
      double value = 3.2;
      CHECK(ck_view_fill(c1, &value, 1) == CK_OK);
      CHECK(ck_view_fill(v, &value, 1) == CK_OK);
      value = 9.9;
      CHECK(ck_view_fill(t, &value, 1) == CK_OK);
      CHECK(ck_view_copy(c0, t, 1) == CK_OK);
      for (size_t i = 0; i < 3; ++i)
      {
        NEAR(D(m, i, 0), 9.9, 0);
        NEAR(D(m, i, 1), 3.2, 0);
        NEAR(vec[i], 3.2, 0);
      }
    }
    else
    {
      CHECK(ck_view_copy(c1, v, 1) == CK_OK);
      for (size_t i = 0; i < 3; ++i)
      {
        NEAR(D(m, i, 0), 0, 0);
        NEAR(D(m, i, 1), vec[i], 0);
      }
      CHECK(ck_view_copy(ov, c1, 1) == CK_OK);
      for (size_t i = 0; i < 3; ++i)
        NEAR(out[i], vec[i], 0);
      CHECK(ck_view_copy(c0, c1, 1) == CK_OK);
      for (size_t i = 0; i < 3; ++i)
        NEAR(D(m, i, 0), vec[i], 0);
      double zero = 0;
      CHECK(ck_view_fill(c0, &zero, 1) == CK_OK);
      c1.writable = 0;
      CHECK(ck_view_copy(c0, c1, 1) == CK_OK);
      for (size_t i = 0; i < 3; ++i)
        NEAR(D(m, i, 0), vec[i], 0);
      CHECK(ck_view_copy(t, c1, 1) == CK_OK);
      c1.writable = 1;
      CHECK(ck_view_fill(c1, &zero, 1) == CK_OK);
      CHECK(ck_view_copy(c1, t, 1) == CK_OK);
      for (size_t i = 0; i < 3; ++i)
        NEAR(D(m, i, 1), vec[i], 0);
    }
    ck_matrix_free(&m);
    return 0;
  }
  if (IS("VectorCapture"))
  {
    NEW(m, 3, 3, width, 0);
    size_t first[] = {2, 0, 1}, second[] = {10, 20, 30};
    ck_view vectors[2];
    CHECK(ck_vector_view(first, 3, sizeof(size_t), 0, &vectors[0]) == CK_OK);
    CHECK(ck_vector_view(second, 3, sizeof(size_t), 0, &vectors[1]) == CK_OK);
    ck_argument args[] = {ck_matrix_argument(&m, 1)};
    ck_array_function f;
    CHECK(ck_array_function_create(args, 1, captured_vectors, vectors, &f) == CK_OK);
    CHECK(ck_array_function_invoke(&f, args, 1) == CK_OK);
    for (size_t i = 0; i < 3; ++i)
    {
      NEAR(D(m, i, 0), 60, 0);
      NEAR(D(m, i, 1), 10, 0);
      NEAR(D(m, i, 2), 0, 0);
    }
    ck_array_function_free(&f);
    ck_matrix_free(&m);
    return 0;
  }
  int kind = IS("ReduceSum")      ? 0
             : IS("ReduceMax")    ? 1
             : IS("ReduceLOr")    ? 2
             : IS("ReduceLAnd")   ? 3
             : IS("ReduceStrict") ? 4
                                  : -1;
  if (kind < 0)
    return 2;
  NEW(m, 3, kind == 4 ? 1 : 2, width, kind == 2 || kind == 3 ? 1 : 0);
  if (kind == 0)
  {
    double vals[] = {1, 2, 3, 4, 5, 6};
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 2; ++j)
        D(m, i, j) = vals[i * 2 + j];
  }
  if (kind == 1)
  {
    double vals[] = {1, 8, 5, 2, 3, 4};
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 2; ++j)
        D(m, i, j) = vals[i * 2 + j];
  }
  for (size_t pass = 0; pass < (kind == 2 || kind == 3 ? 2u : 1u); ++pass)
  {
    double acc = kind == 1 ? -DBL_MAX : 0;
    int boolean = kind == 3;
    if (pass)
    {
      if (kind == 2)
        D(m, 2, 1) = 42;
      else
        D(m, 1, 0) = NAN;
    }
    for (size_t j = 0; j < m.columns; ++j)
    {
      ck_view v;
      CHECK(ck_matrix_column(&m, j, 0, &v) == CK_OK);
      CHECK(ck_view_reduce(v, kind == 4, reducer, kind == 2 || kind == 3 ? (void *)&boolean : (void *)&acc,
                           &kind) == CK_OK);
    }
    if (kind == 0)
      NEAR(acc, 91, 0);
    else if (kind == 1)
      NEAR(acc, 8, 0);
    else if (kind == 2)
      CHECK(boolean == (int)pass);
    else if (kind == 3)
      CHECK(boolean == !pass);
    else
      NEAR(acc, 3, 0);
  }
  ck_matrix_free(&m);
  return 0;
}
int test_matrix(const char *suite, const char *test)
{
  int vector = !strcmp(suite, "VectorMatrix");
  if (!vector && strcmp(suite, "Matrix"))
    return 2;
  for (size_t w = vector ? 1 : 0; w <= (vector ? 4u : 0u); ++w)
  {
    size_t width = IS("SmallVectorMatrix") ? 2 : IS("SmallConstVectorMatrix") ? 4 : w;
    int status = dense_basic(test, width);
    if (status == 2)
      status = dense_functions(test, width);
    if (status == 2)
      status = dense_validation(test, width);
    if (status == 2)
      status = dense_vectors(test, width);
    if (status == 2)
      status = dense_views(test, width);
    if (status)
      return status;
  }
  return 0;
}
