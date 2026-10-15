#include "chem_matrix.h"
#include "internal.h"

ck_status ck_matrix_create(size_t rows, size_t columns, size_t width, size_t element_size,
                           const void *initial, ck_matrix *m)
{
  if (!m || !element_size)
    return CK_INVALID;
  memset(m, 0, sizeof *m);
  if (width && rows > SIZE_MAX - (width - 1))
    return CK_INVALID;
  size_t padded = width ? ((rows + width - 1) / width) * width : rows;
  if (columns && padded > SIZE_MAX / columns)
    return CK_INVALID;
  size_t count = padded * columns;
  if (count > SIZE_MAX / element_size)
    return CK_INVALID;
  m->data = ck_allocate(count, element_size);
  if (!m->data)
    return CK_MEMORY;
  m->rows = rows;
  m->columns = columns;
  m->width = width;
  m->padded_rows = padded;
  m->size = count;
  m->element_size = element_size;
  if (initial)
    ck_matrix_fill(m, initial);
  return CK_OK;
}
void ck_matrix_free(ck_matrix *m)
{
  if (m)
  {
    free(m->data);
    memset(m, 0, sizeof *m);
  }
}
ck_status ck_matrix_clone(const ck_matrix *s, ck_matrix *d)
{
  if (!s || !d || s == d)
    return CK_INVALID;
  ck_status status = ck_matrix_create(s->rows, s->columns, s->width, s->element_size, NULL, d);
  if (!status)
    memcpy(d->data, s->data, s->size * s->element_size);
  return status;
}
size_t ck_matrix_index(const ck_matrix *m, size_t row, size_t col)
{
  if (!m || row >= m->rows || col >= m->columns)
    return SIZE_MAX;
  return m->width ? (row / m->width) * (m->width * m->columns) + col * m->width + row % m->width
                  : row * m->columns + col;
}
const void *ck_matrix_const_at(const ck_matrix *m, size_t row, size_t col)
{
  size_t index = ck_matrix_index(m, row, col);
  return index == SIZE_MAX ? NULL : (const unsigned char *)m->data + index * m->element_size;
}
void *ck_matrix_at(ck_matrix *m, size_t row, size_t col)
{
  return (void *)ck_matrix_const_at(m, row, col);
}
ck_status ck_matrix_set_row(ck_matrix *m, size_t row, const void *v, size_t count)
{
  if (!m || row >= m->rows || count < m->columns || (!v && count))
    return CK_INVALID;
  for (size_t j = 0; j < m->columns; ++j)
    memcpy(ck_matrix_at(m, row, j), (const unsigned char *)v + j * m->element_size, m->element_size);
  return CK_OK;
}
ck_status ck_matrix_get_row(const ck_matrix *m, size_t row, void *v, size_t count)
{
  if (!m || row >= m->rows || count < m->columns || (!v && count))
    return CK_INVALID;
  for (size_t j = 0; j < m->columns; ++j)
    memcpy((unsigned char *)v + j * m->element_size, ck_matrix_const_at(m, row, j), m->element_size);
  return CK_OK;
}
ck_status ck_matrix_from_rows(size_t rows, const void *const *v, const size_t *lengths, size_t width,
                              size_t element_size, ck_matrix *m)
{
  if (rows && (!v || !lengths))
    return CK_INVALID;
  size_t cols = rows ? lengths[0] : 0;
  for (size_t i = 0; i < rows; ++i)
    if (lengths[i] != cols || (!v[i] && cols))
      return CK_INVALID;
  ck_status status = ck_matrix_create(rows, cols, width, element_size, NULL, m);
  if (status)
    return status;
  for (size_t i = 0; i < rows; ++i)
    ck_matrix_set_row(m, i, v[i], cols);
  return CK_OK;
}
void ck_matrix_fill(ck_matrix *m, const void *v)
{
  if (m && v)
    for (size_t i = 0; i < m->size; ++i)
      memcpy((unsigned char *)m->data + i * m->element_size, v, m->element_size);
}
ck_status ck_matrix_axpy(ck_matrix *y, double alpha, const ck_matrix *x)
{
  if (!x || !y || x->rows != y->rows || x->columns != y->columns || x->element_size != sizeof(double) ||
      y->element_size != sizeof(double))
    return CK_INVALID;
  for (size_t j = 0; j < y->columns; ++j)
    for (size_t i = 0; i < y->rows; ++i)
      *(double *)ck_matrix_at(y, i, j) += alpha * *(const double *)ck_matrix_const_at(x, i, j);
  return CK_OK;
}
void ck_matrix_bound(ck_matrix *m, double value, int lower)
{
  if (!m || m->element_size != sizeof(double))
    return;
  double *data = m->data;
  for (size_t i = 0; i < m->size; ++i)
    data[i] = lower ? (data[i] > value ? data[i] : value) : (data[i] < value ? data[i] : value);
}
ck_status ck_matrix_print(const ck_matrix *m, FILE *out)
{
  if (!m || !out || m->element_size != sizeof(double))
    return CK_INVALID;
  for (size_t i = 0; i < m->rows; ++i)
  {
    for (size_t j = 0; j < m->columns; ++j)
      fprintf(out, "%s%g", j ? "," : "", *(const double *)ck_matrix_const_at(m, i, j));
    fputc('\n', out);
  }
  return ferror(out) ? CK_INVALID : CK_OK;
}
ck_status ck_sparse_create(size_t n, size_t blocks, size_t width, ck_sparse_order ordering,
                           const unsigned char *pattern, size_t element_size, const void *initial,
                           ck_sparse *m)
{
  if (!m || !n || n > SIZE_MAX / n || !pattern || (unsigned)ordering > CK_SPARSE_COLUMN)
    return CK_INVALID;
  memset(m, 0, sizeof *m);
  m->order = n;
  m->ordering = ordering;
  for (size_t i = 0; i < n * n; ++i)
    if (pattern[i])
      ++m->nnz;
  m->map = ck_allocate(n * n, sizeof(size_t));
  m->rows = ck_allocate(m->nnz, sizeof(size_t));
  m->columns = ck_allocate(m->nnz, sizeof(size_t));
  if (!m->map || !m->rows || !m->columns)
  {
    ck_sparse_free(m);
    return CK_MEMORY;
  }
  for (size_t i = 0; i < n * n; ++i)
    m->map[i] = SIZE_MAX;
  size_t next = 0;
  for (size_t major = 0; major < n; ++major)
    for (size_t minor = 0; minor < n; ++minor)
    {
      size_t row = ordering == CK_SPARSE_ROW ? major : minor, col = ordering == CK_SPARSE_ROW ? minor : major;
      if (pattern[row * n + col])
      {
        m->map[row * n + col] = next;
        m->rows[next] = row;
        m->columns[next++] = col;
      }
    }
  ck_status status = ck_matrix_create(blocks, m->nnz, width, element_size, initial, &m->values);
  if (status)
    ck_sparse_free(m);
  return status;
}
ck_status ck_sparse_from_entries(size_t n, size_t blocks, size_t width, ck_sparse_order ordering,
                                 const size_t *rows, const size_t *columns, size_t count, size_t element_size,
                                 const void *initial, ck_sparse *m)
{
  if (!n || n > SIZE_MAX / n || (count && (!rows || !columns)))
    return CK_INVALID;
  for (size_t j = 0; j < count; ++j)
    if (rows[j] >= n || columns[j] >= n)
      return CK_INVALID;
  unsigned char *pattern = ck_allocate(n * n, 1);
  if (!pattern)
    return CK_MEMORY;
  for (size_t j = 0; j < count; ++j)
    pattern[rows[j] * n + columns[j]] = 1;
  ck_status status = ck_sparse_create(n, blocks, width, ordering, pattern, element_size, initial, m);
  free(pattern);
  return status;
}
void ck_sparse_free(ck_sparse *m)
{
  if (m)
  {
    free(m->map);
    free(m->rows);
    free(m->columns);
    ck_matrix_free(&m->values);
    memset(m, 0, sizeof *m);
  }
}
size_t ck_sparse_index(const ck_sparse *m, size_t block, size_t row, size_t col)
{
  if (!m || row >= m->order || col >= m->order)
    return SIZE_MAX;
  size_t index = m->map[row * m->order + col];
  return index == SIZE_MAX ? SIZE_MAX : ck_matrix_index(&m->values, block, index);
}
const void *ck_sparse_const_at(const ck_sparse *m, size_t block, size_t row, size_t col)
{
  size_t index = ck_sparse_index(m, block, row, col);
  return index == SIZE_MAX ? NULL : (const unsigned char *)m->values.data + index * m->values.element_size;
}
void *ck_sparse_at(ck_sparse *m, size_t block, size_t row, size_t col)
{
  return (void *)ck_sparse_const_at(m, block, row, col);
}
ck_status ck_sparse_diagonal(ck_sparse *m, double alpha)
{
  if (!m || m->values.element_size != sizeof(double))
    return CK_INVALID;
  for (size_t i = 0; i < m->order; ++i)
    if (m->map[i * m->order + i] != SIZE_MAX)
      for (size_t c = 0; c < m->values.rows; ++c)
        *(double *)ck_sparse_at(m, c, i, i) += alpha;
  return CK_OK;
}
ck_status ck_sparse_print(const ck_sparse *m, FILE *out, int nonzero)
{
  if (!m || !out || m->values.element_size != sizeof(double))
    return CK_INVALID;
  for (size_t block = 0; block < m->values.rows; ++block)
  {
    fprintf(out, "Block %zu\n", block);
    for (size_t i = 0; i < m->order; ++i)
    {
      for (size_t j = 0; j < m->order; ++j)
      {
        const double *v = ck_sparse_const_at(m, block, i, j);
        if (nonzero)
        {
          if (v)
            fprintf(out, "%zu, %zu, %g\n", i, j, *v);
        }
        else
        {
          if (j)
            fputc(',', out);
          if (v)
            fprintf(out, "%g", *v);
          else
            fputc('0', out);
        }
      }
      if (!nonzero)
        fputc('\n', out);
    }
  }
  return ferror(out) ? CK_INVALID : CK_OK;
}
ck_status ck_matrix_column(const ck_matrix *m, size_t col, int writable, ck_view *v)
{
  if (!m || !v || col >= m->columns)
    return CK_INVALID;
  *v = (ck_view){.data = m->data,
                 .rows = m->rows,
                 .capacity = m->padded_rows,
                 .element_size = m->element_size,
                 .writable = writable,
                 .width = m->width,
                 .group_stride = m->width * m->columns,
                 .lane_stride = m->width ? 1 : m->columns,
                 .offset = col * (m->width ? m->width : 1)};
  return CK_OK;
}
ck_status ck_sparse_block_view(const ck_sparse *m, size_t index, int writable, ck_view *v)
{
  return m ? ck_matrix_column(&m->values, index, writable, v) : CK_INVALID;
}
ck_status ck_vector_view(void *data, size_t rows, size_t element_size, int writable, ck_view *v)
{
  if (!v || !element_size || (!data && rows))
    return CK_INVALID;
  *v = (ck_view){.data = data,
                 .rows = rows,
                 .capacity = rows,
                 .element_size = element_size,
                 .lane_stride = 1,
                 .writable = writable};
  return CK_OK;
}
const void *ck_view_const_at(const ck_view *v, size_t row)
{
  if (!v || row >= v->capacity)
    return NULL;
  size_t index =
      v->offset + (v->width ? (row / v->width) * v->group_stride + row % v->width : row * v->lane_stride);
  return (const unsigned char *)v->data + index * v->element_size;
}
void *ck_view_at(ck_view *v, size_t row)
{
  return v && v->writable ? (void *)ck_view_const_at(v, row) : NULL;
}
ck_status ck_view_fill(ck_view d, const void *value, int strict)
{
  if (!d.writable || !value)
    return CK_INVALID;
  for (size_t i = 0; i < (strict ? d.rows : d.capacity); ++i)
    memcpy(ck_view_at(&d, i), value, d.element_size);
  return CK_OK;
}
ck_status ck_view_copy(ck_view d, ck_view s, int strict)
{
  if (!d.writable || d.element_size != s.element_size || d.rows != s.rows)
    return CK_INVALID;
  size_t count = strict ? d.rows : d.capacity;
  if (count > s.capacity)
    return CK_INVALID;
  for (size_t i = 0; i < count; ++i)
    memmove(ck_view_at(&d, i), ck_view_const_at(&s, i), d.element_size);
  return CK_OK;
}
ck_status ck_view_apply(ck_view *v, size_t count, int strict, ck_view_callback fn, void *ctx)
{
  if (!v || !count || !fn)
    return CK_INVALID;
  size_t length = strict ? v[0].rows : v[0].capacity;
  for (size_t j = 1; j < count; ++j)
  {
    if (v[j].rows != v[0].rows)
      return CK_INVALID;
    if (!strict && v[j].capacity < length)
      length = v[j].capacity;
  }
  for (size_t i = 0; i < length; ++i)
  {
    ck_status status = fn(v, count, i, ctx);
    if (status)
      return status;
  }
  return CK_OK;
}
ck_status ck_view_reduce(ck_view v, int strict, ck_reduce_callback fn, void *acc, void *ctx)
{
  if (!fn || !acc)
    return CK_INVALID;
  for (size_t i = 0; i < (strict ? v.rows : v.capacity); ++i)
  {
    ck_status status = fn(ck_view_const_at(&v, i), acc, ctx);
    if (status)
      return status;
  }
  return CK_OK;
}
ck_argument ck_matrix_argument(const ck_matrix *m, int writable)
{
  return (ck_argument){.kind = CK_ARGUMENT_MATRIX, .matrix = m, .writable = writable};
}
ck_argument ck_sparse_argument(const ck_sparse *m, int writable)
{
  return (ck_argument){.kind = CK_ARGUMENT_SPARSE, .sparse = m, .writable = writable};
}
ck_argument ck_vector_argument(ck_view v)
{
  return (ck_argument){.kind = CK_ARGUMENT_VECTOR, .vector = v, .writable = v.writable};
}
static size_t argument_rows(const ck_argument *a)
{
  return a->kind == CK_ARGUMENT_VECTOR   ? a->vector.rows
         : a->kind == CK_ARGUMENT_MATRIX ? a->matrix->rows
                                         : a->sparse->values.rows;
}
static const ck_matrix *argument_matrix(const ck_argument *a)
{
  return a->kind == CK_ARGUMENT_MATRIX ? a->matrix : &a->sparse->values;
}
ck_status ck_array_function_create(const ck_argument *prototypes, size_t count, ck_array_callback callback,
                                   void *context, ck_array_function *f)
{
  if (!f || !prototypes || !count || !callback)
    return CK_INVALID;
  memset(f, 0, sizeof *f);
  size_t width = 0;
  for (size_t j = 0; j < count; ++j)
    if (prototypes[j].kind != CK_ARGUMENT_VECTOR)
    {
      const ck_matrix *m = argument_matrix(&prototypes[j]);
      if (!m)
        return CK_INVALID;
      size_t logical_width = m->width ? m->width : 1;
      if (width && logical_width != width)
        return CK_INVALID;
      width = logical_width;
    }
  f->prototypes = ck_allocate(count, sizeof(ck_argument));
  if (!f->prototypes)
    return CK_MEMORY;
  /* Only metadata is captured; callers may destroy the original objects. */
  for (size_t j = 0; j < count; ++j)
  {
    f->prototypes[j] = prototypes[j];
    if (prototypes[j].kind == CK_ARGUMENT_VECTOR)
      continue;
    const ck_matrix *m = argument_matrix(&prototypes[j]);
    ck_matrix *copy = ck_allocate(1, sizeof *copy);
    if (!copy)
    {
      f->count = j;
      ck_array_function_free(f);
      return CK_MEMORY;
    }
    *copy = *m;
    copy->data = NULL;
    if (prototypes[j].kind == CK_ARGUMENT_MATRIX)
      f->prototypes[j].matrix = copy;
    else
    {
      ck_sparse *s = ck_allocate(1, sizeof *s);
      if (!s)
      {
        free(copy);
        f->count = j;
        ck_array_function_free(f);
        return CK_MEMORY;
      }
      *s = *prototypes[j].sparse;
      s->values = *copy;
      s->map = NULL;
      s->rows = NULL;
      s->columns = NULL;
      free(copy);
      f->prototypes[j].sparse = s;
    }
  }
  f->count = count;
  f->callback = callback;
  f->context = context;
  return CK_OK;
}
ck_status ck_array_function_invoke(const ck_array_function *f, const ck_argument *a, size_t count)
{
  if (!f || !a || count != f->count)
    return CK_INVALID;
  size_t rows = argument_rows(&a[0]), width = SIZE_MAX;
  for (size_t j = 0; j < count; ++j)
  {
    const ck_argument *p = &f->prototypes[j];
    if (a[j].kind != p->kind || argument_rows(&a[j]) != rows || (p->writable && !a[j].writable))
      return CK_INVALID;
    if (a[j].kind == CK_ARGUMENT_VECTOR)
    {
      if (a[j].vector.element_size != p->vector.element_size)
        return CK_INVALID;
      continue;
    }
    const ck_matrix *m = argument_matrix(&a[j]), *proto = argument_matrix(p);
    if (m->columns != proto->columns || m->width != proto->width || m->element_size != proto->element_size)
      return CK_INVALID;
    size_t logical_width = m->width ? m->width : 1;
    if (width == SIZE_MAX)
      width = logical_width;
    else if (width != logical_width)
      return CK_INVALID;
    if (a[j].kind == CK_ARGUMENT_SPARSE &&
        (a[j].sparse->order != p->sparse->order || a[j].sparse->ordering != p->sparse->ordering))
      return CK_INVALID;
  }
  return f->callback(a, count, f->context);
}
void ck_array_function_free(ck_array_function *f)
{
  if (!f)
    return;
  for (size_t j = 0; j < f->count; ++j)
  {
    if (f->prototypes[j].kind == CK_ARGUMENT_MATRIX)
      free((void *)f->prototypes[j].matrix);
    else if (f->prototypes[j].kind == CK_ARGUMENT_SPARSE)
      free((void *)f->prototypes[j].sparse);
  }
  free(f->prototypes);
  memset(f, 0, sizeof *f);
}
