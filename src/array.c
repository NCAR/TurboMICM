#include "internal.h"
#include "chem_array.h"
#include <stdarg.h>
#include <stdio.h>

ChemArray *chem_array_dense(size_t rows, size_t columns, size_t group, double initial) {
  if (!group || rows > SIZE_MAX-(group-1)) return NULL;
  size_t padded = ((rows+group-1)/group)*group;
  if (columns && padded > SIZE_MAX/columns/sizeof(double)) return NULL;
  ChemArray *a = calloc(1,sizeof *a);
  if (!a) return NULL;
  a->rows = rows; a->columns = columns; a->group = group; a->allocated = padded*columns;
  a->values = malloc((a->allocated ? a->allocated : 1)*sizeof(double));
  if (!a->values) { free(a); return NULL; }
  chem_array_fill(a,initial); return a;
}
ChemStatus chem_array_from_rows(const double *const *rows, const size_t *lengths, size_t count, size_t group, ChemArray **result) {
  if (!result || (count && (!rows || !lengths))) return CHEM_INVALID;
  size_t columns = count ? lengths[0] : 0;
  for (size_t i = 0; i < count; ++i) if (lengths[i] != columns || (columns && !rows[i])) return CHEM_SIZE;
  ChemArray *a = chem_array_dense(count,columns,group,0.0); if (!a) return CHEM_MEMORY;
  for (size_t i = 0; i < count; ++i) for (size_t j = 0; j < columns; ++j)
    a->values[chem_array_index(a,i,j)] = rows[i][j];
  *result = a; return CHEM_OK;
}
ChemArray *chem_array_sparse(size_t blocks, size_t n, size_t group, int column_order,
                            const size_t *rows, const size_t *cols, size_t elements, double initial) {
  if ((elements && (!rows || !cols)) || (n && n > SIZE_MAX/n/sizeof(int))) return NULL;
  unsigned char *pattern = calloc(n ? n*n : 1,1);
  if (!pattern) return NULL;
  for (size_t i = 0; i < elements; ++i) {
    if (rows[i] >= n || cols[i] >= n) { free(pattern); return NULL; }
    pattern[rows[i]*n+cols[i]] = 1;
  }
  size_t nnz = 0;
  for (size_t i = 0; i < n*n; ++i) nnz += pattern[i];
  ChemArray *a = chem_array_dense(blocks,nnz,group,initial);
  if (!a) { free(pattern); return NULL; }
  a->order = n; a->column_order = !!column_order; a->index = malloc((n ? n*n : 1)*sizeof(int));
  if (!a->index) { free(pattern); chem_array_free(a); return NULL; }
  for (size_t i = 0; i < n*n; ++i) a->index[i] = -1;
  nnz = 0;
  for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) {
    size_t row = column_order ? j : i, col = column_order ? i : j;
    if (pattern[row*n+col]) a->index[row*n+col] = (int)nnz++;
  }
  free(pattern); return a;
}
void chem_array_free(ChemArray *a) { if (a) { free(a->values); free(a->index); free(a); } }
size_t chem_array_index(const ChemArray *a, size_t row, size_t column) {
  if (!a || row >= a->rows || column >= a->columns) return SIZE_MAX;
  return (row/a->group*a->columns+column)*a->group+row%a->group;
}
size_t chem_sparse_index(const ChemArray *a, size_t block, size_t row, size_t column) {
  if (!a || !a->index || row >= a->order || column >= a->order) return SIZE_MAX;
  int index = a->index[row*a->order+column];
  return index < 0 ? SIZE_MAX : chem_array_index(a,block,(size_t)index);
}
ChemStatus chem_array_get(const ChemArray *a, size_t row, size_t column, double *value) {
  size_t i = chem_array_index(a,row,column); if (i == SIZE_MAX || !value) return CHEM_INVALID;
  *value = a->values[i]; return CHEM_OK;
}
ChemStatus chem_array_set(ChemArray *a, size_t row, size_t column, double value) {
  size_t i = chem_array_index(a,row,column); if (i == SIZE_MAX) return CHEM_INVALID;
  a->values[i] = value; return CHEM_OK;
}
ChemStatus chem_array_row(ChemArray *a, size_t row, double *values, size_t count, int write) {
  if (!a || !values || row >= a->rows) return CHEM_INVALID;
  if (count != a->columns) return CHEM_SIZE;
  for (size_t j = 0; j < count; ++j) {
    size_t i = chem_array_index(a,row,j);
    if (write) a->values[i] = values[j]; else values[j] = a->values[i];
  }
  return CHEM_OK;
}
void chem_array_fill(ChemArray *a, double value) {
  if (a) for (size_t i = 0; i < a->allocated; ++i) a->values[i] = value;
}
void chem_array_bound(ChemArray *a, double value, int maximum) {
  if (a) for (size_t i = 0; i < a->allocated; ++i)
    a->values[i] = maximum ? fmax(a->values[i],value) : fmin(a->values[i],value);
}
static int same_structure(const ChemArray *a, const ChemArray *b) {
  if (!a || !b || a->rows != b->rows || a->columns != b->columns || a->order != b->order) return 0;
  if (a->index && b->index) return !memcmp(a->index,b->index,a->order*a->order*sizeof(int));
  return !a->index && !b->index;
}
ChemStatus chem_array_copy(ChemArray *a, const ChemArray *b) {
  if (!same_structure(a,b)) return CHEM_SIZE;
  for (size_t i = 0; i < a->rows; ++i) for (size_t j = 0; j < a->columns; ++j)
    a->values[chem_array_index(a,i,j)] = b->values[chem_array_index(b,i,j)];
  return CHEM_OK;
}
ChemStatus chem_array_axpy(ChemArray *a, double alpha, const ChemArray *b) {
  if (!same_structure(a,b)) return CHEM_SIZE;
  for (size_t i = 0; i < a->rows; ++i) for (size_t j = 0; j < a->columns; ++j)
    a->values[chem_array_index(a,i,j)] += alpha*b->values[chem_array_index(b,i,j)];
  return CHEM_OK;
}
void chem_array_add_diagonal(ChemArray *a, double value) {
  if (!a || !a->index) return;
  for (size_t i = 0; i < a->order; ++i) {
    int col = a->index[i*a->order+i]; if (col < 0) continue;
    for (size_t row = 0; row < a->rows; ++row) a->values[chem_array_index(a,row,(size_t)col)] += value;
  }
}
static void format_append(char *buffer, size_t capacity, size_t *length, const char *format, ...) {
  va_list args; va_start(args,format);
  char *at = *length < capacity ? buffer+*length : NULL;
  size_t room = *length < capacity ? capacity-*length : 0;
  int n = vsnprintf(at,room,format,args);
  if (n > 0) *length += (size_t)n;
  va_end(args);
}
size_t chem_array_format(const ChemArray *a, int nonzero, char *buffer, size_t capacity) {
  size_t length = 0;
  if (!a || (capacity && !buffer)) return 0;
  if (!a->index) {
    for (size_t i = 0; i < a->rows; ++i) {
      for (size_t j = 0; j < a->columns; ++j) format_append(buffer,capacity,&length,j ? ",%g" : "%g",a->values[chem_array_index(a,i,j)]);
      format_append(buffer,capacity,&length,"\n");
    }
  } else for (size_t block = 0; block < a->rows; ++block) {
    format_append(buffer,capacity,&length,"Block %zu\n",block);
    for (size_t i = 0; i < a->order; ++i) {
      for (size_t j = 0; j < a->order; ++j) {
        size_t id = chem_sparse_index(a,block,i,j);
        if (nonzero) {
          if (id != SIZE_MAX) format_append(buffer,capacity,&length,"%zu, %zu, %g\n",i,j,a->values[id]);
        } else format_append(buffer,capacity,&length,j ? ",%g" : "%g",id == SIZE_MAX ? 0.0 : a->values[id]);
      }
      if (!nonzero) format_append(buffer,capacity,&length,"\n");
    }
  }
  return length;
}

ChemView chem_array_view(ChemArray *a, size_t column, int readonly) {
  ChemView v = {a,NULL,column,0,a ? a->rows : 0,readonly}; return v;
}
ChemView chem_vector_view(double *values, size_t count, int readonly) {
  ChemView v = {NULL,values,0,0,count,readonly}; return v;
}
static double *view_pointer(const ChemView *v, size_t i) {
  size_t row = v->start+i;
  if (v->array) {
    ChemArray *a = v->array;
    size_t padded = a->rows ? ((a->rows+a->group-1)/a->group)*a->group : 0;
    if (row >= padded || v->column >= a->columns) return NULL;
    return a->values+(row/a->group*a->columns+v->column)*a->group+row%a->group;
  }
  return v->vector ? v->vector+row : NULL;
}
ChemStatus chem_view_fill(ChemView *v, double value) {
  if (!v || v->readonly) return CHEM_INVALID;
  for (size_t i = 0; i < v->count; ++i) { double *p = view_pointer(v,i); if (!p) return CHEM_INVALID; *p = value; }
  return CHEM_OK;
}
ChemStatus chem_view_copy(ChemView *to, const ChemView *from) {
  if (!to || !from || to->readonly) return CHEM_INVALID;
  if (to->count != from->count) return CHEM_SIZE;
  for (size_t i = 0; i < to->count; ++i) {
    double *p = view_pointer(to,i), *q = view_pointer(from,i); if (!p || !q) return CHEM_INVALID; *p = *q;
  }
  return CHEM_OK;
}
ChemStatus chem_view_for_each(ChemView *views, size_t count, ChemElementFunction fn, void *context) {
  if (!views || !count || !fn) return CHEM_INVALID;
  for (size_t a = 1; a < count; ++a) if (views[a].count != views[0].count) return CHEM_SIZE;
  double *values = malloc(count*sizeof(double)); if (!values) return CHEM_MEMORY;
  for (size_t i = 0; i < views[0].count; ++i) {
    for (size_t a = 0; a < count; ++a) {
      double *p = view_pointer(&views[a],i);
      if (!p) { free(values); return CHEM_INVALID; } values[a] = *p;
    }
    fn(values,count,context);
    for (size_t a = 0; a < count; ++a) if (!views[a].readonly) *view_pointer(&views[a],i) = values[a];
  }
  free(values); return CHEM_OK;
}
ChemStatus chem_view_reduce(const ChemView *view, int operation, int strict, double *acc) {
  if (!view || !acc || operation < 0 || operation > 3) return CHEM_INVALID;
  size_t n = view->count;
  if (!strict && view->array && !view->start && view->count == view->array->rows)
    n = ((n+view->array->group-1)/view->array->group)*view->array->group;
  for (size_t i = 0; i < n; ++i) {
    double *p = view_pointer(view,i); if (!p) return CHEM_INVALID;
    if (operation == 0) *acc += *p;
    else if (operation == 1) *acc = fmax(*acc,*p);
    else if (operation == 2) *acc = *acc != 0.0 || *p != 0.0;
    else *acc = *acc != 0.0 && *p != 0.0;
  }
  return CHEM_OK;
}
ChemStatus chem_view_reduce_function(const ChemView *view, int strict, double *acc, ChemReduceFunction fn, void *context) {
  if (!view || !acc || !fn) return CHEM_INVALID;
  size_t n = view->count;
  if (!strict && view->array && !view->start && view->count == view->array->rows)
    n = ((n+view->array->group-1)/view->array->group)*view->array->group;
  for (size_t i = 0; i < n; ++i) {
    double *p = view_pointer(view,i); if (!p) return CHEM_INVALID;
    fn(*p,acc,context);
  }
  return CHEM_OK;
}

struct ChemFunction {
  size_t count;
  ChemArgument *signature;
  ChemArray *shapes;
  ChemArrayFunction function;
  void *context;
};
static ChemStatus arguments_valid(const ChemArgument *args, size_t count, int check_rows) {
  if (!args || !count) return CHEM_INVALID;
  size_t rows = args[0].array ? args[0].array->rows : args[0].count, group = 0;
  for (size_t i = 0; i < count; ++i) {
    size_t n = args[i].array ? args[i].array->rows : args[i].count;
    if (check_rows && rows != n) return CHEM_SIZE;
    if (args[i].array) {
      if (group && group != args[i].array->group) return CHEM_INVALID;
      group = args[i].array->group;
    } else if (n && !args[i].vector) return CHEM_INVALID;
  }
  return CHEM_OK;
}
ChemStatus chem_function_create(const ChemArgument *args, size_t count, ChemArrayFunction fn, void *context, ChemFunction **result) {
  if (!result || !fn) return CHEM_INVALID;
  ChemStatus status = arguments_valid(args,count,0); if (status) return status;
  ChemFunction *f = calloc(1,sizeof *f); if (!f) return CHEM_MEMORY;
  f->signature = malloc(count*sizeof(ChemArgument)); f->shapes = calloc(count,sizeof(ChemArray));
  if (!f->signature || !f->shapes) { chem_function_free(f); return CHEM_MEMORY; }
  f->count = count; f->function = fn; f->context = context;
  for (size_t i = 0; i < count; ++i) {
    f->signature[i] = args[i];
    if (args[i].array) { f->shapes[i] = *args[i].array; f->signature[i].array = &f->shapes[i]; }
  }
  *result = f; return CHEM_OK;
}
ChemStatus chem_function_invoke(const ChemFunction *f, ChemArgument *args, size_t count) {
  if (!f || count != f->count) return CHEM_SIZE;
  ChemStatus status = arguments_valid(args,count,1); if (status) return status;
  for (size_t i = 0; i < count; ++i) {
    const ChemArgument *s = &f->signature[i];
    if (!!args[i].array != !!s->array || (!s->readonly && args[i].readonly)) return CHEM_INVALID;
    if (s->array) {
      ChemArray *a = args[i].array, *b = s->array;
      if (a->columns != b->columns || a->order != b->order || a->group != b->group || a->column_order != b->column_order) return CHEM_SIZE;
    }
  }
  f->function(args,count,f->context); return CHEM_OK;
}
void chem_function_free(ChemFunction *f) { if (f) { free(f->signature); free(f->shapes); free(f); } }
