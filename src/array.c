/* SPDX-License-Identifier: Apache-2.0 */
#include "chem_array.h"
#include "internal.h"
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct ch_array {
  size_t records, elements, dimension, group, storage;
  int type, kind, order;
  int *indices;
  void *data;
};
typedef struct {
  size_t elements, dimension, group;
  int type, kind, order;
  int *indices;
  unsigned char read_only;
} ch_array_signature;
struct ch_array_function {
  size_t count, group;
  ch_array_signature *signatures;
  ch_group_function function;
  void *context;
};

static size_t item_size(int type) { return type == CH_ARRAY_INT ? sizeof(int) : sizeof(double); }
int ch_pattern_set(size_t n, unsigned char *pattern, size_t row, size_t col)
{
  if (!pattern || (n && n > SIZE_MAX / n)) return CH_INVALID;
  if (row >= n || col >= n) return CH_SIZE;
  pattern[row*n+col] = 1; return CH_OK;
}
static double read_item(const void *p, size_t i, int type)
{ return type == CH_ARRAY_INT ? ((const int *)p)[i] : ((const double *)p)[i]; }
static int write_item(void *p, size_t i, int type, double v)
{
  if (type == CH_ARRAY_INT) {
    if (!isfinite(v) || v > INT_MAX || v < INT_MIN) return CH_INVALID;
    ((int *)p)[i] = (int)v;
  } else ((double *)p)[i] = v;
  return CH_OK;
}
static size_t position(const ch_array *a, size_t r, size_t e)
{
  return a->group ? (r / a->group) * a->group * a->elements + e * a->group + r % a->group
                  : r * a->elements + e;
}
static ch_array *allocate(size_t records, size_t elements, size_t group, int type, double initial)
{
  if (type != CH_ARRAY_REAL && type != CH_ARRAY_INT) return NULL;
  size_t padded = records;
  if (group && records % group) {
    if (records > SIZE_MAX - (group - records % group)) return NULL;
    padded += group - records % group;
  }
  if (elements && padded > SIZE_MAX / elements) return NULL;
  ch_array *a = calloc(1, sizeof(*a));
  if (!a) return NULL;
  a->records = records; a->elements = elements; a->group = group; a->type = type;
  a->storage = padded * elements;
  a->data = ch_allocate(a->storage, item_size(type));
  if (!a->data || ch_array_fill(a, initial) != CH_OK) { ch_array_free(a); return NULL; }
  return a;
}
ch_array *ch_array_dense(size_t rows, size_t cols, size_t group, int type, double initial)
{ return allocate(rows, cols, group, type, initial); }
ch_array *ch_array_vector(size_t size, int type, double initial)
{
  ch_array *a = allocate(size, 1, 0, type, initial);
  if (a) a->kind = CH_ARRAY_VECTOR;
  return a;
}
ch_array *ch_array_sparse(size_t blocks, size_t n, const unsigned char *pattern,
                          int order, size_t group, int type, double initial)
{
  if ((n && !pattern) || (n && n > SIZE_MAX / n) || (order != 0 && order != 1)) return NULL;
  size_t count = 0;
  for (size_t i = 0; i < n * n; ++i) count += pattern[i] != 0;
  if (count > INT_MAX) return NULL;
  ch_array *a = allocate(blocks, count, group, type, initial);
  if (!a) return NULL;
  a->kind = CH_ARRAY_SPARSE; a->order = order; a->dimension = n;
  a->indices = ch_allocate(n * n, sizeof(int));
  if (!a->indices) { ch_array_free(a); return NULL; }
  for (size_t i = 0; i < n * n; ++i) a->indices[i] = -1;
  int k = 0;
  for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) {
    size_t r = order ? j : i, c = order ? i : j;
    if (pattern[r * n + c]) a->indices[r * n + c] = k++;
  }
  return a;
}
void ch_array_free(ch_array *a) { if (a) { free(a->indices); free(a->data); free(a); } }
ch_array *ch_array_clone(const ch_array *src)
{
  if (!src) return NULL;
  ch_array *a = allocate(src->records, src->elements, src->group, src->type, 0);
  if (!a) return NULL;
  a->kind = src->kind; a->order = src->order; a->dimension = src->dimension;
  if (src->indices) {
    size_t bytes = src->dimension * src->dimension * sizeof(int);
    a->indices = malloc(bytes ? bytes : 1);
    if (!a->indices) { ch_array_free(a); return NULL; }
    memcpy(a->indices, src->indices, bytes);
  }
  memcpy(a->data, src->data, src->storage * item_size(src->type));
  return a;
}
size_t ch_array_records(const ch_array *a) { return a ? a->records : 0; }
size_t ch_array_elements(const ch_array *a) { return a ? a->elements : 0; }
size_t ch_array_storage_size(const ch_array *a) { return a ? a->storage : 0; }
size_t ch_array_group_size(const ch_array *a) { return a ? a->group : 0; }
size_t ch_array_groups(const ch_array *a)
{ return !a || !a->records ? 0 : a->group ? 1 + (a->records - 1) / a->group : a->records; }
size_t ch_array_row_stride(const ch_array *a) { return !a ? 0 : a->group ? 1 : a->elements; }
size_t ch_array_column_stride(const ch_array *a) { return !a ? 0 : a->group ? a->group : 1; }
void *ch_array_data(ch_array *a) { return a ? a->data : NULL; }
int ch_array_element(const ch_array *a, size_t r, size_t c)
{
  if (!a || a->kind != CH_ARRAY_SPARSE) return CH_INVALID;
  if (r >= a->dimension || c >= a->dimension) return CH_SIZE;
  int i = a->indices[r * a->dimension + c];
  return i < 0 ? CH_NOT_FOUND : i;
}
int ch_array_index(const ch_array *a, size_t r, size_t e, size_t *i)
{
  if (!a || !i) return CH_INVALID;
  if (r >= a->records || e >= a->elements) return CH_SIZE;
  *i = position(a, r, e); return CH_OK;
}
int ch_array_get(const ch_array *a, size_t r, size_t e, double *v)
{
  size_t i; int status = ch_array_index(a, r, e, &i);
  if (status || !v) return status ? status : CH_INVALID;
  *v = read_item(a->data, i, a->type); return CH_OK;
}
int ch_array_set(ch_array *a, size_t r, size_t e, double v)
{
  size_t i; int status = ch_array_index(a, r, e, &i);
  return status ? status : write_item(a->data, i, a->type, v);
}
int ch_array_set_row(ch_array *a, size_t r, const double *v, size_t count)
{
  if (!a || (!v && count)) return CH_INVALID;
  if (r >= a->records || count < a->elements) return CH_SIZE;
  for (size_t e = 0; e < a->elements; ++e) {
    int status = ch_array_set(a, r, e, v[e]); if (status) return status;
  }
  return CH_OK;
}
int ch_array_fill(ch_array *a, double value)
{
  if (!a) return CH_INVALID;
  for (size_t i = 0; i < a->storage; ++i) { int s = write_item(a->data, i, a->type, value); if (s) return s; }
  return CH_OK;
}
static int compatible(const ch_array *a, const ch_array *b)
{
  if (!a || !b) return CH_INVALID;
  if (a->records != b->records || a->elements != b->elements || a->kind != b->kind ||
      a->dimension != b->dimension) return CH_SIZE;
  if (a->kind == CH_ARRAY_SPARSE && memcmp(a->indices, b->indices, a->dimension * a->dimension * sizeof(int)))
    return CH_INVALID;
  return CH_OK;
}
int ch_array_copy(ch_array *a, const ch_array *b)
{
  int s = compatible(a, b); if (s) return s;
  if (a->type == b->type && a->group == b->group) memcpy(a->data, b->data, a->storage * item_size(a->type));
  else for (size_t r = 0; r < a->records; ++r) for (size_t e = 0; e < a->elements; ++e) {
    s = ch_array_set(a, r, e, read_item(b->data, position(b,r,e), b->type)); if (s) return s;
  }
  return CH_OK;
}
int ch_array_axpy(ch_array *a, double alpha, const ch_array *b)
{
  int s = compatible(a, b); if (s) return s;
  for (size_t r = 0; r < a->records; ++r) for (size_t e = 0; e < a->elements; ++e) {
    size_t i = position(a,r,e);
    s = write_item(a->data, i, a->type, read_item(a->data,i,a->type) + alpha * read_item(b->data,position(b,r,e),b->type));
    if (s) return s;
  }
  return CH_OK;
}
static int bound(ch_array *a, double v, int maximum)
{
  if (!a) return CH_INVALID;
  for (size_t i = 0; i < a->storage; ++i) {
    double x = read_item(a->data,i,a->type);
    int s = write_item(a->data,i,a->type, maximum ? (x < v ? v : x) : (v < x ? v : x));
    if (s) return s;
  }
  return CH_OK;
}
int ch_array_max(ch_array *a, double v) { return bound(a,v,1); }
int ch_array_min(ch_array *a, double v) { return bound(a,v,0); }
int ch_array_add_diagonal(ch_array *a, double v)
{
  if (!a || a->kind != CH_ARRAY_SPARSE) return CH_INVALID;
  for (size_t i = 0; i < a->dimension; ++i) {
    int e = a->indices[i*a->dimension+i]; if (e < 0) continue;
    for (size_t r = 0; r < a->records; ++r) {
      size_t j = position(a,r,(size_t)e);
      int s = write_item(a->data,j,a->type,read_item(a->data,j,a->type)+v); if (s) return s;
    }
  }
  return CH_OK;
}
static void append(char *out, size_t cap, size_t *used, const char *text)
{
  size_t n = strlen(text);
  if (out && *used < cap) {
    size_t copy = n < cap - *used - 1 ? n : cap - *used - 1;
    memcpy(out + *used,text,copy); out[*used+copy] = '\0';
  }
  *used += n;
}
size_t ch_array_format(const ch_array *a, char *out, size_t cap)
{
  size_t used = 0; char buf[160];
  if (out && cap) *out = '\0';
  if (!a) return 0;
  for (size_t r = 0; r < a->records; ++r) {
    if (a->kind == CH_ARRAY_SPARSE) {
      snprintf(buf,sizeof(buf),"Block %zu\n",r); append(out,cap,&used,buf);
      for (size_t i = 0; i < a->dimension; ++i) for (size_t j = 0; j < a->dimension; ++j) {
        int e = a->indices[i*a->dimension+j]; if (e < 0) continue;
        snprintf(buf,sizeof(buf),"%zu, %zu, %.6g\n",i,j,read_item(a->data,position(a,r,(size_t)e),a->type));
        append(out,cap,&used,buf);
      }
    } else {
      for (size_t e = 0; e < a->elements; ++e) {
        snprintf(buf,sizeof(buf),"%s%.6g",e ? "," : "",read_item(a->data,position(a,r,e),a->type));
        append(out,cap,&used,buf);
      }
      append(out,cap,&used,"\n");
    }
  }
  return used;
}
size_t ch_array_format_dense(const ch_array *a, char *out, size_t cap)
{
  if (!a || a->kind != CH_ARRAY_SPARSE) return ch_array_format(a,out,cap);
  size_t used = 0; char buf[160];
  if (out && cap) *out = '\0';
  for (size_t r = 0; r < a->records; ++r) {
    snprintf(buf,sizeof(buf),"Block %zu\n",r); append(out,cap,&used,buf);
    for (size_t i = 0; i < a->dimension; ++i) {
      for (size_t j = 0; j < a->dimension; ++j) {
        int e = a->indices[i*a->dimension+j];
        double value = e < 0 ? 0 : read_item(a->data,position(a,r,(size_t)e),a->type);
        snprintf(buf,sizeof(buf),"%s%.6g",j ? "," : "",value); append(out,cap,&used,buf);
      }
      append(out,cap,&used,"\n");
    }
  }
  return used;
}
int ch_array_print(const ch_array *a, FILE *stream)
{
  if (!a || !stream) return CH_INVALID;
  size_t n = ch_array_format(a,NULL,0); char *buf = malloc(n+1);
  if (!buf) return CH_MEMORY;
  ch_array_format(a,buf,n+1); int status = fputs(buf,stream) < 0 ? CH_INVALID : CH_OK;
  free(buf); return status;
}
int ch_array_group(ch_array *a, size_t i, int ro, ch_group *g)
{
  if (!a || !g) return CH_INVALID;
  if (i >= ch_array_groups(a)) return CH_SIZE;
  size_t start = a->group ? i*a->group : i;
  size_t capacity = a->group ? a->group : 1;
  *g = (ch_group){ a, start, a->records-start < capacity ? a->records-start : capacity, capacity, ro != 0 };
  return CH_OK;
}
int ch_group_view(const ch_group *g, size_t e, ch_view *v)
{
  if (!g || !g->array || !v) return CH_INVALID;
  ch_array *a = g->array;
  if (e >= a->elements) return CH_SIZE;
  size_t stride = a->group || a->kind == CH_ARRAY_VECTOR ? 1 : a->elements;
  *v = (ch_view){ (unsigned char *)a->data+position(a,g->start,e)*item_size(a->type),
    g->size, a->kind == CH_ARRAY_VECTOR ? g->size : g->capacity, stride, a->type, g->read_only,
    a->kind == CH_ARRAY_SPARSE ? CH_VIEW_BLOCK : a->kind == CH_ARRAY_VECTOR ? CH_VIEW_VECTOR : CH_VIEW_COLUMN };
  return CH_OK;
}
int ch_group_sparse_view(const ch_group *g, size_t r, size_t c, ch_view *v)
{
  if (!g) return CH_INVALID;
  int e = ch_array_element(g->array,r,c);
  return e < 0 ? e : ch_group_view(g,(size_t)e,v);
}
ch_view ch_temporary_view(double *p, size_t n, size_t cap)
{ return (ch_view){p,n,cap,1,CH_ARRAY_REAL,0,CH_VIEW_TEMPORARY}; }
double ch_view_get(const ch_view *v, size_t i)
{ return !v || i >= v->capacity || !v->data ? NAN : read_item(v->data,i*v->stride,v->type); }
int ch_view_set(ch_view *v, size_t i, double x)
{
  if (!v || !v->data) return CH_INVALID;
  if (v->read_only) return CH_FROZEN;
  if (i >= v->capacity) return CH_SIZE;
  return write_item(v->data,i*v->stride,v->type,x);
}
int ch_view_fill(ch_view *v, double x, int strict)
{
  if (!v) return CH_INVALID;
  if (v->read_only) return CH_FROZEN;
  for (size_t i = 0; i < (strict ? v->size : v->capacity); ++i) {
    int s = ch_view_set(v,i,x); if (s) return s;
  }
  return CH_OK;
}
int ch_view_copy(ch_view *a, const ch_view *b, int strict)
{
  if (!a || !b) return CH_INVALID;
  size_t n = strict ? a->size : a->capacity;
  if (a->size != b->size || n > b->capacity) return CH_SIZE;
  if (a->read_only) return CH_FROZEN;
  for (size_t i = 0; i < n; ++i) { int s = ch_view_set(a,i,ch_view_get(b,i)); if (s) return s; }
  return CH_OK;
}
int ch_views_map(size_t count, ch_view *v, ch_element_function fun, void *ctx, int strict)
{
  if (!count || !v || !fun) return CH_INVALID;
  size_t n = strict ? v[0].size : v[0].capacity;
  for (size_t j = 0; j < count; ++j) {
    if (v[j].size != v[0].size) return CH_SIZE;
    if (v[j].capacity < n) n = v[j].capacity; /* unpadded vectors bound a mixed call */
  }
  double *values = ch_allocate(count,sizeof(double)); if (!values) return CH_MEMORY;
  int status = CH_OK;
  for (size_t i = 0; i < n && !status; ++i) {
    for (size_t j = 0; j < count; ++j) values[j] = ch_view_get(v+j,i);
    fun(values,count,ctx);
    for (size_t j = 0; j < count && !status; ++j) if (!v[j].read_only) status = ch_view_set(v+j,i,values[j]);
  }
  free(values); return status;
}
int ch_view_reduce(const ch_view *v, ch_reduce_function fun, double *acc, void *ctx, int strict)
{
  if (!v || !fun || !acc) return CH_INVALID;
  for (size_t i = 0; i < (strict ? v->size : v->capacity); ++i) fun(ch_view_get(v,i),acc,ctx);
  return CH_OK;
}
int ch_array_map(size_t count, ch_array *const *a, const unsigned char *ro, ch_element_function fun, void *ctx)
{
  if (!count || !a || !a[0] || !fun) return CH_INVALID;
  for (size_t j = 1; j < count; ++j) { int s = compatible(a[0],a[j]); if (s) return s; }
  double *v = ch_allocate(count,sizeof(double)); if (!v) return CH_MEMORY;
  int status = CH_OK;
  for (size_t r = 0; r < a[0]->records && !status; ++r) for (size_t e = 0; e < a[0]->elements && !status; ++e) {
    for (size_t j = 0; j < count; ++j) v[j] = read_item(a[j]->data,position(a[j],r,e),a[j]->type);
    fun(v,count,ctx);
    for (size_t j = 0; j < count && !status; ++j) if (!ro || !ro[j]) status = ch_array_set(a[j],r,e,v[j]);
  }
  free(v); return status;
}
void ch_array_function_free(ch_array_function *f)
{
  if (!f) return;
  for (size_t i = 0; i < f->count; ++i) free(f->signatures[i].indices);
  free(f->signatures); free(f);
}
int ch_array_function_create(size_t count, ch_array *const *a, const unsigned char *ro,
                             ch_group_function fun, void *ctx, ch_array_function **result)
{
  if (!count || !a || !fun || !result) return CH_INVALID;
  *result = NULL;
  ch_array_function *f = calloc(1,sizeof(*f)); if (!f) return CH_MEMORY;
  f->signatures = ch_allocate(count,sizeof(*f->signatures));
  if (!f->signatures) { free(f); return CH_MEMORY; }
  f->count = count; f->function = fun; f->context = ctx;
  int found = 0;
  for (size_t i = 0; i < count; ++i) {
    if (!a[i]) { ch_array_function_free(f); return CH_INVALID; }
    if (a[i]->kind != CH_ARRAY_VECTOR) {
      size_t width = a[i]->group > 1 ? a[i]->group : 1;
      if (found && (f->group > 1 ? f->group : 1) != width) { ch_array_function_free(f); return CH_INVALID; }
      f->group = width;
      found = 1;
    }
    ch_array_signature *s = f->signatures+i;
    *s = (ch_array_signature){a[i]->elements,a[i]->dimension,a[i]->group,a[i]->type,a[i]->kind,a[i]->order,NULL,ro ? ro[i] : 0};
    if (a[i]->indices) {
      size_t bytes = a[i]->dimension*a[i]->dimension*sizeof(int);
      s->indices = malloc(bytes ? bytes : 1);
      if (!s->indices) { ch_array_function_free(f); return CH_MEMORY; }
      memcpy(s->indices,a[i]->indices,bytes);
    }
  }
  *result = f; return CH_OK;
}
int ch_array_function_call(const ch_array_function *f, size_t count, ch_array *const *a)
{
  if (!f || !a || count != f->count) return CH_INVALID;
  for (size_t i = 0; i < count; ++i) {
    const ch_array_signature *s = f->signatures+i;
    if (!a[i]) return CH_INVALID;
    if (a[i]->records != a[0]->records || a[i]->elements != s->elements || a[i]->dimension != s->dimension)
      return CH_SIZE;
    if (a[i]->kind != s->kind || a[i]->type != s->type || a[i]->group != s->group || a[i]->order != s->order)
      return CH_INVALID;
    if (s->indices && memcmp(s->indices,a[i]->indices,s->dimension*s->dimension*sizeof(int))) return CH_INVALID;
  }
  size_t rows = a[0]->records, width = f->group ? f->group : rows;
  if (!rows) return CH_OK;
  ch_group *g = ch_allocate(count,sizeof(*g)); if (!g) return CH_MEMORY;
  int status = CH_OK;
  for (size_t start = 0; start < rows && !status; start += width) {
    size_t n = rows-start < width ? rows-start : width;
    for (size_t i = 0; i < count; ++i) g[i] = (ch_group){a[i],start,n,width,f->signatures[i].read_only};
    status = f->function(count,g,f->context);
  }
  free(g); return status;
}
