/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "chem_matrix.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static size_t element_size(ChemMatrixType type)
{ return type == CHEM_MATRIX_REAL ? sizeof(double) : sizeof(int); }
static int valid_value(ChemMatrixType type, double v)
{ return type == CHEM_MATRIX_REAL || (isfinite(v) && v >= INT_MIN && v <= INT_MAX); }
static double read_value(const void *data, ChemMatrixType type, size_t i)
{ return type == CHEM_MATRIX_REAL ? ((const double *)data)[i] : ((const int *)data)[i]; }
static void write_value(void *data, ChemMatrixType type, size_t i, double v)
{
    if (type == CHEM_MATRIX_REAL) ((double *)data)[i] = v;
    else ((int *)data)[i] = (int)v;
}
ChemStatus chem_matrix_create(size_t rows, size_t columns, size_t group,
                              ChemMatrixType type, double initial, ChemMatrix *out)
{
    if (!out || (type != CHEM_MATRIX_REAL && type != CHEM_MATRIX_INT) || !valid_value(type,initial))
        return CHEM_INVALID;
    size_t padded = rows;
    if (group && rows % group) {
        size_t extra = group - rows % group;
        if (rows > SIZE_MAX - extra) return CHEM_MEMORY;
        padded += extra;
    }
    if (columns && padded > SIZE_MAX / columns) return CHEM_MEMORY;
    size_t size = padded * columns, width = element_size(type);
    if (size > SIZE_MAX / width) return CHEM_MEMORY;
    void *data = size ? malloc(size * width) : NULL;
    if (size && !data) return CHEM_MEMORY;
    *out = (ChemMatrix){rows,columns,group,size,type,data};
    return chem_matrix_fill(out,initial);
}
void chem_matrix_free(ChemMatrix *m)
{ if (m) { free(m->data); *m = (ChemMatrix){0}; } }
ChemStatus chem_matrix_copy(const ChemMatrix *m, ChemMatrix *out)
{
    if (!m || !out || m == out) return CHEM_INVALID;
    ChemStatus s = chem_matrix_create(m->rows,m->columns,m->group,m->type,0,out);
    if (s == CHEM_OK && m->size) memcpy(out->data,m->data,m->size * element_size(m->type));
    return s;
}
size_t chem_matrix_index(const ChemMatrix *m, size_t r, size_t c)
{
    return m->group ? (r / m->group) * (m->group * m->columns) + c * m->group + r % m->group
                    : r * m->columns + c;
}
size_t chem_matrix_row_stride(const ChemMatrix *m) { return m->group ? 1 : m->columns; }
size_t chem_matrix_column_stride(const ChemMatrix *m) { return m->group ? m->group : 1; }
size_t chem_matrix_group_count(const ChemMatrix *m)
{ return m->group ? m->rows / m->group + (m->rows % m->group != 0) : m->rows; }
size_t chem_matrix_group_size(const ChemMatrix *m)
{ return (m->group ? m->group : 1) * m->columns; }
const void *chem_matrix_const_at(const ChemMatrix *m, size_t r, size_t c)
{
    if (!m || r >= m->rows || c >= m->columns) return NULL;
    return (const unsigned char *)m->data + chem_matrix_index(m,r,c)*element_size(m->type);
}
void *chem_matrix_at(ChemMatrix *m, size_t r, size_t c)
{
    if (!m || r >= m->rows || c >= m->columns) return NULL;
    return (unsigned char *)m->data + chem_matrix_index(m,r,c)*element_size(m->type);
}
double chem_matrix_get(const ChemMatrix *m, size_t r, size_t c)
{
    const void *p = chem_matrix_const_at(m,r,c);
    return p ? read_value(p,m->type,0) : NAN;
}
ChemStatus chem_matrix_set(ChemMatrix *m, size_t r, size_t c, double v)
{
    void *p = chem_matrix_at(m,r,c);
    if (!p) return CHEM_DIMENSION;
    if (!valid_value(m->type,v)) return CHEM_INVALID;
    write_value(p,m->type,0,v);
    return CHEM_OK;
}
ChemStatus chem_matrix_copy_row(const ChemMatrix *m, size_t r, void *out, size_t count)
{
    if (!m || r >= m->rows || count < m->columns) return CHEM_DIMENSION;
    if (m->columns && !out) return CHEM_INVALID;
    for (size_t c=0;c<m->columns;++c)
        write_value(out,m->type,c,chem_matrix_get(m,r,c));
    return CHEM_OK;
}
ChemStatus chem_matrix_assign_row(ChemMatrix *m, size_t r, const void *in, size_t count)
{
    if (!m || r >= m->rows || count < m->columns) return CHEM_DIMENSION;
    if (m->columns && !in) return CHEM_INVALID;
    for (size_t c=0;c<m->columns;++c)
        write_value(m->data,m->type,chem_matrix_index(m,r,c),read_value(in,m->type,c));
    return CHEM_OK;
}
ChemStatus chem_matrix_from_rows(size_t rows, const size_t *lengths,
    const void *const *values, size_t group, ChemMatrixType type, ChemMatrix *out)
{
    if (!out || (rows && (!lengths || !values))) return CHEM_INVALID;
    size_t columns = rows ? lengths[0] : 0;
    for (size_t r=0;r<rows;++r) {
        if (lengths[r] != columns) return CHEM_DIMENSION;
        if (columns && !values[r]) return CHEM_INVALID;
    }
    ChemStatus s=chem_matrix_create(rows,columns,group,type,0,out);
    if (s != CHEM_OK) return s;
    for (size_t r=0;r<rows;++r) chem_matrix_assign_row(out,r,values[r],columns);
    return CHEM_OK;
}
ChemStatus chem_matrix_fill(ChemMatrix *m, double v)
{
    if (!m || !valid_value(m->type,v)) return CHEM_INVALID;
    for (size_t i=0;i<m->size;++i) write_value(m->data,m->type,i,v);
    return CHEM_OK;
}
static ChemStatus bound(ChemMatrix *m, double v, int maximum)
{
    if (!m || !valid_value(m->type,v)) return CHEM_INVALID;
    for (size_t i=0;i<m->size;++i) {
        double x=read_value(m->data,m->type,i);
        if (maximum ? x < v : v < x) write_value(m->data,m->type,i,v);
    }
    return CHEM_OK;
}
ChemStatus chem_matrix_max(ChemMatrix *m, double v) { return bound(m,v,1); }
ChemStatus chem_matrix_min(ChemMatrix *m, double v) { return bound(m,v,0); }
ChemStatus chem_matrix_for_each(ChemMatrix *m, size_t n,
    const ChemMatrix *const *inputs, ChemMatrixElementFunction fn, void *context)
{
    if (!m || !fn || (n && !inputs)) return CHEM_INVALID;
    for (size_t j=0;j<n;++j)
        if (!inputs[j] || inputs[j]->rows != m->rows || inputs[j]->columns != m->columns)
            return CHEM_DIMENSION;
    if (n > SIZE_MAX / sizeof(double)) return CHEM_MEMORY;
    double *values = n ? malloc(n*sizeof(double)) : NULL;
    if (n && !values) return CHEM_MEMORY;
    ChemStatus status=CHEM_OK;
    for (size_t r=0;r<m->rows && status==CHEM_OK;++r)
        for (size_t c=0;c<m->columns && status==CHEM_OK;++c) {
            for (size_t j=0;j<n;++j) values[j]=chem_matrix_get(inputs[j],r,c);
            status=chem_matrix_set(m,r,c,fn(chem_matrix_get(m,r,c),values,n,context));
        }
    free(values);
    return status;
}
static double axpy(double y, const double *x, size_t n, void *context)
{ (void)n; return y + *(double *)context * x[0]; }
ChemStatus chem_matrix_axpy(ChemMatrix *y, double alpha, const ChemMatrix *x)
{ return chem_matrix_for_each(y,1,&x,axpy,&alpha); }
ChemStatus chem_matrix_print(const ChemMatrix *m, FILE *out)
{
    if (!m || !out) return CHEM_INVALID;
    for (size_t r=0;r<m->rows;++r) {
        for (size_t c=0;c<m->columns;++c)
            if (fprintf(out,c ? ",%g" : "%g",chem_matrix_get(m,r,c)) < 0) return CHEM_INVALID;
        if (fputc('\n',out)==EOF) return CHEM_INVALID;
    }
    return CHEM_OK;
}

ChemArray chem_array_const_matrix(const ChemMatrix *m)
{ return (ChemArray){m->data,NULL,m->rows,m->columns,m->group,m->type,0,1}; }
ChemArray chem_array_matrix(ChemMatrix *m)
{ ChemArray a=chem_array_const_matrix(m);a.writable=m->data;a.read_only=0;return a; }
ChemArray chem_array_const_vector(const void *p,size_t n,ChemMatrixType type)
{ return (ChemArray){p,NULL,n,1,0,type,1,1}; }
ChemArray chem_array_vector(void *p,size_t n,ChemMatrixType type)
{ ChemArray a=chem_array_const_vector(p,n,type);a.writable=p;a.read_only=0;return a; }
ChemStatus chem_array_column(ChemArray a,size_t c,ChemView *out)
{
    if (!out) return CHEM_INVALID;
    if (c>=a.columns) return CHEM_DIMENSION;
    size_t group=a.group?a.group:1;
    size_t padded=a.rows+(a.rows%group?group-a.rows%group:0);
    size_t offset=(a.group?c*group:c)*element_size(a.type);
    *out=(ChemView){a.data?(const unsigned char *)a.data+offset:NULL,
        a.writable?(unsigned char *)a.writable+offset:NULL,a.rows,padded,group,
        a.columns*group,1,a.type};
    return CHEM_OK;
}
ChemView chem_view_const(ChemView v) { v.writable=NULL;return v; }
static size_t view_index(ChemView v,size_t i)
{ return i/v.group*v.block_stride+i%v.group*v.stride; }
double chem_view_get(ChemView v,size_t r)
{ return r<v.length?read_value(v.data,v.type,view_index(v,r)):NAN; }
ChemStatus chem_view_set(ChemView v,size_t r,double value)
{
    if (r>=v.length) return CHEM_DIMENSION;
    if (!v.writable || !valid_value(v.type,value)) return CHEM_INVALID;
    write_value(v.writable,v.type,view_index(v,r),value);
    return CHEM_OK;
}
ChemStatus chem_view_fill(ChemView v,double value)
{
    if ((v.length && !v.writable) || !valid_value(v.type,value)) return CHEM_INVALID;
    for (size_t i=0;i<v.length;++i) write_value(v.writable,v.type,view_index(v,i),value);
    return CHEM_OK;
}
ChemStatus chem_view_copy(ChemView dst,ChemView src)
{
    if (dst.length!=src.length) return CHEM_DIMENSION;
    if (dst.length && (!dst.writable || !src.data)) return CHEM_INVALID;
    for (size_t i=0;i<src.length;++i) {
        ChemStatus s=chem_view_set(dst,i,chem_view_get(src,i));
        if (s!=CHEM_OK) return s;
    }
    return CHEM_OK;
}
ChemStatus chem_views_for_each(size_t n,const ChemView *views,ChemRowFunction fn,void *ctx)
{
    if (!n || !views || !fn) return CHEM_INVALID;
    for (size_t k=0;k<n;++k) {
        if (views[k].length!=views[0].length) return CHEM_DIMENSION;
        if (views[k].length && !views[k].data) return CHEM_INVALID;
    }
    if (n>SIZE_MAX/sizeof(double)) return CHEM_MEMORY;
    double *values=malloc(n*sizeof(double));
    if (!values) return CHEM_MEMORY;
    ChemStatus status=CHEM_OK;
    for (size_t i=0;i<views[0].length && status==CHEM_OK;++i) {
        for (size_t k=0;k<n;++k) values[k]=chem_view_get(views[k],i);
        fn(values,n,ctx);
        for (size_t k=0;k<n && status==CHEM_OK;++k)
            if (views[k].writable) status=chem_view_set(views[k],i,values[k]);
    }
    free(values);
    return status;
}
ChemStatus chem_view_reduce(ChemView v,int strict,ChemReduceFunction fn,void *ctx,double *acc)
{
    if (!fn || !acc || (v.length && !v.data)) return CHEM_INVALID;
    size_t n=strict?v.length:v.padded_length;
    for (size_t i=0;i<n;++i) *acc=fn(read_value(v.data,v.type,view_index(v,i)),*acc,ctx);
    return CHEM_OK;
}
struct ChemArrayFunction {
    size_t count;
    ChemArray *signature;
    unsigned char *mutable;
    ChemArrayCallback callback;
    void *context;
};
ChemStatus chem_array_function_create(size_t n,const ChemArray *args,
    ChemArrayCallback fn,void *ctx,ChemArrayFunction **out)
{
    if (!n || !args || !fn || !out) return CHEM_INVALID;
    if (n>SIZE_MAX/sizeof(ChemArray)) return CHEM_MEMORY;
    ChemArrayFunction *f=calloc(1,sizeof(*f));
    if (!f) return CHEM_MEMORY;
    f->signature=malloc(n*sizeof(ChemArray));f->mutable=malloc(n);
    if (!f->signature || !f->mutable) { chem_array_function_free(f);return CHEM_MEMORY; }
    f->count=n;f->callback=fn;f->context=ctx;
    for (size_t i=0;i<n;++i) {
        f->signature[i]=args[i];
        f->mutable[i]=!args[i].read_only;
        /* Retain schema, never a borrowed data pointer. */
        f->signature[i].data=NULL;f->signature[i].writable=NULL;
    }
    *out=f;
    return CHEM_OK;
}
ChemStatus chem_array_function_call(const ChemArrayFunction *f,size_t n,const ChemArray *args)
{
    if (!f || !args) return CHEM_INVALID;
    if (n!=f->count) return CHEM_DIMENSION;
    for (size_t i=0;i<n;++i) {
        const ChemArray *s=&f->signature[i];
        if (args[i].rows!=args[0].rows || args[i].columns!=s->columns || args[i].vector!=s->vector)
            return CHEM_DIMENSION;
        if (args[i].type!=s->type || (args[i].rows && args[i].columns &&
            (!args[i].data || (f->mutable[i] && !args[i].writable)))) return CHEM_INVALID;
    }
    ChemArray *call=malloc(n*sizeof(ChemArray));
    if (!call) return CHEM_MEMORY;
    for (size_t i=0;i<n;++i) { call[i]=args[i];if (!f->mutable[i]) {call[i].writable=NULL;call[i].read_only=1;} }
    ChemStatus status=f->callback(call,n,f->context);
    free(call);
    return status;
}
void chem_array_function_free(ChemArrayFunction *f)
{ if (f) { free(f->signature);free(f->mutable);free(f); } }
