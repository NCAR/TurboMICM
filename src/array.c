#include "chem_array.h"
#include <stdlib.h>
#include <string.h>

size_t chem_element_size(chem_element_type t)
{ return t==CHEM_DOUBLE ? sizeof(double) : t==CHEM_INT ? sizeof(int) : t==CHEM_INDEX ? sizeof(size_t) : 0; }
static double read_value(const void *data, size_t i, chem_element_type t)
{ return t==CHEM_DOUBLE ? ((const double *)data)[i] : t==CHEM_INT ? ((const int *)data)[i] : (double)((const size_t *)data)[i]; }
static void write_value(void *data, size_t i, chem_element_type t, double x)
{
    if (t==CHEM_DOUBLE) ((double *)data)[i]=x;
    else if (t==CHEM_INT) ((int *)data)[i]=(int)x;
    else ((size_t *)data)[i]=(size_t)x;
}
chem_status chem_array_create(chem_array *a,size_t rows,size_t cols,size_t lanes,chem_element_type t,double initial)
{
    *a=(chem_array){0};
    size_t unit=chem_element_size(t), padded=rows;
    if (!unit) return CHEM_INVALID;
    if (lanes && rows%lanes) {
        if (rows>SIZE_MAX-(lanes-rows%lanes)) return CHEM_SIZE;
        padded+=lanes-rows%lanes;
    }
    if ((cols && padded>SIZE_MAX/cols) || padded*cols>SIZE_MAX/unit) return CHEM_SIZE;
    *a=(chem_array){rows,cols,lanes,padded*cols,t,CHEM_DENSE_ARRAY,NULL};
    if (a->size && !(a->data=malloc(a->size*unit))) { *a=(chem_array){0}; return CHEM_MEMORY; }
    chem_array_fill(a,initial);
    return CHEM_OK;
}
chem_status chem_vector_create(chem_array *a,size_t n,size_t lanes,chem_element_type t,double initial)
{
    chem_status status=chem_array_create(a,n,1,lanes,t,initial);
    if (!status) a->kind=CHEM_VECTOR_ARRAY;
    return status;
}
void chem_array_destroy(chem_array *a) { free(a->data); *a=(chem_array){0}; }
chem_status chem_array_clone(chem_array *d,const chem_array *s)
{
    chem_status status=chem_array_create(d,s->rows,s->columns,s->vector_length,s->type,0);
    if (!status) { d->kind=s->kind; if (s->size) memcpy(d->data,s->data,s->size*chem_element_size(s->type)); }
    return status;
}
chem_status chem_array_import(chem_array *a,size_t rows,const size_t *lengths,const double *const *values,size_t lanes)
{
    size_t cols=rows ? lengths[0] : 0;
    for (size_t i=0;i<rows;++i) if (lengths[i]!=cols) return CHEM_SIZE;
    chem_status status=chem_array_create(a,rows,cols,lanes,CHEM_DOUBLE,0);
    if (!status) for (size_t i=0;i<rows;++i) chem_array_assign_row(a,i,values[i],cols);
    return status;
}
size_t chem_array_index(const chem_array *a,size_t row,size_t col)
{
    size_t v=a->vector_length;
    return v ? (row/v)*v*a->columns+col*v+row%v : row*a->columns+col;
}
size_t chem_array_row_stride(const chem_array *a) { return a->vector_length ? 1 : a->columns; }
size_t chem_array_column_stride(const chem_array *a) { return a->vector_length ? a->vector_length : 1; }
size_t chem_array_groups(const chem_array *a)
{ size_t v=a->vector_length ? a->vector_length : 1; return a->rows/v+(a->rows%v!=0); }
double chem_array_get(const chem_array *a,size_t row,size_t col)
{ return row<a->rows && col<a->columns ? read_value(a->data,chem_array_index(a,row,col),a->type) : 0; }
chem_status chem_array_set(chem_array *a,size_t row,size_t col,double x)
{
    if (row>=a->rows || col>=a->columns) return CHEM_SIZE;
    write_value(a->data,chem_array_index(a,row,col),a->type,x); return CHEM_OK;
}
chem_status chem_array_assign_row(chem_array *a,size_t row,const double *values,size_t n)
{
    if (row>=a->rows || n<a->columns) return CHEM_SIZE;
    for (size_t j=0;j<a->columns;++j) chem_array_set(a,row,j,values[j]);
    return CHEM_OK;
}
chem_status chem_array_export_row(const chem_array *a,size_t row,double *values,size_t n)
{
    if (row>=a->rows || n<a->columns) return CHEM_SIZE;
    for (size_t j=0;j<a->columns;++j) values[j]=chem_array_get(a,row,j);
    return CHEM_OK;
}
void chem_array_fill(chem_array *a,double x)
{ for (size_t i=0;i<a->size;++i) write_value(a->data,i,a->type,x); }
void chem_array_max(chem_array *a,double x)
{ for (size_t i=0;i<a->size;++i) if (read_value(a->data,i,a->type)<x) write_value(a->data,i,a->type,x); }
void chem_array_min(chem_array *a,double x)
{ for (size_t i=0;i<a->size;++i) if (read_value(a->data,i,a->type)>x) write_value(a->data,i,a->type,x); }
chem_status chem_array_axpy(chem_array *y,double alpha,const chem_array *x)
{
    if (y->rows!=x->rows || y->columns!=x->columns || y->vector_length!=x->vector_length) return CHEM_SIZE;
    for (size_t i=0;i<y->size;++i) write_value(y->data,i,y->type,read_value(y->data,i,y->type)+alpha*read_value(x->data,i,x->type));
    return CHEM_OK;
}
chem_status chem_array_print(FILE *f,const chem_array *a)
{
    for (size_t i=0;i<a->rows;++i) {
        for (size_t j=0;j<a->columns;++j) if (fprintf(f,"%s%g",j ? "," : "",chem_array_get(a,i,j))<0) return CHEM_INVALID;
        if (fputc('\n',f)==EOF) return CHEM_INVALID;
    }
    return CHEM_OK;
}
chem_status chem_array_get_group(const chem_array *a,size_t index,int writable,chem_array_group *g)
{
    *g=(chem_array_group){0};
    if (index>=chem_array_groups(a)) return CHEM_SIZE;
    size_t width=a->vector_length ? a->vector_length : 1,start=index*width;
    *g=(chem_array_group){a,index,start,width,a->rows-start<width ? a->rows-start : width,writable};
    return CHEM_OK;
}
chem_span chem_group_column(chem_array_group g,size_t col)
{
    if (!g.array || col>=g.array->columns) return (chem_span){.status=CHEM_SIZE};
    size_t offset=chem_array_index(g.array,g.start,col)*chem_element_size(g.array->type);
    const void *p=(const unsigned char *)g.array->data+offset;
    return (chem_span){p,g.writable ? (void *)p : NULL,chem_array_row_stride(g.array),g.width,g.active,g.array->type,CHEM_COLUMN_VIEW,CHEM_OK};
}
chem_span chem_group_const_column(chem_array_group g,size_t col)
{ g.writable=0; return chem_group_column(g,col); }
chem_span chem_group_vector(chem_array_group g)
{
    chem_span span=chem_group_column(g,0); span.kind=CHEM_VECTOR_VIEW;
    if (g.array->kind!=CHEM_VECTOR_ARRAY) span.status=CHEM_INVALID;
    return span;
}
chem_span chem_span_temporary(double *p,size_t width,size_t active)
{ return (chem_span){p,p,1,width,active,CHEM_DOUBLE,CHEM_TEMPORARY_VIEW,CHEM_OK}; }
chem_span chem_span_vector(const void *p,void *w,size_t n,chem_element_type type)
{ return (chem_span){p,w,1,n,n,type,CHEM_VECTOR_VIEW,CHEM_OK}; }
chem_span chem_span_const(chem_span s) { s.writable=NULL; return s; }
double chem_span_get(chem_span s,size_t i)
{ return !s.status && i<s.count ? read_value(s.data,i*s.stride,s.type) : 0; }
chem_status chem_span_set(chem_span s,size_t i,double x)
{
    if (s.status) return s.status;
    if (!s.writable) return CHEM_INVALID;
    if (i>=s.count) return CHEM_SIZE;
    write_value(s.writable,i*s.stride,s.type,x); return CHEM_OK;
}
chem_status chem_span_fill(chem_span d,double x)
{
    if (d.status) return d.status;
    if (!d.writable) return CHEM_INVALID;
    for (size_t i=0;i<d.count;++i) write_value(d.writable,i*d.stride,d.type,x);
    return CHEM_OK;
}
chem_status chem_span_copy(chem_span d,chem_span s)
{
    if (d.status || s.status) return d.status ? d.status : s.status;
    if (d.count!=s.count) return CHEM_SIZE;
    if (!d.writable) return CHEM_INVALID;
    for (size_t i=0;i<d.count;++i) write_value(d.writable,i*d.stride,d.type,read_value(s.data,i*s.stride,s.type));
    return CHEM_OK;
}
chem_status chem_span_map(chem_span d,size_t n,const chem_span *s,chem_map_function fn,void *ctx)
{
    if (d.status) return d.status;
    if (!d.writable || !fn) return CHEM_INVALID;
    for (size_t j=0;j<n;++j) { if (s[j].status) return s[j].status; if (s[j].count!=d.count) return CHEM_SIZE; }
    double local[8], *values=n<=8 ? local : malloc(n*sizeof(double));
    if (!values) return CHEM_MEMORY;
    for (size_t i=0;i<d.count;++i) {
        for (size_t j=0;j<n;++j) values[j]=chem_span_get(s[j],i);
        write_value(d.writable,i*d.stride,d.type,fn(values,n,ctx));
    }
    if (values!=local) free(values);
    return CHEM_OK;
}
chem_status chem_span_reduce(chem_span s,int strict,double *acc,chem_reduce_function fn,void *ctx)
{
    if (s.status) return s.status;
    if (!acc || !fn) return CHEM_INVALID;
    for (size_t i=0;i<(strict ? s.active : s.count);++i) *acc=fn(chem_span_get(s,i),*acc,ctx);
    return CHEM_OK;
}
chem_array_argument chem_array_mutable(chem_array *a) { return (chem_array_argument){a,1}; }
chem_array_argument chem_array_const(const chem_array *a) { return (chem_array_argument){a,0}; }
chem_status chem_array_function_create(chem_array_function *f,size_t n,const chem_array_argument *args,chem_group_function fn,void *ctx)
{
    *f=(chem_array_function){0};
    if (!n || !args || !fn) return CHEM_INVALID;
    size_t width=args[0].array->vector_length ? args[0].array->vector_length : 1;
    for (size_t i=1;i<n;++i)
        if ((args[i].array->vector_length ? args[i].array->vector_length : 1)!=width) return CHEM_SIZE;
    f->columns=malloc(n*sizeof(size_t)); f->vector_lengths=malloc(n*sizeof(size_t));
    f->types=malloc(n*sizeof(chem_element_type)); f->kinds=malloc(n*sizeof(chem_array_kind)); f->writable=malloc(n*sizeof(int));
    if (!f->columns || !f->vector_lengths || !f->types || !f->kinds || !f->writable) { chem_array_function_destroy(f); return CHEM_MEMORY; }
    f->count=n; f->function=fn; f->context=ctx;
    for (size_t i=0;i<n;++i) {
        f->columns[i]=args[i].array->columns; f->vector_lengths[i]=args[i].array->vector_length;
        f->types[i]=args[i].array->type; f->kinds[i]=args[i].array->kind; f->writable[i]=args[i].writable;
    }
    return CHEM_OK;
}
chem_status chem_array_function_apply(const chem_array_function *f,size_t n,const chem_array_argument *args)
{
    if (n!=f->count || !n || !args) return CHEM_SIZE;
    size_t rows=args[0].array->rows, width=args[0].array->vector_length;
    if (!width) width=1;
    for (size_t i=0;i<n;++i) {
        const chem_array *a=args[i].array;
        if (a->rows!=rows || a->columns!=f->columns[i] || a->vector_length!=f->vector_lengths[i]) return CHEM_SIZE;
        if ((a->vector_length ? a->vector_length : 1)!=width) return CHEM_SIZE;
        if (a->type!=f->types[i] || a->kind!=f->kinds[i] || (f->writable[i] && !args[i].writable)) return CHEM_INVALID;
    }
    chem_array_group local[8], *groups=n<=8 ? local : malloc(n*sizeof(*groups));
    if (!groups) return CHEM_MEMORY;
    chem_status status=CHEM_OK;
    for (size_t start=0,index=0;start<rows;start+=width,++index) {
        for (size_t i=0;i<n;++i) groups[i]=(chem_array_group){args[i].array,index,start,width,rows-start<width ? rows-start : width,f->writable[i]};
        status=f->function(n,groups,f->context); if (status) break;
    }
    if (groups!=local) free(groups);
    return status;
}
void chem_array_function_destroy(chem_array_function *f)
{
    free(f->columns); free(f->vector_lengths); free(f->types); free(f->kinds); free(f->writable);
    *f=(chem_array_function){0};
}
chem_raw_view chem_array_column(const chem_array *a,size_t col,int writable)
{ return (chem_raw_view){a,col,writable,CHEM_RAW_COLUMN_VIEW}; }
chem_span chem_bind_view(chem_raw_view v,chem_array_group g)
{
    if(v.array!=g.array) return (chem_span){.status=CHEM_INVALID};
    g.writable=g.writable && v.writable;
    size_t slot=v.slot;
    if(v.kind==CHEM_RAW_ELEMENT_VIEW) {
        size_t width=g.array->vector_length ? g.array->vector_length : 1;
        if(g.array->kind!=CHEM_SPARSE_ARRAY || slot%width) return (chem_span){.status=CHEM_INVALID};
        slot/=width;
    } else if(v.kind!=CHEM_RAW_COLUMN_VIEW) return (chem_span){.status=CHEM_INVALID};
    chem_span span=chem_group_column(g,slot);
    span.kind=v.kind==CHEM_RAW_ELEMENT_VIEW ? CHEM_ELEMENT_VIEW : CHEM_COLUMN_VIEW;
    return span;
}
