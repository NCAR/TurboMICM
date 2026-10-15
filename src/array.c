#include "chem_array.h"
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct ChemArray {
    size_t rows, columns, group_length, groups, size;
    ChemElementType type;
    void *data;
};
struct ChemArrayFunction {
    size_t count, *columns, group_length;
    ChemElementType *types;
    int *read_only;
    ChemArrayFunctionBody body;
    void *context;
};
static size_t element_size(ChemElementType t)
{
    switch (t) { case CHEM_REAL:return sizeof(double);case CHEM_INTEGER:return sizeof(int);case CHEM_BOOLEAN:return sizeof(_Bool); }
    return 0;
}
ChemStatus chem_array_create(size_t rows, size_t columns, size_t group, ChemElementType type, const void *initial, ChemArray **out)
{
    if (!out) return CHEM_INVALID;
    *out = NULL;
    size_t bytes = element_size(type);
    if (!bytes) return CHEM_INVALID;
    size_t groups = group ? rows / group + (rows % group != 0) : rows;
    if ((group && groups > SIZE_MAX / group) || (columns && (group ? groups * group : rows) > SIZE_MAX / columns)) return CHEM_INVALID;
    size_t size = (group ? groups * group : rows) * columns;
    if (size > SIZE_MAX / bytes) return CHEM_INVALID;
    ChemArray *a = calloc(1, sizeof(*a));
    if (!a) return CHEM_MEMORY;
    a->data = calloc(size ? size : 1, bytes);
    if (!a->data) { free(a); return CHEM_MEMORY; }
    a->rows=rows;a->columns=columns;a->group_length=group;a->groups=groups;a->size=size;a->type=type;
    if (initial) chem_array_fill(a,initial);
    *out=a;
    return CHEM_OK;
}
void chem_array_free(ChemArray *a) { if (a) { free(a->data);free(a); } }
ChemStatus chem_array_clone(const ChemArray *a, ChemArray **out)
{
    if (!a) return CHEM_INVALID;
    ChemStatus status=chem_array_create(a->rows,a->columns,a->group_length,a->type,NULL,out);
    if (status==CHEM_OK) memcpy((*out)->data,a->data,a->size*element_size(a->type));
    return status;
}
ChemStatus chem_array_move(ChemArray **to, ChemArray **from)
{
    if (!to || !from || to==from || *to) return CHEM_INVALID;
    *to=*from;*from=NULL;return CHEM_OK;
}
size_t chem_array_rows(const ChemArray *a) { return a->rows; }
size_t chem_array_columns(const ChemArray *a) { return a->columns; }
size_t chem_array_group_length(const ChemArray *a) { return a->group_length; }
size_t chem_array_groups(const ChemArray *a) { return a->groups; }
size_t chem_array_storage_size(const ChemArray *a) { return a->size; }
size_t chem_array_index(const ChemArray *a,size_t r,size_t c)
{
    if (!a || r>=a->rows || c>=a->columns) return SIZE_MAX;
    size_t g=a->group_length;
    return g ? (r/g)*g*a->columns+c*g+r%g : r*a->columns+c;
}
void *chem_array_data(ChemArray *a) { return a->data; }
const void *chem_array_const_data(const ChemArray *a) { return a->data; }
ChemStatus chem_array_get(const ChemArray *a,size_t r,size_t c,void *v)
{
    size_t i=chem_array_index(a,r,c);
    if (i==SIZE_MAX || !v) return CHEM_INVALID;
    memcpy(v,(const unsigned char *)a->data+i*element_size(a->type),element_size(a->type));return CHEM_OK;
}
ChemStatus chem_array_set(ChemArray *a,size_t r,size_t c,const void *v)
{
    size_t i=chem_array_index(a,r,c);
    if (i==SIZE_MAX || !v) return CHEM_INVALID;
    memcpy((unsigned char *)a->data+i*element_size(a->type),v,element_size(a->type));return CHEM_OK;
}
ChemStatus chem_array_fill(ChemArray *a,const void *v)
{
    if (!a || !v) return CHEM_INVALID;
    size_t bytes=element_size(a->type);
    for (size_t i=0;i<a->size;++i) memcpy((unsigned char *)a->data+i*bytes,v,bytes);
    return CHEM_OK;
}
ChemStatus chem_array_copy(ChemArray *to,const ChemArray *from)
{
    if (!to || !from || to->rows!=from->rows || to->columns!=from->columns || to->type!=from->type) return CHEM_INVALID;
    if (to==from) return CHEM_OK;
    if (to->group_length==from->group_length) {memcpy(to->data,from->data,to->size*element_size(to->type));return CHEM_OK;}
    size_t bytes=element_size(to->type);
    for (size_t r=0;r<to->rows;++r) for (size_t c=0;c<to->columns;++c)
        memcpy((unsigned char *)to->data+chem_array_index(to,r,c)*bytes,
            (const unsigned char *)from->data+chem_array_index(from,r,c)*bytes,bytes);
    return CHEM_OK;
}
ChemStatus chem_array_assign_row(ChemArray *a,size_t row,const void *v,size_t count)
{
    if (!a || row>=a->rows || count<a->columns || (a->columns && !v)) return CHEM_INVALID;
    size_t bytes=element_size(a->type);
    for (size_t j=0;j<a->columns;++j) chem_array_set(a,row,j,(const unsigned char *)v+j*bytes);
    return CHEM_OK;
}
ChemStatus chem_array_read_row(const ChemArray *a,size_t row,void *v,size_t count)
{
    if (!a || row>=a->rows || count<a->columns || (a->columns && !v)) return CHEM_INVALID;
    size_t bytes=element_size(a->type);
    for (size_t j=0;j<a->columns;++j) chem_array_get(a,row,j,(unsigned char *)v+j*bytes);
    return CHEM_OK;
}
ChemStatus chem_array_from_rows(size_t rows,const size_t *lengths,const void *const *data,size_t group,ChemElementType type,ChemArray **out)
{
    if (!out) return CHEM_INVALID;
    *out=NULL;
    if (rows && (!lengths || !data)) return CHEM_INVALID;
    size_t columns=rows?lengths[0]:0;
    for (size_t r=0;r<rows;++r) if (lengths[r]!=columns || (columns && !data[r])) return CHEM_INVALID;
    ChemStatus status=chem_array_create(rows,columns,group,type,NULL,out);
    if (status!=CHEM_OK) return status;
    for (size_t r=0;r<rows;++r) chem_array_assign_row(*out,r,data[r],columns);
    return CHEM_OK;
}
ChemStatus chem_array_axpy(ChemArray *a,double alpha,const ChemArray *x)
{
    if (!a || !x || a->type!=CHEM_REAL || x->type!=CHEM_REAL || a->rows!=x->rows || a->columns!=x->columns) return CHEM_INVALID;
    double *y=a->data;const double *v=x->data;
    if (a->group_length==x->group_length) {
        #pragma omp simd
        for (size_t i=0;i<a->size;++i) y[i]+=alpha*v[i];
    } else for (size_t r=0;r<a->rows;++r) for (size_t c=0;c<a->columns;++c)
        y[chem_array_index(a,r,c)]+=alpha*v[chem_array_index(x,r,c)];
    return CHEM_OK;
}
ChemStatus chem_array_limit(ChemArray *a,double v,int lower)
{
    if (!a || a->type!=CHEM_REAL) return CHEM_INVALID;
    double *data=a->data;
    for (size_t i=0;i<a->size;++i) data[i]=lower?(data[i]>v?data[i]:v):(data[i]<v?data[i]:v);
    return CHEM_OK;
}
size_t chem_array_format(const ChemArray *a,char *buffer,size_t capacity)
{
    if (!a || (capacity && !buffer)) return SIZE_MAX;
    size_t used=0;
    if (capacity) buffer[0]=0;
    if (!a->columns) {
        if (capacity) {
            size_t count=a->rows<capacity-1?a->rows:capacity-1;
            memset(buffer,'\n',count);buffer[count]=0;
        }
        return a->rows;
    }
    for (size_t r=0;r<a->rows;++r) for (size_t c=0;c<a->columns;++c) {
        char part[80];size_t index=chem_array_index(a,r,c);
        int length;
        if (a->type==CHEM_REAL) length=snprintf(part,sizeof(part),"%.6g%c",((double *)a->data)[index],c+1==a->columns?'\n':',');
        else if (a->type==CHEM_INTEGER) length=snprintf(part,sizeof(part),"%d%c",((int *)a->data)[index],c+1==a->columns?'\n':',');
        else length=snprintf(part,sizeof(part),"%d%c",(int)((_Bool *)a->data)[index],c+1==a->columns?'\n':',');
        if (length<0 || used>SIZE_MAX-(size_t)length) return SIZE_MAX;
        for (int j=0;j<length;++j) if (used+(size_t)j+1<capacity) buffer[used+(size_t)j]=part[j];
        used+=(size_t)length;
        if (capacity) buffer[used<capacity?used:capacity-1]=0;
    }
    return used;
}
ChemStatus chem_array_row(const ChemArray *a,size_t row,int ro,ChemView *out)
{
    if (!a || !out || row>=a->rows) return CHEM_INVALID;
    size_t index=a->columns?chem_array_index(a,row,0):0;
    *out=(ChemView){(unsigned char *)a->data+index*element_size(a->type),a->columns,a->columns,
        a->columns?a->columns:1,0,a->group_length?a->group_length:1,a->type,ro};return CHEM_OK;
}
ChemStatus chem_array_column(const ChemArray *a,size_t col,int ro,ChemView *out)
{
    if (!a || !out || col>=a->columns) return CHEM_INVALID;
    size_t g=a->group_length;
    *out=(ChemView){(unsigned char *)a->data+(a->rows?(g?col*g:col):0)*element_size(a->type),a->rows,
        g?a->groups*g:a->rows,g?g:1,g?g*a->columns:a->columns,1,a->type,ro};return CHEM_OK;
}
ChemStatus chem_array_group_column(const ChemArray *a,size_t group,size_t column,int ro,ChemView *out)
{
    if (!a || !out || column>=a->columns) return CHEM_INVALID;
    size_t g=a->group_length?a->group_length:1;
    if (group>=a->groups) return CHEM_INVALID;
    size_t first=group*g, logical=a->rows-first<g?a->rows-first:g;
    size_t index=chem_array_index(a,first,column);
    *out=(ChemView){(unsigned char *)a->data+index*element_size(a->type),logical,g,g,0,1,a->type,ro};return CHEM_OK;
}
ChemView chem_vector_view(void *data,size_t length,ChemElementType type,int ro)
{ return (ChemView){data,length,length,length?length:1,0,1,type,ro}; }
void *chem_view_element(const ChemView *v,size_t i)
{
    if (!v || !v->data || !v->block_length || i>=v->capacity || !element_size(v->type)) return NULL;
    size_t offset=(i/v->block_length)*v->block_stride+(i%v->block_length)*v->element_stride;
    return (unsigned char *)v->data+offset*element_size(v->type);
}
ChemStatus chem_view_fill(ChemView v,const void *value)
{
    if (!value || v.read_only) return CHEM_INVALID;
    for (size_t i=0;i<v.capacity;++i) {
        void *p=chem_view_element(&v,i);if (!p) return CHEM_INVALID;
        memcpy(p,value,element_size(v.type));
    }
    return CHEM_OK;
}
ChemStatus chem_view_copy(ChemView to,ChemView from)
{
    if (to.read_only || to.type!=from.type || to.length!=from.length) return CHEM_INVALID;
    /* A temporary also makes partially overlapping views well defined. */
    size_t bytes=element_size(to.type);
    if (!bytes || to.length>SIZE_MAX/bytes) return CHEM_INVALID;
    unsigned char *temp=malloc(to.length?to.length*bytes:1);
    if (!temp) return CHEM_MEMORY;
    for (size_t i=0;i<to.length;++i) {
        void *p=chem_view_element(&from,i);if (!p) {free(temp);return CHEM_INVALID;}
        memcpy(temp+i*bytes,p,bytes);
    }
    for (size_t i=0;i<to.length;++i) {
        void *p=chem_view_element(&to,i);if (!p) {free(temp);return CHEM_INVALID;}
        memcpy(p,temp+i*bytes,bytes);
    }
    free(temp);return CHEM_OK;
}
ChemStatus chem_views_apply(size_t count,const ChemView *views,size_t writable,ChemElementFunction fn,void *context)
{
    if (!count || !views || !fn || count>sizeof(size_t)*CHAR_BIT) return CHEM_INVALID;
    for (size_t j=0;j<count;++j)
        if (views[j].length!=views[0].length || ((writable>>j)&1 && views[j].read_only)) return CHEM_INVALID;
    void **elements=malloc(count*sizeof(*elements));if (!elements) return CHEM_MEMORY;
    for (size_t i=0;i<views[0].length;++i) {
        for (size_t j=0;j<count;++j) {
            elements[j]=chem_view_element(&views[j],i);
            if (!elements[j]) {free(elements);return CHEM_INVALID;}
        }
        fn(count,elements,context);
    }
    free(elements);return CHEM_OK;
}
ChemStatus chem_views_reduce(size_t count,const ChemView *views,int strict,ChemReductionFunction fn,void *acc,void *context)
{
    if (!count || !views || !fn || !acc || count>SIZE_MAX/sizeof(void *)) return CHEM_INVALID;
    size_t length=strict?views[0].length:views[0].capacity;
    for (size_t j=0;j<count;++j) if (views[j].length!=views[0].length || (!strict && views[j].capacity!=length)) return CHEM_INVALID;
    const void **elements=malloc(count*sizeof(*elements));if (!elements) return CHEM_MEMORY;
    for (size_t i=0;i<length;++i) {
        for (size_t j=0;j<count;++j) {
            elements[j]=chem_view_element(&views[j],i);
            if (!elements[j]) {free(elements);return CHEM_INVALID;}
        }
        fn(count,elements,acc,context);
    }
    free(elements);return CHEM_OK;
}
static int consistent(size_t count,const ChemArrayArgument *args)
{
    if (!count || !args || !args[0].array) return 0;
    for (size_t i=0;i<count;++i) if (!args[i].array || args[i].array->rows!=args[0].array->rows) return 0;
    return 1;
}
ChemStatus chem_array_function_create(size_t count,const ChemArrayArgument *args,ChemArrayFunctionBody body,void *context,ChemArrayFunction **out)
{
    if (!out) return CHEM_INVALID;
    *out=NULL;
    if (!count || !args || !body || count>SIZE_MAX/sizeof(ChemArrayArgument)) return CHEM_INVALID;
    for (size_t i=0;i<count;++i) {
        if (!args[i].array) return CHEM_INVALID;
        size_t first=args[0].array->group_length?args[0].array->group_length:1;
        size_t current=args[i].array->group_length?args[i].array->group_length:1;
        if (first!=current) return CHEM_INVALID;
    }
    ChemArrayFunction *f=calloc(1,sizeof(*f));if (!f) return CHEM_MEMORY;
    f->columns=malloc(count*sizeof(*f->columns));f->types=malloc(count*sizeof(*f->types));f->read_only=malloc(count*sizeof(*f->read_only));
    if (!f->columns || !f->types || !f->read_only) {chem_array_function_free(f);return CHEM_MEMORY;}
    f->count=count;f->group_length=args[0].array->group_length?args[0].array->group_length:1;f->body=body;f->context=context;
    for (size_t i=0;i<count;++i) {f->columns[i]=args[i].array->columns;f->types[i]=args[i].array->type;f->read_only[i]=args[i].read_only;}
    *out=f;return CHEM_OK;
}
ChemStatus chem_array_function_call(const ChemArrayFunction *f,size_t count,const ChemArrayArgument *args)
{
    if (!f || count!=f->count || !consistent(count,args)) return CHEM_INVALID;
    for (size_t i=0;i<count;++i)
        if ((args[i].array->group_length?args[i].array->group_length:1)!=f->group_length || f->columns[i]!=args[i].array->columns || f->types[i]!=args[i].array->type || (!f->read_only[i] && args[i].read_only)) return CHEM_INVALID;
    ChemArrayArgument *arguments=malloc(count*sizeof(*arguments));
    if (!arguments) return CHEM_MEMORY;
    for (size_t i=0;i<count;++i) {
        arguments[i]=args[i];arguments[i].read_only=f->read_only[i] || args[i].read_only;
    }
    ChemStatus status=f->body(count,arguments,f->context);
    free(arguments);return status;
}
void chem_array_function_free(ChemArrayFunction *f)
{ if (f) {free(f->columns);free(f->types);free(f->read_only);free(f);} }
