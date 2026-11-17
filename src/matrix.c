#include "chem.h"
#include "chem_matrix.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t position(const ChemMatrix *m, size_t row, size_t col)
{
    return m->group_size ? (row/m->group_size)*m->group_size*m->columns+col*m->group_size+row%m->group_size :
                           row*m->columns+col;
}

int chem_matrix_init(ChemMatrix *m, size_t rows, size_t columns, size_t group,
                     size_t element_size, const void *value)
{
    if (!m || !element_size) return CHEM_INVALID;
    memset(m,0,sizeof(*m));
    size_t padded=rows;
    if (group && rows%group) {
        if (rows>SIZE_MAX-(group-rows%group)) return CHEM_MEMORY;
        padded+=group-rows%group;
    }
    if (padded && columns>SIZE_MAX/padded) return CHEM_MEMORY;
    size_t slots=padded*columns;
    if (slots && element_size>SIZE_MAX/slots) return CHEM_MEMORY;
    m->data=calloc(slots ? slots : 1,element_size);
    if (!m->data) return CHEM_MEMORY;
    m->rows=rows; m->columns=columns; m->group_size=group; m->element_size=element_size; m->slots=slots;
    if (value) for (size_t i=0;i<slots;++i) memcpy(m->data+i*element_size,value,element_size);
    return CHEM_OK;
}

int chem_matrix_from_rows(ChemMatrix *m, const void *const *rows, const size_t *lengths,
                          size_t count, size_t group, size_t element_size)
{
    if (!m || (count && (!rows || !lengths))) return CHEM_INVALID;
    size_t columns=count ? lengths[0] : 0;
    for (size_t i=0;i<count;++i)
        if (lengths[i]!=columns || (columns && !rows[i])) return CHEM_INVALID;
    int status=chem_matrix_init(m,count,columns,group,element_size,NULL);
    if (status) return status;
    for (size_t i=0;i<count;++i) chem_matrix_set_row(m,i,rows[i],columns);
    return CHEM_OK;
}

void chem_matrix_destroy(ChemMatrix *m)
{ if (m) { free(m->data); memset(m,0,sizeof(*m)); } }

void *chem_matrix_at(ChemMatrix *m, size_t row, size_t column)
{ return m && m->data && row<m->rows && column<m->columns ? m->data+position(m,row,column)*m->element_size : NULL; }
const void *chem_matrix_at_const(const ChemMatrix *m, size_t row, size_t column)
{ return m && m->data && row<m->rows && column<m->columns ? m->data+position(m,row,column)*m->element_size : NULL; }
size_t chem_matrix_row_stride(const ChemMatrix *m) { return m ? (m->group_size ? 1 : m->columns) : 0; }
size_t chem_matrix_column_stride(const ChemMatrix *m) { return m ? (m->group_size ? m->group_size : 1) : 0; }

int chem_matrix_get_row(const ChemMatrix *m, size_t row, void *values, size_t count)
{
    if (!m || row>=m->rows || count<m->columns || (m->columns && !values)) return CHEM_INVALID;
    for (size_t i=0;i<m->columns;++i)
        memcpy((unsigned char *)values+i*m->element_size,chem_matrix_at_const(m,row,i),m->element_size);
    return CHEM_OK;
}
int chem_matrix_set_row(ChemMatrix *m, size_t row, const void *values, size_t count)
{
    if (!m || row>=m->rows || count<m->columns || (m->columns && !values)) return CHEM_INVALID;
    for (size_t i=0;i<m->columns;++i)
        memcpy(chem_matrix_at(m,row,i),(const unsigned char *)values+i*m->element_size,m->element_size);
    return CHEM_OK;
}
int chem_matrix_get_column(const ChemMatrix *m, size_t col, void *values, size_t count)
{
    if (!m || col>=m->columns || count!=m->rows || (m->rows && !values)) return CHEM_INVALID;
    for (size_t i=0;i<m->rows;++i)
        memcpy((unsigned char *)values+i*m->element_size,chem_matrix_at_const(m,i,col),m->element_size);
    return CHEM_OK;
}
int chem_matrix_set_column(ChemMatrix *m, size_t col, const void *values, size_t count)
{
    if (!m || col>=m->columns || count!=m->rows || (m->rows && !values)) return CHEM_INVALID;
    for (size_t i=0;i<m->rows;++i)
        memcpy(chem_matrix_at(m,i,col),(const unsigned char *)values+i*m->element_size,m->element_size);
    return CHEM_OK;
}
int chem_matrix_fill(ChemMatrix *m, const void *value)
{
    if (!m || !m->data || !value) return CHEM_INVALID;
    for (size_t i=0;i<m->slots;++i) memcpy(m->data+i*m->element_size,value,m->element_size);
    return CHEM_OK;
}
int chem_matrix_fill_column(ChemMatrix *m, size_t column, const void *value)
{
    if (!m || !m->data || column>=m->columns || !value) return CHEM_INVALID;
    for (size_t i=0;i<m->rows;++i) memcpy(chem_matrix_at(m,i,column),value,m->element_size);
    return CHEM_OK;
}
int chem_matrix_copy(ChemMatrix *dst, const ChemMatrix *src)
{
    if (!dst || !src || !dst->data || !src->data || dst->rows!=src->rows ||
        dst->columns!=src->columns || dst->element_size!=src->element_size) return CHEM_INVALID;
    for (size_t i=0;i<src->rows;++i) for (size_t j=0;j<src->columns;++j)
        memmove(chem_matrix_at(dst,i,j),chem_matrix_at_const(src,i,j),src->element_size);
    return CHEM_OK;
}
int chem_matrix_copy_column(ChemMatrix *dst, size_t dc, const ChemMatrix *src, size_t sc)
{
    if (!dst || !src || !dst->data || !src->data || dc>=dst->columns || sc>=src->columns ||
        dst->rows!=src->rows || dst->element_size!=src->element_size) return CHEM_INVALID;
    for (size_t i=0;i<src->rows;++i)
        memmove(chem_matrix_at(dst,i,dc),chem_matrix_at_const(src,i,sc),src->element_size);
    return CHEM_OK;
}
int chem_matrix_axpy(ChemMatrix *dst, double alpha, const ChemMatrix *src)
{
    if (!dst || !src || !dst->data || !src->data || dst->element_size!=sizeof(double) ||
        src->element_size!=sizeof(double) || dst->rows!=src->rows || dst->columns!=src->columns) return CHEM_INVALID;
    for (size_t i=0;i<dst->rows;++i) for (size_t j=0;j<dst->columns;++j)
        *(double *)chem_matrix_at(dst,i,j)+=alpha* *(const double *)chem_matrix_at_const(src,i,j);
    return CHEM_OK;
}
int chem_matrix_bound(ChemMatrix *m, double value, int minimum)
{
    if (!m || !m->data || m->element_size!=sizeof(double)) return CHEM_INVALID;
    for (size_t i=0;i<m->slots;++i) {
        double *x=(double *)(void *)(m->data+i*m->element_size);
        if (minimum ? *x<value : *x>value) *x=value;
    }
    return CHEM_OK;
}
int chem_matrix_format(const ChemMatrix *m, char *buffer, size_t capacity)
{
    if (!m || !m->data || !buffer || !capacity || m->element_size!=sizeof(double)) return CHEM_INVALID;
    size_t offset=0; buffer[0]='\0';
    for (size_t i=0;i<m->rows;++i) {
        for (size_t j=0;j<m->columns;++j) {
            int n=snprintf(buffer+offset,capacity-offset,j ? ",%g" : "%g",*(const double *)chem_matrix_at_const(m,i,j));
            if (n<0 || (size_t)n>=capacity-offset) return CHEM_INVALID;
            offset+=(size_t)n;
        }
        if (offset+1>=capacity) return CHEM_INVALID;
        buffer[offset++]='\n'; buffer[offset]='\0';
    }
    return CHEM_OK;
}
int chem_matrix_const_group_column(const ChemMatrix *m, size_t group, size_t column, ChemConstGroupView *view)
{
    if (!m || !view || !m->data || column>=m->columns) return CHEM_INVALID;
    size_t width=m->group_size ? m->group_size : 1;
    if (group>SIZE_MAX/width || group*width>=m->rows) return CHEM_INVALID;
    size_t first=group*width, length=m->rows-first<width ? m->rows-first : width;
    *view=(ChemConstGroupView){chem_matrix_at_const(m,first,column),length,m->element_size,CHEM_GROUPED_DENSE_VIEW};
    return CHEM_OK;
}
int chem_matrix_group_column(ChemMatrix *m, size_t group, size_t column, ChemGroupView *view)
{
    if (!view) return CHEM_INVALID;
    ChemConstGroupView read;
    int status=chem_matrix_const_group_column(m,group,column,&read);
    if (status) return status;
    *view=(ChemGroupView){(unsigned char *)(void *)read.base,read.length,read.element_size,read.category};
    return CHEM_OK;
}
int chem_group_fill(ChemGroupView view, const void *value)
{
    if (!view.base || !view.element_size || !value) return CHEM_INVALID;
    for (size_t i=0;i<view.length;++i) memcpy(view.base+i*view.element_size,value,view.element_size);
    return CHEM_OK;
}
int chem_group_copy(ChemGroupView destination, ChemConstGroupView source)
{
    if (!destination.base || !source.base || !destination.element_size ||
        destination.length!=source.length || destination.element_size!=source.element_size ||
        destination.length>SIZE_MAX/destination.element_size) return CHEM_INVALID;
    memmove(destination.base,source.base,destination.length*destination.element_size);
    return CHEM_OK;
}
int chem_matrix_reduce_column(const ChemMatrix *m, size_t column, int strict,
                              ChemMatrixReducer reducer, void *accumulator, void *context)
{
    if (!m || !m->data || column>=m->columns || !reducer || !accumulator) return CHEM_INVALID;
    size_t n=m->rows;
    if (!strict && m->group_size && n%m->group_size) n+=m->group_size-n%m->group_size;
    for (size_t i=0;i<n;++i) reducer(m->data+position(m,i,column)*m->element_size,accumulator,context);
    return CHEM_OK;
}
int chem_matrix_function_init(ChemMatrixFunction *f, ChemMatrix *const *matrices,
                              size_t count, size_t vectors, ChemMatrixFunctionBody body, void *context)
{
    if (!f || !matrices || !count || !body || count>SIZE_MAX/sizeof(*f->columns)) return CHEM_INVALID;
    memset(f,0,sizeof(*f));
    f->columns=malloc(count*sizeof(*f->columns));
    if (!f->columns) return CHEM_MEMORY;
    for (size_t i=0;i<count;++i) {
        if (!matrices[i] || !matrices[i]->data) { chem_matrix_function_destroy(f); return CHEM_INVALID; }
        size_t group=matrices[i]->group_size ? matrices[i]->group_size : 1;
        if (i && group!=f->group_size) { chem_matrix_function_destroy(f); return CHEM_INVALID; }
        f->group_size=group;
        f->columns[i]=matrices[i]->columns;
    }
    f->matrix_count=count; f->vector_count=vectors; f->body=body; f->context=context;
    return CHEM_OK;
}
void chem_matrix_function_destroy(ChemMatrixFunction *f)
{ if (f) { free(f->columns); memset(f,0,sizeof(*f)); } }
int chem_matrix_function_call(const ChemMatrixFunction *f, ChemMatrix *const *matrices,
                              void *const *vectors, const size_t *lengths)
{
    if (!f || !f->body || !matrices || !f->matrix_count || (f->vector_count && (!vectors || !lengths))) return CHEM_INVALID;
    if (!matrices[0] || !matrices[0]->data) return CHEM_INVALID;
    size_t rows=matrices[0]->rows;
    for (size_t i=0;i<f->matrix_count;++i)
        if (!matrices[i] || !matrices[i]->data || matrices[i]->rows!=rows || matrices[i]->columns!=f->columns[i] ||
            (matrices[i]->group_size ? matrices[i]->group_size : 1)!=f->group_size) return CHEM_INVALID;
    for (size_t i=0;i<f->vector_count;++i)
        if (lengths[i]!=rows || (rows && !vectors[i])) return CHEM_INVALID;
    return f->body(matrices,vectors,f->context);
}
