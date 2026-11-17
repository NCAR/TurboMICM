#include "matrix.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void chem_matrix_destroy(ChemMatrix *m)
{
    if (m) { free(m->data); *m=(ChemMatrix){0}; }
}
ChemMatrixStatus chem_matrix_init(ChemMatrix *m,size_t rows,size_t columns,size_t width,double value)
{
    if (!m) return CHEM_MATRIX_INVALID;
    size_t groups=width ? rows/width+(rows%width!=0) : rows;
    if (width && groups && width>SIZE_MAX/groups) return CHEM_MATRIX_SIZE;
    size_t padded=width ? groups*width : rows;
    if (padded && columns>SIZE_MAX/padded) return CHEM_MATRIX_SIZE;
    size_t size=padded*columns;
    if (size>SIZE_MAX/sizeof(double)) return CHEM_MATRIX_SIZE;
    double *data=size ? malloc(size*sizeof(double)) : NULL;
    if (size && !data) return CHEM_MATRIX_ALLOCATION;
    for (size_t i=0;i<size;++i) data[i]=value;
    chem_matrix_destroy(m);
    *m=(ChemMatrix){rows,columns,width,groups,size,data};
    return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_matrix_from_rows(ChemMatrix *m,size_t rows,size_t width,
                                       const double *const *values,const size_t *lengths)
{
    if (!m || (rows && (!values || !lengths))) return CHEM_MATRIX_INVALID;
    size_t columns=rows ? lengths[0] : 0;
    for (size_t i=0;i<rows;++i) {
        if (lengths[i]!=columns) return CHEM_MATRIX_SIZE;
        if (columns && !values[i]) return CHEM_MATRIX_INVALID;
    }
    ChemMatrix copy={0}; ChemMatrixStatus status=chem_matrix_init(&copy,rows,columns,width,0);
    if (status!=CHEM_MATRIX_OK) return status;
    for (size_t i=0;i<rows;++i) for (size_t j=0;j<columns;++j)
        *chem_matrix_at(&copy,i,j)=values[i][j];
    chem_matrix_destroy(m); *m=copy;
    return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_matrix_clone(ChemMatrix *dst,const ChemMatrix *src)
{
    if (!dst || !src) return CHEM_MATRIX_INVALID;
    if (dst==src) return CHEM_MATRIX_OK;
    ChemMatrix copy={0};
    ChemMatrixStatus status=chem_matrix_init(&copy,src->rows,src->columns,src->width,0);
    if (status!=CHEM_MATRIX_OK) return status;
    if (src->storage_size) memcpy(copy.data,src->data,src->storage_size*sizeof(double));
    chem_matrix_destroy(dst); *dst=copy;
    return CHEM_MATRIX_OK;
}
size_t chem_matrix_index(const ChemMatrix *m,size_t row,size_t column)
{
    return m->width ? (row/m->width)*m->width*m->columns+column*m->width+row%m->width : row*m->columns+column;
}
double *chem_matrix_at(ChemMatrix *m,size_t row,size_t col)
{
    return !m || row>=m->rows || col>=m->columns ? NULL : m->data+chem_matrix_index(m,row,col);
}
const double *chem_matrix_const_at(const ChemMatrix *m,size_t row,size_t col)
{
    return !m || row>=m->rows || col>=m->columns ? NULL : m->data+chem_matrix_index(m,row,col);
}
ChemMatrixStatus chem_matrix_set_row(ChemMatrix *m,size_t row,const double *v,size_t count)
{
    if (!m || row>=m->rows || (m->columns && !v)) return CHEM_MATRIX_INVALID;
    if (count<m->columns) return CHEM_MATRIX_SIZE;
    for (size_t col=0;col<m->columns;++col) *chem_matrix_at(m,row,col)=v[col];
    return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_matrix_get_row(const ChemMatrix *m,size_t row,double *v,size_t count)
{
    if (!m || row>=m->rows || (m->columns && !v)) return CHEM_MATRIX_INVALID;
    if (count<m->columns) return CHEM_MATRIX_SIZE;
    for (size_t col=0;col<m->columns;++col) v[col]=*chem_matrix_const_at(m,row,col);
    return CHEM_MATRIX_OK;
}
void chem_matrix_fill(ChemMatrix *m,double value)
{ for (size_t i=0;i<m->storage_size;++i) m->data[i]=value; }
void chem_matrix_max(ChemMatrix *m,double value)
{ for (size_t i=0;i<m->storage_size;++i) m->data[i]=fmax(m->data[i],value); }
void chem_matrix_min(ChemMatrix *m,double value)
{ for (size_t i=0;i<m->storage_size;++i) m->data[i]=fmin(m->data[i],value); }
ChemMatrixStatus chem_matrix_apply(ChemMatrix *m,const ChemMatrix *const *inputs,size_t count,
                                    ChemMatrixElementFunction function,void *context)
{
    if (!m || !function || (count && !inputs)) return CHEM_MATRIX_INVALID;
    if (count>SIZE_MAX/sizeof(double)) return CHEM_MATRIX_SIZE;
    for (size_t i=0;i<count;++i) {
        if (!inputs[i]) return CHEM_MATRIX_INVALID;
        if (inputs[i]->rows!=m->rows || inputs[i]->columns!=m->columns) return CHEM_MATRIX_SIZE;
    }
    double *values=count ? malloc(count*sizeof(*values)) : NULL;
    if (count && !values) return CHEM_MATRIX_ALLOCATION;
    for (size_t row=0;row<m->rows;++row) for (size_t col=0;col<m->columns;++col) {
        for (size_t i=0;i<count;++i) values[i]=*chem_matrix_const_at(inputs[i],row,col);
        function(chem_matrix_at(m,row,col),values,count,context);
    }
    free(values); return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_matrix_axpy(ChemMatrix *y,double alpha,const ChemMatrix *x)
{
    if (!x || !y) return CHEM_MATRIX_INVALID;
    if (x->rows!=y->rows || x->columns!=y->columns) return CHEM_MATRIX_SIZE;
    for (size_t row=0;row<y->rows;++row) for (size_t col=0;col<y->columns;++col)
        *chem_matrix_at(y,row,col)+=alpha**chem_matrix_const_at(x,row,col);
    return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_matrix_copy_column(ChemMatrix *d,size_t dc,const ChemMatrix *s,size_t sc)
{
    if (!d || !s || dc>=d->columns || sc>=s->columns) return CHEM_MATRIX_INVALID;
    if (d->rows!=s->rows) return CHEM_MATRIX_SIZE;
    for (size_t row=0;row<d->rows;++row) *chem_matrix_at(d,row,dc)=*chem_matrix_const_at(s,row,sc);
    return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_matrix_set_column(ChemMatrix *m,size_t col,const double *v,size_t count)
{
    if (!m || col>=m->columns || (m->rows && !v)) return CHEM_MATRIX_INVALID;
    if (count!=m->rows) return CHEM_MATRIX_SIZE;
    for (size_t row=0;row<m->rows;++row) *chem_matrix_at(m,row,col)=v[row];
    return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_matrix_get_column(const ChemMatrix *m,size_t col,double *v,size_t count)
{
    if (!m || col>=m->columns || (m->rows && !v)) return CHEM_MATRIX_INVALID;
    if (count!=m->rows) return CHEM_MATRIX_SIZE;
    for (size_t row=0;row<m->rows;++row) v[row]=*chem_matrix_const_at(m,row,col);
    return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_matrix_print(const ChemMatrix *m,FILE *stream)
{
    if (!m || !stream) return CHEM_MATRIX_INVALID;
    for (size_t row=0;row<m->rows;++row) {
        for (size_t col=0;col<m->columns;++col)
            if (fprintf(stream,col ? ",%g" : "%g",*chem_matrix_const_at(m,row,col))<0) return CHEM_MATRIX_INVALID;
        if (fputc('\n',stream)==EOF) return CHEM_MATRIX_INVALID;
    }
    return CHEM_MATRIX_OK;
}
