#include "typed_matrix.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static size_t position(size_t row,size_t columns,size_t column,size_t width)
{ return width ? (row/width)*width*columns+column*width+row%width : row*columns+column; }
void chem_typed_matrix_destroy(ChemTypedMatrix *m)
{ if (m) { free(m->data); *m=(ChemTypedMatrix){0}; } }
ChemMatrixStatus chem_typed_matrix_init(ChemTypedMatrix *m,size_t rows,size_t columns,size_t width,
                                        size_t size,const void *value)
{
    if (!m || !size) return CHEM_MATRIX_INVALID;
    size_t groups=width ? rows/width+(rows%width!=0) : rows;
    if (width && groups && width>SIZE_MAX/groups) return CHEM_MATRIX_SIZE;
    size_t padded=width ? groups*width : rows;
    if (padded && columns>SIZE_MAX/padded) return CHEM_MATRIX_SIZE;
    size_t count=padded*columns;
    if (count>SIZE_MAX/size) return CHEM_MATRIX_SIZE;
    unsigned char *data=count ? calloc(count,size) : NULL;
    if (count && !data) return CHEM_MATRIX_ALLOCATION;
    if (value) for (size_t i=0;i<count;++i) memcpy(data+i*size,value,size);
    chem_typed_matrix_destroy(m); *m=(ChemTypedMatrix){rows,columns,width,count,size,data};
    return CHEM_MATRIX_OK;
}
void *chem_typed_matrix_at(ChemTypedMatrix *m,size_t row,size_t column)
{ return !m || row>=m->rows || column>=m->columns ? NULL : (unsigned char *)m->data+position(row,m->columns,column,m->width)*m->element_size; }
const void *chem_typed_matrix_const_at(const ChemTypedMatrix *m,size_t row,size_t column)
{ return !m || row>=m->rows || column>=m->columns ? NULL : (const unsigned char *)m->data+position(row,m->columns,column,m->width)*m->element_size; }
ChemMatrixStatus chem_typed_const_column(const ChemTypedMatrix *m,size_t column,ChemConstRowView *view)
{
    if (!m || !view || column>=m->columns) return CHEM_MATRIX_INVALID;
    *view=(ChemConstRowView){m->data,m->rows,m->element_size,m->columns,column,m->width}; return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_typed_column(ChemTypedMatrix *m,size_t column,ChemRowView *view)
{
    if (!m || !view || column>=m->columns) return CHEM_MATRIX_INVALID;
    *view=(ChemRowView){m->data,m->rows,m->element_size,m->columns,column,m->width}; return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_double_const_column(const ChemMatrix *m,size_t column,ChemConstRowView *view)
{
    if (!m || !view || column>=m->columns) return CHEM_MATRIX_INVALID;
    *view=(ChemConstRowView){m->data,m->rows,sizeof(double),m->columns,column,m->width}; return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_double_column(ChemMatrix *m,size_t column,ChemRowView *view)
{
    if (!m || !view || column>=m->columns) return CHEM_MATRIX_INVALID;
    *view=(ChemRowView){m->data,m->rows,sizeof(double),m->columns,column,m->width}; return CHEM_MATRIX_OK;
}
ChemConstRowView chem_const_vector_view(const void *data,size_t length,size_t size)
{ return (ChemConstRowView){data,length,size,1,0,0}; }
ChemRowView chem_vector_view(void *data,size_t length,size_t size)
{ return (ChemRowView){data,length,size,1,0,0}; }
static int valid(size_t length,size_t size,size_t columns,size_t column,size_t width,const void *data)
{
    if (!size || !columns || column>=columns || (length && !data)) return 0;
    if (!length) return 1;
    size_t groups=width ? (length-1)/width+1 : length;
    if (width && groups>SIZE_MAX/width) return 0;
    size_t padded=width ? groups*width : length;
    return columns<=SIZE_MAX/padded && padded*columns<=SIZE_MAX/size;
}
ChemMatrixStatus chem_typed_rows_apply(size_t rows,const ChemConstRowView *inputs,size_t ni,
                                       const ChemRowView *outputs,size_t no,ChemTypedRowCallback callback,void *context)
{
    if (!callback || (ni && !inputs) || (no && !outputs)) return CHEM_MATRIX_INVALID;
    if (ni>SIZE_MAX/sizeof(void *) || no>SIZE_MAX/sizeof(void *)) return CHEM_MATRIX_SIZE;
    for (size_t i=0;i<ni;++i) {
        const ChemConstRowView *v=inputs+i;
        if (!valid(v->length,v->element_size,v->columns,v->column,v->width,v->data)) return CHEM_MATRIX_INVALID;
        if (v->length!=rows) return CHEM_MATRIX_SIZE;
    }
    for (size_t i=0;i<no;++i) {
        const ChemRowView *v=outputs+i;
        if (!valid(v->length,v->element_size,v->columns,v->column,v->width,v->data)) return CHEM_MATRIX_INVALID;
        if (v->length!=rows) return CHEM_MATRIX_SIZE;
    }
    const void **in=ni ? malloc(ni*sizeof(*in)) : NULL;
    void **out=no ? malloc(no*sizeof(*out)) : NULL;
    if ((ni && !in) || (no && !out)) { free(in); free(out); return CHEM_MATRIX_ALLOCATION; }
    for (size_t row=0;row<rows;++row) {
        for (size_t i=0;i<ni;++i) in[i]=(const unsigned char *)inputs[i].data+
            position(row,inputs[i].columns,inputs[i].column,inputs[i].width)*inputs[i].element_size;
        for (size_t i=0;i<no;++i) out[i]=(unsigned char *)outputs[i].data+
            position(row,outputs[i].columns,outputs[i].column,outputs[i].width)*outputs[i].element_size;
        callback(row,in,out,context);
    }
    free(in); free(out); return CHEM_MATRIX_OK;
}
