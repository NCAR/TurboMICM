#include "matrix.h"
#include <stdint.h>
#include <stdlib.h>

int chem_matrix_init(ChemMatrix *m, size_t rows, size_t columns, size_t width, int grouped, double initial)
{
    if (!m || !width || (!grouped && width!=1)) return 1;
    if (rows > SIZE_MAX-(width-1)) return 1;
    size_t groups=(rows+width-1)/width;
    if (groups && (columns > SIZE_MAX/groups || groups*columns > SIZE_MAX/width || groups*columns*width > SIZE_MAX/sizeof(double))) return 1;
    size_t count=groups*columns*width;
    double *data=count ? malloc(count*sizeof(*data)) : NULL;
    if (count && !data) return 1;
    for (size_t i=0;i<count;++i) data[i]=initial;
    *m=(ChemMatrix){rows,columns,width,groups,count,grouped,data};
    return 0;
}

int chem_matrix_from_rows(ChemMatrix *m, size_t rows, const double *const *values, const size_t *lengths, size_t width, int grouped)
{
    if (rows && (!values || !lengths)) return 1;
    size_t columns=rows ? lengths[0] : 0;
    for (size_t i=0;i<rows;++i) if (lengths[i]!=columns || (columns && !values[i])) return 1;
    if (chem_matrix_init(m,rows,columns,width,grouped,0)) return 1;
    for (size_t i=0;i<rows;++i) chem_matrix_set_row(m,i,values[i],columns);
    return 0;
}

void chem_matrix_free(ChemMatrix *m)
{
    if (!m) return;
    free(m->data); *m=(ChemMatrix){0};
}

static size_t offset(const ChemMatrix *m, size_t r, size_t c)
{
    return (r/m->width)*m->columns*m->width+c*m->width+r%m->width;
}

double *chem_matrix_at(ChemMatrix *m, size_t row, size_t column)
{
    return m && row<m->rows && column<m->columns ? &m->data[offset(m,row,column)] : NULL;
}
const double *chem_matrix_const_at(const ChemMatrix *m, size_t row, size_t column)
{
    return m && row<m->rows && column<m->columns ? &m->data[offset(m,row,column)] : NULL;
}
int chem_matrix_group_column(ChemMatrix *m, size_t group, size_t column, ChemDenseGroup *view)
{
    if (!m || !view || group >= m->groups || column >= m->columns) return 1;
    size_t first = group*m->width, count = m->rows-first;
    if (count > m->width) count=m->width;
    *view=(ChemDenseGroup){chem_matrix_at(m,first,column),count}; return 0;
}
int chem_matrix_const_group_column(const ChemMatrix *m, size_t group, size_t column, ChemConstDenseGroup *view)
{
    if (!m || !view || group >= m->groups || column >= m->columns) return 1;
    size_t first = group*m->width, count = m->rows-first;
    if (count > m->width) count=m->width;
    *view=(ChemConstDenseGroup){chem_matrix_const_at(m,first,column),count}; return 0;
}
size_t chem_matrix_row_stride(const ChemMatrix *m) { return m->grouped ? 1 : m->columns; }
size_t chem_matrix_column_stride(const ChemMatrix *m) { return m->grouped ? m->width : 1; }
int chem_matrix_get_row(const ChemMatrix *m, size_t row, double *values, size_t count)
{
    if (!m || row>=m->rows || count<m->columns || (m->columns && !values)) return 1;
    for(size_t c=0;c<m->columns;++c) values[c]=*chem_matrix_const_at(m,row,c);
    return 0;
}
int chem_matrix_set_row(ChemMatrix *m, size_t row, const double *values, size_t count)
{
    if (!m || row>=m->rows || count<m->columns || (m->columns && !values)) return 1;
    for(size_t c=0;c<m->columns;++c) *chem_matrix_at(m,row,c)=values[c];
    return 0;
}
int chem_matrix_copy(ChemMatrix *out, const ChemMatrix *in)
{
    if (!out || !in || out->rows!=in->rows || out->columns!=in->columns) return 1;
    for(size_t r=0;r<out->rows;++r) for(size_t c=0;c<out->columns;++c) *chem_matrix_at(out,r,c)=*chem_matrix_const_at(in,r,c);
    return 0;
}
void chem_matrix_fill(ChemMatrix *m,double value) { for(size_t i=0;i<m->storage_size;++i) m->data[i]=value; }
void chem_matrix_max(ChemMatrix *m,double value) { for(size_t i=0;i<m->storage_size;++i) if(m->data[i]<value) m->data[i]=value; }
void chem_matrix_min(ChemMatrix *m,double value) { for(size_t i=0;i<m->storage_size;++i) if(m->data[i]>value) m->data[i]=value; }
int chem_matrix_axpy(ChemMatrix *out,double alpha,const ChemMatrix *in)
{
    if(!out || !in || out->rows!=in->rows || out->columns!=in->columns) return 1;
    for(size_t r=0;r<out->rows;++r) for(size_t c=0;c<out->columns;++c) *chem_matrix_at(out,r,c)+=alpha * *chem_matrix_const_at(in,r,c);
    return 0;
}
int chem_matrix_foreach(ChemMatrix *out,const ChemMatrix *a,const ChemMatrix *b,ChemElementFunction fn,void *context)
{
    if(!out || !fn || (a && (out->rows!=a->rows || out->columns!=a->columns)) ||
       (b && (out->rows!=b->rows || out->columns!=b->columns))) return 1;
    for(size_t r=0;r<out->rows;++r) for(size_t c=0;c<out->columns;++c)
        fn(chem_matrix_at(out,r,c),a ? chem_matrix_const_at(a,r,c) : NULL,b ? chem_matrix_const_at(b,r,c) : NULL,context);
    return 0;
}
int chem_matrix_print(const ChemMatrix *m,FILE *file)
{
    if (!m || !file) return 1;
    for(size_t r=0;r<m->rows;++r) {
        for(size_t c=0;c<m->columns;++c) if(fprintf(file,"%s%g",c ? "," : "",*chem_matrix_const_at(m,r,c))<0) return 1;
        if(fputc('\n',file)==EOF) return 1;
    }
    return 0;
}

int chem_matrix_function_init(ChemMatrixFunction *f, ChemMatrix *const *prototypes, size_t matrix_count,
                              size_t vector_count, ChemRowFunction fn, void *context)
{
    if (!f || !fn || (matrix_count && !prototypes) || matrix_count > SIZE_MAX/sizeof(size_t)) return 1;
    size_t *columns=matrix_count ? malloc(matrix_count*sizeof(*columns)) : NULL;
    if (matrix_count && !columns) return 1;
    for (size_t i=0;i<matrix_count;++i) {
        if (!prototypes[i]) { free(columns); return 1; }
        columns[i]=prototypes[i]->columns;
    }
    *f=(ChemMatrixFunction){matrix_count,vector_count,columns,fn,context};
    return 0;
}

void chem_matrix_function_free(ChemMatrixFunction *f)
{
    if (!f) return;
    free(f->columns); *f=(ChemMatrixFunction){0};
}

int chem_matrix_function_call(const ChemMatrixFunction *f, ChemMatrix *const *matrices, size_t matrix_count,
                              ChemVector *const *vectors, size_t vector_count)
{
    if (!f || !f->function || matrix_count!=f->matrix_count || vector_count!=f->vector_count ||
        (matrix_count && !matrices) || (vector_count && !vectors)) return 1;
    size_t rows=0; int have_rows=0;
    for (size_t i=0;i<matrix_count;++i) {
        if (!matrices[i] || matrices[i]->columns!=f->columns[i] || (have_rows && matrices[i]->rows!=rows)) return 1;
        rows=matrices[i]->rows; have_rows=1;
    }
    for (size_t i=0;i<vector_count;++i) {
        if (!vectors[i] || (vectors[i]->size && !vectors[i]->data) || (have_rows && vectors[i]->size!=rows)) return 1;
        rows=vectors[i]->size; have_rows=1;
    }
    for (size_t row=0;row<rows;++row) if (f->function(row,matrices,vectors,f->context)) return 1;
    return 0;
}

int chem_matrix_reduce_column(const ChemMatrix *m,size_t column,int strict,ChemMatrixReducer reducer,void *context)
{
    if (!m || column>=m->columns || !reducer) return 1;
    size_t rows=strict ? m->rows : m->groups*m->width;
    for (size_t r=0;r<rows;++r) reducer(m->data[offset(m,r,column)],context);
    return 0;
}
int chem_matrix_fill_column(ChemMatrix *m,size_t column,double value)
{
    if (!m || column>=m->columns) return 1;
    for (size_t r=0;r<m->rows;++r) *chem_matrix_at(m,r,column)=value;
    return 0;
}
int chem_matrix_copy_column(ChemMatrix *out,size_t to,const ChemMatrix *in,size_t from)
{
    if (!out || !in || out->rows!=in->rows || to>=out->columns || from>=in->columns) return 1;
    for (size_t r=0;r<out->rows;++r) *chem_matrix_at(out,r,to)=*chem_matrix_const_at(in,r,from);
    return 0;
}
int chem_matrix_column_to_vector(const ChemMatrix *m,size_t column,double *out,size_t size)
{
    if (!m || column>=m->columns || size!=m->rows || (size && !out)) return 1;
    for (size_t r=0;r<size;++r) out[r]=*chem_matrix_const_at(m,r,column);
    return 0;
}
int chem_matrix_vector_to_column(ChemMatrix *m,size_t column,const double *values,size_t size)
{
    if (!m || column>=m->columns || size!=m->rows || (size && !values)) return 1;
    for (size_t r=0;r<size;++r) *chem_matrix_at(m,r,column)=values[r];
    return 0;
}

int chem_int_matrix_init(ChemIntMatrix *m,size_t rows,size_t columns,size_t width,int initial)
{
    if (!m || !width || rows>SIZE_MAX-(width-1)) return 1;
    size_t groups=(rows+width-1)/width;
    if (groups && (columns>SIZE_MAX/groups || groups*columns>SIZE_MAX/width || groups*columns*width>SIZE_MAX/sizeof(int))) return 1;
    size_t count=groups*columns*width;
    int *data=count ? malloc(count*sizeof(*data)) : NULL;
    if (count && !data) return 1;
    for(size_t i=0;i<count;++i) data[i]=initial;
    *m=(ChemIntMatrix){rows,columns,width,groups,data};
    return 0;
}
void chem_int_matrix_free(ChemIntMatrix *m)
{
    if (!m) return;
    free(m->data); *m=(ChemIntMatrix){0};
}
int *chem_int_matrix_at(ChemIntMatrix *m,size_t row,size_t column)
{
    return m && row<m->rows && column<m->columns ? &m->data[(row/m->width)*m->columns*m->width+column*m->width+row%m->width] : NULL;
}
int chem_int_matrix_apply_vector(ChemIntMatrix *m,size_t column,const int *values,size_t size,ChemIntRowFunction fn,void *context)
{
    if (!m || column>=m->columns || size!=m->rows || (size && !values) || !fn) return 1;
    for(size_t r=0;r<size;++r) fn(chem_int_matrix_at(m,r,column),values[r],context);
    return 0;
}
