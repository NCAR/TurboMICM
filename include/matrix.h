#ifndef CHEM_MATRIX_H
#define CHEM_MATRIX_H
#include <stddef.h>
#include <stdio.h>

/* width=0 uses row-major storage. A positive width interleaves that many
   rows within each group, padding the last group to the requested width. */
typedef struct {
    size_t rows, columns, width, groups, storage_size;
    double *data;
} ChemMatrix;
typedef enum { CHEM_MATRIX_OK, CHEM_MATRIX_INVALID, CHEM_MATRIX_SIZE, CHEM_MATRIX_ALLOCATION } ChemMatrixStatus;
/* Initialize with {0}; construction and copying replace existing storage only
   on success. Direct element access requires valid row and column indices. */
ChemMatrixStatus chem_matrix_init(ChemMatrix *m,size_t rows,size_t columns,size_t width,double value);
void chem_matrix_destroy(ChemMatrix *m);
ChemMatrixStatus chem_matrix_from_rows(ChemMatrix *m,size_t rows,size_t width,
                                       const double *const *values,const size_t *lengths);
ChemMatrixStatus chem_matrix_clone(ChemMatrix *destination,const ChemMatrix *source);
size_t chem_matrix_index(const ChemMatrix *m,size_t row,size_t column);
double *chem_matrix_at(ChemMatrix *m,size_t row,size_t column);
const double *chem_matrix_const_at(const ChemMatrix *m,size_t row,size_t column);
ChemMatrixStatus chem_matrix_set_row(ChemMatrix *m,size_t row,const double *values,size_t count);
ChemMatrixStatus chem_matrix_get_row(const ChemMatrix *m,size_t row,double *values,size_t count);
void chem_matrix_fill(ChemMatrix *m,double value);
void chem_matrix_max(ChemMatrix *m,double value);
void chem_matrix_min(ChemMatrix *m,double value);
typedef void (*ChemMatrixElementFunction)(double *element,const double *inputs,size_t count,void *context);
ChemMatrixStatus chem_matrix_apply(ChemMatrix *m,const ChemMatrix *const *inputs,size_t count,
                                    ChemMatrixElementFunction function,void *context);
ChemMatrixStatus chem_matrix_axpy(ChemMatrix *y,double alpha,const ChemMatrix *x);
ChemMatrixStatus chem_matrix_copy_column(ChemMatrix *destination,size_t dest_col,
                                          const ChemMatrix *source,size_t src_col);
ChemMatrixStatus chem_matrix_set_column(ChemMatrix *m,size_t col,const double *values,size_t count);
ChemMatrixStatus chem_matrix_get_column(const ChemMatrix *m,size_t col,double *values,size_t count);
/* Column views borrow a matrix or caller-owned vector. All views in an apply
   operation must have the same logical length; callbacks can be reused with
   other views and different row counts. */
typedef struct { ChemMatrix *matrix; size_t column; double *vector; size_t length; } ChemColumnView;
ChemMatrixStatus chem_matrix_column(ChemMatrix *m,size_t col,ChemColumnView *view);
ChemColumnView chem_vector_column(double *data,size_t length);
typedef void (*ChemRowFunction)(size_t row,double *const *values,void *context);
ChemMatrixStatus chem_columns_apply(const ChemColumnView *views,size_t count,ChemRowFunction function,void *context);
/* A padded reduction visits the tail group; a strict reduction visits only
   logical rows. Vector and row-major views have no padding. */
typedef void (*ChemReduceFunction)(double value,void *accumulator,void *context);
ChemMatrixStatus chem_column_reduce(const ChemColumnView *view,int include_padding,
                                     ChemReduceFunction function,void *accumulator,void *context);
ChemMatrixStatus chem_matrix_print(const ChemMatrix *m,FILE *stream);
#endif
