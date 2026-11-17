#ifndef CHEM_TYPED_MATRIX_H
#define CHEM_TYPED_MATRIX_H
#include "matrix.h"

/* Owned matrix of trivially copyable C values of an arbitrary element size. */
typedef struct { size_t rows, columns, width, storage_size, element_size; void *data; } ChemTypedMatrix;
ChemMatrixStatus chem_typed_matrix_init(ChemTypedMatrix *matrix,size_t rows,size_t columns,
                                        size_t width,size_t element_size,const void *initial_value);
void chem_typed_matrix_destroy(ChemTypedMatrix *matrix);
void *chem_typed_matrix_at(ChemTypedMatrix *matrix,size_t row,size_t column);
const void *chem_typed_matrix_const_at(const ChemTypedMatrix *matrix,size_t row,size_t column);

/* Row views carry element size and layout. Read views expose const pointers
   to callbacks; write views expose mutable pointers. Vector views borrow data. */
typedef struct { const void *data; size_t length,element_size,columns,column,width; } ChemConstRowView;
typedef struct { void *data; size_t length,element_size,columns,column,width; } ChemRowView;
ChemMatrixStatus chem_typed_const_column(const ChemTypedMatrix *matrix,size_t column,ChemConstRowView *view);
ChemMatrixStatus chem_typed_column(ChemTypedMatrix *matrix,size_t column,ChemRowView *view);
ChemMatrixStatus chem_double_const_column(const ChemMatrix *matrix,size_t column,ChemConstRowView *view);
ChemMatrixStatus chem_double_column(ChemMatrix *matrix,size_t column,ChemRowView *view);
ChemConstRowView chem_const_vector_view(const void *data,size_t length,size_t element_size);
ChemRowView chem_vector_view(void *data,size_t length,size_t element_size);
typedef void (*ChemTypedRowCallback)(size_t row,const void *const *inputs,void *const *outputs,void *context);
/* Reject invalid or unequal lengths before invoking the callback. */
ChemMatrixStatus chem_typed_rows_apply(size_t rows,const ChemConstRowView *inputs,size_t input_count,
                                       const ChemRowView *outputs,size_t output_count,
                                       ChemTypedRowCallback callback,void *context);
#endif
