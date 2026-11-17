#include "typed_sparse_matrix.h"
#include <stdint.h>
#include <string.h>
void chem_typed_sparse_destroy(ChemTypedSparseMatrix *m)
{
    if (m) { chem_sparse_matrix_destroy(&m->structure); chem_typed_matrix_destroy(&m->values); *m=(ChemTypedSparseMatrix){0}; }
}
ChemMatrixStatus chem_typed_sparse_init(ChemTypedSparseMatrix *m,size_t n,size_t blocks,size_t width,
    ChemSparseOrdering ordering,const ChemSparseEntry *entries,size_t count,size_t size,const void *value)
{
    if (!m) return CHEM_MATRIX_INVALID;
    ChemTypedSparseMatrix copy={0};
    ChemMatrixStatus status=chem_sparse_matrix_init(&copy.structure,n,0,0,ordering,entries,count,0);
    if (status!=CHEM_MATRIX_OK) return status;
    status=chem_typed_matrix_init(&copy.values,blocks,copy.structure.nonzeros,width,size,value);
    if (status!=CHEM_MATRIX_OK) { chem_typed_sparse_destroy(&copy); return status; }
    chem_typed_sparse_destroy(m); *m=copy; return CHEM_MATRIX_OK;
}
ChemMatrixStatus chem_typed_sparse_clone(ChemTypedSparseMatrix *destination,const ChemTypedSparseMatrix *source)
{
    if (!destination || !source || !source->structure.offsets) return CHEM_MATRIX_INVALID;
    if (destination==source) return CHEM_MATRIX_OK;
    ChemTypedSparseMatrix copy={0}; ChemMatrixStatus status=chem_sparse_matrix_clone(&copy.structure,&source->structure);
    if (status!=CHEM_MATRIX_OK) return status;
    status=chem_typed_matrix_init(&copy.values,source->values.rows,source->values.columns,source->values.width,
                                  source->values.element_size,NULL);
    if (status!=CHEM_MATRIX_OK) { chem_typed_sparse_destroy(&copy); return status; }
    if (source->values.storage_size) memcpy(copy.values.data,source->values.data,
                                            source->values.storage_size*source->values.element_size);
    chem_typed_sparse_destroy(destination); *destination=copy; return CHEM_MATRIX_OK;
}
void *chem_typed_sparse_at(ChemTypedSparseMatrix *m,size_t block,size_t row,size_t col)
{
    if (!m) return NULL;
    size_t entry=chem_sparse_matrix_entry(&m->structure,row,col);
    return entry==SIZE_MAX ? NULL : chem_typed_matrix_at(&m->values,block,entry);
}
const void *chem_typed_sparse_const_at(const ChemTypedSparseMatrix *m,size_t block,size_t row,size_t col)
{
    if (!m) return NULL;
    size_t entry=chem_sparse_matrix_entry(&m->structure,row,col);
    return entry==SIZE_MAX ? NULL : chem_typed_matrix_const_at(&m->values,block,entry);
}
