#undef NDEBUG
#include <assert.h>
#include "sparse_matrix.h"
int main(void)
{
    const ChemSparseEntry entries[]={{0,0},{0,1},{0,2},{1,0},{1,1},{1,2},
        {2,0},{2,2},{3,1},{3,3},{4,0},{4,2},{4,4}};
    for (size_t width=0;width<=4;++width) for (int order=0;order<2;++order)
      for (size_t cells=1;cells<=4;++cells) {
        ChemSparseMatrix matrix={0};
        assert(chem_sparse_matrix_init(&matrix,5,cells,width,(ChemSparseOrdering)order,entries,13,0)==CHEM_MATRIX_OK);
        assert(matrix.values.rows==cells && matrix.dimension==5 && matrix.nonzeros==13 && matrix.values.storage_size>=13*cells);
        for (size_t c=0;c<cells;++c) {
            double values[]={-12.2,-24.3*(c+2),-42.3,-.43,-23.4,-83.4/(c+3),
                             -4.74,-6.91,-59.1,-83.4,-78.5,-53.6,-1};
            for (size_t i=0;i<13;++i) *chem_sparse_matrix_at(&matrix,c,entries[i].row,entries[i].column)=values[i];
        }
        chem_sparse_matrix_add_diagonal(&matrix,42.042);
        for (size_t c=0;c<cells;++c) {
            double expected[]={42.042-12.2,-24.3*(c+2),-42.3,-.43,42.042-23.4,-83.4/(c+3),
                               -4.74,42.042-6.91,-59.1,42.042-83.4,-78.5,-53.6,42.042-1};
            for (size_t i=0;i<13;++i)
                assert(*chem_sparse_matrix_const_at(&matrix,c,entries[i].row,entries[i].column)==expected[i]);
        }
        chem_sparse_matrix_destroy(&matrix);
    }
    return 0;
}
