#include <stdio.h>
#include "matrix.h"
#include "rates.h"
#include "typed_sparse_matrix.h"

static int make_state(ChemMatrix *state,size_t width)
{
    const double values[3][3]={{1.1,2.1,3.1},{1.2,2.2,3.2},{1.3,2.3,3.3}};
    if (chem_matrix_init(state,3,3,width,0)!=CHEM_MATRIX_OK) return 0;
    for (size_t i=0;i<3;++i)
        if (chem_matrix_set_column(state,i,values[i],3)!=CHEM_MATRIX_OK) return 0;
    return 1;
}

/* The reference initializes the solver state to Converged, so its inner
   solve loop is skipped; it still updates rates ten times for each state. */
static void update_parameters(void)
{
    const ChemRate rates[3]={{.type=CHEM_RATE_USER,.custom_index=0},
                             {.type=CHEM_RATE_USER,.custom_index=1},
                             {.type=CHEM_RATE_USER,.custom_index=2}};
    const double custom[]={0.04,3e7,1e4};
    const ChemConditions conditions={272.5,101253.3,1e6};
    double calculated[3];
    for (size_t step=0;step<10;++step)
        for (size_t cell=0;cell<3;++cell) chem_rates(3,rates,conditions,custom,calculated);
}

int main(void)
{
    ChemMatrix regular={0},vector={0};
    ChemTypedSparseMatrix sparse={0},sparse_vector={0};
    const ChemSparseEntry entries[]={{0,1},{2,1},{2,3},{3,2}};
    const char *labels[3][4]={{"0.0.1","0.2.1","0.2.3","0.3.2"},
                              {"1.0.1","1.2.1","1.2.3","1.3.2"},
                              {"2.0.1","2.2.1","2.2.3","2.3.2"}};
    const char *empty="";
    int result=1;
    if (!make_state(&regular,0) || !make_state(&vector,3)) goto done;
    for (size_t i=0;i<regular.storage_size;++i) printf("%g\n",regular.data[i]);
    putchar('\n');
    for (size_t i=0;i<vector.storage_size;++i) printf("%g\n",vector.data[i]);
    putchar('\n');
    update_parameters(); update_parameters();
    for (size_t cell=0;cell<3;++cell) {
        printf("Cell %zu\n%10s%20s%20s\n",cell,"Species","Regular","Vectorized");
        for (size_t col=0;col<3;++col) {
            char species[2]={(char)('A'+col),'\0'};
            printf("%10s%20g%20g\n",species,*chem_matrix_at(&regular,cell,col),*chem_matrix_at(&vector,cell,col));
        }
        putchar('\n');
    }
    if (chem_typed_sparse_init(&sparse,4,3,0,CHEM_SPARSE_CSR,entries,4,sizeof(empty),&empty)!=CHEM_MATRIX_OK ||
        chem_typed_sparse_init(&sparse_vector,4,3,3,CHEM_SPARSE_CSR,entries,4,sizeof(empty),&empty)!=CHEM_MATRIX_OK) goto done;
    for (size_t block=0;block<3;++block) for (size_t i=0;i<4;++i) {
        *(const char **)chem_typed_sparse_at(&sparse,block,entries[i].row,entries[i].column)=labels[block][i];
        *(const char **)chem_typed_sparse_at(&sparse_vector,block,entries[i].row,entries[i].column)=labels[block][i];
    }
    puts("\nSparse matrix standard ordering elements");
    for (size_t i=0;i<sparse.values.storage_size;++i) puts(((const char **)sparse.values.data)[i]);
    puts("\nSparse matrix vector ordering elements");
    for (size_t i=0;i<sparse_vector.values.storage_size;++i) puts(((const char **)sparse_vector.values.data)[i]);
    result=0;
done:
    chem_typed_sparse_destroy(&sparse); chem_typed_sparse_destroy(&sparse_vector);
    chem_matrix_destroy(&regular); chem_matrix_destroy(&vector);
    return result;
}
