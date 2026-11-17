#include "constraint_matrix.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

ChemConstraintStatus chem_constraint_set_apply_matrices(ChemConstraintSet set,
                                                        const double *temperatures, size_t cells,
                                                        const ChemMatrix *state, ChemMatrix *forcing,
                                                        ChemSparseMatrix *jacobian)
{
    if (!state || !forcing || !jacobian || state->rows!=cells || forcing->rows!=cells ||
        jacobian->blocks!=cells || state->columns!=forcing->columns ||
        state->columns!=jacobian->order || (cells && !temperatures)) return CHEM_CONSTRAINT_BAD_ROW;
    size_t n=state->columns;
    ChemConstraintStatus status=chem_constraint_set_validate(set,n);
    if (status || !set.count || !cells) return status;
    if (n>SIZE_MAX/n || n*n>SIZE_MAX-2*n || n*n+2*n>SIZE_MAX/sizeof(double)) return CHEM_CONSTRAINT_BAD_ROW;

    /* Reject an incompatible sparse pattern before modifying any matrix. */
    for (size_t i=0;i<set.count;++i) {
        const ChemConstraint *c=&set.items[i];
        if (!chem_sparse_at(jacobian,0,c->row,c->row)) return CHEM_CONSTRAINT_BAD_ROW;
        for (size_t j=0;j<c->reactant_count;++j)
            if (!chem_sparse_at(jacobian,0,c->row,(size_t)c->reactants[j].species)) return CHEM_CONSTRAINT_BAD_ROW;
        for (size_t j=0;j<c->product_count;++j)
            if (!chem_sparse_at(jacobian,0,c->row,(size_t)c->products[j].species)) return CHEM_CONSTRAINT_BAD_ROW;
    }
    double *scratch=malloc((n*n+2*n)*sizeof(*scratch));
    if (!scratch) return CHEM_CONSTRAINT_BAD_ROW;
    double *y=scratch,*f=scratch+n,*j=scratch+2*n;
    for (size_t cell=0;cell<cells;++cell) {
        if (chem_matrix_get_row(state,cell,y,n)) {status=CHEM_CONSTRAINT_BAD_ROW;break;}
        memset(f,0,n*sizeof(*f));memset(j,0,n*n*sizeof(*j));
        status=chem_constraint_set_apply(set,n,1,temperatures+cell,y,f,j);
        if (status) break;
        for (size_t i=0;i<set.count;++i) {
            size_t row=set.items[i].row;
            *chem_matrix_at(forcing,cell,row)=f[row];
            for (size_t col=0;col<n;++col) {
                double *entry=chem_sparse_at(jacobian,cell,row,col);
                if (entry) *entry=j[row*n+col];
            }
        }
    }
    free(scratch);return status;
}
