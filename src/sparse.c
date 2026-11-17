#include "chem.h"
#include "chem_sparse.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { size_t row, col, major, minor; } Entry;
static int entry_compare(const void *a, const void *b)
{
    const Entry *x=a, *y=b;
    if (x->major!=y->major) return x->major<y->major ? -1 : 1;
    return x->minor==y->minor ? 0 : x->minor<y->minor ? -1 : 1;
}
int chem_sparse_init(ChemSparseMatrix *m, size_t order, size_t blocks,
                     const ChemSparseCoordinate *coords, size_t count,
                     ChemSparseOrder ordering, size_t group, size_t element_size, const void *initial)
{
    if (!m || (count && !coords) || (ordering!=CHEM_SPARSE_ROW && ordering!=CHEM_SPARSE_COLUMN)) return CHEM_INVALID;
    memset(m,0,sizeof(*m));
    if (count>SIZE_MAX/sizeof(Entry) || (order && order>SIZE_MAX/order) ||
        order*order>SIZE_MAX/sizeof(*m->lookup)) return CHEM_MEMORY;
    Entry *entries=malloc((count ? count : 1)*sizeof(*entries));
    if (!entries) return CHEM_MEMORY;
    for (size_t i=0;i<count;++i) {
        if (coords[i].row>=order || coords[i].column>=order) { free(entries); return CHEM_INVALID; }
        entries[i]=(Entry){coords[i].row,coords[i].column,
            ordering==CHEM_SPARSE_ROW ? coords[i].row : coords[i].column,
            ordering==CHEM_SPARSE_ROW ? coords[i].column : coords[i].row};
    }
    qsort(entries,count,sizeof(*entries),entry_compare);
    size_t unique=0;
    for (size_t i=0;i<count;++i) if (!i || entries[i].row!=entries[i-1].row || entries[i].col!=entries[i-1].col)
        entries[unique++]=entries[i];
    m->lookup=malloc((order ? order*order : 1)*sizeof(*m->lookup));
    m->coordinates=malloc((unique ? unique : 1)*sizeof(*m->coordinates));
    if (!m->lookup || !m->coordinates) { free(entries); chem_sparse_destroy(m); return CHEM_MEMORY; }
    for (size_t i=0;i<order*order;++i) m->lookup[i]=SIZE_MAX;
    for (size_t i=0;i<unique;++i) {
        m->coordinates[i]=(ChemSparseCoordinate){entries[i].row,entries[i].col};
        m->lookup[entries[i].row*order+entries[i].col]=i;
    }
    free(entries);
    int status=chem_matrix_init(&m->storage,blocks,unique,group,element_size,initial);
    if (status) { chem_sparse_destroy(m); return status; }
    m->order=order; m->nonzeros=unique; m->ordering=ordering;
    return CHEM_OK;
}
void chem_sparse_destroy(ChemSparseMatrix *m)
{
    if (!m) return;
    free(m->lookup); free(m->coordinates); chem_matrix_destroy(&m->storage); memset(m,0,sizeof(*m));
}
size_t chem_sparse_entry(const ChemSparseMatrix *m, size_t row, size_t col)
{ return m && m->lookup && row<m->order && col<m->order ? m->lookup[row*m->order+col] : SIZE_MAX; }
size_t chem_sparse_slot(const ChemSparseMatrix *m, size_t block, size_t row, size_t col)
{
    size_t entry=chem_sparse_entry(m,row,col);
    if (entry==SIZE_MAX || block>=m->storage.rows) return SIZE_MAX;
    const unsigned char *p=chem_matrix_at_const(&m->storage,block,entry);
    return (size_t)(p-m->storage.data)/m->storage.element_size;
}
void *chem_sparse_at(ChemSparseMatrix *m, size_t block, size_t row, size_t col)
{
    size_t entry=chem_sparse_entry(m,row,col);
    return entry==SIZE_MAX ? NULL : chem_matrix_at(&m->storage,block,entry);
}
const void *chem_sparse_at_const(const ChemSparseMatrix *m, size_t block, size_t row, size_t col)
{
    size_t entry=chem_sparse_entry(m,row,col);
    return entry==SIZE_MAX ? NULL : chem_matrix_at_const(&m->storage,block,entry);
}
void *chem_sparse_block_slot(ChemSparseMatrix *m, size_t block, size_t slot)
{
    if (!m) return NULL;
    size_t stride=m->storage.group_size ? m->storage.group_size : 1;
    if (slot%stride || slot/stride>=m->nonzeros) return NULL;
    return chem_matrix_at(&m->storage,block,slot/stride);
}
int chem_sparse_const_group_block(const ChemSparseMatrix *m, size_t group, size_t row, size_t column,
                                  ChemConstGroupView *view)
{
    size_t at=chem_sparse_entry(m,row,column);
    if (at==SIZE_MAX) return CHEM_INVALID;
    int status=chem_matrix_const_group_column(&m->storage,group,at,view);
    if (!status) view->category=CHEM_GROUPED_SPARSE_VIEW;
    return status;
}
int chem_sparse_group_block(ChemSparseMatrix *m, size_t group, size_t row, size_t column, ChemGroupView *view)
{
    size_t at=chem_sparse_entry(m,row,column);
    if (at==SIZE_MAX) return CHEM_INVALID;
    int status=chem_matrix_group_column(&m->storage,group,at,view);
    if (!status) view->category=CHEM_GROUPED_SPARSE_VIEW;
    return status;
}
int chem_sparse_fill(ChemSparseMatrix *m, const void *value)
{ return m ? chem_matrix_fill(&m->storage,value) : CHEM_INVALID; }
int chem_sparse_add_diagonal(ChemSparseMatrix *m, const void *value, ChemNumberType type)
{
    if (!m || !value || (type==CHEM_NUMBER_INT && m->storage.element_size!=sizeof(int)) ||
        (type==CHEM_NUMBER_DOUBLE && m->storage.element_size!=sizeof(double)) ||
        (type!=CHEM_NUMBER_INT && type!=CHEM_NUMBER_DOUBLE)) return CHEM_INVALID;
    for (size_t block=0;block<m->storage.rows;++block) for (size_t i=0;i<m->order;++i) {
        void *at=chem_sparse_at(m,block,i,i); if (!at) continue;
        if (type==CHEM_NUMBER_INT) *(int *)at+=*(const int *)value;
        else *(double *)at+=*(const double *)value;
    }
    return CHEM_OK;
}
int chem_sparse_format(const ChemSparseMatrix *m, char *buffer, size_t capacity, ChemNumberType type, int nonzero)
{
    if (!m || !m->storage.data || !buffer || !capacity ||
        (type==CHEM_NUMBER_INT && m->storage.element_size!=sizeof(int)) ||
        (type==CHEM_NUMBER_DOUBLE && m->storage.element_size!=sizeof(double)) ||
        (type!=CHEM_NUMBER_INT && type!=CHEM_NUMBER_DOUBLE)) return CHEM_INVALID;
    size_t used=0; buffer[0]='\0';
    for (size_t b=0;b<m->storage.rows;++b) {
        int n=snprintf(buffer+used,capacity-used,"Block %zu\n",b);
        if (n<0 || (size_t)n>=capacity-used) return CHEM_INVALID;
        used+=(size_t)n;
        for (size_t i=0;i<m->order;++i) for (size_t j=0;j<m->order;++j) {
            const void *p=chem_sparse_at_const(m,b,i,j);
            if (nonzero && !p) continue;
            if (nonzero) {
                n=type==CHEM_NUMBER_INT ? snprintf(buffer+used,capacity-used,"%zu, %zu, %d\n",i,j,*(const int *)p) :
                   snprintf(buffer+used,capacity-used,"%zu, %zu, %g\n",i,j,*(const double *)p);
            } else {
                n=type==CHEM_NUMBER_INT ? snprintf(buffer+used,capacity-used,j ? ",%d" : "%d",p ? *(const int *)p : 0) :
                   snprintf(buffer+used,capacity-used,j ? ",%g" : "%g",p ? *(const double *)p : 0.0);
            }
            if (n<0 || (size_t)n>=capacity-used) return CHEM_INVALID;
            used+=(size_t)n;
            if (!nonzero && j+1==m->order) {
                if (used+1>=capacity) return CHEM_INVALID;
                buffer[used++]='\n'; buffer[used]='\0';
            }
        }
    }
    return CHEM_OK;
}
