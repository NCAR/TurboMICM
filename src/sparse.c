#include "chem_sparse.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct ChemSparse {
    size_t order, nonzeros, *map, *outer, *inner;
    ChemCompression compression;
    ChemElementType type;
    ChemArray *values;
};
ChemStatus chem_sparse_create(size_t order,size_t blocks,size_t group,ChemCompression compression,ChemElementType type,
    size_t count,const ChemMatrixIndex *entries,const void *initial,ChemSparse **out)
{
    if (!out) return CHEM_INVALID;
    *out=NULL;
    if ((count && !entries) || (compression!=CHEM_CSR && compression!=CHEM_CSC) ||
        (order && order>SIZE_MAX/order/sizeof(size_t)) || order==SIZE_MAX) return CHEM_INVALID;
    ChemSparse *s=calloc(1,sizeof(*s));if (!s) return CHEM_MEMORY;
    s->order=order;s->compression=compression;s->type=type;
    s->map=malloc((order?order*order:1)*sizeof(*s->map));s->outer=calloc(order+1,sizeof(*s->outer));
    if (!s->map || !s->outer) {chem_sparse_free(s);return CHEM_MEMORY;}
    for (size_t i=0;i<order*order;++i) s->map[i]=SIZE_MAX;
    for (size_t i=0;i<count;++i) {
        if (entries[i].row>=order || entries[i].column>=order) {chem_sparse_free(s);return CHEM_INVALID;}
        size_t *slot=&s->map[entries[i].row*order+entries[i].column];
        if (*slot==SIZE_MAX) { *slot=0;++s->nonzeros; }
    }
    s->inner=malloc((s->nonzeros?s->nonzeros:1)*sizeof(*s->inner));
    if (!s->inner) {chem_sparse_free(s);return CHEM_MEMORY;}
    size_t flat=0;
    for (size_t outer=0;outer<order;++outer) {
        s->outer[outer]=flat;
        for (size_t inner=0;inner<order;++inner) {
            size_t row=compression==CHEM_CSR?outer:inner,col=compression==CHEM_CSR?inner:outer;
            if (s->map[row*order+col]!=SIZE_MAX) {s->map[row*order+col]=flat;s->inner[flat++]=inner;}
        }
    }
    s->outer[order]=flat;
    ChemStatus status=chem_array_create(blocks,s->nonzeros,group,type,initial,&s->values);
    if (status!=CHEM_OK) {chem_sparse_free(s);return status;}
    *out=s;return CHEM_OK;
}
void chem_sparse_free(ChemSparse *s)
{ if (s) {free(s->map);free(s->outer);free(s->inner);chem_array_free(s->values);free(s);} }
ChemStatus chem_sparse_clone(const ChemSparse *source,ChemSparse **out)
{
    if (!source || !out) return CHEM_INVALID;
    *out=NULL;
    ChemSparse *s=calloc(1,sizeof(*s));if (!s) return CHEM_MEMORY;
    s->order=source->order;s->nonzeros=source->nonzeros;s->compression=source->compression;s->type=source->type;
    size_t n=s->order;
    s->map=malloc((n?n*n:1)*sizeof(*s->map));s->outer=malloc((n+1)*sizeof(*s->outer));
    s->inner=malloc((s->nonzeros?s->nonzeros:1)*sizeof(*s->inner));
    if (!s->map || !s->outer || !s->inner) {chem_sparse_free(s);return CHEM_MEMORY;}
    memcpy(s->map,source->map,n*n*sizeof(*s->map));memcpy(s->outer,source->outer,(n+1)*sizeof(*s->outer));
    memcpy(s->inner,source->inner,s->nonzeros*sizeof(*s->inner));
    ChemStatus status=chem_array_clone(source->values,&s->values);
    if (status!=CHEM_OK) {chem_sparse_free(s);return status;}
    *out=s;return CHEM_OK;
}
size_t chem_sparse_order(const ChemSparse *s) {return s->order;}
size_t chem_sparse_nonzeros(const ChemSparse *s) {return s->nonzeros;}
ChemCompression chem_sparse_compression(const ChemSparse *s) {return s->compression;}
ChemArray *chem_sparse_values(ChemSparse *s) {return s->values;}
const size_t *chem_sparse_outer_starts(const ChemSparse *s) {return s->outer;}
const size_t *chem_sparse_inner_indices(const ChemSparse *s) {return s->inner;}
size_t chem_sparse_flat_index(const ChemSparse *s,size_t row,size_t col)
{ return s && row<s->order && col<s->order?s->map[row*s->order+col]:SIZE_MAX; }
size_t chem_sparse_index(const ChemSparse *s,size_t block,size_t row,size_t col)
{
    size_t flat=chem_sparse_flat_index(s,row,col);
    return flat==SIZE_MAX?SIZE_MAX:chem_array_index(s->values,block,flat);
}
size_t chem_sparse_single_index(const ChemSparse *s,size_t row,size_t col)
{ return s && chem_array_rows(s->values)==1?chem_sparse_index(s,0,row,col):SIZE_MAX; }
ChemStatus chem_sparse_get(const ChemSparse *s,size_t block,size_t row,size_t col,void *value)
{
    size_t flat=chem_sparse_flat_index(s,row,col);
    return flat==SIZE_MAX?CHEM_INVALID:chem_array_get(s->values,block,flat,value);
}
ChemStatus chem_sparse_set(ChemSparse *s,size_t block,size_t row,size_t col,const void *value)
{
    size_t flat=chem_sparse_flat_index(s,row,col);
    return flat==SIZE_MAX?CHEM_INVALID:chem_array_set(s->values,block,flat,value);
}
size_t chem_sparse_diagonal_indices(const ChemSparse *s,size_t block,size_t *indices,size_t capacity)
{
    if (!s || block>=chem_array_rows(s->values) || (capacity && !indices)) return SIZE_MAX;
    size_t count=0;
    for (size_t i=0;i<s->order;++i) if (s->map[i*s->order+i]!=SIZE_MAX) {
        if (count<capacity) indices[count]=chem_sparse_index(s,block,i,i);
        ++count;
    }
    return count;
}
ChemStatus chem_sparse_add_diagonal(ChemSparse *s,const void *value)
{
    if (!s || !value) return CHEM_INVALID;
    for (size_t block=0;block<chem_array_rows(s->values);++block)
        for (size_t i=0;i<s->order;++i) if (s->map[i*s->order+i]!=SIZE_MAX) {
            union {double d;int i;_Bool b;} v;
            chem_sparse_get(s,block,i,i,&v);
            switch(s->type) {
                case CHEM_REAL:v.d+=*(const double *)value;break;
                case CHEM_INTEGER:v.i+=*(const int *)value;break;
                case CHEM_BOOLEAN:v.b=(_Bool)(v.b+*(const _Bool *)value);break;
            }
            chem_sparse_set(s,block,i,i,&v);
        }
    return CHEM_OK;
}
ChemStatus chem_sparse_block_view(const ChemSparse *s,size_t row,size_t col,int ro,ChemView *out)
{
    size_t flat=chem_sparse_flat_index(s,row,col);
    return flat==SIZE_MAX?CHEM_INVALID:chem_array_column(s->values,flat,ro,out);
}
ChemStatus chem_sparse_index_view(const ChemSparse *s,size_t index,int ro,ChemView *out)
{
    if (!s) return CHEM_INVALID;
    size_t g=chem_array_group_length(s->values);if (!g) g=1;
    if (index%g || index/g>=s->nonzeros) return CHEM_INVALID;
    return chem_array_column(s->values,index/g,ro,out);
}
static int append(char *buffer,size_t capacity,size_t *used,const char *text)
{
    size_t length=strlen(text);
    if (*used>SIZE_MAX-length) return 0;
    if (buffer && *used<capacity) {
        size_t available=capacity-*used-1,n=length<available?length:available;
        memcpy(buffer+*used,text,n);buffer[*used+n]=0;
    }
    *used+=length;return 1;
}
size_t chem_sparse_format(const ChemSparse *s,char *buffer,size_t capacity,int nonzeros_only)
{
    if (!s || (capacity && !buffer)) return SIZE_MAX;
    size_t used=0;if (capacity) buffer[0]=0;
    for (size_t block=0;block<chem_array_rows(s->values);++block) {
        char part[160];snprintf(part,sizeof(part),"Block %zu\n",block);
        if (!append(buffer,capacity,&used,part)) return SIZE_MAX;
        for (size_t r=0;r<s->order;++r) for (size_t c=0;c<s->order;++c) {
            if (nonzeros_only && s->map[r*s->order+c]==SIZE_MAX) continue;
            union {double d;int i;_Bool b;} v={0};
            if (s->map[r*s->order+c]!=SIZE_MAX) chem_sparse_get(s,block,r,c,&v);
            if (nonzeros_only) {
                snprintf(part,sizeof(part),"%zu, %zu, ",r,c);
                if (!append(buffer,capacity,&used,part)) return SIZE_MAX;
            }
            char ending=nonzeros_only || c+1==s->order?'\n':',';
            if (s->type==CHEM_REAL) snprintf(part,sizeof(part),"%.6g%c",v.d,ending);
            else if (s->type==CHEM_INTEGER) snprintf(part,sizeof(part),"%d%c",v.i,ending);
            else snprintf(part,sizeof(part),"%d%c",(int)v.b,ending);
            if (!append(buffer,capacity,&used,part)) return SIZE_MAX;
        }
    }
    return used;
}
