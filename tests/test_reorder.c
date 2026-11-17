#undef NDEBUG
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"
#include "reorder.h"

static size_t fill_count(const unsigned char *pattern,size_t n)
{
    unsigned char *edges=malloc(n*n); assert(edges); memcpy(edges,pattern,n*n);
    for (size_t pivot=0;pivot<n;++pivot)
        for (size_t row=pivot+1;row<n;++row) if (edges[row*n+pivot])
            for (size_t col=pivot+1;col<n;++col) if (edges[pivot*n+col]) edges[row*n+col]=1;
    size_t count=0; for (size_t i=0;i<n*n;++i) count+=edges[i]!=0;
    free(edges); return count;
}
static uint32_t next_random(uint32_t *state)
{
    uint32_t x=*state; x^=x<<13; x^=x>>17; x^=x<<5; return *state=x;
}
int main(void)
{
    enum { n=50 };
    for (size_t width=0;width<=4;++width) {
        ChemMatrix matrix={0}; assert(chem_matrix_init(&matrix,n,n,width,0)==CHEM_MATRIX_OK);
        unsigned char original[n*n],reordered[n*n]; uint32_t random=4723;
        for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) {
            double value=i==j || (next_random(&random)&1);
            *chem_matrix_at(&matrix,i,j)=value; original[i*n+j]=(unsigned char)value;
        }
        ChemPermutation p={0}; assert(chem_permutation_init(&p,n,original)==CHEM_REORDER_OK);
        unsigned char seen[n]; memset(seen,0,sizeof(seen));
        for (size_t i=0;i<n;++i) {
            assert(p.new_to_old[i]<n && !seen[p.new_to_old[i]]);
            seen[p.new_to_old[i]]=1;
            assert(p.old_to_new[p.new_to_old[i]]==i);
            for (size_t j=0;j<n;++j)
                reordered[i*n+j]=(unsigned char)*chem_matrix_at(&matrix,p.new_to_old[i],p.new_to_old[j]);
        }
        assert(fill_count(reordered,n)<fill_count(original,n));
        chem_permutation_destroy(&p); chem_matrix_destroy(&matrix);
    }
    return 0;
}
