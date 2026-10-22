/* Development benchmark for separately identifying numerical kernel costs. */
#include "chem_specialized.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
extern const chem_specialized ts1_specialized;
extern const chem_mechanism ts1_mechanism;
int main(int argc,char **argv)
{
    size_t cells=argc>1 ? (size_t)strtoul(argv[1],NULL,10) : 32;
    size_t repeats=argc>2 ? (size_t)strtoul(argv[2],NULL,10) : 10000;
    chem_model m; chem_state s;
    if (chem_model_create(&m,ts1_specialized.mechanism,CHEM_LU_MOZART)!=CHEM_OK || chem_state_create(&s,&m,cells)!=CHEM_OK) return 1;
    chem_calculate_rates(&s);
    size_t nz=m.lu.nnz*cells,n=m.n_species*cells;
    double *matrix=calloc(nz,sizeof(double));
    if (!matrix) return 1;
    chem_jacobian(&m,cells,s.rates,s.y,matrix);
    for (size_t i=0;i<m.n_species;++i) for(size_t c=0;c<cells;++c) matrix[m.lu.diagonal[i]*cells+c]+=1e5;
    for(int mode=0;mode<2;++mode) {
        chem_kernel forcing=mode ? ts1_specialized.forcing : ts1_mechanism.forcing;
        chem_kernel jacobian=mode ? ts1_specialized.jacobian : ts1_mechanism.jacobian;
        for(int kind=0;kind<4;++kind) {
            double start=omp_get_wtime();
            for(size_t r=0;r<repeats;++r) {
                if(kind==0) { memset(s.work,0,n*sizeof(double)); forcing(&m,cells,s.rates,s.y,s.work); }
                if(kind==1) { memset(s.jacobian,0,nz*sizeof(double)); jacobian(&m,cells,s.rates,s.y,s.jacobian); }
                if(kind==2) { memcpy(s.jacobian,matrix,nz*sizeof(double)); if(mode) ts1_specialized.factor(cells,s.jacobian,s.scratch); else chem_lu_factor(&m.lu,cells,s.jacobian,s.scratch); }
                if(kind==3) { memcpy(s.work,s.y,n*sizeof(double)); if(mode) ts1_specialized.linear(cells,s.jacobian,s.work); else chem_lu_solve(&m.lu,cells,s.jacobian,s.work); }
            }
            printf("mode=%d kind=%d ms=%.3f\n",mode,kind,(omp_get_wtime()-start)*1000);
        }
    }
    free(matrix);chem_state_destroy(&s);chem_model_destroy(&m);return 0;
}
