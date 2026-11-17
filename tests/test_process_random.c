#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include "chem.h"

static uint32_t rng=1234567;
static uint32_t next(void) { rng^=rng<<13; rng^=rng>>17; rng^=rng<<5; return rng; }
static double random_value(void) { return .01 + (next()%1000)/1000.0; }
static void run(size_t cells, size_t nr, size_t n)
{
    ChemReaction *reactions=calloc(nr,sizeof(*reactions));
    size_t *reactants=calloc(3*nr,sizeof(*reactants));
    ChemProduct *products=calloc(10*nr,sizeof(*products));
    double *y=calloc(n,sizeof(*y)), *k=calloc(nr,sizeof(*k));
    double *actual=calloc(n,sizeof(*actual)), *expected=calloc(n,sizeof(*expected));
    double *jac=calloc(n*n,sizeof(*jac)), *expected_jac=calloc(n*n,sizeof(*expected_jac));
    unsigned char *pattern=calloc(n*n,1);
    assert(reactions && reactants && products && y && k && actual && expected && jac && expected_jac && pattern);
    size_t nrct=0, nprd=0;
    for (size_t r=0;r<nr;++r) {
        size_t a=next()%4, b=next()%11;
        reactions[r]=(ChemReaction){nrct,a,nprd,b};
        for (size_t i=0;i<a;++i) reactants[nrct++]=next()%n;
        for (size_t i=0;i<b;++i) products[nprd++]=(ChemProduct){next()%n,1.2};
    }
    ChemMechanism m={n,nr,reactions,reactants,products};
    assert(chem_jacobian_pattern(&m,pattern));
    for (size_t cell=0;cell<cells;++cell) {
        for (size_t i=0;i<n;++i) { y[i]=random_value(); actual[i]=expected[i]=1000; }
        for (size_t i=0;i<nr;++i) k[i]=random_value();
        for (size_t i=0;i<n*n;++i) jac[i]=expected_jac[i]=100;
        chem_add_forcing(&m,y,k,NULL,actual);
        chem_subtract_jacobian(&m,y,k,NULL,jac);
        for (size_t r=0;r<nr;++r) {
            ChemReaction q=reactions[r]; double rate=k[r];
            for (size_t i=0;i<q.reactant_count;++i) rate*=y[reactants[q.reactant_begin+i]];
            for (size_t i=0;i<q.reactant_count;++i) expected[reactants[q.reactant_begin+i]]-=rate;
            for (size_t i=0;i<q.product_count;++i) expected[products[q.product_begin+i].species]+=1.2*rate;
            for (size_t i=0;i<q.reactant_count;++i) {
                size_t col=reactants[q.reactant_begin+i]; double deriv=k[r];
                for (size_t j=0;j<q.reactant_count;++j) if (i!=j) deriv*=y[reactants[q.reactant_begin+j]];
                for (size_t j=0;j<q.reactant_count;++j) expected_jac[reactants[q.reactant_begin+j]*n+col]+=deriv;
                for (size_t j=0;j<q.product_count;++j) expected_jac[products[q.product_begin+j].species*n+col]-=1.2*deriv;
            }
        }
        for (size_t i=0;i<n;++i) assert(fabs(actual[i]-expected[i])<1e-10);
        for (size_t i=0;i<n*n;++i) {
            assert(fabs(jac[i]-expected_jac[i])<1e-10);
            if (!pattern[i]) assert(jac[i]==100);
        }
    }
    free(reactions); free(reactants); free(products); free(y); free(k); free(actual); free(expected);
    free(jac); free(expected_jac); free(pattern);
}
int main(void)
{
    run(200,50,40); run(300,30,20); run(400,100,80);
    return 0;
}
