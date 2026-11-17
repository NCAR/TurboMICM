#include "rosenbrock.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void)
{
    const double atol[] = {1e-3,1e-3,1e-3,1e-3,1e-3};
    const double rtol = 1e-6;
    for (size_t cells=1; cells<=5; ++cells) {
        double old[25], next[25], error[25], expected=0;
        for (size_t i=0; i<cells; ++i)
            for (size_t j=0; j<5; ++j) {
                size_t k=i*5+j;
                old[k]=7.7*i+j; next[k]=-13.9/(j+1)-i; error[k]=81.57/(i+7)/(j+3);
                double scale=atol[j]+rtol*fmax(fabs(old[k]),fabs(next[k]));
                expected += error[k]*error[k]/(scale*scale);
            }
        expected=fmax(sqrt(expected/(cells*5)),1e-10);
        double actual=chem_rosenbrock_normalized_error(cells,5,old,next,error,atol,rtol);
        assert(fabs(actual-expected)/fmax(fabs(actual),fabs(expected)) <= 1e-11);
    }
    for (size_t cells=2;cells<=5;++cells) {
        double old[15],next[15],error[15],expected=0;
        for (size_t i=0;i<cells;++i) for (size_t j=0;j<3;++j) {
            size_t k=3*i+j;
            old[k]=1.0+i+.1*j; next[k]=.8+.5*i+.2*j; error[k]=.01*(1+i+j);
            double scale=atol[j]+rtol*fmax(fabs(old[k]),fabs(next[k]));
            expected+=error[k]*error[k]/(scale*scale);
        }
        expected=fmax(sqrt(expected/(3*cells)),1e-10);
        for (size_t method=0;method<5;++method) {
            double actual=chem_rosenbrock_normalized_error(cells,3,old,next,error,atol,rtol);
            assert(fabs(actual-expected)<=1e-12);
        }
    }
    puts("Rosenbrock normalized error PASS");
    return 0;
}
