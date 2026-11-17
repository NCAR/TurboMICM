#include "chem.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void)
{
    /* A + 2 B -> 3 C, and C -> A. */
    const size_t reactants[] = {0, 1, 1, 2};
    const ChemProduct products[] = {{2, 3}, {0, 1}};
    const ChemReaction reactions[] = {{0, 3, 0, 1}, {3, 1, 1, 1}};
    const ChemMechanism m = {3, 2, reactions, reactants, products};
    const double y[] = {2, 4, 5}, k[] = {0.5, 2};
    double f[3], jac[9];
    chem_rhs(&m, y, k, f);
    assert(f[0] == -6 && f[1] == -32 && f[2] == 38);
    chem_jacobian(&m, y, k, jac);
    for (size_t col = 0; col < 3; ++col) {
        double x[] = {2, 4, 5}, fp[3], fm[3];
        const double h = 1e-5;
        x[col] += h; chem_rhs(&m, x, k, fp);
        x[col] -= 2*h; chem_rhs(&m, x, k, fm);
        for (size_t row = 0; row < 3; ++row)
            assert(fabs(jac[row*3+col] - (fp[row]-fm[row])/(2*h)) < 1e-8);
    }
    double zero[] = {2, 0, 5};
    chem_jacobian(&m, zero, k, jac);
    assert(jac[1] == 0 && jac[4] == 0 && jac[7] == 0);
    puts("mass-action kernel PASS");
    return 0;
}
