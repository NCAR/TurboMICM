#include "chem.h"
#include <string.h>

void chem_rhs(const ChemMechanism *m, const double *restrict y,
              const double *restrict k, double *restrict f)
{
    memset(f, 0, m->species_count * sizeof(*f));
    for (size_t r = 0; r < m->reaction_count; ++r) {
        const ChemReaction *reaction = &m->reactions[r];
        double rate = k[r];
        for (size_t i = 0; i < reaction->reactant_count; ++i)
            rate *= y[m->reactants[reaction->reactant_begin + i]];
        for (size_t i = 0; i < reaction->reactant_count; ++i)
            f[m->reactants[reaction->reactant_begin + i]] -= rate;
        for (size_t i = 0; i < reaction->product_count; ++i) {
            const ChemProduct *p = &m->products[reaction->product_begin + i];
            f[p->species] += p->yield * rate;
        }
    }
}

void chem_jacobian(const ChemMechanism *m, const double *restrict y,
                   const double *restrict k, double *restrict jac)
{
    size_t n = m->species_count;
    memset(jac, 0, n * n * sizeof(*jac));
    for (size_t r = 0; r < m->reaction_count; ++r) {
        const ChemReaction *reaction = &m->reactions[r];
        for (size_t d = 0; d < reaction->reactant_count; ++d) {
            size_t col = m->reactants[reaction->reactant_begin + d];
            double derivative = k[r];
            /* Do not divide by y[col]: a reactant may be zero. */
            for (size_t i = 0; i < reaction->reactant_count; ++i)
                if (i != d) derivative *= y[m->reactants[reaction->reactant_begin + i]];
            for (size_t i = 0; i < reaction->reactant_count; ++i)
                jac[m->reactants[reaction->reactant_begin + i] * n + col] -= derivative;
            for (size_t i = 0; i < reaction->product_count; ++i) {
                const ChemProduct *p = &m->products[reaction->product_begin + i];
                jac[p->species * n + col] += p->yield * derivative;
            }
        }
    }
}
