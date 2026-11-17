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

void chem_add_forcing(const ChemMechanism *m, const double *y, const double *k,
                      const unsigned char *algebraic, double *forcing)
{
    for (size_t r = 0; r < m->reaction_count; ++r) {
        const ChemReaction *reaction = &m->reactions[r];
        double rate = k[r];
        for (size_t i = 0; i < reaction->reactant_count; ++i)
            rate *= y[m->reactants[reaction->reactant_begin + i]];
        for (size_t i = 0; i < reaction->reactant_count; ++i) {
            size_t row = m->reactants[reaction->reactant_begin + i];
            if (!algebraic || !algebraic[row]) forcing[row] -= rate;
        }
        for (size_t i = 0; i < reaction->product_count; ++i) {
            const ChemProduct *p = &m->products[reaction->product_begin + i];
            if (!algebraic || !algebraic[p->species]) forcing[p->species] += p->yield * rate;
        }
    }
}

void chem_subtract_jacobian(const ChemMechanism *m, const double *y, const double *k,
                            const unsigned char *algebraic, double *jac)
{
    size_t n = m->species_count;
    for (size_t r = 0; r < m->reaction_count; ++r) {
        const ChemReaction *reaction = &m->reactions[r];
        for (size_t d = 0; d < reaction->reactant_count; ++d) {
            size_t col = m->reactants[reaction->reactant_begin + d];
            double derivative = k[r];
            for (size_t i = 0; i < reaction->reactant_count; ++i)
                if (i != d) derivative *= y[m->reactants[reaction->reactant_begin + i]];
            for (size_t i = 0; i < reaction->reactant_count; ++i) {
                size_t row = m->reactants[reaction->reactant_begin + i];
                if (!algebraic || !algebraic[row]) jac[row * n + col] += derivative;
            }
            for (size_t i = 0; i < reaction->product_count; ++i) {
                const ChemProduct *p = &m->products[reaction->product_begin + i];
                if (!algebraic || !algebraic[p->species]) jac[p->species * n + col] -= p->yield * derivative;
            }
        }
    }
}

int chem_used_species(const ChemMechanism *m, unsigned char *used)
{
    if (!m || !used || (m->reaction_count && !m->reactions)) return 0;
    for (size_t r=0;r<m->reaction_count;++r) {
        const ChemReaction *reaction=&m->reactions[r];
        if ((reaction->reactant_count && !m->reactants) || (reaction->product_count && !m->products)) return 0;
        for (size_t i=0;i<reaction->reactant_count;++i) {
            size_t row=m->reactants[reaction->reactant_begin+i];
            if (row>=m->species_count) return 0;
            used[row]=1;
        }
        for (size_t i=0;i<reaction->product_count;++i) {
            size_t row=m->products[reaction->product_begin+i].species;
            if (row>=m->species_count) return 0;
            used[row]=1;
        }
    }
    return 1;
}

int chem_jacobian_pattern(const ChemMechanism *m, unsigned char *pattern)
{
    if (!m || !pattern || (m->reaction_count && !m->reactions)) return 0;
    const size_t n=m->species_count;
    if (n && n > (size_t)-1/n) return 0;
    for (size_t r=0;r<m->reaction_count;++r) {
        const ChemReaction *reaction=&m->reactions[r];
        if ((reaction->reactant_count && !m->reactants) || (reaction->product_count && !m->products)) return 0;
        for (size_t c=0;c<reaction->reactant_count;++c) {
            size_t col=m->reactants[reaction->reactant_begin+c];
            if (col>=n) return 0;
            for (size_t i=0;i<reaction->reactant_count;++i) {
                size_t row=m->reactants[reaction->reactant_begin+i];
                if (row>=n) return 0;
                pattern[row*n+col]=1;
            }
            for (size_t i=0;i<reaction->product_count;++i) {
                size_t row=m->products[reaction->product_begin+i].species;
                if (row>=n) return 0;
                pattern[row*n+col]=1;
            }
        }
    }
    return 1;
}

size_t chem_species_catalog_used_count(const ChemSpeciesCatalog *catalog)
{
    size_t count=0;
    if (!catalog || !catalog->used) return 0;
    for (size_t i=0;i<catalog->count;++i) count+=catalog->used[i]!=0;
    return count;
}
int chem_species_catalog_is_used(const ChemSpeciesCatalog *catalog, const char *name)
{
    if (!catalog || !catalog->names || !catalog->used || !name) return 0;
    for (size_t i=0;i<catalog->count;++i)
        if (catalog->names[i] && !strcmp(name,catalog->names[i])) return catalog->used[i]!=0;
    return 0;
}
