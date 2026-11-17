#ifndef CHEM_H
#define CHEM_H
#include <stddef.h>

/* Repeated reactants represent integer reaction orders.  Products carry yields. */
typedef struct { size_t species; double yield; } ChemProduct;
typedef struct {
    size_t reactant_begin, reactant_count;
    size_t product_begin, product_count;
} ChemReaction;
typedef struct {
    size_t species_count, reaction_count;
    const ChemReaction *reactions;
    const size_t *reactants;
    const ChemProduct *products;
} ChemMechanism;

/* k has one entry per reaction. Output arrays must not alias inputs. */
void chem_rhs(const ChemMechanism *m, const double *restrict y,
              const double *restrict k, double *restrict f);
/* Dense row-major Jacobian; its storage has species_count squared entries. */
void chem_jacobian(const ChemMechanism *m, const double *restrict y,
                   const double *restrict k, double *restrict jac);
/* Add kinetics to existing outputs, optionally leaving algebraic rows untouched.
   The Jacobian operation adds -df/dy, matching the solver convention. */
void chem_add_forcing(const ChemMechanism *m, const double *y, const double *k,
                      const unsigned char *algebraic, double *forcing);
void chem_subtract_jacobian(const ChemMechanism *m, const double *y, const double *k,
                            const unsigned char *algebraic, double *jac);
#endif
