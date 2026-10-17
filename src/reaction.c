/* Reactions (MICM's ChemicalReaction and its builder). */
#include "internal.h"

tmc_reaction tmc_reaction_create(void)
{
  tmc_reaction r;
  memset(&r, 0, sizeof r);
  return r;
}

tmc_reaction tmc_reaction_copy(const tmc_reaction *src)
{
  tmc_reaction r = *src;
  r.reactants = NULL;
  r.products = NULL;
  if (src->n_reactants)
  {
    r.reactants = tmc_malloc((size_t)src->n_reactants * sizeof *r.reactants);
    memcpy(r.reactants, src->reactants, (size_t)src->n_reactants * sizeof *r.reactants);
  }
  if (src->n_products)
  {
    r.products = tmc_malloc((size_t)src->n_products * sizeof *r.products);
    memcpy(r.products, src->products, (size_t)src->n_products * sizeof *r.products);
  }
  return r;
}

void tmc_reaction_free(tmc_reaction *r)
{
  free(r->reactants);
  free(r->products);
  r->reactants = NULL;
  r->products = NULL;
  r->n_reactants = r->n_products = 0;
}

void tmc_reaction_add_reactant(tmc_reaction *r, const tmc_species *s)
{
  r->reactants = tmc_realloc(r->reactants, (size_t)(r->n_reactants + 1) * sizeof *r->reactants);
  r->reactants[r->n_reactants++] = *s;
}

void tmc_reaction_add_product(tmc_reaction *r, const tmc_species *s, double yield)
{
  r->products = tmc_realloc(r->products, (size_t)(r->n_products + 1) * sizeof *r->products);
  r->products[r->n_products].species = *s;
  r->products[r->n_products].coefficient = yield;
  ++r->n_products;
}

void tmc_reaction_set_rate_constant(tmc_reaction *r, tmc_rate_constant rc)
{
  r->rate_constant = rc;
  r->has_rate_constant = true;
}

void tmc_reaction_set_phase(tmc_reaction *r, const tmc_phase *p)
{
  tmc_strlcpy(r->phase, p->name, sizeof r->phase);
}

int tmc_reaction_validate(const tmc_reaction *r, tmc_error *err)
{
  if (!r->has_rate_constant)
    return tmc_error_set(err, TMC_ERROR_PROCESS, TMC_PROCESS_RATE_CONSTANT_IS_NOT_SET, "A rate constant is required for a reaction");
  return TMC_OK;
}

void tmc_reactions_free(tmc_reaction *reactions, int n)
{
  for (int i = 0; i < n; ++i)
    tmc_reaction_free(&reactions[i]);
  free(reactions);
}
