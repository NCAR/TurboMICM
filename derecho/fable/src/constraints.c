/* Equilibrium and linear constraints (MICM's EquilibriumConstraint and
 * LinearConstraint): construction, their residuals and Jacobian terms over a
 * state, and the equilibrium constant's temperature dependence. */
#include "internal.h"

#include <math.h>

#define TMC_GAS_CONSTANT 8.31446261815324

static tmc_stoich_species *copy_stoich(const tmc_stoich_species *src, int n)
{
  if (n <= 0)
    return NULL;
  tmc_stoich_species *p = tmc_malloc((size_t)n * sizeof *p);
  memcpy(p, src, (size_t)n * sizeof *p);
  return p;
}

int tmc_equilibrium_constraint(tmc_constraint *out, const char *name, const tmc_species *algebraic_species, const tmc_stoich_species *reactants, int n_reactants, const tmc_stoich_species *products, int n_products, tmc_vant_hoff vant_hoff, tmc_error *err)
{
  memset(out, 0, sizeof *out);
  if (n_reactants <= 0)
    return tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_EMPTY_REACTANTS, "Equilibrium constraint requires at least one reactant");
  if (n_products <= 0)
    return tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_EMPTY_PRODUCTS, "Equilibrium constraint requires at least one product");
  for (int i = 0; i < n_reactants; ++i)
    if (reactants[i].coefficient <= 0)
      return tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_INVALID_STOICHIOMETRY, "Stoichiometric coefficients must be positive");
  for (int i = 0; i < n_products; ++i)
    if (products[i].coefficient <= 0)
      return tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_INVALID_STOICHIOMETRY, "Stoichiometric coefficients must be positive");
  if (vant_hoff.K_HLC_ref <= 0)
    return tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_INVALID_EQUILIBRIUM_CONSTANT, "Henry's Law constant (K_HLC_ref) must be positive");
  out->type = TMC_CONSTRAINT_EQUILIBRIUM;
  tmc_strlcpy(out->name, name, sizeof out->name);
  out->algebraic_species = *algebraic_species;
  out->n_reactants = n_reactants;
  out->reactants = copy_stoich(reactants, n_reactants);
  out->n_products = n_products;
  out->products = copy_stoich(products, n_products);
  out->vant_hoff = vant_hoff;
  if (out->vant_hoff.R == 0.0)
    out->vant_hoff.R = TMC_GAS_CONSTANT;
  if (out->vant_hoff.T_ref == 0.0)
    out->vant_hoff.T_ref = 298.15;
  return TMC_OK;
}

int tmc_linear_constraint(tmc_constraint *out, const char *name, const tmc_species *algebraic_species, const tmc_stoich_species *terms, int n_terms, double constant, tmc_error *err)
{
  memset(out, 0, sizeof *out);
  if (n_terms <= 0)
    return tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_EMPTY_REACTANTS, "Linear constraint requires at least one term");
  out->type = TMC_CONSTRAINT_LINEAR;
  tmc_strlcpy(out->name, name, sizeof out->name);
  out->algebraic_species = *algebraic_species;
  out->n_terms = n_terms;
  out->terms = copy_stoich(terms, n_terms);
  out->constant = constant;
  return TMC_OK;
}

tmc_constraint tmc_constraint_copy(const tmc_constraint *c)
{
  tmc_constraint r = *c;
  r.reactants = copy_stoich(c->reactants, c->n_reactants);
  r.products = copy_stoich(c->products, c->n_products);
  r.terms = copy_stoich(c->terms, c->n_terms);
  return r;
}

void tmc_constraint_free(tmc_constraint *c)
{
  free(c->reactants);
  free(c->products);
  free(c->terms);
  c->reactants = c->products = c->terms = NULL;
  c->n_reactants = c->n_products = c->n_terms = 0;
}

int tmc_constraint_dependencies(const tmc_constraint *c, char (*names)[TMC_NAME_MAX], int max_names)
{
  int n = 0;
  if (c->type == TMC_CONSTRAINT_EQUILIBRIUM)
  {
    for (int i = 0; i < c->n_reactants; ++i, ++n)
      if (names && n < max_names)
        memcpy(names[n], c->reactants[i].species.name, TMC_NAME_MAX);
    for (int i = 0; i < c->n_products; ++i, ++n)
      if (names && n < max_names)
        memcpy(names[n], c->products[i].species.name, TMC_NAME_MAX);
  }
  else
  {
    for (int i = 0; i < c->n_terms; ++i, ++n)
      if (names && n < max_names)
        memcpy(names[n], c->terms[i].species.name, TMC_NAME_MAX);
  }
  return n;
}

int tmc_constraint_parameter_names(const tmc_constraint *c, char (*names)[TMC_NAME_MAX], int max_names)
{
  if (c->type == TMC_CONSTRAINT_EQUILIBRIUM)
  {
    if (names && max_names > 0)
      memcpy(names[0], c->name, TMC_NAME_MAX);
    return 1;
  }
  return 0;
}

double tmc_vant_hoff_constant(tmc_vant_hoff p, double temperature)
{
  double R = p.R == 0.0 ? TMC_GAS_CONSTANT : p.R;
  double T_ref = p.T_ref == 0.0 ? 298.15 : p.T_ref;
  return p.K_HLC_ref * exp((p.delta_H / R) * (1.0 / temperature - 1.0 / T_ref));
}

/* ---- evaluation over a state (plans built by the builder) --------------- */

static double clamp_pow(double conc, double stoich)
{
  return pow(conc > 0.0 ? conc : 0.0, stoich);
}

/* The constraints' state parameters (equilibrium constants) from the conditions. */
void tmc_constraints_update_parameters(const tmc_solver *s, const tmc_conditions *conditions, int n_cells, tmc_dense *params)
{
  for (int k = 0; k < s->n_constraints; ++k)
  {
    const tmc_constraint *c = &s->constraints[k];
    const tmc_constraint_plan *plan = &s->constraint_plans[k];
    if (c->type != TMC_CONSTRAINT_EQUILIBRIUM)
      continue;
    for (int cell = 0; cell < n_cells; ++cell)
      tmc_dense_set(params, cell, plan->parameter_index, tmc_vant_hoff_constant(c->vant_hoff, conditions[cell].temperature));
  }
}

/* forcing[row] = G(y) for every constraint (assigned, as MICM does). */
void tmc_constraints_add_forcing(const tmc_solver *s, const tmc_dense *variables, const tmc_dense *params, tmc_dense *forcing)
{
  const int L = variables->L, n_vars = s->n_vars, n_params = params->n_cols;
  for (int k = 0; k < s->n_constraints; ++k)
  {
    const tmc_constraint *c = &s->constraints[k];
    const tmc_constraint_plan *plan = &s->constraint_plans[k];
    for (int g = 0; g < variables->n_groups; ++g)
    {
      const double *y = variables->data + (size_t)g * n_vars * L;
      double *f = forcing->data + (size_t)g * n_vars * L;
      const double *p = params->data + (size_t)g * n_params * L;
      for (int lane = 0; lane < L; ++lane)
      {
        double value;
        if (c->type == TMC_CONSTRAINT_EQUILIBRIUM)
        {
          double rp = p[(size_t)plan->parameter_index * L + lane], pp = 1.0;
          for (int i = 0; i < c->n_reactants; ++i)
            rp *= clamp_pow(y[(size_t)plan->dependencies[i] * L + lane], c->reactants[i].coefficient);
          for (int i = 0; i < c->n_products; ++i)
            pp *= clamp_pow(y[(size_t)plan->dependencies[c->n_reactants + i] * L + lane], c->products[i].coefficient);
          value = rp - pp;
        }
        else
        {
          double sum = 0.0;
          for (int i = 0; i < c->n_terms; ++i)
            sum += c->terms[i].coefficient * y[(size_t)plan->dependencies[i] * L + lane];
          value = sum - c->constant;
        }
        f[(size_t)plan->row * L + lane] = value;
      }
    }
  }
}

/* jacobian -= dG/dy for every constraint, into the LU pattern's elements. */
void tmc_constraints_subtract_jacobian(const tmc_solver *s, const tmc_dense *variables, const tmc_dense *params, tmc_sparse *jacobian)
{
  const int L = variables->L, n_vars = s->n_vars, n_params = params->n_cols, nnz = jacobian->nnz;
  for (int k = 0; k < s->n_constraints; ++k)
  {
    const tmc_constraint *c = &s->constraints[k];
    const tmc_constraint_plan *plan = &s->constraint_plans[k];
    for (int g = 0; g < variables->n_groups; ++g)
    {
      const double *y = variables->data + (size_t)g * n_vars * L;
      double *J = jacobian->data + (size_t)g * nnz * L;
      const double *p = params->data + (size_t)g * n_params * L;
      for (int lane = 0; lane < L; ++lane)
      {
        if (c->type == TMC_CONSTRAINT_EQUILIBRIUM)
        {
          const double K_eq = p[(size_t)plan->parameter_index * L + lane];
          for (int i = 0; i < c->n_reactants; ++i)
          {
            const double stoich_i = c->reactants[i].coefficient;
            const double conc = y[(size_t)plan->dependencies[i] * L + lane];
            double prod = K_eq;
            for (int j = 0; j < c->n_reactants; ++j)
              if (j != i)
                prod *= clamp_pow(y[(size_t)plan->dependencies[j] * L + lane], c->reactants[j].coefficient);
            double partial;
            if (stoich_i == 1.0)
              partial = prod;
            else if (conc > 0.0)
              partial = stoich_i * prod * pow(conc, stoich_i - 1.0);
            else
              partial = 0.0;
            J[(size_t)plan->flat_ids[i] * L + lane] -= partial;
          }
          for (int i = 0; i < c->n_products; ++i)
          {
            const double stoich_i = c->products[i].coefficient;
            const double conc = y[(size_t)plan->dependencies[c->n_reactants + i] * L + lane];
            double prod = 1.0;
            for (int j = 0; j < c->n_products; ++j)
              if (j != i)
                prod *= clamp_pow(y[(size_t)plan->dependencies[c->n_reactants + j] * L + lane], c->products[j].coefficient);
            double partial;
            if (stoich_i == 1.0)
              partial = prod;
            else if (conc > 0.0)
              partial = stoich_i * prod * pow(conc, stoich_i - 1.0);
            else
              partial = 0.0;
            J[(size_t)plan->flat_ids[c->n_reactants + i] * L + lane] += partial;
          }
        }
        else
        {
          for (int i = 0; i < c->n_terms; ++i)
            J[(size_t)plan->flat_ids[i] * L + lane] -= c->terms[i].coefficient;
        }
      }
    }
  }
}

/* Yerror[row] = Ynew[row] - Y[row] for the algebraic rows (MICM's SetAlgebraicErrors). */
void tmc_set_algebraic_errors(const tmc_solver *s, tmc_dense *yerror, const tmc_dense *y, const tmc_dense *ynew)
{
  const int L = y->L, n_vars = s->n_vars;
  for (int g = 0; g < y->n_groups; ++g)
    for (int a = 0; a < s->n_algebraic; ++a)
    {
      size_t base = ((size_t)g * n_vars + (size_t)s->algebraic_ids[a]) * L;
      for (int lane = 0; lane < L; ++lane)
        yerror->data[base + lane] = ynew->data[base + lane] - y->data[base + lane];
    }
}
