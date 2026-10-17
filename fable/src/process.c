/* Rate constants of a state and the generic forcing and Jacobian kernels
 * (MICM's ReactionRateConstantStore and ProcessSet). */
#include "internal.h"

void tmc_calculate_rate_constants(const tmc_solver *s, tmc_state *st)
{
  const int L = st->L;
  for (int cell = 0; cell < st->n_cells; ++cell)
  {
    const tmc_conditions *c = &st->conditions[cell];
    for (int r = 0; r < s->n_reactions; ++r)
    {
      const tmc_process_plan *p = &s->processes[r];
      double custom[2] = { 0.0, 0.0 };
      if (p->custom_index >= 0)
      {
        int n = tmc_rate_constant_custom_parameter_count(&s->reactions[r].rate_constant);
        for (int i = 0; i < n; ++i)
          custom[i] = tmc_dense_get(st->custom_rate_parameters, cell, p->custom_index + i);
      }
      double k = tmc_rate_constant_calculate(&s->reactions[r].rate_constant, c, custom);
      /* MICM multiplies parameterized (third-body) reactants into the rate constant */
      for (int i = 0; i < p->n_parameterized; ++i)
        k *= tmc_species_parameterized_value(&p->parameterized[i], c);
      tmc_dense_set(st->rate_constants, cell, r, k);
    }
  }
  (void)L;
}

void tmc_solver_update_state_parameters(const tmc_solver *s, tmc_state *st)
{
  tmc_calculate_rate_constants(s, st);
}

/* forcing += F(y), one reaction at a time: rate = k * prod(reactants); the
 * reactants lose the rate, the products gain yield * rate. */
void tmc_generic_forcing(const tmc_solver *s, const tmc_dense *rate_constants, const tmc_dense *variables, tmc_dense *forcing)
{
  const int L = variables->L, n_vars = s->n_vars, n_rxn = s->n_reactions;
  double *rate = tmc_malloc((size_t)L * sizeof(double));
  for (int g = 0; g < variables->n_groups; ++g)
  {
    const double *restrict y = variables->data + (size_t)g * (size_t)n_vars * L;
    const double *restrict rc = rate_constants->data + (size_t)g * (size_t)n_rxn * L;
    double *restrict f = forcing->data + (size_t)g * (size_t)n_vars * L;
    for (int r = 0; r < n_rxn; ++r)
    {
      const tmc_process_plan *p = &s->processes[r];
      const double *k = rc + (size_t)r * L;
      for (int c = 0; c < L; ++c)
        rate[c] = k[c];
      for (int i = 0; i < p->n_reactants; ++i)
      {
        const double *yi = y + (size_t)p->reactants[i] * L;
#pragma omp simd
        for (int c = 0; c < L; ++c)
          rate[c] *= yi[c];
      }
      for (int i = 0; i < p->n_reactants; ++i)
      {
        double *fi = f + (size_t)p->reactants[i] * L;
#pragma omp simd
        for (int c = 0; c < L; ++c)
          fi[c] -= rate[c];
      }
      for (int i = 0; i < p->n_products; ++i)
      {
        double *fi = f + (size_t)p->products[i] * L;
        const double yield = p->yields[i];
#pragma omp simd
        for (int c = 0; c < L; ++c)
          fi[c] += yield * rate[c];
      }
    }
  }
  free(rate);
}

/* jacobian -= dF/dy, one (independent variable, reaction) job at a time, in
 * MICM's order: d = k * prod(other reactants); reactant rows gain d, product
 * rows lose yield * d. */
void tmc_generic_jacobian(const tmc_solver *s, const tmc_dense *rate_constants, const tmc_dense *variables, tmc_sparse *jacobian)
{
  const int L = variables->L, n_vars = s->n_vars, n_rxn = s->n_reactions, nnz = jacobian->nnz;
  double *d = tmc_malloc((size_t)L * sizeof(double));
  for (int g = 0; g < variables->n_groups; ++g)
  {
    const double *restrict y = variables->data + (size_t)g * (size_t)n_vars * L;
    const double *restrict rc = rate_constants->data + (size_t)g * (size_t)n_rxn * L;
    double *restrict J = jacobian->data + (size_t)g * (size_t)nnz * L;
    const int *dep = s->job_dependent_ids;
    const double *yield = s->job_yields;
    const int *flat = s->job_flat_ids;
    for (int q = 0; q < s->n_jobs; ++q)
    {
      const double *k = rc + (size_t)s->job_process[q] * L;
      for (int c = 0; c < L; ++c)
        d[c] = k[c];
      const int nd = s->job_n_dependent[q], np = s->job_n_products[q];
      for (int i = 0; i < nd; ++i)
      {
        const double *yi = y + (size_t)dep[i] * L;
#pragma omp simd
        for (int c = 0; c < L; ++c)
          d[c] *= yi[c];
      }
      for (int i = 0; i < nd; ++i, ++flat)
      {
        double *t = J + (size_t)*flat * L;
#pragma omp simd
        for (int c = 0; c < L; ++c)
          t[c] += d[c];
      }
      {
        double *t = J + (size_t)*flat * L;
#pragma omp simd
        for (int c = 0; c < L; ++c)
          t[c] += d[c];
        ++flat;
      }
      for (int i = 0; i < np; ++i, ++flat)
      {
        double *t = J + (size_t)*flat * L;
        const double yv = yield[i];
#pragma omp simd
        for (int c = 0; c < L; ++c)
          t[c] -= yv * d[c];
      }
      dep += nd;
      yield += np;
    }
  }
  free(d);
}

void tmc_solver_add_forcing_terms(const tmc_solver *s, const tmc_state *st, const tmc_dense *variables, tmc_dense *forcing)
{
  tmc_generic_forcing(s, st->rate_constants, variables, forcing);
}

void tmc_solver_subtract_jacobian_terms(const tmc_solver *s, const tmc_state *st, const tmc_dense *variables, tmc_sparse *jacobian)
{
  tmc_generic_jacobian(s, st->rate_constants, variables, jacobian);
}
