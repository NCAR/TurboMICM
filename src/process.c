/* Translated operation ordering: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

void chem_generic_forcing(const ChemModel *m, size_t nc, size_t s, const double *y, const double *k, double *f) {
  for (size_t r = 0; r < m->nreactions; ++r) {
    const ChemReaction *rx = &m->reactions[r];
    for (size_t c = 0; c < nc; ++c) {
      double rate = k[r*s+c];
      for (size_t j = 0; j < rx->nreactants; ++j) {
        size_t v = m->species[rx->reactants[j]].variable;
        if (v != SIZE_MAX) rate *= y[v*s+c];
      }
      for (size_t j = 0; j < rx->nreactants; ++j) {
        size_t v = m->species[rx->reactants[j]].variable;
        if (v != SIZE_MAX && !m->algebraic[v]) f[v*s+c] -= rate;
      }
      for (size_t j = 0; j < rx->nproducts; ++j) {
        size_t v = m->species[rx->products[j]].variable;
        if (v != SIZE_MAX && !m->algebraic[v]) f[v*s+c] += rx->yields[j]*rate;
      }
    }
  }
}
void chem_generic_jacobian(const ChemModel *m, size_t nc, size_t s, const double *y, const double *k, double *jac) {
  size_t n = m->nvariables;
  for (size_t z = 0; z < m->njobs; ++z) {
    const ChemJacJob *job = &m->jobs[z];
    const ChemReaction *rx = &m->reactions[job->reaction];
    for (size_t c = 0; c < nc; ++c) {
      double rate = k[job->reaction*s+c];
      for (size_t j = 0; j < rx->nreactants; ++j) if (j != job->omitted) {
        size_t v = m->species[rx->reactants[j]].variable;
        if (v != SIZE_MAX) rate *= y[v*s+c];
      }
      /* Reference emits dependent reactants first, then the independent reactant. */
      for (size_t j = 0; j < rx->nreactants; ++j) if (j != job->omitted) {
        size_t v = m->species[rx->reactants[j]].variable;
        if (v != SIZE_MAX && !m->algebraic[v]) jac[(size_t)m->index[v*n+job->column]*s+c] += rate;
      }
      if (!m->algebraic[job->column]) jac[m->diagonal[job->column]*s+c] += rate;
      for (size_t j = 0; j < rx->nproducts; ++j) {
        size_t v = m->species[rx->products[j]].variable;
        if (v != SIZE_MAX && !m->algebraic[v]) jac[(size_t)m->index[v*n+job->column]*s+c] -= rx->yields[j]*rate;
      }
    }
  }
}
void chem_forcing(const ChemState *s, const double *y, double *f) {
  const ChemModel *m = s->model;
  memset(f, 0, m->nvariables*s->stride*sizeof(double));
  m->force_kernel(m, s->cells, s->stride, y, s->rates, f);
  if (m->external_force) m->external_force(m, s, y, f, m->external_context);
  if (m->nconstraints) chem_constraints(s, y, f, NULL);
}
void chem_jacobian(const ChemState *s, const double *y, double *j) {
  const ChemModel *m = s->model;
  memset(j, 0, m->nnz*s->stride*sizeof(double));
  m->jacobian_kernel(m, s->cells, s->stride, y, s->rates, j);
  if (m->external_jacobian) m->external_jacobian(m, s, y, j, m->external_context);
  if (m->nconstraints) chem_constraints(s, y, NULL, j);
}

void chem_constraints(const ChemState *s, const double *y, double *f, double *jac) {
  const ChemModel *m = s->model;
  size_t n = m->nvariables, stride = s->stride;
  for (size_t z = 0; z < m->nconstraints; ++z) {
    const ChemConstraint *q = &m->constraints[z];
    for (size_t c = 0; c < s->cells; ++c) {
      double value = 0.0;
      if (q->type == CHEM_LINEAR) {
        for (size_t i = 0; i < q->nterms; ++i) {
          value += q->coefficients[i]*y[q->species[i]*stride+c];
          if (jac) jac[(size_t)m->index[q->algebraic*n+q->species[i]]*stride+c] -= q->coefficients[i];
        }
        value -= q->constant;
      } else if (q->type == CHEM_EQUILIBRIUM) {
        double rp = q->constant, pp = 1.0;
        for (size_t i = 0; i < q->nterms; ++i) {
          double v = y[q->species[i]*stride+c];
          double factor = pow(v > 0.0 ? v : 0.0, q->coefficients[i]);
          if (i < q->nreactants) rp *= factor; else pp *= factor;
        }
        value = rp-pp;
        if (jac) for (size_t i = 0; i < q->nterms; ++i) {
          double prod = i < q->nreactants ? q->constant : 1.0;
          size_t begin = i < q->nreactants ? 0 : q->nreactants;
          size_t end = i < q->nreactants ? q->nreactants : q->nterms;
          for (size_t j = begin; j < end; ++j) if (j != i) {
            double v = y[q->species[j]*stride+c];
            prod *= pow(v > 0.0 ? v : 0.0, q->coefficients[j]);
          }
          double v = y[q->species[i]*stride+c];
          double partial = v > 0.0 ? q->coefficients[i]*prod*pow(v, q->coefficients[i]-1.0)
            : (q->coefficients[i] == 1.0 ? prod : 0.0);
          size_t at = (size_t)m->index[q->algebraic*n+q->species[i]]*stride+c;
          if (i < q->nreactants) jac[at] -= partial; else jac[at] += partial;
        }
      } else {
        double *v = malloc((n+m->ncustom+n)*sizeof(double));
        if (!v) { if (f) f[q->algebraic*stride+c] = NAN; continue; }
        double *u = v+n, *gradient = u+m->ncustom;
        for (size_t i = 0; i < n; ++i) v[i] = y[i*stride+c];
        for (size_t i = 0; i < m->ncustom; ++i) u[i] = s->custom[i*stride+c];
        if (q->residual) value = q->residual(v, u, q->context);
        if (jac && q->jacobian) {
          memset(gradient, 0, n*sizeof(double)); q->jacobian(v, u, gradient, q->context);
          for (size_t i = 0; i < q->nterms; ++i) {
            size_t col = q->species[i]; jac[(size_t)m->index[q->algebraic*n+col]*stride+c] -= gradient[col];
          }
        }
        free(v);
      }
      if (f) f[q->algebraic*stride+c] = value;
    }
  }
}
