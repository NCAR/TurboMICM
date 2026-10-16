/* SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

ch_state *ch_state_create(const ch_model *m, size_t cells) {
  if (!m || !m->finalized || !m->nv || !cells || cells > SIZE_MAX/m->nv/6) return NULL;
  if ((m->nparam && cells > SIZE_MAX/m->nparam) ||
      (m->nr && cells > SIZE_MAX/m->nr) || cells > SIZE_MAX/m->sparse->nnz) return NULL;
  ch_state *s = calloc(1,sizeof(*s));
  if (!s) return NULL;
  s->model = m; s->cells = cells; s->rtol = 1e-6;
  size_t n = m->nv*cells;
  s->y = ch_allocate(n,sizeof(double)); s->custom = ch_allocate(m->nparam*cells,sizeof(double));
  s->rates = ch_allocate(m->nr*cells,sizeof(double)); s->atol = ch_allocate(m->nv,sizeof(double));
  s->conditions = ch_allocate(cells,sizeof(ch_conditions));
  s->new_y = ch_allocate(n,sizeof(double)); s->forcing = ch_allocate(n,sizeof(double));
  s->error = ch_allocate(n,sizeof(double)); s->stages = ch_allocate(6*n,sizeof(double));
  s->inverse_diagonal = ch_allocate(n,sizeof(double));
  s->jacobian = ch_allocate(m->sparse->nnz*cells,sizeof(double)); s->scratch = ch_allocate(cells,sizeof(double));
  if (!s->y || !s->custom || !s->rates || !s->atol || !s->conditions || !s->new_y ||
      !s->forcing || !s->error || !s->stages || !s->jacobian || !s->scratch || !s->inverse_diagonal) {
    ch_state_free(s); return NULL;
  }
  for (size_t i = 0; i < m->nv; ++i) s->atol[i] = m->species[m->order[i]].atol;
  for (size_t i = 0; i < m->nparam; ++i)
    for (size_t c = 0; c < cells; ++c) s->custom[i*cells+c] = m->parameters[i].initial;
  return s;
}
ch_state *ch_state_clone(const ch_state *s) {
  if (!s) return NULL;
  ch_state *d = ch_state_create(s->model,s->cells);
  if (!d) return NULL;
  ch_state_copy(d,s);
  return d;
}
int ch_state_copy(ch_state *d, const ch_state *s) {
  if (!d || !s) return CH_INVALID;
  if (d->model != s->model || d->cells != s->cells) return CH_SIZE;
  if (d == s) return CH_OK;
  memcpy(d->y,s->y,s->model->nv*s->cells*sizeof(double));
  memcpy(d->custom,s->custom,s->model->nparam*s->cells*sizeof(double));
  memcpy(d->rates,s->rates,s->model->nr*s->cells*sizeof(double));
  memcpy(d->atol,s->atol,s->model->nv*sizeof(double));
  memcpy(d->conditions,s->conditions,s->cells*sizeof(ch_conditions)); d->rtol = s->rtol;
  d->lu_algorithm = s->lu_algorithm; d->factor = s->factor; d->linear_solve = s->linear_solve;
  d->linear_context = s->linear_context;
  return CH_OK;
}
void ch_state_move(ch_state **d, ch_state **s) {
  if (!d || !s || d == s) return;
  if (*d != *s) ch_state_free(*d);
  *d = *s; *s = NULL;
}
void ch_state_free(ch_state *s) {
  if (!s) return;
  free(s->y); free(s->custom); free(s->rates); free(s->atol); free(s->conditions);
  free(s->inverse_diagonal);
  free(s->new_y); free(s->forcing); free(s->error); free(s->stages); free(s->jacobian); free(s->scratch); free(s);
}
size_t ch_state_cells(const ch_state *s) { return s ? s->cells : 0; }
double *ch_state_values(ch_state *s) { return s ? s->y : NULL; }
double *ch_state_parameters(ch_state *s) { return s ? s->custom : NULL; }
double *ch_state_rates(ch_state *s) { return s ? s->rates : NULL; }
ch_conditions *ch_state_conditions(ch_state *s) { return s ? s->conditions : NULL; }
double ch_state_relative_tolerance(const ch_state *s) { return s ? s->rtol : NAN; }
const double *ch_state_absolute_tolerances(const ch_state *s) { return s ? s->atol : NULL; }

int ch_state_set_lu(ch_state *s, int algorithm) {
  if (!s || (algorithm != 0 && algorithm != 1)) return CH_INVALID;
  s->lu_algorithm = algorithm; s->factor = NULL; s->linear_solve = NULL;
  return CH_OK;
}
int ch_state_set_linear_solver(ch_state *s, ch_factor_function factor, ch_linear_function solve, void *context) {
  if (!s || !factor || !solve) return CH_INVALID;
  s->factor = factor; s->linear_solve = solve; s->linear_context = context; return CH_OK;
}
int ch_factor_state(ch_state *s) {
  if (s->factor) return s->factor(s->model->sparse,s->cells,s->jacobian,s->scratch,s->linear_context);
  return ch_lu_factor_work(s->model->sparse,s->cells,s->jacobian,s->lu_algorithm,s->scratch,s->inverse_diagonal);
}
void ch_linear_state(ch_state *s, double *x) {
  if (s->linear_solve) s->linear_solve(s->model->sparse,s->cells,s->jacobian,x,s->linear_context);
  else ch_lu_solve_work(s->model->sparse,s->cells,s->jacobian,x,s->inverse_diagonal);
}
void ch_alpha_minus_jacobian(const ch_state *s, double alpha, double *a) {
  for (size_t i = 0; i < s->model->nv; ++i) if (!s->model->algebraic[i]) {
    double *restrict diag = a+(size_t)s->model->sparse->diagonal[i]*s->cells;
    #pragma omp simd
    for (size_t c = 0; c < s->cells; ++c) diag[c] += alpha;
  }
}

static int set_column(ch_state *s, int id, double *dest, const double *v, size_t n) {
  if (id < 0) return id;
  if (!v) return CH_INVALID;
  if (n != 1 && n != s->cells) return CH_SIZE;
  for (size_t c = 0; c < s->cells; ++c) dest[(size_t)id*s->cells+c] = v[n == 1 ? 0 : c];
  return CH_OK;
}
int ch_state_set(ch_state *s, const char *name, const double *v, size_t n) {
  return s ? set_column(s,ch_model_species_index(s->model,name),s->y,v,n) : CH_INVALID;
}
int ch_state_set_parameter(ch_state *s, const char *name, const double *v, size_t n) {
  return s ? set_column(s,ch_model_parameter_index(s->model,name),s->custom,v,n) : CH_INVALID;
}
int ch_state_set_exact(ch_state *s, const char *name, const double *v, size_t n) {
  if (!s) return CH_INVALID;
  if (n != s->cells) return CH_SIZE;
  return ch_state_set(s,name,v,n);
}
int ch_state_set_parameter_exact(ch_state *s, const char *name, const double *v, size_t n) {
  if (!s) return CH_INVALID;
  if (n != s->cells) return CH_SIZE;
  return ch_state_set_parameter(s,name,v,n);
}
int ch_state_set_parameter_matrix(ch_state *s, const double *v, size_t rows, size_t cols) {
  if (!s || !v) return CH_INVALID;
  if (rows != s->cells || cols != s->model->nparam) return CH_SIZE;
  for (size_t c = 0; c < rows; ++c) for (size_t p = 0; p < cols; ++p) s->custom[p*rows+c] = v[c*cols+p];
  return CH_OK;
}
int ch_state_tolerances(ch_state *s, double rtol, const double *atol, size_t n) {
  if (!s || !atol || !isfinite(rtol) || rtol < 0) return CH_INVALID;
  if (n != 1 && n != s->model->nv) return CH_SIZE;
  for (size_t i = 0; i < n; ++i) if (!isfinite(atol[i]) || atol[i] < 0 || (atol[i] == 0 && rtol == 0)) return CH_INVALID;
  s->rtol = rtol;
  for (size_t i = 0; i < s->model->nv; ++i) s->atol[i] = atol[n == 1 ? 0 : i];
  return CH_OK;
}
int ch_update_rates(ch_state *s) {
  if (!s) return CH_INVALID;
  const ch_model *m = s->model;
  for (size_t i = 0; i < m->nexternal; ++i) if (m->externals[i].update)
    m->externals[i].update(s,m->externals[i].context);
  ch_update_constraints(s);
  double *custom = ch_allocate(m->nparam,sizeof(double));
  if (!custom) return CH_MEMORY;
  for (size_t c = 0; c < s->cells; ++c) {
    for (size_t i = 0; i < m->nparam; ++i) custom[i] = s->custom[i*s->cells+c];
    for (size_t ri = 0; ri < m->nr; ++ri) {
      const ch_reaction *r = &m->reactions[ri];
      double value = ch_rate_calculate(&r->rate,&s->conditions[c],custom);
      for (size_t i = 0; i < r->nparameterized; ++i)
        value *= ch_species_parameter_value(m->species[r->parameterized[i]].spec,&s->conditions[c]);
      s->rates[ri*s->cells+c] = value;
    }
  }
  free(custom); return CH_OK;
}

void ch_force(const ch_state *s, const double *y, double *f) {
  const ch_model *m = s->model;
  size_t nc = s->cells;
  double *restrict rate = s->scratch;
  memset(f,0,m->nv*nc*sizeof(double));
  for (size_t ri = 0; ri < m->nr; ++ri) {
    const ch_reaction *r = &m->reactions[ri];
    memcpy(rate,s->rates+ri*nc,nc*sizeof(double));
    for (size_t j = 0; j < r->nr; ++j) {
      const double *restrict reactant = y + (size_t)r->reactants[j]*nc;
      #pragma omp simd
      for (size_t c = 0; c < nc; ++c) rate[c] *= reactant[c];
    }
    for (size_t j = 0; j < r->nr; ++j) {
      if (m->algebraic[r->reactants[j]]) continue;
      double *restrict dst = f + (size_t)r->reactants[j]*nc;
      #pragma omp simd
      for (size_t c = 0; c < nc; ++c) dst[c] -= rate[c];
    }
    for (size_t j = 0; j < r->np; ++j) {
      if (m->algebraic[r->products[j]]) continue;
      double *restrict dst = f + (size_t)r->products[j]*nc;
      double yield = r->yields[j];
      #pragma omp simd
      for (size_t c = 0; c < nc; ++c) dst[c] += yield*rate[c];
    }
  }
  for (size_t i = 0; i < m->nexternal; ++i) if (m->externals[i].force)
    m->externals[i].force(s,y,f,m->externals[i].context);
  if (m->nexternal) for (size_t i = 0; i < m->nv; ++i) if (m->algebraic[i])
    memset(f+i*nc,0,nc*sizeof(double));
}
void ch_jacobian(const ch_state *s, const double *y, double *a) {
  const ch_model *m = s->model;
  size_t nc = s->cells;
  double *restrict rate = s->scratch;
  memset(a,0,m->sparse->nnz*nc*sizeof(double));
  for (size_t i = 0; i < m->nderiv; ++i) {
    const ch_derivative *d = &m->derivatives[i];
    memcpy(rate,s->rates+(size_t)d->reaction*nc,nc*sizeof(double));
    for (size_t j = 0; j < d->nr; ++j) {
      const double *restrict reactant = y + (size_t)d->reactants[j]*nc;
      #pragma omp simd
      for (size_t c = 0; c < nc; ++c) rate[c] *= reactant[c];
    }
    for (size_t j = 0; j <= d->nr; ++j) {
      int row = j == d->nr ? d->independent : d->reactants[j];
      if (m->algebraic[row]) continue;
      double *restrict dst = a + (size_t)d->destinations[j]*nc;
      #pragma omp simd
      for (size_t c = 0; c < nc; ++c) dst[c] += rate[c];
    }
    const ch_reaction *r = &m->reactions[d->reaction];
    for (size_t j = 0; j < d->np; ++j) {
      if (m->algebraic[r->products[j]]) continue;
      double *restrict dst = a + (size_t)d->destinations[d->nr+1+j]*nc;
      double yield = d->yields[j];
      #pragma omp simd
      for (size_t c = 0; c < nc; ++c) dst[c] -= yield*rate[c];
    }
  }
  for (size_t i = 0; i < m->nexternal; ++i) if (m->externals[i].jacobian)
    m->externals[i].jacobian(s,y,a,m->externals[i].context);
  if (m->nexternal) for (size_t i = 0; i < m->nv; ++i) if (m->algebraic[i])
    memset(a+(size_t)m->sparse->row_start[i]*nc,0,
           (size_t)(m->sparse->row_start[i+1]-m->sparse->row_start[i])*nc*sizeof(double));
}
