/* Constraint formulas and weighted Newton projection ported from MICM.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

static int add_constraint(ch_model *m, const char *name, int row, size_t nr, const int *r,
                          const double *cr, size_t np, const int *p, const double *cp,
                          double constant, double enthalpy, int kind) {
  if (!m || !name || row < 0 || (size_t)row >= m->ns || m->species[row].third_body ||
      (nr && (!r || !cr)) || (np && (!p || !cp))) return CH_INVALID;
  if (m->finalized) return CH_FROZEN;
  if (kind && (!nr || !np || constant <= 0)) return CH_INVALID;
  for (size_t i = 0; i < nr+np; ++i) {
    int id = i < nr ? r[i] : p[i-nr];
    double v = i < nr ? cr[i] : cp[i-nr];
    if (id < 0 || (size_t)id >= m->ns || m->species[id].third_body) return CH_NOT_FOUND;
    if (!isfinite(v) || (kind && v <= 0)) return CH_INVALID;
  }
  for (size_t i = 0; i < m->nconstraint; ++i) if (m->constraints[i].row == row) return CH_INVALID;
  ch_constraint c = {0};
  c.name = ch_string_copy(name); c.row = row; c.kind = kind; c.nr = nr; c.np = np;
  c.constant = constant; c.enthalpy = enthalpy; c.custom_index = -1;
  c.ids = ch_allocate(nr+np,sizeof(int)); c.flat = ch_allocate(nr+np,sizeof(int));
  c.coefficients = ch_allocate(nr+np,sizeof(double));
  if (!c.name || !c.ids || !c.flat || !c.coefficients) goto fail;
  for (size_t i = 0; i < nr+np; ++i) {
    c.ids[i] = i < nr ? r[i] : p[i-nr];
    c.coefficients[i] = i < nr ? cr[i] : cp[i-nr];
  }
  if (kind) {
    c.custom_index = ch_model_add_parameter(m,name,constant);
    if (c.custom_index < 0) goto fail;
  }
  ch_constraint *list = realloc(m->constraints,(m->nconstraint+1)*sizeof(c));
  if (!list) goto fail;
  m->constraints = list; m->constraints[m->nconstraint++] = c; return CH_OK;
fail:
  free(c.name); free(c.ids); free(c.flat); free(c.coefficients); return CH_MEMORY;
}
int ch_model_add_linear(ch_model *m, const char *name, int row, size_t n, const int *ids,
                        const double *coeff, double constant) {
  return add_constraint(m,name,row,n,ids,coeff,0,NULL,NULL,constant,0,0);
}
int ch_model_add_equilibrium(ch_model *m, const char *name, int row, size_t nr, const int *r,
                             const double *cr, size_t np, const int *p, const double *cp,
                             double constant, double enthalpy) {
  return add_constraint(m,name,row,nr,r,cr,np,p,cp,constant,enthalpy,1);
}
size_t ch_model_constraint_count(const ch_model *m) { return m ? m->nalgebraic : 0; }
int ch_model_is_algebraic(const ch_model *m, size_t i) { return m && m->finalized && i < m->nv ? m->algebraic[i] : CH_INVALID; }
const ch_model *ch_state_model(const ch_state *s) { return s ? s->model : NULL; }

int ch_model_add_external(ch_model *m, const ch_external *e) {
  if (!m || !e || (e->elements && (!e->rows || !e->columns)) ||
      (e->algebraic_count && !e->algebraic_species)) return CH_INVALID;
  if (m->finalized) return CH_FROZEN;
  for (size_t i = 0; i < e->elements; ++i) {
    if (e->rows[i] < 0 || e->columns[i] < 0 || (size_t)e->rows[i] >= m->ns ||
        (size_t)e->columns[i] >= m->ns || m->species[e->rows[i]].third_body ||
        m->species[e->columns[i]].third_body) return CH_INVALID;
  }
  for (size_t i = 0; i < e->algebraic_count; ++i)
    if (e->algebraic_species[i] < 0 || (size_t)e->algebraic_species[i] >= m->ns ||
        m->species[e->algebraic_species[i]].third_body) return CH_INVALID;
  int *rows = ch_allocate(e->elements,sizeof(int)), *cols = ch_allocate(e->elements,sizeof(int));
  int *alg = ch_allocate(e->algebraic_count,sizeof(int));
  if (!rows || !cols || !alg) { free(rows); free(cols); free(alg); return CH_MEMORY; }
  if (e->elements) { memcpy(rows,e->rows,e->elements*sizeof(int)); memcpy(cols,e->columns,e->elements*sizeof(int)); }
  if (e->algebraic_count) memcpy(alg,e->algebraic_species,e->algebraic_count*sizeof(int));
  ch_external *list = realloc(m->externals,(m->nexternal+1)*sizeof(*e));
  if (!list) { free(rows); free(cols); free(alg); return CH_MEMORY; }
  m->externals = list; ch_external *dest = &m->externals[m->nexternal++];
  *dest = *e; dest->rows = rows; dest->columns = cols; dest->algebraic_species = alg;
  return CH_OK;
}
void ch_update_constraints(ch_state *s) {
  for (size_t i = 0; i < s->model->nconstraint; ++i) {
    const ch_constraint *c = &s->model->constraints[i];
    if (!c->kind) continue;
    for (size_t k = 0; k < s->cells; ++k)
      s->custom[(size_t)c->custom_index*s->cells+k] = c->constant *
        exp((c->enthalpy/CH_GAS_CONSTANT)*(1.0/s->conditions[k].temperature-1.0/298.15));
  }
}
void ch_constraint_force(const ch_state *s, const double *y, double *f) {
  size_t nc = s->cells;
  for (size_t i = 0; i < s->model->nconstraint; ++i) {
    const ch_constraint *c = &s->model->constraints[i];
    for (size_t k = 0; k < nc; ++k) {
      if (!c->kind) {
        double value = -c->constant;
        for (size_t j = 0; j < c->nr; ++j) value += c->coefficients[j]*y[(size_t)c->ids[j]*nc+k];
        f[(size_t)c->row*nc+k] = value;
      } else {
        double left = s->custom[(size_t)c->custom_index*nc+k], right = 1.0;
        for (size_t j = 0; j < c->nr; ++j) left *= pow(fmax(y[(size_t)c->ids[j]*nc+k],0.0),c->coefficients[j]);
        for (size_t j = c->nr; j < c->nr+c->np; ++j) right *= pow(fmax(y[(size_t)c->ids[j]*nc+k],0.0),c->coefficients[j]);
        f[(size_t)c->row*nc+k] = left-right;
      }
    }
  }
  for (size_t i = 0; i < s->model->nexternal; ++i) {
    const ch_external *e = &s->model->externals[i];
    if (e->algebraic_count && e->constraint_force) e->constraint_force(s,y,f,e->context);
  }
}
void ch_constraint_jacobian(const ch_state *s, const double *y, double *a) {
  size_t nc = s->cells;
  for (size_t i = 0; i < s->model->nconstraint; ++i) {
    const ch_constraint *c = &s->model->constraints[i];
    for (size_t j = 0; j < c->nr+c->np; ++j) {
      double coeff = c->coefficients[j];
      for (size_t k = 0; k < nc; ++k) {
        double partial = coeff;
        if (c->kind) {
          size_t first = j < c->nr ? 0 : c->nr, last = j < c->nr ? c->nr : c->nr+c->np;
          partial = j < c->nr ? s->custom[(size_t)c->custom_index*nc+k] : -1.0;
          for (size_t q = first; q < last; ++q) if (q != j)
            partial *= pow(fmax(y[(size_t)c->ids[q]*nc+k],0.0),c->coefficients[q]);
          if (coeff != 1.0) {
            double conc = y[(size_t)c->ids[j]*nc+k];
            partial = conc > 0 ? coeff*partial*pow(conc,coeff-1.0) : 0.0;
          }
        }
        a[(size_t)c->flat[j]*nc+k] -= partial;
      }
    }
  }
  for (size_t i = 0; i < s->model->nexternal; ++i) {
    const ch_external *e = &s->model->externals[i];
    if (e->algebraic_count && e->constraint_jacobian) e->constraint_jacobian(s,y,a,e->context);
  }
}

static double correction_norm(const ch_state *s, const double *y, const double *delta) {
  double norm = 0;
  for (size_t i = 0; i < s->model->nv; ++i) if (s->model->algebraic[i]) {
    for (size_t c = 0; c < s->cells; ++c) {
      size_t k = i*s->cells+c;
      double scale = s->atol[i]+s->rtol*fmax(fabs(y[k]),fabs(y[k]+delta[k]));
      double q = scale > 0 ? fabs(delta[k])/scale : (delta[k] == 0 ? 0 : INFINITY);
      if (q > norm) norm = q;
    }
  }
  return norm;
}
static int check_algebraic(const ch_state *s, const double *values, double *maximum) {
  *maximum = 0;
  int status = CH_OK;
  for (size_t i = 0; i < s->model->nv; ++i) if (s->model->algebraic[i]) {
    for (size_t c = 0; c < s->cells; ++c) {
      double v = values[i*s->cells+c];
      if (isnan(v)) return CH_NAN;
      if (isinf(v)) status = CH_INF;
      if (fabs(v) > *maximum) *maximum = fabs(v);
    }
  }
  return status;
}
int ch_initialize_constraints(ch_state *s, const ch_rosenbrock_parameters *p, ch_stats *stats) {
  if (!s || !p || !stats) return CH_INVALID;
  const ch_model *m = s->model;
  if (!m->nalgebraic) return CH_OK;
  size_t n = m->nv*s->cells;
  double *delta = s->forcing, *original = s->stages, maximum;
  int modified = 0, status = CH_CONSTRAINT_INIT;
  for (int update = 0; update <= p->init_max_iterations; ++update) {
    memset(delta,0,n*sizeof(double)); ch_constraint_force(s,s->y,delta);
    status = check_algebraic(s,delta,&maximum); if (status) goto fail;
    ++stats->constraint_init_iterations;
    if (maximum == 0) return CH_OK;
    memset(s->jacobian,0,m->sparse->nnz*s->cells*sizeof(double));
    ch_constraint_jacobian(s,s->y,s->jacobian); ++stats->jacobian_updates;
    for (size_t i = 0; i < m->nv; ++i) if (!m->algebraic[i])
      for (size_t c = 0; c < s->cells; ++c) s->jacobian[(size_t)m->sparse->diagonal[i]*s->cells+c] += 1.0;
    status = ch_factor_state(s); if (status) goto fail;
    ++stats->decompositions;
    ch_linear_state(s,delta); ++stats->solves;
    status = check_algebraic(s,delta,&maximum); if (status) goto fail;
    double current = correction_norm(s,s->y,delta);
    if (current <= p->init_tolerance) return CH_OK;
    if (update == p->init_max_iterations) break;
    if (!modified) { memcpy(original,s->y,n*sizeof(double)); modified = 1; }
    if (p->init_max_backtracks > 0) {
      double step = 1.0;
      int accepted = 0;
      for (int back = 0; back <= p->init_max_backtracks; ++back, step *= p->init_backtrack_factor) {
        memcpy(s->new_y,s->y,n*sizeof(double));
        for (size_t i = 0; i < m->nv; ++i) if (m->algebraic[i])
          for (size_t c = 0; c < s->cells; ++c) s->new_y[i*s->cells+c] += step*delta[i*s->cells+c];
        if (check_algebraic(s,s->new_y,&maximum)) continue;
        memset(s->error,0,n*sizeof(double)); ch_constraint_force(s,s->new_y,s->error);
        if (check_algebraic(s,s->error,&maximum)) continue;
        if (maximum == 0) { accepted = 1; break; }
        ch_linear_state(s,s->error); ++stats->solves;
        if (check_algebraic(s,s->error,&maximum)) continue;
        double next = correction_norm(s,s->new_y,s->error);
        if (next <= p->init_tolerance || next < (1.0-p->init_sufficient_decrease*step)*current) { accepted = 1; break; }
      }
      if (!accepted) break;
      memcpy(s->y,s->new_y,n*sizeof(double));
    } else {
      for (size_t i = 0; i < m->nv; ++i) if (m->algebraic[i])
        for (size_t c = 0; c < s->cells; ++c) s->y[i*s->cells+c] += delta[i*s->cells+c];
    }
  }
  status = CH_CONSTRAINT_INIT;
fail:
  if (modified) memcpy(s->y,original,n*sizeof(double));
  return status;
}
void ch_post_solve(ch_state *s) {
  for (size_t i = 0; i < s->model->nv; ++i) if (!s->model->algebraic[i]) {
    for (size_t c = 0; c < s->cells; ++c) {
      double *v = &s->y[i*s->cells+c];
      if (*v < 0) *v = 0;
    }
  }
}
