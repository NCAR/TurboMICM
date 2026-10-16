/* SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <limits.h>

void *ch_allocate(size_t n, size_t size) {
  if (size && n > SIZE_MAX / size) return NULL;
  return calloc(n ? n : 1, size);
}
char *ch_string_copy(const char *s) {
  if (!s) return NULL;
  size_t n = strlen(s)+1;
  char *copy = malloc(n);
  if (copy) memcpy(copy,s,n);
  return copy;
}
const char *ch_status_string(int status) {
  switch (status) {
    case CH_OK: return "Converged";
    case CH_INVALID: return "Invalid argument";
    case CH_MEMORY: return "Allocation failed";
    case CH_NOT_FOUND: return "Unknown name";
    case CH_SIZE: return "Dimension mismatch";
    case CH_FROZEN: return "Model already finalized";
    case CH_SINGULAR: return "Singular matrix";
    case CH_MAX_STEPS: return "Convergence exceeded maximum steps";
    case CH_SMALL_STEP: return "Step size too small";
    case CH_NAN: return "NaN detected";
    case CH_INF: return "Inf detected";
    case CH_CONSTRAINT_INIT: return "Constraint Initialization Failed";
    default: return "Unknown status";
  }
}

ch_model *ch_model_create(void) { return calloc(1,sizeof(ch_model)); }
void ch_model_free(ch_model *m) {
  if (!m) return;
  for (size_t i = 0; i < m->ns; ++i) free(m->species[i].name);
  for (size_t i = 0; i < m->nparam; ++i) free(m->parameters[i].name);
  for (size_t i = 0; i < m->nr; ++i) {
    free(m->reactions[i].reactants); free(m->reactions[i].products); free(m->reactions[i].yields);
  }
  for (size_t i = 0; i < m->nderiv; ++i) {
    free(m->derivatives[i].reactants); free(m->derivatives[i].destinations); free(m->derivatives[i].yields);
  }
  for (size_t i = 0; i < m->nconstraint; ++i) {
    free(m->constraints[i].name); free(m->constraints[i].ids);
    free(m->constraints[i].flat); free(m->constraints[i].coefficients);
  }
  for (size_t i = 0; i < m->nexternal; ++i) {
    free((void *)m->externals[i].rows); free((void *)m->externals[i].columns);
    free((void *)m->externals[i].algebraic_species);
  }
  free(m->constraints); free(m->externals);
  free(m->species); free(m->parameters); free(m->reactions); free(m->derivatives);
  free(m->order); free(m->algebraic); ch_sparse_free(m->sparse); free(m);
}

int ch_model_add_species(ch_model *m, const char *name, double atol, int third_body) {
  if (!m || !name || !name[0] || !isfinite(atol) || atol < 0 || m->ns >= INT_MAX) return CH_INVALID;
  if (m->finalized) return CH_FROZEN;
  for (size_t i = 0; i < m->ns; ++i) if (!strcmp(m->species[i].name,name)) return CH_INVALID;
  ch_species *s = realloc(m->species,(m->ns+1)*sizeof(*s));
  if (!s) return CH_MEMORY;
  m->species = s;
  s = &m->species[m->ns];
  *s = (ch_species){0};
  s->name = ch_string_copy(name);
  if (!s->name) return CH_MEMORY;
  s->atol = atol; s->third_body = third_body; s->index = -1;
  return (int)m->ns++;
}
int ch_model_set_species_property(ch_model *m, int species, const char *key, double value) {
  if (!m || species < 0 || (size_t)species >= m->ns || !key) return CH_INVALID;
  if (m->finalized) return CH_FROZEN;
  ch_species *s = &m->species[species];
  if (!strcmp(key,"molecular weight [kg mol-1]")) { s->molecular_weight = value; s->has_weight = 1; }
  else if (!strcmp(key,"diffusion coefficient [m2 s-1]")) { s->diffusion = value; s->has_diffusion = 1; }
  else if (!strcmp(key,"absolute tolerance")) s->atol = value;
  else return CH_NOT_FOUND;
  return CH_OK;
}
int ch_model_add_parameter(ch_model *m, const char *name, double initial) {
  if (!m || !name || !name[0] || m->nparam >= INT_MAX) return CH_INVALID;
  if (m->finalized) return CH_FROZEN;
  for (size_t i = 0; i < m->nparam; ++i) if (!strcmp(m->parameters[i].name,name)) return (int)i;
  ch_parameter *p = realloc(m->parameters,(m->nparam+1)*sizeof(*p));
  if (!p) return CH_MEMORY;
  m->parameters = p;
  p[m->nparam] = (ch_parameter){ch_string_copy(name),initial};
  if (!p[m->nparam].name) return CH_MEMORY;
  return (int)m->nparam++;
}
int ch_model_add_reaction(ch_model *m, size_t nr, const int *reactants, size_t np,
                          const int *products, const double *yields, const ch_rate *rate) {
  if (!m || !rate || (nr && !reactants) || (np && (!products || !yields))) return CH_INVALID;
  if (m->finalized) return CH_FROZEN;
  for (size_t i = 0; i < nr; ++i) if (reactants[i] < 0 || (size_t)reactants[i] >= m->ns) return CH_NOT_FOUND;
  for (size_t i = 0; i < np; ++i) if (products[i] < 0 || (size_t)products[i] >= m->ns) return CH_NOT_FOUND;
  if (rate->type < CH_ARRHENIUS || rate->type > CH_LAMBDA ||
      (rate->type == CH_TAYLOR && (rate->coefficient_count < 1 || rate->coefficient_count > 16))) return CH_INVALID;
  if ((rate->type == CH_USER || rate->type == CH_SURFACE) &&
      (rate->custom_index < 0 || (size_t)rate->custom_index + (rate->type == CH_SURFACE) >= m->nparam)) return CH_NOT_FOUND;
  ch_reaction r = {0};
  r.rate = *rate;
  if (rate->type == CH_SURFACE) {
    if (nr != 1) return CH_INVALID;
    ch_species *s = &m->species[reactants[0]];
    if (!s->has_weight || !s->has_diffusion) return CH_NOT_FOUND;
    r.rate.p[0] = s->diffusion;
    r.rate.p[1] = 8.0 * CH_GAS_CONSTANT / (3.14159265358979323846 * s->molecular_weight);
  }
  ch_rate_prepare(&r.rate);
  r.nr = nr; r.np = np;
  r.reactants = ch_allocate(nr,sizeof(int)); r.products = ch_allocate(np,sizeof(int));
  r.yields = ch_allocate(np,sizeof(double));
  if (!r.reactants || !r.products || !r.yields) goto memory_error;
  if (nr) memcpy(r.reactants,reactants,nr*sizeof(int));
  if (np) { memcpy(r.products,products,np*sizeof(int)); memcpy(r.yields,yields,np*sizeof(double)); }
  ch_reaction *list = realloc(m->reactions,(m->nr+1)*sizeof(*list));
  if (!list) goto memory_error;
  m->reactions = list;
  m->reactions[m->nr] = r;
  return (int)m->nr++;
memory_error:
  free(r.reactants); free(r.products); free(r.yields); return CH_MEMORY;
}

int ch_model_finalize(ch_model *m, int reorder, int ignore_unused) {
  if (!m || !m->ns) return CH_INVALID;
  if (m->finalized) return CH_FROZEN;
  size_t n = 0;
  for (size_t i = 0; i < m->ns; ++i) if (!m->species[i].third_body) m->species[i].index = (int)n++;
  if (!n) return CH_INVALID;
  m->nv = n;
  unsigned char *p = ch_allocate(n*n,1), *used = ch_allocate(m->ns,1);
  int *perm = ch_allocate(n,sizeof(int)), *original = ch_allocate(n,sizeof(int));
  if (!p || !used || !perm || !original) { free(p); free(used); free(perm); free(original); return CH_MEMORY; }
  int status = CH_OK;
  for (size_t i = 0; i < m->ns; ++i) if (m->species[i].index >= 0) original[m->species[i].index] = (int)i;
  for (size_t ri = 0; ri < m->nr; ++ri) {
    ch_reaction *r = &m->reactions[ri];
    for (size_t j = 0; j < r->nr; ++j) used[r->reactants[j]] = 1;
    for (size_t j = 0; j < r->np; ++j) used[r->products[j]] = 1;
    for (size_t j = 0; j < r->nr; ++j) {
      int col = m->species[r->reactants[j]].index;
      if (col < 0) continue;
      for (size_t k = 0; k < r->nr; ++k) {
        int row = m->species[r->reactants[k]].index;
        if (row >= 0) p[(size_t)row*n+col] = 1;
      }
      for (size_t k = 0; k < r->np; ++k) {
        int row = m->species[r->products[k]].index;
        if (row >= 0) p[(size_t)row*n+col] = 1;
      }
    }
  }
  for (size_t i = 0; i < m->nconstraint; ++i) {
    ch_constraint *c = &m->constraints[i];
    used[c->row] = 1;
    for (size_t j = 0; j < c->nr+c->np; ++j) used[c->ids[j]] = 1;
  }
  for (size_t i = 0; i < m->nexternal; ++i) {
    ch_external *e = &m->externals[i];
    for (size_t j = 0; j < e->elements; ++j) {
      int row = m->species[e->rows[j]].index, col = m->species[e->columns[j]].index;
      p[(size_t)row*n+col] = 1; used[e->rows[j]] = used[e->columns[j]] = 1;
    }
    for (size_t j = 0; j < e->algebraic_count; ++j) used[e->algebraic_species[j]] = 1;
  }
  if (!ignore_unused) for (size_t i = 0; i < m->ns; ++i) if (!used[i]) { status = CH_NOT_FOUND; goto done; }
  if (reorder) { status = ch_markowitz(n,p,perm); if (status) goto done; }
  else for (size_t i = 0; i < n; ++i) perm[i] = (int)i;
  m->order = ch_allocate(n,sizeof(int)); m->algebraic = ch_allocate(n,1);
  if (!m->order || !m->algebraic) { status = CH_MEMORY; goto done; }
  for (size_t i = 0; i < n; ++i) {
    m->order[i] = original[perm[i]];
    m->species[m->order[i]].index = (int)i;
  }
  memset(p,0,n*n);
  size_t jobs = 0;
  for (size_t ri = 0; ri < m->nr; ++ri) {
    ch_reaction *r = &m->reactions[ri];
    size_t nr = 0, np = 0;
    for (size_t j = 0; j < r->nr; ++j) {
      int id = m->species[r->reactants[j]].index;
      if (id >= 0) r->reactants[nr++] = id;
      else ++r->nparameterized;
    }
    for (size_t j = 0; j < r->np; ++j) {
      int id = m->species[r->products[j]].index;
      if (id >= 0) { r->products[np] = id; r->yields[np++] = r->yields[j]; }
    }
    r->nr = nr; r->np = np; jobs += nr;
    for (size_t j = 0; j < nr; ++j) {
      int col = r->reactants[j];
      for (size_t k = 0; k < nr; ++k) p[(size_t)r->reactants[k]*n+col] = 1;
      for (size_t k = 0; k < np; ++k) p[(size_t)r->products[k]*n+col] = 1;
    }
  }
  for (size_t i = 0; i < m->nconstraint; ++i) {
    ch_constraint *c = &m->constraints[i];
    c->row = m->species[c->row].index;
    m->algebraic[c->row] = 1;
    for (size_t j = 0; j < c->nr+c->np; ++j) {
      c->ids[j] = m->species[c->ids[j]].index;
      p[(size_t)c->row*n+c->ids[j]] = 1;
    }
  }
  for (size_t i = 0; i < m->nexternal; ++i) {
    ch_external *e = &m->externals[i];
    for (size_t j = 0; j < e->elements; ++j) {
      int row = m->species[e->rows[j]].index, col = m->species[e->columns[j]].index;
      p[(size_t)row*n+col] = 1;
    }
    for (size_t j = 0; j < e->algebraic_count; ++j) m->algebraic[m->species[e->algebraic_species[j]].index] = 1;
  }
  for (size_t i = 0; i < n; ++i) m->nalgebraic += m->algebraic[i] != 0;
  m->sparse = ch_sparse_create(n,p);
  m->derivatives = ch_allocate(jobs,sizeof(ch_derivative));
  if (!m->sparse || !m->derivatives) { status = CH_MEMORY; goto done; }
  for (size_t i = 0; i < m->nconstraint; ++i) {
    ch_constraint *c = &m->constraints[i];
    for (size_t j = 0; j < c->nr+c->np; ++j) c->flat[j] = m->sparse->index[(size_t)c->row*n+c->ids[j]];
  }
  for (size_t col = 0; col < n; ++col) for (size_t ri = 0; ri < m->nr; ++ri) {
    ch_reaction *r = &m->reactions[ri];
    for (size_t which = 0; which < r->nr; ++which) if (r->reactants[which] == (int)col) {
      ch_derivative *d = &m->derivatives[m->nderiv++];
      d->reaction = (int)ri; d->independent = (int)col; d->nr = r->nr-1; d->np = r->np;
      d->reactants = ch_allocate(d->nr,sizeof(int));
      d->destinations = ch_allocate(r->nr+r->np,sizeof(int));
      d->yields = ch_allocate(d->np,sizeof(double));
      if (!d->reactants || !d->destinations || !d->yields) { status = CH_MEMORY; goto done; }
      size_t dep = 0; int skipped = 0;
      for (size_t j = 0; j < r->nr; ++j) {
        if (r->reactants[j] == (int)col && !skipped) { skipped = 1; continue; }
        d->reactants[dep] = r->reactants[j];
        d->destinations[dep++] = m->sparse->index[(size_t)r->reactants[j]*n+col];
      }
      d->destinations[d->nr] = m->sparse->index[col*n+col];
      for (size_t j = 0; j < r->np; ++j) {
        d->destinations[d->nr+1+j] = m->sparse->index[(size_t)r->products[j]*n+col];
        d->yields[j] = r->yields[j];
      }
    }
  }
  m->finalized = 1;
  for (size_t i = 0; i < m->nexternal; ++i) if (m->externals[i].bind) m->externals[i].bind(m,m->externals[i].context);
done:
  free(p); free(used); free(perm); free(original); return status;
}

size_t ch_model_species_count(const ch_model *m) { return m ? m->nv : 0; }
size_t ch_model_reaction_count(const ch_model *m) { return m ? m->nr : 0; }
size_t ch_model_parameter_count(const ch_model *m) { return m ? m->nparam : 0; }
const char *ch_model_species_name(const ch_model *m, size_t i) { return m && i < m->nv && m->order ? m->species[m->order[i]].name : NULL; }
const char *ch_model_parameter_name(const ch_model *m, size_t i) { return m && i < m->nparam ? m->parameters[i].name : NULL; }
int ch_model_species_index(const ch_model *m, const char *name) {
  if (!m || !name) return CH_INVALID;
  for (size_t i = 0; i < m->ns; ++i) if (!strcmp(name,m->species[i].name)) return m->species[i].index;
  return CH_NOT_FOUND;
}
int ch_model_parameter_index(const ch_model *m, const char *name) {
  if (!m || !name) return CH_INVALID;
  for (size_t i = 0; i < m->nparam; ++i) if (!strcmp(name,m->parameters[i].name)) return (int)i;
  return CH_NOT_FOUND;
}
ch_rate *ch_model_rate(ch_model *m, size_t i) { return m && i < m->nr ? &m->reactions[i].rate : NULL; }
size_t ch_model_jacobian_size(const ch_model *m) { return m ? ch_sparse_size(m->sparse) : 0; }
int ch_model_jacobian_index(const ch_model *m, size_t r, size_t c) { return m ? ch_sparse_index(m->sparse,r,c) : CH_INVALID; }
