#include "internal.h"
#include "chem_system.h"
#include <stdio.h>

char *chem_strdup(const char *s) {
  if (!s) return NULL;
  size_t n = strlen(s) + 1;
  char *p = malloc(n);
  if (p) memcpy(p, s, n);
  return p;
}
int chem_resize(void **p, size_t n, size_t size) {
  if (size && n > SIZE_MAX / size) return 0;
  void *q = realloc(*p, n ? n * size : 1);
  if (!q) return 0;
  *p = q;
  return 1;
}
const char *chem_status_string(ChemStatus status) {
  static const char *names[] = { "Converged", "Invalid argument", "Allocation failed", "Unknown name",
    "Incorrect size", "Convergence exceeded maximum steps", "Step size too small", "NaN detected",
    "Inf detected", "Constraint Initialization Failed", "Accepting unconverged integration" };
  return (unsigned)status < sizeof names / sizeof *names ? names[status] : "Unknown status";
}
ChemModel *chem_model_create(void) {
  ChemModel *m = calloc(1, sizeof *m);
  if (m) m->ignore_unused = 1;
  return m;
}
void chem_model_free(ChemModel *m) {
  if (!m) return;
  for (size_t i = 0; i < m->nspecies; ++i) chem_species_destroy(&m->species[i]);
  for (size_t i = 0; i < m->nreactions; ++i) {
    free(m->reactions[i].reactants); free(m->reactions[i].products);
    free(m->reactions[i].yields); free(m->reactions[i].rate.label);
  }
  for (size_t i = 0; i < m->nconstraints; ++i) {
    free(m->constraints[i].species); free(m->constraints[i].coefficients);
    free(m->constraints[i].name);
  }
  for (size_t i = 0; i < m->ncustom; ++i) free(m->custom_names[i]);
  free(m->species); free(m->reactions); free(m->constraints);
  free(m->custom_names); free(m->custom_defaults); free(m->variable_species);
  free(m->index); free(m->diagonal); free(m->algebraic);
  free(m->lower); free(m->upper); free(m->lower_ids); free(m->upper_ids);
  free(m->lower_columns); free(m->upper_columns); free(m->scale); free(m->scale_ids);
  free(m->updates); free(m->update_ops); free(m->jobs); free(m->declared_sparsity); free(m->external_models); free(m);
}
size_t chem_species_id(const ChemModel *m, const char *name) {
  if (m && name) for (size_t i = 0; i < m->nspecies; ++i)
    if (!strcmp(m->species[i].name, name)) return i;
  return SIZE_MAX;
}
size_t chem_variable_id(const ChemModel *m, const char *name) {
  size_t i = chem_species_id(m, name);
  return i == SIZE_MAX ? SIZE_MAX : m->species[i].variable;
}
size_t chem_custom_id(const ChemModel *m, const char *name) {
  if (m && name) for (size_t i = 0; i < m->ncustom; ++i)
    if (!strcmp(m->custom_names[i], name)) return i;
  return SIZE_MAX;
}
ChemRate *chem_model_lambda(ChemModel *m, const char *label) {
  if (!m || !label) return NULL;
  for (size_t i = 0; i < m->nreactions; ++i) {
    ChemRate *r = &m->reactions[i].rate;
    if (r->type == CHEM_LAMBDA && r->label && !strcmp(r->label,label)) return r;
  }
  return NULL;
}
size_t chem_rate_offset(const ChemModel *m, ChemRateType type) {
  size_t n = 0;
  if (m) for (size_t i = 0; i < m->nreactions; ++i) n += m->reactions[i].rate.type < type;
  return n;
}
ChemStatus chem_model_species(ChemModel *m, const char *name, double initial, int parameterized, size_t *id) {
  if (!m || !name || !*name || m->finalized || parameterized < 0 || parameterized > 2) return CHEM_INVALID;
  if (chem_species_id(m, name) != SIZE_MAX) return CHEM_INVALID;
  char *label = chem_strdup(name);
  if (!label || !chem_resize((void **)&m->species, m->nspecies + 1, sizeof *m->species)) { free(label); return CHEM_MEMORY; }
  ChemSpecies s = {0};
  s.name = label; s.initial = initial; s.parameterized = parameterized;
  s.absolute_tolerance = 1e-3; s.molecular_weight = s.diffusion = NAN;
  s.variable = parameterized ? SIZE_MAX : m->nvariables++;
  if (id) *id = m->nspecies;
  m->species[m->nspecies++] = s;
  return CHEM_OK;
}
ChemStatus chem_model_reaction(ChemModel *m, const size_t *reactants, size_t nr,
                              const size_t *products, const double *yields, size_t np, ChemRate rate) {
  if (!m || m->finalized || (nr && !reactants) || (np && (!products || !yields)) ||
      rate.type < CHEM_ARRHENIUS || rate.type > CHEM_LAMBDA || rate.ncoefficients > 16) return CHEM_INVALID;
  for (size_t j = 0; j < nr; ++j) if (reactants[j] >= m->nspecies) return CHEM_UNKNOWN;
  for (size_t j = 0; j < np; ++j) if (products[j] >= m->nspecies) return CHEM_UNKNOWN;
  if (rate.type == CHEM_SURFACE && (isnan(rate.diffusion) || isnan(rate.molecular_weight))) return CHEM_INVALID;
  ChemReaction r = {0};
  r.nreactants = nr; r.nproducts = np; r.rate = rate; r.rate.label = chem_strdup(rate.label);
  r.reactants = malloc((nr ? nr : 1) * sizeof *r.reactants);
  r.products = malloc((np ? np : 1) * sizeof *r.products);
  r.yields = malloc((np ? np : 1) * sizeof *r.yields);
  if (!r.reactants || !r.products || !r.yields || (rate.label && !r.rate.label) ||
      !chem_resize((void **)&m->reactions, m->nreactions + 1, sizeof r)) {
    free(r.reactants); free(r.products); free(r.yields); free(r.rate.label); return CHEM_MEMORY;
  }
  if (nr) memcpy(r.reactants, reactants, nr * sizeof *reactants);
  if (np) { memcpy(r.products, products, np * sizeof *products); memcpy(r.yields, yields, np * sizeof *yields); }
  m->reactions[m->nreactions++] = r;
  return CHEM_OK;
}
ChemStatus chem_model_constraint(ChemModel *m, const ChemConstraint *c) {
  if (!m || !c || m->finalized || c->algebraic >= m->nvariables || !c->nterms || !c->species || !c->coefficients)
    return CHEM_INVALID;
  if (c->type == CHEM_EQUILIBRIUM && (!(c->constant > 0.0) || !c->nreactants || c->nreactants >= c->nterms)) return CHEM_INVALID;
  for (size_t i = 0; i < c->nterms; ++i) {
    if (c->species[i] >= m->nvariables) return CHEM_UNKNOWN;
    if (c->type == CHEM_EQUILIBRIUM && c->coefficients[i] <= 0.0) return CHEM_INVALID;
  }
  for (size_t i = 0; i < m->nconstraints; ++i) if (c->algebraic == m->constraints[i].algebraic) return CHEM_INVALID;
  ChemConstraint d = *c;
  d.name = chem_strdup(c->name); d.custom_index = SIZE_MAX;
  d.species = malloc(c->nterms * sizeof *d.species);
  d.coefficients = malloc(c->nterms * sizeof *d.coefficients);
  if (!d.species || !d.coefficients || (c->name && !d.name) || !chem_resize((void **)&m->constraints, m->nconstraints + 1, sizeof d)) {
    free(d.species); free(d.coefficients); free(d.name); return CHEM_MEMORY;
  }
  memcpy(d.species, c->species, c->nterms * sizeof *d.species);
  memcpy(d.coefficients, c->coefficients, c->nterms * sizeof *d.coefficients);
  m->constraints[m->nconstraints++] = d;
  return CHEM_OK;
}
ChemStatus chem_constraint_sparsity(const ChemModel *m, unsigned char *pattern, size_t n) {
  if (!m || !pattern || !n || n != m->nvariables || n > SIZE_MAX/n) return CHEM_INVALID;
  memset(pattern,0,n*n);
  for (size_t i = 0; i < m->nconstraints; ++i) {
    const ChemConstraint *q = &m->constraints[i];
    pattern[q->algebraic*n+q->algebraic] = 1;
    for (size_t j = 0; j < q->nterms; ++j) pattern[q->algebraic*n+q->species[j]] = 1;
  }
  return CHEM_OK;
}
static ChemStatus custom_parameter(ChemModel *m, const char *label, size_t *index) {
  *index = chem_custom_id(m, label);
  if (*index != SIZE_MAX) return CHEM_OK;
  char *s = chem_strdup(label);
  if (!s || !chem_resize((void **)&m->custom_names, m->ncustom + 1, sizeof *m->custom_names) ||
      !chem_resize((void **)&m->custom_defaults, m->ncustom + 1, sizeof(double))) { free(s); return CHEM_MEMORY; }
  *index = m->ncustom; m->custom_names[m->ncustom] = s; m->custom_defaults[m->ncustom++] = 0.0;
  return CHEM_OK;
}
ChemStatus chem_model_custom(ChemModel *m, const char *label, double initial, size_t *index) {
  if (!m || m->finalized || !label || !*label || !index || !isfinite(initial)) return CHEM_INVALID;
  ChemStatus status = custom_parameter(m,label,index);
  if (!status) m->custom_defaults[*index] = initial;
  return status;
}
ChemStatus chem_model_external(ChemModel *m, ChemExternal external) {
  if (!m || m->finalized) return CHEM_INVALID;
  if (!chem_resize((void **)&m->external_models,m->nexternal_models+1,sizeof external)) return CHEM_MEMORY;
  m->external_models[m->nexternal_models++] = external;
  return CHEM_OK;
}
ChemStatus chem_model_sparsity(ChemModel *m, const unsigned char *pattern, size_t n) {
  if (!m || !pattern || m->finalized || n != m->nvariables || !n || n > SIZE_MAX/n) return CHEM_INVALID;
  unsigned char *copy = malloc(n*n); if (!copy) return CHEM_MEMORY;
  memcpy(copy,pattern,n*n); free(m->declared_sparsity); m->declared_sparsity = copy; return CHEM_OK;
}

ChemStatus chem_model_finalize(ChemModel *m) {
  if (!m || !m->nvariables || m->finalized) return CHEM_INVALID;
  if (!m->ignore_unused && !m->external_force && !m->nexternal_models) {
    for (size_t s = 0; s < m->nspecies; ++s) {
      int used = 0;
      for (size_t i = 0; i < m->nreactions; ++i) {
        const ChemReaction *r = &m->reactions[i];
        for (size_t j = 0; j < r->nreactants; ++j) used |= r->reactants[j] == s;
        for (size_t j = 0; j < r->nproducts; ++j) used |= r->products[j] == s;
      }
      for (size_t i = 0; i < m->nconstraints; ++i) {
        const ChemConstraint *q = &m->constraints[i];
        for (size_t j = 0; j < q->nterms; ++j) used |= q->species[j] == m->species[s].variable;
        used |= q->algebraic == m->species[s].variable;
      }
      if (!used) return CHEM_UNKNOWN;
    }
  }
  /* Insertion sort is stable and runs only when a mechanism is built. */
  for (size_t i = 1; i < m->nreactions; ++i) {
    ChemReaction r = m->reactions[i]; size_t j = i;
    while (j && m->reactions[j-1].rate.type > r.rate.type) { m->reactions[j] = m->reactions[j-1]; --j; }
    m->reactions[j] = r;
  }
  for (size_t r = 0; r < m->nreactions; ++r) {
    ChemRate *p = &m->reactions[r].rate;
    if (p->type == CHEM_BRANCHED) {
      p->k0 = 2.0e-22 * 6.02214076e23 * 1.0e-6 * exp(p->n);
      double a = p->k0*(2.45e19/6.02214076e23*1.0e6), b = 0.43*pow(293.0/298.0,-8.0);
      double value = a/(1.0+a/b)*pow(0.41,1.0/(1.0+pow(log10(a/b),2.0)));
      p->z = value*(1.0-p->a0)/p->a0;
    }
    if (p->type == CHEM_SURFACE) p->mean_free_speed_factor = 8.0*CHEM_GAS_CONSTANT/(CHEM_PI*p->molecular_weight);
    if (p->type == CHEM_USER || p->type == CHEM_SURFACE) {
      if (!p->label) return CHEM_INVALID;
      if (p->type == CHEM_USER) {
        ChemStatus code = custom_parameter(m, p->label, &p->custom_index); if (code) return code;
      } else {
        static const char *suffix[] = { ".effective radius [m]", ".particle number concentration [# m-3]" };
        for (size_t i = 0; i < 2; ++i) {
          size_t index, n = strlen(p->label) + strlen(suffix[i]) + 1;
          char *label = malloc(n); if (!label) return CHEM_MEMORY;
          snprintf(label, n, "%s%s", p->label, suffix[i]);
          ChemStatus code = custom_parameter(m, label, &index); free(label); if (code) return code;
          if (!i) p->custom_index = index;
          else if (index != p->custom_index + 1) return CHEM_INVALID;
        }
      }
    }
  }
  for (size_t i = 0; i < m->nconstraints; ++i) {
    ChemConstraint *q = &m->constraints[i];
    if (q->type == CHEM_EQUILIBRIUM && q->name) {
      ChemStatus code = custom_parameter(m,q->name,&q->custom_index); if (code) return code;
      m->custom_defaults[q->custom_index] = q->constant;
    }
  }
  const size_t n = m->nvariables;
  if (n > SIZE_MAX / n || n * n > SIZE_MAX / sizeof(int)) return CHEM_MEMORY;
  unsigned char *pattern = calloc(n*n, 1);
  m->algebraic = calloc(n, 1); m->variable_species = malloc(n * sizeof(size_t));
  m->index = malloc(n*n * sizeof(int)); m->diagonal = malloc(n * sizeof(size_t));
  m->lower = calloc(n, sizeof(ChemRange)); m->upper = calloc(n, sizeof(ChemRange));
  m->scale = calloc(n, sizeof(ChemRange)); m->updates = calloc(n, sizeof(ChemRange));
  if (!pattern || !m->algebraic || !m->variable_species || !m->index || !m->diagonal || !m->lower ||
      !m->upper || !m->scale || !m->updates) { free(pattern); return CHEM_MEMORY; }
  if (m->declared_sparsity) memcpy(pattern,m->declared_sparsity,n*n);
  for (size_t i = 0; i < m->nspecies; ++i) if (!m->species[i].parameterized) m->variable_species[m->species[i].variable] = i;
  for (size_t i = 0; i < m->nconstraints; ++i) m->algebraic[m->constraints[i].algebraic] = 1;
  for (size_t i = 0; i < n; ++i) pattern[i*n+i] = 1;
  size_t njobs = 0;
  for (size_t col = 0; col < n; ++col) for (size_t r = 0; r < m->nreactions; ++r) {
    ChemReaction *rx = &m->reactions[r];
    for (size_t k = 0; k < rx->nreactants; ++k) {
      if (m->species[rx->reactants[k]].variable != col) continue;
      if (!chem_resize((void **)&m->jobs, njobs + 1, sizeof *m->jobs)) { free(pattern); return CHEM_MEMORY; }
      m->jobs[njobs++] = (ChemJacJob){r, col, k, rx->nreactants, rx->nproducts};
      for (size_t j = 0; j < rx->nreactants; ++j) {
        size_t row = m->species[rx->reactants[j]].variable;
        if (row != SIZE_MAX && !m->algebraic[row]) pattern[row*n+col] = 1;
      }
      for (size_t j = 0; j < rx->nproducts; ++j) {
        size_t row = m->species[rx->products[j]].variable;
        if (row != SIZE_MAX && !m->algebraic[row]) pattern[row*n+col] = 1;
      }
    }
  }
  m->njobs = njobs;
  for (size_t c = 0; c < m->nconstraints; ++c) {
    ChemConstraint *q = &m->constraints[c];
    for (size_t j = 0; j < q->nterms; ++j) pattern[q->algebraic*n+q->species[j]] = 1;
  }
  /* External kernels declare conservative dense structure, including algebraic rows. */
  if ((m->external_jacobian || m->nexternal_models) && !m->declared_sparsity) memset(pattern, 1, n*n);
  for (size_t i = 0; i < n; ++i) for (size_t k = i+1; k < n; ++k) if (pattern[i*n+k])
    for (size_t j = i+1; j < n; ++j) if (pattern[j*n+i]) pattern[j*n+k] = 1;
  for (size_t i = 0; i < n*n; ++i) m->index[i] = pattern[i] ? (int)m->nnz++ : -1;
  size_t nl = 0, nu = 0, ns = 0, no = 0;
  for (size_t i = 0; i < n; ++i) {
    m->diagonal[i] = (size_t)m->index[i*n+i];
    for (size_t j = 0; j < i; ++j) if (pattern[i*n+j]) ++nl;
    for (size_t j = i+1; j < n; ++j) if (pattern[i*n+j]) ++nu;
    for (size_t j = i+1; j < n; ++j) if (pattern[j*n+i]) ++ns;
    for (size_t k = i+1; k < n; ++k) if (pattern[i*n+k])
      for (size_t j = i+1; j < n; ++j) if (pattern[j*n+i]) ++no;
  }
  m->lower_ids = malloc((nl+1)*sizeof(size_t)); m->lower_columns = malloc((nl+1)*sizeof(size_t));
  m->upper_ids = malloc((nu+1)*sizeof(size_t)); m->upper_columns = malloc((nu+1)*sizeof(size_t));
  m->scale_ids = malloc((ns+1)*sizeof(size_t)); m->update_ops = malloc((no+1)*sizeof(ChemTriple));
  if (!m->lower_ids || !m->lower_columns || !m->upper_ids || !m->upper_columns || !m->scale_ids || !m->update_ops) {
    free(pattern); return CHEM_MEMORY;
  }
  nl = nu = ns = no = 0;
  for (size_t i = 0; i < n; ++i) {
    m->lower[i].first = nl; m->upper[i].first = nu;
    m->scale[i].first = ns; m->updates[i].first = no;
    for (size_t j = 0; j < i; ++j) if (pattern[i*n+j]) {
      m->lower_ids[nl] = (size_t)m->index[i*n+j]; m->lower_columns[nl++] = j; ++m->lower[i].count;
    }
    for (size_t j = i+1; j < n; ++j) if (pattern[i*n+j]) {
      m->upper_ids[nu] = (size_t)m->index[i*n+j]; m->upper_columns[nu++] = j; ++m->upper[i].count;
    }
    for (size_t j = i+1; j < n; ++j) if (pattern[j*n+i]) {
      m->scale_ids[ns++] = (size_t)m->index[j*n+i]; ++m->scale[i].count;
    }
    for (size_t k = i+1; k < n; ++k) if (pattern[i*n+k])
      for (size_t j = i+1; j < n; ++j) if (pattern[j*n+i]) {
        m->update_ops[no++] = (ChemTriple){(size_t)m->index[j*n+k], (size_t)m->index[j*n+i], (size_t)m->index[i*n+k]};
        ++m->updates[i].count;
      }
  }
  free(pattern);
  m->force_kernel = chem_generic_forcing; m->jacobian_kernel = chem_generic_jacobian;
  m->finalized = 1;
  if (m->setup) m->setup(m,m->setup_context);
  for (size_t i = 0; i < m->nexternal_models; ++i) {
    const ChemExternal *e = m->external_models+i;
    if (e->setup) e->setup(m,e->context);
  }
  return CHEM_OK;
}
