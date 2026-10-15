#include "internal.h"
#include <stdio.h>

void *ck_allocate(size_t n, size_t size)
{
  if (size && n > SIZE_MAX / size) return NULL;
  return calloc(n ? n : 1, size ? size : 1);
}
char *ck_string_copy(const char *s)
{
  if (!s) return NULL;
  size_t n = strlen(s) + 1;
  char *copy = ck_allocate(n, 1);
  if (copy) memcpy(copy, s, n);
  return copy;
}
static char *custom_label(const char *name, const char *suffix)
{
  if (!name) return NULL;
  size_t a = strlen(name), b = strlen(suffix);
  char *s = ck_allocate(a + b + 1, 1);
  if (s) { memcpy(s, name, a); memcpy(s + a, suffix, b + 1); }
  return s;
}
void ck_model_free(ck_model *m)
{
  if (!m) return;
  for (size_t r = 0; r < m->nr; ++r) {
    if (m->reactants) free(m->reactants[r]);
    if (m->products) free(m->products[r]);
  }
  for (size_t j = 0; j < m->njobs; ++j) {
    free(m->jobs[j].reactants); free(m->jobs[j].products); free(m->jobs[j].targets);
  }
  for (size_t c = 0; c < m->nc; ++c) if (m->custom_names) free(m->custom_names[c]);
  for (size_t k = 0; k < m->nconstraints; ++k) {
    free(m->constraints[k].reactants); free(m->constraints[k].products);
  }
  free(m->constraints); free(m->algebraic);
  free(m->species_to_variable); free(m->variable_to_species); free(m->custom_offset); free(m->reaction_order);
  free(m->nreactants); free(m->reactants); free(m->nproducts); free(m->products); free(m->rates);
  free(m->custom_names); free(m->map); free(m->diagonal); free(m->row_start); free(m->columns);
  free(m->jobs); free(m->ops); free(m->divide_start); free(m->divides); free(m->update_start); free(m);
}

static ck_status symbolic(ck_model *m, const unsigned char *supplied_pattern)
{
  size_t n = m->n;
  unsigned char *pattern = ck_allocate(n * n, 1);
  if (!pattern) return CK_MEMORY;
  if (supplied_pattern) memcpy(pattern, supplied_pattern, n * n);
  for (size_t i = 0; i < n; ++i) pattern[i * n + i] = 1;
  for (size_t r = 0; r < m->nr; ++r) {
    for (size_t k = 0; k < m->nreactants[r]; ++k) {
      size_t col = m->reactants[r][k];
      for (size_t j = 0; j < m->nreactants[r]; ++j) pattern[m->reactants[r][j] * n + col] = 1;
      for (size_t j = 0; j < m->nproducts[r]; ++j) pattern[m->products[r][j].species * n + col] = 1;
    }
  }
  for (size_t k = 0; k < m->nconstraints; ++k) {
    ck_bound_constraint *b = &m->constraints[k];
    for (size_t j = 0; j < b->description.nreactants; ++j) pattern[b->row * n + b->reactants[j].species] = 1;
    for (size_t j = 0; j < b->description.nproducts; ++j) pattern[b->row * n + b->products[j].species] = 1;
  }
  for (size_t i = 0; i < n; ++i)
    for (size_t k = i + 1; k < n; ++k) if (pattern[i * n + k])
      for (size_t j = i + 1; j < n; ++j) if (pattern[j * n + i]) pattern[j * n + k] = 1;
  for (size_t i = 0; i < n * n; ++i) if (pattern[i]) ++m->nnz;
  m->map = ck_allocate(n * n, sizeof *m->map);
  m->diagonal = ck_allocate(n, sizeof *m->diagonal);
  m->row_start = ck_allocate(n + 1, sizeof *m->row_start);
  m->columns = ck_allocate(m->nnz, sizeof *m->columns);
  m->divide_start = ck_allocate(n + 1, sizeof *m->divide_start);
  m->divides = ck_allocate(m->nnz, sizeof *m->divides);
  m->update_start = ck_allocate(n + 1, sizeof *m->update_start);
  if (!m->map || !m->diagonal || !m->row_start || !m->columns || !m->divide_start || !m->divides || !m->update_start) {
    free(pattern); return CK_MEMORY;
  }
  size_t off = 0, divides = 0;
  for (size_t i = 0; i < n; ++i) {
    m->row_start[i] = off;
    for (size_t j = 0; j < n; ++j) {
      m->map[i * n + j] = pattern[i * n + j] ? off : SIZE_MAX;
      if (pattern[i * n + j]) m->columns[off++] = j;
    }
    m->diagonal[i] = m->map[i * n + i];
  }
  m->row_start[n] = off;
  for (size_t i = 0; i < n; ++i) {
    m->divide_start[i] = divides;
    m->update_start[i] = m->nops;
    for (size_t j = i + 1; j < n; ++j) if (pattern[j * n + i]) m->divides[divides++] = m->map[j * n + i];
    for (size_t k = i + 1; k < n; ++k) if (pattern[i * n + k])
      for (size_t j = i + 1; j < n; ++j) if (pattern[j * n + i]) ++m->nops;
  }
  m->divide_start[n] = divides; m->update_start[n] = m->nops;
  m->ops = ck_allocate(m->nops, sizeof *m->ops);
  if (!m->ops) { free(pattern); return CK_MEMORY; }
  off = 0;
  for (size_t i = 0; i < n; ++i)
    for (size_t k = i + 1; k < n; ++k) if (pattern[i * n + k])
      for (size_t j = i + 1; j < n; ++j) if (pattern[j * n + i])
        m->ops[off++] = (ck_operation){m->map[j * n + k], m->map[j * n + i], m->map[i * n + k]};
  free(pattern);
  return CK_OK;
}
static ck_status derivatives(ck_model *m)
{
  for (size_t r = 0; r < m->nr; ++r) m->njobs += m->nreactants[r];
  m->jobs = ck_allocate(m->njobs, sizeof *m->jobs);
  if (!m->jobs) return CK_MEMORY;
  size_t job = 0;
  for (size_t col = 0; col < m->n; ++col) {
    for (size_t r = 0; r < m->nr; ++r) {
      for (size_t occurrence = 0; occurrence < m->nreactants[r]; ++occurrence) {
        if (m->reactants[r][occurrence] != col) continue;
        ck_derivative *d = &m->jobs[job++];
        d->reaction = r; d->independent = col;
        d->nreactants = m->nreactants[r] - 1; d->nproducts = m->nproducts[r];
        d->reactants = ck_allocate(d->nreactants, sizeof *d->reactants);
        d->products = ck_allocate(d->nproducts, sizeof *d->products);
        d->targets = ck_allocate(d->nreactants + 1 + d->nproducts, sizeof *d->targets);
        if (!d->reactants || !d->products || !d->targets) return CK_MEMORY;
        size_t k = 0; int found = 0;
        for (size_t j = 0; j < m->nreactants[r]; ++j) {
          size_t s = m->reactants[r][j];
          if (s == col && !found) { found = 1; continue; }
          d->reactants[k] = s; d->targets[k++] = m->map[s * m->n + col];
        }
        d->targets[k++] = m->map[col * m->n + col];
        for (size_t j = 0; j < d->nproducts; ++j) {
          d->products[j] = m->products[r][j];
          d->targets[k++] = m->map[d->products[j].species * m->n + col];
        }
      }
    }
  }
  return CK_OK;
}
static int rate_priority(ck_rate_type type)
{
  static const int priorities[] = {0,1,2,3,5,4,4,6,7,8,9,10};
  return (unsigned)type < sizeof priorities / sizeof *priorities ? priorities[type] : 11;
}
ck_status ck_model_create(const ck_mechanism *desc, ck_model **out)
{ return ck_model_create_with_order(desc, 1, out); }
ck_status ck_model_create_with_order(const ck_mechanism *desc, int reorder, ck_model **out)
{
  if (!out) return CK_INVALID;
  *out = NULL;
  if (!desc || !desc->nspecies || !desc->species || (desc->nreactions && !desc->reactions) || (desc->nconstraints && !desc->constraints)) return CK_INVALID;
  ck_model *m = ck_allocate(1, sizeof *m);
  if (!m) return CK_MEMORY;
  m->description = desc; m->nr = desc->nreactions;
  m->species_to_variable = ck_allocate(desc->nspecies, sizeof(size_t));
  m->variable_to_species = ck_allocate(desc->nspecies, sizeof(size_t));
  m->reaction_order = ck_allocate(m->nr, sizeof(size_t));
  m->custom_offset = ck_allocate(m->nr, sizeof(size_t));
  m->nreactants = ck_allocate(m->nr, sizeof(size_t)); m->nproducts = ck_allocate(m->nr, sizeof(size_t));
  m->reactants = ck_allocate(m->nr, sizeof(size_t *)); m->products = ck_allocate(m->nr, sizeof(ck_product *));
  m->rates = ck_allocate(m->nr, sizeof(ck_rate));
  m->custom_names = ck_allocate(2 * m->nr + desc->nconstraints, sizeof(char *));
  ck_status error = CK_MEMORY;
  if (!m->species_to_variable || !m->variable_to_species || !m->reaction_order || !m->custom_offset || !m->nreactants ||
      !m->nproducts || !m->reactants || !m->products || !m->rates || !m->custom_names) goto fail;
  error = CK_INVALID;
  for (size_t s = 0; s < desc->nspecies; ++s) {
    if (!desc->species[s].name) goto fail;
    for (size_t j = 0; j < s; ++j) if (!strcmp(desc->species[j].name, desc->species[s].name)) goto fail;
    m->species_to_variable[s] = SIZE_MAX;
    if (!desc->species[s].parameterized) {
      m->species_to_variable[s] = m->n; m->variable_to_species[m->n++] = s;
    }
  }
  if (!m->n || m->n > SIZE_MAX / m->n) goto fail;
  if (reorder) {
    unsigned char *pattern = ck_allocate(m->n * m->n, 1);
    size_t *permutation = ck_allocate(m->n, sizeof(size_t));
    size_t *old_variables = ck_allocate(m->n, sizeof(size_t));
    if (!pattern || !permutation || !old_variables) {
      free(pattern); free(permutation); free(old_variables); error = CK_MEMORY; goto fail;
    }
    for (size_t r = 0; r < m->nr; ++r) {
      const ck_reaction *rx = &desc->reactions[r];
      for (size_t k = 0; k < rx->nreactants; ++k) {
        if (rx->reactants[k] >= desc->nspecies) { free(pattern); free(permutation); free(old_variables); goto fail; }
        size_t col = m->species_to_variable[rx->reactants[k]];
        if (col == SIZE_MAX) continue;
        for (size_t j = 0; j < rx->nreactants + rx->nproducts; ++j) {
          size_t sp = j < rx->nreactants ? rx->reactants[j] : rx->products[j - rx->nreactants].species;
          if (sp >= desc->nspecies) { free(pattern); free(permutation); free(old_variables); goto fail; }
          size_t row = m->species_to_variable[sp];
          if (row != SIZE_MAX) pattern[row * m->n + col] = 1;
        }
      }
    }
    error = ck_markowitz_order(m->n, pattern, permutation);
    memcpy(old_variables, m->variable_to_species, m->n * sizeof(size_t));
    if (!error) for (size_t i = 0; i < m->n; ++i) {
      size_t sp = old_variables[permutation[i]];
      m->variable_to_species[i] = sp; m->species_to_variable[sp] = i;
    }
    free(pattern); free(permutation); free(old_variables);
    if (error) goto fail;
    error = CK_INVALID;
  }
  for (size_t r = 0; r < m->nr; ++r) {
    m->reaction_order[r] = r;
    size_t j = r;
    while (j > 0 && rate_priority(desc->reactions[m->reaction_order[j - 1]].rate.type) > rate_priority(desc->reactions[r].rate.type)) {
      m->reaction_order[j] = m->reaction_order[j - 1]; --j;
    }
    m->reaction_order[j] = r;
  }
  for (size_t r = 0; r < m->nr; ++r) {
    const ck_reaction *rx = &desc->reactions[m->reaction_order[r]];
    if ((rx->nreactants && !rx->reactants) || (rx->nproducts && !rx->products)) goto fail;
    m->reactants[r] = ck_allocate(rx->nreactants, sizeof(size_t));
    m->products[r] = ck_allocate(rx->nproducts, sizeof(ck_product));
    if (!m->reactants[r] || !m->products[r]) { error = CK_MEMORY; goto fail; }
    for (size_t j = 0; j < rx->nreactants; ++j) {
      size_t s = rx->reactants[j];
      if (s >= desc->nspecies) goto fail;
      if (m->species_to_variable[s] != SIZE_MAX) m->reactants[r][m->nreactants[r]++] = m->species_to_variable[s];
    }
    for (size_t j = 0; j < rx->nproducts; ++j) {
      ck_product p = rx->products[j];
      if (p.species >= desc->nspecies) goto fail;
      if (m->species_to_variable[p.species] == SIZE_MAX) continue;
      p.species = m->species_to_variable[p.species]; m->products[r][m->nproducts[r]++] = p;
    }
    m->rates[r] = rx->rate;
    error = ck_rate_prepare(&m->rates[r]);
    if (error) goto fail;
    m->custom_offset[r] = m->nc;
    if (rx->rate.type == CK_USER) m->custom_names[m->nc++] = custom_label(rx->rate.label, "");
    else if (rx->rate.type == CK_SURFACE) {
      m->custom_names[m->nc++] = custom_label(rx->rate.label, ".effective radius [m]");
      m->custom_names[m->nc++] = custom_label(rx->rate.label, ".particle number concentration [# m-3]");
    }
    for (size_t c = m->custom_offset[r]; c < m->nc; ++c) if (!m->custom_names[c]) { error = CK_INVALID; goto fail; }
    error = CK_INVALID;
  }
  m->algebraic = ck_allocate(m->n, 1);
  m->constraints = ck_allocate(desc->nconstraints, sizeof *m->constraints);
  if (!m->algebraic || !m->constraints) { error = CK_MEMORY; goto fail; }
  for (size_t k = 0; k < desc->nconstraints; ++k) {
    const ck_constraint *constraint = &desc->constraints[k];
    error = ck_constraint_validate(constraint); if (error) goto fail;
    error = CK_UNKNOWN;
    if (constraint->algebraic_species >= desc->nspecies) goto fail;
    size_t row = m->species_to_variable[constraint->algebraic_species];
    if (row == SIZE_MAX) goto fail;
    error = CK_INVALID; if (m->algebraic[row]) goto fail;
    ck_bound_constraint *b = &m->constraints[k];
    b->description = *constraint; b->row = row; m->algebraic[row] = 1;
    b->reactants = ck_allocate(constraint->nreactants, sizeof(ck_product));
    b->products = ck_allocate(constraint->nproducts, sizeof(ck_product));
    ++m->nconstraints;
    if (!b->reactants || !b->products) { error = CK_MEMORY; goto fail; }
    for (size_t j = 0; j < constraint->nreactants + constraint->nproducts; ++j) {
      ck_product term = j < constraint->nreactants ? constraint->reactants[j] : constraint->products[j - constraint->nreactants];
      error = CK_UNKNOWN;
      if (term.species >= desc->nspecies || m->species_to_variable[term.species] == SIZE_MAX) goto fail;
      term.species = m->species_to_variable[term.species];
      if (j < constraint->nreactants) b->reactants[j] = term; else b->products[j - constraint->nreactants] = term;
    }
    b->parameter = SIZE_MAX;
    if (constraint->type == CK_EQUILIBRIUM_CONSTRAINT) {
      char *name = ck_string_copy(constraint->name);
      if (!name) { error = CK_MEMORY; goto fail; }
      size_t suffix = 0; int duplicate;
      do {
        duplicate = 0;
        for (size_t j = 0; j < m->nc; ++j) if (!strcmp(m->custom_names[j], name)) { duplicate = 1; break; }
        if (duplicate) {
          free(name); name = ck_allocate(strlen(constraint->name) + 32, 1);
          if (!name) { error = CK_MEMORY; goto fail; }
          snprintf(name, strlen(constraint->name) + 32, "%s_%zu", constraint->name, ++suffix);
        }
      } while (duplicate);
      b->parameter = m->nc; m->custom_names[m->nc++] = name;
    }
  }
  error = symbolic(m, NULL); if (error) goto fail;
  error = derivatives(m); if (error) goto fail;
  *out = m;
  return CK_OK;
fail:
  ck_model_free(m); return error;
}
size_t ck_model_nspecies(const ck_model *m) { return m ? m->n : 0; }
size_t ck_model_nreactions(const ck_model *m) { return m ? m->nr : 0; }
size_t ck_model_ncustom(const ck_model *m) { return m ? m->nc : 0; }
const char *ck_species_name(const ck_model *m, size_t i) { return m && m->description && i < m->n ? m->description->species[m->variable_to_species[i]].name : NULL; }
const char *ck_custom_name(const ck_model *m, size_t i) { return m && i < m->nc ? m->custom_names[i] : NULL; }
ck_status ck_species_index(const ck_model *m, const char *name, size_t *index)
{
  if (!m || !name || !index) return CK_INVALID;
  for (size_t i = 0; i < m->n; ++i) if (!strcmp(ck_species_name(m, i), name)) { *index = i; return CK_OK; }
  return CK_UNKNOWN;
}
ck_status ck_custom_index(const ck_model *m, const char *name, size_t *index)
{
  if (!m || !name || !index) return CK_INVALID;
  for (size_t i = 0; i < m->nc; ++i) if (!strcmp(m->custom_names[i], name)) { *index = i; return CK_OK; }
  return CK_UNKNOWN;
}
size_t ck_jacobian_size(const ck_model *m) { return m ? m->nnz : 0; }
size_t ck_jacobian_index(const ck_model *m, size_t r, size_t c) { return m && r < m->n && c < m->n ? m->map[r * m->n + c] : SIZE_MAX; }

ck_status ck_markowitz_order(size_t n, const unsigned char *pattern, size_t *permutation)
{
  if (!n || n > SIZE_MAX / n || !pattern || !permutation) return CK_INVALID;
  unsigned char *a = ck_allocate(n * n, 1), *alive = ck_allocate(n, 1);
  if (!a || !alive) { free(a); free(alive); return CK_MEMORY; }
  memcpy(a, pattern, n * n); memset(alive, 1, n);
  for (size_t step = 0; step < n; ++step) {
    size_t pivot = n, best = SIZE_MAX;
    for (size_t v = 0; v < n; ++v) if (alive[v]) {
      size_t rows = 0, cols = 0;
      for (size_t j = 0; j < n; ++j) if (alive[j]) { rows += a[v * n + j] != 0; cols += a[j * n + v] != 0; }
      /* Unsigned arithmetic also matches the reference's isolated-node rule. */
      size_t cost = (rows - 1) * (cols - 1);
      if (pivot == n || cost < best) { best = cost; pivot = v; }
    }
    permutation[step] = pivot; alive[pivot] = 0;
    for (size_t i = 0; i < n; ++i) if (alive[i] && a[i * n + pivot])
      for (size_t j = 0; j < n; ++j) if (alive[j] && a[pivot * n + j]) a[i * n + j] = 1;
  }
  free(a); free(alive); return CK_OK;
}

ck_status ck_linear_model_create(size_t n, const unsigned char *pattern, ck_model **out)
{
  if (!out) return CK_INVALID;
  *out = NULL;
  if (!n || n > SIZE_MAX / n || !pattern) return CK_INVALID;
  ck_model *m = ck_allocate(1, sizeof *m);
  if (!m) return CK_MEMORY;
  m->n = n;
  ck_status error = symbolic(m, pattern);
  if (error) { ck_model_free(m); return error; }
  *out = m; return CK_OK;
}

size_t ck_model_nconstraints(const ck_model *m) { return m ? m->nconstraints : 0; }
int ck_is_algebraic(const ck_model *m, size_t variable) { return m && m->algebraic && variable < m->n ? m->algebraic[variable] : 0; }
