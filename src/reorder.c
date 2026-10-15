/* Diagonal Markowitz ordering, translated from the reference algorithm.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

ChemStatus chem_markowitz(size_t n, const unsigned char *source, size_t *perm) {
  if (!source || !perm || !n || n > SIZE_MAX/n) return CHEM_INVALID;
  unsigned char *pattern = malloc(n*n), *alive = malloc(n);
  if (!pattern || !alive) { free(pattern); free(alive); return CHEM_MEMORY; }
  memcpy(pattern, source, n*n); memset(alive, 1, n);
  for (size_t i = 0; i < n; ++i) pattern[i*n+i] = 1;
  for (size_t step = 0; step < n; ++step) {
    size_t best = n, cost = SIZE_MAX;
    for (size_t v = 0; v < n; ++v) if (alive[v]) {
      size_t rows = 0, columns = 0;
      for (size_t j = 0; j < n; ++j) if (alive[j]) { rows += !!pattern[v*n+j]; columns += !!pattern[j*n+v]; }
      size_t candidate = (rows-1)*(columns-1);
      if (best == n || candidate < cost) { best = v; cost = candidate; }
    }
    perm[step] = best; alive[best] = 0;
    for (size_t row = 0; row < n; ++row) if (alive[row] && pattern[row*n+best])
      for (size_t col = 0; col < n; ++col) if (alive[col] && pattern[best*n+col]) pattern[row*n+col] = 1;
  }
  free(pattern); free(alive); return CHEM_OK;
}
ChemStatus chem_model_reorder(ChemModel *m) {
  if (!m || m->finalized || !m->nvariables) return CHEM_INVALID;
  size_t n = m->nvariables;
  if (n > SIZE_MAX/n || n > SIZE_MAX/sizeof(size_t)) return CHEM_MEMORY;
  unsigned char *pattern = calloc(n*n, 1);
  size_t *perm = malloc(2*n*sizeof(size_t));
  if (!pattern || !perm) { free(pattern); free(perm); return CHEM_MEMORY; }
  for (size_t i = 0; i < n; ++i) pattern[i*n+i] = 1;
  for (size_t r = 0; r < m->nreactions; ++r) {
    ChemReaction *q = &m->reactions[r];
    for (size_t i = 0; i < q->nreactants; ++i) {
      size_t col = m->species[q->reactants[i]].variable;
      if (col == SIZE_MAX) continue;
      for (size_t j = 0; j < q->nreactants; ++j) {
        size_t row = m->species[q->reactants[j]].variable; if (row != SIZE_MAX) pattern[row*n+col] = 1;
      }
      for (size_t j = 0; j < q->nproducts; ++j) {
        size_t row = m->species[q->products[j]].variable; if (row != SIZE_MAX) pattern[row*n+col] = 1;
      }
    }
  }
  for (size_t i = 0; i < m->nconstraints; ++i) {
    ChemConstraint *q = &m->constraints[i];
    for (size_t j = 0; j < q->nterms; ++j) pattern[q->algebraic*n+q->species[j]] = 1;
  }
  if (m->external_jacobian) memset(pattern, 1, n*n);
  ChemStatus status = chem_markowitz(n, pattern, perm);
  if (!status) {
    size_t *inverse = perm+n;
    for (size_t i = 0; i < n; ++i) inverse[perm[i]] = i;
    for (size_t i = 0; i < m->nspecies; ++i) if (!m->species[i].parameterized) m->species[i].variable = inverse[m->species[i].variable];
    for (size_t i = 0; i < m->nconstraints; ++i) {
      ChemConstraint *q = &m->constraints[i]; q->algebraic = inverse[q->algebraic];
      for (size_t j = 0; j < q->nterms; ++j) q->species[j] = inverse[q->species[j]];
    }
  }
  free(pattern); free(perm); return status;
}
