/* Sparse Mozart and Doolittle algorithms ported from MICM.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <limits.h>

int ch_markowitz(size_t n, const unsigned char *pattern, int *permutation) {
  if ((!pattern || !permutation) && n) return CH_INVALID;
  if (n > INT_MAX || (n && n > SIZE_MAX/n)) return CH_SIZE;
  unsigned char *a = ch_allocate(n*n, 1), *alive = ch_allocate(n, 1);
  if (!a || !alive) { free(a); free(alive); return CH_MEMORY; }
  memcpy(a, pattern, n*n);
  memset(alive, 1, n);
  for (size_t s = 0; s < n; ++s) {
    size_t pivot = n, best = SIZE_MAX;
    for (size_t v = 0; v < n; ++v) if (alive[v]) {
      size_t rows = 0, cols = 0;
      for (size_t j = 0; j < n; ++j) if (alive[j]) {
        rows += a[v*n+j] != 0;
        cols += a[j*n+v] != 0;
      }
      size_t cost = (rows - 1) * (cols - 1);
      if (pivot == n || cost < best) { pivot = v; best = cost; }
    }
    permutation[s] = (int)pivot;
    alive[pivot] = 0;
    for (size_t i = 0; i < n; ++i) if (alive[i] && a[i*n+pivot])
      for (size_t j = 0; j < n; ++j) if (alive[j] && a[pivot*n+j]) a[i*n+j] = 1;
  }
  free(a); free(alive);
  return CH_OK;
}

ch_sparse *ch_sparse_create(size_t n, const unsigned char *pattern) {
  if (n > INT_MAX || (n && !pattern) || (n && n > SIZE_MAX/n)) return NULL;
  ch_sparse *s = calloc(1, sizeof(*s));
  unsigned char *p = ch_allocate(n*n, 1);
  if (!s || !p) { free(s); free(p); return NULL; }
  s->n = n;
  memcpy(p, pattern, n*n);
  for (size_t i = 0; i < n; ++i) p[i*n+i] = 1;
  /* Symbolic elimination is independent of the numerical initial values. */
  for (size_t i = 0; i < n; ++i)
    for (size_t k = i+1; k < n; ++k) if (p[i*n+k])
      for (size_t j = i+1; j < n; ++j) if (p[j*n+i]) p[j*n+k] = 1;
  s->index = ch_allocate(n*n, sizeof(int));
  s->row_start = ch_allocate(n+1, sizeof(int));
  s->diagonal = ch_allocate(n, sizeof(int));
  if (!s->index || !s->row_start || !s->diagonal) goto fail;
  for (size_t i = 0; i < n; ++i) {
    s->row_start[i] = (int)s->nnz;
    for (size_t j = 0; j < n; ++j)
      if (p[i*n+j]) {
        if (s->nnz == INT_MAX) goto fail;
        s->index[i*n+j] = (int)s->nnz++;
      } else s->index[i*n+j] = -1;
    s->diagonal[i] = s->index[i*n+i];
  }
  s->row_start[n] = (int)s->nnz;
  s->columns = ch_allocate(s->nnz, sizeof(int));
  if (!s->columns) goto fail;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      if (p[i*n+j]) s->columns[s->index[i*n+j]] = (int)j;
  /* One inverse per pivot, lower-column scaling, then outer-product updates. */
  size_t count = n;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i+1; j < n; ++j) if (p[j*n+i]) ++count;
    for (size_t k = i+1; k < n; ++k) if (p[i*n+k])
      for (size_t j = i+1; j < n; ++j) if (p[j*n+i]) ++count;
  }
  s->ops = ch_allocate(count, sizeof(ch_lu_op));
  if (!s->ops) goto fail;
  for (size_t i = 0; i < n; ++i) {
    s->ops[s->nops++] = (ch_lu_op){s->diagonal[i], (int)i, 0, 0};
    for (size_t j = i+1; j < n; ++j) if (p[j*n+i])
      s->ops[s->nops++] = (ch_lu_op){s->index[j*n+i], 0, 0, 1};
    for (size_t k = i+1; k < n; ++k) if (p[i*n+k])
      for (size_t j = i+1; j < n; ++j) if (p[j*n+i])
        s->ops[s->nops++] = (ch_lu_op){s->index[j*n+k],s->index[j*n+i],s->index[i*n+k],2};
  }
  free(p);
  return s;
fail:
  free(p); ch_sparse_free(s); return NULL;
}

void ch_sparse_free(ch_sparse *s) {
  if (!s) return;
  free(s->index); free(s->row_start); free(s->columns); free(s->diagonal); free(s->ops); free(s);
}
size_t ch_sparse_size(const ch_sparse *s) { return s ? s->nnz : 0; }
size_t ch_sparse_dimension(const ch_sparse *s) { return s ? s->n : 0; }
const int *ch_sparse_row_offsets(const ch_sparse *s) { return s ? s->row_start : NULL; }
const int *ch_sparse_column_indices(const ch_sparse *s) { return s ? s->columns : NULL; }
const int *ch_sparse_diagonal_indices(const ch_sparse *s) { return s ? s->diagonal : NULL; }
int ch_sparse_index(const ch_sparse *s, size_t row, size_t col) {
  return s && row < s->n && col < s->n ? s->index[row*s->n+col] : CH_INVALID;
}

void ch_sparse_multiply(const ch_sparse *s, size_t cells, const double *a, const double *x, double *output) {
  if (!s || !a || !x || !output) return;
  memset(output,0,s->n*cells*sizeof(double));
  for (size_t i = 0; i < s->n; ++i) {
    double *restrict dst = output+i*cells;
    for (int j = s->row_start[i]; j < s->row_start[i+1]; ++j) {
      const double *restrict coefficient = a+(size_t)j*cells;
      const double *restrict value = x+(size_t)s->columns[j]*cells;
      #pragma omp simd
      for (size_t c = 0; c < cells; ++c) dst[c] += coefficient[c]*value[c];
    }
  }
}

int ch_lu_factor(const ch_sparse *s, size_t cells, double *a, int algorithm) {
  if (!s || !a || !cells || (algorithm != 0 && algorithm != 1)) return CH_INVALID;
  double *inverse = ch_allocate(cells, sizeof(double));
  if (!inverse) return CH_MEMORY;
  int result = ch_lu_factor_work(s,cells,a,algorithm,inverse,NULL);
  free(inverse);
  return result;
}

int ch_lu_factor_copy(const ch_sparse *s, size_t cells, const double *a, double *lu, int algorithm) {
  if (!s || !a || !lu) return CH_INVALID;
  if (a != lu) memcpy(lu,a,s->nnz*cells*sizeof(double));
  return ch_lu_factor(s,cells,lu,algorithm);
}

int ch_lu_factor_work(const ch_sparse *s, size_t cells, double *a, int algorithm,
                      double *inverse, double *inverse_diagonal) {
  if (algorithm == 0) {
    for (size_t i = 0; i < s->nops; ++i) {
      const ch_lu_op *op = &s->ops[i];
      double *restrict dst = a + (size_t)op->destination*cells;
      if (op->kind == 0) {
        if (inverse_diagonal) inverse = inverse_diagonal+(size_t)op->left*cells;
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) inverse[c] = 1.0/dst[c];
      } else if (op->kind == 1) {
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) dst[c] *= inverse[c];
      } else {
        const double *restrict left = a + (size_t)op->left*cells;
        const double *restrict right = a + (size_t)op->right*cells;
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) dst[c] -= left[c]*right[c];
      }
    }
  } else {
    size_t n = s->n;
    for (size_t i = 0; i < n; ++i) {
      if (inverse_diagonal) inverse = inverse_diagonal+i*cells;
      for (size_t k = i; k < n; ++k) {
        int dst = s->index[i*n+k];
        if (dst < 0) continue;
        for (size_t j = 0; j < i; ++j) {
          int l = s->index[i*n+j], r = s->index[j*n+k];
          if (l < 0 || r < 0) continue;
          #pragma omp simd
          for (size_t c = 0; c < cells; ++c)
            a[(size_t)dst*cells+c] -= a[(size_t)l*cells+c]*a[(size_t)r*cells+c];
        }
      }
      #pragma omp simd
      for (size_t c = 0; c < cells; ++c) inverse[c] = 1.0/a[(size_t)s->diagonal[i]*cells+c];
      for (size_t k = i+1; k < n; ++k) {
        int dst = s->index[k*n+i];
        if (dst < 0) continue;
        for (size_t j = 0; j < i; ++j) {
          int l = s->index[k*n+j], r = s->index[j*n+i];
          if (l < 0 || r < 0) continue;
          #pragma omp simd
          for (size_t c = 0; c < cells; ++c)
            a[(size_t)dst*cells+c] -= a[(size_t)l*cells+c]*a[(size_t)r*cells+c];
        }
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) a[(size_t)dst*cells+c] *= inverse[c];
      }
    }
  }
  return CH_OK;
}

void ch_lu_solve(const ch_sparse *s, size_t cells, const double *a, double *x) {
  ch_lu_solve_work(s,cells,a,x,NULL);
}

void ch_lu_solve_work(const ch_sparse *s, size_t cells, const double *a, double *x,
                       const double *inverse_diagonal) {
  for (size_t i = 0; i < s->n; ++i) {
    double *restrict dst = x + i*cells;
    for (int j = s->row_start[i]; j < s->diagonal[i]; ++j) {
      const double *restrict l = a + (size_t)j*cells;
      const double *restrict v = x + (size_t)s->columns[j]*cells;
      #pragma omp simd
      for (size_t c = 0; c < cells; ++c) dst[c] -= l[c]*v[c];
    }
  }
  for (size_t i = s->n; i-- > 0;) {
    double *restrict dst = x + i*cells;
    for (int j = s->diagonal[i]+1; j < s->row_start[i+1]; ++j) {
      const double *restrict u = a + (size_t)j*cells;
      const double *restrict v = x + (size_t)s->columns[j]*cells;
      #pragma omp simd
      for (size_t c = 0; c < cells; ++c) dst[c] -= u[c]*v[c];
    }
    if (inverse_diagonal) {
      const double *restrict inverse = inverse_diagonal+i*cells;
      #pragma omp simd
      for (size_t c = 0; c < cells; ++c) dst[c] *= inverse[c];
    } else {
      const double *restrict d = a + (size_t)s->diagonal[i]*cells;
      #pragma omp simd
      for (size_t c = 0; c < cells; ++c) dst[c] /= d[c];
    }
  }
}
