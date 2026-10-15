#include "internal.h"

/* Every operation on a species/reaction vector traverses contiguous cells.
   The reaction and derivative order matches the reference, including duplicate
   reactants and products; cancellations are deliberately not folded. */
void ck_forcing(const ck_state *s, const double *restrict y, double *restrict f)
{
  const ck_model *m = s->model;
  const size_t cells = s->cells;
  memset(f, 0, m->n * cells * sizeof(double));
  for (size_t r = 0; r < m->nr; ++r) {
    const double *restrict k = s->rates + r * cells;
    /* Cell-local rate avoids a temporary allocation or shared workspace. */
    #pragma omp simd
    for (size_t c = 0; c < cells; ++c) {
      double rate = k[c];
      for (size_t j = 0; j < m->nreactants[r]; ++j) rate *= y[m->reactants[r][j] * cells + c];
      for (size_t j = 0; j < m->nreactants[r]; ++j) f[m->reactants[r][j] * cells + c] -= rate;
      for (size_t j = 0; j < m->nproducts[r]; ++j) f[m->products[r][j].species * cells + c] += m->products[r][j].yield * rate;
    }
  }
}
void ck_jacobian(const ck_state *s, const double *restrict y, double *restrict a)
{
  const ck_model *m = s->model;
  const size_t cells = s->cells;
  memset(a, 0, m->nnz * cells * sizeof(double));
  for (size_t job = 0; job < m->njobs; ++job) {
    const ck_derivative *d = &m->jobs[job];
    const double *restrict k = s->rates + d->reaction * cells;
    #pragma omp simd
    for (size_t c = 0; c < cells; ++c) {
      double rate = k[c];
      for (size_t j = 0; j < d->nreactants; ++j) rate *= y[d->reactants[j] * cells + c];
      for (size_t j = 0; j <= d->nreactants; ++j) a[d->targets[j] * cells + c] += rate;
      for (size_t j = 0; j < d->nproducts; ++j) a[d->targets[d->nreactants + 1 + j] * cells + c] -= d->products[j].yield * rate;
    }
  }
}
void ck_factor(const ck_model *m, size_t cells, double *restrict a, ck_lu_algorithm algorithm)
{
  if (algorithm == CK_LU_MOZART) {
    for (size_t i = 0; i < m->n; ++i) {
      const double *restrict diagonal = a + m->diagonal[i] * cells;
      /* Divide all below-diagonal values by the same reciprocal, as Mozart. */
      for (size_t j = m->divide_start[i]; j < m->divide_start[i + 1]; ++j) {
        double *restrict v = a + m->divides[j] * cells;
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) v[c] *= 1.0 / diagonal[c];
      }
      for (size_t j = m->update_start[i]; j < m->update_start[i + 1]; ++j) {
        ck_operation op = m->ops[j];
        double *restrict target = a + op.target * cells;
        const double *restrict left = a + op.left * cells, *restrict right = a + op.right * cells;
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) target[c] -= left[c] * right[c];
      }
    }
  } else {
    /* Doolittle updates each upper row, then each lower column. */
    const size_t n = m->n;
    for (size_t i = 0; i < n; ++i) {
      for (size_t p = m->row_start[i]; p < m->row_start[i + 1]; ++p) {
        size_t k = m->columns[p];
        if (k < i) continue;
        double *restrict v = a + p * cells;
        for (size_t j = 0; j < i; ++j) {
          size_t l = m->map[i * n + j], u = m->map[j * n + k];
          if (l == SIZE_MAX || u == SIZE_MAX) continue;
          const double *restrict left = a + l * cells, *restrict right = a + u * cells;
          #pragma omp simd
          for (size_t c = 0; c < cells; ++c) v[c] -= left[c] * right[c];
        }
      }
      for (size_t k = i + 1; k < n; ++k) {
        size_t off = m->map[k * n + i];
        if (off == SIZE_MAX) continue;
        double *restrict v = a + off * cells;
        for (size_t j = 0; j < i; ++j) {
          size_t l = m->map[k * n + j], u = m->map[j * n + i];
          if (l == SIZE_MAX || u == SIZE_MAX) continue;
          const double *restrict left = a + l * cells, *restrict right = a + u * cells;
          #pragma omp simd
          for (size_t c = 0; c < cells; ++c) v[c] -= left[c] * right[c];
        }
        const double *restrict diag = a + m->diagonal[i] * cells;
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) v[c] /= diag[c];
      }
    }
  }
}
void ck_linear_solve(const ck_model *m, size_t cells, const double *restrict a, double *restrict x)
{
  for (size_t i = 0; i < m->n; ++i) {
    double *restrict xi = x + i * cells;
    for (size_t p = m->row_start[i]; p < m->diagonal[i]; ++p) {
      const double *restrict v = a + p * cells, *restrict xj = x + m->columns[p] * cells;
      #pragma omp simd
      for (size_t c = 0; c < cells; ++c) xi[c] -= v[c] * xj[c];
    }
  }
  for (size_t i = m->n; i-- > 0;) {
    double *restrict xi = x + i * cells;
    for (size_t p = m->diagonal[i] + 1; p < m->row_start[i + 1]; ++p) {
      const double *restrict v = a + p * cells, *restrict xj = x + m->columns[p] * cells;
      #pragma omp simd
      for (size_t c = 0; c < cells; ++c) xi[c] -= v[c] * xj[c];
    }
    const double *restrict diag = a + m->diagonal[i] * cells;
    #pragma omp simd
    for (size_t c = 0; c < cells; ++c) xi[c] /= diag[c];
  }
}
