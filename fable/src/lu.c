/* Sparse LU: diagonal Markowitz ordering, symbolic fill, MICM's Mozart and
 * Doolittle in-place decompositions as replayed index lists, and the
 * triangular solves. All routines act on every block of a vector-ordered
 * sparse matrix, with the lane (cell) loop innermost. */
#include "internal.h"

#include <math.h>

/* ---- Greedy minimum-fill ordering ---------------------------------------- */

void tmc_minimum_fill_reorder(int n, const unsigned char *pattern, int *permutation)
{
  unsigned char *live = tmc_malloc((size_t)n * (size_t)n);
  memcpy(live, pattern, (size_t)n * (size_t)n);
  unsigned char *alive = tmc_malloc((size_t)n);
  int *cols = tmc_malloc((size_t)n * sizeof(int));
  int *ins = tmc_malloc((size_t)n * sizeof(int));
  memset(alive, 1, (size_t)n);
  for (int step = 0; step < n; ++step)
  {
    int pivot = -1;
    unsigned long long best_fill = 0, best_cost = 0;
    for (int v = 0; v < n; ++v)
    {
      if (!alive[v])
        continue;
      int nc = 0, ni = 0;
      for (int c = 0; c < n; ++c)
        if (c != v && alive[c] && live[(size_t)v * n + c])
          cols[nc++] = c;
      for (int i = 0; i < n; ++i)
        if (i != v && alive[i] && live[(size_t)i * n + v])
          ins[ni++] = i;
      unsigned long long fill = 0;
      for (int a = 0; a < ni; ++a)
        for (int b = 0; b < nc; ++b)
          if (!live[(size_t)ins[a] * n + cols[b]])
            ++fill;
      unsigned long long cost = (unsigned long long)nc * (unsigned long long)ni;
      if (pivot < 0 || fill < best_fill || (fill == best_fill && cost < best_cost))
      {
        best_fill = fill;
        best_cost = cost;
        pivot = v;
      }
    }
    permutation[step] = pivot;
    alive[pivot] = 0;
    for (int i = 0; i < n; ++i)
    {
      if (i == pivot || !alive[i] || !live[(size_t)i * n + pivot])
        continue;
      for (int c = 0; c < n; ++c)
        if (c != pivot && alive[c] && live[(size_t)pivot * n + c])
          live[(size_t)i * n + c] = 1;
    }
  }
  free(ins);
  free(cols);
  free(alive);
  free(live);
}

/* ---- Diagonal Markowitz ordering (MICM's DiagonalMarkowitzReorder) ------ */

void tmc_markowitz_reorder(int n, const unsigned char *pattern, int *permutation)
{
  /* Live pattern as a dense byte matrix, degrees maintained exactly. */
  unsigned char *live = tmc_malloc((size_t)n * (size_t)n);
  memcpy(live, pattern, (size_t)n * (size_t)n);
  int *row_deg = tmc_calloc((size_t)n, sizeof(int));
  int *col_deg = tmc_calloc((size_t)n, sizeof(int));
  unsigned char *alive = tmc_malloc((size_t)n);
  for (int i = 0; i < n; ++i)
  {
    alive[i] = 1;
    for (int j = 0; j < n; ++j)
      if (live[(size_t)i * n + j])
      {
        ++row_deg[i];
        ++col_deg[j];
      }
  }
  for (int step = 0; step < n; ++step)
  {
    int pivot = -1;
    unsigned long long best = 0;
    for (int v = 0; v < n; ++v)
    {
      if (!alive[v])
        continue;
      /* MICM computes the cost in std::size_t, so a degree of zero wraps around */
      unsigned long long cost = ((unsigned long long)row_deg[v] - 1ULL) * ((unsigned long long)col_deg[v] - 1ULL);
      if (pivot < 0 || cost < best)
      {
        best = cost;
        pivot = v;
      }
    }
    permutation[step] = pivot;
    alive[pivot] = 0;
    /* fill: every live in-neighbor i (i -> pivot) gets an edge to every live out-neighbor c */
    for (int i = 0; i < n; ++i)
    {
      if (i == pivot || !alive[i] || !live[(size_t)i * n + pivot])
        continue;
      for (int c = 0; c < n; ++c)
      {
        if (c == pivot || !alive[c] || !live[(size_t)pivot * n + c])
          continue;
        if (!live[(size_t)i * n + c])
        {
          live[(size_t)i * n + c] = 1;
          ++row_deg[i];
          ++col_deg[c];
        }
      }
    }
    for (int c = 0; c < n; ++c)
      if (c != pivot && alive[c] && live[(size_t)pivot * n + c])
      {
        live[(size_t)pivot * n + c] = 0;
        --col_deg[c];
      }
    for (int i = 0; i < n; ++i)
      if (i != pivot && alive[i] && live[(size_t)i * n + pivot])
      {
        live[(size_t)i * n + pivot] = 0;
        --row_deg[i];
      }
  }
  free(live);
  free(row_deg);
  free(col_deg);
  free(alive);
}

/* ---- Symbolic factorization ------------------------------------------- */

tmc_lu *tmc_lu_create(const tmc_sparse *A, tmc_lu_algorithm algorithm, tmc_error *err)
{
  int n = A->n;
  unsigned char *p = tmc_calloc((size_t)n * (size_t)n, 1);
  for (int e = 0; e < A->nnz; ++e)
    p[(size_t)A->row_of[e] * n + A->col_of[e]] = 1;
  for (int i = 0; i < n; ++i)
  {
    if (!p[(size_t)i * n + i])
    {
      free(p);
      tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_ZERO_ELEMENT_ACCESS, "Diagonal element (%d, %d) is zero in LU decomposition", i, i);
      return NULL;
    }
    for (int k = i + 1; k < n; ++k)
    {
      if (!p[(size_t)i * n + k])
        continue;
      for (int j = i + 1; j < n; ++j)
        if (p[(size_t)j * n + i])
          p[(size_t)j * n + k] = 1;
    }
  }
  tmc_sparse_builder *b = tmc_sparse_builder_create(n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (p[(size_t)i * n + j])
        tmc_sparse_builder_with_element(b, i, j, NULL);
  tmc_sparse *LU = tmc_sparse_build(b, A->ordering, A->L, true);
  tmc_sparse_builder_free(b);
  free(p);

  tmc_lu *lu = tmc_calloc(1, sizeof *lu);
  lu->algorithm = algorithm;
  lu->n = n;
  lu->pattern = LU;
  int nnz = LU->nnz;
#define ID(i, j) (LU->lookup[(size_t)(i) * n + (j)])

  /* Mozart plan */
  lu->mz_diag = tmc_malloc((size_t)n * sizeof(int));
  lu->mz_nji = tmc_calloc((size_t)n, sizeof(int));
  lu->mz_nki = tmc_calloc((size_t)n, sizeof(int));
  tmc_ivec aji = { 0 }, aik = { 0 }, njk = { 0 }, ajk = { 0 }, ajk_aji = { 0 };
  for (int i = 0; i < n; ++i)
  {
    lu->mz_diag[i] = ID(i, i);
    for (int j = i + 1; j < n; ++j)
      if (ID(j, i) >= 0)
      {
        tmc_ivec_push(&aji, ID(j, i));
        ++lu->mz_nji[i];
      }
    for (int k = i + 1; k < n; ++k)
    {
      if (ID(i, k) < 0)
        continue;
      int count = 0;
      for (int j = i + 1; j < n; ++j)
        if (ID(j, i) >= 0)
        {
          tmc_ivec_push(&ajk, ID(j, k));
          tmc_ivec_push(&ajk_aji, ID(j, i));
          ++count;
        }
      tmc_ivec_push(&aik, ID(i, k));
      tmc_ivec_push(&njk, count);
      ++lu->mz_nki[i];
    }
  }
  lu->mz_aji = aji.v;
  lu->mz_aik = aik.v;
  lu->mz_njk = njk.v;
  lu->mz_ajk = ajk.v;
  lu->mz_ajk_aji = ajk_aji.v;

  /* Doolittle plan: per row i, the U entries (k >= i) then the L entries (k > i in column i) */
  lu->dl_n_upper = tmc_calloc((size_t)n, sizeof(int));
  lu->dl_n_lower = tmc_calloc((size_t)n, sizeof(int));
  tmc_ivec elem = { 0 }, nterms = { 0 }, terms = { 0 }, diag_of = { 0 };
  for (int i = 0; i < n; ++i)
  {
    for (int k = i; k < n; ++k)
    {
      if (ID(i, k) < 0)
        continue;
      int count = 0;
      for (int j = 0; j < i; ++j)
        if (ID(i, j) >= 0 && ID(j, k) >= 0)
        {
          tmc_ivec_push(&terms, ID(i, j));
          tmc_ivec_push(&terms, ID(j, k));
          ++count;
        }
      tmc_ivec_push(&elem, ID(i, k));
      tmc_ivec_push(&nterms, count);
      tmc_ivec_push(&diag_of, -1);
      ++lu->dl_n_upper[i];
    }
    for (int k = i + 1; k < n; ++k)
    {
      if (ID(k, i) < 0)
        continue;
      int count = 0;
      for (int j = 0; j < i; ++j)
        if (ID(k, j) >= 0 && ID(j, i) >= 0)
        {
          tmc_ivec_push(&terms, ID(k, j));
          tmc_ivec_push(&terms, ID(j, i));
          ++count;
        }
      tmc_ivec_push(&elem, ID(k, i));
      tmc_ivec_push(&nterms, count);
      tmc_ivec_push(&diag_of, ID(i, i));
      ++lu->dl_n_lower[i];
    }
  }
  lu->dl_elem = elem.v;
  lu->dl_nterms = nterms.v;
  lu->dl_terms = terms.v;
  lu->dl_diag_of_col = diag_of.v;
  (void)nnz;

  /* Solve plan */
  lu->fw_count = tmc_calloc((size_t)n, sizeof(int));
  lu->bw_count = tmc_calloc((size_t)n, sizeof(int));
  lu->bw_diag = tmc_malloc((size_t)n * sizeof(int));
  tmc_ivec fe = { 0 }, fc = { 0 }, be = { 0 }, bc = { 0 };
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < i; ++j)
      if (ID(i, j) >= 0)
      {
        tmc_ivec_push(&fe, ID(i, j));
        tmc_ivec_push(&fc, j);
        ++lu->fw_count[i];
      }
  for (int r = 0; r < n; ++r)
  {
    int i = n - 1 - r;
    for (int j = i + 1; j < n; ++j)
      if (ID(i, j) >= 0)
      {
        tmc_ivec_push(&be, ID(i, j));
        tmc_ivec_push(&bc, j);
        ++lu->bw_count[r];
      }
    lu->bw_diag[r] = ID(i, i);
  }
  lu->fw_elem = fe.v;
  lu->fw_col = fc.v;
  lu->bw_elem = be.v;
  lu->bw_col = bc.v;
  lu->n_lower = fe.n;
  lu->n_upper = be.n;
#undef ID
  return lu;
}

void tmc_lu_free(tmc_lu *lu)
{
  if (!lu)
    return;
  tmc_sparse_free(lu->pattern);
  free(lu->mz_diag);
  free(lu->mz_nji);
  free(lu->mz_nki);
  free(lu->mz_aji);
  free(lu->mz_aik);
  free(lu->mz_njk);
  free(lu->mz_ajk);
  free(lu->mz_ajk_aji);
  free(lu->dl_n_lower);
  free(lu->dl_n_upper);
  free(lu->dl_elem);
  free(lu->dl_nterms);
  free(lu->dl_terms);
  free(lu->dl_diag_of_col);
  free(lu->fw_count);
  free(lu->fw_elem);
  free(lu->fw_col);
  free(lu->bw_count);
  free(lu->bw_elem);
  free(lu->bw_col);
  free(lu->bw_diag);
  free(lu);
}

tmc_lu_algorithm tmc_lu_algorithm_of(const tmc_lu *lu)
{
  return lu->algorithm;
}

tmc_sparse *tmc_lu_matrix(const tmc_lu *lu, int n_blocks, double initial_value)
{
  return tmc_sparse_create_like(lu->pattern, n_blocks, initial_value);
}

static tmc_sparse *triangle(const tmc_lu *lu, int n_blocks, double initial_value, bool lower)
{
  const tmc_sparse *P = lu->pattern;
  tmc_sparse_builder *b = tmc_sparse_builder_create(P->n);
  for (int e = 0; e < P->nnz; ++e)
  {
    int i = P->row_of[e], j = P->col_of[e];
    if ((lower && j <= i) || (!lower && j >= i))
      tmc_sparse_builder_with_element(b, i, j, NULL);
  }
  tmc_sparse_builder_set_number_of_blocks(b, n_blocks);
  tmc_sparse_builder_set_initial_value(b, initial_value);
  tmc_sparse *m = tmc_sparse_build(b, P->ordering, P->L, false);
  tmc_sparse_builder_free(b);
  return m;
}

tmc_sparse *tmc_lu_lower_matrix(const tmc_lu *lu, int n_blocks, double initial_value)
{
  return triangle(lu, n_blocks, initial_value, true);
}

tmc_sparse *tmc_lu_upper_matrix(const tmc_lu *lu, int n_blocks, double initial_value)
{
  return triangle(lu, n_blocks, initial_value, false);
}

/* ---- Decompositions ---------------------------------------------------- */

void tmc_lu_decompose_mozart_in_place(const tmc_lu *lu, tmc_sparse *ALU, double *reciprocals)
{
  const int L = ALU->L, n = lu->n;
  const size_t group_stride = (size_t)ALU->nnz * (size_t)L;
  double *inv = tmc_malloc((size_t)L * sizeof(double));
  for (int g = 0; g < ALU->n_groups; ++g)
  {
    double *restrict a = ALU->data + (size_t)g * group_stride;
    double *recip = reciprocals ? reciprocals + (size_t)g * (size_t)n * (size_t)L : NULL;
    const int *aji = lu->mz_aji;
    const int *aik = lu->mz_aik;
    const int *njk = lu->mz_njk;
    const int *ajk = lu->mz_ajk;
    const int *ajk_aji = lu->mz_ajk_aji;
    for (int i = 0; i < n; ++i)
    {
      const double *d = a + (size_t)lu->mz_diag[i] * L;
      for (int c = 0; c < L; ++c)
        inv[c] = 1.0 / d[c];
      if (recip)
        for (int c = 0; c < L; ++c)
          recip[(size_t)i * L + c] = inv[c];
      for (int q = 0; q < lu->mz_nji[i]; ++q, ++aji)
      {
        double *t = a + (size_t)*aji * L;
#pragma omp simd
        for (int c = 0; c < L; ++c)
          t[c] *= inv[c];
      }
      for (int q = 0; q < lu->mz_nki[i]; ++q, ++aik, ++njk)
      {
        const double *r = a + (size_t)*aik * L;
        for (int m = 0; m < *njk; ++m, ++ajk, ++ajk_aji)
        {
          double *t = a + (size_t)*ajk * L;
          const double *l = a + (size_t)*ajk_aji * L;
#pragma omp simd
          for (int c = 0; c < L; ++c)
            t[c] -= l[c] * r[c];
        }
      }
    }
  }
  free(inv);
}

void tmc_lu_decompose_doolittle_in_place(const tmc_lu *lu, tmc_sparse *ALU, double *reciprocals)
{
  const int L = ALU->L, n = lu->n;
  const size_t group_stride = (size_t)ALU->nnz * (size_t)L;
  for (int g = 0; g < ALU->n_groups; ++g)
  {
    double *restrict a = ALU->data + (size_t)g * group_stride;
    double *recip = reciprocals ? reciprocals + (size_t)g * (size_t)n * (size_t)L : NULL;
    const int *elem = lu->dl_elem;
    const int *nterms = lu->dl_nterms;
    const int *terms = lu->dl_terms;
    const int *diag_of = lu->dl_diag_of_col;
    for (int i = 0; i < n; ++i)
    {
      int count = lu->dl_n_upper[i] + lu->dl_n_lower[i];
      for (int q = 0; q < count; ++q, ++elem, ++nterms, ++diag_of)
      {
        double *t = a + (size_t)*elem * L;
        for (int m = 0; m < *nterms; ++m, terms += 2)
        {
          const double *x = a + (size_t)terms[0] * L;
          const double *y = a + (size_t)terms[1] * L;
#pragma omp simd
          for (int c = 0; c < L; ++c)
            t[c] -= x[c] * y[c];
        }
        if (*diag_of >= 0)
        {
          const double *d = a + (size_t)*diag_of * L;
#pragma omp simd
          for (int c = 0; c < L; ++c)
            t[c] /= d[c];
        }
      }
      if (recip)
      {
        const double *d = a + (size_t)lu->mz_diag[i] * L;
        for (int c = 0; c < L; ++c)
          recip[(size_t)i * L + c] = 1.0 / d[c];
      }
    }
  }
}

void tmc_lu_decompose_in_place(const tmc_lu *lu, tmc_sparse *ALU, double *reciprocals)
{
  if (lu->algorithm == TMC_LU_DOOLITTLE)
    tmc_lu_decompose_doolittle_in_place(lu, ALU, reciprocals);
  else
    tmc_lu_decompose_mozart_in_place(lu, ALU, reciprocals);
}

void tmc_lu_decompose(const tmc_lu *lu, const tmc_sparse *A, tmc_sparse *Lm, tmc_sparse *U)
{
  /* Factor a copy of A in the combined pattern, then split it. */
  tmc_sparse *ALU = tmc_lu_matrix(lu, A->n_blocks, 0.0);
  for (int b = 0; b < A->n_blocks; ++b)
    for (int e = 0; e < A->nnz; ++e)
      tmc_sparse_set(ALU, b, A->row_of[e], A->col_of[e], A->data[tmc_sparse_data_index(A, b, e)]);
  tmc_lu_decompose_in_place(lu, ALU, NULL);
  for (int b = 0; b < A->n_blocks; ++b)
    for (int e = 0; e < ALU->nnz; ++e)
    {
      int i = ALU->row_of[e], j = ALU->col_of[e];
      double v = ALU->data[tmc_sparse_data_index(ALU, b, e)];
      if (j < i)
        tmc_sparse_set(Lm, b, i, j, v);
      else
        tmc_sparse_set(U, b, i, j, v);
      if (i == j)
        tmc_sparse_set(Lm, b, i, i, 1.0);
    }
  tmc_sparse_free(ALU);
}

/* ---- Solves ----------------------------------------------------------- */

void tmc_lu_solve_in_place(const tmc_lu *lu, const tmc_sparse *ALU, const double *reciprocals, tmc_dense *x)
{
  const int L = ALU->L, n = lu->n;
  const size_t group_stride = (size_t)ALU->nnz * (size_t)L;
  const size_t x_stride = (size_t)n * (size_t)L;
  for (int g = 0; g < ALU->n_groups; ++g)
  {
    const double *restrict a = ALU->data + (size_t)g * group_stride;
    double *restrict xg = x->data + (size_t)g * x_stride;
    const double *recip = reciprocals ? reciprocals + (size_t)g * x_stride : NULL;
    const int *fe = lu->fw_elem, *fc = lu->fw_col;
    for (int i = 0; i < n; ++i)
    {
      double *xi = xg + (size_t)i * L;
      for (int q = 0; q < lu->fw_count[i]; ++q, ++fe, ++fc)
      {
        const double *l = a + (size_t)*fe * L;
        const double *xj = xg + (size_t)*fc * L;
#pragma omp simd
        for (int c = 0; c < L; ++c)
          xi[c] -= l[c] * xj[c];
      }
    }
    const int *be = lu->bw_elem, *bc = lu->bw_col;
    for (int r = 0; r < n; ++r)
    {
      int i = n - 1 - r;
      double *xi = xg + (size_t)i * L;
      for (int q = 0; q < lu->bw_count[r]; ++q, ++be, ++bc)
      {
        const double *u = a + (size_t)*be * L;
        const double *xj = xg + (size_t)*bc * L;
#pragma omp simd
        for (int c = 0; c < L; ++c)
          xi[c] -= u[c] * xj[c];
      }
      if (recip)
      {
        const double *d = recip + (size_t)i * L;
#pragma omp simd
        for (int c = 0; c < L; ++c)
          xi[c] *= d[c];
      }
      else
      {
        const double *d = a + (size_t)lu->bw_diag[r] * L;
#pragma omp simd
        for (int c = 0; c < L; ++c)
          xi[c] /= d[c];
      }
    }
  }
}

void tmc_lu_solve(const tmc_lu *lu, const tmc_sparse *Lm, const tmc_sparse *U, tmc_dense *x)
{
  const int n = lu->n;
  for (int b = 0; b < Lm->n_blocks; ++b)
  {
    for (int i = 0; i < n; ++i)
    {
      double v = tmc_dense_get(x, b, i);
      for (int j = 0; j < i; ++j)
        if (!tmc_sparse_is_zero(Lm, i, j))
          v -= tmc_sparse_get(Lm, b, i, j) * tmc_dense_get(x, b, j);
      tmc_dense_set(x, b, i, v);
    }
    for (int i = n - 1; i >= 0; --i)
    {
      double v = tmc_dense_get(x, b, i);
      for (int j = i + 1; j < n; ++j)
        if (!tmc_sparse_is_zero(U, i, j))
          v -= tmc_sparse_get(U, b, i, j) * tmc_dense_get(x, b, j);
      tmc_dense_set(x, b, i, v / tmc_sparse_get(U, b, i, i));
    }
  }
}

bool tmc_lu_is_singular(const tmc_lu *lu, const tmc_sparse *ALU, int block)
{
  for (int i = 0; i < lu->n; ++i)
  {
    double d = ALU->data[tmc_sparse_data_index(ALU, block, lu->mz_diag[i])];
    if (d == 0.0 || !isfinite(d))
      return true;
  }
  return false;
}
