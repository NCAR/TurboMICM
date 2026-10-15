/* Translated sparse LU algorithms: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

void chem_sparse_multiply(const ChemModel *m, size_t cells, size_t stride, const double *a,
                          const double *x, double *product) {
  memset(product,0,m->nvariables*stride*sizeof(double));
  for (size_t i = 0; i < m->nvariables; ++i) for (size_t j = 0; j < m->nvariables; ++j) {
    int at = m->index[i*m->nvariables+j];
    if (at < 0) continue;
    #pragma omp simd
    for (size_t c = 0; c < cells; ++c) product[i*stride+c] += a[(size_t)at*stride+c]*x[j*stride+c];
  }
}

void chem_lu_mozart(const ChemModel *m, size_t nc, size_t s, double *a) {
  double *inverse = malloc(nc*sizeof(double));
  if (!inverse) { for (size_t c = 0; c < nc; ++c) a[m->diagonal[0]*s+c] = NAN; return; }
  for (size_t i = 0; i < m->nvariables; ++i) {
    /* The reciprocal is shared by the column, exactly as in the reference. */
    #pragma omp simd
    for (size_t c = 0; c < nc; ++c) inverse[c] = 1.0/a[m->diagonal[i]*s+c];
    ChemRange range = m->scale[i];
    for (size_t j = range.first; j < range.first+range.count; ++j) {
      double *restrict target = a+m->scale_ids[j]*s;
      #pragma omp simd
      for (size_t c = 0; c < nc; ++c) target[c] *= inverse[c];
    }
    range = m->updates[i];
    for (size_t j = range.first; j < range.first+range.count; ++j) {
      ChemTriple op = m->update_ops[j];
      double *restrict target = a+op.target*s;
      const double *restrict left = a+op.left*s, *restrict right = a+op.right*s;
      #pragma omp simd
      for (size_t c = 0; c < nc; ++c) target[c] -= left[c]*right[c];
    }
  }
  free(inverse);
}
void chem_lu_doolittle(const ChemModel *m, size_t nc, size_t s, double *a) {
  size_t n = m->nvariables;
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = i; k < n; ++k) {
      int ik = m->index[i*n+k]; if (ik < 0) continue;
      for (size_t j = 0; j < i; ++j) {
        int ij = m->index[i*n+j], jk = m->index[j*n+k]; if (ij < 0 || jk < 0) continue;
        #pragma omp simd
        for (size_t c = 0; c < nc; ++c) a[(size_t)ik*s+c] -= a[(size_t)ij*s+c]*a[(size_t)jk*s+c];
      }
    }
    for (size_t k = i+1; k < n; ++k) {
      int ki = m->index[k*n+i]; if (ki < 0) continue;
      for (size_t j = 0; j < i; ++j) {
        int kj = m->index[k*n+j], ji = m->index[j*n+i]; if (kj < 0 || ji < 0) continue;
        #pragma omp simd
        for (size_t c = 0; c < nc; ++c) a[(size_t)ki*s+c] -= a[(size_t)kj*s+c]*a[(size_t)ji*s+c];
      }
      #pragma omp simd
      for (size_t c = 0; c < nc; ++c) a[(size_t)ki*s+c] /= a[m->diagonal[i]*s+c];
    }
  }
}
void chem_linear_solve(const ChemModel *m, size_t nc, size_t s, const double *a, double *x) {
  for (size_t i = 0; i < m->nvariables; ++i) {
    ChemRange range = m->lower[i];
    double *restrict xi = x+i*s;
    for (size_t j = range.first; j < range.first+range.count; ++j) {
      const double *restrict lij = a+m->lower_ids[j]*s, *restrict xj = x+m->lower_columns[j]*s;
      #pragma omp simd
      for (size_t c = 0; c < nc; ++c) xi[c] -= lij[c]*xj[c];
    }
  }
  for (size_t i = m->nvariables; i-- > 0;) {
    ChemRange range = m->upper[i];
    double *restrict xi = x+i*s;
    for (size_t j = range.first; j < range.first+range.count; ++j) {
      const double *restrict uij = a+m->upper_ids[j]*s, *restrict xj = x+m->upper_columns[j]*s;
      #pragma omp simd
      for (size_t c = 0; c < nc; ++c) xi[c] -= uij[c]*xj[c];
    }
    const double *restrict diagonal = a+m->diagonal[i]*s;
    #pragma omp simd
    for (size_t c = 0; c < nc; ++c) xi[c] /= diagonal[c];
  }
}
ChemStatus chem_lu_separate(const ChemModel *m, size_t nc, size_t s, const double *a, double *lower, double *upper, int doolittle) {
  if (!m || !a || !lower || !upper || nc > s || m->nnz > SIZE_MAX/s/sizeof(double)) return CHEM_INVALID;
  size_t n = m->nvariables;
  double *work = calloc(m->nnz*s,sizeof(double)); if (!work) return CHEM_MEMORY;
  for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) {
    int at = m->index[i*n+j];
    if (at < 0 || (m->declared_sparsity && !m->declared_sparsity[i*n+j])) continue;
    for (size_t c = 0; c < nc; ++c) work[(size_t)at*s+c] = a[(size_t)at*s+c];
  }
  if (doolittle) chem_lu_doolittle(m,nc,s,work); else chem_lu_mozart(m,nc,s,work);
  memset(lower,0,m->nnz*s*sizeof(double)); memset(upper,0,m->nnz*s*sizeof(double));
  for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) {
    int at = m->index[i*n+j]; if (at < 0) continue;
    for (size_t c = 0; c < nc; ++c) {
      if (i >= j) lower[(size_t)at*s+c] = i == j ? 1.0 : work[(size_t)at*s+c];
      if (i <= j) upper[(size_t)at*s+c] = work[(size_t)at*s+c];
    }
  }
  free(work); return CHEM_OK;
}
void chem_linear_solve_separate(const ChemModel *m, size_t nc, size_t s, const double *lower, const double *upper, double *x) {
  for (size_t i = 0; i < m->nvariables; ++i) {
    ChemRange range = m->lower[i];
    for (size_t j = range.first; j < range.first+range.count; ++j)
      for (size_t c = 0; c < nc; ++c) x[i*s+c] -= lower[m->lower_ids[j]*s+c]*x[m->lower_columns[j]*s+c];
  }
  for (size_t i = m->nvariables; i-- > 0;) {
    ChemRange range = m->upper[i];
    for (size_t j = range.first; j < range.first+range.count; ++j)
      for (size_t c = 0; c < nc; ++c) x[i*s+c] -= upper[m->upper_ids[j]*s+c]*x[m->upper_columns[j]*s+c];
    for (size_t c = 0; c < nc; ++c) x[i*s+c] /= upper[m->diagonal[i]*s+c];
  }
}
