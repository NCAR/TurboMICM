/* Translated numerical verification: Copyright (C) 2024-2026 University Corporation for Atmospheric Research.
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include "chem_verify.h"

ChemStatus chem_finite_difference(const ChemArray *v, size_t n, double perturbation,
                                  ChemArrayForcing function, void *context, ChemArray **result) {
  if (!v || !result || !function || n != v->columns || v->index || !(perturbation > 0.0) || n > SIZE_MAX/n) return CHEM_INVALID;
  ChemArray *r = chem_array_dense(v->rows,n*n,v->group,0);
  ChemArray *plus = chem_array_dense(v->rows,n,v->group,0), *minus = chem_array_dense(v->rows,n,v->group,0);
  ChemArray *fp = chem_array_dense(v->rows,n,v->group,0), *fm = chem_array_dense(v->rows,n,v->group,0);
  ChemArray *base = chem_array_dense(v->rows,n,v->group,0);
  if (!r || !plus || !minus || !fp || !fm || !base) {
    chem_array_free(r); chem_array_free(plus); chem_array_free(minus);
    chem_array_free(fp); chem_array_free(fm); chem_array_free(base); return CHEM_MEMORY;
  }
  for (size_t col = 0; col < n; ++col) {
    chem_array_copy(plus,v); chem_array_copy(minus,v);
    for (size_t b = 0; b < v->rows; ++b) {
      size_t at = chem_array_index(v,b,col);
      double x = v->values[at], h = perturbation*fmax(1.0,fabs(x));
      plus->values[at] = x+h; minus->values[at] = x-h;
    }
    chem_array_fill(fp,0); chem_array_fill(fm,0);
    function(plus,fp,context); function(minus,fm,context);
    for (size_t b = 0; b < v->rows; ++b) {
      double x = v->values[chem_array_index(v,b,col)], h = perturbation*fmax(1.0,fabs(x));
      int central = x-h >= 0;
      if (!central) { chem_array_fill(base,0); function(v,base,context); }
      for (size_t row = 0; row < n; ++row) {
        size_t at = chem_array_index(v,b,row);
        r->values[chem_array_index(r,b,row*n+col)] = (fp->values[at]-(central ? fm->values[at] : base->values[at]))/(central ? 2*h : h);
      }
    }
  }
  chem_array_free(plus); chem_array_free(minus); chem_array_free(fp); chem_array_free(fm); chem_array_free(base);
  *result = r; return CHEM_OK;
}
ChemJacobianComparison chem_compare_jacobian(const ChemArray *a, const ChemArray *fd, size_t n, double atol, double rtol) {
  ChemJacobianComparison r = {0}; r.passed = 1;
  if (!a || !fd || !a->index || a->order != n || fd->rows != a->rows || fd->columns != n*n) { r.passed = 0; return r; }
  for (size_t b = 0; b < a->rows; ++b) for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) {
    size_t at = chem_sparse_index(a,b,i,j); if (at == SIZE_MAX) continue;
    double analytical = -a->values[at], f = fd->values[chem_array_index(fd,b,i*n+j)];
    double error = fabs(analytical-f), scale = fmax(fabs(analytical),fabs(f));
    if (error > r.max_abs_error) {
      r.max_abs_error = error; r.worst_block = b; r.worst_row = i; r.worst_col = j; r.worst_analytical = analytical; r.worst_fd = f;
    }
    if (scale > 0 && error/scale > r.max_rel_error) r.max_rel_error = error/scale;
    if (!isfinite(error) || error > atol+rtol*scale) r.passed = 0;
  }
  return r;
}
ChemJacobianComparison chem_check_sparsity(const ChemArray *a, const ChemArray *fd, size_t n, double threshold) {
  ChemJacobianComparison r = {0}; r.passed = 1;
  if (!a || !fd || !a->index || a->order != n || fd->rows != a->rows || fd->columns != n*n) { r.passed = 0; return r; }
  for (size_t b = 0; b < a->rows; ++b) for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) {
    if (a->index[i*n+j] >= 0) continue;
    double f = fabs(fd->values[chem_array_index(fd,b,i*n+j)]);
    if (!isfinite(f) || f > threshold) {
      r.passed = 0;
      if (f > r.max_abs_error) { r.max_abs_error = f; r.worst_block = b; r.worst_row = i; r.worst_col = j; r.worst_fd = f; }
    }
  }
  return r;
}
