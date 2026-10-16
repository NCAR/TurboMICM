/* SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

const char *ch_version(void) { return "1.0.0.0"; }
int ch_version_component(int index)
{
  const int parts[] = {CH_VERSION_MAJOR,CH_VERSION_MINOR,CH_VERSION_PATCH,CH_VERSION_TWEAK};
  return index >= 0 && index < 4 ? parts[index] : CH_INVALID;
}
int ch_model_reaction_info(const ch_model *m, size_t index, ch_reaction_info *info)
{
  if (!m || !info) return CH_INVALID;
  if (index >= m->nr) return CH_SIZE;
  const ch_reaction *r = m->reactions+index;
  *info = (ch_reaction_info){r->nr,r->np,r->nparameterized,r->reactants,r->products,r->parameterized,r->yields,&r->rate};
  return CH_OK;
}
static int type_order(int type)
{
  const int order[] = {0,1,2,3,5,4,9,8,6,7,10};
  return type >= CH_ARRHENIUS && type <= CH_LAMBDA ? order[type] : -1;
}
int ch_model_rate_group(const ch_model *m, int type, size_t *offset, size_t *count)
{
  int rank = type_order(type);
  if (!m || !offset || !count || rank < 0) return CH_INVALID;
  *offset = *count = 0;
  for (size_t i = 0; i < m->nr; ++i) {
    int r = type_order(m->reactions[i].rate.type);
    *offset += r < rank; *count += r == rank;
  }
  return CH_OK;
}
const ch_rate *ch_model_rate_in_group(const ch_model *m, int type, size_t index)
{
  if (!m) return NULL;
  for (size_t i = 0; i < m->nr; ++i) if (m->reactions[i].rate.type == type)
    if (!index--) return &m->reactions[i].rate;
  return NULL;
}
int ch_model_species_used(const ch_model *m, size_t id)
{
  if (!m || id >= m->ns) return CH_NOT_FOUND;
  if (m->finalized) return m->species[id].used;
  for (size_t i = 0; i < m->nr; ++i) {
    const ch_reaction *r = m->reactions+i;
    for (size_t j = 0; j < r->nr; ++j) if ((size_t)r->reactants[j] == id) return 1;
    for (size_t j = 0; j < r->np; ++j) if ((size_t)r->products[j] == id) return 1;
  }
  return 0;
}
int ch_model_kinetic_pattern(const ch_model *m, unsigned char *pattern, size_t size)
{
  if (!m || !m->finalized || !pattern) return CH_INVALID;
  if (size != m->nv*m->nv) return CH_SIZE;
  memset(pattern,0,size);
  for (size_t i = 0; i < m->nr; ++i) {
    const ch_reaction *r = m->reactions+i;
    for (size_t j = 0; j < r->nr; ++j) {
      for (size_t k = 0; k < r->nr; ++k) pattern[(size_t)r->reactants[k]*m->nv+r->reactants[j]] = 1;
      for (size_t k = 0; k < r->np; ++k) pattern[(size_t)r->products[k]*m->nv+r->reactants[j]] = 1;
    }
  }
  return CH_OK;
}
int ch_add_force(const ch_state *s, const double *y, double *f)
{
  if (!s || !y || !f) return CH_INVALID;
  size_t n = s->model->nv*s->cells;
  double *tmp = ch_allocate(n,sizeof(double)); if (!tmp) return CH_MEMORY;
  ch_force(s,y,tmp);
  for (size_t i = 0; i < n; ++i) f[i] += tmp[i];
  free(tmp); return CH_OK;
}
int ch_add_negative_jacobian(const ch_state *s, const double *y, double *a)
{
  if (!s || !y || !a) return CH_INVALID;
  size_t n = s->model->sparse->nnz*s->cells;
  double *tmp = ch_allocate(n,sizeof(double)); if (!tmp) return CH_MEMORY;
  ch_jacobian(s,y,tmp);
  for (size_t i = 0; i < n; ++i) a[i] += tmp[i];
  free(tmp); return CH_OK;
}
/* Central differences and one-sided fallback follow MICM's verification helper.
 * Evaluate f(y) once when needed rather than once for each cell near zero. */
int ch_finite_difference(size_t n, size_t nc, const double *y, double step,
                         ch_forcing_function fun, void *ctx, double *jac)
{
  if (!n || !nc || !y || !jac || !fun || !isfinite(step) || step <= 0) return CH_INVALID;
  if (n > SIZE_MAX/n || nc > SIZE_MAX/(n*n) || nc > SIZE_MAX/n/5) return CH_SIZE;
  size_t size = n*nc;
  double *storage = ch_allocate(5*size,sizeof(double)); if (!storage) return CH_MEMORY;
  double *plus = storage, *minus = plus+size, *fp = minus+size, *fm = fp+size, *base = fm+size;
  int has_base = 0;
  for (size_t col = 0; col < n; ++col) {
    memcpy(plus,y,size*sizeof(double)); memcpy(minus,y,size*sizeof(double));
    memset(fp,0,size*sizeof(double)); memset(fm,0,size*sizeof(double));
    for (size_t c = 0; c < nc; ++c) {
      double x = y[col*nc+c], h = step*fmax(1,fabs(x));
      plus[col*nc+c] = x+h; minus[col*nc+c] = x-h;
    }
    fun(n,nc,plus,fp,ctx); fun(n,nc,minus,fm,ctx);
    for (size_t c = 0; c < nc; ++c) {
      double x = y[col*nc+c], h = step*fmax(1,fabs(x));
      int central = x-h >= 0;
      if (!central && !has_base) { fun(n,nc,y,base,ctx); has_base = 1; }
      for (size_t row = 0; row < n; ++row)
        jac[(row*n+col)*nc+c] = central ? (fp[row*nc+c]-fm[row*nc+c])/(2*h) : (fp[row*nc+c]-base[row*nc+c])/h;
    }
  }
  free(storage); return CH_OK;
}
ch_jacobian_comparison ch_compare_jacobian(size_t n, size_t nc, const unsigned char *pattern,
                                           const double *a, const double *fd, double atol, double rtol)
{
  ch_jacobian_comparison result = {0}; result.passed = 1;
  if (!pattern || !a || !fd || atol < 0 || rtol < 0) { result.passed = 0; return result; }
  for (size_t c = 0; c < nc; ++c) for (size_t r = 0; r < n; ++r) for (size_t j = 0; j < n; ++j) {
    if (!pattern[r*n+j]) continue;
    size_t i = (r*n+j)*nc+c;
    double analytical = -a[i], finite_difference = fd[i], delta = fabs(analytical-finite_difference);
    double scale = fmax(fabs(analytical),fabs(finite_difference));
    double relative = scale > 0 ? delta/scale : 0;
    if (delta > result.max_absolute_error) {
      result.max_absolute_error = delta; result.worst_cell = c; result.worst_row = r; result.worst_column = j;
      result.worst_analytical = analytical; result.worst_finite_difference = finite_difference;
    }
    if (relative > result.max_relative_error) result.max_relative_error = relative;
    if (!isfinite(delta) || delta > atol+rtol*scale) result.passed = 0;
  }
  return result;
}
ch_jacobian_comparison ch_check_sparsity(size_t n, size_t nc, const unsigned char *pattern,
                                       const double *fd, double threshold)
{
  ch_jacobian_comparison result = {0}; result.passed = 1;
  if (!pattern || !fd || threshold < 0) { result.passed = 0; return result; }
  for (size_t c = 0; c < nc; ++c) for (size_t r = 0; r < n; ++r) for (size_t j = 0; j < n; ++j) {
    if (pattern[r*n+j]) continue;
    double value = fabs(fd[(r*n+j)*nc+c]);
    if (!isfinite(value)) result.passed = 0;
    if (value > threshold) {
      result.passed = 0;
      if (value > result.max_absolute_error) {
        result.max_absolute_error = value; result.worst_cell = c; result.worst_row = r; result.worst_column = j;
        result.worst_finite_difference = value;
      }
    }
  }
  return result;
}
