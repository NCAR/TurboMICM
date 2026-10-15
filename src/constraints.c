/* Linear and Van't Hoff equilibrium constraints, translated from MICM. */
#include "internal.h"

ck_status ck_constraint_validate(const ck_constraint *p)
{
  if (!p || !p->name || !p->nreactants || !p->reactants) return CK_INVALID;
  if (p->type == CK_LINEAR_CONSTRAINT) return CK_OK;
  if (p->type != CK_EQUILIBRIUM_CONSTRAINT || !p->nproducts || !p->products || !(p->equilibrium_ref > 0)) return CK_INVALID;
  for (size_t j = 0; j < p->nreactants; ++j) if (!(p->reactants[j].yield > 0)) return CK_INVALID;
  for (size_t j = 0; j < p->nproducts; ++j) if (!(p->products[j].yield > 0)) return CK_INVALID;
  return CK_OK;
}
static double product(const double *y, size_t cells, size_t cell, const ck_product *terms, size_t count, size_t exclude, double initial)
{
  for (size_t j = 0; j < count; ++j) if (j != exclude) {
    double conc = y[terms[j].species * cells + cell];
    initial *= pow(conc > 0 ? conc : 0, terms[j].yield);
  }
  return initial;
}
void ck_constraint_forcing(const ck_state *s, const double *y, double *f)
{
  const ck_model *m = s->model;
  for (size_t k = 0; k < m->nconstraints; ++k) {
    const ck_bound_constraint *b = &m->constraints[k];
    const ck_constraint *p = &b->description;
    for (size_t c = 0; c < s->cells; ++c) {
      double value;
      if (p->type == CK_LINEAR_CONSTRAINT) {
        value = 0;
        for (size_t j = 0; j < p->nreactants; ++j) value += b->reactants[j].yield * y[b->reactants[j].species * s->cells + c];
        value -= p->constant;
      } else {
        double rp = product(y, s->cells, c, b->reactants, p->nreactants, SIZE_MAX, s->custom[b->parameter * s->cells + c]);
        double pp = product(y, s->cells, c, b->products, p->nproducts, SIZE_MAX, 1);
        value = rp - pp;
      }
      f[b->row * s->cells + c] = value;
    }
  }
}
void ck_constraint_jacobian(const ck_state *s, const double *y, double *a)
{
  const ck_model *m = s->model;
  for (size_t k = 0; k < m->nconstraints; ++k) {
    const ck_bound_constraint *b = &m->constraints[k];
    const ck_constraint *p = &b->description;
    for (size_t side = 0; side < (p->type == CK_LINEAR_CONSTRAINT ? 1u : 2u); ++side) {
      const ck_product *terms = side == 0 ? b->reactants : b->products;
      size_t count = side == 0 ? p->nreactants : p->nproducts;
      for (size_t j = 0; j < count; ++j) {
        size_t off = m->map[b->row * m->n + terms[j].species];
        for (size_t c = 0; c < s->cells; ++c) {
          double partial = terms[j].yield;
          if (p->type == CK_EQUILIBRIUM_CONSTRAINT) {
            double prod = product(y, s->cells, c, terms, count, j, side == 0 ? s->custom[b->parameter * s->cells + c] : 1);
            double conc = y[terms[j].species * s->cells + c];
            partial = terms[j].yield == 1 ? prod : conc > 0 ? terms[j].yield * prod * pow(conc, terms[j].yield - 1) : 0;
          }
          if (side == 0) a[off * s->cells + c] -= partial;
          else a[off * s->cells + c] += partial;
        }
      }
    }
  }
}
static ck_status finite_algebraic(const ck_state *s, const double *v, double *maximum)
{
  *maximum = 0;
  int inf = 0;
  for (size_t j = 0; j < s->model->n; ++j) if (s->model->algebraic[j])
    for (size_t c = 0; c < s->cells; ++c) {
      double val = v[j * s->cells + c];
      if (isnan(val)) return CK_NAN;
      if (isinf(val)) inf = 1;
      if (fabs(val) > *maximum) *maximum = fabs(val);
    }
  return inf ? CK_INF : CK_OK;
}
static double weighted_correction(const ck_state *s, const double *y, const double *delta)
{
  double maximum = 0;
  for (size_t j = 0; j < s->model->n; ++j) if (s->model->algebraic[j])
    for (size_t c = 0; c < s->cells; ++c) {
      size_t i = j * s->cells + c;
      double corrected = y[i] + delta[i];
      double ymax = fabs(y[i]) > fabs(corrected) ? fabs(y[i]) : fabs(corrected);
      double scale = s->absolute_tolerance[j] + s->relative_tolerance * ymax;
      double q = scale > 0 ? fabs(delta[i]) / scale : delta[i] == 0 ? 0 : INFINITY;
      if (q > maximum) maximum = q;
    }
  return maximum;
}
static void apply_correction(const ck_state *s, double *candidate, const double *y, const double *delta, double step)
{
  memcpy(candidate, y, s->model->n * s->cells * sizeof(double));
  for (size_t j = 0; j < s->model->n; ++j) if (s->model->algebraic[j])
    for (size_t c = 0; c < s->cells; ++c) candidate[j * s->cells + c] += step * delta[j * s->cells + c];
}
ck_status ck_initialize_constraints(ck_state *s, const ck_rosenbrock_parameters *p, ck_stats *stats)
{
  if (!s || !p || !stats) return CK_INVALID;
  const ck_model *m = s->model;
  if (!m->nconstraints) return CK_OK;
  ck_workspace *w = s->workspace;
  size_t size = m->n * s->cells;
  double *delta = w->initial, *original = w->stage[0];
  int modified = 0;
  ck_status status = CK_CONSTRAINT_FAILED;
  for (size_t update = 0; update <= p->constraint_init_max_iterations; ++update) {
    memset(delta, 0, size * sizeof(double)); ck_constraint_forcing(s, s->y, delta);
    double maximum;
    status = finite_algebraic(s, delta, &maximum); if (status) goto fail;
    ++stats->constraint_init_iterations;
    if (maximum == 0) return CK_OK;
    memset(w->jacobian, 0, m->nnz * s->cells * sizeof(double));
    ck_constraint_jacobian(s, s->y, w->jacobian); ++stats->jacobian_updates;
    for (size_t j = 0; j < m->n; ++j) if (!m->algebraic[j])
      for (size_t c = 0; c < s->cells; ++c) w->jacobian[m->diagonal[j] * s->cells + c] = 1;
    ck_factor(m, s->cells, w->jacobian, p->lu); ++stats->decompositions;
    ck_linear_solve(m, s->cells, w->jacobian, delta); ++stats->solves;
    status = finite_algebraic(s, delta, &maximum); if (status) goto fail;
    double current_norm = weighted_correction(s, s->y, delta);
    if (current_norm <= p->constraint_init_tolerance) return CK_OK;
    status = CK_CONSTRAINT_FAILED;
    if (update == p->constraint_init_max_iterations) goto fail;
    if (!modified) { memcpy(original, s->y, size * sizeof(double)); modified = 1; }
    double step = 1;
    for (size_t backtrack = 0;; ++backtrack, step *= p->constraint_init_backtrack_factor) {
      if (p->constraint_init_max_backtracks && backtrack > p->constraint_init_max_backtracks) { status = CK_CONSTRAINT_FAILED; goto fail; }
      apply_correction(s, w->next, s->y, delta, step);
      status = finite_algebraic(s, w->next, &maximum);
      if (status) { if (!p->constraint_init_max_backtracks) goto fail; continue; }
      if (!p->constraint_init_max_backtracks) break;
      memset(w->error, 0, size * sizeof(double)); ck_constraint_forcing(s, w->next, w->error);
      if (finite_algebraic(s, w->error, &maximum)) continue;
      if (maximum == 0) break;
      ck_linear_solve(m, s->cells, w->jacobian, w->error); ++stats->solves;
      if (finite_algebraic(s, w->error, &maximum)) continue;
      double norm = weighted_correction(s, w->next, w->error);
      if (norm <= p->constraint_init_tolerance || norm < (1 - p->constraint_init_sufficient_decrease * step) * current_norm) break;
    }
    memcpy(s->y, w->next, size * sizeof(double));
  }
  status = CK_CONSTRAINT_FAILED;
fail:
  if (modified) memcpy(s->y, original, size * sizeof(double));
  return status;
}
