#include "internal.h"

ChemState *chem_state_create(const ChemModel *m, size_t cells) {
  if (!m || !m->finalized || !cells || m->nvariables > SIZE_MAX / cells / (10*sizeof(double)) ||
      m->nnz > SIZE_MAX / cells / sizeof(double)) return NULL;
  ChemState *s = calloc(1, sizeof *s);
  if (!s) return NULL;
  s->model = m; s->cells = s->stride = cells; s->relative_tolerance = 1e-6;
  const size_t ny = m->nvariables * cells;
  s->variables = calloc(ny, sizeof(double));
  s->custom = calloc((m->ncustom ? m->ncustom : 1)*cells, sizeof(double));
  s->rates = calloc((m->nreactions ? m->nreactions : 1)*cells, sizeof(double));
  s->conditions = calloc(cells, sizeof(ChemConditions));
  s->absolute_tolerance = malloc(m->nvariables*sizeof(double));
  s->jacobian = calloc(m->nnz*cells, sizeof(double));
  s->work = calloc(ny*9, sizeof(double));
  if (!s->variables || !s->custom || !s->rates || !s->conditions || !s->absolute_tolerance || !s->jacobian || !s->work) {
    chem_state_free(s); return NULL;
  }
  s->new_y = s->work; s->forcing = s->work+ny; s->error = s->work+2*ny;
  for (size_t i = 0; i < 6; ++i) s->stages[i] = s->work+(i+3)*ny;
  for (size_t i = 0; i < m->nvariables; ++i) {
    const ChemSpecies *p = &m->species[m->variable_species[i]];
    s->absolute_tolerance[i] = p->absolute_tolerance;
    for (size_t c = 0; c < cells; ++c) s->variables[i*cells+c] = p->initial;
  }
  for (size_t i = 0; i < m->ncustom; ++i)
    for (size_t c = 0; c < cells; ++c) s->custom[i*cells+c] = m->custom_defaults[i];
  for (size_t c = 0; c < cells; ++c) {
    s->conditions[c].temperature = m->default_temperature;
    s->conditions[c].pressure = m->default_pressure;
    chem_ideal_air_density(&s->conditions[c]);
  }
  return s;
}
void chem_state_free(ChemState *s) {
  if (!s) return;
  free(s->variables); free(s->custom); free(s->rates); free(s->conditions);
  free(s->absolute_tolerance); free(s->jacobian); free(s->work); free(s);
}
ChemState *chem_state_copy(const ChemState *s) {
  if (!s) return NULL;
  ChemState *t = chem_state_create(s->model, s->cells);
  if (!t) return NULL;
  size_t n = s->cells; const ChemModel *m = s->model;
  memcpy(t->variables, s->variables, n*m->nvariables*sizeof(double));
  memcpy(t->custom, s->custom, n*m->ncustom*sizeof(double));
  memcpy(t->rates, s->rates, n*m->nreactions*sizeof(double));
  memcpy(t->conditions, s->conditions, n*sizeof(ChemConditions));
  memcpy(t->absolute_tolerance, s->absolute_tolerance, m->nvariables*sizeof(double));
  t->relative_tolerance = s->relative_tolerance;
  return t;
}
static ChemStatus set_values(ChemState *s, size_t index, double *array, const double *v, size_t count) {
  if (!s || !v) return CHEM_INVALID;
  if (index == SIZE_MAX) return CHEM_UNKNOWN;
  if (count != 1 && count != s->cells) return CHEM_SIZE;
  for (size_t c = 0; c < s->cells; ++c) array[index*s->stride+c] = v[count == 1 ? 0 : c];
  return CHEM_OK;
}
ChemStatus chem_state_set(ChemState *s, const char *name, const double *v, size_t count) {
  if (!s) return CHEM_INVALID;
  return set_values(s, chem_variable_id(s->model, name), s->variables, v, count);
}
ChemStatus chem_state_set_custom(ChemState *s, const char *name, const double *v, size_t count) {
  if (!s) return CHEM_INVALID;
  return set_values(s, chem_custom_id(s->model, name), s->custom, v, count);
}
ChemStatus chem_state_tolerances(ChemState *s, double relative, const double *absolute, size_t count) {
  if (!s || !(relative >= 0.0) || !absolute) return CHEM_INVALID;
  if (count != 1 && count != s->model->nvariables) return CHEM_SIZE;
  for (size_t i = 0; i < count; ++i) if (!(absolute[i] >= 0.0)) return CHEM_INVALID;
  s->relative_tolerance = relative;
  for (size_t i = 0; i < s->model->nvariables; ++i) s->absolute_tolerance[i] = absolute[count == 1 ? 0 : i];
  return CHEM_OK;
}
void chem_calculate_rates(ChemState *s) {
  const ChemModel *m = s->model;
  double *u = malloc((m->ncustom+1)*sizeof(double));
  if (!u) { for (size_t i = 0; i < m->nreactions*s->stride; ++i) s->rates[i] = NAN; return; }
  for (size_t c = 0; c < s->cells; ++c) {
    for (size_t i = 0; i < m->ncustom; ++i) u[i] = s->custom[i*s->stride+c];
    for (size_t r = 0; r < m->nreactions; ++r) {
      const ChemReaction *rx = &m->reactions[r];
      double k = chem_rate_value(&rx->rate, &s->conditions[c], u);
      for (size_t j = 0; j < rx->nreactants; ++j) {
        const ChemSpecies *p = &m->species[rx->reactants[j]];
        if (p->parameterized == 1) k *= s->conditions[c].air_density;
        else if (p->parameterized == 2 && p->parameterization) k *= p->parameterization(&s->conditions[c], u, p->context);
      }
      s->rates[r*s->stride+c] = k;
    }
  }
  free(u);
}

void chem_axpy(size_t n, double *restrict x, double a, const double *restrict y) {
  #pragma omp simd
  for (size_t i = 0; i < n; ++i) x[i] += a*y[i];
}
