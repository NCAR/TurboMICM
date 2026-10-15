#include "internal.h"

void ck_state_free(ck_state *s)
{
  if (!s)
    return;
  ck_workspace *w = s->workspace;
  if (w)
  {
    free(w->initial);
    free(w->next);
    free(w->error);
    free(w->jacobian);
    free(w->scratch);
    for (size_t j = 0; j < 6; ++j)
      free(w->stage[j]);
    free(w);
  }
  free(s->y);
  free(s->rates);
  free(s->custom);
  free(s->conditions);
  free(s->absolute_tolerance);
  free(s);
}
ck_status ck_state_create(const ck_model *m, size_t cells, ck_state **out)
{
  if (!out)
    return CK_INVALID;
  *out = NULL;
  if (!m || !m->description || !cells || cells > SIZE_MAX / m->n || cells > SIZE_MAX / m->nnz ||
      (m->nr && cells > SIZE_MAX / m->nr) || (m->nc && cells > SIZE_MAX / m->nc))
    return CK_INVALID;
  ck_state *s = ck_allocate(1, sizeof *s);
  if (!s)
    return CK_MEMORY;
  s->model = m;
  s->cells = cells;
  s->relative_tolerance = 1e-6;
  size_t size = cells * m->n;
  s->y = ck_allocate(size, sizeof(double));
  s->rates = ck_allocate(cells * m->nr, sizeof(double));
  s->custom = ck_allocate(cells * m->nc, sizeof(double));
  s->conditions = ck_allocate(cells, sizeof(ck_conditions));
  s->absolute_tolerance = ck_allocate(m->n, sizeof(double));
  ck_workspace *w = ck_allocate(1, sizeof *w);
  s->workspace = w;
  if (!s->y || !s->rates || !s->custom || !s->conditions || !s->absolute_tolerance || !w)
    goto fail;
  w->scratch = ck_allocate(cells, sizeof(double));
  w->initial = ck_allocate(size, sizeof(double));
  w->next = ck_allocate(size, sizeof(double));
  w->error = ck_allocate(size, sizeof(double));
  w->jacobian = ck_allocate(cells * m->nnz, sizeof(double));
  if (!w->scratch || !w->initial || !w->next || !w->error || !w->jacobian)
    goto fail;
  for (size_t j = 0; j < 6; ++j)
  {
    w->stage[j] = ck_allocate(size, sizeof(double));
    if (!w->stage[j])
      goto fail;
  }
  for (size_t j = 0; j < m->n; ++j)
  {
    const ck_species *sp = &m->description->species[m->variable_to_species[j]];
    s->absolute_tolerance[j] = sp->absolute_tolerance > 0 ? sp->absolute_tolerance : 1e-3;
    for (size_t c = 0; c < cells; ++c)
      s->y[j * cells + c] = sp->initial;
  }
  for (size_t c = 0; c < cells; ++c)
  {
    s->conditions[c].temperature = m->description->temperature;
    s->conditions[c].pressure = m->description->pressure;
    if (s->conditions[c].temperature > 0)
      ck_ideal_air_density(&s->conditions[c]);
  }
  for (size_t j = 0; j < m->description->ncustom_defaults; ++j)
  {
    size_t index;
    if (ck_custom_index(m, m->description->custom_defaults[j].name, &index))
    {
      ck_state_free(s);
      return CK_UNKNOWN;
    }
    for (size_t c = 0; c < cells; ++c)
      s->custom[index * cells + c] = m->description->custom_defaults[j].value;
  }
  *out = s;
  return CK_OK;
fail:
  ck_state_free(s);
  return CK_MEMORY;
}
ck_status ck_state_clone(const ck_state *s, ck_state **out)
{
  if (!s)
    return CK_INVALID;
  ck_status status = ck_state_create(s->model, s->cells, out);
  if (status)
    return status;
  ck_state *d = *out;
  memcpy(d->y, s->y, s->cells * s->model->n * sizeof(double));
  memcpy(d->rates, s->rates, s->cells * s->model->nr * sizeof(double));
  memcpy(d->custom, s->custom, s->cells * s->model->nc * sizeof(double));
  memcpy(d->conditions, s->conditions, s->cells * sizeof(ck_conditions));
  memcpy(d->absolute_tolerance, s->absolute_tolerance, s->model->n * sizeof(double));
  d->relative_tolerance = s->relative_tolerance;
  return CK_OK;
}
ck_status ck_state_set(ck_state *s, const char *name, const double *v, size_t count)
{
  if (!s || !v || count != s->cells)
    return CK_INVALID;
  size_t i;
  ck_status error = ck_species_index(s->model, name, &i);
  if (error)
    return error;
  for (size_t c = 0; c < s->cells; ++c)
    s->y[i * s->cells + c] = v[count == 1 ? 0 : c];
  return CK_OK;
}
ck_status ck_state_set_custom(ck_state *s, const char *name, const double *v, size_t count)
{
  if (!s || !v || count != s->cells)
    return CK_INVALID;
  size_t i;
  ck_status error = ck_custom_index(s->model, name, &i);
  if (error)
    return error;
  for (size_t c = 0; c < s->cells; ++c)
    s->custom[i * s->cells + c] = v[count == 1 ? 0 : c];
  return CK_OK;
}
ck_status ck_state_tolerances(ck_state *s, double relative, const double *absolute, size_t count)
{
  if (!s || !absolute || count != s->model->n || !(relative >= 0) || !isfinite(relative))
    return CK_INVALID;
  for (size_t i = 0; i < count; ++i)
    if (!(absolute[i] > 0) || !isfinite(absolute[i]))
      return CK_INVALID;
  s->relative_tolerance = relative;
  memcpy(s->absolute_tolerance, absolute, count * sizeof(double));
  return CK_OK;
}
ck_status ck_update_rates(ck_state *s)
{
  if (!s)
    return CK_INVALID;
  const ck_model *m = s->model;
  for (size_t r = 0; r < m->nr; ++r)
  {
    const ck_rate *rate = &m->rates[r];
    size_t off = m->custom_offset[r];
    for (size_t c = 0; c < s->cells; ++c)
    {
      double custom[2] = {0, 0};
      if (rate->type == CK_USER || rate->type == CK_SURFACE)
        custom[0] = s->custom[off * s->cells + c];
      if (rate->type == CK_SURFACE)
        custom[1] = s->custom[(off + 1) * s->cells + c];
      double k = ck_rate_calculate(rate, &s->conditions[c], custom);
      const ck_reaction *rx = &m->description->reactions[m->reaction_order[r]];
      for (size_t j = 0; j < rx->nreactants; ++j)
      {
        const ck_species *sp = &m->description->species[rx->reactants[j]];
        if (sp->parameterized == 1)
          k *= s->conditions[c].air_density;
        else if (sp->parameterized == 2 && sp->function)
          k *= sp->function(&s->conditions[c], NULL, sp->context);
      }
      s->rates[r * s->cells + c] = k;
    }
  }
  for (size_t k = 0; k < m->nbuiltinconstraints; ++k)
  {
    const ck_bound_constraint *b = &m->constraints[k];
    if (b->description.type != CK_EQUILIBRIUM_CONSTRAINT)
      continue;
    const ck_constraint *p = &b->description;
    double gas = p->gas_constant == 0 ? CK_GAS_CONSTANT : p->gas_constant;
    double tref = p->reference_temperature == 0 ? 298.15 : p->reference_temperature;
    for (size_t c = 0; c < s->cells; ++c)
      s->custom[b->parameter * s->cells + c] =
          p->equilibrium_ref * exp((p->enthalpy / gas) * (1.0 / s->conditions[c].temperature - 1.0 / tref));
  }

  for (size_t k = 0; k < m->nexternal; ++k)
  {
    const ck_external_model *e = &m->external[k];
    if (e->update)
    {
      ck_status status = e->update(s, e->context);
      if (status)
        return status;
    }
    if (e->nalgebraic && e->update_constraints)
    {
      ck_status status = e->update_constraints(s, e->context);
      if (status)
        return status;
    }
  }
  return CK_OK;
}
ck_status ck_state_set_many(ck_state *s, const ck_named_values *v, size_t count, int custom)
{
  if (!s || (!v && count))
    return CK_INVALID;
  for (size_t j = 0; j < count; ++j)
  {
    ck_status error = custom ? ck_state_set_custom(s, v[j].name, v[j].values, v[j].count)
                             : ck_state_set(s, v[j].name, v[j].values, v[j].count);
    if (error)
      return error;
  }
  return CK_OK;
}
ck_status ck_state_set_custom_rows(ck_state *s, const double *v, size_t rows, size_t columns)
{
  if (!s || !v || rows != s->cells || columns != s->model->nc)
    return CK_INVALID;
  for (size_t j = 0; j < columns; ++j)
    for (size_t c = 0; c < rows; ++c)
      s->custom[j * rows + c] = v[c * columns + j];
  return CK_OK;
}
ck_status ck_state_move(ck_state **source, ck_state **destination)
{
  if (!source || !destination || source == destination || *destination)
    return CK_INVALID;
  *destination = *source;
  *source = NULL;
  return CK_OK;
}
