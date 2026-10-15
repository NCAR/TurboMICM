#include "internal.h"

ck_status ck_copy_description(const ck_mechanism *d, ck_model *m)
{
  size_t nspecies = d->nspecies;
  for (size_t k = 0; k < d->nexternal; ++k)
  {
    const ck_external_model *e = &d->external[k];
    if ((e->nspecies && !e->species) || (e->ncustom && !e->custom_names) ||
        (e->nalgebraic && !e->algebraic_names) || (e->njacobian && !e->jacobian) ||
        (e->clone_context && !e->destroy_context))
      return CK_INVALID;
    if (e->nspecies > SIZE_MAX - nspecies || e->ncustom > SIZE_MAX - m->extra_custom)
      return CK_INVALID;
    nspecies += e->nspecies;
    m->extra_custom += e->ncustom;
  }
  m->owned_description = ck_allocate(1, sizeof(ck_mechanism));
  m->owned_species = ck_allocate(nspecies, sizeof(ck_species));
  m->external = ck_allocate(d->nexternal, sizeof(ck_external_model));
  if (!m->owned_description || !m->owned_species || !m->external)
    return CK_MEMORY;
  *m->owned_description = *d;
  m->description = m->owned_description;
  m->owned_description->nspecies = nspecies;
  m->owned_description->species = m->owned_species;
  if (d->nspecies)
    memcpy(m->owned_species, d->species, d->nspecies * sizeof(ck_species));
  size_t offset = d->nspecies;
  for (size_t k = 0; k < d->nexternal; ++k)
  {
    m->external[k] = d->external[k];
    ck_external_model *e = &m->external[k];
    if (e->clone_context)
    {
      e->context = e->clone_context(d->external[k].context);
      if (!e->context)
        return CK_MEMORY;
    }
    ++m->nexternal;
    if (e->nspecies)
      memcpy(m->owned_species + offset, e->species, e->nspecies * sizeof(ck_species));
    offset += e->nspecies;
  }
  return CK_OK;
}
ck_status ck_external_pattern(const ck_model *m, unsigned char *pattern, int constraints)
{
  for (size_t k = 0; k < m->nexternal; ++k)
  {
    const ck_external_model *e = &m->external[k];
    for (size_t j = 0; j < e->njacobian; ++j)
    {
      const ck_named_jacobian_entry *entry = &e->jacobian[j];
      if (entry->constraint && (!constraints || !e->nalgebraic))
        continue;
      size_t row, col;
      if (ck_species_index(m, entry->row, &row) || ck_species_index(m, entry->column, &col))
        return CK_UNKNOWN;
      pattern[row * m->n + col] = 1;
    }
  }
  return CK_OK;
}
ck_status ck_external_bind(ck_model *m)
{
  for (size_t k = 0; k < m->nexternal; ++k)
  {
    const ck_external_model *e = &m->external[k];
    for (size_t j = 0; j < e->ncustom; ++j)
    {
      size_t index;
      if (!e->custom_names[j])
        return CK_INVALID;
      if (ck_custom_index(m, e->custom_names[j], &index) == CK_UNKNOWN)
      {
        char *name = ck_string_copy(e->custom_names[j]);
        if (!name)
          return CK_MEMORY;
        m->custom_names[m->nc++] = name;
      }
    }
    for (size_t j = 0; j < e->nalgebraic; ++j)
    {
      size_t row;
      if (ck_species_index(m, e->algebraic_names[j], &row))
        return CK_UNKNOWN;
      if (m->algebraic[row])
        return CK_INVALID;
      m->algebraic[row] = 1;
      ++m->nconstraints;
    }
  }
  return CK_OK;
}
ck_status ck_external_finalize(ck_model *m)
{
  for (size_t k = 0; k < m->nexternal; ++k)
    if (m->external[k].finalize)
    {
      ck_status status = m->external[k].finalize(m, m->external[k].context);
      if (status)
        return status;
    }
  return CK_OK;
}
ck_status ck_external_initialize(ck_state *s)
{
  if (!s)
    return CK_INVALID;
  for (size_t k = 0; k < s->model->nexternal; ++k)
  {
    const ck_external_model *e = &s->model->external[k];
    if (e->nalgebraic && e->initialize_constraints)
    {
      ck_status status = e->initialize_constraints(s, e->context);
      if (status)
        return status;
    }
  }
  return CK_OK;
}
