#include "internal.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void *array(size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size) return NULL;
    return calloc(n ? n : 1, size);
}

static ChemState *create_state(const ChemMechanism *m, size_t cells, const ChemState *model)
{
    if (!m || !cells || !m->species_count || !m->species ||
        (m->reaction_count && (!m->reactions || !m->reactants || !m->products)) ||
        (m->parameter_count && (!m->parameter_names || !m->parameter_defaults))) return NULL;
    if (m->species_count > SIZE_MAX / cells || m->reaction_count > SIZE_MAX / cells ||
        m->parameter_count > SIZE_MAX / cells || m->species_count * cells > SIZE_MAX / 6) return NULL;
    ChemState *s = calloc(1, sizeof(*s));
    if (!s) return NULL;
    s->mechanism = m;
    s->cells = cells;
    if (model) {
        s->nonzeros = model->nonzeros;
        s->index = model->index; s->columns = model->columns; s->row_start = model->row_start;
        s->pivots = model->pivots; s->lower = model->lower;
        s->upper = model->upper; s->updates = model->updates;
        s->symbolic_references = model->symbolic_references;
        atomic_fetch_add(s->symbolic_references, 1);
    } else {
        s->symbolic_references = array(1, sizeof(*s->symbolic_references));
        if (!s->symbolic_references) { free(s); return NULL; }
        atomic_init(s->symbolic_references, 1);
        if (chem_make_symbolic(s) != CHEM_OK) { chem_state_free(s); return NULL; }
    }
    size_t ny = m->species_count * cells;
    s->values = array(ny, sizeof(double));
    s->parameters = array(m->parameter_count * cells, sizeof(double));
    s->temperatures = array(cells, sizeof(double));
    s->pressures = array(cells, sizeof(double));
    s->rates = array(m->reaction_count * cells, sizeof(double));
    s->scratch = array(cells, sizeof(double));
    s->forcing = array(ny, sizeof(double));
    s->candidate = array(ny, sizeof(double));
    s->error = array(ny, sizeof(double));
    s->stages = array(6 * ny, sizeof(double));
    if (!s->values || !s->parameters || !s->temperatures || !s->pressures ||
        !s->rates || !s->scratch || !s->forcing || !s->candidate || !s->error || !s->stages ||
        s->nonzeros > SIZE_MAX / cells ||
        !(s->matrix = array(s->nonzeros * cells, sizeof(double)))) {
        chem_state_free(s);
        return NULL;
    }
    for (size_t i = 0; i < m->species_count; ++i)
        for (size_t c = 0; c < cells; ++c) s->values[i*cells+c] = m->species[i].initial;
    for (size_t i = 0; i < m->parameter_count; ++i)
        for (size_t c = 0; c < cells; ++c) s->parameters[i*cells+c] = m->parameter_defaults[i];
    for (size_t c = 0; c < cells; ++c) {
        s->temperatures[c] = m->initial_temperature;
        s->pressures[c] = m->initial_pressure;
    }
    return s;
}

ChemState *chem_state_create(const ChemMechanism *m, size_t cells)
{ return create_state(m, cells, NULL); }

ChemState *chem_state_create_like(const ChemState *s, size_t cells)
{ return s ? create_state(s->mechanism, cells, s) : NULL; }

void chem_state_free(ChemState *s)
{
    if (!s) return;
    free(s->values); free(s->parameters); free(s->temperatures); free(s->pressures);
    free(s->rates); free(s->matrix); free(s->scratch); free(s->forcing);
    free(s->candidate); free(s->error); free(s->stages);
    if (s->symbolic_references && atomic_fetch_sub(s->symbolic_references, 1) == 1) {
        free(s->index); free(s->columns); free(s->row_start); free(s->pivots);
        free(s->lower); free(s->upper); free(s->updates); free(s->symbolic_references);
    }
    free(s);
}

size_t chem_cell_count(const ChemState *s) { return s ? s->cells : 0; }
double *chem_species_values(ChemState *s, size_t i)
{ return s && i < s->mechanism->species_count ? s->values + i*s->cells : NULL; }
double *chem_parameter_values(ChemState *s, size_t i)
{ return s && i < s->mechanism->parameter_count ? s->parameters + i*s->cells : NULL; }
double *chem_temperatures(ChemState *s) { return s ? s->temperatures : NULL; }
double *chem_pressures(ChemState *s) { return s ? s->pressures : NULL; }

int chem_find_species(const ChemMechanism *m, const char *name)
{
    if (!m || !name) return -1;
    for (size_t i = 0; i < m->species_count; ++i)
        if (!strcmp(name, m->species[i].name)) return (int)i;
    return -1;
}
int chem_find_parameter(const ChemMechanism *m, const char *name)
{
    if (!m || !name) return -1;
    for (size_t i = 0; i < m->parameter_count; ++i)
        if (!strcmp(name, m->parameter_names[i])) return (int)i;
    return -1;
}

const char *chem_error_string(int code)
{
    switch (code) {
    case CHEM_OK: return "success";
    case CHEM_INVALID: return "invalid argument or mechanism";
    case CHEM_MEMORY: return "out of memory";
    case CHEM_SINGULAR: return "singular linear system";
    case CHEM_TOO_MANY_STEPS: return "maximum internal steps exceeded";
    case CHEM_STEP_TOO_SMALL: return "internal step size too small";
    case CHEM_NONFINITE: return "non-finite numerical result";
    default: return "unknown error";
    }
}
