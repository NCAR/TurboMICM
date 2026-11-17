#include "internal.h"

#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void *array(size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size) return NULL;
    return calloc(n ? n : 1, size);
}

static ChemState *create_state(const ChemMechanism *m, size_t cells, const ChemConstraint *constraints,
                               size_t constraint_count, const ChemState *model)
{
    if (!m || !cells || !m->species_count || !m->species ||
        (m->reaction_count && !m->reactions) ||
        (m->parameter_count && (!m->parameter_names || !m->parameter_defaults))) return NULL;
    if (m->species_count > SIZE_MAX / cells || m->reaction_count > SIZE_MAX / cells ||
        m->parameter_count > SIZE_MAX / cells || constraint_count > SIZE_MAX / cells ||
        m->species_count * cells > SIZE_MAX / 6 ||
        (!model && chem_validate_constraints(m, constraints, constraint_count) != CHEM_OK)) return NULL;
    ChemState *s = calloc(1, sizeof(*s));
    if (!s) return NULL;
    s->mechanism = m;
    s->cells = cells;
    if (model) {
        s->constraints = model->constraints; s->constraint_count = model->constraint_count;
        s->nonzeros = model->nonzeros;
        s->index = model->index; s->columns = model->columns; s->row_start = model->row_start;
        s->pivots = model->pivots; s->lower = model->lower;
        s->upper = model->upper; s->updates = model->updates; s->algebraic = model->algebraic;
        s->symbolic_references = model->symbolic_references;
        atomic_fetch_add(s->symbolic_references, 1);
    } else {
        s->constraints = constraints; s->constraint_count = constraint_count;
        s->symbolic_references = array(1, sizeof(*s->symbolic_references));
        if (!s->symbolic_references) { free(s); return NULL; }
        atomic_init(s->symbolic_references, 1);
        s->algebraic = array(m->species_count, sizeof(*s->algebraic));
        if (!s->algebraic) { chem_state_free(s); return NULL; }
        for (size_t i=0; i<constraint_count; ++i) s->algebraic[constraints[i].algebraic_species] = 1;
        if (chem_make_symbolic(s) != CHEM_OK) { chem_state_free(s); return NULL; }
    }
    size_t ny = m->species_count * cells;
    s->values = array(ny, sizeof(double));
    s->parameters = array(m->parameter_count * cells, sizeof(double));
    s->temperatures = array(cells, sizeof(double));
    s->pressures = array(cells, sizeof(double));
    s->air_densities = array(cells, sizeof(double));
    s->explicit_air_density = array(cells, sizeof(*s->explicit_air_density));
    s->absolute_tolerances = array(m->species_count, sizeof(double));
    s->constraint_constants = array(s->constraint_count * cells, sizeof(double));
    s->constraint_gradient = array(m->species_count, sizeof(double));
    s->rates = array(m->reaction_count * cells, sizeof(double));
    s->rate_functions = array(m->reaction_count, sizeof(*s->rate_functions));
    s->rate_contexts = array(m->reaction_count, sizeof(*s->rate_contexts));
    s->scratch = array(cells, sizeof(double));
    s->forcing = array(ny, sizeof(double));
    s->candidate = array(ny, sizeof(double));
    s->error = array(ny, sizeof(double));
    s->stages = array(6 * ny, sizeof(double));
    if (!s->values || !s->parameters || !s->temperatures || !s->pressures ||
        !s->air_densities || !s->explicit_air_density || !s->absolute_tolerances ||
        !s->constraint_constants || !s->constraint_gradient ||
        !s->rates || !s->rate_functions || !s->rate_contexts || !s->scratch || !s->forcing ||
        !s->candidate || !s->error || !s->stages ||
        s->nonzeros > SIZE_MAX / cells ||
        !(s->matrix = array(s->nonzeros * cells, sizeof(double)))) {
        chem_state_free(s);
        return NULL;
    }
    for (size_t i = 0; i < m->species_count; ++i) {
        s->absolute_tolerances[i] = m->species[i].absolute_tolerance;
        for (size_t c = 0; c < cells; ++c) s->values[i*cells+c] = m->species[i].initial;
    }
    for (size_t i = 0; i < m->parameter_count; ++i)
        for (size_t c = 0; c < cells; ++c) s->parameters[i*cells+c] = m->parameter_defaults[i];
    for (size_t c = 0; c < cells; ++c) {
        s->temperatures[c] = m->initial_temperature;
        s->pressures[c] = m->initial_pressure;
    }
    for (size_t i = 0; i < m->reaction_count; ++i) {
        s->rate_functions[i] = m->reactions[i].rate.function;
        s->rate_contexts[i] = m->reactions[i].rate.context;
    }
    return s;
}

ChemState *chem_state_create(const ChemMechanism *m, size_t cells)
{ return create_state(m, cells, NULL, 0, NULL); }

ChemState *chem_state_create_constrained(const ChemMechanism *m, const ChemConstraint *constraints,
                                         size_t count, size_t cells)
{ return create_state(m, cells, constraints, count, NULL); }

ChemState *chem_state_create_like(const ChemState *s, size_t cells)
{ return s ? create_state(s->mechanism, cells, NULL, s->constraint_count, s) : NULL; }

ChemState *chem_state_clone(const ChemState *s)
{
    if (!s) return NULL;
    ChemState *copy=chem_state_create_like(s,s->cells);
    if (!copy) return NULL;
    size_t cells=s->cells, species=s->mechanism->species_count, parameters=s->mechanism->parameter_count;
    size_t reactions=s->mechanism->reaction_count;
    memcpy(copy->values,s->values,species*cells*sizeof(*s->values));
    memcpy(copy->parameters,s->parameters,parameters*cells*sizeof(*s->parameters));
    memcpy(copy->temperatures,s->temperatures,cells*sizeof(*s->temperatures));
    memcpy(copy->pressures,s->pressures,cells*sizeof(*s->pressures));
    memcpy(copy->air_densities,s->air_densities,cells*sizeof(*s->air_densities));
    memcpy(copy->explicit_air_density,s->explicit_air_density,cells*sizeof(*s->explicit_air_density));
    memcpy(copy->absolute_tolerances,s->absolute_tolerances,species*sizeof(*s->absolute_tolerances));
    memcpy(copy->constraint_constants,s->constraint_constants,s->constraint_count*cells*sizeof(*s->constraint_constants));
    memcpy(copy->rates,s->rates,reactions*cells*sizeof(*s->rates));
    memcpy(copy->rate_functions,s->rate_functions,reactions*sizeof(*s->rate_functions));
    memcpy(copy->rate_contexts,s->rate_contexts,reactions*sizeof(*s->rate_contexts));
    copy->linear_solver=s->linear_solver; copy->linear_solver_context=s->linear_solver_context;
    return copy;
}

void chem_state_free(ChemState *s)
{
    if (!s) return;
    free(s->values); free(s->parameters); free(s->temperatures); free(s->pressures);
    free(s->air_densities); free(s->explicit_air_density); free(s->absolute_tolerances);
    free(s->rates); free(s->rate_functions); free(s->rate_contexts);
    free(s->matrix); free(s->scratch); free(s->forcing);
    free(s->candidate); free(s->error); free(s->stages); free(s->constraint_constants); free(s->constraint_gradient);
    if (s->symbolic_references && atomic_fetch_sub(s->symbolic_references, 1) == 1) {
        free(s->index); free(s->columns); free(s->row_start); free(s->pivots);
        free(s->lower); free(s->upper); free(s->updates); free(s->algebraic); free(s->symbolic_references);
    }
    free(s);
}

size_t chem_cell_count(const ChemState *s) { return s ? s->cells : 0; }
double *chem_species_values(ChemState *s, size_t i)
{ return s && i < s->mechanism->species_count ? s->values + i*s->cells : NULL; }
const double *chem_species_values_const(const ChemState *s, size_t i)
{ return s && i < s->mechanism->species_count ? s->values + i*s->cells : NULL; }
double *chem_parameter_values(ChemState *s, size_t i)
{ return s && i < s->mechanism->parameter_count ? s->parameters + i*s->cells : NULL; }
double *chem_temperatures(ChemState *s) { return s ? s->temperatures : NULL; }
double *chem_pressures(ChemState *s) { return s ? s->pressures : NULL; }
double *chem_absolute_tolerances(ChemState *s) { return s ? s->absolute_tolerances : NULL; }
const double *chem_rate_constants(const ChemState *s, size_t i)
{ return s && i < s->mechanism->reaction_count ? s->rates+i*s->cells : NULL; }
double *chem_constraint_constants(ChemState *s, size_t i)
{ return s && i < s->constraint_count ? s->constraint_constants+i*s->cells : NULL; }
int chem_find_constraint(const ChemState *s, const char *name)
{
    if (!s || !name) return -1;
    for (size_t i=0;i<s->constraint_count;++i)
        if (s->constraints[i].name && !strcmp(s->constraints[i].name,name)) return (int)i;
    return -1;
}
int chem_is_algebraic(const ChemState *s, size_t species)
{ return s && species<s->mechanism->species_count ? !!s->algebraic[species] : -1; }
int chem_set_air_density(ChemState *s, size_t cell, double density)
{
    if (!s || cell >= s->cells || !isfinite(density) || density < 0) return CHEM_INVALID;
    s->air_densities[cell] = density; s->explicit_air_density[cell] = 1;
    return CHEM_OK;
}
int chem_use_ideal_air_density(ChemState *s, size_t cell)
{
    if (!s || cell >= s->cells) return CHEM_INVALID;
    s->explicit_air_density[cell] = 0;
    return CHEM_OK;
}

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

static int set_named(ChemState *s, const ChemNamedValues *items, size_t count, int parameters)
{
    if (!s || (count && !items)) return CHEM_INVALID;
    for (size_t i=0;i<count;++i) {
        int index=parameters ? chem_find_parameter(s->mechanism,items[i].name) :
                               chem_find_species(s->mechanism,items[i].name);
        if (index<0 || !items[i].values || items[i].count!=s->cells) return CHEM_INVALID;
    }
    for (size_t i=0;i<count;++i) {
        int index=parameters ? chem_find_parameter(s->mechanism,items[i].name) :
                               chem_find_species(s->mechanism,items[i].name);
        double *target=parameters ? s->parameters+(size_t)index*s->cells : s->values+(size_t)index*s->cells;
        memmove(target,items[i].values,s->cells*sizeof(*target));
    }
    return CHEM_OK;
}

int chem_set_concentrations(ChemState *s, const ChemNamedValues *items, size_t count)
{ return set_named(s,items,count,0); }

int chem_set_parameters(ChemState *s, const ChemNamedValues *items, size_t count)
{ return set_named(s,items,count,1); }

int chem_set_parameters_row_major(ChemState *s, const double *values, size_t cells, size_t parameters)
{
    if (!s || cells!=s->cells || parameters!=s->mechanism->parameter_count ||
        (cells && parameters && !values)) return CHEM_INVALID;
    for (size_t c=0;c<cells;++c) for (size_t i=0;i<parameters;++i)
        s->parameters[i*cells+c]=values[c*parameters+i];
    return CHEM_OK;
}

int chem_override_rate_function(ChemState *s, const char *label, ChemRateFunction function, void *context)
{
    if (!s || !label || !function) return CHEM_INVALID;
    for (size_t i = 0; i < s->mechanism->reaction_count; ++i) {
        const ChemRate *rate = &s->mechanism->reactions[i].rate;
        if (rate->kind == CHEM_CALLBACK && rate->label && !strcmp(rate->label, label)) {
            s->rate_functions[i] = function; s->rate_contexts[i] = context;
            return CHEM_OK;
        }
    }
    return CHEM_INVALID;
}

int chem_set_linear_solver(ChemState *s, const ChemLinearSolver *solver, void *context)
{
    if (!s || (solver && (!solver->factor || !solver->solve))) return CHEM_INVALID;
    s->linear_solver = solver ? *solver : (ChemLinearSolver){0};
    s->linear_solver_context = context;
    return CHEM_OK;
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
    case CHEM_MISSING_DIFFUSION: return "surface diffusion coefficient is missing";
    case CHEM_MISSING_MOLECULAR_WEIGHT: return "surface molecular weight is missing";
    case CHEM_UNCONVERGED: return "backward Euler accepted an unconverged integration";
    case CHEM_CONSTRAINT_INITIALIZATION_FAILED: return "Constraint Initialization Failed";
    default: return "unknown error";
    }
}
