#include "chem_builder.h"
#include "internal.h"
#include <stdlib.h>
#include <string.h>

/* Tables are owned here; strings, equations, reaction arrays, and model contexts
   are borrowed. Keeping this owner with the mechanism permits builder reuse. */
typedef struct {
    ChemMechanismData data;
    const char **species_names, **parameter_names;
    double *initial_variables, *initial_parameters, *absolute_tolerances;
    ChemExternalModel *models;
    ChemConstraint *constraints;
} BuiltData;
static void free_built_data(void *pointer)
{
    BuiltData *b = pointer;
    if (!b) return;
    free(b->species_names); free(b->parameter_names); free(b->initial_variables);
    free(b->initial_parameters); free(b->absolute_tolerances); free(b->models); free(b->constraints); free(b);
}
static size_t find_name(size_t count, const char *const *names, const char *name)
{
    if (name && names) for (size_t i = 0; i < count; ++i) if (names[i] && !strcmp(names[i], name)) return i;
    return SIZE_MAX;
}
size_t chem_data_species_index(const ChemMechanismData *d, const char *name)
{ return d ? find_name(d->species_count, d->species_names, name) : SIZE_MAX; }
size_t chem_data_parameter_index(const ChemMechanismData *d, const char *name)
{ return d ? find_name(d->parameter_count, d->parameter_names, name) : SIZE_MAX; }
static int accumulate(size_t *count, size_t addition)
{
    if (addition > SIZE_MAX - *count) return 0;
    *count += addition;
    return 1;
}
static int all_used(const ChemBuilder *builder, const ChemMechanism *m)
{
    const ChemMechanismData *d = m->data;
    for (size_t i = 0; i < d->species_count; ++i) {
        int used = chem_reaction_species_used(m, d->species_names[i]);
        for (size_t j = 0; !used && j < d->constraint_count; ++j) {
            const ChemConstraint *c = &d->constraints[j];
            used |= (size_t)c->algebraic_species == i;
            for (size_t k = 0; k < c->reactant_count; ++k) used |= (size_t)c->reactants[k].species == i;
            for (size_t k = 0; k < c->product_count; ++k) used |= (size_t)c->products[k].species == i;
        }
        for (size_t j = 0; !used && j < builder->model_count; ++j) {
            const ChemModelDefinition *model = &builder->models[j];
            used |= find_name(model->variable_count, model->variable_names, d->species_names[i]) != SIZE_MAX;
            used |= find_name(model->used_species_count, model->used_species_names, d->species_names[i]) != SIZE_MAX;
        }
        if (!used) return 0;
    }
    return 1;
}
ChemStatus chem_builder_build(const ChemBuilder *builder, ChemMechanism **out)
{
    if (!out) return CHEM_INVALID;
    *out = NULL;
    if (!builder || !builder->base.species_count || !builder->base.species_names
        || (builder->model_count && !builder->models)) return CHEM_INVALID;
    const ChemMechanismData *base = &builder->base;
    size_t n = base->species_count, np = base->parameter_count, nc = base->constraint_count, nm = base->external_model_count;
    if ((np && !base->parameter_names) || (nc && !base->constraints) || (nm && !base->external_models)) return CHEM_INVALID;
    for (size_t i = 0; i < builder->model_count; ++i) {
        const ChemModelDefinition *model = &builder->models[i];
        if ((model->variable_count && !model->variable_names) || (model->parameter_count && !model->parameter_names)
            || (model->constraint_count && !model->configure) || (model->used_species_count && !model->used_species_names)
            || !accumulate(&n, model->variable_count) || !accumulate(&np, model->parameter_count)
            || !accumulate(&nc, model->constraint_count) || !accumulate(&nm, 1)) return CHEM_INVALID;
    }
    if (n > SIZE_MAX / sizeof(double) || n > SIZE_MAX / sizeof(char *) || np > SIZE_MAX / sizeof(double)
        || np > SIZE_MAX / sizeof(char *) || nc > SIZE_MAX / sizeof(ChemConstraint)
        || nm > SIZE_MAX / sizeof(ChemExternalModel)) return CHEM_INVALID;
    BuiltData *b = calloc(1, sizeof(*b));
    if (!b) return CHEM_MEMORY;
    b->data = *base;
    b->species_names = calloc(n, sizeof(*b->species_names));
    b->parameter_names = calloc(np ? np : 1, sizeof(*b->parameter_names));
    b->initial_variables = calloc(n, sizeof(double));
    b->absolute_tolerances = malloc(n * sizeof(double));
    b->initial_parameters = calloc(np ? np : 1, sizeof(double));
    b->constraints = calloc(nc ? nc : 1, sizeof(*b->constraints));
    b->models = calloc(nm ? nm : 1, sizeof(*b->models));
    if (!b->species_names || !b->parameter_names || !b->initial_variables || !b->absolute_tolerances
        || !b->initial_parameters || !b->constraints || !b->models) { free_built_data(b); return CHEM_MEMORY; }
    for (size_t i = 0; i < n; ++i) b->absolute_tolerances[i] = 1e-3;
    for (size_t i = 0; i < base->species_count; ++i) {
        b->species_names[i] = base->species_names[i];
        b->initial_variables[i] = base->initial_concentrations ? base->initial_concentrations[i] : 0;
        if (base->absolute_tolerances) b->absolute_tolerances[i] = base->absolute_tolerances[i];
    }
    for (size_t i = 0; i < base->parameter_count; ++i) {
        b->parameter_names[i] = base->parameter_names[i];
        b->initial_parameters[i] = base->initial_parameters ? base->initial_parameters[i] : 0;
    }
    if (base->constraint_count) memcpy(b->constraints, base->constraints, base->constraint_count * sizeof(*b->constraints));
    if (base->external_model_count) memcpy(b->models, base->external_models, base->external_model_count * sizeof(*b->models));
    size_t variable = base->species_count, parameter = base->parameter_count;
    for (size_t i = 0; i < builder->model_count; ++i) {
        const ChemModelDefinition *model = &builder->models[i];
        for (size_t j = 0; j < model->variable_count; ++j, ++variable) {
            b->species_names[variable] = model->variable_names[j];
            b->initial_variables[variable] = model->initial_variables ? model->initial_variables[j] : 0;
            if (model->absolute_tolerances) b->absolute_tolerances[variable] = model->absolute_tolerances[j];
        }
        for (size_t j = 0; j < model->parameter_count; ++j, ++parameter) {
            b->parameter_names[parameter] = model->parameter_names[j];
            b->initial_parameters[parameter] = model->initial_parameters ? model->initial_parameters[j] : 0;
        }
    }
    b->data.species_count = n; b->data.parameter_count = np; b->data.constraint_count = nc; b->data.external_model_count = nm;
    b->data.species_names = b->species_names; b->data.parameter_names = b->parameter_names;
    b->data.initial_concentrations = b->initial_variables; b->data.absolute_tolerances = b->absolute_tolerances;
    b->data.initial_parameters = b->initial_parameters; b->data.constraints = b->constraints; b->data.external_models = b->models;
    size_t constraint = base->constraint_count;
    for (size_t i = 0; i < builder->model_count; ++i) {
        const ChemModelDefinition *model = &builder->models[i];
        ChemExternalModel *process = &b->models[base->external_model_count + i];
        process->context = model->context;
        if (model->configure) {
            ChemStatus status = model->configure(&b->data, process, b->constraints + constraint, model->context);
            if (status != CHEM_OK) { free_built_data(b); return status; }
        }
        constraint += model->constraint_count;
    }
    ChemStatus status = chem_mechanism_create(&b->data, out);
    if (status != CHEM_OK) { free_built_data(b); return status; }
    if (!builder->ignore_unused_species && !all_used(builder, *out)) {
        chem_mechanism_free(*out); *out = NULL; free_built_data(b); return CHEM_INVALID;
    }
    (*out)->owned_descriptor = b;
    (*out)->free_descriptor = free_built_data;
    return CHEM_OK;
}
