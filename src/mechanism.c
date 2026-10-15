#include "internal.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int valid_data(const ChemMechanismData *d)
{
    if (!d || !d->species_count || !d->species_names || (d->reaction_count && !d->reactions)
        || (d->parameter_count && !d->parameter_names) || d->species_count > SIZE_MAX / d->species_count)
        return 0;
    for (size_t i = 0; i < d->species_count; ++i) {
        if (!d->species_names[i]) return 0;
        for (size_t j = 0; j < i; ++j)
            if (!strcmp(d->species_names[i], d->species_names[j])) return 0;
    }
    for (size_t i = 0; i < d->reaction_count; ++i) {
        const ChemReaction *r = &d->reactions[i];
        if ((r->reactant_count && !r->reactants) || (r->product_count && !r->products)) return 0;
        for (size_t j = 0; j < r->reactant_count; ++j)
            if (r->reactants[j] < 0 || (size_t)r->reactants[j] >= d->species_count) return 0;
        for (size_t j = 0; j < r->product_count; ++j)
            if (r->products[j].species < 0 || (size_t)r->products[j].species >= d->species_count
                || !isfinite(r->products[j].yield)) return 0;
        const ChemRate *rate = &r->rate;
        if (rate->type < CHEM_ARRHENIUS || rate->type > CHEM_CALLBACK || rate->third_body_count < 0
            || rate->coefficient_count < 0 || rate->coefficient_count > 19) return 0;
        if (rate->type == CHEM_USER || rate->type == CHEM_SURFACE) {
            size_t needed = rate->type == CHEM_SURFACE ? 2 : 1;
            if (rate->parameter < 0 || (size_t)rate->parameter + needed > d->parameter_count) return 0;
        }
        if (rate->type == CHEM_SURFACE && (!(rate->p[0] > 0) || !(rate->p[1] > 0))) return 0;
        if (rate->multiplier_count && !rate->multipliers) return 0;
    }
    if ((d->constraint_count && !d->constraints) || (d->external_model_count && !d->external_models)) return 0;
    for (size_t i = 0; i < d->constraint_count; ++i) {
        const ChemConstraint *c = &d->constraints[i];
        if (c->type < CHEM_LINEAR_CONSTRAINT || c->type > CHEM_CALLBACK_CONSTRAINT
            || c->algebraic_species < 0 || (size_t)c->algebraic_species >= d->species_count
            || !c->reactant_count || !c->reactants || (c->product_count && !c->products)) return 0;
        for (size_t j = 0; j < i; ++j) if (d->constraints[j].algebraic_species == c->algebraic_species) return 0;
        for (size_t j = 0; j < c->reactant_count + c->product_count; ++j) {
            ChemProduct term = j < c->reactant_count ? c->reactants[j] : c->products[j - c->reactant_count];
            if (term.species < 0 || (size_t)term.species >= d->species_count || !isfinite(term.yield)) return 0;
            if (c->type == CHEM_EQUILIBRIUM_CONSTRAINT && !(term.yield > 0)) return 0;
        }
        if (c->type == CHEM_EQUILIBRIUM_CONSTRAINT && (!c->product_count || !(c->equilibrium_constant > 0)
            || !(c->gas_constant > 0) || !(c->reference_temperature > 0)
            || (c->parameter >= 0 && (size_t)c->parameter >= d->parameter_count))) return 0;
        if (c->type == CHEM_CALLBACK_CONSTRAINT && !c->function) return 0;
    }
    for (size_t i = 0; i < d->external_model_count; ++i) {
        const ChemExternalModel *model = &d->external_models[i];
        if (model->nonzero_count && !model->nonzeros) return 0;
        for (size_t j = 0; j < model->nonzero_count; ++j)
            if (model->nonzeros[j].row >= d->species_count || model->nonzeros[j].column >= d->species_count) return 0;
    }
    return 1;
}

ChemStatus chem_mechanism_create(const ChemMechanismData *d, ChemMechanism **out)
{
    if (!out) return CHEM_INVALID;
    *out = NULL;
    if (!valid_data(d)) return CHEM_INVALID;
    const size_t n = d->species_count;
    ChemMechanism *m = calloc(1, sizeof(*m));
    unsigned char *pattern = calloc(n * n, 1);
    if (!m || !pattern) { free(m); free(pattern); return CHEM_MEMORY; }
    m->data = d;
    m->rate_specs = calloc(d->reaction_count + 1, sizeof(*m->rate_specs));
    if (!m->rate_specs) goto memory_error;
    for (size_t i = 0; i < d->reaction_count; ++i) m->rate_specs[i] = d->reactions[i].rate;
    for (size_t i = 0; i < n; ++i) pattern[i * n + i] = 1;
    for (size_t r = 0; r < d->reaction_count; ++r) {
        const ChemReaction *rx = &d->reactions[r];
        m->derivative_count += rx->reactant_count;
        for (size_t j = 0; j < rx->reactant_count; ++j) {
            size_t col = (size_t)rx->reactants[j];
            for (size_t k = 0; k < rx->reactant_count; ++k) pattern[rx->reactants[k] * n + col] = 1;
            for (size_t k = 0; k < rx->product_count; ++k) pattern[rx->products[k].species * n + col] = 1;
        }
    }
    m->algebraic = calloc(n, 1);
    m->constraint_pattern = calloc(n * n, 1);
    if (!m->algebraic || !m->constraint_pattern) goto memory_error;
    for (size_t i = 0; i < d->constraint_count; ++i) {
        const ChemConstraint *c = &d->constraints[i];
        size_t row = (size_t)c->algebraic_species;
        m->algebraic[row] = 1;
        m->constraint_pattern[row * n + row] = 1;
        for (size_t j = 0; j < c->reactant_count; ++j) pattern[row * n + c->reactants[j].species] = m->constraint_pattern[row * n + c->reactants[j].species] = 1;
        for (size_t j = 0; j < c->product_count; ++j) pattern[row * n + c->products[j].species] = m->constraint_pattern[row * n + c->products[j].species] = 1;
    }
    for (size_t i = 0; i < d->external_model_count; ++i)
        for (size_t j = 0; j < d->external_models[i].nonzero_count; ++j) {
            ChemMatrixIndex ij = d->external_models[i].nonzeros[j];
            pattern[ij.row * n + ij.column] = 1;
        }
    for (size_t i = 0; i < n; ++i)
        for (size_t k = i + 1; k < n; ++k) if (pattern[i * n + k])
            for (size_t j = i + 1; j < n; ++j) if (pattern[j * n + i]) pattern[j * n + k] = 1;
    for (size_t i = 0; i < n * n; ++i) m->nnz += pattern[i] != 0;
    m->index = malloc(n * n * sizeof(*m->index));
    m->row_start = calloc(n + 1, sizeof(*m->row_start));
    m->entries = malloc(m->nnz * sizeof(*m->entries));
    m->diagonal = malloc(n * sizeof(*m->diagonal));
    m->lower_start = calloc(n + 1, sizeof(*m->lower_start));
    m->lower = malloc(m->nnz * sizeof(*m->lower));
    m->operation_start = calloc(n + 1, sizeof(*m->operation_start));
    m->derivatives = calloc(m->derivative_count + 1, sizeof(*m->derivatives));
    if (!m->index || !m->row_start || !m->entries || !m->diagonal || !m->lower_start
        || !m->lower || !m->operation_start || !m->derivatives) goto memory_error;
    size_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        m->row_start[i] = count;
        for (size_t j = 0; j < n; ++j) {
            m->index[i * n + j] = pattern[i * n + j] ? count : SIZE_MAX;
            if (pattern[i * n + j]) { m->entries[count] = (ChemEntry){count, j}; ++count; }
        }
        m->diagonal[i] = m->index[i * n + i];
    }
    m->row_start[n] = count;
    size_t lower_count = 0, operations = 0;
    for (size_t i = 0; i < n; ++i) {
        m->lower_start[i] = lower_count;
        m->operation_start[i] = operations;
        for (size_t j = i + 1; j < n; ++j)
            if (pattern[j * n + i]) m->lower[lower_count++] = m->index[j * n + i];
        for (size_t k = i + 1; k < n; ++k) if (pattern[i * n + k])
            for (size_t j = i + 1; j < n; ++j) if (pattern[j * n + i]) ++operations;
    }
    m->lower_start[n] = lower_count;
    m->operation_start[n] = operations;
    m->operations = malloc((operations + 1) * sizeof(*m->operations));
    if (!m->operations) goto memory_error;
    count = 0;
    for (size_t i = 0; i < n; ++i)
        for (size_t k = i + 1; k < n; ++k) if (pattern[i * n + k])
            for (size_t j = i + 1; j < n; ++j) if (pattern[j * n + i])
                m->operations[count++] = (ChemElimination){m->index[j * n + k], m->index[j * n + i], m->index[i * n + k]};
    count = 0;
    /* Match the reference's stable independent-variable, reaction ordering.
       Repeated reactants produce repeated derivative jobs, each skipping the
       first matching occurrence, and hence the correct multiplicity. */
    for (size_t col = 0; col < n; ++col)
        for (size_t r = 0; r < d->reaction_count; ++r) {
            const ChemReaction *rx = &d->reactions[r];
            size_t first = SIZE_MAX;
            for (size_t j = 0; j < rx->reactant_count; ++j)
                if ((size_t)rx->reactants[j] == col) {
                    if (first == SIZE_MAX) first = j;
                    m->derivatives[count++] = (ChemDerivative){r, col, first};
                }
        }
    free(pattern);
    *out = m;
    return CHEM_OK;
memory_error:
    free(pattern);
    chem_mechanism_free(m);
    return CHEM_MEMORY;
}

void chem_mechanism_free(ChemMechanism *m)
{
    if (!m) return;
    free(m->index); free(m->row_start); free(m->entries); free(m->diagonal);
    free(m->lower_start); free(m->lower); free(m->operation_start);
    free(m->rate_specs); free(m->operations); free(m->derivatives); free(m->algebraic); free(m->constraint_pattern); free(m);
}
const ChemMechanismData *chem_mechanism_data(const ChemMechanism *m) { return m->data; }
size_t chem_jacobian_size(const ChemMechanism *m) { return m->nnz; }
size_t chem_jacobian_index(const ChemMechanism *m, size_t r, size_t c)
{
    size_t n = m->data->species_count;
    return r < n && c < n ? m->index[r * n + c] : SIZE_MAX;
}
const char *chem_status_string(ChemStatus s)
{
    static const char *const names[] = {"converged", "invalid argument", "allocation failed", "maximum steps exceeded",
        "step size too small", "nonfinite value", "singular matrix", "accepting unconverged integration", "Constraint Initialization Failed",
        "property not found", "invalid property type"};
    return s >= CHEM_OK && s <= CHEM_WRONG_TYPE ? names[s] : "unknown status";
}

int chem_is_algebraic(const ChemMechanism *m, size_t species)
{ return species < m->data->species_count ? m->algebraic[species] : 0; }
int chem_constraint_nonzero(const ChemMechanism *m, size_t row, size_t col)
{
    size_t n = m->data->species_count;
    return row < n && col < n ? m->constraint_pattern[row * n + col] : 0;
}

ChemRate *chem_rate_by_name(ChemMechanism *m, const char *name)
{
    if (!m || !name) return NULL;
    for (size_t i = 0; i < m->data->reaction_count; ++i)
        if (m->rate_specs[i].name && !strcmp(m->rate_specs[i].name, name)) return &m->rate_specs[i];
    return NULL;
}

static ChemRateFamily rate_family(ChemRateType t)
{
    static const ChemRateFamily families[] = {
        CHEM_RATE_ARRHENIUS, CHEM_RATE_TROE, CHEM_RATE_TERNARY, CHEM_RATE_JPL19,
        CHEM_RATE_TUNNELING, CHEM_RATE_BRANCHED, CHEM_RATE_BRANCHED,
        CHEM_RATE_TAYLOR, CHEM_RATE_REVERSIBLE, CHEM_RATE_USER, CHEM_RATE_SURFACE, CHEM_RATE_CALLBACK
    };
    return families[t];
}
size_t chem_rate_family_count(const ChemMechanism *m, ChemRateFamily family)
{
    size_t count = 0;
    if (m) for (size_t i = 0; i < m->data->reaction_count; ++i) count += rate_family(m->rate_specs[i].type) == family;
    return count;
}
size_t chem_rate_family_offset(const ChemMechanism *m, ChemRateFamily family)
{
    size_t count = 0;
    if (m) for (size_t i = 0; i < m->data->reaction_count; ++i) count += rate_family(m->rate_specs[i].type) < family;
    return count;
}
size_t chem_rate_family_index(const ChemMechanism *m, ChemRateFamily family, size_t ordinal)
{
    if (m) for (size_t i = 0; i < m->data->reaction_count; ++i)
        if (rate_family(m->rate_specs[i].type) == family && ordinal-- == 0) return i;
    return SIZE_MAX;
}
const ChemRate *chem_rate_specification(const ChemMechanism *m, size_t i)
{ return m && i < m->data->reaction_count ? &m->rate_specs[i] : NULL; }

const char *chem_version(void) { return "3.14.0"; }
unsigned chem_version_major(void) { return 3; }
unsigned chem_version_minor(void) { return 14; }
unsigned chem_version_patch(void) { return 0; }
unsigned chem_version_tweak(void) { return 0; }
