#include "internal.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

ChemStatus chem_state_create(const ChemMechanism *m, size_t cells, ChemState **out)
{
    if (!out) return CHEM_INVALID;
    *out = NULL;
    if (!m || !cells) return CHEM_INVALID;
    const ChemMechanismData *d = m->data;
    size_t n = d->species_count;
    if (cells > SIZE_MAX / sizeof(double) / (n + m->nnz + d->reaction_count + d->parameter_count + 8)) return CHEM_INVALID;
    ChemState *s = calloc(1, sizeof(*s));
    if (!s) return CHEM_MEMORY;
    s->mechanism = m; s->cells = cells; s->rtol = 1e-6; s->linear = chem_linear_mozart();
    s->y = calloc(n * cells, sizeof(double));
    s->custom = calloc((d->parameter_count + 1) * cells, sizeof(double));
    s->temperature = malloc(cells * sizeof(double));
    s->pressure = malloc(cells * sizeof(double));
    s->air_density = calloc(cells, sizeof(double));
    s->rates = calloc((d->reaction_count + 1) * cells, sizeof(double));
    s->atol = malloc(n * sizeof(double));
    s->jacobian = calloc(m->nnz * cells, sizeof(double));
    s->forcing = calloc(n * cells, sizeof(double));
    s->candidate = calloc(n * cells, sizeof(double));
    s->error = calloc(n * cells, sizeof(double));
    s->stages = calloc(6 * n * cells, sizeof(double));
    s->scratch = calloc(cells + d->parameter_count + 2 * n + 1, sizeof(double));
    if (!s->y || !s->custom || !s->temperature || !s->pressure || !s->air_density || !s->rates || !s->atol || !s->jacobian
        || !s->forcing || !s->candidate || !s->error || !s->stages || !s->scratch) {
        chem_state_free(s); return CHEM_MEMORY;
    }
    for (size_t i = 0; i < n; ++i) {
        s->atol[i] = d->absolute_tolerances ? d->absolute_tolerances[i] : 1e-3;
        for (size_t c = 0; c < cells; ++c) s->y[i * cells + c] = d->initial_concentrations ? d->initial_concentrations[i] : 0;
    }
    for (size_t i = 0; i < d->parameter_count; ++i)
        for (size_t c = 0; c < cells; ++c) s->custom[i * cells + c] = d->initial_parameters ? d->initial_parameters[i] : 0;
    for (size_t c = 0; c < cells; ++c) { s->temperature[c] = d->temperature; s->pressure[c] = d->pressure; }
    *out = s;
    return CHEM_OK;
}
void chem_state_free(ChemState *s)
{
    if (!s) return;
    free(s->y); free(s->custom); free(s->temperature); free(s->pressure); free(s->air_density); free(s->rates);
    free(s->atol); free(s->jacobian); free(s->forcing); free(s->candidate); free(s->error); free(s->stages); free(s->scratch); free(s);
}
size_t chem_cell_count(const ChemState *s) { return s->cells; }
double *chem_concentrations(ChemState *s) { return s->y; }
double *chem_parameters(ChemState *s) { return s->custom; }
double *chem_temperatures(ChemState *s) { return s->temperature; }
double *chem_pressures(ChemState *s) { return s->pressure; }
double *chem_air_densities(ChemState *s) { return s->air_density; }
ChemStatus chem_calculate_air_density(ChemState *s)
{
    if (!s) return CHEM_INVALID;
    for (size_t c = 0; c < s->cells; ++c) {
        if (!(s->temperature[c] > 0) || !isfinite(s->temperature[c]) || !isfinite(s->pressure[c])) return CHEM_INVALID;
        s->air_density[c] = s->pressure[c] / ((1.380649e-23 * 6.02214076e23) * s->temperature[c]);
    }
    return CHEM_OK;
}
double *chem_absolute_tolerances(ChemState *s) { return s->atol; }
void chem_relative_tolerance(ChemState *s, double rtol) { s->rtol = rtol; }
const double *chem_rate_constants(const ChemState *s) { return s->rates; }
ChemStatus chem_update_rates(ChemState *s)
{
    const ChemMechanismData *d = s->mechanism->data;
    size_t cells = s->cells;
    for (size_t c = 0; c < cells; ++c)
        if (!(s->temperature[c] > 0) || !(s->pressure[c] >= 0)
            || !isfinite(s->temperature[c]) || !isfinite(s->pressure[c])) return CHEM_INVALID;
    for (size_t i = 0; i < d->external_model_count; ++i)
        if (d->external_models[i].update) d->external_models[i].update(s, d->external_models[i].context);
    for (size_t i = 0; i < d->constraint_count; ++i) {
        const ChemConstraint *c = &d->constraints[i];
        if (c->type != CHEM_EQUILIBRIUM_CONSTRAINT || c->parameter < 0) continue;
        for (size_t j = 0; j < cells; ++j)
            s->custom[(size_t)c->parameter * cells + j] = c->equilibrium_constant *
                exp((c->enthalpy / c->gas_constant) * (1 / s->temperature[j] - 1 / c->reference_temperature));
    }
    for (size_t i = 0; i < d->reaction_count; ++i) {
        ChemRate rate = s->mechanism->rate_specs[i];
        if ((rate.type == CHEM_USER || rate.type == CHEM_SURFACE) &&
            (rate.parameter < 0 || (size_t)rate.parameter + (rate.type == CHEM_SURFACE ? 2 : 1) > d->parameter_count)) return CHEM_INVALID;
        int base = rate.parameter;
        if (rate.type == CHEM_USER || rate.type == CHEM_SURFACE) rate.parameter = 0;
        for (size_t c = 0; c < cells; ++c) {
            double local[2] = {0, 0};
            const double *custom = local;
            if (rate.type == CHEM_USER || rate.type == CHEM_SURFACE) {
                local[0] = s->custom[(size_t)base * cells + c];
                if (rate.type == CHEM_SURFACE) local[1] = s->custom[((size_t)base + 1) * cells + c];
            } else if (rate.type == CHEM_CALLBACK) {
                for (size_t j = 0; j < d->parameter_count; ++j) s->scratch[j] = s->custom[j * cells + c];
                custom = s->scratch;
            }
            double air = s->air_density[c];
            double value = chem_rate_evaluate(&rate, s->temperature[c], s->pressure[c], air, custom);
            if (!isfinite(value)) return CHEM_NONFINITE;
            s->rates[i * cells + c] = value;
        }
    }
    return CHEM_OK;
}

void chem_post_solve_clamp(ChemState *s)
{
    size_t n = s->mechanism->data->species_count, nc = s->cells;
    for (size_t i = 0; i < n; ++i) if (!s->mechanism->algebraic[i]) {
        double *restrict values = s->y + i * nc;
        #pragma omp simd
        for (size_t c = 0; c < nc; ++c) values[c] = values[c] > 0 ? values[c] : 0;
    }
}
