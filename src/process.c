#include "internal.h"
#include <stdlib.h>
#include <string.h>

/* Each elementary operation visits contiguous cells. SIMD never reorders the
   accumulation of different reactions into the same species or matrix entry. */
void chem_forcing(const ChemState *s, const double *y, double *f)
{
    memset(f, 0, s->mechanism->data->species_count * s->cells * sizeof(double));
    chem_add_forcing(s, y, f);
}
void chem_add_forcing(const ChemState *s, const double *y, double *f)
{
    const ChemMechanismData *d = s->mechanism->data;
    const size_t nc = s->cells;
    double *restrict work = s->scratch;
    for (size_t r = 0; r < d->reaction_count; ++r) {
        const ChemReaction *rx = &d->reactions[r];
        memcpy(work, s->rates + r * nc, nc * sizeof(double));
        for (size_t j = 0; j < rx->reactant_count; ++j) {
            const double *restrict conc = y + (size_t)rx->reactants[j] * nc;
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) work[c] *= conc[c];
        }
        for (size_t j = 0; j < rx->reactant_count; ++j) {
            if (s->mechanism->algebraic[rx->reactants[j]]) continue;
            double *restrict dst = f + (size_t)rx->reactants[j] * nc;
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) dst[c] -= work[c];
        }
        for (size_t j = 0; j < rx->product_count; ++j) {
            if (s->mechanism->algebraic[rx->products[j].species]) continue;
            double *restrict dst = f + (size_t)rx->products[j].species * nc;
            double yield = rx->products[j].yield;
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) dst[c] += yield * work[c];
        }
    }
    for (size_t i = 0; i < d->external_model_count; ++i)
        if (d->external_models[i].forcing)
            d->external_models[i].forcing((ChemState *)s, y, f, d->external_models[i].context);
}
void chem_negative_jacobian(const ChemState *s, const double *y, double *jacobian)
{
    memset(jacobian, 0, s->mechanism->nnz * s->cells * sizeof(double));
    chem_subtract_jacobian(s, y, jacobian);
}
void chem_subtract_jacobian(const ChemState *s, const double *y, double *jacobian)
{
    const ChemMechanism *m = s->mechanism;
    const size_t nc = s->cells, n = m->data->species_count;
    double *restrict work = s->scratch;
    for (size_t job = 0; job < m->derivative_count; ++job) {
        const ChemDerivative *deriv = &m->derivatives[job];
        const ChemReaction *rx = &m->data->reactions[deriv->reaction];
        size_t col = deriv->independent;
        memcpy(work, s->rates + deriv->reaction * nc, nc * sizeof(double));
        for (size_t j = 0; j < rx->reactant_count; ++j) if (j != deriv->skipped) {
            const double *restrict conc = y + (size_t)rx->reactants[j] * nc;
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) work[c] *= conc[c];
        }
        for (size_t j = 0; j < rx->reactant_count; ++j) if (j != deriv->skipped) {
            if (m->algebraic[rx->reactants[j]]) continue;
            double *restrict dst = jacobian + m->index[(size_t)rx->reactants[j] * n + col] * nc;
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) dst[c] += work[c];
        }
        double *restrict diagonal = jacobian + m->diagonal[col] * nc;
        if (!m->algebraic[col]) {
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) diagonal[c] += work[c];
        }
        for (size_t j = 0; j < rx->product_count; ++j) {
            if (m->algebraic[rx->products[j].species]) continue;
            double *restrict dst = jacobian + m->index[(size_t)rx->products[j].species * n + col] * nc;
            double yield = rx->products[j].yield;
            #pragma omp simd
            for (size_t c = 0; c < nc; ++c) dst[c] -= yield * work[c];
        }
    }
    for (size_t i = 0; i < m->data->external_model_count; ++i) {
        const ChemExternalModel *model = &m->data->external_models[i];
        if (model->negative_jacobian) model->negative_jacobian((ChemState *)s, y, jacobian, model->context);
    }
}

int chem_reaction_species_used(const ChemMechanism *m, const char *name)
{
    if (!m || !name) return 0;
    const ChemMechanismData *d = m->data;
    for (size_t i = 0; i < d->reaction_count; ++i) {
        const ChemReaction *r = &d->reactions[i];
        for (size_t j = 0; j < r->reactant_count; ++j)
            if (!strcmp(name, d->species_names[r->reactants[j]])) return 1;
        for (size_t j = 0; j < r->product_count; ++j)
            if (!strcmp(name, d->species_names[r->products[j].species])) return 1;
        const ChemRate *rate = &m->rate_specs[i];
        if (rate->multiplier_names) for (size_t j = 0; j < rate->multiplier_count; ++j)
            if (rate->multiplier_names[j] && !strcmp(name, rate->multiplier_names[j])) return 1;
    }
    return 0;
}
int chem_kinetic_nonzero(const ChemMechanism *m, size_t row, size_t column)
{
    if (!m || row >= m->data->species_count || column >= m->data->species_count) return 0;
    for (size_t i = 0; i < m->data->reaction_count; ++i) {
        const ChemReaction *r = &m->data->reactions[i];
        int independent = 0;
        for (size_t j = 0; j < r->reactant_count; ++j) independent |= (size_t)r->reactants[j] == column;
        if (!independent) continue;
        for (size_t j = 0; j < r->reactant_count; ++j) if ((size_t)r->reactants[j] == row) return 1;
        for (size_t j = 0; j < r->product_count; ++j) if ((size_t)r->products[j].species == row) return 1;
    }
    return 0;
}
