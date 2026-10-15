#include "internal.h"
#include <stdlib.h>
#include <string.h>

/* Each elementary operation visits contiguous cells. SIMD never reorders the
   accumulation of different reactions into the same species or matrix entry. */
void chem_forcing(const ChemState *s, const double *y, double *f)
{
    const ChemMechanismData *d = s->mechanism->data;
    const size_t nc = s->cells;
    double *restrict work = s->scratch;
    memset(f, 0, d->species_count * nc * sizeof(double));
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
    const ChemMechanism *m = s->mechanism;
    const size_t nc = s->cells, n = m->data->species_count;
    double *restrict work = s->scratch;
    memset(jacobian, 0, m->nnz * nc * sizeof(double));
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
