/* Constraint residuals and derivatives ported from the reference.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <math.h>
#include <string.h>

ChemConstraint chem_constraint_default(ChemConstraintType type)
{
    return (ChemConstraint){.type=type, .algebraic_species=-1, .equilibrium_constant=1,
        .gas_constant=1.380649e-23*6.02214076e23, .reference_temperature=298.15, .parameter=-1};
}
static double term_product(const ChemProduct *terms, size_t count, const double *y, size_t nc, size_t cell,
                           double initial, size_t skip)
{
    double result = initial;
    for (size_t j = 0; j < count; ++j) if (j != skip) {
        double conc = y[(size_t)terms[j].species * nc + cell];
        result *= pow(conc > 0 ? conc : 0, terms[j].yield);
    }
    return result;
}
static double value(const ChemState *s, const ChemConstraint *c, const double *y, size_t cell, double *gradient)
{
    size_t nc = s->cells, n = s->mechanism->data->species_count;
    if (gradient) memset(gradient, 0, n * sizeof(double));
    if (c->type == CHEM_CALLBACK_CONSTRAINT) {
        double *local_y = s->scratch + nc, *local_custom = local_y + n;
        for (size_t j = 0; j < n; ++j) local_y[j] = y[j * nc + cell];
        for (size_t j = 0; j < s->mechanism->data->parameter_count; ++j) local_custom[j] = s->custom[j * nc + cell];
        return c->function(s->temperature[cell], local_y, local_custom, gradient, c->context);
    }
    if (c->type == CHEM_LINEAR_CONSTRAINT) {
        double residual = -c->constant;
        for (size_t j = 0; j < c->reactant_count; ++j) {
            ChemProduct term = c->reactants[j];
            residual += term.yield * y[(size_t)term.species * nc + cell];
            if (gradient) gradient[term.species] += term.yield;
        }
        return residual;
    }
    double keq = c->parameter >= 0 ? s->custom[(size_t)c->parameter * nc + cell] : c->equilibrium_constant *
        exp((c->enthalpy / c->gas_constant) * (1 / s->temperature[cell] - 1 / c->reference_temperature));
    double rp = term_product(c->reactants, c->reactant_count, y, nc, cell, keq, SIZE_MAX);
    double pp = term_product(c->products, c->product_count, y, nc, cell, 1, SIZE_MAX);
    if (gradient) {
        for (int side = 0; side < 2; ++side) {
            size_t count = side ? c->product_count : c->reactant_count;
            const ChemProduct *terms = side ? c->products : c->reactants;
            for (size_t i = 0; i < count; ++i) {
                double product = term_product(terms, count, y, nc, cell, side ? 1 : keq, i);
                double conc = y[(size_t)terms[i].species * nc + cell], stoich = terms[i].yield;
                double partial = stoich == 1 ? product : conc > 0 ? stoich * product * pow(conc, stoich - 1) : 0;
                gradient[terms[i].species] += side ? -partial : partial;
            }
        }
    }
    return rp - pp;
}
void chem_constraint_forcing(const ChemState *s, const double *y, double *forcing)
{
    const ChemMechanismData *d = s->mechanism->data;
    for (size_t i = 0; i < d->constraint_count; ++i) {
        const ChemConstraint *c = &d->constraints[i];
        for (size_t cell = 0; cell < s->cells; ++cell)
            forcing[(size_t)c->algebraic_species * s->cells + cell] = value(s, c, y, cell, NULL);
    }
}
void chem_constraint_negative_jacobian(const ChemState *s, const double *y, double *jacobian)
{
    const ChemMechanism *m = s->mechanism;
    const ChemMechanismData *d = m->data;
    size_t n = d->species_count;
    double *gradient = s->scratch + s->cells + n + d->parameter_count;
    for (size_t i = 0; i < d->constraint_count; ++i) {
        const ChemConstraint *c = &d->constraints[i];
        for (size_t cell = 0; cell < s->cells; ++cell) {
            value(s, c, y, cell, gradient);
            for (size_t col = 0; col < n; ++col) {
                size_t index = m->index[(size_t)c->algebraic_species * n + col];
                if (index != SIZE_MAX) jacobian[index * s->cells + cell] -= gradient[col];
            }
        }
    }
}

static int finite_algebraic(const ChemState *s, const double *v, double *maximum)
{
    *maximum = 0;
    for (size_t i = 0; i < s->mechanism->data->species_count; ++i) if (s->mechanism->algebraic[i])
        for (size_t c = 0; c < s->cells; ++c) {
            double a = fabs(v[i * s->cells + c]);
            if (!isfinite(a)) return 0;
            if (a > *maximum) *maximum = a;
        }
    return 1;
}
static double correction_norm(const ChemState *s, const double *y, const double *delta)
{
    double maximum = 0;
    for (size_t i = 0; i < s->mechanism->data->species_count; ++i) if (s->mechanism->algebraic[i])
        for (size_t c = 0; c < s->cells; ++c) {
            size_t j = i * s->cells + c;
            double zmax = fmax(fabs(y[j]), fabs(y[j] + delta[j]));
            double scale = s->atol[i] + s->rtol * zmax;
            double q = scale > 0 ? fabs(delta[j]) / scale : delta[j] == 0 ? 0 : INFINITY;
            if (q > maximum) maximum = q;
        }
    return maximum;
}
static void apply_correction(const ChemState *s, double *y, const double *delta, double step)
{
    for (size_t i = 0; i < s->mechanism->data->species_count; ++i) if (s->mechanism->algebraic[i])
        for (size_t c = 0; c < s->cells; ++c) {
            size_t j = i * s->cells + c;
            y[j] += step * delta[j];
        }
}
ChemStatus chem_initialize_constraints(ChemState *s, const ChemRosenbrock *p, ChemStats *stats)
{
    if (!s || !p || !stats) return CHEM_INVALID;
    const ChemMechanism *m = s->mechanism;
    if (!m->data->constraint_count) return CHEM_OK;
    size_t nc = s->cells, n = m->data->species_count, bytes = n * nc * sizeof(double);
    double *delta = s->forcing, *original = s->stages;
    int modified = 0;
    ChemStatus status = CHEM_INITIALIZATION_FAILED;
    for (size_t update = 0; update <= p->constraint_init_max_iterations; ++update) {
        memset(delta, 0, bytes);
        chem_constraint_forcing(s, s->y, delta);
        double maximum;
        if (!finite_algebraic(s, delta, &maximum)) { status = CHEM_NONFINITE; goto rollback; }
        ++stats->constraint_init_iterations;
        if (maximum == 0) return CHEM_OK;
        memset(s->jacobian, 0, m->nnz * nc * sizeof(double));
        chem_constraint_negative_jacobian(s, s->y, s->jacobian);
        ++stats->jacobian_updates;
        for (size_t i = 0; i < n; ++i) if (!m->algebraic[i])
            for (size_t c = 0; c < nc; ++c) s->jacobian[m->diagonal[i] * nc + c] += 1;
        status = s->linear.factor(m, nc, s->jacobian);
        ++stats->decompositions;
        if (status != CHEM_OK) { status = CHEM_INITIALIZATION_FAILED; goto rollback; }
        s->linear.solve(m, nc, s->jacobian, delta);
        ++stats->solves;
        if (!finite_algebraic(s, delta, &maximum)) { status = CHEM_NONFINITE; goto rollback; }
        double norm = correction_norm(s, s->y, delta);
        if (norm <= p->constraint_init_tolerance) return CHEM_OK;
        status = CHEM_INITIALIZATION_FAILED;
        if (update == p->constraint_init_max_iterations) break;
        if (!modified) { memcpy(original, s->y, bytes); modified = 1; }
        if (p->constraint_init_max_backtracks > 0) {
            double step = 1;
            size_t backtrack = 0;
            for (;; ++backtrack, step *= p->constraint_init_backtrack_factor) {
                if (backtrack > p->constraint_init_max_backtracks) goto rollback;
                memcpy(s->candidate, s->y, bytes);
                apply_correction(s, s->candidate, delta, step);
                if (!finite_algebraic(s, s->candidate, &maximum)) continue;
                memset(s->error, 0, bytes);
                chem_constraint_forcing(s, s->candidate, s->error);
                if (!finite_algebraic(s, s->error, &maximum)) continue;
                if (maximum == 0) break;
                s->linear.solve(m, nc, s->jacobian, s->error);
                ++stats->solves;
                if (!finite_algebraic(s, s->error, &maximum)) continue;
                double trial_norm = correction_norm(s, s->candidate, s->error);
                if (trial_norm <= p->constraint_init_tolerance || trial_norm <
                    (1 - p->constraint_init_sufficient_decrease * step) * norm) break;
            }
            memcpy(s->y, s->candidate, bytes);
        } else apply_correction(s, s->y, delta, 1);
    }
rollback:
    if (modified) memcpy(s->y, original, bytes);
    return status;
}
