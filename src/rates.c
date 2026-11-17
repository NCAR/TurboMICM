#include "internal.h"

#include <math.h>
#include <string.h>

static const double gas_constant = 1.380649e-23 * 6.02214076e23;
static const double pi = 3.14159265358979323846;

static double falloff(const ChemRate *p, double temperature, double density, double scale)
{
    double k0 = p->k0_a * exp(p->k0_c / temperature) * pow(temperature / 300.0, p->k0_b);
    double kinf = p->kinf_a * exp(p->kinf_c / temperature) * pow(temperature / 300.0, p->kinf_b);
    double ratio = k0 * density / kinf;
    return k0 * scale / (1.0 + ratio) * pow(p->fc, p->n / (p->n + pow(log10(ratio), 2.0)));
}

static double rate_value(const ChemRate *p, double temperature, double pressure,
                         double density, const double *parameters, size_t cells, size_t cell)
{
    double k0, kinf, kint, ratio, f, a, b, z, ref, poly, power;
    switch (p->kind) {
    case CHEM_ARRHENIUS:
        return p->a * exp(p->c / temperature) * pow(temperature / p->d, p->b) * (1.0 + p->e * pressure);
    case CHEM_TROE: return falloff(p, temperature, density, density);
    case CHEM_TERNARY: return falloff(p, temperature, density, 1.0);
    case CHEM_TERNARY_JPL19:
        k0 = p->k0_a * exp(p->k0_c / temperature) * pow(temperature / p->k0_d, p->k0_b);
        kinf = p->kinf_a * exp(p->kinf_c / temperature) * pow(temperature / p->kinf_d, p->kinf_b);
        kint = p->kint_a * exp(p->kint_c / temperature) * pow(temperature / p->kint_d, p->kint_b);
        ratio = k0 * density / kinf;
        f = k0 * density / (1.0 + ratio) * pow(p->fc, p->n / (p->n + pow(log10(ratio), 2.0)));
        return f + kint * (1.0 - f / kinf);
    case CHEM_BRANCHED:
        k0 = 2.0e-22 * 6.02214076e23 * 1.0e-6 * exp(p->branch_n);
        ref = 2.45e19 / 6.02214076e23 * 1.0e6;
        a = k0 * ref;
        b = 0.43 * pow(293.0 / 298.0, -8.0);
        z = a / (1.0 + a/b) * pow(0.41, 1.0 / (1.0 + pow(log10(a/b), 2.0))) * (1.0-p->branch_a0) / p->branch_a0;
        a = k0 * density;
        b = 0.43 * pow(temperature / 298.0, -8.0);
        f = a / (1.0 + a/b) * pow(0.41, 1.0 / (1.0 + pow(log10(a/b), 2.0)));
        return p->branch_x * exp(-p->branch_y / temperature) *
            (p->nitrate_branch ? f / (f+z) : z / (z+f));
    case CHEM_TUNNELING:
        return p->a * exp(-p->b / temperature + p->c / (temperature*temperature*temperature));
    case CHEM_TAYLOR:
        poly = 0; power = 1;
        for (size_t i = 0; i < p->coefficient_count; ++i, power *= temperature) poly += p->coefficients[i] * power;
        return poly * p->a * exp(p->c / temperature) * pow(temperature / p->d, p->b) * (1.0+p->e*pressure);
    case CHEM_REVERSIBLE: return p->a * exp(p->c / temperature) * p->b;
    case CHEM_USER: return parameters[(size_t)p->parameter*cells+cell] * p->a;
    case CHEM_SURFACE:
        a = parameters[(size_t)p->parameter*cells+cell];
        b = parameters[((size_t)p->parameter+1)*cells+cell];
        f = sqrt(8.0 * gas_constant / (pi * p->molecular_weight) * temperature);
        return 4.0 * b * pi * a * a / (a / p->diffusion + 4.0 / (f * p->probability));
    case CHEM_CALLBACK: return p->function(temperature, pressure, density, p->context);
    }
    return NAN;
}

ChemRate chem_rate_defaults(ChemRateKind kind)
{
    ChemRate p = {0};
    p.kind = kind;
    p.parameter = -1;
    switch (kind) {
    case CHEM_ARRHENIUS: case CHEM_TAYLOR:
        p.a = 1; p.d = 300;
        if (kind == CHEM_TAYLOR) { p.coefficient_count = 1; p.coefficients[0] = 1; }
        break;
    case CHEM_TROE: case CHEM_TERNARY: case CHEM_TERNARY_JPL19:
        p.k0_a = p.kinf_a = p.kint_a = 1; p.fc = .6; p.n = 1;
        p.k0_d = p.kinf_d = p.kint_d = 298;
        break;
    case CHEM_TUNNELING: case CHEM_REVERSIBLE: case CHEM_USER: p.a = 1; break;
    case CHEM_SURFACE: p.probability = 1; break;
    case CHEM_BRANCHED: case CHEM_CALLBACK: break;
    }
    return p;
}

static int validate_rate(const ChemRate *p, const double *parameters, size_t parameter_count)
{
    if (!p || p->kind < CHEM_ARRHENIUS || p->kind > CHEM_CALLBACK) return CHEM_INVALID;
    if ((p->kind == CHEM_USER || p->kind == CHEM_SURFACE) &&
        (!parameters || p->parameter < 0 || (size_t)p->parameter >= parameter_count ||
         (p->kind == CHEM_SURFACE && (size_t)p->parameter+1 >= parameter_count))) return CHEM_INVALID;
    if (p->kind == CHEM_SURFACE && p->diffusion <= 0) return CHEM_MISSING_DIFFUSION;
    if (p->kind == CHEM_SURFACE && p->molecular_weight <= 0) return CHEM_MISSING_MOLECULAR_WEIGHT;
    if (p->kind == CHEM_TAYLOR && p->coefficient_count > 16) return CHEM_INVALID;
    if (p->kind == CHEM_CALLBACK && !p->function) return CHEM_INVALID;
    return CHEM_OK;
}

int chem_rate_calculate(const ChemRate *rate, double temperature, double pressure,
                        double density, const double *parameters, size_t parameter_count, double *result)
{
    if (!result) return CHEM_INVALID;
    int status = validate_rate(rate, parameters, parameter_count);
    if (status) return status;
    *result = rate_value(rate, temperature, pressure, density, parameters, 1, 0);
    return isfinite(*result) ? CHEM_OK : CHEM_NONFINITE;
}

int chem_update_rates(ChemState *s)
{
    if (!s) return CHEM_INVALID;
    const ChemMechanism *m = s->mechanism;
    size_t cells = s->cells;
    for (size_t r = 0; r < m->reaction_count; ++r) {
        const ChemReaction *rx = &m->reactions[r];
        ChemRate configured = rx->rate;
        if (configured.kind == CHEM_CALLBACK) {
            configured.function = s->rate_functions[r];
            configured.context = s->rate_contexts[r];
        }
        const ChemRate *p = &configured;
        int status = validate_rate(p, s->parameters, m->parameter_count);
        if (status) return status;
        for (size_t c = 0; c < cells; ++c) {
            double t = s->temperatures[c], p_air = s->pressures[c];
            double density = s->explicit_air_density[c] ? s->air_densities[c] : p_air / (gas_constant * t);
            double value = rate_value(p, t, p_air, density, s->parameters, cells, c);
            for (unsigned i = 0; i < rx->air_density_power; ++i) value *= density;
            if (!isfinite(value)) return CHEM_NONFINITE;
            s->rates[r*cells+c] = value;
        }
    }
    return chem_update_constraint_constants(s);
}

void chem_forcing(const ChemState *s, const double *y, double *out, double *scratch)
{
    const ChemMechanism *m = s->mechanism;
    size_t cells = s->cells;
    memset(out, 0, m->species_count*cells*sizeof(*out));
    for (size_t r = 0; r < m->reaction_count; ++r) {
        const ChemReaction *rx = &m->reactions[r];
        const double *rate = s->rates+r*cells;
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) scratch[c] = rate[c];
        for (size_t k = 0; k < rx->reactant_count; ++k) {
            const double *v = y+(size_t)m->reactants[rx->reactant_start+k]*cells;
            #pragma omp simd
            for (size_t c = 0; c < cells; ++c) scratch[c] *= v[c];
        }
        for (size_t k = 0; k < rx->reactant_count; ++k) {
            int row=m->reactants[rx->reactant_start+k];
            if (s->algebraic[row]) continue;
            double *v = out+(size_t)row*cells;
            #pragma omp simd
            for (size_t c = 0; c < cells; ++c) v[c] -= scratch[c];
        }
        for (size_t k = 0; k < rx->product_count; ++k) {
            ChemProduct product = m->products[rx->product_start+k];
            if (s->algebraic[product.species]) continue;
            double *v = out+(size_t)product.species*cells, yield = product.yield;
            #pragma omp simd
            for (size_t c = 0; c < cells; ++c) v[c] += yield*scratch[c];
        }
    }
    chem_constraint_residuals(s,y,out);
}

void chem_negative_jacobian(ChemState *s, const double *y)
{
    const ChemMechanism *m = s->mechanism;
    size_t cells = s->cells, n = m->species_count;
    double *jac = s->matrix, *dr = s->scratch;
    memset(jac, 0, s->nonzeros*cells*sizeof(*jac));
    for (size_t r = 0; r < m->reaction_count; ++r) {
        const ChemReaction *rx = &m->reactions[r];
        const double *rate = s->rates+r*cells;
        for (size_t k = 0; k < rx->reactant_count; ++k) {
            size_t col = (size_t)m->reactants[rx->reactant_start+k];
            #pragma omp simd
            for (size_t c = 0; c < cells; ++c) dr[c] = rate[c];
            for (size_t j = 0; j < rx->reactant_count; ++j) if (j != k) {
                const double *v = y+(size_t)m->reactants[rx->reactant_start+j]*cells;
                #pragma omp simd
                for (size_t c = 0; c < cells; ++c) dr[c] *= v[c];
            }
            for (size_t j = 0; j < rx->reactant_count; ++j) {
                size_t row = (size_t)m->reactants[rx->reactant_start+j];
                if (s->algebraic[row]) continue;
                double *v = jac+(size_t)s->index[row*n+col]*cells;
                #pragma omp simd
                for (size_t c = 0; c < cells; ++c) v[c] += dr[c];
            }
            for (size_t j = 0; j < rx->product_count; ++j) {
                ChemProduct product = m->products[rx->product_start+j];
                if (s->algebraic[product.species]) continue;
                double *v = jac+(size_t)s->index[(size_t)product.species*n+col]*cells, yield = product.yield;
                #pragma omp simd
                for (size_t c = 0; c < cells; ++c) v[c] -= yield*dr[c];
            }
        }
    }
    chem_constraint_jacobian(s,y);
}

int chem_evaluate(const ChemState *s, size_t cell, double *forcing, double *jacobian)
{
    if (!s || cell >= s->cells || !forcing || !jacobian) return CHEM_INVALID;
    const ChemMechanism *m = s->mechanism;
    size_t n = m->species_count, cells = s->cells;
    memset(forcing, 0, n*sizeof(*forcing));
    memset(jacobian, 0, n*n*sizeof(*jacobian));
    for (size_t r = 0; r < m->reaction_count; ++r) {
        const ChemReaction *rx = &m->reactions[r];
        double rate = s->rates[r*cells+cell], flux = rate;
        for (size_t k = 0; k < rx->reactant_count; ++k)
            flux *= s->values[(size_t)m->reactants[rx->reactant_start+k]*cells+cell];
        for (size_t k = 0; k < rx->reactant_count; ++k) {
            int row=m->reactants[rx->reactant_start+k];
            if (!s->algebraic[row]) forcing[row] -= flux;
        }
        for (size_t k = 0; k < rx->product_count; ++k) {
            ChemProduct p = m->products[rx->product_start+k];
            if (!s->algebraic[p.species]) forcing[p.species] += p.yield*flux;
        }
        for (size_t k = 0; k < rx->reactant_count; ++k) {
            size_t col = (size_t)m->reactants[rx->reactant_start+k];
            double derivative = rate;
            for (size_t j = 0; j < rx->reactant_count; ++j) if (j != k)
                derivative *= s->values[(size_t)m->reactants[rx->reactant_start+j]*cells+cell];
            for (size_t j = 0; j < rx->reactant_count; ++j) {
                int row=m->reactants[rx->reactant_start+j];
                if (!s->algebraic[row]) jacobian[(size_t)row*n+col] += derivative;
            }
            for (size_t j = 0; j < rx->product_count; ++j) {
                ChemProduct p = m->products[rx->product_start+j];
                if (!s->algebraic[p.species]) jacobian[(size_t)p.species*n+col] -= p.yield*derivative;
            }
        }
    }
    for (size_t i=0; i<s->constraint_count; ++i) {
        const ChemConstraint *constraint=&s->constraints[i];
        size_t row=(size_t)constraint->algebraic_species;
        forcing[row]=chem_constraint_value(constraint,s->values,cells,cell,s->constraint_constants[i*cells+cell],
                                           jacobian+row*n,n);
        for (size_t col=0; col<n; ++col) jacobian[row*n+col] = -jacobian[row*n+col];
    }
    return CHEM_OK;
}
