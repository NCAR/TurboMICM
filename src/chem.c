#include "chem.h"

#include <math.h>
#include <string.h>

static const double gas_constant = (1.380649e-23 * 6.02214076e23);
static const double avogadro = 6.02214076e23;
static const double pi = 3.14159265358979323846;

static double branch_factor(double a, double b)
{
    double ratio = a / b;
    double log_ratio = log10(ratio);
    return a / (1.0 + ratio) * pow(0.41, 1.0 / (1.0 + log_ratio * log_ratio));
}

int chem_rate_constants(const ChemMechanism *m, ChemConditions c, const double *parameters, double *rates)
{
    if (!m || !rates || (m->parameter_count && !parameters) || !(c.temperature > 0) || !(c.pressure >= 0)) return 1;
    double t = c.temperature, air = c.pressure / (gas_constant * t);
    for (size_t i = 0; i < m->reaction_count; ++i) {
        const ChemReaction *r = &m->reactions[i];
        const double *p = r->p;
        double k = 0.0;
        switch (r->type) {
        case CHEM_ARRHENIUS:
            k = p[0] * exp(p[2] / t) * pow(t / p[3], p[1]) * (1.0 + p[4] * c.pressure);
            break;
        case CHEM_TROE:
        case CHEM_TERNARY: {
            double k0 = p[0] * exp(p[2] / t) * pow(t / 300.0, p[1]);
            double ki = p[3] * exp(p[5] / t) * pow(t / 300.0, p[4]);
            double ratio = k0 * air / ki;
            double log_ratio = log10(ratio);
            k = k0 * (r->type == CHEM_TROE ? air : 1.0) / (1.0 + ratio)
                * pow(p[6], p[7] / (p[7] + log_ratio * log_ratio));
            break;
        }
        case CHEM_TUNNELING:
            k = p[0] * exp(-p[1] / t + p[2] / (t * t * t));
            break;
        case CHEM_USER:
            k = parameters[r->parameter] * p[0];
            break;
        case CHEM_SURFACE: {
            double radius = parameters[r->parameter], count = parameters[r->parameter + 1];
            double speed = sqrt(8.0 * gas_constant * t / (pi * p[1]));
            k = 4.0 * count * pi * radius * radius / (radius / p[0] + 4.0 / (speed * p[2]));
            break;
        }
        case CHEM_BRANCH_ALKOXY:
        case CHEM_BRANCH_NITRATE: {
            double k0 = 2.0e-22 * avogadro * 1.0e-6 * exp(p[3]);
            double ref = 2.45e19 / avogadro * 1.0e6;
            double z = branch_factor(k0 * ref, 0.43 * pow(293.0 / 298.0, -8.0)) * (1.0 - p[2]) / p[2];
            double a = branch_factor(k0 * air, 0.43 * pow(t / 298.0, -8.0));
            double pre = p[0] * exp(-p[1] / t);
            k = pre * (r->type == CHEM_BRANCH_ALKOXY ? z / (z + a) : a / (a + z));
            break;
        }
        default: return 1;
        }
        if (!isfinite(k)) return 1;
        rates[i] = k;
    }
    return 0;
}

void chem_forcing(const ChemMechanism *m, ChemConditions c, const double *rates, const double *y, double *f)
{
    double air = c.pressure / (gas_constant * c.temperature);
    memset(f, 0, m->species_count * sizeof(*f));
    for (size_t i = 0; i < m->reaction_count; ++i) {
        const ChemReaction *r = &m->reactions[i];
        double rate = rates[i];
        for (size_t q = 0; q < r->reactant_count; ++q) {
            int s = m->reactants[r->reactant_begin + q].species;
            rate *= s < 0 ? air : y[s];
        }
        for (size_t q = 0; q < r->reactant_count; ++q) {
            int s = m->reactants[r->reactant_begin + q].species;
            if (s >= 0) f[s] -= rate;
        }
        for (size_t q = 0; q < r->product_count; ++q) {
            ChemTerm term = m->products[r->product_begin + q];
            f[term.species] += term.yield * rate;
        }
    }
}

void chem_jacobian(const ChemMechanism *m, ChemConditions c, const double *rates, const double *y, double *j)
{
    size_t n = m->species_count;
    double air = c.pressure / (gas_constant * c.temperature);
    memset(j, 0, n * n * sizeof(*j));
    for (size_t i = 0; i < m->reaction_count; ++i) {
        const ChemReaction *r = &m->reactions[i];
        for (size_t q = 0; q < r->reactant_count; ++q) {
            int col = m->reactants[r->reactant_begin + q].species;
            if (col < 0) continue;
            double derivative = rates[i];
            for (size_t k = 0; k < r->reactant_count; ++k) {
                if (k == q) continue;
                int s = m->reactants[r->reactant_begin + k].species;
                derivative *= s < 0 ? air : y[s];
            }
            for (size_t k = 0; k < r->reactant_count; ++k) {
                int row = m->reactants[r->reactant_begin + k].species;
                if (row >= 0) j[(size_t)row*n + (size_t)col] -= derivative;
            }
            for (size_t k = 0; k < r->product_count; ++k) {
                ChemTerm term = m->products[r->product_begin + k];
                j[(size_t)term.species*n + (size_t)col] += term.yield * derivative;
            }
        }
    }
}
