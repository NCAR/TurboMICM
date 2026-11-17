#include "chem.h"
#include "internal.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char *chem_version(void) { return "0.1.0"; }
unsigned chem_version_major(void) { return 0; }
unsigned chem_version_minor(void) { return 1; }
unsigned chem_version_patch(void) { return 0; }
unsigned chem_version_tweak(void) { return 0; }

static const double gas_constant = (1.380649e-23 * 6.02214076e23);
static const double avogadro = 6.02214076e23;
static const double pi = 3.14159265358979323846;

ChemValidation chem_validate_mechanism(const ChemMechanism *m, char *message, size_t capacity)
{
    if (message && capacity) message[0] = '\0';
    if (!m || !m->species_count || !m->species_names || (m->reaction_count && !m->reactions)) return CHEM_INVALID_MECHANISM;
    for (size_t i = 0; i < m->reaction_count; ++i) {
        const ChemReaction *r = &m->reactions[i];
        if (r->type == CHEM_SURFACE) {
            int index = r->reactant_count && m->reactants ? m->reactants[r->reactant_begin].species : -1;
            const char *name = index >= 0 && (size_t)index < m->species_count ? m->species_names[index] : "?";
            if (!(r->p[0] > 0)) {
                if (message && capacity) snprintf(message, capacity, "Diffusion coefficient for species '%s' is not defined", name);
                return CHEM_MISSING_DIFFUSION;
            }
            if (!(r->p[1] > 0)) {
                if (message && capacity) snprintf(message, capacity, "Species '%s' molecular weight is not defined", name);
                return CHEM_MISSING_MOLECULAR_WEIGHT;
            }
        }
    }
    return CHEM_VALID;
}

ptrdiff_t chem_parameter_index(const ChemMechanism *m, const char *name)
{
    if (!m || !name || (m->parameter_count && !m->parameter_names)) return -1;
    for (size_t i = 0; i < m->parameter_count; ++i)
        if (!strcmp(m->parameter_names[i], name)) return (ptrdiff_t)i;
    return -1;
}

static double density(ChemConditions c)
{
    return c.has_air_density ? c.air_density : c.pressure / (gas_constant * c.temperature);
}

double chem_parameterized_value(ChemParameterization p, ChemConditions c)
{
    return p.constant + p.temperature*c.temperature + p.pressure*c.pressure +
           (p.air_density == 0 ? 0 : p.air_density*density(c));
}

static double reactant_value(const ChemMechanism *m, int species, ChemConditions c, const double *y)
{
    if (species >= 0) return y[species];
    if (species == -1) return density(c);
    size_t index = (size_t)(-(species+2));
    return index < m->parameterized_species_count && m->parameterized_species ?
        chem_parameterized_value(m->parameterized_species[index],c) : NAN;
}

static double branch_factor(double a, double b)
{
    double ratio = a / b;
    double log_ratio = log10(ratio);
    return a / (1.0 + ratio) * pow(0.41, 1.0 / (1.0 + log_ratio * log_ratio));
}

int chem_rate_constants(const ChemMechanism *m, ChemConditions c, const double *parameters, double *rates)
{
    return chem_rate_constants_with_callbacks(m, c, parameters, NULL, rates);
}

int chem_rate_constants_with_callbacks(const ChemMechanism *m, ChemConditions c, const double *parameters,
                                       const ChemCallback *callbacks, double *rates)
{
    if (!m || !rates || (m->parameter_count && !parameters) || !(c.temperature > 0) || !(c.pressure >= 0)) return 1;
    double t = c.temperature, air = density(c);
    if (air < 0 || !isfinite(air)) return 1;
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
        case CHEM_TERNARY_JPL19: {
            double k0 = p[0]*exp(p[2]/t)*pow(t/p[3],p[1]);
            double ki = p[4]*exp(p[6]/t)*pow(t/p[7],p[5]);
            double kint = p[8]*exp(p[10]/t)*pow(t/p[11],p[9]);
            double ratio = k0*air/ki, log_ratio = log10(ratio);
            double kf = k0*air/(1.0+ratio)*pow(p[12],p[13]/(p[13]+log_ratio*log_ratio));
            k = kf + kint*(1.0-kf/ki);
            break;
        }
        case CHEM_TUNNELING:
            k = p[0] * exp(-p[1] / t + p[2] / (t * t * t));
            break;
        case CHEM_REVERSIBLE:
            k = p[0] * exp(p[1] / t) * p[2];
            break;
        case CHEM_TAYLOR: {
            int count = (int)p[5];
            if (count < 1 || count > 16 || p[5] != count) return 1;
            double polynomial = 0.0, power = 1.0;
            for (int j = 0; j < count; ++j, power *= t) polynomial += p[6+j] * power;
            k = polynomial * p[0] * exp(p[2]/t) * pow(t/p[3], p[1]) * (1.0+p[4]*c.pressure);
            break;
        }
        case CHEM_CALLBACK:
            if (!callbacks || r->parameter < 0 || (size_t)r->parameter >= m->parameter_count || !callbacks[r->parameter].function) return 1;
            k = callbacks[r->parameter].function(c, callbacks[r->parameter].context);
            break;
        case CHEM_USER:
            if (r->parameter < 0 || (size_t)r->parameter >= m->parameter_count) return 1;
            k = parameters[r->parameter] * p[0];
            break;
        case CHEM_SURFACE: {
            if (r->parameter < 0 || (size_t)r->parameter >= m->parameter_count || m->parameter_count-(size_t)r->parameter < 2) return 1;
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
        for (size_t q = 0; q < r->reactant_count; ++q) {
            int species = m->reactants[r->reactant_begin+q].species;
            if (species < 0) k *= reactant_value(m, species, c, NULL);
        }
        if (!isfinite(k)) return 1;
        rates[i] = k;
    }
    return 0;
}

void chem_forcing(const ChemMechanism *m, ChemConditions c, const double *rates, const double *y, double *f)
{
    (void)c; /* parameterized reactants are incorporated in rates */
    memset(f, 0, m->species_count * sizeof(*f));
    for (size_t i = 0; i < m->reaction_count; ++i) {
        const ChemReaction *r = &m->reactions[i];
        double rate = rates[i];
        for (size_t q = 0; q < r->reactant_count; ++q) {
            int s = m->reactants[r->reactant_begin + q].species;
            if (s >= 0) rate *= y[s];
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
    memset(j, 0, m->species_count * m->species_count * sizeof(*j));
    chem_jacobian_add(m,c,rates,y,j);
}

void chem_jacobian_add(const ChemMechanism *m, ChemConditions c, const double *rates, const double *y, double *j)
{
    (void)c;
    size_t n = m->species_count;
    for (size_t i = 0; i < m->reaction_count; ++i) {
        const ChemReaction *r = &m->reactions[i];
        for (size_t q = 0; q < r->reactant_count; ++q) {
            int col = m->reactants[r->reactant_begin + q].species;
            if (col < 0) continue;
            double derivative = rates[i];
            for (size_t k = 0; k < r->reactant_count; ++k) {
                if (k == q) continue;
                int s = m->reactants[r->reactant_begin + k].species;
                if (s >= 0) derivative *= y[s];
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

int chem_process_add_forcing(const ChemMechanism *m,ChemConditions c,const double *rates,const double *y,
                             const unsigned char *algebraic,double *forcing)
{
    if(!m || !rates || !y || !forcing || m->species_count>SIZE_MAX/sizeof(double)) return 1;
    double *temp=malloc((m->species_count?m->species_count:1)*sizeof(double)); if(!temp) return 1;
    chem_forcing(m,c,rates,y,temp);
    for(size_t i=0;i<m->species_count;++i) if(!algebraic || !algebraic[i]) forcing[i]+=temp[i];
    free(temp); return 0;
}
int chem_process_subtract_jacobian(const ChemMechanism *m,ChemConditions c,const double *rates,const double *y,
                                   const unsigned char *algebraic,double *jacobian)
{
    if(!m || !rates || !y || !jacobian || (m->species_count && m->species_count>SIZE_MAX/m->species_count)) return 1;
    size_t n=m->species_count, count=n*n; if(count>SIZE_MAX/sizeof(double)) return 1;
    double *temp=malloc((count?count:1)*sizeof(double)); if(!temp) return 1;
    chem_jacobian(m,c,rates,y,temp);
    for(size_t i=0;i<n;++i) if(!algebraic || !algebraic[i]) for(size_t j=0;j<n;++j) jacobian[i*n+j]-=temp[i*n+j];
    free(temp); return 0;
}
int chem_process_jacobian_pattern(const ChemMechanism *m,unsigned char *pattern)
{
    if(!m || !pattern || (m->species_count && m->species_count>SIZE_MAX/m->species_count)) return 1;
    size_t n=m->species_count; memset(pattern,0,n*n);
    for(size_t i=0;i<m->reaction_count;++i) {
        const ChemReaction *r=&m->reactions[i];
        for(size_t j=0;j<r->reactant_count;++j) {
            int col=m->reactants[r->reactant_begin+j].species;
            if(col<0) continue;
            if((size_t)col>=n) return 1;
            for(size_t k=0;k<r->reactant_count;++k) {
                int row=m->reactants[r->reactant_begin+k].species;
                if(row>=0) { if((size_t)row>=n) return 1; pattern[(size_t)row*n+(size_t)col]=1; }
            }
            for(size_t k=0;k<r->product_count;++k) {
                int row=m->products[r->product_begin+k].species;
                if(row<0 || (size_t)row>=n) return 1;
                pattern[(size_t)row*n+(size_t)col]=1;
            }
        }
    }
    return 0;
}
