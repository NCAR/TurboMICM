#ifndef RATES_H
#define RATES_H
#include <stddef.h>

typedef enum {
    CHEM_RATE_ARRHENIUS, CHEM_RATE_TROE, CHEM_RATE_TERNARY,
    CHEM_RATE_TUNNELING, CHEM_RATE_ALKOXY, CHEM_RATE_NITRATE,
    CHEM_RATE_SURFACE, CHEM_RATE_USER, CHEM_RATE_REVERSIBLE,
    CHEM_RATE_TAYLOR, CHEM_RATE_TERNARY_JPL19, CHEM_RATE_CALLBACK
} ChemRateType;
typedef struct {
    ChemRateType type;
    double p[21];
    size_t custom_index; /* Taylor uses this as its coefficient count. */
} ChemRate;
typedef struct { double temperature, pressure, air_density; } ChemConditions;
typedef struct { double constant, temperature, pressure, air_density; } ChemParameterizedSpecies;
typedef struct { size_t begin, count; } ChemRateMultiplier;
typedef double (*ChemRateCallback)(ChemConditions conditions, void *context);
typedef struct { ChemRateCallback function; void *context; } ChemRateBinding;
/* custom entries are addressed by the generated custom parameter labels. */
double chem_rate(const ChemRate *rate, ChemConditions conditions, const double *custom);
void chem_rates(size_t count, const ChemRate *rates, ChemConditions conditions,
                const double *custom, double *out);
/* Parameterized reactants are excluded from state and multiply their rate.
   Descriptors refer to a shared array of linear functions of conditions. */
void chem_rates_parameterized(size_t count, const ChemRate *rates,
    const ChemRateMultiplier *multipliers, const ChemParameterizedSpecies *functions,
    ChemConditions conditions, const double *custom, double *out);
/* Callback rate custom_index selects a binding; unavailable bindings yield NaN. */
void chem_rates_parameterized_with_callbacks(size_t count, const ChemRate *rates,
    const ChemRateMultiplier *multipliers, const ChemParameterizedSpecies *functions,
    const ChemRateBinding *bindings, size_t binding_count,
    ChemConditions conditions, const double *custom, double *out);
/* Optional third-body powers multiply the computed rate by air density. */
void chem_rates_with_air(size_t count, const ChemRate *rates, const unsigned *air_orders,
                         ChemConditions conditions, const double *custom, double *out);
#endif
