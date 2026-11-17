#ifndef RATES_H
#define RATES_H
#include <stddef.h>

typedef enum {
    CHEM_RATE_ARRHENIUS, CHEM_RATE_TROE, CHEM_RATE_TERNARY,
    CHEM_RATE_TUNNELING, CHEM_RATE_ALKOXY, CHEM_RATE_NITRATE,
    CHEM_RATE_SURFACE, CHEM_RATE_USER, CHEM_RATE_REVERSIBLE,
    CHEM_RATE_TAYLOR, CHEM_RATE_TERNARY_JPL19
} ChemRateType;
typedef struct {
    ChemRateType type;
    double p[21];
    size_t custom_index; /* Taylor uses this as its coefficient count. */
} ChemRate;
typedef struct { double temperature, pressure, air_density; } ChemConditions;
/* custom entries are addressed by the generated custom parameter labels. */
double chem_rate(const ChemRate *rate, ChemConditions conditions, const double *custom);
void chem_rates(size_t count, const ChemRate *rates, ChemConditions conditions,
                const double *custom, double *out);
/* Optional third-body powers multiply the computed rate by air density. */
void chem_rates_with_air(size_t count, const ChemRate *rates, const unsigned *air_orders,
                         ChemConditions conditions, const double *custom, double *out);
#endif
