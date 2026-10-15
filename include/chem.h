#ifndef CHEM_H
#define CHEM_H
#include <stddef.h>

/* All rates use concentrations in mol m^-3, temperature K, pressure Pa. */
typedef enum {
    CHEM_ARRHENIUS, CHEM_TROE, CHEM_TERNARY, CHEM_TERNARY_JPL19,
    CHEM_TUNNELING, CHEM_BRANCHED_ALKOXY, CHEM_BRANCHED_NITRATE,
    CHEM_TAYLOR, CHEM_REVERSIBLE, CHEM_USER, CHEM_SURFACE, CHEM_CALLBACK
} ChemRateType;
typedef double (*ChemRateFunction)(double temperature, double pressure,
                                   double air_density, const double *custom, void *context);
typedef struct {
    ChemRateType type;
    /* Arrhenius/Taylor: A B C D E, then polynomial coefficients.
       Troe/ternary: k0 A B C, kinf A B C, Fc N.
       JPL19: k0 A B C D, kinf A B C D, kint A B C D, Fc N.
       Tunnel: A B C. Branched: X Y a0 n. Reversible: A C kr.
       User: scaling. Surface: diffusion, molecular weight, probability. */
    double p[24];
    int parameter, coefficient_count, third_body_count;
    ChemRateFunction function;
    void *context;
} ChemRate;
typedef struct { int species; double yield; } ChemProduct;
typedef struct {
    size_t reactant_count, product_count;
    const int *reactants;
    const ChemProduct *products;
    ChemRate rate;
} ChemReaction;
typedef struct {
    size_t species_count, reaction_count, parameter_count;
    const char *const *species_names;
    const char *const *parameter_names;
    const double *initial_concentrations, *initial_parameters, *absolute_tolerances;
    double temperature, pressure;
    const ChemReaction *reactions;
} ChemMechanismData;
typedef struct ChemMechanism ChemMechanism;
typedef struct ChemState ChemState;
typedef enum { CHEM_OK, CHEM_INVALID, CHEM_MEMORY, CHEM_MAX_STEPS,
               CHEM_SMALL_STEP, CHEM_NONFINITE, CHEM_SINGULAR } ChemStatus;
typedef struct {
    size_t internal_steps, accepted, rejected, function_calls,
           jacobian_updates, decompositions, solves;
    double final_time;
} ChemStats;
typedef struct {
    int stages, new_function[6];
    size_t max_steps;
    double a[15], c[15], m[6], e[6], alpha[6], gamma[6], order;
    double factor_min, factor_max, rejection_factor, safety;
    double h_min, h_max, h_start;
} ChemRosenbrock;

/* Data and all pointer members must outlive the mechanism. Construction validates
   indices and precomputes the sparse fill pattern. Separate states are independent. */
ChemStatus chem_mechanism_create(const ChemMechanismData *, ChemMechanism **);
void chem_mechanism_free(ChemMechanism *);
const ChemMechanismData *chem_mechanism_data(const ChemMechanism *);
size_t chem_jacobian_size(const ChemMechanism *);
ChemStatus chem_state_create(const ChemMechanism *, size_t cells, ChemState **);
void chem_state_free(ChemState *);
/* Species-major and parameter-major arrays, with cells contiguous. */
double *chem_concentrations(ChemState *);
double *chem_parameters(ChemState *);
double *chem_temperatures(ChemState *);
double *chem_pressures(ChemState *);
double *chem_absolute_tolerances(ChemState *);
void chem_relative_tolerance(ChemState *, double);
ChemStatus chem_update_rates(ChemState *);
const double *chem_rate_constants(const ChemState *);
void chem_forcing(const ChemState *, const double *y, double *forcing);
void chem_negative_jacobian(const ChemState *, const double *y, double *jacobian);
/* Diagnostic row/column to storage index, SIZE_MAX if structural zero. */
size_t chem_jacobian_index(const ChemMechanism *, size_t row, size_t column);
ChemRosenbrock chem_rosenbrock_parameters(int stages); /* 2, 3, 4, -4 (DA), 6 (DA) */
ChemStatus chem_solve(ChemState *, const ChemRosenbrock *, double dt, ChemStats *);
const char *chem_status_string(ChemStatus);
double chem_rate_evaluate(const ChemRate *, double t, double p, double air, const double *custom);
#endif
