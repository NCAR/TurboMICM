/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * C translation and independent CPU implementation. */
#ifndef CHEM_H
#define CHEM_H
#include <stddef.h>

#define CHEM_AVOGADRO 6.02214076e23
#define CHEM_GAS_CONSTANT (1.380649e-23 * CHEM_AVOGADRO)
#define CHEM_PI 3.14159265358979323846

typedef enum {
    CHEM_OK, CHEM_INVALID, CHEM_MEMORY, CHEM_NOT_FOUND, CHEM_DIMENSION,
    CHEM_SINGULAR, CHEM_MAX_STEPS, CHEM_STEP_TOO_SMALL, CHEM_CONSTRAINT_FAILURE, CHEM_UNCONVERGED, CHEM_NAN, CHEM_INF,
    CHEM_EMPTY_REACTANTS, CHEM_EMPTY_PRODUCTS, CHEM_BAD_STOICHIOMETRY, CHEM_BAD_EQUILIBRIUM
} ChemStatus;
const char *chem_status_string(ChemStatus status);

typedef struct { double temperature, pressure, air_density; } ChemConditions;
void chem_ideal_air_density(ChemConditions *conditions);

typedef enum {
    CHEM_ARRHENIUS, CHEM_TROE, CHEM_TERNARY, CHEM_TERNARY_JPL19,
    CHEM_BRANCHED, CHEM_TUNNELING, CHEM_TAYLOR, CHEM_REVERSIBLE,
    CHEM_USER, CHEM_SURFACE, CHEM_CALLBACK, CHEM_RATE_TYPE_COUNT
} ChemRateType;
typedef double (*ChemRateFunction)(const ChemConditions *, void *);
typedef struct { double A, B, C, D, E; } ChemArrhenius;
typedef struct {
    double k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C, Fc, N;
} ChemFalloff;
typedef struct {
    double k0_A, k0_B, k0_C, k0_D, kinf_A, kinf_B, kinf_C, kinf_D;
    double kint_A, kint_B, kint_C, kint_D, Fc, N;
} ChemJpl19;
typedef struct { int branch; double X, Y, a0; size_t n; double k0, z; } ChemBranched;
typedef struct { double A, B, C; } ChemTunneling;
typedef struct { ChemArrhenius arrhenius; double coefficients[16]; size_t count; } ChemTaylor;
typedef struct { double A, C, reverse_rate; } ChemReversible;
typedef struct { double scaling; } ChemUser;
typedef struct { double diffusion, molecular_weight, probability, speed_factor; } ChemSurface;
typedef struct {
    ChemRateType type;
    const char *label;
    size_t custom_index;
    union {
        ChemArrhenius arrhenius;
        ChemFalloff falloff;
        ChemJpl19 jpl19;
        ChemBranched branched;
        ChemTunneling tunneling;
        ChemTaylor taylor;
        ChemReversible reversible;
        ChemUser user;
        ChemSurface surface;
        struct { ChemRateFunction function; void *context; } callback;
    } data;
} ChemRate;
ChemRate chem_rate_default(ChemRateType type);
ChemStatus chem_rate_prepare(ChemRate *rate);
double chem_rate_evaluate(const ChemRate *rate, const ChemConditions *conditions, const double *custom);

/* Species with a parameter function are excluded from the solution vector.
 * Their values multiply the rate constant once per reactant occurrence. */
typedef struct {
    const char *name;
    double absolute_tolerance;
    double molecular_weight, diffusion;
    ChemRateFunction parameter;
    void *parameter_context;
} ChemSpecies;
typedef struct { size_t species; double yield; } ChemProduct;
typedef struct {
    size_t reactant_count, product_count;
    const size_t *reactants;
    const ChemProduct *products;
    ChemRate rate;
} ChemReaction;
typedef struct ChemMechanism ChemMechanism;
ChemStatus chem_mechanism_create(size_t species_count, const ChemSpecies *species,
                                 size_t reaction_count, const ChemReaction *reactions,
                                 ChemMechanism **out);
void chem_mechanism_free(ChemMechanism *mechanism);
size_t chem_species_count(const ChemMechanism *mechanism);
size_t chem_reaction_count(const ChemMechanism *mechanism);
size_t chem_custom_count(const ChemMechanism *mechanism);
const char *chem_species_name(const ChemMechanism *mechanism, size_t variable);
const char *chem_custom_name(const ChemMechanism *mechanism, size_t parameter);
const ChemRate *chem_rate_definition(const ChemMechanism *,size_t sorted_reaction);
size_t chem_parameterized_reactant_count(const ChemMechanism *,size_t sorted_reaction);
double chem_reaction_multiplier(const ChemMechanism *,size_t sorted_reaction,const ChemConditions *);
size_t chem_rate_offset(const ChemMechanism *mechanism, ChemRateType type);
ChemStatus chem_species_index(const ChemMechanism *, const char *, size_t *);
ChemStatus chem_custom_index(const ChemMechanism *, const char *, size_t *);
ChemStatus chem_set_rate_callback(ChemMechanism *, const char *, ChemRateFunction, void *);

/* All fields are column-major: field[column * cells + cell]. This layout
 * gives the cell loop contiguous storage for OpenMP SIMD. */
typedef struct {
    const ChemMechanism *mechanism;
    size_t cells;
    ChemConditions *conditions;
    double *values, *custom, *rates, *absolute_tolerance;
    double *scratch; /* private cell workspace; one state per concurrent solve */
    double relative_tolerance;
} ChemState;
ChemStatus chem_state_create(const ChemMechanism *, size_t cells, ChemState **out);
ChemStatus chem_state_copy(const ChemState *, ChemState **out);
void chem_state_free(ChemState *);
ChemStatus chem_state_assign(ChemState **destination,const ChemState *source);
void chem_state_move(ChemState **destination,ChemState **source);
/* Named setters require one value per cell. Fill functions explicitly broadcast. */
ChemStatus chem_state_fill_concentration(ChemState *,const char *,double value);
ChemStatus chem_state_fill_custom(ChemState *,const char *,double value);
typedef struct { const char *name;const double *values;size_t count; } ChemNamedValues;
ChemStatus chem_state_set_concentrations(ChemState *,size_t count,const ChemNamedValues *);
ChemStatus chem_state_set_custom_parameters(ChemState *,size_t count,const ChemNamedValues *);
/* Matrix input is row-major [cell*parameters+parameter]. */
ChemStatus chem_state_set_custom_matrix(ChemState *,size_t cells,size_t parameters,const double *);
double *chem_state_concentration_data(ChemState *,const char *);
const double *chem_state_concentration_const_data(const ChemState *,const char *);
ChemStatus chem_state_set_concentration(ChemState *, const char *, const double *, size_t count);
ChemStatus chem_state_set_custom(ChemState *, const char *, const double *, size_t count);
void chem_update_rates(ChemState *);
/* Kinetic accumulation leaves existing values and algebraic rows intact. */
void chem_add_kinetic_forcing(const ChemState *,const double *values,double *forcing);
void chem_subtract_kinetic_jacobian(const ChemState *,const double *values,double *jacobian);
int chem_reactions_use_species(const ChemMechanism *,const char *name);
int chem_kinetic_has_element(const ChemMechanism *,size_t row,size_t column);
void chem_forcing(const ChemState *, const double *values, double *forcing);
/* Dense column-major blocks indexed [(row * species + column) * cells + cell]. */
void chem_jacobian(const ChemState *, const double *values, double *jacobian);
double chem_parameter_air_density(const ChemConditions *, void *);
#endif
