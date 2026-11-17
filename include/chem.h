#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>

/* All mechanism data are immutable and may be shared between solver states.
 * A caller may initialize these structs in C or link a generated definition. */
typedef enum {
    CHEM_ARRHENIUS, CHEM_TROE, CHEM_TERNARY, CHEM_TERNARY_JPL19,
    CHEM_BRANCHED, CHEM_TUNNELING, CHEM_TAYLOR, CHEM_REVERSIBLE,
    CHEM_USER, CHEM_SURFACE
} ChemRateKind;

typedef struct {
    const char *name;
    double initial, absolute_tolerance, molecular_weight, diffusion;
} ChemSpecies;

typedef struct { int species; double yield; } ChemProduct;

typedef struct {
    ChemRateKind kind;
    int parameter;                 /* index, or -1; surface uses index and index+1 */
    double a, b, c, d, e;
    double k0_a, k0_b, k0_c, k0_d;
    double kinf_a, kinf_b, kinf_c, kinf_d;
    double kint_a, kint_b, kint_c, kint_d;
    double fc, n;
    double branch_x, branch_y, branch_a0, branch_n;
    int nitrate_branch;
    double diffusion, molecular_weight, probability;
    size_t coefficient_count;
    double coefficients[16];
} ChemRate;

typedef struct {
    size_t reactant_start, reactant_count, product_start, product_count;
    unsigned air_density_power;
    ChemRate rate;
} ChemReaction;

typedef struct {
    size_t species_count, reaction_count, parameter_count;
    const ChemSpecies *species;
    const ChemReaction *reactions;
    const int *reactants;
    const ChemProduct *products;
    const char *const *parameter_names;
    const double *parameter_defaults;
    double initial_temperature, initial_pressure;
} ChemMechanism;

typedef struct {
    unsigned long internal_steps, accepted, rejected, function_calls;
    unsigned long jacobian_updates, decompositions, solves;
    double final_time;
} ChemStats;

/* The arrays a and c store the strictly lower triangular matrices in row order.
 * The maximum of six stages matches the supported Rosenbrock parameter sets. */
typedef struct {
    size_t stages;
    double a[15], c[15], m[6], e[6], gamma[6];
    int new_function[6];
    double order, relative_tolerance, h_min, h_max, h_start;
    double factor_min, factor_max, safety, rejection_factor;
    unsigned max_steps;
} ChemRosenbrock;

typedef struct ChemState ChemState;

ChemRosenbrock chem_rosenbrock_three_stage(void);
ChemState *chem_state_create(const ChemMechanism *mechanism, size_t cells);
/* Share the immutable symbolic linear-solver plan with an existing state. */
ChemState *chem_state_create_like(const ChemState *state, size_t cells);
void chem_state_free(ChemState *state);
size_t chem_cell_count(const ChemState *state);
/* Each returned array has chem_cell_count() entries. */
double *chem_species_values(ChemState *state, size_t species);
double *chem_parameter_values(ChemState *state, size_t parameter);
double *chem_temperatures(ChemState *state);
double *chem_pressures(ChemState *state);
int chem_find_species(const ChemMechanism *mechanism, const char *name);
int chem_find_parameter(const ChemMechanism *mechanism, const char *name);
int chem_update_rates(ChemState *state);
/* Negative Jacobian, in dense row-major form for one cell, and forcing, for
 * diagnostics and tests. The solver itself uses sparse storage. */
int chem_evaluate(const ChemState *state, size_t cell, double *forcing, double *negative_jacobian);
int chem_solve(ChemState *state, const ChemRosenbrock *parameters, double dt, ChemStats *stats);
const char *chem_error_string(int code);

enum { CHEM_OK = 0, CHEM_INVALID = 1, CHEM_MEMORY = 2, CHEM_SINGULAR = 3,
       CHEM_TOO_MANY_STEPS = 4, CHEM_STEP_TOO_SMALL = 5, CHEM_NONFINITE = 6 };

#endif
