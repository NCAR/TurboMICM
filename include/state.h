#ifndef CHEM_STATE_H
#define CHEM_STATE_H
#include <stddef.h>
#include "rates.h"

typedef enum {
    CHEM_STATE_OK, CHEM_STATE_INVALID, CHEM_STATE_SIZE, CHEM_STATE_UNKNOWN_NAME,
    CHEM_STATE_DUPLICATE_NAME, CHEM_STATE_NO_MEMORY
} ChemStateStatus;

typedef struct {
    size_t cells, species_count, parameter_count, rate_count;
    char **species_names, **parameter_names;
    ChemConditions *conditions;
    double *concentrations, *parameters, *rates;
    double relative_tolerance, *absolute_tolerances;
} ChemState;

typedef struct {
    size_t cells, species_count, parameter_count, rate_count;
    const char *const *species_names, *const *parameter_names;
    double relative_tolerance;
    const double *absolute_tolerances; /* NULL initializes to zero. */
} ChemStateConfig;

/* Start with ChemState s = {0}; init and clone replace an initialized state
   only after success. Names and all numeric arrays are owned by the state. */
ChemStateStatus chem_state_init(ChemState *state, const ChemStateConfig *config);
ChemStateStatus chem_state_clone(ChemState *destination, const ChemState *source);
void chem_state_move(ChemState *destination, ChemState *source);
void chem_state_destroy(ChemState *state);
/* A named column update accepts exactly one value per cell. */
ChemStateStatus chem_state_set_concentration(ChemState *state, const char *name,
                                             const double *values, size_t count);
ChemStateStatus chem_state_set_parameter(ChemState *state, const char *name,
                                         const double *values, size_t count);
typedef struct { const char *name; const double *values; size_t count; } ChemStateColumn;
/* Validate all columns before changing any values. */
ChemStateStatus chem_state_set_concentrations(ChemState *state, const ChemStateColumn *columns, size_t count);
ChemStateStatus chem_state_set_parameters(ChemState *state, const ChemStateColumn *columns, size_t count);
/* Direct access for one cell; returns NULL for an invalid name or cell. */
double *chem_state_concentration(ChemState *state, size_t cell, const char *name);
/* Copy a cell-major parameter table after validating both dimensions. */
ChemStateStatus chem_state_set_parameter_table(ChemState *state, const double *values,
                                               size_t cells, size_t parameters);
#endif
