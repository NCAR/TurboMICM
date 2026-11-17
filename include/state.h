#ifndef CHEM_STATE_H
#define CHEM_STATE_H
#include "chem.h"
#include "matrix.h"

typedef struct {
    size_t rate_count, variable_count, parameter_count;
    const char *const *variable_names, *const *parameter_names;
    double relative_tolerance;
    const double *absolute_tolerances;
} ChemStateParameters;
/* Name arrays and their strings are borrowed and must outlive the state and
 * any copies; matrices, conditions, and tolerance values are owned. Initialize
 * destinations to zero before init, copy, or move. */
typedef struct {
    size_t cells;
    ChemStateParameters parameters;
    ChemMatrix variables, custom_parameters, rate_constants;
    ChemConditions *conditions;
    double *absolute_tolerances;
} ChemState;
typedef struct { const char *name; const double *values; size_t count; } ChemNamedValues;
int chem_state_init(ChemState *, const ChemStateParameters *, size_t cells, size_t width);
void chem_state_free(ChemState *);
int chem_state_copy(ChemState *destination, const ChemState *source);
void chem_state_move(ChemState *destination, ChemState *source);
ptrdiff_t chem_state_variable_index(const ChemState *, const char *);
ptrdiff_t chem_state_parameter_index(const ChemState *, const char *);
double *chem_state_variable(ChemState *, size_t cell, const char *name);
const double *chem_state_const_variable(const ChemState *, size_t cell, const char *name);
int chem_state_set_concentration(ChemState *, const char *, const double *values, size_t count);
int chem_state_set_concentrations(ChemState *, const ChemNamedValues *, size_t count);
int chem_state_set_parameter(ChemState *, const char *, const double *values, size_t count);
int chem_state_set_parameters(ChemState *, const ChemNamedValues *, size_t count);
int chem_state_set_parameter_rows(ChemState *, const double *const *rows, const size_t *sizes, size_t count);
int chem_state_update_rate_constants(ChemState *, const ChemMechanism *, const ChemCallback *callbacks);
int chem_state_apply_phase_tolerances(ChemState *, ChemPhase);
int chem_state_update_constraint_parameters(ChemState *, ChemConstraintSet);
#endif
