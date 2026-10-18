/*
 * tmc support: functions the ported MICM tests need that tmc.h lacks. Each
 * mirrors a MICM function named in its comment; none changes the behavior of
 * the functions in tmc.h.
 */
#ifndef TMC_SUPPORT_H
#define TMC_SUPPORT_H

#include "tmc.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Lambda rate constants by label ----------------------------------------
 * MICM's Solver::GetLambdaRateConstantByName: the solver's own copy of the
 * lambda parameters with this label, which the caller may change (a new
 * function takes effect at the next tmc_solver_update_state_parameters).
 * Returns NULL and sets TMC_ERROR_SOLVER / TMC_SOLVER_RATE_CONSTANT_NOT_FOUND
 * when no reaction has a lambda rate constant with this label. */
tmc_lambda_parameters *tmc_support_solver_lambda_rate_constant(tmc_solver *s, const char *label, tmc_error *err);

/* ---- Unique names in a given order -----------------------------------------
 * MICM's System::UniqueNames(f): name i of the result is order(names, n, i,
 * context), where names are the system's n unique state-variable names in
 * declaration order (tmc_system_unique_names); a NULL order keeps that order.
 * Writes up to max_names names and returns n. */
typedef const char *(*tmc_support_name_order)(const char (*names)[TMC_NAME_MAX], int n_names, int i, void *context);
int tmc_support_system_unique_names_ordered(const tmc_system *s, tmc_support_name_order order, void *context,
                                            char (*names)[TMC_NAME_MAX], int max_names);

#ifdef __cplusplus
}
#endif

#endif /* TMC_SUPPORT_H */
