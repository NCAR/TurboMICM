/*
 * Additions to the external-model API that MICM's external-model tests need
 * (test/integration/test_external_model_constraints.cpp):
 *
 * - name maps built by hand, as those tests build std::unordered_map index
 *   maps and hand them to a model's callbacks without a solver;
 * - the constraint half of the forcing and the Jacobian (MICM's
 *   ConstraintBundle::AddForcingTerms and SubtractJacobianTerms), which the
 *   tests call after the rate half (tmc_solver_add_forcing_terms and
 *   tmc_solver_subtract_jacobian_terms).
 */
#ifndef TMC_EXTERNAL_EXT_H
#define TMC_EXTERNAL_EXT_H

#include "tmc.h"

#ifdef __cplusplus
extern "C" {
#endif

/* A name map in which names[i] has index i. Free it with tmc_name_map_free. */
tmc_name_map *tmc_name_map_create(const char *const *names, int n);
void tmc_name_map_free(tmc_name_map *map);

/* forcing[row] = G(y) for every algebraic row: the built-in constraints, then
 * the active external models' residuals, with the state's parameters. The
 * rows are assigned, as MICM does. */
void tmc_solver_add_constraint_terms(const tmc_solver *s, const tmc_state *state, const tmc_dense *variables, tmc_dense *forcing);
/* jacobian -= dG/dy for the same constraints, in the LU pattern. */
void tmc_solver_subtract_constraint_jacobian_terms(const tmc_solver *s, const tmc_state *state, const tmc_dense *variables, tmc_sparse *jacobian);

#ifdef __cplusplus
}
#endif

#endif /* TMC_EXTERNAL_EXT_H */
