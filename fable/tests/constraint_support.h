/* Support for the ports of MICM's constraint tests (constraint_set,
 * equilibrium, linear, constraint_initialization, and the DAE integration
 * tests): systems and reactions in a line, the C equivalent of MICM's
 * stand-alone ConstraintSet, its residual and Jacobian evaluators, and MICM's
 * finite-difference Jacobian checks (micm/util/jacobian_verification.hpp). */
#ifndef TMC_TEST_CONSTRAINT_SUPPORT_H
#define TMC_TEST_CONSTRAINT_SUPPORT_H

#include "harness.h"
#include "tmc_constraint_ext.h"

/* TMC_ASSERT_TRUE for tests that own allocations: records the failure and
 * jumps to the test's `done:` label, which frees them. */
#define CST_REQUIRE(cond)                                                          \
  do                                                                               \
  {                                                                                \
    if (!(cond))                                                                   \
    {                                                                              \
      tmc_test_fail(__FILE__, __LINE__, "required: %s", #cond);                    \
      goto done;                                                                   \
    }                                                                              \
  } while (0)

/* CST_REQUIRE that a build succeeded, recording the builder's message. */
#define CST_REQUIRE_BUILT(solver, err)                                             \
  do                                                                               \
  {                                                                                \
    if (!(solver))                                                                 \
    {                                                                              \
      tmc_test_fail(__FILE__, __LINE__, "build failed: %s", (err).message);        \
      goto done;                                                                   \
    }                                                                              \
  } while (0)

/* EXPECT_THROW of a MicmException: the call returned the error category and code. */
#define CST_EXPECT_ERROR(rc, err, category_, code_)                                \
  do                                                                               \
  {                                                                                \
    int cst_rc_ = (rc);                                                            \
    if (cst_rc_ != (category_) || (err).category != (category_) || (err).code != (code_)) \
      tmc_test_fail(__FILE__, __LINE__, "expected error %d/%d, got %d (%d/%d: %s)", (int)(category_), (int)(code_), cst_rc_, (int)(err).category, (err).code, (err).message); \
  } while (0)

/* ---- systems, reactions, constraints ------------------------------------ */

/* StoichSpecies(Species(name), coefficient) */
tmc_stoich_species cst_stoich(const char *name, double coefficient);
/* A system whose gas phase holds the species, in order. */
tmc_system cst_system(const tmc_species *species, int n);
/* A reaction; MICM's ChemicalReactionBuilder with reactants, products, and a rate constant. */
tmc_reaction cst_reaction(const tmc_species *reactants, int n_reactants, const tmc_stoich_species *products, int n_products, tmc_rate_constant rate_constant);
/* ArrheniusRateConstantParameters{ .A_ = A, .B_ = 0, .C_ = 0 } (D and E at their defaults) */
tmc_rate_constant cst_arrhenius(double A);
/* Equilibrium constraints with reference temperature and gas constant at
 * MICM's defaults; the constraint must be valid (failure is recorded). */
tmc_constraint cst_equilibrium(const char *name, const tmc_species *algebraic, const tmc_stoich_species *reactants, int n_reactants, const tmc_stoich_species *products, int n_products, double K_HLC_ref, double delta_H);
tmc_constraint cst_linear(const char *name, const tmc_species *algebraic, const tmc_stoich_species *terms, int n_terms, double constant);

/* A solver from the system, reactions, and constraints, as MICM's
 * CpuSolverBuilder builds one: vector width L (1 for MICM's Matrix and
 * standard sparse ordering, L for VectorMatrix<L>), the builder's default
 * Doolittle LU, and the reordering given. NULL with err set on failure. */
tmc_solver *cst_build(tmc_rosenbrock_parameters parameters, const tmc_system *system, const tmc_reaction *reactions, int n_reactions, const tmc_constraint *constraints, int n_constraints, int L, bool reorder_state, tmc_error *err);

/* K_HLC_ref exp((delta_H / R) (1 / T - 1 / 298.15)), R = k_B N_A as in
 * MICM's constants.hpp (the tests' ComputeEquilibriumConstant). */
double cst_equilibrium_constant(double K_HLC_ref, double delta_H, double temperature);

/* ---- MICM's stand-alone ConstraintSet ------------------------------------ */

/* ConstraintSet(constraints, variable_map): the C library holds a constraint
 * set inside a solver, so this builds one from the species in variable-map
 * order (no kinetics, no reordering, and no unused-species check, which the
 * constraint set does not make) with vector width L. The species must keep
 * the given indices; NULL with err set on failure. */
tmc_solver *cst_constraint_set(const char *const *species, int n_species, const tmc_constraint *constraints, int n_constraints, int L, tmc_error *err);

/* ConstraintSet::NonZeroJacobianElements, newly allocated (*n pairs). */
tmc_index_pair *cst_nonzero_elements(const tmc_solver *s, int *n);
bool cst_has_pair(const tmc_index_pair *pairs, int n, int row, int col);

/* SetConstraintFunctions(state_parameter_indices): the column of a
 * constraint parameter in the solver's parameter map (-1 if absent). */
int cst_parameter_index(const tmc_solver *s, const char *label);
/* A state-parameter matrix (n_cells x the solver's parameter count) at a value. */
tmc_dense *cst_parameters(const tmc_solver *s, int n_cells, double value);
/* A Jacobian with the solver's sparsity pattern, n_blocks blocks, zero. */
tmc_sparse *cst_jacobian(const tmc_solver *s, int n_blocks);
/* jacobian[block][row][col]; a structurally zero element records a failure
 * (MICM throws on access) and reads as NaN. */
double cst_jacobian_at(const tmc_sparse *jacobian, int block, int row, int col, const char *file, int line);
#define CST_JAC(jacobian, block, row, col) cst_jacobian_at((jacobian), (block), (row), (col), __FILE__, __LINE__)

/* ConstraintSet::AddForcingTerms and SubtractJacobianTerms over the solver's
 * constraints: the residuals G are assigned to the algebraic rows of forcing,
 * and dG/dy is subtracted from the algebraic rows of the Jacobian. */
void cst_add_forcing_terms(const tmc_solver *s, const tmc_dense *variables, const tmc_dense *parameters, tmc_dense *forcing);
void cst_subtract_jacobian_terms(const tmc_solver *s, const tmc_dense *variables, const tmc_dense *parameters, tmc_sparse *jacobian);

/* ---- finite-difference Jacobians (micm/util/jacobian_verification.hpp) --- */

typedef void (*cst_forcing_function)(void *context, const tmc_dense *variables, tmc_dense *forcing);

/* cst_add_forcing_terms as a cst_forcing_function; the context is a cst_forcing_context. */
typedef struct cst_forcing_context
{
  const tmc_solver *solver;
  const tmc_dense *parameters;
} cst_forcing_context;
void cst_constraint_forcing(void *context, const tmc_dense *variables, tmc_dense *forcing);

typedef struct cst_jacobian_comparison
{
  bool passed;
  double max_abs_error, max_rel_error;
  int worst_block, worst_row, worst_col;
  double worst_analytical, worst_fd;
} cst_jacobian_comparison;

/* MICM's defaults in double precision */
#define CST_FD_PERTURBATION 1.0e-8
#define CST_FD_ATOL 1.0e-7
#define CST_FD_RTOL 1.0e-7
#define CST_FD_SPARSITY_THRESHOLD 1.0e-6

/* FiniteDifferenceJacobian: element [block][row * n_species + col] is
 * d forcing_row / d x_col by central differences, or forward differences
 * where x - h would be negative; h = perturbation max(1, |x|). */
tmc_dense *cst_finite_difference_jacobian(cst_forcing_function forcing, void *context, const tmc_dense *variables, int n_species, double perturbation);
/* CompareJacobianToFiniteDifference: over the analytical Jacobian's pattern,
 * |-J - fd| <= atol + rtol max(|J|, |fd|) (the analytical Jacobian holds -df/dx). */
cst_jacobian_comparison cst_compare_jacobian_to_finite_difference(const tmc_sparse *analytical, const tmc_dense *fd, int n_species, double atol, double rtol);
/* CheckJacobianSparsityCompleteness: outside the pattern, |fd| <= threshold. */
cst_jacobian_comparison cst_check_jacobian_sparsity_completeness(const tmc_sparse *analytical, const tmc_dense *fd, int n_species, double threshold);

#endif
