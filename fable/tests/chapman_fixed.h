/* MICM's hard-coded Chapman solver (test/regression/RosenbrockChapman/
 * chapman_ode_solver.hpp) and the Chapman mechanism its regression tests build
 * for the general solver (util.hpp), for the regression_test_* ports.
 *
 * The hard-coded solver orders its species M, Ar, CO2, H2O, N2, O1D, O, O2, O3
 * and its rate constants N2_O1D_1, O1D_O2_1, O_O3_1, M_O_O2_1 (Arrhenius), then
 * O2_1, O3_1, O3_2 (the photolysis rates jO2, jO3a, jO3b). Its Jacobian has 23
 * elements, numbered column by column. */
#ifndef TMC_TEST_CHAPMAN_FIXED_H
#define TMC_TEST_CHAPMAN_FIXED_H

#include "harness.h"

enum
{
  CF_N_SPECIES = 9,
  CF_N_RATE_CONSTANTS = 7,
  CF_N_PHOTOLYSIS = 3,
  CF_N_ELEMENTS = 23
};

extern const char *const cf_species_names[CF_N_SPECIES];
extern const char *const cf_photolysis_names[CF_N_PHOTOLYSIS];
/* (row, column) of each of the 23 elements, in species indices */
extern const int cf_elements[CF_N_ELEMENTS][2];

/* ChapmanODESolver::Force: the forcing of each species */
void cf_force(const double *rate_constants, const double *number_densities, double *force);
/* ChapmanODESolver::DforceDy: d force / d y in the 23-element numbering */
void cf_dforce_dy(const double *rate_constants, const double *number_densities, double *jacobian);
/* ChapmanODESolver::UpdateState: the rate constants from one cell's
 * conditions and photolysis rates (the Arrhenius ones through the library's
 * rate-constant evaluation, as MICM's go through its shared CalculateArrhenius) */
void cf_update_state(double temperature, double pressure, const double *photolysis_rates, double *rate_constants);

typedef struct cf_result
{
  double y[CF_N_SPECIES];
  double final_time;
  long function_calls, jacobian_updates, number_of_steps, accepted, rejected, decompositions, solves;
} cf_result;

/* ChapmanODESolver::Solve: its three-stage Rosenbrock integration of one cell
 * from time_start to time_end, with its settings (relative tolerance 1e-4,
 * absolute 1e-3, at most 100 steps, first step 1e-6 s). */
cf_result cf_solve(double time_start, double time_end, const double *y, const double *rate_constants);

/* util.hpp's GetChapmanSolver: the Chapman mechanism (phase M, Ar, CO2, H2O,
 * N2, O1D, O, O2, O3; reactions jO2, jO3a, jO3b, then the four Arrhenius
 * reactions) built on the given builder, which is freed. Returns NULL after
 * recording a failure. */
tmc_solver *cf_chapman_solver(tmc_builder *builder);

#endif
