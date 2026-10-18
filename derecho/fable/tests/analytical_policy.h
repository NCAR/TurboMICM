/* MICM test/integration/analytical_policy.hpp and analytical_surface_rxn_policy.hpp.
 *
 * MICM instantiates each analytical test for a list of builder policies: the
 * five Rosenbrock parameter sets on the standard layout, the three-stage method
 * on VectorMatrix<1..4>, each with the Doolittle or the Mozart LU, in CSR or
 * CSC order, in place or not, and backward Euler likewise. In the library the
 * layout is the width L, the LU algorithm a run-time choice, and the sparse
 * order and in-place factorization fixed (CSR, in place), so each MICM builder
 * maps onto an ap_config, and the CSC and in-place builders onto the same
 * configurations as their CSR counterparts. */
#ifndef TMC_TEST_ANALYTICAL_POLICY_H
#define TMC_TEST_ANALYTICAL_POLICY_H

#include "harness.h"

typedef struct ap_config
{
  const char *name; /* MICM's builder, for failure messages */
  tmc_method method;
  int stages; /* Rosenbrock: 2, 3, 4, 44 (four-stage differential algebraic), 6 (six-stage differential algebraic) */
  int L;      /* MICM's Matrix is 1; VectorMatrix<L> is L */
  tmc_lu_algorithm lu;
} ap_config;

/* A builder for the configuration; the caller frees it (tmc_builder_free) or
 * hands it to one of the tests below, which free it. */
tmc_builder *ap_create_builder(const ap_config *c);

/* MICM's default tolerances for each test */
#define AP_TROE_TOL 1e-10
#define AP_STIFF_TROE_TOL 1e-5
#define AP_PHOTOLYSIS_TOL 2e-6
#define AP_STIFF_PHOTOLYSIS_TOL 2e-5
#define AP_TERNARY_TOL 1e-8
#define AP_STIFF_TERNARY_TOL 1e-6
#define AP_TUNNELING_TOL 1e-8
#define AP_STIFF_TUNNELING_TOL 1e-6
#define AP_ARRHENIUS_TOL 1e-9
#define AP_STIFF_ARRHENIUS_TOL 1e-6
#define AP_BRANCHED_TOL 1e-13
#define AP_STIFF_BRANCHED_TOL 1e-6
#define AP_SURFACE_RXN_TOL 1e-8
#define AP_ROBERTSON_TOL 1e-8
#define AP_OREGONATOR_TOL 1e-6
#define AP_OREGONATOR_SOLVER_RTOL 1e-9
#define AP_HIRES_TOL 1e-8
#define AP_E5_TOL 1e-8

/* One call of each MICM TestAnalytical* function on one configuration. The
 * A -> B -> C tests take an absolute tolerance on the concentrations; the
 * surface reaction, Robertson, Oregonator, and E5 tests a relative one; HIRES
 * an absolute one. Failures are recorded through the harness with the
 * configuration's name. */
void ap_troe(const ap_config *c, double abs_tol);
void ap_stiff_troe(const ap_config *c, double abs_tol);
void ap_photolysis(const ap_config *c, double abs_tol);
void ap_stiff_photolysis(const ap_config *c, double abs_tol);
void ap_ternary_chemical_activation(const ap_config *c, double abs_tol);
void ap_stiff_ternary_chemical_activation(const ap_config *c, double abs_tol);
void ap_tunneling(const ap_config *c, double abs_tol);
void ap_stiff_tunneling(const ap_config *c, double abs_tol);
void ap_arrhenius(const ap_config *c, double abs_tol);
void ap_stiff_arrhenius(const ap_config *c, double abs_tol);
void ap_branched(const ap_config *c, double abs_tol);
void ap_stiff_branched(const ap_config *c, double abs_tol);
void ap_surface_rxn(const ap_config *c, double tol);
void ap_robertson(const ap_config *c, double rel_tol);
/* substeps_per_output is 1 for a solver that controls its own step size;
 * backward Euler needs many (MICM uses 18000). */
void ap_oregonator(const ap_config *c, double rel_tol, int substeps_per_output, double solver_rtol);
void ap_hires(const ap_config *c, double abs_tol);
void ap_e5(const ap_config *c, double rel_tol);

#endif
