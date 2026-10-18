/* MICM test/integration/test_analytical_backward_euler.cpp.
 *
 * One library configuration stands for each group of MICM builders that differ
 * only in what the library fixes (see test_analytical_rosenbrock_integration.c):
 *   backward_euler (L = 1, the library's default LU, Mozart) also stands for
 *       backard_euler_vector_1 [MICM's spelling], backward_euler_vector_mozart_1,
 *       backward_euler_vector_mozart_csc_1, backward_euler_vector_mozart_in_place_1;
 *   backward_euler_vector_L (L = 2, 3, 4) for backward_euler_vector_mozart_L,
 *       backward_euler_vector_mozart_csc_L, backward_euler_vector_mozart_in_place_L;
 *   backward_euler_vector_doolittle_L (L = 1 .. 4) for
 *       backward_euler_vector_doolittle_csc_L and backward_euler_vector_doolittle_in_place_L. */
#include "analytical_policy.h"

#define BE TMC_METHOD_BACKWARD_EULER
#define MOZART TMC_LU_MOZART
#define DOOLITTLE TMC_LU_DOOLITTLE

static const ap_config backward_euler = { "backward_euler", BE, 0, 1, MOZART };
static const ap_config backward_euler_vector_2 = { "backward_euler_vector_2", BE, 0, 2, MOZART };
static const ap_config backward_euler_vector_3 = { "backward_euler_vector_3", BE, 0, 3, MOZART };
static const ap_config backward_euler_vector_4 = { "backward_euler_vector_4", BE, 0, 4, MOZART };
static const ap_config backward_euler_vector_doolittle_1 = { "backward_euler_vector_doolittle_1", BE, 0, 1, DOOLITTLE };
static const ap_config backward_euler_vector_doolittle_2 = { "backward_euler_vector_doolittle_2", BE, 0, 2, DOOLITTLE };
static const ap_config backward_euler_vector_doolittle_3 = { "backward_euler_vector_doolittle_3", BE, 0, 3, DOOLITTLE };
static const ap_config backward_euler_vector_doolittle_4 = { "backward_euler_vector_doolittle_4", BE, 0, 4, DOOLITTLE };

/* backward_euler and backard_euler_vector_1 .. 4 */
static const ap_config *const vector_configs[] = { &backward_euler, &backward_euler_vector_2, &backward_euler_vector_3, &backward_euler_vector_4, NULL };

/* every distinct configuration MICM's lists reach */
static const ap_config *const all_configs[] = {
  &backward_euler, &backward_euler_vector_2, &backward_euler_vector_3, &backward_euler_vector_4,
  &backward_euler_vector_doolittle_1, &backward_euler_vector_doolittle_2, &backward_euler_vector_doolittle_3, &backward_euler_vector_doolittle_4,
  NULL
};

TMC_TEST(troe)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_troe(*c, 1e-6);
  return TMC_TEST_PASS;
}

TMC_TEST(troe_stiff)
{
  for (const ap_config *const *c = vector_configs; *c; ++c)
    ap_stiff_troe(*c, AP_STIFF_TROE_TOL);
  return TMC_TEST_PASS;
}

TMC_TEST(photolysis)
{
  for (const ap_config *const *c = vector_configs; *c; ++c)
    ap_photolysis(*c, 1e-3);
  return TMC_TEST_PASS;
}

TMC_TEST(photolysis_stiff)
{
  for (const ap_config *const *c = vector_configs; *c; ++c)
    ap_stiff_photolysis(*c, 1e-3);
  return TMC_TEST_PASS;
}

TMC_TEST(ternary_chemical_activation)
{
  for (const ap_config *const *c = vector_configs; *c; ++c)
    ap_ternary_chemical_activation(*c, 1e-5);
  return TMC_TEST_PASS;
}

TMC_TEST(ternary_chemical_activation_stiff)
{
  for (const ap_config *const *c = vector_configs; *c; ++c)
    ap_stiff_ternary_chemical_activation(*c, 1e-2);
  return TMC_TEST_PASS;
}

TMC_TEST(tunneling)
{
  for (const ap_config *const *c = vector_configs; *c; ++c)
    ap_tunneling(*c, 1e-3);
  return TMC_TEST_PASS;
}

TMC_TEST(tunneling_stiff)
{
  for (const ap_config *const *c = vector_configs; *c; ++c)
    ap_stiff_tunneling(*c, 1e-3);
  return TMC_TEST_PASS;
}

TMC_TEST(arrhenius)
{
  for (const ap_config *const *c = vector_configs; *c; ++c)
    ap_arrhenius(*c, 1e-3);
  return TMC_TEST_PASS;
}

TMC_TEST(arrhenius_stiff)
{
  for (const ap_config *const *c = vector_configs; *c; ++c)
    ap_stiff_arrhenius(*c, 1e-3);
  return TMC_TEST_PASS;
}

TMC_TEST(branched)
{
  for (const ap_config *const *c = vector_configs; *c; ++c)
    ap_branched(*c, 1e-5);
  return TMC_TEST_PASS;
}

TMC_TEST(branched_stiff)
{
  for (const ap_config *const *c = vector_configs; *c; ++c)
    ap_stiff_branched(*c, 1e-2);
  return TMC_TEST_PASS;
}

TMC_TEST(surface_rxn)
{
  ap_surface_rxn(&backward_euler, 0.05);
  return TMC_TEST_PASS;
}

TMC_TEST(hires)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_hires(*c, 1e-1);
  return TMC_TEST_PASS;
}

TMC_TEST(oregonator)
{
  /* One large step (H = 30 tau) makes Newton converge to the wrong attractor;
   * 18000 sub-steps per output interval let backward Euler track the limit
   * cycle with first-order error, about 3e-3 relative. */
  const int substeps = 18000;
  const double tolerance = 5e-3;
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_oregonator(*c, tolerance, substeps, AP_OREGONATOR_SOLVER_RTOL);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.Troe", troe),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.TroeSuperStiffButAnalytical", troe_stiff),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.Photolysis", photolysis),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.PhotolysisSuperStiffButAnalytical", photolysis_stiff),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.TernaryChemicalActivation", ternary_chemical_activation),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.TernaryChemicalActivationSuperStiffButAnalytical", ternary_chemical_activation_stiff),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.Tunneling", tunneling),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.TunnelingSuperStiffButAnalytical", tunneling_stiff),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.Arrhenius", arrhenius),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.ArrheniusSuperStiffButAnalytical", arrhenius_stiff),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.Branched", branched),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.BranchedSuperStiffButAnalytical", branched_stiff),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.SurfaceRxn", surface_rxn),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.HIRES", hires),
  TMC_CASE("analytical_backward_euler:AnalyticalExamples.Oregonator", oregonator),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_analytical_backward_euler(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
