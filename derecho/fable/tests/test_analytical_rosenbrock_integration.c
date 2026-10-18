/* MICM test/integration/test_analytical_rosenbrock.cpp.
 *
 * MICM runs each analytical test on up to 35 builders. One library
 * configuration stands for each group of builders that differ only in what the
 * library fixes (CSR against CSC order, in place against separate L and U):
 *   rosenbrock_2stage, _3stage, _4stage, _4stage_da, _6stage_da
 *       the five parameter sets on MICM's standard layout (L = 1) with the
 *       library's default LU (Mozart);
 *   rosenbrock_3stage (L = 1) also stands for rosenbrock_vector_1,
 *       rosenbrock_standard_mozart, rosenbrock_vector_mozart_1,
 *       rosenbrock_vector_mozart_csc_1, rosenbrock_vector_mozart_in_place_1;
 *   rosenbrock_vector_L (L = 2, 3, 4) for rosenbrock_vector_mozart_L,
 *       rosenbrock_vector_mozart_csc_L, rosenbrock_vector_mozart_in_place_L;
 *   rosenbrock_vector_doolittle_L (L = 1 .. 4) for rosenbrock_standard_doolittle
 *       (L = 1), rosenbrock_vector_doolittle_csc_L, and
 *       rosenbrock_vector_doolittle_in_place_L.
 * The tolerances are the ones MICM passes for each builder. */
#include "analytical_policy.h"

#define R TMC_METHOD_ROSENBROCK
#define MOZART TMC_LU_MOZART
#define DOOLITTLE TMC_LU_DOOLITTLE

static const ap_config rosenbrock_2stage = { "rosenbrock_2stage", R, 2, 1, MOZART };
static const ap_config rosenbrock_3stage = { "rosenbrock_3stage", R, 3, 1, MOZART };
static const ap_config rosenbrock_4stage = { "rosenbrock_4stage", R, 4, 1, MOZART };
static const ap_config rosenbrock_4stage_da = { "rosenbrock_4stage_da", R, 44, 1, MOZART };
static const ap_config rosenbrock_6stage_da = { "rosenbrock_6stage_da", R, 6, 1, MOZART };
static const ap_config rosenbrock_vector_2 = { "rosenbrock_vector_2", R, 3, 2, MOZART };
static const ap_config rosenbrock_vector_3 = { "rosenbrock_vector_3", R, 3, 3, MOZART };
static const ap_config rosenbrock_vector_4 = { "rosenbrock_vector_4", R, 3, 4, MOZART };
static const ap_config rosenbrock_vector_doolittle_1 = { "rosenbrock_vector_doolittle_1", R, 3, 1, DOOLITTLE };
static const ap_config rosenbrock_vector_doolittle_2 = { "rosenbrock_vector_doolittle_2", R, 3, 2, DOOLITTLE };
static const ap_config rosenbrock_vector_doolittle_3 = { "rosenbrock_vector_doolittle_3", R, 3, 3, DOOLITTLE };
static const ap_config rosenbrock_vector_doolittle_4 = { "rosenbrock_vector_doolittle_4", R, 3, 4, DOOLITTLE };

/* every distinct configuration MICM's lists reach */
static const ap_config *const all_configs[] = {
  &rosenbrock_2stage,   &rosenbrock_3stage,   &rosenbrock_4stage,   &rosenbrock_4stage_da,
  &rosenbrock_6stage_da, &rosenbrock_vector_2, &rosenbrock_vector_3, &rosenbrock_vector_4,
  &rosenbrock_vector_doolittle_1, &rosenbrock_vector_doolittle_2, &rosenbrock_vector_doolittle_3, &rosenbrock_vector_doolittle_4,
  NULL
};

/* the five parameter sets on the standard layout */
static const ap_config *const parameter_sets[] = {
  &rosenbrock_2stage, &rosenbrock_3stage, &rosenbrock_4stage, &rosenbrock_4stage_da, &rosenbrock_6stage_da, NULL
};

/* MICM gives the two-stage method a looser tolerance in several tests */
static double two_stage_or(const ap_config *c, double two_stage, double others)
{
  return c->stages == 2 ? two_stage : others;
}

TMC_TEST(troe)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_troe(*c, AP_TROE_TOL);
  return TMC_TEST_PASS;
}

TMC_TEST(troe_stiff)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_stiff_troe(*c, AP_STIFF_TROE_TOL);
  return TMC_TEST_PASS;
}

TMC_TEST(photolysis)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_photolysis(*c, AP_PHOTOLYSIS_TOL);
  return TMC_TEST_PASS;
}

TMC_TEST(photolysis_stiff)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_stiff_photolysis(*c, AP_STIFF_PHOTOLYSIS_TOL);
  return TMC_TEST_PASS;
}

TMC_TEST(ternary_chemical_activation)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_ternary_chemical_activation(*c, AP_TERNARY_TOL);
  return TMC_TEST_PASS;
}

TMC_TEST(ternary_chemical_activation_stiff)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_stiff_ternary_chemical_activation(*c, 2e-3);
  return TMC_TEST_PASS;
}

TMC_TEST(tunneling)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_tunneling(*c, two_stage_or(*c, 2e-5, AP_TUNNELING_TOL));
  return TMC_TEST_PASS;
}

TMC_TEST(tunneling_stiff)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_stiff_tunneling(*c, 1e-4);
  return TMC_TEST_PASS;
}

TMC_TEST(arrhenius)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_arrhenius(*c, two_stage_or(*c, 4e-6, AP_ARRHENIUS_TOL));
  return TMC_TEST_PASS;
}

TMC_TEST(arrhenius_stiff)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
  {
    const double tol = (*c)->stages == 2 ? 1e-4 : (*c)->stages == 6 ? 1e-5 : 2e-5;
    ap_stiff_arrhenius(*c, tol);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(branched)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_branched(*c, two_stage_or(*c, 1e-10, AP_BRANCHED_TOL));
  return TMC_TEST_PASS;
}

TMC_TEST(branched_stiff)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_stiff_branched(*c, 2e-3);
  return TMC_TEST_PASS;
}

TMC_TEST(surface_rxn)
{
  static const double tolerances[5] = { 1e-2, 1e-5, 1e-6, 1e-5, 1e-7 };
  for (int i = 0; parameter_sets[i]; ++i)
    ap_surface_rxn(parameter_sets[i], tolerances[i]);
  return TMC_TEST_PASS;
}

TMC_TEST(robertson)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_robertson(*c, 1e-6);
  return TMC_TEST_PASS;
}

TMC_TEST(e5)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_e5(*c, 1e-10);
  return TMC_TEST_PASS;
}

TMC_TEST(oregonator)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
    ap_oregonator(*c, 1e-6, 1, AP_OREGONATOR_SOLVER_RTOL);
  return TMC_TEST_PASS;
}

TMC_TEST(hires)
{
  for (const ap_config *const *c = all_configs; *c; ++c)
  {
    const double tol = (*c)->stages == 2 || (*c)->stages == 44 || (*c)->stages == 6 ? 1e-6 : 1e-7;
    ap_hires(*c, tol);
  }
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.Troe", troe),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.TroeSuperStiffButAnalytical", troe_stiff),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.Photolysis", photolysis),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.PhotolysisSuperStiffButAnalytical", photolysis_stiff),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.TernaryChemicalActivation", ternary_chemical_activation),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.TernaryChemicalActivationSuperStiffButAnalytical", ternary_chemical_activation_stiff),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.Tunneling", tunneling),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.TunnelingSuperStiffButAnalytical", tunneling_stiff),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.Arrhenius", arrhenius),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.ArrheniusSuperStiffButAnalytical", arrhenius_stiff),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.Branched", branched),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.BranchedSuperStiffButAnalytical", branched_stiff),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.SurfaceRxn", surface_rxn),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.Robertson", robertson),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.E5", e5),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.Oregonator", oregonator),
  TMC_CASE("analytical_rosenbrock_integration:AnalyticalExamples.HIRES", hires),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_analytical_rosenbrock_integration(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
