/* MICM test/unit/process/rate_constant/: the ten rate-constant executables
 * (arrhenius, branched, lambda, reversible, surface, taylor_series,
 * ternary_chemical_activation, troe, tunneling, user_defined).
 *
 * MICM's tests call its Calculate* functions on parameter structs; here each
 * parameter struct becomes a tmc_rate_constant and the value comes from
 * tmc_rate_constant_calculate. The expected values are computed in the test
 * with the formulas and tolerances MICM's tests use (a tolerance of 0 is
 * MICM's exact equality in double precision). */
#include "harness.h"

#define AVOGADRO_CONSTANT 6.02214076e23
#define BOLTZMANN_CONSTANT 1.380649e-23
#define GAS_CONSTANT (BOLTZMANN_CONSTANT * AVOGADRO_CONSTANT)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double calculate(tmc_rate_constant rc, tmc_conditions conditions, const double *custom)
{
  return tmc_rate_constant_calculate(&rc, &conditions, custom);
}

/* ---- arrhenius_rate_constant ---- */

static void arrhenius_cases(void)
{
  const double tolerance = 1e-13;
  tmc_conditions conditions = { .temperature = 301.24 };

  tmc_arrhenius_parameters zero = tmc_arrhenius_defaults();
  double k = calculate(tmc_rate_arrhenius(zero), conditions, NULL);
  double expected = 1;
  TMC_EXPECT_NEAR(k, expected, tolerance * expected);

  tmc_arrhenius_parameters parameters = tmc_arrhenius_defaults();
  parameters.A = expected;
  k = calculate(tmc_rate_arrhenius(parameters), conditions, NULL);
  TMC_EXPECT_NEAR(k, expected, tolerance * expected);

  parameters.A = 2.2e-10;
  k = calculate(tmc_rate_arrhenius(parameters), conditions, NULL);
  expected = 2.2e-10;
  TMC_EXPECT_NEAR(k, expected, tolerance * expected);

  /* O + HO2 -> OH + O2 */
  parameters.A = 3e-11;
  parameters.C = -200;
  k = calculate(tmc_rate_arrhenius(parameters), conditions, NULL);
  expected = 3e-11 * exp(-200 / 301.24);
  TMC_EXPECT_NEAR(k, expected, tolerance * expected);

  /* OH + HCl -> H2O + Cl */
  parameters.A = 2.6e-12;
  parameters.C = -350;
  k = calculate(tmc_rate_arrhenius(parameters), conditions, NULL);
  expected = 2.6e-12 * exp(-350 / 301.24);
  TMC_EXPECT_NEAR(k, expected, tolerance * expected);
}

TMC_TEST(arrhenius_calculate_with_system)
{
  arrhenius_cases();
  return TMC_TEST_PASS;
}

TMC_TEST(arrhenius_calculate_with_prescribed_arguments)
{
  /* MICM's two Arrhenius cases differ only in how A = 1 is set; both set it
   * as a parameter here. */
  arrhenius_cases();
  return TMC_TEST_PASS;
}

/* ---- branched_rate_constant ---- */

static double branched_expected(tmc_branch branch, double temperature)
{
  double air_dens_n_cm3 = 42.2 * AVOGADRO_CONSTANT * 1.0e-6;
  double a = 2.0e-22 * exp(2) * 2.45e19;
  double b = 0.43 * pow((293.0 / 298.0), -8.0);
  double z = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
  a = 2.0e-22 * exp(2) * air_dens_n_cm3;
  b = 0.43 * pow((temperature / 298.0), -8.0);
  double A = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2)));
  return branch == TMC_BRANCH_ALKOXY ? 1.2 * exp(-204.3 / temperature) * (z / (z + A))
                                     : 1.2 * exp(-204.3 / temperature) * (A / (z + A));
}

static void branched_case(tmc_branch branch)
{
  double temperature = 301.24;
  tmc_conditions conditions = { .temperature = temperature, .air_density = 42.2 };
  tmc_branched_parameters params = { .branch = branch, .X = 1.2, .Y = 204.3, .a0 = 1.0e-3, .n = 2 };
  /* MICM derives k0 and z from n and a0 when it builds its rate-constant
   * store (the test's ComputeDerivedFields); the library derives them inside
   * the calculation, so the parameters are passed as given. */
  double k = calculate(tmc_rate_branched(params), conditions, NULL);
  double expected = branched_expected(branch, temperature);
  TMC_EXPECT_NEAR(k, expected, 1e-13 * expected);
}

TMC_TEST(branched_alkoxy)
{
  branched_case(TMC_BRANCH_ALKOXY);
  return TMC_TEST_PASS;
}

TMC_TEST(branched_nitrate)
{
  branched_case(TMC_BRANCH_NITRATE);
  return TMC_TEST_PASS;
}

/* ---- lambda_rate_constant ---- */

static double twice_temperature(const tmc_conditions *conditions, void *context)
{
  (void)context;
  return 2.0 * conditions->temperature;
}

TMC_TEST(lambda_calculate_with_system)
{
  tmc_rate_constant rc = tmc_rate_lambda("test lambda", twice_temperature, NULL);
  TMC_EXPECT_STREQ(rc.p.lambda.label, "test lambda");
  tmc_conditions conditions = { .temperature = 300.0 };
  double k = calculate(rc, conditions, NULL);
  double expected = 600.0;
  TMC_EXPECT_NEAR(k, expected, 1e-13 * expected);
  return TMC_TEST_PASS;
}

/* ---- reversible_rate_constant ---- */

TMC_TEST(reversible_default_constructor)
{
  tmc_reversible_parameters params = tmc_reversible_defaults();
  TMC_EXPECT_EQ(params.A, 1.0);
  TMC_EXPECT_EQ(params.C, 0.0);
  TMC_EXPECT_EQ(params.k_r, 0.0);
  tmc_conditions conditions = { .temperature = 298.15 };
  /* default: A = 1, C = 0, k_r = 0, so k = 1 * exp(0 / T) * 0 = 0 */
  double k = calculate(tmc_rate_reversible(params), conditions, NULL);
  TMC_EXPECT_NEAR(k, 0.0, 1e-13);
  return TMC_TEST_PASS;
}

TMC_TEST(reversible_calculate_with_system)
{
  tmc_conditions conditions = { .temperature = 301.24 };
  tmc_reversible_parameters parameters = tmc_reversible_defaults();
  parameters.A = 1.14e-2;
  parameters.C = 2300.0;
  parameters.k_r = 0.32;
  double k = calculate(tmc_rate_reversible(parameters), conditions, NULL);
  double K_eq = 1.14e-2 * exp(2300.0 / 301.24);
  double expected = K_eq * 0.32;
  TMC_EXPECT_NEAR(k, expected, 1e-13 * expected);
  return TMC_TEST_PASS;
}

TMC_TEST(reversible_calculate_with_prescribed_arguments)
{
  double temperature = 298.15;
  tmc_reversible_parameters parameters = tmc_reversible_defaults();
  parameters.A = 1.14e-2;
  parameters.C = 2300.0;
  parameters.k_r = 0.32;
  tmc_conditions conditions = { .temperature = temperature };
  double k = calculate(tmc_rate_reversible(parameters), conditions, NULL);
  double K_eq = 1.14e-2 * exp(2300.0 / temperature);
  double expected = K_eq * 0.32;
  TMC_EXPECT_NEAR(k, expected, 1e-13 * expected);
  return TMC_TEST_PASS;
}

TMC_TEST(reversible_temperature_dependence)
{
  tmc_reversible_parameters parameters = tmc_reversible_defaults();
  parameters.A = 1.5;
  parameters.C = 2000.0;
  parameters.k_r = 0.5;
  const double temperatures[3] = { 273.15, 298.15, 323.15 };
  for (int i = 0; i < 3; ++i)
  {
    tmc_conditions conditions = { .temperature = temperatures[i] };
    double k = calculate(tmc_rate_reversible(parameters), conditions, NULL);
    double expected = 1.5 * exp(2000.0 / temperatures[i]) * 0.5;
    TMC_EXPECT_NEAR(k, expected, 1e-13 * expected);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(reversible_zero_reverse_rate)
{
  tmc_reversible_parameters parameters = tmc_reversible_defaults();
  parameters.A = 5.0;
  parameters.C = 1000.0;
  parameters.k_r = 0.0;
  tmc_conditions conditions = { .temperature = 298.15 };
  double k = calculate(tmc_rate_reversible(parameters), conditions, NULL);
  TMC_EXPECT_NEAR(k, 0.0, 1e-13);
  return TMC_TEST_PASS;
}

TMC_TEST(reversible_equilibrium_constant_only)
{
  tmc_reversible_parameters parameters = tmc_reversible_defaults();
  parameters.A = 3.5;
  parameters.C = 1800.0;
  parameters.k_r = 1.0;
  double temperature = 298.15;
  tmc_conditions conditions = { .temperature = temperature };
  double k = calculate(tmc_rate_reversible(parameters), conditions, NULL);
  double K_eq = 3.5 * exp(1800.0 / temperature);
  TMC_EXPECT_NEAR(k, K_eq, 1e-13 * K_eq);
  return TMC_TEST_PASS;
}

/* ---- surface_rate_constant ---- */

/* The phase species MICM's tests construct as PhaseSpecies(foo, D). */
static tmc_phase_species phase_species(const tmc_species *species, bool with_diffusion, double diffusion_coefficient)
{
  tmc_phase phase = tmc_phase_create("gas");
  if (with_diffusion)
    tmc_phase_add_species_with_diffusion(&phase, species, diffusion_coefficient);
  else
    tmc_phase_add_species(&phase, species);
  tmc_phase_species ps = phase.species[0];
  tmc_phase_free(&phase);
  return ps;
}

static void surface_case(double probability)
{
  tmc_species foo = tmc_species_create_with_property("foo", TMC_PROPERTY_MOLECULAR_WEIGHT, 0.025);
  double foo_diffusion_coefficient = 2.3e2;
  tmc_phase_species foo_gas_species = phase_species(&foo, true, foo_diffusion_coefficient);
  tmc_rate_constant rc = tmc_rate_surface("foo", &foo_gas_species, probability);
  TMC_EXPECT_EQ(tmc_rate_constant_custom_parameter_count(&rc), 2);
  char labels[2][TMC_NAME_MAX];
  TMC_EXPECT_EQ(tmc_rate_constant_custom_parameter_labels(&rc, labels), 2);
  TMC_EXPECT_STREQ(labels[0], "foo.effective radius [m]");
  TMC_EXPECT_STREQ(labels[1], "foo.particle number concentration [# m-3]");

  tmc_conditions conditions = { .temperature = 273.65 };
  const double custom_params[2] = { 1.0e-7, 2.5e6 }; /* effective radius [m], particle number concentration [# m-3] */
  double k = calculate(rc, conditions, custom_params);
  double expected = 4.0 * 2.5e6 * M_PI * pow(1.0e-7, 2) /
                    (1.0e-7 / 2.3e2 + 4.0 / (sqrt(8.0 * GAS_CONSTANT * 273.65 / (M_PI * 0.025)) * probability));
  TMC_EXPECT_NEAR(k, expected, 1e-13 * expected);
}

TMC_TEST(surface_default_probability)
{
  /* MICM's default reaction probability is 1 (SurfaceRateConstantParameters);
   * the C constructor takes the probability, so the default is passed. */
  surface_case(1.0);
  return TMC_TEST_PASS;
}

TMC_TEST(surface_specified_probability)
{
  surface_case(0.74);
  return TMC_TEST_PASS;
}

/* A solver for foo -> (nothing) with a surface rate constant labelled "foo"
 * on foo's gas-phase species. MICM's solver builder throws while building its
 * rate-constant store (ReactionRateConstantStore::BuildFrom) when the species
 * lacks a diffusion coefficient or a molecular weight; tmc_builder_build
 * returns the same error. */
static tmc_solver *build_surface_solver(const tmc_species *foo, bool with_diffusion, double diffusion_coefficient, tmc_error *err)
{
  tmc_phase gas = tmc_phase_create("gas");
  if (with_diffusion)
    tmc_phase_add_species_with_diffusion(&gas, foo, diffusion_coefficient);
  else
    tmc_phase_add_species(&gas, foo);
  tmc_system system = tmc_system_create(&gas);
  tmc_reaction r = tmc_reaction_create();
  tmc_reaction_add_reactant(&r, foo);
  tmc_reaction_set_rate_constant(&r, tmc_rate_surface("foo", &gas.species[0], 1.0));
  tmc_reaction_set_phase(&r, &gas);
  tmc_builder *b = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, &r, 1);
  tmc_solver *s = tmc_builder_build(b, err);
  tmc_builder_free(b);
  tmc_reaction_free(&r);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return s;
}

TMC_TEST(surface_diffusion_coefficient_is_missing)
{
  tmc_species foo = tmc_species_create_with_property("foo", TMC_PROPERTY_MOLECULAR_WEIGHT, 0.025);
  tmc_phase_species foo_gas_species = phase_species(&foo, false, 0.0); /* no diffusion coefficient */
  TMC_EXPECT_FALSE(foo_gas_species.has_diffusion_coefficient);
  tmc_error err;
  tmc_solver *s = build_surface_solver(&foo, false, 0.0, &err);
  TMC_EXPECT_TRUE(s == NULL);
  TMC_EXPECT_EQ(err.category, TMC_ERROR_SPECIES);
  TMC_EXPECT_EQ(err.code, TMC_SPECIES_PROPERTY_NOT_FOUND);
  TMC_EXPECT_TRUE(strstr(err.message, "Diffusion coefficient for species 'foo' is not defined") != NULL);
  tmc_solver_free(s);
  return TMC_TEST_PASS;
}

TMC_TEST(surface_molecular_weight_is_missing)
{
  tmc_species foo = tmc_species_create("foo");
  TMC_EXPECT_FALSE(tmc_species_has_property(&foo, TMC_PROPERTY_MOLECULAR_WEIGHT));
  tmc_error err;
  tmc_solver *s = build_surface_solver(&foo, true, 2.3e2, &err);
  TMC_EXPECT_TRUE(s == NULL);
  TMC_EXPECT_EQ(err.category, TMC_ERROR_SPECIES);
  TMC_EXPECT_EQ(err.code, TMC_SPECIES_PROPERTY_NOT_FOUND);
  /* MICM's message: "Species: 'foo' Property: 'molecular weight [kg mol-1]'" */
  TMC_EXPECT_TRUE(strstr(err.message, "foo") != NULL);
  TMC_EXPECT_TRUE(strstr(err.message, "molecular weight") != NULL);
  tmc_solver_free(s);
  return TMC_TEST_PASS;
}

/* ---- taylor_series_rate_constant ---- */

static void taylor_series_cases(void)
{
  const double tolerance = 1e-13;
  tmc_conditions conditions = { .temperature = 301.24 };

  /* default parameters: A = 1, one coefficient equal to 1, so k = 1 */
  tmc_taylor_series_parameters params = tmc_taylor_series_defaults();
  double k = calculate(tmc_rate_taylor_series(params), conditions, NULL);
  double expected = 1;
  TMC_EXPECT_NEAR(k, expected, tolerance * expected);

  params.A = 1;
  k = calculate(tmc_rate_taylor_series(params), conditions, NULL);
  TMC_EXPECT_NEAR(k, expected, tolerance * expected);

  params.A = 2.2e-10;
  k = calculate(tmc_rate_taylor_series(params), conditions, NULL);
  expected = 2.2e-10;
  TMC_EXPECT_NEAR(k, expected, tolerance * expected);

  /* O + HO2 -> OH + O2 */
  params.A = 3e-11;
  params.C = -200;
  params.coefficients[0] = 12.5;
  params.coefficients[1] = 1.3e-2;
  params.coefficients[2] = 5.2e-4;
  params.n_coefficients = 3;
  k = calculate(tmc_rate_taylor_series(params), conditions, NULL);
  expected = 3e-11 * exp(-200 / 301.24) * (12.5 + 1.3e-2 * 301.24 + 5.2e-4 * pow(301.24, 2));
  TMC_EXPECT_NEAR(k, expected, tolerance * expected);

  /* OH + HCl -> H2O + Cl */
  params.A = 2.6e-12;
  params.C = -350;
  params.coefficients[0] = 1.0;
  params.coefficients[1] = 4.3e-1;
  params.coefficients[2] = 7.3e-3;
  params.n_coefficients = 3;
  k = calculate(tmc_rate_taylor_series(params), conditions, NULL);
  expected = 2.6e-12 * exp(-350 / 301.24) * (1.0 + 4.3e-1 * 301.24 + 7.3e-3 * pow(301.24, 2));
  TMC_EXPECT_NEAR(k, expected, tolerance * expected);
}

TMC_TEST(taylor_series_calculate_with_system)
{
  taylor_series_cases();
  return TMC_TEST_PASS;
}

TMC_TEST(taylor_series_calculate_with_prescribed_arguments)
{
  /* identical to CalculateWithSystem in MICM */
  taylor_series_cases();
  return TMC_TEST_PASS;
}

/* ---- ternary_chemical_activation_rate_constant (exact in double) ---- */

TMC_TEST(ternary_minimal_arguments)
{
  tmc_conditions conditions = { .temperature = 301.24, .air_density = 42.2 };
  tmc_falloff_parameters ternary_params = tmc_falloff_defaults();
  ternary_params.k0_A = 1.0;
  ternary_params.kinf_A = 1.0;
  double k = calculate(tmc_rate_ternary_chemical_activation(ternary_params), conditions, NULL);
  double k0 = 1.0;
  double kinf = 1.0;
  double expected = k0 / (1.0 + k0 * 42.2 / kinf) * pow(0.6, 1.0 / (1 + pow(log10(k0 * 42.2 / kinf), 2)));
  TMC_EXPECT_NEAR(k, expected, 0.0);
  return TMC_TEST_PASS;
}

TMC_TEST(ternary_all_arguments)
{
  double temperature = 301.24;
  tmc_conditions conditions = { .temperature = temperature, .air_density = 42.2 };
  tmc_falloff_parameters params = { .k0_A = 1.2, .k0_B = 2.3, .k0_C = 302.3, .kinf_A = 2.6, .kinf_B = -3.1, .kinf_C = 402.1, .Fc = 0.9, .N = 1.2 };
  double k = calculate(tmc_rate_ternary_chemical_activation(params), conditions, NULL);
  double k0 = 1.2 * exp(302.3 / temperature) * pow(temperature / 300.0, 2.3);
  double kinf = 2.6 * exp(402.1 / temperature) * pow(temperature / 300.0, -3.1);
  double expected = k0 / (1.0 + k0 * 42.2 / kinf) * pow(0.9, 1.0 / (1.0 + 1.0 / 1.2 * pow(log10(k0 * 42.2 / kinf), 2)));
  TMC_EXPECT_NEAR(k, expected, 0.0);
  return TMC_TEST_PASS;
}

TMC_TEST(ternary_jpl19_all_arguments)
{
  double temperature = 301.24;
  double air_density = 42.2;
  tmc_ternary_jpl19_parameters params = { .k0_A = 1.2,   .k0_B = 2.3,  .k0_C = 302.3,  .k0_D = 298.0,  .kinf_A = 2.6,
                                          .kinf_B = -3.1, .kinf_C = 402.1, .kinf_D = 310.0, .kint_A = 0.7, .kint_B = 0.5,
                                          .kint_C = 120.0, .kint_D = 290.0, .Fc = 0.9,      .N = 1.2 };
  tmc_conditions conditions = { .temperature = temperature, .air_density = air_density };
  double k = calculate(tmc_rate_ternary_jpl19(params), conditions, NULL);
  double k0 = 1.2 * exp(302.3 / temperature) * pow(temperature / 298.0, 2.3);
  double kinf = 2.6 * exp(402.1 / temperature) * pow(temperature / 310.0, -3.1);
  double kint = 0.7 * exp(120.0 / temperature) * pow(temperature / 290.0, 0.5);
  double ratio = k0 * air_density / kinf;
  double k_f = k0 * air_density / (1.0 + ratio) * pow(0.9, 1.2 / (1.2 + pow(log10(ratio), 2)));
  double expected = k_f + kint * (1.0 - k_f / kinf);
  TMC_EXPECT_NEAR(k, expected, 0.0);
  return TMC_TEST_PASS;
}

TMC_TEST(ternary_jpl19_zero_kint_reduces_to_troe)
{
  double temperature = 250.0;
  double air_density = 30.0;
  tmc_ternary_jpl19_parameters params = tmc_ternary_jpl19_defaults();
  params.k0_A = 3.0e-3;
  params.k0_B = -2.0;
  params.kinf_A = 1.5;
  params.kinf_B = -1.0;
  params.kint_A = 0.0;
  tmc_conditions conditions = { .temperature = temperature, .air_density = air_density };
  double k = calculate(tmc_rate_ternary_jpl19(params), conditions, NULL);
  /* with k_int = 0 the result is the Troe falloff expression with D = 298 K */
  double k0 = 3.0e-3 * pow(temperature / 298.0, -2.0);
  double kinf = 1.5 * pow(temperature / 298.0, -1.0);
  double ratio = k0 * air_density / kinf;
  double expected = k0 * air_density / (1.0 + ratio) * pow(0.6, 1.0 / (1.0 + pow(log10(ratio), 2)));
  TMC_EXPECT_NEAR(k, expected, 1e-30);
  return TMC_TEST_PASS;
}

/* ---- troe_rate_constant (exact in double, except ExampleBC) ---- */

TMC_TEST(troe_minimal_arguments)
{
  tmc_conditions conditions = { .temperature = 301.24, .air_density = 42.2 };
  tmc_falloff_parameters troe_params = tmc_falloff_defaults();
  troe_params.k0_A = 1.0;
  troe_params.kinf_A = 1.0;
  double k = calculate(tmc_rate_troe(troe_params), conditions, NULL);
  double k0 = 1.0;
  double kinf = 1.0;
  double expected = 42.2 * k0 / (1.0 + 42.2 * k0 / kinf) * pow(0.6, 1.0 / (1 + pow(log10(42.2 * k0 / kinf), 2)));
  TMC_EXPECT_NEAR(k, expected, 0.0);
  return TMC_TEST_PASS;
}

TMC_TEST(troe_all_arguments)
{
  double temperature = 301.24;
  tmc_conditions conditions = { .temperature = temperature, .air_density = 42.2 };
  tmc_falloff_parameters params = { .k0_A = 1.2, .k0_B = 2.3, .k0_C = 302.3, .kinf_A = 2.6, .kinf_B = -3.1, .kinf_C = 402.1, .Fc = 0.9, .N = 1.2 };
  double k = calculate(tmc_rate_troe(params), conditions, NULL);
  double k0 = 1.2 * exp(302.3 / temperature) * pow(temperature / 300.0, 2.3);
  double kinf = 2.6 * exp(402.1 / temperature) * pow(temperature / 300.0, -3.1);
  double expected = 42.2 * k0 / (1.0 + 42.2 * k0 / kinf) * pow(0.9, 1.0 / (1.0 + 1.0 / 1.2 * pow(log10(42.2 * k0 / kinf), 2)));
  TMC_EXPECT_NEAR(k, expected, 0.0);
  return TMC_TEST_PASS;
}

TMC_TEST(troe_analytical_example_ab)
{
  tmc_conditions conditions = { .temperature = 301.24, .air_density = 42.2 };
  tmc_falloff_parameters troe_params = tmc_falloff_defaults();
  troe_params.k0_A = 4.0e-18;
  double k = calculate(tmc_rate_troe(troe_params), conditions, NULL);
  double k_0 = 4.0e-18;
  double k_inf = 1;
  double k1 = k_0 * 42.2 / (1.0 + k_0 * 42.2 / k_inf) * pow(0.6, 1.0 / (1.0 + 1.0 / 1.0 * pow(log10(k_0 * 42.2 / k_inf), 2)));
  TMC_EXPECT_NEAR(k, k1, 0.0); /* MICM: EXPECT_EQ */
  return TMC_TEST_PASS;
}

TMC_TEST(troe_analytical_example_bc)
{
  tmc_conditions conditions = { .temperature = 301.24, .air_density = 42.2 };
  tmc_falloff_parameters troe_params = { .k0_A = 1.2e-12, .k0_B = 167.0, .k0_C = 3.0, .kinf_A = 136.0, .kinf_B = 5.0, .kinf_C = 24.0, .Fc = 0.9, .N = 0.8 };
  double k = calculate(tmc_rate_troe(troe_params), conditions, NULL);
  double k_0 = 1.2e-12 * exp(3.0 / 301.24) * pow(301.24 / 300.0, 167.0);
  double k_inf = 136.0 * exp(24.0 / 301.24) * pow(301.24 / 300.0, 5.0);
  double k1 = k_0 * 42.2 / (1.0 + k_0 * 42.2 / k_inf) * pow(0.9, 1.0 / (1.0 + (1.0 / 0.8) * pow(log10(k_0 * 42.2 / k_inf), 2)));
  /* MICM: relative error at most 1e-14 */
  TMC_EXPECT_REL(k, k1, 1.e-14);
  return TMC_TEST_PASS;
}

/* ---- tunneling_rate_constant ---- */

TMC_TEST(tunneling_minimal_arguments)
{
  tmc_conditions conditions = { .temperature = 301.24 };
  double k = calculate(tmc_rate_tunneling(tmc_tunneling_defaults()), conditions, NULL);
  double expected = 1.0;
  TMC_EXPECT_NEAR(k, expected, 1e-13 * expected);
  return TMC_TEST_PASS;
}

TMC_TEST(tunneling_all_arguments)
{
  double temperature = 301.24;
  tmc_conditions conditions = { .temperature = temperature };
  tmc_tunneling_parameters tunneling_params = { .A = 1.2, .B = 2.3, .C = 302.3 };
  double k = calculate(tmc_rate_tunneling(tunneling_params), conditions, NULL);
  double expected = 1.2 * exp(-2.3 / temperature) * exp(302.3 / pow(temperature, 3));
  TMC_EXPECT_NEAR(k, expected, 1e-13 * expected);
  return TMC_TEST_PASS;
}

/* ---- user_defined_rate_constant (exact in double) ---- */

static double user_defined(double scaling_factor, double custom_param)
{
  tmc_rate_constant rc = tmc_rate_user_defined("label", scaling_factor);
  tmc_conditions conditions = { 0 };
  const double custom_params[1] = { custom_param };
  return calculate(rc, conditions, custom_params);
}

TMC_TEST(user_defined_calculate_with_system)
{
  double k = user_defined(1.0, 0.5);
  TMC_EXPECT_EQ(k, 0.5);
  return TMC_TEST_PASS;
}

TMC_TEST(user_defined_constructor_with_rate)
{
  double k = user_defined(1.0, 1.1);
  TMC_EXPECT_NEAR(k, 1.1, 0.0);
  return TMC_TEST_PASS;
}

TMC_TEST(user_defined_constructor_with_rate_and_name)
{
  tmc_rate_constant rc = tmc_rate_user_defined("a name", 1.0);
  TMC_EXPECT_STREQ(rc.p.user_defined.label, "a name");
  TMC_EXPECT_EQ(tmc_rate_constant_custom_parameter_count(&rc), 1);
  double k = user_defined(1.0, 1.1);
  TMC_EXPECT_NEAR(k, 1.1, 0.0);
  return TMC_TEST_PASS;
}

TMC_TEST(user_defined_custom_scaling_factor)
{
  double k = user_defined(2.0, 1.2);
  TMC_EXPECT_NEAR(k, 1.2 * 2.0, 0.0);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithSystem", arrhenius_calculate_with_system),
  TMC_CASE("arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithPrescribedArguments", arrhenius_calculate_with_prescribed_arguments),
  TMC_CASE("branched_rate_constant:BranchedRateConstant.CalculateAlkoxyBranchWithAllArguments", branched_alkoxy),
  TMC_CASE("branched_rate_constant:BranchedRateConstant.CalculateNitrateBranchWithAllArguments", branched_nitrate),
  TMC_CASE("lambda_rate_constant:LambdaRateConstant.CalculateWithSystem", lambda_calculate_with_system),
  TMC_CASE("reversible_rate_constant:ReversibleRateConstant.DefaultConstructor", reversible_default_constructor),
  TMC_CASE("reversible_rate_constant:ReversibleRateConstant.CalculateWithSystem", reversible_calculate_with_system),
  TMC_CASE("reversible_rate_constant:ReversibleRateConstant.CalculateWithPrescribedArguments", reversible_calculate_with_prescribed_arguments),
  TMC_CASE("reversible_rate_constant:ReversibleRateConstant.TemperatureDependence", reversible_temperature_dependence),
  TMC_CASE("reversible_rate_constant:ReversibleRateConstant.ZeroReverseRate", reversible_zero_reverse_rate),
  TMC_CASE("reversible_rate_constant:ReversibleRateConstant.EquilibriumConstantOnly", reversible_equilibrium_constant_only),
  TMC_CASE("surface_rate_constant:SurfaceRateConstant.CalculateDefaultProbability", surface_default_probability),
  TMC_CASE("surface_rate_constant:SurfaceRateConstant.CalculateSpecifiedProbability", surface_specified_probability),
  TMC_CASE("surface_rate_constant:SurfaceRateConstant.DiffusionCoefficientIsMissing", surface_diffusion_coefficient_is_missing),
  TMC_CASE("surface_rate_constant:SurfaceRateConstant.MolecularWeightIsMissing", surface_molecular_weight_is_missing),
  TMC_CASE("taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithSystem", taylor_series_calculate_with_system),
  TMC_CASE("taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithPrescribedArguments", taylor_series_calculate_with_prescribed_arguments),
  TMC_CASE("ternary_chemical_activation_rate_constant:TernaryChemicalActivationRateConstant.CalculateWithMinimalArguments", ternary_minimal_arguments),
  TMC_CASE("ternary_chemical_activation_rate_constant:TernaryChemicalActivationRateConstant.CalculateWithAllArguments", ternary_all_arguments),
  TMC_CASE("ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.CalculateWithAllArguments", ternary_jpl19_all_arguments),
  TMC_CASE("ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.ZeroKintReducesToTroe", ternary_jpl19_zero_kint_reduces_to_troe),
  TMC_CASE("troe_rate_constant:TroeRateConstant.CalculateWithMinimalArguments", troe_minimal_arguments),
  TMC_CASE("troe_rate_constant:TroeRateConstant.CalculateWithAllArguments", troe_all_arguments),
  TMC_CASE("troe_rate_constant:TroeRateConstant.AnalyticalTroeExampleAB", troe_analytical_example_ab),
  TMC_CASE("troe_rate_constant:TroeRateConstant.AnalyticalTroeExampleBC", troe_analytical_example_bc),
  TMC_CASE("tunneling_rate_constant:TunnelingRateConstant.CalculateWithMinimalArguments", tunneling_minimal_arguments),
  TMC_CASE("tunneling_rate_constant:TunnelingRateConstant.CalculateWithAllArguments", tunneling_all_arguments),
  TMC_CASE("user_defined_rate_constant:UserDefinedRateConstant.CalculateWithSystem", user_defined_calculate_with_system),
  TMC_CASE("user_defined_rate_constant:UserDefinedRateConstant.ConstructorWithRate", user_defined_constructor_with_rate),
  TMC_CASE("user_defined_rate_constant:UserDefinedRateConstant.ConstructorWithRateAndName", user_defined_constructor_with_rate_and_name),
  TMC_CASE("user_defined_rate_constant:UserDefinedRateConstant.CustomScalingFactor", user_defined_custom_scaling_factor),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_rate_constants(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
