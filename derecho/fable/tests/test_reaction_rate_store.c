/* MICM test/unit/process/test_reaction_rate_store.cpp
 *
 * MICM's ReactionRateConstantStore holds each rate-constant type's parameters
 * in one array, in the order of the reactions sorted by type, with an offset
 * per type into the state's rate constants, precomputed fields (branched k0
 * and z, the surface mean-free-speed factor), and custom-parameter indices.
 * The C library has no separate store: the built solver keeps the reactions
 * sorted the same way (tmc_solver_reactions), labels the custom parameters in
 * that order (tmc_solver_custom_parameter_label), and computes the rate
 * constants of a state in those columns (tmc_solver_update_state_parameters),
 * deriving the precomputed fields inside the calculation. Each case checks
 * the store's fields through these: the per-type counts and offsets from the
 * solver's reaction order, the parameters from the solver's reactions, and
 * the indices and derived fields from the rate constants they produce. */
#include "harness.h"
#include "tmc_support.h"

#define AVOGADRO_CONSTANT 6.02214076e23
#define BOLTZMANN_CONSTANT 1.380649e-23
#define GAS_CONSTANT (BOLTZMANN_CONSTANT * AVOGADRO_CONSTANT)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- helpers ---- */

static tmc_reaction reaction(const tmc_species *reactant, const tmc_species *product, tmc_rate_constant rc, const tmc_phase *gas)
{
  tmc_reaction r = tmc_reaction_create();
  tmc_reaction_add_reactant(&r, reactant);
  tmc_reaction_add_product(&r, product, 1.0);
  tmc_reaction_set_rate_constant(&r, rc);
  tmc_reaction_set_phase(&r, gas);
  return r;
}

/* a solver over a gas phase with MICM's dense layout (Matrix: one cell per block) */
static tmc_solver *build(const tmc_phase *gas, const tmc_reaction *reactions, int n, tmc_error *err)
{
  tmc_system system = tmc_system_create(gas);
  tmc_builder *b = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, reactions, n);
  tmc_builder_set_width(b, 1);
  tmc_solver *s = tmc_builder_build(b, err);
  tmc_builder_free(b);
  tmc_system_free(&system);
  return s;
}

static void free_all(tmc_reaction *reactions, int n, tmc_phase *gas, tmc_solver *s)
{
  for (int i = 0; i < n; ++i)
    tmc_reaction_free(&reactions[i]);
  tmc_phase_free(gas);
  tmc_solver_free(s);
}

static int count_type(const tmc_solver *s, tmc_rate_type type)
{
  int n, count = 0;
  const tmc_reaction *r = tmc_solver_reactions(s, &n);
  for (int i = 0; i < n; ++i)
    count += r[i].rate_constant.type == type;
  return count;
}

/* MICM's offsets: the first column of a type follows the columns of every
 * earlier type */
static int type_offset(const tmc_solver *s, tmc_rate_type type)
{
  int n, offset = 0;
  const tmc_reaction *r = tmc_solver_reactions(s, &n);
  for (int i = 0; i < n; ++i)
    offset += r[i].rate_constant.type < type;
  return offset;
}

/* each type occupies one contiguous block of columns */
static bool types_contiguous(const tmc_solver *s)
{
  int n;
  const tmc_reaction *r = tmc_solver_reactions(s, &n);
  for (int i = 1; i < n; ++i)
    if (r[i].rate_constant.type < r[i - 1].rate_constant.type)
      return false;
  return true;
}

/* MICM's CalculateSurfaceOne with the mean-free-speed factor BuildFrom precomputes */
static double surface_rate(double mw, double diffusion_coefficient, double probability, double temperature, double radius, double number_concentration)
{
  double mean_free_speed_factor = 8.0 * GAS_CONSTANT / (M_PI * mw);
  double mean_free_speed = sqrt(mean_free_speed_factor * temperature);
  return 4.0 * number_concentration * M_PI * radius * radius / (radius / diffusion_coefficient + 4.0 / (mean_free_speed * probability));
}

/* ---- Offset arithmetic ---- */

TMC_TEST(offsets_are_contiguous_cumulative_sizes)
{
  tmc_species a = tmc_species_create("a"), b = tmc_species_create("b");
  tmc_species c = tmc_species_create_with_property("c", TMC_PROPERTY_MOLECULAR_WEIGHT, 0.025);
  double c_diff = 1.0e-5;
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &a);
  tmc_phase_add_species(&gas, &b);
  tmc_phase_add_species_with_diffusion(&gas, &c, c_diff);
  const tmc_phase_species *gas_c = &gas.species[2];

  tmc_falloff_parameters troe = tmc_falloff_defaults();
  troe.k0_A = 1.0;
  troe.kinf_A = 1.0;
  /* MICM's test sorts the processes by type (SortByTypeOrder) before
   * BuildFrom; the solver builder sorts them itself, so they are given here
   * out of order */
  tmc_reaction procs[6] = {
    reaction(&a, &b, tmc_rate_user_defined("ud1", 1.0), &gas),
    reaction(&a, &b, tmc_rate_arrhenius(tmc_arrhenius_defaults()), &gas),
    reaction(&c, &b, tmc_rate_surface("surf1", gas_c, 1.0), &gas),
    reaction(&a, &b, tmc_rate_tunneling(tmc_tunneling_defaults()), &gas),
    reaction(&a, &b, tmc_rate_arrhenius(tmc_arrhenius_defaults()), &gas),
    reaction(&a, &b, tmc_rate_troe(troe), &gas),
  };
  tmc_error err;
  tmc_solver *s = build(&gas, procs, 6, &err);
  TMC_ASSERT_TRUE(s != NULL);

  /* Arrhenius 2, Troe 1, ternary 0, branched 0, tunneling 1, Taylor 0,
   * reversible 0, user-defined 1, surface 1 */
  TMC_EXPECT_EQ(count_type(s, TMC_RATE_ARRHENIUS), 2);
  TMC_EXPECT_EQ(count_type(s, TMC_RATE_TROE), 1);
  TMC_EXPECT_EQ(count_type(s, TMC_RATE_TERNARY_CHEMICAL_ACTIVATION), 0);
  TMC_EXPECT_EQ(count_type(s, TMC_RATE_BRANCHED), 0);
  TMC_EXPECT_EQ(count_type(s, TMC_RATE_TUNNELING), 1);
  TMC_EXPECT_EQ(count_type(s, TMC_RATE_TAYLOR_SERIES), 0);
  TMC_EXPECT_EQ(count_type(s, TMC_RATE_REVERSIBLE), 0);
  TMC_EXPECT_EQ(count_type(s, TMC_RATE_USER_DEFINED), 1);
  TMC_EXPECT_EQ(count_type(s, TMC_RATE_SURFACE), 1);
  TMC_EXPECT_TRUE(types_contiguous(s));

  /* cumulative offsets */
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_TROE), 2);
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_TERNARY_CHEMICAL_ACTIVATION), 3);
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_TERNARY_JPL19), 3);
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_BRANCHED), 3);
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_TUNNELING), 3);
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_TAYLOR_SERIES), 4);
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_REVERSIBLE), 4);
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_USER_DEFINED), 4);
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_SURFACE), 5);
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_LAMBDA), 6);

  /* the total (lambda offset plus lambda count) is the number of reactions,
   * one rate-constant column each */
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_LAMBDA) + count_type(s, TMC_RATE_LAMBDA), 6);
  TMC_EXPECT_EQ(tmc_solver_num_reactions(s), 6);

  /* and the state's rate constants follow the offsets */
  tmc_state *state = tmc_solver_create_state(s, 1);
  TMC_EXPECT_EQ(tmc_state_rate_constants(state)->n_cols, 6);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0] = (tmc_conditions){ .temperature = 287.4, .pressure = 101325.0, .air_density = 42.2 };
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "ud1", 0.5, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "surf1.effective radius [m]", 1.0e-7, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "surf1.particle number concentration [# m-3]", 2.5e6, &err), TMC_OK);
  tmc_solver_update_state_parameters(s, state);
  const tmc_dense *k = tmc_state_rate_constants(state);
  double troe_value = tmc_rate_constant_calculate(&procs[5].rate_constant, &conditions[0], NULL);
  TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(k, 0, 0), 1.0);
  TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(k, 0, 1), 1.0);
  TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(k, 0, type_offset(s, TMC_RATE_TROE)), troe_value);
  TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(k, 0, type_offset(s, TMC_RATE_TUNNELING)), 1.0);
  TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(k, 0, type_offset(s, TMC_RATE_USER_DEFINED)), 0.5);
  TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(k, 0, type_offset(s, TMC_RATE_SURFACE)), surface_rate(0.025, c_diff, 1.0, 287.4, 1.0e-7, 2.5e6));
  TMC_EXPECT_TRUE(troe_value != 1.0);
  tmc_state_free(state);
  free_all(procs, 6, &gas, s);
  return TMC_TEST_PASS;
}

/* ---- Single Arrhenius: parameters preserved ---- */

TMC_TEST(arrhenius_parameters_preserved)
{
  tmc_species a = tmc_species_create("a"), b = tmc_species_create("b");
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &a);
  tmc_phase_add_species(&gas, &b);
  tmc_arrhenius_parameters params = { .A = 2.15e-4, .B = 1.2, .C = 110.0, .D = 300.0, .E = 0.0 };
  tmc_reaction procs[1] = { reaction(&a, &b, tmc_rate_arrhenius(params), &gas) };
  tmc_error err;
  tmc_solver *s = build(&gas, procs, 1, &err);
  TMC_ASSERT_TRUE(s != NULL);
  TMC_ASSERT_TRUE(count_type(s, TMC_RATE_ARRHENIUS) == 1);
  const tmc_arrhenius_parameters *stored = &tmc_solver_reactions(s, NULL)[0].rate_constant.p.arrhenius;
  TMC_EXPECT_DOUBLE_EQ(stored->A, 2.15e-4);
  TMC_EXPECT_DOUBLE_EQ(stored->B, 1.2);
  TMC_EXPECT_DOUBLE_EQ(stored->C, 110.0);
  TMC_EXPECT_DOUBLE_EQ(stored->D, 300.0);

  /* and the state's rate constant is computed from them */
  tmc_state *state = tmc_solver_create_state(s, 1);
  tmc_state_conditions(state)[0] = (tmc_conditions){ .temperature = 287.4, .pressure = 101325.0 };
  tmc_solver_update_state_parameters(s, state);
  double expected = 2.15e-4 * exp(110.0 / 287.4) * pow(287.4 / 300.0, 1.2) * (1.0 + 0.0 * 101325.0);
  TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(tmc_state_rate_constants(state), 0, 0), expected);
  tmc_state_free(state);
  free_all(procs, 1, &gas, s);
  return TMC_TEST_PASS;
}

/* ---- Branched: derived fields k0 and z ---- */

TMC_TEST(branched_derived_fields_computed)
{
  tmc_species a = tmc_species_create("a"), b = tmc_species_create("b");
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &a);
  tmc_phase_add_species(&gas, &b);
  tmc_branched_parameters params = { .branch = TMC_BRANCH_ALKOXY, .X = 1.0, .Y = 0.0, .a0 = 0.5, .n = 3 };
  tmc_reaction procs[1] = { reaction(&a, &b, tmc_rate_branched(params), &gas) };
  tmc_error err;
  tmc_solver *s = build(&gas, procs, 1, &err);
  TMC_ASSERT_TRUE(s != NULL);
  TMC_ASSERT_TRUE(count_type(s, TMC_RATE_BRANCHED) == 1);

  /* k0 = 2e-22 N_A 1e-6 exp(n); z = A(293 K, [M]_ref) (1 - a0) / a0 */
  double expected_k0 = 2.0e-22 * AVOGADRO_CONSTANT * 1.0e-6 * exp(3.0);
  double air_ref = 2.45e19 / AVOGADRO_CONSTANT * 1.0e6;
  double a_val = expected_k0 * air_ref;
  double b_val = 0.43 * pow(293.0 / 298.0, -8.0);
  double A_val = a_val / (1.0 + a_val / b_val) * pow(0.41, 1.0 / (1.0 + pow(log10(a_val / b_val), 2.0)));
  double expected_z = A_val * (1.0 - 0.5) / 0.5;

  /* the solver derives k0 and z inside the calculation; the alkoxy rate
   * X exp(-Y / T) z / (z + A(T, M)) with these k0 and z must match it at
   * several conditions (at the reference conditions it is 1 - a0) */
  const tmc_conditions conditions[3] = { { .temperature = 293.0, .air_density = air_ref },
                                         { .temperature = 301.24, .air_density = 42.2 },
                                         { .temperature = 250.0, .air_density = 15.0 } };
  tmc_state *state = tmc_solver_create_state(s, 3);
  for (int c = 0; c < 3; ++c)
    tmc_state_conditions(state)[c] = conditions[c];
  tmc_solver_update_state_parameters(s, state);
  for (int c = 0; c < 3; ++c)
  {
    double T = conditions[c].temperature, M = conditions[c].air_density;
    double a = expected_k0 * M;
    double bb = 0.43 * pow(T / 298.0, -8.0);
    double A = a / (1.0 + a / bb) * pow(0.41, 1.0 / (1.0 + pow(log10(a / bb), 2.0)));
    double expected = 1.0 * exp(-0.0 / T) * (expected_z / (expected_z + A));
    TMC_EXPECT_NEAR(tmc_dense_get(tmc_state_rate_constants(state), c, 0), expected, 1.0e-10 * fabs(expected));
  }
  TMC_EXPECT_NEAR(tmc_dense_get(tmc_state_rate_constants(state), 0, 0), 0.5, 1.0e-10);
  tmc_state_free(state);
  free_all(procs, 1, &gas, s);
  return TMC_TEST_PASS;
}

/* ---- UserDefined: custom parameter indices ---- */

TMC_TEST(user_defined_custom_param_index)
{
  tmc_species a = tmc_species_create("a"), b = tmc_species_create("b");
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &a);
  tmc_phase_add_species(&gas, &b);
  tmc_reaction procs[2] = { reaction(&a, &b, tmc_rate_user_defined("rate1", 2.0), &gas),
                            reaction(&a, &b, tmc_rate_user_defined("rate2", 0.5), &gas) };
  tmc_error err;
  tmc_solver *s = build(&gas, procs, 2, &err);
  TMC_ASSERT_TRUE(s != NULL);
  TMC_ASSERT_TRUE(count_type(s, TMC_RATE_USER_DEFINED) == 2);
  const tmc_reaction *sorted = tmc_solver_reactions(s, NULL);

  /* custom parameter 0 is rate1's, 1 is rate2's; the scaling factors are kept */
  TMC_ASSERT_TRUE(tmc_solver_num_custom_parameters(s) == 2);
  TMC_EXPECT_STREQ(tmc_solver_custom_parameter_label(s, 0), "rate1");
  TMC_EXPECT_STREQ(tmc_solver_custom_parameter_label(s, 1), "rate2");
  TMC_EXPECT_STREQ(sorted[0].rate_constant.p.user_defined.label, "rate1");
  TMC_EXPECT_DOUBLE_EQ(sorted[0].rate_constant.p.user_defined.scaling_factor, 2.0);
  TMC_EXPECT_STREQ(sorted[1].rate_constant.p.user_defined.label, "rate2");
  TMC_EXPECT_DOUBLE_EQ(sorted[1].rate_constant.p.user_defined.scaling_factor, 0.5);

  /* each rate constant reads its own parameter */
  tmc_state *state = tmc_solver_create_state(s, 1);
  TMC_EXPECT_EQ(tmc_state_custom_parameter_index(state, "rate1"), 0);
  TMC_EXPECT_EQ(tmc_state_custom_parameter_index(state, "rate2"), 1);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "rate1", 3.0, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "rate2", 7.0, &err), TMC_OK);
  tmc_solver_update_state_parameters(s, state);
  TMC_EXPECT_EQ(tmc_dense_get(tmc_state_rate_constants(state), 0, 0), 3.0 * 2.0);
  TMC_EXPECT_EQ(tmc_dense_get(tmc_state_rate_constants(state), 0, 1), 7.0 * 0.5);
  tmc_state_free(state);
  free_all(procs, 2, &gas, s);
  return TMC_TEST_PASS;
}

/* ---- Surface: data fields and custom parameter index ---- */

TMC_TEST(surface_data_fields_and_custom_param_index)
{
  double mw = 0.025, diff_coeff = 2.3e2, prob = 0.74;
  tmc_species c = tmc_species_create_with_property("c", TMC_PROPERTY_MOLECULAR_WEIGHT, mw);
  tmc_species b = tmc_species_create("b");
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species_with_diffusion(&gas, &c, diff_coeff);
  tmc_phase_add_species(&gas, &b);
  tmc_reaction procs[1] = { reaction(&c, &b, tmc_rate_surface("surf", &gas.species[0], prob), &gas) };
  tmc_error err;
  tmc_solver *s = build(&gas, procs, 1, &err);
  TMC_ASSERT_TRUE(s != NULL);
  TMC_ASSERT_TRUE(count_type(s, TMC_RATE_SURFACE) == 1);

  /* the fields BuildFrom stores: diffusion coefficient, reaction probability,
   * and (from the molecular weight) the mean-free-speed factor */
  const tmc_surface_parameters *stored = &tmc_solver_reactions(s, NULL)[0].rate_constant.p.surface;
  TMC_EXPECT_TRUE(stored->species.has_diffusion_coefficient);
  TMC_EXPECT_DOUBLE_EQ(stored->species.diffusion_coefficient, diff_coeff);
  TMC_EXPECT_DOUBLE_EQ(stored->reaction_probability, prob);
  double stored_mw = 0.0;
  TMC_EXPECT_EQ(tmc_species_get_double(&stored->species.species, TMC_PROPERTY_MOLECULAR_WEIGHT, &stored_mw, &err), TMC_OK);
  TMC_EXPECT_DOUBLE_EQ(stored_mw, mw);

  /* custom parameter base index 0: radius at 0, number concentration at 1 */
  TMC_ASSERT_TRUE(tmc_solver_num_custom_parameters(s) == 2);
  TMC_EXPECT_STREQ(tmc_solver_custom_parameter_label(s, 0), "surf.effective radius [m]");
  TMC_EXPECT_STREQ(tmc_solver_custom_parameter_label(s, 1), "surf.particle number concentration [# m-3]");

  /* the rate constant uses the factor 8 R / (pi Mw) to within 4 ULP */
  tmc_state *state = tmc_solver_create_state(s, 1);
  tmc_state_conditions(state)[0] = (tmc_conditions){ .temperature = 273.65 };
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "surf.effective radius [m]", 1.0e-7, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "surf.particle number concentration [# m-3]", 2.5e6, &err), TMC_OK);
  tmc_solver_update_state_parameters(s, state);
  TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(tmc_state_rate_constants(state), 0, 0), surface_rate(mw, diff_coeff, prob, 273.65, 1.0e-7, 2.5e6));
  tmc_state_free(state);
  free_all(procs, 1, &gas, s);
  return TMC_TEST_PASS;
}

TMC_TEST(surface_custom_param_index_after_user_defined)
{
  double mw = 0.025, diff_coeff = 1.0e-5;
  tmc_species a = tmc_species_create("a"), b = tmc_species_create("b");
  tmc_species c = tmc_species_create_with_property("c", TMC_PROPERTY_MOLECULAR_WEIGHT, mw);
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &a);
  tmc_phase_add_species(&gas, &b);
  tmc_phase_add_species_with_diffusion(&gas, &c, diff_coeff);
  tmc_reaction procs[2] = { reaction(&a, &b, tmc_rate_user_defined("ud", 1.0), &gas),
                            reaction(&c, &b, tmc_rate_surface("surf", &gas.species[2], 1.0), &gas) };
  tmc_error err;
  tmc_solver *s = build(&gas, procs, 2, &err);
  TMC_ASSERT_TRUE(s != NULL);
  TMC_ASSERT_TRUE(count_type(s, TMC_RATE_USER_DEFINED) == 1);
  TMC_ASSERT_TRUE(count_type(s, TMC_RATE_SURFACE) == 1);

  /* the user-defined parameter takes slot 0; the surface takes 1 (radius) and 2 */
  TMC_ASSERT_TRUE(tmc_solver_num_custom_parameters(s) == 3);
  TMC_EXPECT_STREQ(tmc_solver_custom_parameter_label(s, 0), "ud");
  TMC_EXPECT_STREQ(tmc_solver_custom_parameter_label(s, 1), "surf.effective radius [m]");
  TMC_EXPECT_STREQ(tmc_solver_custom_parameter_label(s, 2), "surf.particle number concentration [# m-3]");

  tmc_state *state = tmc_solver_create_state(s, 1);
  tmc_state_conditions(state)[0] = (tmc_conditions){ .temperature = 273.65 };
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "ud", 2.0, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "surf.effective radius [m]", 3.0e-7, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "surf.particle number concentration [# m-3]", 4.0e5, &err), TMC_OK);
  tmc_solver_update_state_parameters(s, state);
  const tmc_dense *k = tmc_state_rate_constants(state);
  TMC_EXPECT_EQ(tmc_dense_get(k, 0, type_offset(s, TMC_RATE_USER_DEFINED)), 2.0);
  TMC_EXPECT_DOUBLE_EQ(tmc_dense_get(k, 0, type_offset(s, TMC_RATE_SURFACE)), surface_rate(mw, diff_coeff, 1.0, 273.65, 3.0e-7, 4.0e5));
  tmc_state_free(state);
  free_all(procs, 2, &gas, s);
  return TMC_TEST_PASS;
}

/* ---- Lambda: entries and their rate-constant index ---- */

static double milli_temperature(const tmc_conditions *conditions, void *context)
{
  (void)context;
  return conditions->temperature * 1.0e-3;
}

TMC_TEST(lambda_entries_rc_index)
{
  tmc_species a = tmc_species_create("a"), b = tmc_species_create("b");
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &a);
  tmc_phase_add_species(&gas, &b);
  tmc_reaction procs[2] = { reaction(&a, &b, tmc_rate_arrhenius(tmc_arrhenius_defaults()), &gas),
                            reaction(&a, &b, tmc_rate_lambda("lam", milli_temperature, NULL), &gas) };
  tmc_error err;
  tmc_solver *s = build(&gas, procs, 2, &err);
  TMC_ASSERT_TRUE(s != NULL);
  TMC_EXPECT_EQ(count_type(s, TMC_RATE_ARRHENIUS), 1);
  TMC_ASSERT_TRUE(count_type(s, TMC_RATE_LAMBDA) == 1);

  /* the lambda is the second process: rate-constant index 1 */
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_LAMBDA), 1);
  TMC_EXPECT_EQ(tmc_solver_reactions(s, NULL)[1].rate_constant.type, TMC_RATE_LAMBDA);
  /* its source is the solver's copy of the parameters, and the function works */
  tmc_lambda_parameters *source = tmc_support_solver_lambda_rate_constant(s, "lam", &err);
  TMC_ASSERT_TRUE(source != NULL && source->function != NULL);
  tmc_conditions cond = { .temperature = 300.0 };
  TMC_EXPECT_NEAR(source->function(&cond, source->context), 0.3, 1.0e-14);
  /* and the state's rate constant 1 is the lambda's value */
  tmc_state *state = tmc_solver_create_state(s, 1);
  tmc_state_conditions(state)[0] = cond;
  tmc_solver_update_state_parameters(s, state);
  TMC_EXPECT_NEAR(tmc_dense_get(tmc_state_rate_constants(state), 0, 1), 0.3, 1.0e-14);
  TMC_EXPECT_EQ(tmc_dense_get(tmc_state_rate_constants(state), 0, 0), 1.0);
  tmc_state_free(state);
  free_all(procs, 2, &gas, s);
  return TMC_TEST_PASS;
}

/* ---- Parameterized multipliers ---- */

TMC_TEST(parameterized_multipliers)
{
  tmc_species a = tmc_species_create("a"), b = tmc_species_create("b");
  tmc_species_set_parameterization(&b, 0.0, 0.0, 0.0, 2.0); /* MICM: { .c_rho_ = 2.0 } */
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &a);
  tmc_phase_add_species(&gas, &b);
  tmc_reaction procs[2] = { reaction(&a, &b, tmc_rate_arrhenius(tmc_arrhenius_defaults()), &gas),   /* no parameterized reactant */
                            reaction(&b, &a, tmc_rate_arrhenius(tmc_arrhenius_defaults()), &gas) }; /* b is parameterized */
  tmc_error err;
  tmc_solver *s = build(&gas, procs, 2, &err);
  TMC_ASSERT_TRUE(s != NULL);
  TMC_EXPECT_EQ(tmc_solver_num_variables(s), 1); /* b holds no state variable */

  /* only the second reaction (rate-constant index 1) is multiplied, by
   * 2.0 * air density = 10 at an air density of 5; Arrhenius() is 1 */
  tmc_state *state = tmc_solver_create_state(s, 1);
  tmc_state_conditions(state)[0] = (tmc_conditions){ .temperature = 300.0, .air_density = 5.0 };
  tmc_solver_update_state_parameters(s, state);
  TMC_EXPECT_NEAR(tmc_dense_get(tmc_state_rate_constants(state), 0, 0), 1.0, 1.0e-14);
  TMC_EXPECT_NEAR(tmc_dense_get(tmc_state_rate_constants(state), 0, 1), 10.0, 1.0e-14);
  tmc_state_free(state);
  free_all(procs, 2, &gas, s);
  return TMC_TEST_PASS;
}

/* ---- Total rate constants for a mixed mechanism ---- */

static double one(const tmc_conditions *conditions, void *context)
{
  (void)conditions;
  (void)context;
  return 1.0;
}

TMC_TEST(total_rate_constants_matches_process_count)
{
  tmc_species a = tmc_species_create("a"), b = tmc_species_create("b");
  tmc_species c = tmc_species_create_with_property("c", TMC_PROPERTY_MOLECULAR_WEIGHT, 0.025);
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &a);
  tmc_phase_add_species(&gas, &b);
  tmc_phase_add_species_with_diffusion(&gas, &c, 1.0e-5);
  tmc_falloff_parameters troe = tmc_falloff_defaults();
  troe.k0_A = 1.0;
  troe.kinf_A = 1.0;
  tmc_reaction procs[5] = {
    reaction(&a, &b, tmc_rate_arrhenius(tmc_arrhenius_defaults()), &gas),
    reaction(&a, &b, tmc_rate_troe(troe), &gas),
    reaction(&a, &b, tmc_rate_user_defined("ud", 1.0), &gas),
    reaction(&c, &b, tmc_rate_surface("surf", &gas.species[2], 1.0), &gas),
    reaction(&a, &b, tmc_rate_lambda("lam", one, NULL), &gas),
  };
  tmc_error err;
  tmc_solver *s = build(&gas, procs, 5, &err);
  TMC_ASSERT_TRUE(s != NULL);
  TMC_EXPECT_EQ(type_offset(s, TMC_RATE_LAMBDA) + count_type(s, TMC_RATE_LAMBDA), 5);
  TMC_EXPECT_EQ(tmc_solver_num_reactions(s), 5);
  tmc_state *state = tmc_solver_create_state(s, 1);
  TMC_EXPECT_EQ(tmc_state_rate_constants(state)->n_cols, 5);
  tmc_state_free(state);
  free_all(procs, 5, &gas, s);
  return TMC_TEST_PASS;
}

/* ---- Errors: a surface species without a diffusion coefficient or molecular weight ---- */

TMC_TEST(surface_missing_diffusion_coefficient_throws)
{
  tmc_species c = tmc_species_create_with_property("c", TMC_PROPERTY_MOLECULAR_WEIGHT, 0.025);
  tmc_species b = tmc_species_create("b");
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &c); /* no diffusion coefficient */
  tmc_phase_add_species(&gas, &b);
  TMC_EXPECT_FALSE(gas.species[0].has_diffusion_coefficient);
  tmc_reaction procs[1] = { reaction(&c, &b, tmc_rate_surface("surf", &gas.species[0], 1.0), &gas) };
  tmc_error err;
  tmc_solver *s = build(&gas, procs, 1, &err);
  TMC_EXPECT_TRUE(s == NULL);
  TMC_EXPECT_EQ(err.category, TMC_ERROR_SPECIES);
  TMC_EXPECT_EQ(err.code, TMC_SPECIES_PROPERTY_NOT_FOUND);
  free_all(procs, 1, &gas, s);
  return TMC_TEST_PASS;
}

TMC_TEST(surface_missing_molecular_weight_throws)
{
  tmc_species c = tmc_species_create("c"); /* no properties at all */
  tmc_species b = tmc_species_create("b");
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species_with_diffusion(&gas, &c, 1.0e-5);
  tmc_phase_add_species(&gas, &b);
  tmc_reaction procs[1] = { reaction(&c, &b, tmc_rate_surface("surf", &gas.species[0], 1.0), &gas) };
  tmc_error err;
  tmc_solver *s = build(&gas, procs, 1, &err);
  TMC_EXPECT_TRUE(s == NULL);
  TMC_EXPECT_EQ(err.category, TMC_ERROR_SPECIES);
  TMC_EXPECT_EQ(err.code, TMC_SPECIES_PROPERTY_NOT_FOUND);
  free_all(procs, 1, &gas, s);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("reaction_rate_store:ReactionRateConstantStore.OffsetsAreContiguousCumulativeSizes", offsets_are_contiguous_cumulative_sizes),
  TMC_CASE("reaction_rate_store:ReactionRateConstantStore.ArrheniusParametersPreserved", arrhenius_parameters_preserved),
  TMC_CASE("reaction_rate_store:ReactionRateConstantStore.BranchedDerivedFieldsComputed", branched_derived_fields_computed),
  TMC_CASE("reaction_rate_store:ReactionRateConstantStore.UserDefinedCustomParamIndex", user_defined_custom_param_index),
  TMC_CASE("reaction_rate_store:ReactionRateConstantStore.SurfaceDataFieldsAndCustomParamIndex", surface_data_fields_and_custom_param_index),
  TMC_CASE("reaction_rate_store:ReactionRateConstantStore.SurfaceCustomParamIndexAfterUserDefined", surface_custom_param_index_after_user_defined),
  TMC_CASE("reaction_rate_store:ReactionRateConstantStore.LambdaEntriesRcIndex", lambda_entries_rc_index),
  TMC_CASE("reaction_rate_store:ReactionRateConstantStore.ParameterizedMultipliers", parameterized_multipliers),
  TMC_CASE("reaction_rate_store:ReactionRateConstantStore.TotalRateConstantsMatchesProcessCount", total_rate_constants_matches_process_count),
  TMC_CASE("reaction_rate_store:ReactionRateConstantStore.SurfaceMissingDiffusionCoefficientThrows", surface_missing_diffusion_coefficient_throws),
  TMC_CASE("reaction_rate_store:ReactionRateConstantStore.SurfaceMissingMolecularWeightThrows", surface_missing_molecular_weight_throws),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_reaction_rate_store(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
