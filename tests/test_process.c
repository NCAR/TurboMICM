/* MICM test/unit/process/test_process.cpp
 *
 * MICM builds a ReactionRateConstantStore and a State by hand and checks the
 * rate constants CalculateRateConstants writes. Here the store and state come
 * from a built solver: tmc_solver_update_state_parameters writes the state's
 * rate constants, one column per reaction in the solver's order (stably
 * sorted by rate-constant type, as MICM's SolverBuilder sorts them). */
#include "harness.h"

#define AVOGADRO_CONSTANT 6.02214076e23
#define BOLTZMANN_CONSTANT 1.380649e-23
#define GAS_CONSTANT (BOLTZMANN_CONSTANT * AVOGADRO_CONSTANT)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Lognormal values from a fixed-seed generator (MICM: std::lognormal_distribution
 * with std::default_random_engine); the checks compare against values computed
 * from the same draws. */
typedef struct
{
  unsigned long long state;
} rng;

static double uniform(rng *g)
{
  g->state = g->state * 6364136223846793005ULL + 1442695040888963407ULL;
  return ((double)(g->state >> 11) + 0.5) / 9007199254740992.0;
}

static double lognormal(rng *g, double mean, double sigma)
{
  double u1 = uniform(g), u2 = uniform(g);
  return exp(mean + sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2));
}

/* MICM's store offsets: the first rate constant of a type follows those of
 * the earlier types. */
static int type_offset(const tmc_solver *s, tmc_rate_type type)
{
  int n;
  const tmc_reaction *r = tmc_solver_reactions(s, &n);
  int offset = 0;
  for (int i = 0; i < n; ++i)
    if (r[i].rate_constant.type < type)
      ++offset;
  return offset;
}

static void process_update_state(int L, int number_of_grid_cells)
{
  tmc_species foo = tmc_species_create_with_property("foo", TMC_PROPERTY_MOLECULAR_WEIGHT, 0.025);
  tmc_species bar = tmc_species_create("bar");
  tmc_species_set_parameterization(&bar, 0.0, 0.0, 0.0, 0.82); /* MICM: { .c_rho_ = 0.82 } */

  double foo_diff_coeff = 2.3e2;
  tmc_phase gas_phase = tmc_phase_create("gas");
  tmc_phase_add_species_with_diffusion(&gas_phase, &foo, foo_diff_coeff);
  tmc_phase_add_species(&gas_phase, &bar);
  const tmc_phase_species *gas_foo = &gas_phase.species[0];

  tmc_arrhenius_parameters rc1_params = tmc_arrhenius_defaults();
  rc1_params.A = 12.2;
  rc1_params.C = 300.0;

  tmc_reaction r[3] = { tmc_reaction_create(), tmc_reaction_create(), tmc_reaction_create() };
  tmc_reaction_add_reactant(&r[0], &foo);
  tmc_reaction_add_reactant(&r[0], &bar);
  tmc_reaction_set_rate_constant(&r[0], tmc_rate_arrhenius(rc1_params));
  tmc_reaction_add_reactant(&r[1], &foo);
  tmc_reaction_set_rate_constant(&r[1], tmc_rate_surface("foo_surf", gas_foo, 1.0));
  tmc_reaction_add_reactant(&r[2], &bar);
  tmc_reaction_set_rate_constant(&r[2], tmc_rate_user_defined("bar_user", 1.0));
  for (int i = 0; i < 3; ++i)
    tmc_reaction_set_phase(&r[i], &gas_phase);

  tmc_system system = tmc_system_create(&gas_phase);
  tmc_builder *b = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(b, &system);
  tmc_builder_set_reactions(b, r, 3);
  tmc_builder_set_width(b, L);
  tmc_error err;
  tmc_solver *s = tmc_builder_build(b, &err);
  TMC_EXPECT_TRUE(s != NULL);
  if (s)
  {
    /* sorted as MICM's SolverBuilder sorts: Arrhenius, user-defined, surface */
    int n_rxn;
    const tmc_reaction *sorted = tmc_solver_reactions(s, &n_rxn);
    TMC_EXPECT_EQ(n_rxn, 3);
    const int arr = type_offset(s, TMC_RATE_ARRHENIUS);
    const int ud = type_offset(s, TMC_RATE_USER_DEFINED);
    const int surf = type_offset(s, TMC_RATE_SURFACE);
    TMC_EXPECT_EQ(arr, 0);
    TMC_EXPECT_EQ(ud, 1);
    TMC_EXPECT_EQ(surf, 2);
    TMC_EXPECT_EQ(sorted[arr].rate_constant.type, TMC_RATE_ARRHENIUS);
    TMC_EXPECT_EQ(sorted[ud].rate_constant.type, TMC_RATE_USER_DEFINED);
    TMC_EXPECT_EQ(sorted[surf].rate_constant.type, TMC_RATE_SURFACE);
    /* custom parameter labels in sorted order */
    TMC_EXPECT_EQ(tmc_solver_num_custom_parameters(s), 3);
    if (tmc_solver_num_custom_parameters(s) == 3)
    {
      TMC_EXPECT_STREQ(tmc_solver_custom_parameter_label(s, 0), "bar_user");
      TMC_EXPECT_STREQ(tmc_solver_custom_parameter_label(s, 1), "foo_surf.effective radius [m]");
      TMC_EXPECT_STREQ(tmc_solver_custom_parameter_label(s, 2), "foo_surf.particle number concentration [# m-3]");
    }

    tmc_state *state = tmc_solver_create_state(s, number_of_grid_cells);
    TMC_EXPECT_EQ(tmc_state_width(state), L);
    rng g = { 1 };
    double *user_rate = malloc(sizeof(double) * (size_t)number_of_grid_cells);
    double *radius = malloc(sizeof(double) * (size_t)number_of_grid_cells);
    double *num_conc = malloc(sizeof(double) * (size_t)number_of_grid_cells);
    tmc_conditions *conditions = tmc_state_conditions(state);
    for (int i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
    {
      conditions[i_cell].temperature = lognormal(&g, 0.0, 0.01) * 285.0;
      conditions[i_cell].pressure = lognormal(&g, 0.0, 0.01) * 101100.0;
      conditions[i_cell].air_density = lognormal(&g, 0.0, 0.01) * 10.0;
      user_rate[i_cell] = lognormal(&g, 0.0, 0.01) * 1.0e-2;
      radius[i_cell] = lognormal(&g, 0.0, 0.01) * 1.0e-8;
      num_conc[i_cell] = lognormal(&g, 0.0, 0.01) * 1.0e5;
    }
    TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameters(state, "bar_user", user_rate, number_of_grid_cells, &err), TMC_OK);
    TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameters(state, "foo_surf.effective radius [m]", radius, number_of_grid_cells, &err), TMC_OK);
    TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameters(state, "foo_surf.particle number concentration [# m-3]", num_conc,
                                                       number_of_grid_cells, &err),
                  TMC_OK);

    tmc_solver_update_state_parameters(s, state);

    const tmc_dense *rate_constants = tmc_state_rate_constants(state);
    const tmc_dense *params = tmc_state_custom_rate_parameters(state);
    for (int i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
    {
      const tmc_conditions *cond = &conditions[i_cell];
      double ur = tmc_dense_get(params, i_cell, tmc_state_custom_parameter_index(state, "bar_user"));
      double rad = tmc_dense_get(params, i_cell, tmc_state_custom_parameter_index(state, "foo_surf.effective radius [m]"));
      double nc = tmc_dense_get(params, i_cell, tmc_state_custom_parameter_index(state, "foo_surf.particle number concentration [# m-3]"));
      TMC_EXPECT_EQ(ur, user_rate[i_cell]);
      TMC_EXPECT_EQ(rad, radius[i_cell]);
      TMC_EXPECT_EQ(nc, num_conc[i_cell]);

      /* r1 (Arrhenius); bar is parameterized: air density * 0.82 */
      double expected_arr = 12.2 * exp(300.0 / cond->temperature) * pow(cond->temperature / 300.0, 0.0) * (1.0 + 0.0 * cond->pressure);
      expected_arr *= (cond->air_density * 0.82);
      TMC_EXPECT_NEAR(tmc_dense_get(rate_constants, i_cell, arr), expected_arr, 1.0e-12 * expected_arr);

      /* r3 (user-defined); bar is parameterized */
      double expected_ud = ur * (cond->air_density * 0.82);
      TMC_EXPECT_NEAR(tmc_dense_get(rate_constants, i_cell, ud), expected_ud, 1.0e-12 * expected_ud);

      /* r2 (surface); foo is not parameterized */
      double mean_free_speed = sqrt(8.0 * GAS_CONSTANT / (M_PI * 0.025) * cond->temperature);
      double expected_surf = 4.0 * nc * M_PI * rad * rad / (rad / foo_diff_coeff + 4.0 / mean_free_speed);
      TMC_EXPECT_NEAR(tmc_dense_get(rate_constants, i_cell, surf), expected_surf, 1.0e-10 * expected_surf);
    }
    free(user_rate);
    free(radius);
    free(num_conc);
    tmc_state_free(state);
  }
  tmc_solver_free(s);
  tmc_builder_free(b);
  tmc_system_free(&system);
  for (int i = 0; i < 3; ++i)
    tmc_reaction_free(&r[i]);
  tmc_phase_free(&gas_phase);
}

TMC_TEST(process_matrix)
{
  process_update_state(1, 5);
  return TMC_TEST_PASS;
}

TMC_TEST(process_vector_matrix)
{
  for (int L = 1; L <= 4; ++L)
    process_update_state(L, 5);
  return TMC_TEST_PASS;
}

TMC_TEST(builds_chemical_reaction)
{
  tmc_species O3 = tmc_species_create("O3");
  tmc_species NO = tmc_species_create("NO");
  tmc_species NO2 = tmc_species_create("NO2");
  tmc_species O2 = tmc_species_create("O2");
  tmc_phase gas_phase = tmc_phase_create("gas");
  tmc_phase_add_species(&gas_phase, &O3);
  tmc_phase_add_species(&gas_phase, &NO);
  tmc_phase_add_species(&gas_phase, &NO2);
  tmc_phase_add_species(&gas_phase, &O2);

  tmc_reaction reaction = tmc_reaction_create();
  tmc_species reactant_O3 = tmc_species_create("O3"), reactant_NO = tmc_species_create("NO");
  tmc_species product_NO2 = tmc_species_create("NO2"), product_O2 = tmc_species_create("O2");
  tmc_reaction_add_reactant(&reaction, &reactant_O3);
  tmc_reaction_add_reactant(&reaction, &reactant_NO);
  tmc_reaction_add_product(&reaction, &product_NO2, 1.0);
  tmc_reaction_add_product(&reaction, &product_O2, 1.0);
  tmc_reaction_set_rate_constant(&reaction, tmc_rate_arrhenius(tmc_arrhenius_defaults()));
  tmc_reaction_set_phase(&reaction, &gas_phase);
  /* MICM's ChemicalReactionBuilder::Build validates the reaction */
  tmc_error err;
  TMC_EXPECT_EQ(tmc_reaction_validate(&reaction, &err), TMC_OK);

  TMC_EXPECT_EQ(reaction.n_reactants, 2);
  TMC_EXPECT_EQ(reaction.n_products, 2);
  TMC_EXPECT_STREQ(reaction.products[0].species.name, "NO2");
  TMC_EXPECT_EQ(reaction.products[0].coefficient, 1.0);
  TMC_EXPECT_STREQ(reaction.phase, "gas");
  TMC_EXPECT_TRUE(reaction.has_rate_constant);
  TMC_EXPECT_EQ(reaction.rate_constant.type, TMC_RATE_ARRHENIUS);

  /* and Build refuses a reaction without a rate constant */
  tmc_reaction no_rate = tmc_reaction_create();
  tmc_reaction_add_reactant(&no_rate, &reactant_O3);
  TMC_EXPECT_EQ(tmc_reaction_validate(&no_rate, &err), TMC_ERROR_PROCESS);
  TMC_EXPECT_EQ(err.code, TMC_PROCESS_RATE_CONSTANT_IS_NOT_SET);

  tmc_reaction_free(&no_rate);
  tmc_reaction_free(&reaction);
  tmc_phase_free(&gas_phase);
  return TMC_TEST_PASS;
}

TMC_TEST(chemical_reaction_copy_assignment)
{
  tmc_species foo = tmc_species_create_with_property("foo", TMC_PROPERTY_MOLECULAR_WEIGHT, 0.025);
  tmc_species_set_double(&foo, TMC_PROPERTY_DIFFUSION_COEFFICIENT, 2.3e2);
  tmc_species bar = tmc_species_create("bar");
  tmc_species_set_parameterization(&bar, 0.0, 0.0, 0.0, 0.82);
  tmc_phase gas_phase = tmc_phase_create("gas");
  tmc_phase_add_species(&gas_phase, &foo);
  tmc_phase_add_species(&gas_phase, &bar);

  tmc_arrhenius_parameters params = tmc_arrhenius_defaults();
  params.A = 12.2;
  params.C = 300.0;
  tmc_reaction reaction = tmc_reaction_create();
  tmc_reaction_add_reactant(&reaction, &foo);
  tmc_reaction_add_reactant(&reaction, &bar);
  tmc_reaction_set_rate_constant(&reaction, tmc_rate_arrhenius(params));
  tmc_reaction_set_phase(&reaction, &gas_phase);

  tmc_reaction copy = tmc_reaction_copy(&reaction);
  TMC_EXPECT_STREQ(copy.reactants[0].name, reaction.reactants[0].name);
  TMC_EXPECT_EQ(copy.n_reactants, reaction.n_reactants);
  TMC_EXPECT_EQ(copy.n_products, reaction.n_products);
  TMC_EXPECT_STREQ(copy.phase, reaction.phase);
  /* value semantics: the rate constant is copied by value */
  TMC_EXPECT_EQ(copy.rate_constant.type, TMC_RATE_ARRHENIUS);
  TMC_EXPECT_EQ(reaction.rate_constant.type, TMC_RATE_ARRHENIUS);
  TMC_EXPECT_EQ(copy.rate_constant.p.arrhenius.A, reaction.rate_constant.p.arrhenius.A);
  TMC_EXPECT_TRUE(copy.reactants[1].parameterized);
  TMC_EXPECT_EQ(copy.reactants[1].c_rho, 0.82);

  /* the copy owns its own arrays */
  TMC_EXPECT_TRUE(copy.reactants != reaction.reactants);
  copy.reactants[0].name[0] = 'X';
  copy.rate_constant.p.arrhenius.A = 1.0;
  TMC_EXPECT_STREQ(reaction.reactants[0].name, "foo");
  TMC_EXPECT_EQ(reaction.rate_constant.p.arrhenius.A, 12.2);

  tmc_reaction_free(&copy);
  tmc_reaction_free(&reaction);
  tmc_phase_free(&gas_phase);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("process:Process.Matrix", process_matrix),
  TMC_CASE("process:Process.VectorMatrix", process_vector_matrix),
  TMC_CASE("process:Process.BuildsChemicalReaction", builds_chemical_reaction),
  TMC_CASE("process:Process.ChemicalReactionCopyAssignmentSucceeds", chemical_reaction_copy_assignment),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_process(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
