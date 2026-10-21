#include "chemistry.h"

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>

static double test_lambda(double temperature, double pressure,
                          double air_density, const double *parameters) {
  (void)parameters;
  (void)pressure;
  (void)air_density;
  return 2.0 * temperature;
}

static double test_lambda_override(double temperature, double pressure,
                                   double air_density,
                                   const double *parameters) {
  (void)pressure;
  (void)air_density;
  (void)parameters;
  return 2.0e-3 * temperature;
}

static int near(double actual, double expected) {
  return fabs(actual - expected) <= 1e-13 * fabs(expected) + 1e-30;
}

static int tutorial_print_match(double actual, double expected) {
  if (expected == 0.0) return actual == 0.0;
  const double tolerance = 0.5 * pow(10.0, floor(log10(expected)) - 2.0);
  return fabs(actual - expected) <= tolerance;
}

static int check_solver_workspace_argument_errors(void) {
  chem_solver_workspace workspace = {0};
  return chem_solver_workspace_init(&workspace,NULL) != 0 &&
         chem_solver_workspace_init(NULL,NULL) != 0 ? 0 : 1;
}

static int check_mechanism_species_usage(void) {
  const chem_species species[] = {
      {.name="A"}, {.name="B"}, {.name="C"}, {.name="D"}};
  const chem_term reactants[] = {{.species=0,.coefficient=1.0}};
  const chem_term products[] = {{.species=1,.coefficient=1.0}};
  const chem_reaction reaction = {.kind=CHEM_USER_DEFINED,
      .reactants=reactants,.reactant_count=1,
      .products=products,.product_count=1};
  const chem_mechanism mechanism = {.species=species,.species_count=4,
      .reactions=&reaction,.reaction_count=1};
  unsigned char used[4];
  if (chem_mechanism_validate_species_usage(&mechanism,used,4) != 2 ||
      used[0] != 1 || used[1] != 1 || used[2] != 0 || used[3] != 0 ||
      chem_mechanism_validate_species_usage(&mechanism,used,3) == 0)
    return 1;
  const chem_term c_reactant[] = {{.species=0,.coefficient=1.0},
                                  {.species=2,.coefficient=1.0}};
  const chem_term cd_products[] = {{.species=1,.coefficient=1.0},
                                   {.species=3,.coefficient=1.0}};
  const chem_reaction all_species_reaction = {.kind=CHEM_USER_DEFINED,
      .reactants=c_reactant,.reactant_count=2,
      .products=cd_products,.product_count=2};
  const chem_mechanism all_used = {.species=species,.species_count=4,
      .reactions=&all_species_reaction,.reaction_count=1};
  if (chem_mechanism_validate_species_usage(&all_used,used,4) != 0 ||
      used[0] != 1 || used[1] != 1 || used[2] != 1 || used[3] != 1)
    return 2;
  return 0;
}

static int check_rate_constants_no_user_tutorial(void) {
  const chem_species species[] = {
      {.name="A"}, {.name="B"},
      {.name="C", .molecular_weight=0.025, .diffusion=230.0},
      {.name="D"}, {.name="E"}, {.name="F"}, {.name="G"}};
  const chem_term a[] = {{.species=0, .coefficient=1.0}};
  const chem_term b[] = {{.species=1, .coefficient=1.0}};
  const chem_term c[] = {{.species=2, .coefficient=1.0}};
  const chem_term d[] = {{.species=3, .coefficient=1.0}};
  const chem_term e_e[] = {
      {.species=4, .coefficient=1.0}, {.species=4, .coefficient=1.0}};
  const chem_term f[] = {{.species=5, .coefficient=1.0}};
  const chem_term product_b[] = {{.species=1, .coefficient=1.0}};
  const chem_term product_a[] = {{.species=0, .coefficient=1.0}};
  const chem_term product_c[] = {{.species=2, .coefficient=1.0}};
  const chem_term product_d[] = {{.species=3, .coefficient=1.0}};
  const chem_term product_e[] = {{.species=4, .coefficient=1.0}};
  const chem_term product_f2[] = {{.species=5, .coefficient=2.0}};
  const chem_term product_g[] = {{.species=6, .coefficient=1.0}};
  const double conversion = 6.02214076e17;
  const double reference_a =
      2.0e-22 * 6.02214076e23 * 1.0e-6 * exp(2.0) * (2.45e19/conversion);
  const double reference_b = 0.43 * pow(293.0/298.0, -8.0);
  const double reference_ratio = reference_a/reference_b;
  const double av = reference_a/(1.0+reference_ratio) *
      pow(0.41, 1.0/(1.0+pow(log10(reference_ratio), 2.0)));
  const double branch_z = av * (1.0-1.0e-3)/1.0e-3;
  const chem_reaction reactions[] = {
      {.kind=CHEM_ARRHENIUS, .reactants=a, .reactant_count=1,
       .products=product_b, .product_count=1,
       .p={2.15e-4, 0.0, 110.0, 300.0, 0.0}},
      {.kind=CHEM_BRANCHED_ALKOXY, .reactants=b, .reactant_count=1,
       .products=product_c, .product_count=1,
       .p={1.2, 204.3, 1.0e-3, 2.0, branch_z}},
      {.kind=CHEM_BRANCHED_NITRATE, .reactants=b, .reactant_count=1,
       .products=product_d, .product_count=1,
       .p={1.2, 204.3, 1.0e-3, 2.0, branch_z}},
      {.kind=CHEM_SURFACE, .reactants=c, .reactant_count=1,
       .products=product_e, .product_count=1, .surface_species=2,
       .p={0.90, 0, 0, 0, 0, 0, 0, 0}},
      {.kind=CHEM_TERNARY_ACTIVATION, .reactants=d, .reactant_count=1,
       .products=product_f2, .product_count=1,
       .p={1.2, 2.3, 302.3, 2.6/conversion, -3.1, 402.1, 0.9, 1.2}},
      {.kind=CHEM_TROE, .reactants=e_e, .reactant_count=2,
       .products=product_g, .product_count=1,
       .p={1.8e-30*conversion*conversion, -3.0, 0.0,
           2.8e-11*conversion, 0.0, 0.0, 0.6, 1.0}},
      {.kind=CHEM_TUNNELING, .reactants=f, .reactant_count=1,
       .products=product_g, .product_count=1, .p={1.2, 2.3, 302.3}},
  };
  const chem_mechanism mechanism = {
      .species=species, .species_count=7,
      .reactions=reactions, .reaction_count=7};
  const double rate_parameters[] = {1.0e-7, 2.5e6};
  const double expected[10][7] = {
      {8.54e-1, 4.57e-4, 1.44e-1, 1.55e-4, 0.0, 1.23e-22, 6.44e-4},
      {7.30e-1, 3.90e-4, 2.65e-1, 2.89e-4, 8.18e-7, 2.28e-22, 2.44e-3},
      {6.23e-1, 3.33e-4, 3.66e-1, 4.02e-4, 9.66e-7, 3.18e-22, 5.20e-3},
      {5.32e-1, 2.85e-4, 4.49e-1, 5.00e-4, 1.07e-6, 3.95e-22, 8.77e-3},
      {4.55e-1, 2.43e-4, 5.18e-1, 5.83e-4, 1.15e-6, 4.61e-22, 1.30e-2},
      {3.88e-1, 2.08e-4, 5.75e-1, 6.54e-4, 1.21e-6, 5.17e-22, 1.78e-2},
      {3.32e-1, 1.77e-4, 6.21e-1, 7.14e-4, 1.26e-6, 5.65e-22, 2.30e-2},
      {2.83e-1, 1.52e-4, 6.59e-1, 7.66e-4, 1.30e-6, 6.06e-22, 2.86e-2},
      {2.42e-1, 1.29e-4, 6.88e-1, 8.10e-4, 1.33e-6, 6.41e-22, 3.45e-2},
      {2.07e-1, 1.11e-4, 7.11e-1, 8.48e-4, 1.35e-6, 6.71e-22, 4.06e-2}};
  double state[] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double temperature = 287.45, pressure = 101319.9;
  const double density = pressure/(8.31446261815324*temperature);
  chem_solver_workspace workspace = {0};
  if (chem_solver_workspace_init(&workspace, &mechanism) != 0) return 1;
  int result = 0;
  for (size_t step = 0; step < 10 && !result; ++step) {
    result = chem_rosenbrock3_solve_with_workspace(
        &mechanism, &workspace, state, temperature, pressure, density,
        rate_parameters, 500.0, NULL, NULL);
    for (size_t i = 0; i < 7 && !result; ++i)
      if (!tutorial_print_match(state[i], expected[step][i])) result = 1;
  }
  chem_solver_workspace_destroy(&workspace);
  if (result) return result;

  /* The user-defined variant adds changing photolysis and fixed emission/loss
     rates to the same seven base reactions. */
  chem_reaction user_reactions[10];
  memcpy(user_reactions, reactions, sizeof(reactions));
  user_reactions[7] = (chem_reaction){
      .kind=CHEM_USER_DEFINED, .label="my photolysis rate",
      .reactants=c, .reactant_count=1, .products=product_g,
      .product_count=1, .p={1, 0, 0, 0, 0, 0, 0, 2}};
  user_reactions[8] = (chem_reaction){
      .kind=CHEM_USER_DEFINED, .label="my emission rate",
      .reactants=NULL, .reactant_count=0, .products=product_a,
      .product_count=1, .p={1, 0, 0, 0, 0, 0, 0, 3}};
  user_reactions[9] = (chem_reaction){
      .kind=CHEM_USER_DEFINED, .label="my loss rate",
      .reactants=b, .reactant_count=1, .products=NULL,
      .product_count=0, .p={1, 0, 0, 0, 0, 0, 0, 4}};
  const char *custom_labels[] = {
      "C.effective radius [m]", "C.particle number concentration [# m-3]",
      "my photolysis rate", "my emission rate", "my loss rate"};
  const chem_mechanism user_mechanism = {
      .species=species, .species_count=7,
      .reactions=user_reactions, .reaction_count=10,
      .custom_labels=custom_labels, .custom_parameter_count=5};
  const chem_solver_options options = chem_solver_default_options();
  double user_state[] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double photo_rate = 1.0e-10;
  chem_solver_workspace user_workspace = {0};
  if (chem_solver_workspace_init(&user_workspace, &user_mechanism) != 0) return 2;
  for (size_t step = 0; step < 10 && !result; ++step) {
    const double user_rates[] = {1.0e-7, 2.5e6, photo_rate, 1.0e-20, 1.0e-23};
    result = chem_rosenbrock3_solve_with_workspace(
        &user_mechanism, &user_workspace, user_state, temperature, pressure,
        density, user_rates, 500.0, &options, NULL);
    for (size_t i = 0; i < 7 && !result; ++i)
      if (!tutorial_print_match(user_state[i], expected[step][i])) result = 1;
    photo_rate *= 1.5;
  }
  chem_solver_workspace_destroy(&user_workspace);
  return result;
}

static int check_lambda_rate_constant_position(void) {
  const chem_term reactant[] = {{.species = 0, .coefficient = 1.0}};
  const chem_term product[] = {{.species = 1, .coefficient = 1.0}};
  const chem_reaction reactions[] = {
      {.kind = CHEM_ARRHENIUS,
       .reactants = reactant, .reactant_count = 1,
       .products = product, .product_count = 1,
       .p = {2.0, 0.0, 0.0, 300.0, 0.0}},
      {.kind = CHEM_LAMBDA,
       .reactants = reactant, .reactant_count = 1,
       .products = product, .product_count = 1,
       .rate_function = test_lambda}};
  const chem_species species[] = {{.name = "a"}, {.name = "b"}};
  const chem_mechanism mechanism = {
      .species = species, .species_count = 2,
      .reactions = reactions, .reaction_count = 2};
  double constants[2] = {0.0, 0.0};
  if (chem_rate_constants(&mechanism, 300.0, 101325.0, 1.0, NULL,
                          constants) != 0)
    return 1;
  return near(constants[0], 2.0) && near(constants[1], 600.0) ? 0 : 1;
}

static void simple_forcing(const double *x, size_t blocks, double *f) {
  for (size_t b = 0; b < blocks; ++b) {
    const double a = x[2*b], c = x[2*b + 1];
    f[2*b] = a * c;
    f[2*b + 1] = a * a - c;
  }
}

static void finite_difference_jacobian(const double *x, size_t blocks,
                                       double *jacobian) {
  const double epsilon = 1e-8;
  for (size_t col = 0; col < 2; ++col) {
    double plus[6], minus[6], f_plus[6], f_minus[6], f_base[6];
    memcpy(plus, x, blocks * 2 * sizeof(double));
    memcpy(minus, x, blocks * 2 * sizeof(double));
    for (size_t b = 0; b < blocks; ++b) {
      const double value = x[2*b + col];
      const double h = epsilon * fmax(1.0, fabs(value));
      plus[2*b + col] = value + h;
      minus[2*b + col] = value - h;
    }
    simple_forcing(plus, blocks, f_plus);
    simple_forcing(minus, blocks, f_minus);
    simple_forcing(x, blocks, f_base);
    for (size_t b = 0; b < blocks; ++b) {
      const double value = x[2*b + col];
      const double h = epsilon * fmax(1.0, fabs(value));
      for (size_t row = 0; row < 2; ++row) {
        jacobian[4*b + 2*row + col] = value - h >= 0.0
            ? (f_plus[2*b + row] - f_minus[2*b + row]) / (2.0*h)
            : (f_plus[2*b + row] - f_base[2*b + row]) / h;
      }
    }
  }
}

static int check_jacobian_verification(void) {
  const double points[] = {3.0, 5.0, 0.5, 4.0, 10.0, 0.1};
  const double expected[] = {5.0, 3.0, 6.0, -1.0,
                             4.0, 0.5, 1.0, -1.0,
                             0.1, 10.0, 20.0, -1.0};
  double fd[12] = {0};
  finite_difference_jacobian(points, 3, fd);
  for (size_t i = 0; i < 12; ++i)
    if (fabs(fd[i] - expected[i]) > 1e-5) return 1;

  /* Compare MICM's negative-Jacobian convention with finite differences. */
  const double negative_jacobian[] = {-5.0, -3.0, -6.0, 1.0};
  for (size_t i = 0; i < 4; ++i)
    if (fabs(-negative_jacobian[i] - fd[i]) > 1e-7 + 1e-7*fabs(fd[i])) return 2;
  double wrong[4] = {-999.0, -3.0, -6.0, 1.0};
  size_t worst = 0;
  for (size_t i = 1; i < 4; ++i)
    if (fabs(-wrong[i] - fd[i]) > fabs(-wrong[worst] - fd[worst])) worst = i;
  if (fabs(-wrong[worst] - fd[worst]) <= 1e-7 + 1e-7*fabs(fd[worst]) ||
      worst != 0) return 3;

  /* A full sparsity pattern is complete; omitting (1,0) is detected. */
  unsigned char full_pattern[] = {1, 1, 1, 1};
  unsigned char incomplete_pattern[] = {1, 1, 0, 1};
  for (size_t i = 0; i < 4; ++i)
    if (!full_pattern[i] && fabs(fd[i]) > 1e-5) return 4;
  size_t missing = 4;
  for (size_t i = 0; i < 4; ++i)
    if (!incomplete_pattern[i] && fabs(fd[i]) > 1e-5) missing = i;
  if (missing != 2) return 5;

  const double near_zero[] = {0.0, 1.0};
  double zero_fd[4] = {0};
  finite_difference_jacobian(near_zero, 1, zero_fd);
  const double zero_expected[] = {1.0, 0.0, 0.0, -1.0};
  for (size_t i = 0; i < 4; ++i)
    if (fabs(zero_fd[i] - zero_expected[i]) > 1e-5) return 6;
  return 0;
}

static int check_version_api(void) {
  return strcmp(chem_version(), "0.1.0+0") == 0 &&
         chem_version_major() == 0 && chem_version_minor() == 1 &&
         chem_version_patch() == 0 && chem_version_tweak() == 0 ? 0 : 1;
}

static int check_one_solver_many_thread_states(void) {
  const chem_species species[] = {{.name = "A"}, {.name = "B"}, {.name = "C"}};
  const chem_term a[] = {{.species = 0, .coefficient = 1.0}};
  const chem_term b[] = {{.species = 1, .coefficient = 1.0}};
  const chem_term bb[] = {{.species = 1, .coefficient = 2.0}};
  const chem_term bc[] = {{.species = 1, .coefficient = 1.0},
                          {.species = 2, .coefficient = 1.0}};
  const chem_term product_b[] = {{.species = 1, .coefficient = 1.0}};
  const chem_term product_bc[] = {{.species = 1, .coefficient = 1.0},
                                  {.species = 2, .coefficient = 1.0}};
  const chem_term product_c[] = {{.species = 2, .coefficient = 1.0}};
  const chem_term product_ac[] = {{.species = 0, .coefficient = 1.0},
                                  {.species = 2, .coefficient = 1.0}};
  const chem_reaction reactions[] = {
      {.kind = CHEM_USER_DEFINED, .reactants = a, .reactant_count = 1,
       .products = product_b, .product_count = 1, .p = {1,0,0,0,0,0,0,0}},
      {.kind = CHEM_USER_DEFINED, .reactants = bb, .reactant_count = 1,
       .products = product_bc, .product_count = 2, .p = {1,0,0,0,0,0,0,1}},
      {.kind = CHEM_USER_DEFINED, .reactants = bc, .reactant_count = 2,
       .products = product_ac, .product_count = 2, .p = {1,0,0,0,0,0,0,2}}};
  (void)b;
  (void)product_c;
  const chem_mechanism mechanism = {
      .species = species, .species_count = 3,
      .reactions = reactions, .reaction_count = 3};
  double states[8][3] = {{0}};
  int statuses[8] = {0};
  const double rates[] = {0.04, 3.0e7, 1.0e4};
  #pragma omp parallel num_threads(8)
  {
    const int thread = omp_get_thread_num();
    states[thread][0] = 1.0;
    chem_solver_workspace workspace = {0};
    statuses[thread] = chem_solver_workspace_init(&workspace, &mechanism);
    if (!statuses[thread]) {
      for (size_t step = 0; step < 10 && !statuses[thread]; ++step)
        statuses[thread] = chem_rosenbrock3_solve_with_workspace(
            &mechanism, &workspace, states[thread], 272.5, 101253.3, 1.0e6,
            rates, 200.0, NULL, NULL);
    }
    chem_solver_workspace_destroy(&workspace);
  }
  for (size_t thread = 0; thread < 8; ++thread) {
    if (statuses[thread]) return 1;
    for (size_t species_index = 0; species_index < 3; ++species_index)
      if (states[thread][species_index] != states[0][species_index]) return 2;
  }
  return 0;
}

static int check_state_view(void) {
  chem_state_view empty_state = {0};
  if (empty_state.cell_count != 0 || empty_state.species_names != NULL ||
      empty_state.concentrations != NULL) return 13;
  const char *species[] = {"foo", "bar", "baz", "quz"};
  const char *parameters[] = {"O1", "O2", "O3", "AAA", "BBB"};
  double values[12] = {0};
  double custom[15] = {0};
  chem_state_view state;
  if (chem_state_view_init(&state, species, 4, parameters, 5, 3, values,
                           custom)) return 1;
  if (state.cell_count != 3 || state.species_count != 4 ||
      state.custom_parameter_count != 5 || state.concentrations != values ||
      state.custom_parameters != custom) return 14;

  const double bar[] = {12.0, 42.0, 35.2};
  if (chem_state_set_concentration(&state, "bar", bar, 3)) return 2;
  for (size_t cell = 0; cell < 3; ++cell) {
    double observed = 0.0;
    if (chem_state_get_concentration(&state, "bar", cell, &observed) ||
        observed != bar[cell]) return 3;
  }
  const double foo[] = {0.9, 0.99, 0.999};
  const double baz[] = {0.3, 0.33, 0.333};
  const double quz[] = {0.4, 0.44, 0.444};
  const chem_state_series concentrations[] = {
      {.label = "quz", .values = quz, .value_count = 3},
      {.label = "foo", .values = foo, .value_count = 3},
      {.label = "baz", .values = baz, .value_count = 3}};
  if (chem_state_set_concentrations(&state, concentrations, 3)) return 4;
  const double expected[] = {0.9, 12.0, 0.3, 0.4,
                             0.99, 42.0, 0.33, 0.44,
                             0.999, 35.2, 0.333, 0.444};
  for (size_t i = 0; i < 12; ++i)
    if (values[i] != expected[i]) return 5;
  const double invalid[] = {0.1};
  const chem_state_series unknown_species[] = {
      {.label = "FUU", .values = invalid, .value_count = 1}};
  if (chem_state_set_concentrations(&state, unknown_species, 1) == 0 ||
      chem_state_set_concentration(&state, "missing", invalid, 3) == 0 ||
      chem_state_set_concentration(&state, "foo", invalid, 1) == 0) return 6;

  const double o1[] = {0.1, 0.11, 0.111};
  const double o2[] = {0.5, 0.55, 0.555};
  const double o3[] = {0.3, 0.33, 0.333};
  const double aaa[] = {0.5, 0.55, 0.555};
  const double bbb[] = {0.7, 0.77, 0.777};
  const chem_state_series custom_series[] = {
      {.label = "O3", .values = o3, .value_count = 3},
      {.label = "O1", .values = o1, .value_count = 3},
      {.label = "BBB", .values = bbb, .value_count = 3},
      {.label = "O2", .values = o2, .value_count = 3},
      {.label = "AAA", .values = aaa, .value_count = 3}};
  if (chem_state_set_custom_parameters(&state, custom_series, 5)) return 7;
  const double expected_custom[] = {0.1,0.5,0.3,0.5,0.7,
                                   0.11,0.55,0.33,0.55,0.77,
                                   0.111,0.555,0.333,0.555,0.777};
  for (size_t i = 0; i < 15; ++i)
    if (custom[i] != expected_custom[i]) return 8;
  const chem_state_series bad_custom_size[] = {
      {.label = "O1", .values = invalid, .value_count = 1}};
  const chem_state_series bad_custom_label[] = {
      {.label = "CCC", .values = o1, .value_count = 3}};
  if (chem_state_set_custom_parameters(&state, bad_custom_size, 1) == 0 ||
      chem_state_set_custom_parameters(&state, bad_custom_label, 1) == 0 ||
      chem_state_set_custom_parameter(&state, "CCC", o1, 3) == 0) return 9;

  /* The caller-owned flat matrix is the C equivalent of MICM's unchecked bulk setter. */
  if (chem_state_set_custom_parameter_matrix(&state, expected_custom, 3, 5) ||
      chem_state_set_custom_parameter_matrix(&state, expected_custom, 2, 5) == 0 ||
      chem_state_set_custom_parameter_matrix(&state, expected_custom, 3, 4) == 0)
    return 10;

  double one_cell_values[4] = {0};
  double one_cell_custom[5] = {0};
  chem_state_view one_cell;
  if (chem_state_view_init(&one_cell, species, 4, parameters, 5, 1,
                           one_cell_values, one_cell_custom)) return 11;
  const double scalar[] = {324.2};
  const double param[] = {42.3};
  if (chem_state_set_concentration(&one_cell, "bar", scalar, 1) ||
      chem_state_set_custom_parameter(&one_cell, "O2", param, 1) ||
      one_cell_values[1] != 324.2 || one_cell_custom[1] != 42.3) return 12;
  return 0;
}

static int check_process_rate_constant_layouts(void) {
  const chem_species species[] = {
      {.name = "foo", .molecular_weight = 0.025, .diffusion = 230.0},
      {.name = "bar", .density_parameterized = 1, .density_factor = 0.82}};
  const chem_term foo[] = {{.species = 0, .coefficient = 1.0}};
  const chem_term bar[] = {{.species = 1, .coefficient = 1.0}};
  const chem_term foo_bar[] = {{.species = 0, .coefficient = 1.0},
                               {.species = 1, .coefficient = 1.0}};
  const chem_reaction reactions[] = {
      {.kind = CHEM_ARRHENIUS, .reactants = foo_bar, .reactant_count = 2,
       .p = {12.2, 0.0, 300.0, 300.0, 0.0}},
      {.kind = CHEM_USER_DEFINED, .reactants = bar, .reactant_count = 1,
       .p = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      {.kind = CHEM_SURFACE, .reactants = foo, .reactant_count = 1,
       .surface_species = 0, .p = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}}};
  const chem_mechanism mechanism = {
      .species = species, .species_count = 2,
      .reactions = reactions, .reaction_count = 3};
  const double r = 8.31446261815324;
  double scalar_results[5][3] = {{0}};
  for (size_t cell = 0; cell < 5; ++cell) {
    const double temperature = 280.0 + 2.0 * (double)cell;
    const double pressure = 99000.0 + 100.0 * (double)cell;
    const double density = 6.0 + 0.1 * (double)cell;
    const double user_rate = 0.01 + 0.001 * (double)cell;
    const double radius = 1.0e-8 * (1.0 + 0.1 * (double)cell);
    const double number = 1.0e5 * (1.0 + 0.05 * (double)cell);
    const double custom[] = {user_rate, radius, number};
    double constants[3] = {0};
    if (chem_rate_constants(&mechanism, temperature, pressure, density,
                            custom, constants)) return 1;
    const double expected_arrhenius = 12.2 * exp(300.0 / temperature) *
                                      (density * 0.82);
    const double expected_user = user_rate * (density * 0.82);
    const double speed = sqrt(8.0 * r * temperature /
                              (3.14159265358979323846 * 0.025));
    const double expected_surface =
        4.0 * number * 3.14159265358979323846 * radius * radius /
        (radius / 230.0 + 4.0 / speed);
    if (!near(constants[0], expected_arrhenius) ||
        !near(constants[1], expected_user) ||
        fabs(constants[2] - expected_surface) > 1e-10 * expected_surface)
      return 2;
    memcpy(scalar_results[cell], constants, sizeof(constants));
  }

  /* The C API's cell-major arrays are valid for scalar and grouped cell widths. */
  for (size_t width = 1; width <= 4; ++width) {
    size_t offset = 0;
    while (offset < 5) {
      const size_t count = width < 5 - offset ? width : 5 - offset;
      double grouped[4][3] = {{0}};
      for (size_t local = 0; local < count; ++local) {
        const size_t cell = offset + local;
        const double temperature = 280.0 + 2.0 * (double)cell;
        const double pressure = 99000.0 + 100.0 * (double)cell;
        const double density = 6.0 + 0.1 * (double)cell;
        const double custom[] = {0.01 + 0.001 * (double)cell,
                                 1.0e-8 * (1.0 + 0.1 * (double)cell),
                                 1.0e5 * (1.0 + 0.05 * (double)cell)};
        if (chem_rate_constants(&mechanism, temperature, pressure, density,
                                custom, grouped[local])) return 3;
      }
      for (size_t local = 0; local < count; ++local)
        for (size_t rate = 0; rate < 3; ++rate)
          if (!isfinite(grouped[local][rate]) ||
              !near(grouped[local][rate], scalar_results[offset + local][rate]))
            return 4;
      offset += count;
    }
  }
  return 0;
}

static int check_process_set_finite_difference_jacobian(void) {
  const chem_species species[] = {{.name = "A"}, {.name = "B"}, {.name = "C"}};
  const chem_term ab[] = {{.species = 0, .coefficient = 1.0},
                          {.species = 1, .coefficient = 1.0}};
  const chem_term c[] = {{.species = 2, .coefficient = 1.0}};
  const chem_term two_c[] = {{.species = 2, .coefficient = 2.0}};
  const chem_term a[] = {{.species = 0, .coefficient = 1.0}};
  const chem_reaction reactions[] = {
      {.kind = CHEM_ARRHENIUS, .reactants = ab, .reactant_count = 2,
       .products = two_c, .product_count = 1},
      {.kind = CHEM_ARRHENIUS, .reactants = c, .reactant_count = 1,
       .products = a, .product_count = 1}};
  const chem_mechanism mechanism = {
      .species = species, .species_count = 3,
      .reactions = reactions, .reaction_count = 2};
  const double constants[] = {1.3, 0.8};
  const double states[][3] = {{1.2, 0.7, 0.5}, {0.4, 1.1, 0.9}};
  for (size_t cell = 0; cell < 2; ++cell) {
    double jacobian[9] = {0};
    if (chem_jacobian_with_constants(&mechanism, states[cell], constants,
                                     1.0e6, jacobian)) return 1;
    for (size_t col = 0; col < 3; ++col) {
      double plus[3], minus[3], f_plus[3], f_minus[3];
      memcpy(plus, states[cell], sizeof(plus));
      memcpy(minus, states[cell], sizeof(minus));
      const double h = 1.0e-6 * fmax(1.0, fabs(states[cell][col]));
      plus[col] += h;
      minus[col] -= h;
      if (chem_forcing_with_constants(&mechanism, plus, constants, 1.0e6,
                                      f_plus) ||
          chem_forcing_with_constants(&mechanism, minus, constants, 1.0e6,
                                      f_minus)) return 2;
      for (size_t row = 0; row < 3; ++row) {
        const double finite_difference = (f_plus[row] - f_minus[row]) / (2*h);
        if (fabs(jacobian[row*3 + col] - finite_difference) >
            1.0e-8 * fmax(1.0, fabs(finite_difference))) return 3;
      }
    }
  }
  return 0;
}

static int check_random_process_set(size_t cells, size_t reaction_count,
                                    size_t species_count) {
  chem_species *species = calloc(species_count, sizeof(*species));
  chem_reaction *reactions = calloc(reaction_count, sizeof(*reactions));
  chem_term *reactants = calloc(reaction_count * 3, sizeof(*reactants));
  chem_term *products = calloc(reaction_count * 10, sizeof(*products));
  double *concentration = calloc(species_count, sizeof(*concentration));
  double *constants = calloc(reaction_count, sizeof(*constants));
  double *forcing = calloc(species_count, sizeof(*forcing));
  if (!species || !reactions || !reactants || !products || !concentration ||
      !constants || !forcing) {
    free(species); free(reactions); free(reactants); free(products);
    free(concentration); free(constants); free(forcing);
    return 1;
  }
  for (size_t r = 0; r < reaction_count; ++r) {
    const size_t n_reactants = r % 4;
    const size_t n_products = r % 11;
    reactions[r].kind = CHEM_ARRHENIUS;
    reactions[r].reactants = reactants + 3*r;
    reactions[r].reactant_count = n_reactants;
    reactions[r].products = products + 10*r;
    reactions[r].product_count = n_products;
    for (size_t i = 0; i < n_reactants; ++i)
      reactants[3*r + i] = (chem_term){
          .species = (17*r + 31*i + 7) % species_count,
          .coefficient = 1.0};
    for (size_t i = 0; i < n_products; ++i)
      products[10*r + i] = (chem_term){
          .species = (19*r + 23*i + 11) % species_count,
          .coefficient = 1.2};
  }
  const chem_mechanism mechanism = {
      .species = species, .species_count = species_count,
      .reactions = reactions, .reaction_count = reaction_count};
  int result = 0;
  for (size_t cell = 0; cell < cells && !result; ++cell) {
    for (size_t i = 0; i < species_count; ++i)
      concentration[i] = 0.2 + 0.001 * (double)((i*37 + cell*13) % 300);
    for (size_t r = 0; r < reaction_count; ++r)
      constants[r] = 0.001 + 0.00001 * (double)((r*43 + cell*17) % 1000);
    result = chem_forcing_with_constants(&mechanism, concentration, constants,
                                         1.0e6, forcing);
    for (size_t i = 0; i < species_count && !result; ++i)
      if (!isfinite(forcing[i])) result = 2;
  }
  free(species); free(reactions); free(reactants); free(products);
  free(concentration); free(constants); free(forcing);
  return result;
}

static int check_random_process_sets(void) {
  return check_random_process_set(200, 50, 40) ||
         check_random_process_set(300, 30, 20) ||
         check_random_process_set(400, 100, 80);
}

static int check_process_set_forcing_and_jacobian(void) {
  const chem_species species[] = {
      {.name = "foo"}, {.name = "bar"}, {.name = "baz"}, {.name = "quz"},
      {.name = "quuz"}, {.name = "corge"},
      {.name = "qux", .density_parameterized = 1, .density_factor = 0.72}};
  const chem_term foo_baz[] = {{.species = 0, .coefficient = 1.0},
                               {.species = 2, .coefficient = 1.0}};
  const chem_term bar_qux[] = {{.species = 1, .coefficient = 1.0},
                               {.species = 6, .coefficient = 1.0}};
  const chem_term quz[] = {{.species = 3, .coefficient = 1.0}};
  const chem_term baz_qux[] = {{.species = 2, .coefficient = 1.0},
                               {.species = 6, .coefficient = 1.0}};
  const chem_term bar_quuz[] = {{.species = 1, .coefficient = 1.0},
                                {.species = 4, .coefficient = 2.4}};
  const chem_term foo_quz[] = {{.species = 0, .coefficient = 1.0},
                               {.species = 3, .coefficient = 1.4}};
  const chem_term bar_quz[] = {{.species = 1, .coefficient = 1.0},
                               {.species = 3, .coefficient = 2.5}};
  const chem_reaction reactions[] = {
      {.kind = CHEM_ARRHENIUS, .reactants = foo_baz, .reactant_count = 2,
       .products = bar_quuz, .product_count = 2},
      {.kind = CHEM_ARRHENIUS, .reactants = bar_qux, .reactant_count = 2,
       .products = foo_quz, .product_count = 2},
      {.kind = CHEM_ARRHENIUS, .reactants = quz, .reactant_count = 1,
       .products = NULL, .product_count = 0},
      {.kind = CHEM_ARRHENIUS, .reactants = baz_qux, .reactant_count = 2,
       .products = bar_quz, .product_count = 2}};
  const chem_mechanism mechanism = {
      .species = species, .species_count = 7,
      .reactions = reactions, .reaction_count = 4};
  const double y[2][7] = {{0.1, 0.2, 0.3, 0.4, 0.5, 0.0, 0.0},
                          {1.1, 1.2, 1.3, 1.4, 1.5, 0.0, 0.0}};
  const double rates[2][4] = {{10.0, 20.0*70.0*0.72, 30.0,
                               40.0*70.0*0.72},
                              {110.0, 120.0*80.0*0.72, 130.0,
                               140.0*80.0*0.72}};
  const double expected[2][6] = {
      {1000.0 - 10.0*0.1*0.3 + 20.0*70.0*0.72*0.2,
       1000.0 + 10.0*0.1*0.3 - 20.0*0.2*70.0*0.72 + 40.0*70.0*0.72*0.3,
       1000.0 - 10.0*0.1*0.3 - 40.0*70.0*0.72*0.3,
       1000.0 + 20.0*70.0*0.72*0.2*1.4 - 30.0*0.4 + 40.0*70.0*0.72*2.5*0.3,
       1000.0 + 10.0*0.1*0.3*2.4, 1000.0},
      {1000.0 - 110.0*1.1*1.3 + 120.0*80.0*0.72*1.2,
       1000.0 + 110.0*1.1*1.3 - 120.0*1.2*80.0*0.72 + 140.0*80.0*0.72*1.3,
       1000.0 - 110.0*1.1*1.3 - 140.0*80.0*0.72*1.3,
       1000.0 + 120.0*80.0*0.72*1.2*1.4 - 130.0*1.4 + 140.0*80.0*0.72*2.5*1.3,
       1000.0 + 110.0*1.1*1.3*2.4, 1000.0}};
  for (size_t cell = 0; cell < 2; ++cell) {
    double forcing[7] = {0};
    if (chem_forcing_with_constants(&mechanism, y[cell], rates[cell],
                                    70.0 + 10.0*cell, forcing)) return 1;
    for (size_t row = 0; row < 6; ++row)
      if (!near(1000.0 + forcing[row], expected[cell][row])) return 2;
    double negative_jacobian[49] = {0};
    if (chem_negative_jacobian_with_constants(&mechanism, y[cell], rates[cell],
                                              negative_jacobian)) return 3;
    for (size_t col = 0; col < 7; ++col) {
      double plus[7], minus[7], f_plus[7], f_minus[7];
      memcpy(plus, y[cell], sizeof(plus));
      memcpy(minus, y[cell], sizeof(minus));
      const double h = 1.0e-7 * fmax(1.0, fabs(y[cell][col]));
      plus[col] += h;
      minus[col] -= h;
      if (chem_forcing_with_constants(&mechanism, plus, rates[cell], 70.0,
                                      f_plus) ||
          chem_forcing_with_constants(&mechanism, minus, rates[cell], 70.0,
                                      f_minus)) return 4;
      for (size_t row = 0; row < 7; ++row) {
        const double fd = (f_plus[row] - f_minus[row]) / (2*h);
        if (fabs(negative_jacobian[row*7 + col] + fd) >
            1.0e-7 * fmax(1.0, fabs(fd))) return 5;
      }
    }
  }
  return 0;
}

static double terminator_photolysis(size_t cell, size_t cells, double lon) {
  const double pi = 3.14159265358979323846;
  const double lat = pi / 180.0 * ((double)cell * (90.0 / (double)cells));
  const double lat_center = pi * 20.0 / 180.0;
  const double lon_center = pi * 300.0 / 180.0;
  return fmax(0.0, sin(lat) * sin(lat_center) +
      cos(lat) * cos(lat_center) * cos(lon - lon_center));
}

static int check_terminator_solver(size_t cells, int grouped) {
  const chem_species species[] = {{.name = "Cl2"}, {.name = "Cl"}};
  const chem_term cl2[] = {{.species = 0, .coefficient = 1.0}};
  const chem_term two_cl[] = {{.species = 1, .coefficient = 2.0}};
  const chem_reaction reactions[] = {
      {.kind = CHEM_USER_DEFINED, .label = "toy_k1",
       .reactants = cl2, .reactant_count = 1,
       .products = two_cl, .product_count = 1,
       .p = {1.0, 0, 0, 0, 0, 0, 0, 0}},
      {.kind = CHEM_ARRHENIUS, .reactants = (const chem_term[]){{1, 2.0}},
       .reactant_count = 1, .products = cl2, .product_count = 1,
       .p = {1.0, 0, 0, 300.0, 0}}};
  const chem_mechanism mechanism = {
      .species = species, .species_count = 2,
      .reactions = reactions, .reaction_count = 2,
      .custom_parameter_count = 1};
  double states[4][2] = {{0}};
  double rates[4] = {0};
  double temperatures[4], pressures[4], densities[4];
  double initial[4][2];
  for (size_t cell = 0; cell < cells; ++cell) {
    initial[cell][0] = 1.0e-6 * (1.0 + 0.1 * (double)cell);
    initial[cell][1] = 1.0e-10 * (1.0 + 0.2 * (double)cell);
    temperatures[cell] = 298.0;
    pressures[cell] = 101300.0;
    densities[cell] = 42.0;
  }
  chem_solver_workspace workspace = {0};
  const int init_status = grouped
      ? chem_solver_workspace_init_group(&workspace, &mechanism, cells)
      : chem_solver_workspace_init(&workspace, &mechanism);
  if (init_status) return 1;
  chem_solver_options options = chem_solver_default_options();
  options.relative_tolerance = 1.0e-8;
  options.absolute_tolerance = 1.0e-20;
  options.max_steps = 100000;
  const size_t steps = (size_t)floor(2.0 * 3.14159265358979323846 / 0.3);
  int result = 0;
  for (size_t step = 0; step < steps && !result; ++step) {
    const double lon = (double)step * 0.3;
    for (size_t cell = 0; cell < cells; ++cell) {
      states[cell][0] = initial[cell][0];
      states[cell][1] = initial[cell][1];
      rates[cell] = terminator_photolysis(cell, cells, lon);
    }
    if (grouped) {
      result = chem_rosenbrock3_solve_group_with_workspace(
          &mechanism, &workspace, &states[0][0], temperatures, pressures,
          densities, rates, cells, 30.0, &options, NULL);
    } else {
      for (size_t cell = 0; cell < cells && !result; ++cell)
        result = chem_rosenbrock3_solve_with_workspace(
            &mechanism, &workspace, states[cell], temperatures[cell],
            pressures[cell], densities[cell], &rates[cell], 30.0, &options,
            NULL);
    }
    for (size_t cell = 0; cell < cells && !result; ++cell) {
      const double k1 = rates[cell], k2 = 1.0, dt = 30.0;
      const double r = k1 / (4.0 * k2);
      const double cl_i = initial[cell][1], cl2_i = initial[cell][0];
      const double cly = cl_i + 2.0 * cl2_i;
      const double det = sqrt(r*r + 2.0*r*cly);
      const double e = exp(-4.0*k2*det*dt);
      const double l = det*k2*dt > 1.0e-16 ? (1.0-e)/det/dt : 4.0*k2;
      const double cl_f = -l*(cl_i-det+r)*(cl_i+det+r) /
                          (1.0+e+dt*l*(cl_i+r));
      const double expected_cl = cl_i + dt*cl_f;
      const double expected_cl2 = cl2_i - dt*cl_f/2.0;
      if (fabs(states[cell][1] - expected_cl) >
              fabs(expected_cl)*1.0e-6 + 1.0e-14 ||
          fabs(states[cell][0] - expected_cl2) >
              fabs(expected_cl2)*1.0e-6 + 1.0e-14) result = 2;
    }
  }
  chem_solver_workspace_destroy(&workspace);
  return result;
}

static int check_terminator_cases(void) {
  for (size_t cells = 1; cells <= 4; ++cells)
    if (check_terminator_solver(cells, 0) ||
        check_terminator_solver(cells, 1)) return 1;
  return 0;
}

static int check_chapman_solver_comparisons(const chem_mechanism *mechanism) {
  const double temperatures[]={284.19,215.02,299.31};
  const double pressures[]={101245.0,100789.2,101398.0};
  const double photons[3][3]={{1e-4,1e-5,1e-6},
      {3.2e-4,7.3e-5,3.2e-6},{5.2e-4,8.2e-5,4.6e-6}};
  const double initial[3][9]={
      {0.8,0.2,8.84,0.1,1.0,0.392,32.9,0.01,0.0169},
      {1.2,0.4,7.2,0.3,1.0,0.39,30.0,0.02,0.014},
      {0.4,0.1,9.3,0.5,1.0,0.4,34.0,0.01,0.02}};
  const double gas_constant=8.31446261815324;
  const chem_solver_options options={.relative_tolerance=1e-6,
      .absolute_tolerance=1e-12,.max_steps=100000};
  double reference[3][9];
  const chem_rosenbrock_parameters baseline=
      chem_rosenbrock_builtin_parameters(CHEM_ROSENBROCK_THREE_STAGE);
  chem_solver_workspace baseline_workspace={0};
  if (chem_solver_workspace_init_group(&baseline_workspace,mechanism,3)) return 1;
  for (size_t cell=0;cell<3;++cell) {
    memcpy(reference[cell],initial[cell],sizeof(reference[cell]));
    const double density=pressures[cell]/(gas_constant*temperatures[cell]);
    chem_solver_stats stats={0};
    if (chem_rosenbrock_solve_with_parameters(mechanism,&baseline_workspace,
          reference[cell],temperatures[cell],pressures[cell],density,photons[cell],
          500,&baseline,&options,&stats) || stats.accepted==0) {
      chem_solver_workspace_destroy(&baseline_workspace); return 2;
    }
  }
  chem_solver_workspace_destroy(&baseline_workspace);

  const chem_rosenbrock_parameter_set schemes[]={CHEM_ROSENBROCK_TWO_STAGE,
      CHEM_ROSENBROCK_THREE_STAGE,CHEM_ROSENBROCK_FOUR_STAGE,
      CHEM_ROSENBROCK_FOUR_STAGE_DAE,CHEM_ROSENBROCK_SIX_STAGE_DAE};
  const double tolerances[]={1e-2,1e-4,1e-4,1e-4,1e-4};
  for (size_t scheme=0;scheme<5;++scheme) {
    chem_solver_workspace workspace={0};
    if (chem_solver_workspace_init_group(&workspace,mechanism,3)) return 3;
    const chem_rosenbrock_parameters parameters=
        chem_rosenbrock_builtin_parameters(schemes[scheme]);
    for (size_t cell=0;cell<3;++cell) {
      double y[9]; memcpy(y,initial[cell],sizeof(y));
      const double density=pressures[cell]/(gas_constant*temperatures[cell]);
      chem_solver_stats stats={0};
      if (chem_rosenbrock_solve_with_parameters(mechanism,&workspace,y,
            temperatures[cell],pressures[cell],density,photons[cell],500,
            &parameters,&options,&stats) || stats.accepted==0) {
        chem_solver_workspace_destroy(&workspace); return 4;
      }
      const double abs_tol=1e-12;
      for (size_t variable=0;variable<9;++variable)
        if (fabs(y[variable]-reference[cell][variable])>
            tolerances[scheme]*(fabs(y[variable])+fabs(reference[cell][variable]))+
            abs_tol) {
          chem_solver_workspace_destroy(&workspace); return 5;
        }
    }
    chem_solver_workspace_destroy(&workspace);
  }

  /* VectorSolve instantiates the same solver at lane widths one through four. */
  for (size_t width=1;width<=4;++width) {
    chem_solver_workspace workspace={0};
    if (chem_solver_workspace_init_group(&workspace,mechanism,width)) return 6;
    for (size_t cell=0;cell<3;++cell) {
      double y[9]; memcpy(y,initial[cell],sizeof(y));
      const double density=pressures[cell]/(gas_constant*temperatures[cell]);
      if (chem_rosenbrock3_solve_with_workspace(mechanism,&workspace,y,
            temperatures[cell],pressures[cell],density,photons[cell],500,
            &options,NULL)) {
        chem_solver_workspace_destroy(&workspace); return 7;
      }
      for (size_t variable=0;variable<9;++variable)
        if (fabs(y[variable]-reference[cell][variable])>
            1e-4*(fabs(y[variable])+fabs(reference[cell][variable]))+1e-12) {
          chem_solver_workspace_destroy(&workspace); return 8;
        }
    }
    chem_solver_workspace_destroy(&workspace);
  }
  return 0;
}

static int check_chapman_integration(void) {
  const chem_species species[] = {
      {.name = "O"}, {.name = "O1D"}, {.name = "O2"}, {.name = "O3"},
      {.name = "M"}, {.name = "Ar"}, {.name = "N2"}, {.name = "H2O"},
      {.name = "CO2"}};
  const chem_term o1d_n2[] = {{.species=1,.coefficient=1}, {.species=6,.coefficient=1}};
  const chem_term o1d_o2[] = {{.species=1,.coefficient=1}, {.species=2,.coefficient=1}};
  const chem_term o_o3[] = {{.species=0,.coefficient=1}, {.species=3,.coefficient=1}};
  const chem_term o_o2_m[] = {{.species=0,.coefficient=1}, {.species=2,.coefficient=1}, {.species=4,.coefficient=1}};
  const chem_term o2[] = {{.species=2,.coefficient=1}};
  const chem_term o3[] = {{.species=3,.coefficient=1}};
  const chem_term o_n2[] = {{.species=0,.coefficient=1}, {.species=6,.coefficient=1}};
  const chem_term o_o2[] = {{.species=0,.coefficient=1}, {.species=2,.coefficient=1}};
  const chem_term two_o2[] = {{.species=2,.coefficient=2}};
  const chem_term o3_m[] = {{.species=3,.coefficient=1}, {.species=4,.coefficient=1}};
  const chem_term two_o[] = {{.species=0,.coefficient=2}};
  const chem_term o1d_o2_product[] = {{.species=1,.coefficient=1}, {.species=2,.coefficient=1}};
  const chem_term o_o2_product[] = {{.species=0,.coefficient=1}, {.species=2,.coefficient=1}};
  const chem_reaction reactions[] = {
      {.kind=CHEM_ARRHENIUS,.reactants=o1d_n2,.reactant_count=2,.products=o_n2,.product_count=2,.p={2.15e-11,0,110,300,0}},
      {.kind=CHEM_ARRHENIUS,.reactants=o1d_o2,.reactant_count=2,.products=o_o2,.product_count=2,.p={3.3e-11,0,55,300,0}},
      {.kind=CHEM_ARRHENIUS,.reactants=o_o3,.reactant_count=2,.products=two_o2,.product_count=1,.p={8e-12,0,-2060,300,0}},
      {.kind=CHEM_ARRHENIUS,.reactants=o_o2_m,.reactant_count=3,.products=o3_m,.product_count=2,.p={6e-34,2.4,0,300,0}},
      {.kind=CHEM_USER_DEFINED,.label="jO2",.reactants=o2,.reactant_count=1,.products=two_o,.product_count=1,.p={1,0,0,0,0,0,0,0}},
      {.kind=CHEM_USER_DEFINED,.label="jO3a",.reactants=o3,.reactant_count=1,.products=o1d_o2_product,.product_count=2,.p={1,0,0,0,0,0,0,1}},
      {.kind=CHEM_USER_DEFINED,.label="jO3b",.reactants=o3,.reactant_count=1,.products=o_o2_product,.product_count=2,.p={1,0,0,0,0,0,0,2}}};
  const chem_mechanism mechanism = {
      .species=species,.species_count=9,.reactions=reactions,.reaction_count=7,
      .custom_labels=(const char *const[]){"jO2","jO3a","jO3b"},
      .custom_parameter_count=3};
  const char *expected_species_names[]={"O","O1D","O2","O3","M",
                                       "Ar","N2","H2O","CO2"};
  const chem_rosenbrock_parameters default_parameters=
      chem_rosenbrock_builtin_parameters(CHEM_ROSENBROCK_THREE_STAGE);
  if (default_parameters.stages!=3 || mechanism.species_count!=9 ||
      mechanism.reaction_count!=7 || mechanism.custom_parameter_count!=3)
    return 12;
  for (size_t i=0;i<mechanism.species_count;++i)
    if (strcmp(species[i].name,expected_species_names[i])!=0) return 13;
  double unit_state[9],unit_rates[7],unit_forcing[9];
  for (size_t i=0;i<9;++i) unit_state[i]=1.0;
  for (size_t i=0;i<7;++i) unit_rates[i]=1.0;
  if (chem_forcing_with_constants(&mechanism,unit_state,unit_rates,0.0,
                                  unit_forcing)) return 14;
  const double expected_unit_forcing[]={3,-1,2,-2,0,0,0,0,0};
  for (size_t i=0;i<9;++i)
    if (unit_forcing[i]!=expected_unit_forcing[i]) return 15;
  double small_state[9],small_rates[7],small_forcing[9];
  for (size_t i=0;i<9;++i) small_state[i]=5e-6;
  for (size_t i=0;i<7;++i) small_rates[i]=3e-8;
  if (chem_forcing_with_constants(&mechanism,small_state,small_rates,6e-14,
                                  small_forcing)) return 16;
  const double small_o1d=small_rates[5]*small_state[3]-small_rates[0]*
      small_state[6]*small_state[1]-small_rates[1]*small_state[1]*small_state[2];
  const double small_o=2*small_rates[4]*small_state[2]+small_rates[6]*small_state[3]+
      small_rates[0]*small_state[6]*small_state[1]+small_rates[1]*small_state[1]*
      small_state[2]-small_rates[2]*small_state[0]*small_state[3]-
      small_rates[3]*small_state[4]*small_state[0]*small_state[2];
  const double small_o2=-small_rates[4]*small_state[2]+small_rates[5]*small_state[3]+
      small_rates[6]*small_state[3]+2*small_rates[2]*small_state[0]*small_state[3]-
      small_rates[3]*small_state[4]*small_state[0]*small_state[2];
  const double small_o3=-small_rates[5]*small_state[3]-small_rates[6]*small_state[3]-
      small_rates[2]*small_state[0]*small_state[3]+small_rates[3]*small_state[4]*
      small_state[0]*small_state[2];
  const double expected_small_forcing[]={small_o,small_o1d,small_o2,small_o3,0,0,0,0,0};
  for (size_t i=0;i<9;++i)
    if (fabs(small_forcing[i]-expected_small_forcing[i])>
        1e-12*fmax(1e-30,fabs(expected_small_forcing[i]))) return 17;
  double negative_jacobian[81];
  if (chem_negative_jacobian_with_constants(&mechanism,small_state,small_rates,
                                             negative_jacobian)) return 18;
  const double direction[]={0.3,-0.2,0.4,-0.1,0.5,-0.6,0.7,0.8,-0.9};
  double plus[9],minus[9],fplus[9],fminus[9];
  for (size_t i=0;i<9;++i) {
    plus[i]=small_state[i]+1e-10*direction[i];
    minus[i]=small_state[i]-1e-10*direction[i];
  }
  if (chem_forcing_with_constants(&mechanism,plus,small_rates,6e-14,fplus) ||
      chem_forcing_with_constants(&mechanism,minus,small_rates,6e-14,fminus))
    return 19;
  for (size_t row=0;row<9;++row) {
    double product=0;
    for (size_t col=0;col<9;++col)
      product-=negative_jacobian[row*9+col]*direction[col];
    const double finite_difference=(fplus[row]-fminus[row])/(2e-10);
    if (fabs(product-finite_difference)>1e-10*fmax(1e-20,fabs(finite_difference)))
      return 20;
  }
  const double cell_temperatures[] = {284.19,215.02,299.31};
  const double cell_pressures[] = {101245.0,100789.2,101398.0};
  const double cell_photolysis[3][3] = {
      {1.0e-4,1.0e-5,1.0e-6},
      {3.2e-4,7.3e-5,3.2e-6},
      {5.2e-4,8.2e-5,4.6e-6}};
  const double gas_constant = 8.31446261815324;
  for (size_t width=1;width<=4;++width) {
    double packed[4*9]={0};
    chem_vector_matrix_view packed_view;
    if (chem_vector_matrix_view_init(&packed_view,packed,3,9,width)) return 6;
    for (size_t cell=0;cell<3;++cell)
      for (size_t col=0;col<9;++col)
        if (chem_vector_matrix_view_set(&packed_view,cell,col,
              1.0+0.1*(double)(9*cell+col))) return 7;
    for (size_t cell=0;cell<3;++cell) {
      const double temperature=cell_temperatures[cell];
      const double pressure=cell_pressures[cell];
      const double density=pressure/(gas_constant*temperature);
      double custom[3],rate[7],state_values[9],forcing[9];
      memcpy(custom,cell_photolysis[cell],sizeof(custom));
      if (chem_vector_matrix_view_copy_row_to_vector(
              &packed_view,cell,state_values,9) ||
          chem_rate_constants(&mechanism,temperature,pressure,density,custom,rate))
        return 8;
      const double expected_rates[] = {
          2.15e-11*exp(110.0/temperature),
          3.3e-11*exp(55.0/temperature),
          8e-12*exp(-2060.0/temperature),
          6e-34*pow(temperature/300.0,2.4),
          custom[0],custom[1],custom[2]};
      for (size_t r=0;r<7;++r)
        if (fabs(rate[r]-expected_rates[r])>1e-12*fabs(expected_rates[r])+1e-30)
          return 9;
      if (chem_forcing_with_constants(&mechanism,state_values,rate,density,forcing))
        return 10;
      const double o1d=rate[5]*state_values[3]-rate[0]*state_values[6]*state_values[1]-
          rate[1]*state_values[1]*state_values[2];
      const double oxygen=2.0*rate[4]*state_values[2]+rate[6]*state_values[3]+
          rate[0]*state_values[6]*state_values[1]+rate[1]*state_values[1]*state_values[2]-
          rate[2]*state_values[0]*state_values[3]-rate[3]*state_values[4]*
          state_values[0]*state_values[2];
      const double oxygen2=-rate[4]*state_values[2]+rate[5]*state_values[3]+
          rate[6]*state_values[3]+2.0*rate[2]*state_values[0]*state_values[3]-
          rate[3]*state_values[4]*state_values[0]*state_values[2];
      const double ozone=-rate[5]*state_values[3]-rate[6]*state_values[3]-
          rate[2]*state_values[0]*state_values[3]+rate[3]*state_values[4]*
          state_values[0]*state_values[2];
      const double expected_forcing[9] = {oxygen,o1d,oxygen2,ozone,0,0,0,0,0};
      for (size_t i=0;i<9;++i)
        if (fabs(forcing[i]-expected_forcing[i])>1e-12*fabs(expected_forcing[i])+1e-30)
          return 11;
    }
  }
  double state[] = {0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3};
  double initial[9];
  memcpy(initial, state, sizeof(state));
  const double photo_rates[] = {0.1,0.2,0.3};
  const double oxygen_initial = state[0] + state[1] + 2.0*state[2] + 3.0*state[3];
  chem_solver_workspace workspace = {0};
  if (chem_solver_workspace_init(&workspace, &mechanism)) return 1;
  chem_solver_options options = chem_solver_default_options();
  options.max_steps = 100000;
  const double chapman_initial[]={0,0,8.84,0,1,0.392,32.9,0,0.0169};
  const double chapman_expected[]={0.00176789,1.8039e-41,8.83912,
      4.5031e-36,1,0.392,32.9,0,0.0169};
  const double chapman_scales[]={1.0,0.1,10.0};
  const double chapman_photos[]={1e-4,1e-5,1e-6};
  for (size_t sample=0;sample<3;++sample) {
    double one_step[9];
    for (size_t i=0;i<9;++i)
      one_step[i]=chapman_initial[i]*chapman_scales[sample];
    const int one_step_status=chem_rosenbrock3_solve_with_workspace(
        &mechanism,&workspace,one_step,273.15,100000.0,2.7e19,
        chapman_photos,1.0,&options,NULL);
    if (one_step_status) {
      chem_solver_workspace_destroy(&workspace);
      return 21;
    }
    for (size_t i=0;i<9;++i)
      if (fabs(one_step[i]-chapman_expected[i]*chapman_scales[sample])>1e-4) {
        chem_solver_workspace_destroy(&workspace);
        return 22;
      }
  }
  double random_densities[9];
  for (size_t i=0;i<9;++i)
    random_densities[i]=((double)((i*37+11)%101)+0.5)/101.0;
  if (chem_rosenbrock3_solve_with_workspace(&mechanism,&workspace,
          random_densities,273.15,100000.0,2.7e19,chapman_photos,1.0,
          &options,NULL)) {
    chem_solver_workspace_destroy(&workspace);
    return 23;
  }
  for (size_t i=0;i<9;++i)
    if (!isfinite(random_densities[i]) || random_densities[i]<0.0) {
      chem_solver_workspace_destroy(&workspace);
      return 24;
    }
  const double temperature = 2.0, pressure = 3.0;
  const double density = pressure/(8.31446261815324*temperature);
  int result = 0;
  for (size_t step = 0; step < 100 && !result; ++step) {
    result = chem_rosenbrock3_solve_with_workspace(
        &mechanism,&workspace,state,temperature,pressure,density,photo_rates,
        30.0,&options,NULL);
    if (!result) {
      for (size_t i=0;i<9;++i)
        if (!isfinite(state[i]) || state[i]<0.0) result=2;
      const double oxygen = state[0]+state[1]+2.0*state[2]+3.0*state[3];
      if (fabs(oxygen-oxygen_initial)>1.0e-8*oxygen_initial) result=3;
      if (state[4]!=initial[4] || state[5]!=initial[5] ||
          state[6]!=initial[6] || state[7]!=initial[7] ||
          state[8]!=initial[8]) result=4;
    }
  }
  chem_solver_workspace_destroy(&workspace);
  if (result) return result;
  if (!(state[2] < 0.5*initial[2] && state[0] > initial[0] &&
        state[0] > 0.5*oxygen_initial)) return 5;
  return check_chapman_solver_comparisons(&mechanism);
}

static int check_integrated_reaction_rates(void) {
  const chem_species species[] = {{.name="a"}, {.name="b"}, {.name="c"},
                                  {.name="irr_1"}, {.name="irr_2"}};
  const chem_term a[] = {{.species=0,.coefficient=1.0}};
  const chem_term b[] = {{.species=1,.coefficient=1.0}};
  const chem_term b_irr1[] = {{.species=1,.coefficient=1.0},
                              {.species=3,.coefficient=1.0}};
  const chem_term c_irr2[] = {{.species=2,.coefficient=1.0},
                              {.species=4,.coefficient=1.0}};
  const chem_reaction reactions[] = {
      {.kind=CHEM_ARRHENIUS,.reactants=a,.reactant_count=1,
       .products=b_irr1,.product_count=2,.p={2.15e-11,0,110,300,0}},
      {.kind=CHEM_USER_DEFINED,.label="r2",.reactants=b,.reactant_count=1,
       .products=c_irr2,.product_count=2,.p={1,0,0,0,0,0,0,0}}};
  const chem_mechanism mechanism = {.species=species,.species_count=5,
      .reactions=reactions,.reaction_count=2,.custom_parameter_count=1};
  double state[] = {1.0,0,0,0,0};
  chem_solver_workspace workspace = {0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  const double temperature=273.0, pressure=1000.0;
  const double density=pressure/(8.31446261815324*temperature);
  double previous_irr1=0.0, previous_irr2=0.0;
  int result=0;
  for (size_t step=0;step<100 && !result;++step) {
    const double rate=(step>50)?0.0:1.0;
    result=chem_rosenbrock3_solve_with_workspace(
        &mechanism,&workspace,state,temperature,pressure,density,&rate,30.0,
        NULL,NULL);
    if (!result && (state[3]<previous_irr1 || state[4]<previous_irr2)) result=2;
    previous_irr1=state[3];
    previous_irr2=state[4];
  }
  chem_solver_workspace_destroy(&workspace);
  return result;
}

static int check_flow_tube_reference(int backward_euler) {
  static const double expected[2][10][4] = {
      {{0,0,8e-8,2e-5},
       {1.51776e-11,7.43233e-12,7.99153e-8,1.99999e-5},
       {3.01881e-11,1.44946e-11,7.98306e-8,1.99998e-5},
       {4.50332e-11,2.12049e-11,7.97460e-8,1.99997e-5},
       {5.97144e-11,2.75804e-11,7.96616e-8,1.99997e-5},
       {7.42334e-11,3.36373e-11,7.95772e-8,1.99996e-5},
       {8.85919e-11,3.93913e-11,7.94929e-8,1.99995e-5},
       {1.02791e-10,4.48570e-11,7.94087e-8,1.99994e-5},
       {1.16834e-10,5.00486e-11,7.93245e-8,1.99993e-5},
       {1.30720e-10,5.49794e-11,7.92405e-8,1.99992e-5}},
      {{0,0,8e-8,2e-5},
       {1.50953e-11,7.26014e-12,7.99153e-8,1.99999e-5},
       {3.00252e-11,1.41668e-11,7.98307e-8,1.99998e-5},
       {4.47911e-11,2.07369e-11,7.97462e-8,1.99997e-5},
       {5.93949e-11,2.69865e-11,7.96617e-8,1.99997e-5},
       {7.38381e-11,3.29307e-11,7.95774e-8,1.99996e-5},
       {8.81223e-11,3.85842e-11,7.94931e-8,1.99995e-5},
       {1.02249e-10,4.39608e-11,7.94090e-8,1.99994e-5},
       {1.16220e-10,4.90738e-11,7.93249e-8,1.99993e-5},
       {1.30037e-10,5.39355e-11,7.92409e-8,1.99992e-5}}};
  const chem_species species[] = {{.name="SOA1"}, {.name="SOA2"},
      {.name="a-pinene"}, {.name="O3"}};
  const chem_term soa1[] = {{.species=0,.coefficient=1.0}};
  const chem_term soa2[] = {{.species=1,.coefficient=1.0}};
  const chem_term oxidation[] = {{.species=2,.coefficient=1.0},
                                  {.species=3,.coefficient=1.0}};
  const chem_term soa_products[] = {{.species=0,.coefficient=0.18},
                                     {.species=1,.coefficient=0.09}};
  const chem_reaction reactions[] = {
      {.kind=CHEM_USER_DEFINED,.label="r1",.reactants=soa1,.reactant_count=1,
       .p={1,0,0,0,0,0,0,0}},
      {.kind=CHEM_USER_DEFINED,.label="r2",.reactants=soa2,.reactant_count=1,
       .p={1,0,0,0,0,0,0,0}},
      {.kind=CHEM_ARRHENIUS,.reactants=oxidation,.reactant_count=2,
       .products=soa_products,.product_count=2,
       .p={8.8e-17*6.02214076e17,0,0,298.15,0}}};
  const char *labels[] = {"r1","r2"};
  const chem_mechanism mechanism = {.species=species,.species_count=4,
      .reactions=reactions,.reaction_count=3,.custom_labels=labels,
      .custom_parameter_count=2};
  const double rates[] = {0.01,0.05};
  const double temperature=298.15, pressure=101325.0;
  const double density=pressure/(8.31446261815324*temperature);
  double state[] = {0.0,0.0,8e-8,2e-5};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  int status=0;
  for (size_t step=1;step<10 && !status;++step) {
    if (backward_euler) {
      const chem_backward_euler_options options=chem_backward_euler_default_options();
      status=chem_backward_euler_solve_with_workspace(
          &mechanism,&workspace,state,temperature,pressure,density,rates,1.0,
          &options,NULL);
    } else {
      const chem_solver_options options=chem_solver_default_options();
      status=chem_rosenbrock3_solve_with_workspace(
          &mechanism,&workspace,state,temperature,pressure,density,rates,1.0,
          &options,NULL);
    }
    for (size_t i=0;i<4 && !status;++i)
      if (fabs(state[i]-expected[backward_euler][step][i])>1e-10) status=2;
  }
  chem_solver_workspace_destroy(&workspace);
  return status;
}

static double stub_temperature_rate(double temperature, double pressure,
                                   double air_density,
                                   const double *parameters) {
  (void)pressure; (void)air_density; (void)parameters;
  return 0.005*temperature;
}

static int check_stub_aerosol_mechanism(int backward_euler) {
  const char *names[] = {
      "FO2", "BAR", "STUB1.MODE1.QUUX.BAZ", "STUB1.MODE1.QUUX.QUX",
      "STUB1.MODE2.QUUX.BAZ", "STUB1.MODE2.QUUX.QUX",
      "STUB1.MODE2.CORGE.FO2", "STUB1.MODE2.CORGE.BAZ",
      "STUB1.MODE2.CORGE.QUX", "STUB2.MODE1.NUMBER",
      "STUB2.MODE1.QUUX.BAZ", "STUB2.MODE1.QUUX.QUX",
      "STUB2.MODE2.NUMBER", "STUB2.MODE2.CORGE.FO2",
      "STUB2.MODE2.CORGE.BAZ", "STUB2.MODE2.CORGE.QUX",
      "STUB2.MODE3.NUMBER", "STUB2.MODE3.QUUX.BAZ",
      "STUB2.MODE3.QUUX.QUX", "STUB2.MODE3.CORGE.FO2",
      "STUB2.MODE3.CORGE.BAZ", "STUB2.MODE3.CORGE.QUX"};
  chem_species species[22];
  for (size_t i=0;i<22;++i) species[i]=(chem_species){.name=names[i]};
  const chem_term fo2[]={ {0,1.0} }, fo2_dest[]={ {6,1.0} };
  const chem_term baz_source[]={ {2,1.0} }, baz_dest[]={ {4,1.0} };
  const chem_term stub2_fo2[]={ {13,1.0} }, stub2_baz[]={ {14,1.0} };
  const chem_term mode3_baz[]={ {17,1.0} }, mode3_qux[]={ {18,1.0} };
  const char *custom_labels[]={"STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT"};
  const chem_reaction reactions[]={
      {.kind=CHEM_ARRHENIUS,.reactants=fo2,.reactant_count=1,
       .products=fo2_dest,.product_count=1,.p={0.001,0,0,0,0}},
      {.kind=CHEM_ARRHENIUS,.reactants=baz_source,.reactant_count=1,
       .products=baz_dest,.product_count=1,.p={0.002,0,0,0,0}},
      {.kind=CHEM_USER_DEFINED,.label=custom_labels[0],.reactants=stub2_fo2,
       .reactant_count=1,.products=stub2_baz,.product_count=1,
       .p={1,0,0,0,0,0,0,0}},
      {.kind=CHEM_LAMBDA,.reactants=mode3_baz,.reactant_count=1,
       .products=mode3_qux,.product_count=1,.rate_function=stub_temperature_rate}};
  const chem_mechanism mechanism={.species=species,.species_count=22,
      .reactions=reactions,.reaction_count=4,.custom_labels=custom_labels,
      .custom_parameter_count=1};
  double state_values[3*22]={0};
  double custom_values[3]={0};
  chem_state_view state;
  if (chem_state_view_init(&state,names,22,custom_labels,1,3,state_values,
                           custom_values)) return 1;
  const double fo2_values[]={1.0,1.5,2.0};
  const double bar_values[]={2.0,2.5,3.0};
  const double baz_values[]={0.5,0.7,0.9};
  const double mode2_fo2[]={0.4,0.6,0.8};
  const chem_state_series series[]={
      {"FO2",fo2_values,3},{"BAR",bar_values,3},
      {"STUB1.MODE1.QUUX.BAZ",baz_values,3},
      {"STUB2.MODE2.CORGE.FO2",mode2_fo2,3},
      {"STUB2.MODE3.QUUX.BAZ",(const double[]){0.2,0.3,0.4},3}};
  if (chem_state_set_concentrations(&state,series,5)) return 2;
  const double rates[]={0.01};
  double forcing[22], jacobian[22*22];
  if (chem_forcing(&mechanism,state_values,275.0,1e5,rates,forcing) ||
      chem_jacobian(&mechanism,state_values,275.0,1e5,rates,jacobian)) return 3;
  if (fabs(forcing[0]+0.001)>1e-14 || fabs(forcing[6]-0.001)>1e-14 ||
      fabs(forcing[2]+0.001)>1e-14 || fabs(forcing[4]-0.001)>1e-14 ||
      fabs(forcing[13]+0.004)>1e-14 || fabs(forcing[14]-0.004)>1e-14 ||
      fabs(forcing[17]+0.275)>1e-14 || fabs(forcing[18]-0.275)>1e-14 ||
      fabs(jacobian[0*22+0]+0.001)>1e-14 ||
      fabs(jacobian[6*22+0]-0.001)>1e-14 ||
      fabs(jacobian[2*22+2]+0.002)>1e-14 ||
      fabs(jacobian[4*22+2]-0.002)>1e-14) return 4;
  const double number_values[]={1000.0,1100.0,1200.0};
  if (chem_state_set_concentration(&state,"STUB2.MODE1.NUMBER",
                                   number_values,3)) return 5;
  for (size_t i=0;i<22;++i) {
    double observed=0.0;
    if (chem_state_get_concentration(&state,names[i],0,&observed)) return 6;
  }
  /* Solve three cells independently; unused aerosol state stays unchanged. */
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 7;
  const double temperatures[]={275.0,285.0,295.0};
  const double stub2_rates[]={0.01,0.015,0.02};
  int status=0;
  for (size_t cell=0;cell<3 && !status;++cell) {
    double y[22]={0};
    y[0]=fo2_values[cell]; y[1]=bar_values[cell]; y[2]=baz_values[cell];
    y[4]=0.0; y[6]=0.8+0.2*(double)cell;
    y[13]=mode2_fo2[cell]; y[14]=0.1+0.05*(double)cell;
    y[17]=0.2+0.1*(double)cell; y[18]=0.3+0.1*(double)cell;
    const double cell_rates[]={stub2_rates[cell]};
    if (backward_euler) {
      const chem_backward_euler_options options=chem_backward_euler_default_options();
      status=chem_backward_euler_solve_with_workspace(&mechanism,&workspace,y,
          temperatures[cell],101325.0,1.0e5,cell_rates,10.0,&options,NULL);
    } else {
      const chem_solver_options options=chem_solver_default_options();
      status=chem_rosenbrock3_solve_with_workspace(&mechanism,&workspace,y,
          temperatures[cell],101325.0,1.0e5,cell_rates,10.0,&options,NULL);
    }
    const double dfo2=fo2_values[cell]*(1.0-exp(-0.001*10.0));
    const double dbaz=baz_values[cell]*(1.0-exp(-0.002*10.0));
    const double dstub2=mode2_fo2[cell]*(1.0-exp(-stub2_rates[cell]*10.0));
    const double dmode3=(0.2+0.1*(double)cell)*
        (1.0-exp(-temperatures[cell]*0.005*10.0));
    const double tol=backward_euler?1.0e-1:5.0e-5;
    if (!status && (fabs(y[0]-(fo2_values[cell]-dfo2))>tol*fo2_values[cell] ||
        y[1]!=bar_values[cell] ||
        fabs(y[2]-(baz_values[cell]-dbaz))>tol*baz_values[cell] ||
        fabs(y[4]-dbaz)>tol*baz_values[cell] ||
        fabs(y[6]-(0.8+0.2*(double)cell+dfo2))>tol*fo2_values[cell] ||
        fabs(y[13]-(mode2_fo2[cell]-dstub2))>tol*mode2_fo2[cell] ||
        fabs(y[14]-(0.1+0.05*(double)cell+dstub2))>tol*mode2_fo2[cell] ||
        fabs(y[17]-(0.2+0.1*(double)cell-dmode3))>tol*(0.2+0.1*(double)cell) ||
        fabs(y[18]-(0.3+0.1*(double)cell+dmode3))>tol*(0.2+0.1*(double)cell))) status=8;
  }
  if (backward_euler && !status) {
    const double zero_rate[]={0.0};
    const double cell_fo2[]={1.0,1.5,2.0};
    const double cell_baz[]={0.5,0.7,0.9};
    for (size_t cell=0;cell<3 && !status;++cell) {
      double y[22]={0};
      y[0]=cell_fo2[cell]; y[1]=2.0+0.5*(double)cell;
      y[2]=cell_baz[cell]; y[6]=0.8+0.2*(double)cell;
      const chem_backward_euler_options options=chem_backward_euler_default_options();
      status=chem_backward_euler_solve_with_workspace(&mechanism,&workspace,y,
          298.15,101325.0,1.0e5,zero_rate,10.0,&options,NULL);
      const double dfo2=cell_fo2[cell]*(1.0-exp(-0.001*10.0));
      const double dbaz=cell_baz[cell]*(1.0-exp(-0.002*10.0));
      if (!status && (fabs(y[0]-(cell_fo2[cell]-dfo2))>2e-4*cell_fo2[cell] ||
          y[1]!=2.0+0.5*(double)cell ||
          fabs(y[2]-(cell_baz[cell]-dbaz))>2e-4*cell_baz[cell] ||
          fabs(y[4]-dbaz)>2e-4*cell_baz[cell] ||
          fabs(y[6]-(0.8+0.2*(double)cell+dfo2))>2e-4*cell_fo2[cell])) status=9;
    }
  }
  chem_solver_workspace_destroy(&workspace);
  return status;
}

static int check_lambda_callback_override(void) {
  const chem_species species[] = {{.name="A"}, {.name="B"}};
  const chem_term a[] = {{.species=0,.coefficient=1.0}};
  const chem_term b[] = {{.species=1,.coefficient=1.0}};
  chem_reaction reaction = {.kind=CHEM_LAMBDA,.label="lambda_rc",
      .reactants=a,.reactant_count=1,.products=b,.product_count=1,
      .rate_function=test_lambda};
  const chem_mechanism mechanism = {.species=species,.species_count=2,
      .reactions=&reaction,.reaction_count=1};
  reaction.rate_function = test_lambda_override;
  double rate = 0.0;
  if (chem_rate_constants(&mechanism,300.0,101253.3,1.0e6,NULL,&rate) ||
      !near(rate,0.6)) return 1;
  double state[] = {1.0,0.0};
  chem_solver_workspace workspace = {0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 2;
  const int status = chem_rosenbrock3_solve_with_workspace(
      &mechanism,&workspace,state,300.0,101253.3,1.0e6,NULL,1.0,NULL,NULL);
  chem_solver_workspace_destroy(&workspace);
  return status == 0 && state[0] < 1.0 && state[1] > 0.0 ? 0 : 3;
}

static int check_phase_species_metadata(void) {
  chem_species without_diffusion = {.name="CO2"};
  if (strcmp(without_diffusion.name,"CO2") != 0 ||
      without_diffusion.diffusion != 0.0) return 1;
  chem_species with_diffusion = {.name="CO2",.diffusion=1.23e-5};
  if (strcmp(with_diffusion.name,"CO2") != 0 ||
      !near(with_diffusion.diffusion,1.23e-5)) return 2;
  without_diffusion.diffusion = 2.5e-6;
  if (!near(without_diffusion.diffusion,2.5e-6)) return 3;
  const chem_species third_body = {.name="M",.third_body=1};
  if (!chem_species_is_parameterized(&third_body) ||
      !near(chem_species_parameterize(&third_body,42.4),42.4)) return 4;
  return 0;
}

static int check_constraint_api(void) {
  const chem_constraint_term ab_terms[] = {
      {.species=0,.coefficient=1.0}, {.species=1,.coefficient=1.0}};
  const chem_constraint linear = {.kind=CHEM_CONSTRAINT_LINEAR,
      .name="A_B_conservation",.algebraic_species=1,.terms=ab_terms,
      .term_count=2,.constant=1.0};
  const chem_constraint_set linear_set = {
      .constraints=&linear,.constraint_count=1,.species_count=2};
  if (chem_constraint_set_validate(&linear_set)) return 1;
  unsigned char mask[4] = {0}, pattern[16] = {0};
  if (chem_constraint_set_algebraic_mask(&linear_set,mask) ||
      mask[0] || mask[1] != 1 ||
      chem_constraint_set_jacobian_pattern(&linear_set,pattern) ||
      pattern[2] != 1 || pattern[3] != 1 || pattern[0] || pattern[1]) return 2;
  const double states[] = {0.3,0.7, 0.5,0.6};
  const double temperatures[] = {298.15,298.15};
  double forcing[4] = {0};
  if (chem_constraint_set_add_forcing(&linear_set,states,temperatures,2,forcing) ||
      !near(forcing[1],0.0) || !near(forcing[3],0.1)) return 3;
  double jacobian[8] = {0};
  if (chem_constraint_set_add_negative_jacobian(
          &linear_set,states,temperatures,2,jacobian) ||
      !near(jacobian[2],-1.0) || !near(jacobian[3],-1.0) ||
      !near(jacobian[6],-1.0) || !near(jacobian[7],-1.0)) return 4;

  const chem_constraint_term weighted_terms[] = {
      {.species=0,.coefficient=2.0}, {.species=1,.coefficient=3.0},
      {.species=2,.coefficient=-1.0}};
  const chem_constraint weighted = {.kind=CHEM_CONSTRAINT_LINEAR,
      .name="weighted_sum",.algebraic_species=2,.terms=weighted_terms,
      .term_count=3,.constant=5.0};
  if (strcmp(weighted.name,"weighted_sum") != 0 ||
      weighted.algebraic_species != 2 || weighted.term_count != 3 ||
      weighted.constant != 5.0 || weighted.terms[0].coefficient != 2.0 ||
      weighted.terms[1].coefficient != 3.0 ||
      weighted.terms[2].coefficient != -1.0) return 51;
  const chem_constraint_set weighted_set = {
      .constraints=&weighted,.constraint_count=1,.species_count=3};
  const double weighted_state[] = {1.0,2.0,3.0};
  double weighted_forcing[3] = {0};
  double weighted_jacobian[9] = {0};
  if (chem_constraint_set_add_forcing(&weighted_set,weighted_state,NULL,1,
                                      weighted_forcing) ||
      !near(weighted_forcing[2],0.0) ||
      chem_constraint_set_add_negative_jacobian(
          &weighted_set,weighted_state,NULL,1,weighted_jacobian) ||
      !near(weighted_jacobian[6],-2.0) ||
      !near(weighted_jacobian[7],-3.0) ||
      !near(weighted_jacobian[8],1.0)) return 5;
  const double weighted_state_2[] = {9.0,0.2,8.0};
  double weighted_jacobian_2[9]={0};
  if (chem_constraint_set_add_negative_jacobian(
          &weighted_set,weighted_state_2,NULL,1,weighted_jacobian_2)) return 52;
  for (size_t i=0;i<9;++i)
    if (weighted_jacobian[i] != weighted_jacobian_2[i]) return 53;

  const chem_constraint_term conservation_terms[] = {
      {.species=0,.coefficient=1.0}, {.species=1,.coefficient=1.0},
      {.species=2,.coefficient=1.0}};
  const chem_constraint conservation = {.kind=CHEM_CONSTRAINT_LINEAR,
      .name="conservation",.algebraic_species=2,.terms=conservation_terms,
      .term_count=3,.constant=1.0};
  const chem_constraint_set conservation_set = {.constraints=&conservation,
      .constraint_count=1,.species_count=3};
  const double conservation_state[] = {0.2,0.3,0.5};
  double conservation_forcing[3]={0};
  if (chem_constraint_set_add_forcing(&conservation_set,conservation_state,
                                      NULL,1,conservation_forcing) ||
      fabs(conservation_forcing[2])>1.0e-14) return 54;
  const double conservation_other[] = {0.4,0.3,0.2};
  double conservation_jacobian[9]={0};
  if (chem_constraint_set_add_negative_jacobian(
          &conservation_set,conservation_other,NULL,1,conservation_jacobian)) return 60;
  for (size_t col=0;col<3;++col) {
    double plus[3],minus[3],fp[3]={0},fm[3]={0};
    memcpy(plus,conservation_other,sizeof(plus));
    memcpy(minus,conservation_other,sizeof(minus));
    plus[col]+=1.0e-7; minus[col]-=1.0e-7;
    if (chem_constraint_set_add_forcing(&conservation_set,plus,NULL,1,fp) ||
        chem_constraint_set_add_forcing(&conservation_set,minus,NULL,1,fm) ||
        fabs(conservation_jacobian[6+col]+(fp[2]-fm[2])/2.0e-7)>1.0e-7) return 61;
  }
  const chem_constraint_term zero_terms[] = {
      {.species=0,.coefficient=1.0}, {.species=1,.coefficient=-1.0}};
  const chem_constraint zero_linear = {.kind=CHEM_CONSTRAINT_LINEAR,
      .name="equal",.algebraic_species=1,.terms=zero_terms,.term_count=2,
      .constant=0.0};
  const chem_constraint_set zero_set = {.constraints=&zero_linear,
      .constraint_count=1,.species_count=2};
  const double equal_state[] = {0.4,0.4};
  double zero_forcing[2]={0};
  if (zero_linear.constant != 0.0 ||
      chem_constraint_set_add_forcing(&zero_set,equal_state,NULL,1,zero_forcing) ||
      fabs(zero_forcing[1])>1.0e-14) return 55;
  const chem_constraint_term fractional_terms[] = {
      {.species=0,.coefficient=0.5}, {.species=1,.coefficient=1.5}};
  const chem_constraint fractional = {.kind=CHEM_CONSTRAINT_LINEAR,
      .name="fractional",.algebraic_species=1,.terms=fractional_terms,
      .term_count=2,.constant=2.0};
  const chem_constraint_set fractional_set = {.constraints=&fractional,
      .constraint_count=1,.species_count=2};
  const double fractional_state[] = {1.0,1.0};
  double fractional_forcing[2]={0},fractional_jac[4]={0};
  if (fractional.terms[0].coefficient!=0.5 ||
      fractional.terms[1].coefficient!=1.5 ||
      chem_constraint_set_add_forcing(&fractional_set,fractional_state,NULL,1,
                                      fractional_forcing) ||
      fabs(fractional_forcing[1])>1.0e-14 ||
      chem_constraint_set_add_negative_jacobian(&fractional_set,
          fractional_state,NULL,1,fractional_jac) ||
      !near(fractional_jac[2],-0.5) || !near(fractional_jac[3],-1.5)) return 56;
  for (size_t col=0;col<2;++col) {
    double plus[2]={1.0,1.0},minus[2]={1.0,1.0};
    double fp[2]={0},fm[2]={0};
    plus[col]+=1.0e-7; minus[col]-=1.0e-7;
    if (chem_constraint_set_add_forcing(&fractional_set,plus,NULL,1,fp) ||
        chem_constraint_set_add_forcing(&fractional_set,minus,NULL,1,fm) ||
        fabs(fractional_jac[2+col]+(fp[1]-fm[1])/2.0e-7)>1.0e-7) return 57;
  }
  for (size_t col=0;col<3;++col) {
    double plus[3]={1.0,2.0,3.0},minus[3]={1.0,2.0,3.0};
    double fp[3]={0},fm[3]={0};
    plus[col]+=1.0e-7; minus[col]-=1.0e-7;
    if (chem_constraint_set_add_forcing(&weighted_set,plus,NULL,1,fp) ||
        chem_constraint_set_add_forcing(&weighted_set,minus,NULL,1,fm) ||
        fabs(weighted_jacobian[6+col]+(fp[2]-fm[2])/2.0e-7)>1.0e-7) return 62;
  }

  /* Validate multiple constraints, shared dependencies, and batched rows. */
  const chem_constraint_term second_terms[] = {
      {.species=0,.coefficient=1.0}, {.species=2,.coefficient=1.0}};
  const chem_constraint constraints[] = {
      {.kind=CHEM_CONSTRAINT_LINEAR,.name="first",.algebraic_species=1,
       .terms=ab_terms,.term_count=2,.constant=1.0},
      {.kind=CHEM_CONSTRAINT_LINEAR,.name="second",.algebraic_species=3,
       .terms=second_terms,.term_count=2,.constant=2.0}};
  const chem_constraint_set set = {.constraints=constraints,
      .constraint_count=2,.species_count=4};
  unsigned char set_mask[4], set_pattern[16];
  if (chem_constraint_set_validate(&set) ||
      chem_constraint_set_algebraic_mask(&set,set_mask) ||
      !set_mask[1] || !set_mask[3] || set_mask[0] || set_mask[2] ||
      chem_constraint_set_jacobian_pattern(&set,set_pattern) ||
      !set_pattern[4] || !set_pattern[5] || !set_pattern[12] ||
      !set_pattern[14]) return 6;
  for (size_t cells=1;cells<=4;++cells) {
    double y[16]={0}, f[16]={0};
    for (size_t cell=0;cell<cells;++cell) {
      y[4*cell]=0.1*(double)(cell+1);
      y[4*cell+1]=0.8;
      y[4*cell+2]=0.2*(double)(cell+1);
      y[4*cell+3]=0.5;
    }
    if (chem_constraint_set_add_forcing(&set,y,NULL,cells,f)) return 7;
    for (size_t cell=0;cell<cells;++cell) {
      if (!near(f[4*cell+1],y[4*cell]+y[4*cell+1]-1.0) ||
          !near(f[4*cell+3],y[4*cell]+y[4*cell+2]-2.0)) return 8;
    }
  }

  /* Equilibrium residual and analytic derivatives for ordinary and complex stoichiometry. */
  const chem_constraint_term eq_reactants[] = {
      {.species=0,.coefficient=1.0}, {.species=1,.coefficient=1.0}};
  const chem_constraint_term eq_products[] = {{.species=2,.coefficient=1.0}};
  const chem_constraint equilibrium = {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
      .name="A_B_equilibrium",.algebraic_species=2,
      .reactants=eq_reactants,.reactant_count=2,
      .products=eq_products,.product_count=1,
      .equilibrium_constant_ref=1000.0,.enthalpy_change=0.0,
      .reference_temperature=298.15};
  if (strcmp(equilibrium.name,"A_B_equilibrium") != 0 ||
      equilibrium.algebraic_species != 2 || equilibrium.reactant_count != 2 ||
      equilibrium.product_count != 1) return 58;
  const chem_constraint_set equilibrium_set = {.constraints=&equilibrium,
      .constraint_count=1,.species_count=3};
  const double equilibrium_state[] = {0.2,0.4,80.0};
  double equilibrium_forcing[3]={0}, equilibrium_jacobian[9]={0};
  if (chem_constraint_set_validate(&equilibrium_set)) return 91;
  if (chem_constraint_set_add_forcing(&equilibrium_set,equilibrium_state,
                                      temperatures,1,equilibrium_forcing)) return 92;
  if (fabs(equilibrium_forcing[2]) > 1.0e-10) return 93;
  if (chem_constraint_set_add_negative_jacobian(
          &equilibrium_set,equilibrium_state,temperatures,1,
          equilibrium_jacobian)) return 94;
  if (!near(equilibrium_jacobian[6],-400.0)) return 95;
  if (!near(equilibrium_jacobian[7],-200.0)) return 96;
  if (!near(equilibrium_jacobian[8],1.0)) return 97;
  const chem_constraint_term complex_reactant[] = {{.species=0,.coefficient=2.0}};
  const chem_constraint_term complex_products[] = {
      {.species=1,.coefficient=1.0}, {.species=2,.coefficient=1.0}};
  const chem_constraint complex_eq = {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
      .algebraic_species=1,.reactants=complex_reactant,.reactant_count=1,
      .products=complex_products,.product_count=2,
      .equilibrium_constant_ref=100.0,.enthalpy_change=0.0,
      .reference_temperature=298.15};
  const chem_constraint_set complex_set = {.constraints=&complex_eq,
      .constraint_count=1,.species_count=3};
  const double complex_state[] = {0.2,0.4,1.0};
  double complex_forcing[3]={0}, complex_jacobian[9]={0};
  if (chem_constraint_set_add_forcing(&complex_set,complex_state,NULL,1,
                                      complex_forcing) ||
      !near(complex_forcing[1],3.6) ||
      chem_constraint_set_add_negative_jacobian(
          &complex_set,complex_state,NULL,1,complex_jacobian) ||
      !near(complex_jacobian[3],-40.0) ||
      !near(complex_jacobian[4],1.0) ||
      !near(complex_jacobian[5],0.4)) return 10;
  /* Every analytic equilibrium derivative agrees with a finite difference. */
  for (size_t column=0;column<3;++column) {
    double plus[3],minus[3],fp[3]={0},fm[3]={0};
    memcpy(plus,complex_state,sizeof(plus));
    memcpy(minus,complex_state,sizeof(minus));
    const double h=1e-7;
    plus[column]+=h; minus[column]-=h;
    if (chem_constraint_set_add_forcing(&complex_set,plus,NULL,1,fp) ||
        chem_constraint_set_add_forcing(&complex_set,minus,NULL,1,fm)) return 11;
    const double fd=(fp[1]-fm[1])/(2*h);
    if (fabs(complex_jacobian[3+column]+fd)>1e-6*fmax(1.0,fabs(fd))) return 12;
  }

  const chem_constraint_term simple_reactants[] = {{.species=0,.coefficient=1.0}};
  const chem_constraint_term simple_products[] = {{.species=1,.coefficient=1.0}};
  const chem_constraint simple_equilibrium = {
      .kind=CHEM_CONSTRAINT_EQUILIBRIUM,.name="A_B_simple",
      .algebraic_species=1,.reactants=simple_reactants,.reactant_count=1,
      .products=simple_products,.product_count=1,
      .equilibrium_constant_ref=10.0,.enthalpy_change=-2400.0,
      .reference_temperature=298.15};
  const chem_constraint_set simple_set = {.constraints=&simple_equilibrium,
      .constraint_count=1,.species_count=2};
  const double simple_state[] = {0.2,2.0};
  double simple_forcing[2]={0}, simple_jacobian[4]={0};
  if (chem_constraint_set_add_forcing(&simple_set,simple_state,NULL,1,
                                      simple_forcing) ||
      fabs(simple_forcing[1])>1.0e-14 ||
      chem_constraint_set_add_negative_jacobian(&simple_set,simple_state,NULL,
                                                  1,simple_jacobian) ||
      !near(simple_jacobian[2],-10.0) || !near(simple_jacobian[3],1.0)) return 59;
  for (size_t col=0;col<2;++col) {
    double plus[2]={0.2,2.0},minus[2]={0.2,2.0};
    double fp[2]={0},fm[2]={0};
    plus[col]+=1.0e-7; minus[col]-=1.0e-7;
    if (chem_constraint_set_add_forcing(&simple_set,plus,NULL,1,fp) ||
        chem_constraint_set_add_forcing(&simple_set,minus,NULL,1,fm) ||
        fabs(simple_jacobian[2+col]+(fp[1]-fm[1])/2.0e-7)>1.0e-6) return 63;
  }

  const chem_constraint_set empty_set = {
      .constraints=NULL,.constraint_count=0,.species_count=3};
  unsigned char empty_mask[3]={1,1,1}, empty_pattern[9];
  if (chem_constraint_set_validate(&empty_set) ||
      chem_constraint_set_algebraic_mask(&empty_set,empty_mask) ||
      empty_mask[0] || empty_mask[1] || empty_mask[2] ||
      chem_constraint_set_jacobian_pattern(&empty_set,empty_pattern)) return 13;
  for (size_t i=0;i<9;++i) if (empty_pattern[i]) return 14;

  const chem_constraint_term invalid_term[] = {{.species=4,.coefficient=1.0}};
  const chem_constraint bad_linear = {.kind=CHEM_CONSTRAINT_LINEAR,
      .algebraic_species=0,.terms=invalid_term,.term_count=1};
  const chem_constraint_set bad_set = {.constraints=&bad_linear,
      .constraint_count=1,.species_count=3};
  const chem_constraint invalid_k = {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
      .algebraic_species=1,.reactants=eq_reactants,.reactant_count=2,
      .products=eq_products,.product_count=1,.equilibrium_constant_ref=0.0,
      .reference_temperature=298.15};
  const chem_constraint_set invalid_k_set = {.constraints=&invalid_k,
      .constraint_count=1,.species_count=3};
  const chem_constraint invalid_negative_k = {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
      .algebraic_species=1,.reactants=eq_reactants,.reactant_count=2,
      .products=eq_products,.product_count=1,.equilibrium_constant_ref=-1.0,
      .reference_temperature=298.15};
  const chem_constraint_set invalid_negative_k_set = {
      .constraints=&invalid_negative_k,.constraint_count=1,.species_count=3};
  const chem_constraint empty_reactants = {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
      .algebraic_species=1,.products=eq_products,.product_count=1,
      .equilibrium_constant_ref=1.0,.reference_temperature=298.15};
  const chem_constraint_set empty_reactants_set = {.constraints=&empty_reactants,
      .constraint_count=1,.species_count=3};
  const chem_constraint empty_products = {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
      .algebraic_species=1,.reactants=eq_reactants,.reactant_count=2,
      .equilibrium_constant_ref=1.0,.reference_temperature=298.15};
  const chem_constraint_set empty_products_set = {.constraints=&empty_products,
      .constraint_count=1,.species_count=3};
  const chem_constraint_term invalid_stoich[] = {{.species=0,.coefficient=0.0}};
  const chem_constraint bad_stoich = {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
      .algebraic_species=1,.reactants=invalid_stoich,.reactant_count=1,
      .products=eq_products,.product_count=1,.equilibrium_constant_ref=1.0,
      .reference_temperature=298.15};
  const chem_constraint_set bad_stoich_set = {.constraints=&bad_stoich,
      .constraint_count=1,.species_count=3};
  if (chem_constraint_set_validate(&bad_set)==0 ||
      chem_constraint_set_validate(&invalid_k_set)==0 ||
      chem_constraint_set_validate(&invalid_negative_k_set)==0 ||
      chem_constraint_set_validate(&empty_reactants_set)==0 ||
      chem_constraint_set_validate(&empty_products_set)==0 ||
      chem_constraint_set_validate(&bad_stoich_set)==0) return 15;
  return 0;
}

static int check_analytical_user_rate_chain(void) {
  const chem_species species[] = {{.name="A"},{.name="B"},{.name="C"}};
  const chem_term a[] = {{.species=0,.coefficient=1.0}};
  const chem_term b[] = {{.species=1,.coefficient=1.0}};
  const chem_term product_b[] = {{.species=1,.coefficient=1.0}};
  const chem_term product_c[] = {{.species=2,.coefficient=1.0}};
  const chem_reaction reactions[] = {
      {.kind=CHEM_USER_DEFINED,.reactants=a,.reactant_count=1,
       .products=product_b,.product_count=1,.p={1,0,0,0,0,0,0,0}},
      {.kind=CHEM_USER_DEFINED,.reactants=b,.reactant_count=1,
       .products=product_c,.product_count=1,.p={1,0,0,0,0,0,0,1}}};
  const char *labels[] = {"photoA","photoB"};
  const chem_mechanism mechanism = {.species=species,.species_count=3,
      .reactions=reactions,.reaction_count=2,.custom_labels=labels,
      .custom_parameter_count=2};
  const double rate_parameters[] = {2e-3,3e-3};
  double be_states[] = {1,0,0, 2.0/3.0,0,0, 1.0/3.0,0,0};
  double ros_states[9];
  memcpy(ros_states,be_states,sizeof(be_states));
  chem_solver_workspace be_workspace = {0}, ros_workspace = {0};
  if (chem_solver_workspace_init(&be_workspace,&mechanism) ||
      chem_solver_workspace_init(&ros_workspace,&mechanism)) {
    chem_solver_workspace_destroy(&be_workspace);
    chem_solver_workspace_destroy(&ros_workspace);
    return 1;
  }
  const chem_backward_euler_options be_options =
      chem_backward_euler_default_options();
  const chem_solver_options ros_options = chem_solver_default_options();
  for (size_t step = 1; step < 1000; ++step) {
    for (size_t cell = 0; cell < 3; ++cell) {
      chem_solver_stats stats;
      const int be_status = chem_backward_euler_solve_with_workspace(
          &mechanism,&be_workspace,be_states+3*cell,298.15,101325.0,40.9,
          rate_parameters,1.0,&be_options,&stats);
      const int ros_status = chem_rosenbrock3_solve_with_workspace(
          &mechanism,&ros_workspace,ros_states+3*cell,298.15,101325.0,40.9,
          rate_parameters,1.0,&ros_options,NULL);
      const double initial_a = 1.0-(double)cell/3.0;
      const double time = (double)step;
      const double exp1 = exp(-rate_parameters[0]*time);
      const double exp2 = exp(-rate_parameters[1]*time);
      const double expected[] = {
          initial_a*exp1,
          initial_a*rate_parameters[0]/(rate_parameters[1]-rate_parameters[0])*
              (exp1-exp2),
          initial_a*(1.0+(rate_parameters[0]*exp2-rate_parameters[1]*exp1)/
              (rate_parameters[1]-rate_parameters[0]))};
      if (be_status || ros_status || stats.accepted == 0) {
        chem_solver_workspace_destroy(&be_workspace);
        chem_solver_workspace_destroy(&ros_workspace);
        return 2;
      }
      for (size_t i = 0; i < 3; ++i)
        if (fabs(be_states[3*cell+i]-expected[i]) > 1e-3 ||
            fabs(ros_states[3*cell+i]-expected[i]) > 1e-3) {
          chem_solver_workspace_destroy(&be_workspace);
          chem_solver_workspace_destroy(&ros_workspace);
          return 3;
        }
    }
  }
  chem_solver_workspace_destroy(&be_workspace);
  chem_solver_workspace_destroy(&ros_workspace);
  return 0;
}

static int check_analytical_mechanism_chain(const chem_mechanism *mechanism,
                                           double backward_tolerance,
                                           double rosenbrock_tolerance) {
  const double temperatures[] = {272.5,254.7,312.6};
  const double pressures[] = {101253.3,100672.5,101319.8};
  const double gas_constant = 8.31446261815324;
  double be_states[] = {1,0,0, 2.0/3.0,0,0, 1.0/3.0,0,0};
  double ros_states[9];
  memcpy(ros_states,be_states,sizeof(be_states));
  chem_solver_workspace be_workspace = {0}, ros_workspace = {0};
  if (chem_solver_workspace_init(&be_workspace,mechanism) ||
      chem_solver_workspace_init(&ros_workspace,mechanism)) {
    chem_solver_workspace_destroy(&be_workspace);
    chem_solver_workspace_destroy(&ros_workspace);
    return 1;
  }
  double rates[3][2];
  for (size_t cell = 0; cell < 3; ++cell) {
    const double density = pressures[cell]/(gas_constant*temperatures[cell]);
    if (chem_rate_constants(mechanism,temperatures[cell],pressures[cell],
                            density,NULL,rates[cell])) {
      chem_solver_workspace_destroy(&be_workspace);
      chem_solver_workspace_destroy(&ros_workspace);
      return 2;
    }
  }
  const chem_backward_euler_options be_options =
      chem_backward_euler_default_options();
  const chem_solver_options ros_options = chem_solver_default_options();
  for (size_t step = 1; step < 1000; ++step)
    for (size_t cell = 0; cell < 3; ++cell) {
      const double density = pressures[cell]/(gas_constant*temperatures[cell]);
      chem_solver_stats stats;
      const int be_status = chem_backward_euler_solve_with_workspace(
          mechanism,&be_workspace,be_states+3*cell,temperatures[cell],
          pressures[cell],density,NULL,1.0,&be_options,&stats);
      const int ros_status = chem_rosenbrock3_solve_with_workspace(
          mechanism,&ros_workspace,ros_states+3*cell,temperatures[cell],
          pressures[cell],density,NULL,1.0,&ros_options,NULL);
      const double initial_a = 1.0-(double)cell/3.0;
      const double exp1 = exp(-rates[cell][0]*(double)step);
      const double exp2 = exp(-rates[cell][1]*(double)step);
      const double expected[] = {
          initial_a*exp1,
          initial_a*rates[cell][0]/(rates[cell][1]-rates[cell][0])*
              (exp1-exp2),
          initial_a*(1.0+(rates[cell][0]*exp2-rates[cell][1]*exp1)/
              (rates[cell][1]-rates[cell][0]))};
      if (be_status || ros_status || stats.accepted == 0) {
        chem_solver_workspace_destroy(&be_workspace);
        chem_solver_workspace_destroy(&ros_workspace);
        return 3;
      }
      for (size_t i = 0; i < 3; ++i)
        if (fabs(be_states[3*cell+i]-expected[i]) > backward_tolerance ||
            fabs(ros_states[3*cell+i]-expected[i]) > rosenbrock_tolerance) {
          chem_solver_workspace_destroy(&be_workspace);
          chem_solver_workspace_destroy(&ros_workspace);
          return 4;
        }
    }
  chem_solver_workspace_destroy(&be_workspace);
  chem_solver_workspace_destroy(&ros_workspace);
  return 0;
}

static int check_superstiff_mechanism(int variant, int backward_euler) {
  const chem_species species[]={{.name="A1"},{.name="A2"},
                                {.name="B"},{.name="C"}};
  const chem_term a1[]={{0,1}},a2[]={{1,1}},b[]={{2,1}};
  const chem_term to_b[]={{2,1}},to_c[]={{3,1}},to_a2[]={{1,1}},to_a1[]={{0,1}};
  chem_reaction reactions[]={
    {.reactants=a1,.reactant_count=1,.products=to_b,.product_count=1},
    {.reactants=a2,.reactant_count=1,.products=to_b,.product_count=1},
    {.reactants=b,.reactant_count=1,.products=to_c,.product_count=1},
    {.kind=CHEM_ARRHENIUS,.reactants=a1,.reactant_count=1,.products=to_a2,
      .product_count=1,.p={4e10,0,0,300,0}},
    {.kind=CHEM_ARRHENIUS,.reactants=a2,.reactant_count=1,.products=to_a1,
      .product_count=1,.p={3.6e10,0,0,300,0}}};
  double tolerance=1e-5;
  if (variant==0) {
    reactions[0].kind=reactions[1].kind=reactions[2].kind=CHEM_TROE;
    memcpy(reactions[0].p,(double[]){4e-11,0,0,1,0,0,0.6,1},8*sizeof(double));
    memcpy(reactions[1].p,reactions[0].p,8*sizeof(double));
    memcpy(reactions[2].p,(double[]){1.2e-3,1.6,3,136,5,24,0.9,0.8},8*sizeof(double));
  } else if (variant==1) {
    reactions[0].kind=reactions[1].kind=reactions[2].kind=CHEM_TERNARY_ACTIVATION;
    memcpy(reactions[0].p,(double[]){4e-5,0,0,1,0,0,0.6,1},8*sizeof(double));
    memcpy(reactions[1].p,reactions[0].p,8*sizeof(double));
    memcpy(reactions[2].p,(double[]){1.2e-3,1.6,3,136,5,24,0.9,0.8},8*sizeof(double));
    tolerance=backward_euler?1e-2:2e-3;
  } else if (variant==2) {
    reactions[0].kind=reactions[1].kind=reactions[2].kind=CHEM_TUNNELING;
    memcpy(reactions[0].p,(double[]){4e-3,0,0},3*sizeof(double));
    memcpy(reactions[1].p,reactions[0].p,3*sizeof(double));
    memcpy(reactions[2].p,(double[]){1.2e-4,167,1e8},3*sizeof(double));
    tolerance=backward_euler?1e-3:1e-4;
  } else if (variant==3) {
    reactions[0].kind=reactions[1].kind=reactions[2].kind=CHEM_ARRHENIUS;
    memcpy(reactions[0].p,(double[]){4e-3,0,50,300,0},5*sizeof(double));
    memcpy(reactions[1].p,reactions[0].p,5*sizeof(double));
    memcpy(reactions[2].p,(double[]){1.2e-4,1.6,75,50,0.5},5*sizeof(double));
    tolerance=backward_euler?1e-3:2e-5;
  } else {
    const double ref_a=2e-22*6.02214076e23*1e-6*exp(2.0)*2.45e19;
    const double ref_b=0.43*pow(293.0/298.0,-8.0);
    const double ref_av=ref_a/(1.0+ref_a/ref_b)*
        pow(0.41,1.0/(1.0+pow(log10(ref_a/ref_b),2.0)));
    const double z=ref_av*(1.0-1e-3)/1e-3;
    reactions[0].kind=reactions[1].kind=CHEM_BRANCHED_ALKOXY;
    reactions[2].kind=CHEM_BRANCHED_NITRATE;
    memcpy(reactions[0].p,(double[]){1e-4,204.3,1e-3,2.0,z},5*sizeof(double));
    memcpy(reactions[1].p,reactions[0].p,5*sizeof(double));
    memcpy(reactions[2].p,reactions[0].p,5*sizeof(double));
    tolerance=backward_euler?1e-2:2e-3;
  }
  const chem_mechanism mechanism={.species=species,.species_count=4,
      .reactions=reactions,.reaction_count=5};
  const double temperatures[]={272.5,254.7,312.6};
  const double pressures[]={101253.3,100672.5,101319.8};
  const double gas_constant=8.31446261815324;
  double rates[3][5];
  for (size_t c=0;c<3;++c) {
    const double density=pressures[c]/(gas_constant*temperatures[c]);
    if (chem_rate_constants(&mechanism,temperatures[c],pressures[c],density,
                            NULL,rates[c])) return 1;
  }
  const chem_rosenbrock_parameter_set schemes[]={CHEM_ROSENBROCK_TWO_STAGE,
      CHEM_ROSENBROCK_THREE_STAGE,CHEM_ROSENBROCK_FOUR_STAGE,
      CHEM_ROSENBROCK_FOUR_STAGE_DAE,CHEM_ROSENBROCK_SIX_STAGE_DAE};
  const size_t scheme_count=backward_euler?1:5;
  for (size_t scheme=0;scheme<scheme_count;++scheme) {
    double y[12]={0};
    for (size_t cell=0;cell<3;++cell) {
      const double total=1.0-(double)cell/3.0;
      y[cell*4]=0.5*total; y[cell*4+1]=0.5*total;
    }
    chem_solver_workspace workspace={0};
    if (chem_solver_workspace_init_group(&workspace,&mechanism,3)) return 2;
    const chem_solver_options ros_options=chem_solver_default_options();
    const chem_backward_euler_options be_options=chem_backward_euler_default_options();
    const chem_rosenbrock_parameters parameters=
        chem_rosenbrock_builtin_parameters(schemes[scheme]);
    for (size_t step=1;step<1000;++step) {
      for (size_t cell=0;cell<3;++cell) {
        const double density=pressures[cell]/(gas_constant*temperatures[cell]);
        chem_solver_stats stats={0};
        const int status=backward_euler
          ? chem_backward_euler_solve_with_workspace(&mechanism,&workspace,
                y+4*cell,temperatures[cell],pressures[cell],density,NULL,1.0,
                &be_options,&stats)
          : chem_rosenbrock_solve_with_parameters(&mechanism,&workspace,
                y+4*cell,temperatures[cell],pressures[cell],density,NULL,1.0,
                &parameters,&ros_options,&stats);
        const double initial=1.0-(double)cell/3.0;
        const double k1=rates[cell][0],k2=rates[cell][2],t=(double)step;
        const double e1=exp(-k1*t),e2=exp(-k2*t);
        const double expected[]={initial*e1,
          initial*k1/(k2-k1)*(e1-e2),
          initial*(1.0+(k1*e2-k2*e1)/(k2-k1))};
        if (status || stats.accepted==0 ||
            fabs(y[4*cell]+y[4*cell+1]-expected[0])>tolerance ||
            fabs(y[4*cell+2]-expected[1])>tolerance ||
            fabs(y[4*cell+3]-expected[2])>tolerance) {
          chem_solver_workspace_destroy(&workspace);
          return 3;
        }
      }
    }
    chem_solver_workspace_destroy(&workspace);
  }
  return 0;
}

static int check_backward_euler(void) {
  const double absolute_tolerances[] = {1e-6,1e-6,1e-6};
  double residual_data[12] = {0}, state_data[12] = {0};
  chem_const_matrix_view residual, state;
  int converged = 0;
  if (chem_const_matrix_view_init(&residual,residual_data,4,3,3,1) ||
      chem_const_matrix_view_init(&state,state_data,4,3,3,1) ||
      chem_backward_euler_is_converged(&residual,&state,absolute_tolerances,
                                       3,1e-3,1e-6,&converged) || !converged)
    return 1;
  residual_data[1] = 1e-5;
  if (chem_backward_euler_is_converged(&residual,&state,absolute_tolerances,
                                       3,1e-3,1e-6,&converged) || converged)
    return 2;
  residual_data[1] = 1e-5;
  if (chem_backward_euler_is_converged(&residual,&state,absolute_tolerances,
                                       3,1e-3,1e-4,&converged) || !converged)
    return 3;
  residual_data[11] = 1e-3;
  if (chem_backward_euler_is_converged(&residual,&state,absolute_tolerances,
                                       3,1e-3,1e-4,&converged) || converged)
    return 4;
  state_data[11] = 10.0;
  if (chem_backward_euler_is_converged(&residual,&state,absolute_tolerances,
                                       3,1e-3,1e-4,&converged) || !converged)
    return 5;
  residual_data[11] = 1e-1;
  if (chem_backward_euler_is_converged(&residual,&state,absolute_tolerances,
                                       3,1e-3,1e-4,&converged) || converged)
    return 6;
  const double relaxed_tolerances[] = {1e-6,1e-6,1.0};
  if (chem_backward_euler_is_converged(&residual,&state,relaxed_tolerances,
                                       3,1e-3,1e-4,&converged) || !converged)
    return 7;
  for (size_t width = 1; width <= 5; ++width) {
    double packed_residual[20] = {0}, packed_state[20] = {0};
    chem_vector_matrix_view residual_mut, state_mut;
    chem_const_vector_matrix_view residual_view, state_view;
    if (chem_vector_matrix_view_init(&residual_mut,packed_residual,4,3,width) ||
        chem_vector_matrix_view_init(&state_mut,packed_state,4,3,width)) return 8;
    for (size_t row = 0; row < 4; ++row)
      for (size_t column = 0; column < 3; ++column)
        if (chem_vector_matrix_view_set(&residual_mut,row,column,
              row == 0 && column == 1 ? 1e-5 : 0.0) ||
            chem_vector_matrix_view_set(&state_mut,row,column,0.0)) return 9;
    if (chem_const_vector_matrix_view_init(&residual_view,packed_residual,4,3,width) ||
        chem_const_vector_matrix_view_init(&state_view,packed_state,4,3,width) ||
        chem_backward_euler_vector_is_converged(
            &residual_view,&state_view,absolute_tolerances,3,1e-3,1e-6,
            &converged) || converged) return 10;
  }

  const chem_species species[] = {{.name="A"},{.name="B"},{.name="C"}};
  const chem_term a[] = {{.species=0,.coefficient=1.0}};
  const chem_term b[] = {{.species=1,.coefficient=1.0}};
  const chem_term product_b[] = {{.species=1,.coefficient=1.0}};
  const chem_term product_c[] = {{.species=2,.coefficient=1.0}};
  const chem_reaction reactions[] = {
      {.kind=CHEM_USER_DEFINED,.reactants=a,.reactant_count=1,
       .products=product_b,.product_count=1,.p={1,0,0,0,0,0,0,0}},
      {.kind=CHEM_USER_DEFINED,.reactants=b,.reactant_count=1,
       .products=product_c,.product_count=1,.p={1,0,0,0,0,0,0,1}}};
  const char *labels[] = {"k1","k2"};
  const chem_mechanism mechanism = {.species=species,.species_count=3,
      .reactions=reactions,.reaction_count=2,.custom_labels=labels,
      .custom_parameter_count=2};
  const double rates[] = {0.2,0.1};
  double concentrations[] = {1.0,0.0,0.0};
  chem_solver_workspace workspace = {0};
  if (chem_solver_workspace_init_group(&workspace,&mechanism,4)) return 11;
  chem_backward_euler_options options = chem_backward_euler_default_options();
  options.absolute_tolerances = absolute_tolerances;
  options.absolute_tolerance_count = 3;
  chem_solver_stats stats;
  const int status = chem_backward_euler_solve_with_workspace(
      &mechanism,&workspace,concentrations,272.5,101253.3,1e6,rates,1.0,
      &options,&stats);
  chem_solver_workspace_destroy(&workspace);
  if (status || stats.accepted == 0 || stats.solves == 0 ||
      fabs(concentrations[0]+concentrations[1]+concentrations[2]-1.0)>1e-12 ||
      concentrations[0] >= 1.0 || concentrations[2] <= 0.0) return 12;

  const chem_species troe_species[] = {{.name="A"},{.name="B"},{.name="C"}};
  const chem_term troe_a[] = {{.species=0,.coefficient=1.0}};
  const chem_term troe_b[] = {{.species=1,.coefficient=1.0}};
  const chem_term troe_product_b[] = {{.species=1,.coefficient=1.0}};
  const chem_term troe_product_c[] = {{.species=2,.coefficient=1.0}};
  const chem_reaction troe_reactions[] = {
      {.kind=CHEM_TROE,.reactants=troe_a,.reactant_count=1,
       .products=troe_product_b,.product_count=1,
       .p={4.0e-11,0,0,1,0,0,0.6,1}},
      {.kind=CHEM_TROE,.reactants=troe_b,.reactant_count=1,
       .products=troe_product_c,.product_count=1,
       .p={1.2e-3,1.6,3,136,5,24,0.9,0.8}}};
  const chem_mechanism troe_mechanism = {.species=troe_species,
      .species_count=3,.reactions=troe_reactions,.reaction_count=2};
  const double temperatures[] = {272.5,254.7,312.6};
  const double pressures[] = {101253.3,100672.5,101319.8};
  const double gas_constant = 8.31446261815324;
  double troe_states[] = {1.0,0.0,0.0, 2.0/3.0,0.0,0.0, 1.0/3.0,0.0,0.0};
  double rosenbrock_states[] = {1.0,0.0,0.0, 2.0/3.0,0.0,0.0, 1.0/3.0,0.0,0.0};
  double troe_rates[3][2];
  if (chem_solver_workspace_init(&workspace,&troe_mechanism)) return 13;
  chem_solver_workspace rosenbrock_workspace = {0};
  if (chem_solver_workspace_init(&rosenbrock_workspace,&troe_mechanism)) {
    chem_solver_workspace_destroy(&workspace);
    return 13;
  }
  const chem_solver_options rosenbrock_options = chem_solver_default_options();
  for (size_t cell = 0; cell < 3; ++cell) {
    const double density = pressures[cell]/(gas_constant*temperatures[cell]);
    if (chem_rate_constants(&troe_mechanism,temperatures[cell],pressures[cell],
                            density,NULL,troe_rates[cell])) return 14;
  }
  for (size_t step = 1; step < 1000; ++step) {
    for (size_t cell = 0; cell < 3; ++cell) {
      const double density = pressures[cell]/(gas_constant*temperatures[cell]);
      chem_solver_stats one = {0};
      const int troe_status = chem_backward_euler_solve_with_workspace(
          &troe_mechanism,&workspace,troe_states+3*cell,temperatures[cell],
          pressures[cell],density,NULL,1.0,NULL,&one);
      if (troe_status || !one.accepted) {
        chem_solver_workspace_destroy(&workspace);
        chem_solver_workspace_destroy(&rosenbrock_workspace);
        return 15;
      }
      const int rosenbrock_status = chem_rosenbrock3_solve_with_workspace(
          &troe_mechanism,&rosenbrock_workspace,rosenbrock_states+3*cell,
          temperatures[cell],pressures[cell],density,NULL,1.0,
          &rosenbrock_options,NULL);
      if (rosenbrock_status) {
        chem_solver_workspace_destroy(&workspace);
        chem_solver_workspace_destroy(&rosenbrock_workspace);
        return 17;
      }
      const double initial_a = 1.0-(double)cell/3.0;
      const double k1 = troe_rates[cell][0], k2 = troe_rates[cell][1];
      const double time = (double)step;
      const double exp1 = exp(-k1*time), exp2 = exp(-k2*time);
      const double expected_a = initial_a*exp1;
      const double expected_b = initial_a*k1/(k2-k1)*(exp1-exp2);
      const double expected_c = initial_a*(1.0+(k1*exp2-k2*exp1)/(k2-k1));
      if (fabs(troe_states[3*cell]-expected_a)>1e-6 ||
          fabs(troe_states[3*cell+1]-expected_b)>1e-6 ||
          fabs(troe_states[3*cell+2]-expected_c)>1e-6) {
        chem_solver_workspace_destroy(&workspace);
        chem_solver_workspace_destroy(&rosenbrock_workspace);
        return 16;
      }
      if (fabs(rosenbrock_states[3*cell]-expected_a)>1e-10 ||
          fabs(rosenbrock_states[3*cell+1]-expected_b)>1e-10 ||
          fabs(rosenbrock_states[3*cell+2]-expected_c)>1e-10) {
        chem_solver_workspace_destroy(&workspace);
        chem_solver_workspace_destroy(&rosenbrock_workspace);
        return 16;
      }
    }
  }
  chem_solver_workspace_destroy(&workspace);
  chem_solver_workspace_destroy(&rosenbrock_workspace);
  if (check_analytical_user_rate_chain()) return 18;
  const chem_species chain_species[] = {{.name="A"},{.name="B"},{.name="C"}};
  const chem_term chain_a[] = {{.species=0,.coefficient=1.0}};
  const chem_term chain_b[] = {{.species=1,.coefficient=1.0}};
  const chem_term chain_product_b[] = {{.species=1,.coefficient=1.0}};
  const chem_term chain_product_c[] = {{.species=2,.coefficient=1.0}};
  chem_reaction chain_reactions[] = {
      {.kind=CHEM_TUNNELING,.reactants=chain_a,.reactant_count=1,
       .products=chain_product_b,.product_count=1,.p={4e-3,0,0}},
      {.kind=CHEM_TUNNELING,.reactants=chain_b,.reactant_count=1,
       .products=chain_product_c,.product_count=1,.p={1.2e-4,167,1e8}}};
  chem_mechanism chain_mechanism = {.species=chain_species,.species_count=3,
      .reactions=chain_reactions,.reaction_count=2};
  if (check_analytical_mechanism_chain(&chain_mechanism,1e-3,1e-8)) return 19;
  chain_reactions[0] = (chem_reaction){.kind=CHEM_ARRHENIUS,
      .reactants=chain_a,.reactant_count=1,.products=chain_product_b,
      .product_count=1,.p={4e-3,0,50,300,0}};
  chain_reactions[1] = (chem_reaction){.kind=CHEM_ARRHENIUS,
      .reactants=chain_b,.reactant_count=1,.products=chain_product_c,
      .product_count=1,.p={1.2e-4,7,75,50,0.5}};
  if (check_analytical_mechanism_chain(&chain_mechanism,1e-3,1e-9)) return 20;
  const double branch_reference_a =
      2e-22*6.02214076e23*1e-6*exp(2.0)*2.45e19;
  const double branch_reference_b = 0.43*pow(293.0/298.0,-8.0);
  const double branch_reference_ratio = branch_reference_a/branch_reference_b;
  const double branch_reference_av = branch_reference_a/(1.0+branch_reference_ratio)*
      pow(0.41,1.0/(1.0+pow(log10(branch_reference_ratio),2.0)));
  const double branch_z = branch_reference_av*(1.0-1e-3)/1e-3;
  chain_reactions[0] = (chem_reaction){.kind=CHEM_BRANCHED_ALKOXY,
      .reactants=chain_a,.reactant_count=1,.products=chain_product_b,
      .product_count=1,.p={1e-4,204.3,1e-3,2.0,branch_z}};
  chain_reactions[1] = (chem_reaction){.kind=CHEM_BRANCHED_NITRATE,
      .reactants=chain_b,.reactant_count=1,.products=chain_product_c,
      .product_count=1,.p={1e-4,204.3,1e-3,2.0,branch_z}};
  if (check_analytical_mechanism_chain(&chain_mechanism,1e-5,1e-13)) return 21;
  chain_reactions[0] = (chem_reaction){.kind=CHEM_TERNARY_ACTIVATION,
      .reactants=chain_a,.reactant_count=1,.products=chain_product_b,
      .product_count=1,.p={4e-5,0,0,1,0,0,0.6,1}};
  chain_reactions[1] = (chem_reaction){.kind=CHEM_TERNARY_ACTIVATION,
      .reactants=chain_b,.reactant_count=1,.products=chain_product_c,
      .product_count=1,.p={1.2e-3,1.6,3,136,5,24,0.9,0.8}};
  if (check_analytical_mechanism_chain(&chain_mechanism,1e-5,1e-8)) return 22;
  const chem_species surface_species[] = {
      {.name="foo",.molecular_weight=0.04607,.diffusion=0.95e-5},
      {.name="bar"},{.name="baz"}};
  const chem_term surface_reactants[] = {{.species=0,.coefficient=1.0}};
  const chem_term surface_products[] = {{.species=1,.coefficient=1.0},
                                        {.species=2,.coefficient=0.4}};
  const chem_reaction surface_reaction = {.kind=CHEM_SURFACE,
      .reactants=surface_reactants,.reactant_count=1,
      .products=surface_products,.product_count=2,.surface_species=0,
      .p={0.02,0,0,0,0,0,0,0}};
  const char *surface_labels[] = {"radius","number"};
  const chem_mechanism surface_mechanism = {.species=surface_species,
      .species_count=3,.reactions=&surface_reaction,.reaction_count=1,
      .custom_labels=surface_labels,.custom_parameter_count=2};
  const double temperature=272.5, pressure=101253.3;
  const double radius=1e-6/2.0*exp(5.0*log(0.1)*log(0.1)/2.0);
  const double number_concentration=6.0/(3.14159265358979323846*pow(1e-6,3.0)*
      exp(9.0/2.0*log(0.1)*log(0.1)))*(2e-3/1000.0+3e-3/1000.0);
  const double surface_parameters[] = {radius,number_concentration};
  const double mean_speed=sqrt(8.0*gas_constant/(3.14159265358979323846*0.04607)*
                               temperature);
  const double surface_k=4.0*number_concentration*3.14159265358979323846*
      radius*radius/(radius/0.95e-5+4.0/(mean_speed*0.02));
  double surface_rate;
  if (chem_rate_constants(&surface_mechanism,temperature,pressure,
        pressure/(gas_constant*temperature),surface_parameters,&surface_rate) ||
      fabs(surface_rate-surface_k)>1e-12*fabs(surface_k)) return 23;
  chem_solver_workspace surface_be_workspace={0},surface_ros_workspace={0};
  if (chem_solver_workspace_init(&surface_be_workspace,&surface_mechanism) ||
      chem_solver_workspace_init(&surface_ros_workspace,&surface_mechanism)) {
    chem_solver_workspace_destroy(&surface_be_workspace);
    chem_solver_workspace_destroy(&surface_ros_workspace);
    return 24;
  }
  double surface_be_state[]={1.0,0.0,0.0};
  double surface_ros_state[]={1.0,0.0,0.0};
  const double surface_dt=0.1/surface_k;
  const chem_backward_euler_options surface_be_options=
      chem_backward_euler_default_options();
  const chem_solver_options surface_ros_options=chem_solver_default_options();
  for (size_t step=1;step<=10;++step) {
    chem_solver_stats surface_stats;
    const int be_status=chem_backward_euler_solve_with_workspace(
        &surface_mechanism,&surface_be_workspace,surface_be_state,temperature,
        pressure,pressure/(gas_constant*temperature),surface_parameters,
        surface_dt,&surface_be_options,&surface_stats);
    const int ros_status=chem_rosenbrock3_solve_with_workspace(
        &surface_mechanism,&surface_ros_workspace,surface_ros_state,temperature,
        pressure,pressure/(gas_constant*temperature),surface_parameters,
        surface_dt,&surface_ros_options,NULL);
    const double expected_foo=exp(-surface_k*surface_dt*(double)step);
    const double expected_bar=1.0-expected_foo;
    const double expected_baz=0.4*(1.0-expected_foo);
    if (be_status || ros_status || !surface_stats.accepted ||
        fabs(surface_be_state[0]-expected_foo)>0.05 ||
        fabs(surface_be_state[1]-expected_bar)>0.05 ||
        fabs(surface_be_state[2]-expected_baz)>0.05 ||
        fabs(surface_ros_state[0]-expected_foo)>1e-5 ||
        fabs(surface_ros_state[1]-expected_bar)>1e-5 ||
        fabs(surface_ros_state[2]-expected_baz)>1e-5) {
      chem_solver_workspace_destroy(&surface_be_workspace);
      chem_solver_workspace_destroy(&surface_ros_workspace);
      return 25;
    }
  }
  chem_solver_workspace_destroy(&surface_be_workspace);
  chem_solver_workspace_destroy(&surface_ros_workspace);
  return 0;
}

static int check_constrained_solver(void) {
  if (strcmp(chem_solver_status_string(
                 CHEM_SOLVER_STATUS_CONSTRAINT_INITIALIZATION_FAILED),
             "Constraint Initialization Failed") ||
      strcmp(chem_solver_status_string((chem_solver_status)99),"Unknown"))
    return 16;
  const chem_species species[] = {{.name="X",.initial_concentration=0.0}};
  const chem_mechanism mechanism = {.species=species,.species_count=1};
  const chem_constraint_term terms[] = {{.species=0,.coefficient=1.0}};
  const chem_constraint constraint = {.kind=CHEM_CONSTRAINT_LINEAR,
      .algebraic_species=0,.terms=terms,.term_count=1,.constant=1.0};
  const chem_constraint_set constraints = {.constraints=&constraint,
      .constraint_count=1,.species_count=1};
  chem_solver_workspace workspace = {0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  double state[] = {0.0};
  chem_solver_stats stats = {0};
  const chem_rosenbrock_parameters parameters =
      chem_rosenbrock_builtin_parameters(CHEM_ROSENBROCK_FOUR_STAGE_DAE);
  const int status = chem_rosenbrock_solve_constrained_with_parameters(
      &mechanism,&workspace,state,298.15,101325.0,40.0,NULL,1.0,
      &parameters,&constraints,NULL,&stats);
  if (status || fabs(state[0]-1.0)>1.0e-12 || stats.accepted==0 ||
      stats.constraint_init_iterations == 0) return 2;
  const chem_constraint_term scaled_terms[] = {{.species=0,.coefficient=1.0e-9}};
  const chem_constraint scaled_constraint = {.kind=CHEM_CONSTRAINT_LINEAR,
      .algebraic_species=0,.terms=scaled_terms,.term_count=1,.constant=1.0e-9};
  const chem_constraint_set scaled_set = {.constraints=&scaled_constraint,
      .constraint_count=1,.species_count=1};
  double scaled_state[] = {0.0};
  const int scaled_status = chem_rosenbrock_solve_constrained_with_parameters(
      &mechanism,&workspace,scaled_state,298.15,101325.0,40.0,NULL,1.0,
      &parameters,&scaled_set,NULL,&stats);
  chem_solver_workspace_destroy(&workspace);
  if (scaled_status || fabs(scaled_state[0]-1.0)>1.0e-12) return 3;

  const chem_species pair_species[] = {
      {.name="B",.initial_concentration=0.5},
      {.name="C",.initial_concentration=1000.0}};
  const chem_mechanism pair_mechanism = {.species=pair_species,.species_count=2};
  const chem_constraint_term reactants[] = {{.species=0,.coefficient=1.0}};
  const chem_constraint_term products[] = {{.species=1,.coefficient=1.0}};
  const chem_constraint equilibrium = {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
      .algebraic_species=1,.reactants=reactants,.reactant_count=1,
      .products=products,.product_count=1,.equilibrium_constant_ref=2.0,
      .reference_temperature=298.15};
  const chem_constraint_set equilibrium_set = {.constraints=&equilibrium,
      .constraint_count=1,.species_count=2};
  chem_solver_workspace pair_workspace = {0};
  if (chem_solver_workspace_init(&pair_workspace,&pair_mechanism)) return 3;
  double pair_state[] = {0.5,1000.0};
  const int pair_status = chem_rosenbrock_solve_constrained_with_parameters(
      &pair_mechanism,&pair_workspace,pair_state,298.15,101325.0,40.0,NULL,0.1,
      NULL,&equilibrium_set,NULL,&stats);
  if (pair_status || fabs(pair_state[0]-0.5)>1.0e-12 ||
      fabs(pair_state[1]-1.0)>1.0e-12 || stats.constraint_init_iterations == 0)
    return 4;
  const double consistent_state[] = {0.5,1.0};
  pair_state[0] = consistent_state[0]; pair_state[1] = consistent_state[1];
  const int consistent_status = chem_rosenbrock_solve_constrained_with_parameters(
      &pair_mechanism,&pair_workspace,pair_state,298.15,101325.0,40.0,NULL,0.1,
      NULL,&equilibrium_set,NULL,&stats);
  if (consistent_status || fabs(pair_state[0]-0.5)>1.0e-12 ||
      fabs(pair_state[1]-1.0)>1.0e-12 || stats.constraint_init_iterations > 1)
    return 5;
  pair_state[1] = 8.0;
  const int repeat_status = chem_rosenbrock_solve_constrained_with_parameters(
      &pair_mechanism,&pair_workspace,pair_state,298.15,101325.0,40.0,NULL,0.1,
      NULL,&equilibrium_set,NULL,&stats);
  chem_solver_workspace_destroy(&pair_workspace);
  if (repeat_status || fabs(pair_state[1]-1.0)>1.0e-12 ||
      stats.constraint_init_iterations == 0) return 6;

  const chem_constraint_term weighted_terms[] = {
      {.species=0,.coefficient=1.0},{.species=1,.coefficient=-1.0}};
  const chem_constraint weighted_constraint = {.kind=CHEM_CONSTRAINT_LINEAR,
      .algebraic_species=1,.terms=weighted_terms,.term_count=2};
  const chem_constraint_set weighted_set = {.constraints=&weighted_constraint,
      .constraint_count=1,.species_count=2};
  chem_solver_workspace weighted_workspace = {0};
  if (chem_solver_workspace_init(&weighted_workspace,&pair_mechanism)) return 7;
  const double deviation = 5.0e-8;
  double weighted_state[] = {1.0,1.0+deviation};
  const chem_solver_options loose = {
      .relative_tolerance=1.0e-6, .absolute_tolerance=1.0e-8,
      .max_steps=1000};
  const int loose_status = chem_rosenbrock_solve_constrained_with_parameters(
      &pair_mechanism,&weighted_workspace,weighted_state,298.15,101325.0,
      40.0,NULL,0.1,NULL,&weighted_set,&loose,&stats);
  if (loose_status || weighted_state[1] != 1.0+deviation ||
      stats.constraint_init_iterations != 1) return 8;
  weighted_state[0] = 1.0; weighted_state[1] = 1.0+deviation;
  const chem_solver_options tight = {
      .relative_tolerance=1.0e-8, .absolute_tolerance=1.0e-12,
      .max_steps=1000};
  const int tight_status = chem_rosenbrock_solve_constrained_with_parameters(
      &pair_mechanism,&weighted_workspace,weighted_state,298.15,101325.0,
      40.0,NULL,0.1,NULL,&weighted_set,&tight,&stats);
  chem_solver_workspace_destroy(&weighted_workspace);
  if (tight_status || fabs(weighted_state[1]-1.0)>1.0e-12 ||
      stats.constraint_init_iterations != 2) return 9;

  const chem_constraint_term squared_reactants[] = {{.species=0,.coefficient=1.0}};
  const chem_constraint_term squared_products[] = {{.species=1,.coefficient=2.0}};
  const chem_constraint squared_constraint = {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
      .algebraic_species=1,.reactants=squared_reactants,.reactant_count=1,
      .products=squared_products,.product_count=1,.equilibrium_constant_ref=2.0,
      .reference_temperature=298.15};
  const chem_constraint_set squared_set = {.constraints=&squared_constraint,
      .constraint_count=1,.species_count=2};
  chem_solver_workspace squared_workspace = {0};
  if (chem_solver_workspace_init(&squared_workspace,&pair_mechanism)) return 10;
  double squared_state[] = {0.5,1.0e-6};
  const int squared_status = chem_rosenbrock_solve_constrained_with_parameters(
      &pair_mechanism,&squared_workspace,squared_state,298.15,101325.0,40.0,
      NULL,0.1,NULL,&squared_set,NULL,&stats);
  chem_solver_workspace_destroy(&squared_workspace);
  if (squared_status || fabs(squared_state[1]-1.0)>1.0e-3 ||
      stats.constraint_init_iterations <= 1) return 11;

  const chem_species singular_species[] = {{.name="X"},{.name="Y"}};
  const chem_mechanism singular_mechanism = {
      .species=singular_species,.species_count=2};
  const chem_constraint_term dependent_only[] = {{.species=1,.coefficient=1.0}};
  const chem_constraint singular_constraint = {.kind=CHEM_CONSTRAINT_LINEAR,
      .algebraic_species=0,.terms=dependent_only,.term_count=1,.constant=1.0};
  const chem_constraint_set singular_set = {.constraints=&singular_constraint,
      .constraint_count=1,.species_count=2};
  chem_solver_workspace singular_workspace = {0};
  if (chem_solver_workspace_init(&singular_workspace,&singular_mechanism)) return 12;
  double singular_state[] = {7.0,0.0};
  const int singular_status = chem_rosenbrock_solve_constrained_with_parameters(
      &singular_mechanism,&singular_workspace,singular_state,298.15,101325.0,
      40.0,NULL,0.1,NULL,&singular_set,NULL,NULL);
  chem_solver_workspace_destroy(&singular_workspace);
  if (singular_status == 0 || singular_state[0] != 7.0 || singular_state[1] != 0.0)
    return 13;

  chem_solver_workspace group_workspace = {0};
  if (chem_solver_workspace_init_group(&group_workspace,&mechanism,3)) return 14;
  double group_states[] = {0.0,2.0,-1.0};
  const double group_temperatures[] = {298.15,300.0,280.0};
  const double group_pressures[] = {101325.0,101325.0,101325.0};
  const double group_densities[] = {40.0,40.0,40.0};
  const int group_status = chem_rosenbrock_solve_constrained_group_with_parameters(
      &mechanism,&group_workspace,group_states,group_temperatures,
      group_pressures,group_densities,NULL,3,0.1,NULL,&constraints,NULL,&stats);
  chem_solver_workspace_destroy(&group_workspace);
  if (group_status || fabs(group_states[0]-1.0)>1.0e-12 ||
      fabs(group_states[1]-1.0)>1.0e-12 || fabs(group_states[2]-1.0)>1.0e-12 ||
      stats.constraint_init_iterations < 3) return 15;
  return 0;
}

static int check_equilibrium_dae_case(int mode) {
  chem_species species[4]={{.name="A"},{.name="B"},{.name="C"},{.name="D"}};
  chem_term reactant={.species=0,.coefficient=1.0};
  chem_term reaction_products={.species=1,.coefficient=1.0};
  chem_reaction reaction={.kind=CHEM_ARRHENIUS,.reactants=&reactant,
      .reactant_count=1,.products=&reaction_products,.product_count=1,
      .p={1.0,0,0,300,0}};
  chem_constraint_term eq1_reactants[2],eq1_products[1],eq2_products[1];
  chem_constraint constraints[2]={0};
  size_t constraint_count=1, species_count=3, steps=100;
  double state[4]={1.0,0.0,0.0,0.0};
  if (mode==1) {
    /* Same mechanism with algebraic C first in the caller's ordering. */
    species[0].name="C"; species[1].name="A"; species[2].name="B";
    state[0]=0.0; state[1]=1.0; state[2]=0.0;
    reactant.species=1; reaction_products.species=2;
    eq1_reactants[0]=(chem_constraint_term){.species=2,.coefficient=1.0};
    eq1_products[0]=(chem_constraint_term){.species=0,.coefficient=1.0};
    constraints[0]=(chem_constraint){.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
        .name="B_C_eq",.algebraic_species=0,.reactants=eq1_reactants,
        .reactant_count=1,.products=eq1_products,.product_count=1,
        .equilibrium_constant_ref=2.0,.enthalpy_change=-2400.0,
        .reference_temperature=298.15};
  } else if (mode==2) {
    species_count=4; constraint_count=2;
    eq1_reactants[0]=(chem_constraint_term){.species=1,.coefficient=1.0};
    eq1_products[0]=(chem_constraint_term){.species=2,.coefficient=1.0};
    eq2_products[0]=(chem_constraint_term){.species=3,.coefficient=1.0};
    constraints[0]=(chem_constraint){.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
        .name="B_C_eq",.algebraic_species=2,.reactants=eq1_reactants,
        .reactant_count=1,.products=eq1_products,.product_count=1,
        .equilibrium_constant_ref=3.0,.enthalpy_change=-2400.0,
        .reference_temperature=298.15};
    constraints[1]=constraints[0];
    constraints[1].name="B_D_eq"; constraints[1].algebraic_species=3;
    constraints[1].products=eq2_products;
    constraints[1].equilibrium_constant_ref=5.0;
    state[0]=1.0; state[1]=0.0;
  } else if (mode==3) {
    reaction.reactants=&(chem_term){.species=2,.coefficient=1.0};
    reaction_products.species=0;
    reactant.species=2; reaction.reactants=&reactant;
    reaction_products.species=0;
    reaction.p[0]=0.5;
    eq1_reactants[0]=(chem_constraint_term){.species=0,.coefficient=2.0};
    eq1_products[0]=(chem_constraint_term){.species=1,.coefficient=1.0};
    constraints[0]=(chem_constraint){.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
        .name="A2_B_eq",.algebraic_species=1,.reactants=eq1_reactants,
        .reactant_count=1,.products=eq1_products,.product_count=1,
        .equilibrium_constant_ref=10.0,.enthalpy_change=-2400.0,
        .reference_temperature=298.15};
    state[0]=0.1; state[1]=0.1; state[2]=1.0; steps=50;
  } else {
    eq1_reactants[0]=(chem_constraint_term){.species=1,.coefficient=1.0};
    eq1_products[0]=(chem_constraint_term){.species=2,.coefficient=1.0};
    constraints[0]=(chem_constraint){.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
        .name="B_C_eq",.algebraic_species=2,.reactants=eq1_reactants,
        .reactant_count=1,.products=eq1_products,.product_count=1,
        .equilibrium_constant_ref=2.0,.enthalpy_change=-2400.0,
        .reference_temperature=298.15};
  }
  reaction.reactants=&reactant;
  const chem_mechanism mechanism={.species=species,.species_count=species_count,
      .reactions=&reaction,.reaction_count=1};
  const chem_constraint_set set={.constraints=constraints,
      .constraint_count=constraint_count,.species_count=species_count};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  const chem_rosenbrock_parameters parameters=
      chem_rosenbrock_builtin_parameters(CHEM_ROSENBROCK_THREE_STAGE);
  int status=0;
  for (size_t i=0;i<steps && !status;++i) {
    status=chem_rosenbrock_solve_constrained_with_parameters(&mechanism,
        &workspace,state,mode==1?270.0:mode==2?300.0:mode==3?298.0:270.0,
        101325.0,40.0,NULL,0.001,&parameters,&set,NULL,NULL);
    double residual[4]={0};
    if (!status && chem_constraint_set_add_forcing(&set,state,
            &(double){mode==1?270.0:mode==2?300.0:mode==3?298.0:270.0},
            1,residual)) status=2;
    if (!status)
      for (size_t j=0;j<constraint_count;++j)
        if (fabs(residual[constraints[j].algebraic_species])>1e-6) {
          status=3;
        }
  }
  chem_solver_workspace_destroy(&workspace);
  if (!status && mode!=3 && fabs(state[0]+state[1]-1.0)>0.01 && mode!=1)
    status=4;
  return status;
}

static int check_robertson_terminator_dae(void) {
  const chem_species species[]={
      {.name="Cl2"},{.name="Cl"},{.name="A"},{.name="B"},{.name="C"}};
  const chem_term cl2[]={ {0,1.0} };
  const chem_term two_cl[]={ {1,1.0},{1,1.0} };
  const chem_term a[]={ {2,1.0} }, two_b[]={ {3,1.0},{3,1.0} },
      bc[]={ {3,1.0},{4,1.0} };
  const chem_term two_cl_product[]={ {1,2.0} }, cl2_product[]={ {0,1.0} };
  const chem_term b_product[]={ {3,1.0} }, bc_product[]={ {3,1.0},{4,1.0} };
  const chem_term bc_to_ac[]={ {2,1.0},{4,1.0} };
  const chem_reaction reactions[]={
      {.kind=CHEM_USER_DEFINED,.label="terminator_k1",.reactants=cl2,
       .reactant_count=1,.products=two_cl_product,.product_count=1,
       .p={1,0,0,0,0,0,0,3}},
      {.kind=CHEM_ARRHENIUS,.reactants=two_cl,.reactant_count=2,
       .products=cl2_product,.product_count=1,.p={1,0,0,300,0}},
      {.kind=CHEM_USER_DEFINED,.label="robertson_r1",.reactants=a,
       .reactant_count=1,.products=b_product,.product_count=1,
       .p={1,0,0,0,0,0,0,0}},
      {.kind=CHEM_USER_DEFINED,.label="robertson_r2",.reactants=two_b,
       .reactant_count=2,.products=bc_product,.product_count=2,
       .p={1,0,0,0,0,0,0,1}},
      {.kind=CHEM_USER_DEFINED,.label="robertson_r3",.reactants=bc,
       .reactant_count=2,.products=bc_to_ac,.product_count=2,
       .p={1,0,0,0,0,0,0,2}}};
  const char *labels[]={"robertson_r1","robertson_r2","robertson_r3",
                        "terminator_k1"};
  const chem_mechanism mechanism={.species=species,.species_count=5,
      .reactions=reactions,.reaction_count=5,.custom_labels=labels,
      .custom_parameter_count=4};
  const chem_constraint_term conservation_terms[]={
      {2,1.0},{3,1.0},{4,1.0}};
  const chem_constraint conservation={.kind=CHEM_CONSTRAINT_LINEAR,
      .name="mass_conservation",.algebraic_species=4,
      .terms=conservation_terms,.term_count=3,.constant=1.0};
  const chem_constraint_set constraints={.constraints=&conservation,
      .constraint_count=1,.species_count=5};
  const chem_rosenbrock_parameters parameters=
      chem_rosenbrock_builtin_parameters(CHEM_ROSENBROCK_FOUR_STAGE_DAE);
  const double atol[]={1e-20,1e-20,1e-3,1e-3,1e-3};
  const chem_solver_options options={.relative_tolerance=1e-8,
      .absolute_tolerance=1e-3,.absolute_tolerances=atol,
      .absolute_tolerance_count=5,.max_steps=1000,.initial_step=0.0};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  double state[]={1.8e-10,1.2e-6,1.0,0.0,0.0};
  const double chlorine_initial=2.0*state[0]+state[1];
  const double rates[]={0.04,3.0e7,1.0e4,0.0};
  double duration=1.0;
  int status=0;
  for (size_t interval=0;interval<12 && !status;++interval) {
    status=chem_rosenbrock_solve_constrained_with_parameters(&mechanism,
        &workspace,state,298.0,101300.0,42.0,rates,duration,&parameters,
        &constraints,&options,NULL);
    double residual[5]={0};
    if (!status && chem_constraint_set_add_forcing(&constraints,state,NULL,1,
                                                    residual)) status=2;
    if (!status && (fabs(residual[4])>1e-10 ||
        fabs(state[2]+state[3]+state[4]-1.0)>1e-10 ||
        fabs((2.0*state[0]+state[1])-chlorine_initial)>1e-6*chlorine_initial))
      status=3;
    for (size_t i=0;i<5 && !status;++i)
      if (!isfinite(state[i]) || state[i]<0.0) status=4;
    if (interval==0 && !status &&
        (state[2]>=0.99 || state[3]<=0.0 || state[4]<=1e-3)) status=5;
    duration*=10.0;
  }
  chem_solver_workspace_destroy(&workspace);
  if (!status && (state[4]<=0.95 || state[2]>=0.05 ||
      state[0]<=0.4*chlorine_initial || state[1]>=0.1*chlorine_initial)) status=6;
  return status;
}

static int check_constrained_overshoot_order(chem_rosenbrock_parameter_set scheme) {
  const chem_species species[]={{.name="A"},{.name="B"},{.name="C"}};
  const chem_term reactant[]={{0,1.0}}, product[]={{1,1.0}};
  const chem_reaction reaction={.kind=CHEM_ARRHENIUS,.reactants=reactant,
      .reactant_count=1,.products=product,.product_count=1,
      .p={1.0e4,0,0,298,0}};
  const chem_mechanism mechanism={.species=species,.species_count=3,
      .reactions=&reaction,.reaction_count=1};
  const chem_constraint_term terms[]={{0,1.0},{1,1.0},{2,1.0}};
  const chem_constraint constraint={.kind=CHEM_CONSTRAINT_LINEAR,
      .name="mass_conservation",.algebraic_species=2,.terms=terms,
      .term_count=3,.constant=1.0e-6};
  const chem_constraint_set set={.constraints=&constraint,
      .constraint_count=1,.species_count=3};
  const chem_rosenbrock_parameters parameters=chem_rosenbrock_builtin_parameters(scheme);
  const double atol[]={1.0e-12,1.0e-12,1.0e-12};
  const chem_solver_options options={.relative_tolerance=1.0e-6,
      .absolute_tolerance=1.0e-12,.absolute_tolerances=atol,
      .absolute_tolerance_count=3,.max_steps=10000};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  double y[]={0.9e-6,0.0,0.1e-6};
  int status=chem_rosenbrock_solve_constrained_with_parameters(&mechanism,
      &workspace,y,298.0,101325.0,40.0,NULL,30.0,&parameters,
      &set,&options,NULL);
  const double sum=y[0]+y[1]+y[2];
  if (!status && (fabs(sum-1.0e-6)>1.0e-12 || y[2]<-1.0e-18)) {
    chem_solver_workspace_destroy(&workspace);
    return 3;
  }
  if (!status && (y[0]>=1.0e-12 || fabs(y[1]-0.9e-6)>1.0e-10 ||
                  fabs(y[2]-0.1e-6)>1.0e-10)) status=4;
  chem_solver_workspace_destroy(&workspace);
  return status;
}

static int check_all_constrained_overshoot_orders(void) {
  const chem_rosenbrock_parameter_set schemes[]={
      CHEM_ROSENBROCK_THREE_STAGE,CHEM_ROSENBROCK_FOUR_STAGE,
      CHEM_ROSENBROCK_FOUR_STAGE_DAE,CHEM_ROSENBROCK_SIX_STAGE_DAE};
  for (size_t i=0;i<sizeof(schemes)/sizeof(schemes[0]);++i)
    if (check_constrained_overshoot_order(schemes[i])) return (int)i+1;
  return 0;
}

static int check_equilibrium_conservation_overshoot(void) {
  const chem_species species[]={{.name="A_gas"},{.name="A_aq"},{.name="P"}};
  const chem_term reactant[]={{1,1.0}}, product[]={{2,1.0}};
  const chem_reaction reaction={.kind=CHEM_ARRHENIUS,.reactants=reactant,
      .reactant_count=1,.products=product,.product_count=1,.p={1.0e3,0,0,298,0}};
  const chem_mechanism mechanism={.species=species,.species_count=3,
      .reactions=&reaction,.reaction_count=1};
  const chem_constraint_term eq_r[]={{0,1.0}}, eq_p[]={{1,1.0}};
  const chem_constraint_term mass[]={{0,1.0},{1,1.0},{2,1.0}};
  const chem_constraint constraints[]={
    {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,.name="gas_aq_eq",
      .algebraic_species=1,.reactants=eq_r,.reactant_count=1,
      .products=eq_p,.product_count=1,.equilibrium_constant_ref=10.0,
      .reference_temperature=298.0},
    {.kind=CHEM_CONSTRAINT_LINEAR,.name="mass_conservation",
      .algebraic_species=0,.terms=mass,.term_count=3,.constant=1.0e-6}};
  const chem_constraint_set set={.constraints=constraints,
      .constraint_count=2,.species_count=3};
  const chem_rosenbrock_parameters parameters=
      chem_rosenbrock_builtin_parameters(CHEM_ROSENBROCK_THREE_STAGE);
  const double atol[]={1.0e-8,1.0e-8,1.0e-12};
  const chem_solver_options options={.relative_tolerance=1.0e-6,
      .absolute_tolerance=1.0e-12,.absolute_tolerances=atol,
      .absolute_tolerance_count=3,.max_steps=10000};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  double y[]={1.0e-6/11.0,10.0e-6/11.0,0.0};
  int status=chem_rosenbrock_solve_constrained_with_parameters(&mechanism,
      &workspace,y,298.0,101325.0,40.0,NULL,30.0,&parameters,&set,&options,NULL);
  if (!status && (fabs(y[0]+y[1]+y[2]-1.0e-6)>1.0e-12 ||
      y[0]<-1.0e-18 || y[1]<-1.0e-18 || fabs(y[2]-1.0e-6)>1.0e-8)) status=2;
  chem_solver_workspace_destroy(&workspace);
  return status;
}

static int check_process_algebraic_mask(void) {
  const chem_species species[] = {
      {.name="A"},{.name="B"},{.name="C"},{.name="D"}};
  const chem_term reactants[] = {{.species=0,.coefficient=1.0},
                                  {.species=1,.coefficient=1.0}};
  const chem_term products[] = {{.species=2,.coefficient=1.0},
                                 {.species=3,.coefficient=1.0}};
  const chem_reaction reaction = {.kind=CHEM_ARRHENIUS,.reactants=reactants,
      .reactant_count=2,.products=products,.product_count=2};
  const chem_mechanism mechanism = {.species=species,.species_count=4,
      .reactions=&reaction,.reaction_count=1};
  const unsigned char mask[] = {0,1,0,1};
  const double rates[] = {10.0};
  const double states[][4] = {{1.0,2.0,3.0,4.0},{5.0,6.0,7.0,8.0}};
  const double expected[][4] = {{980.0,1000.0,1020.0,1000.0},
                                 {700.0,1000.0,1300.0,1000.0}};
  for (size_t cell = 0; cell < 2; ++cell) {
    double forcing[] = {1000.0,1000.0,1000.0,1000.0};
    if (chem_forcing_add_with_algebraic_mask(&mechanism,states[cell],rates,
          40.0,mask,forcing)) return 1;
    for (size_t row = 0; row < 4; ++row)
      if (fabs(forcing[row]-expected[cell][row])>1.0e-12) return 2;
    double jacobian[] = {500,500,500,500,500,500,500,500,
                         500,500,500,500,500,500,500,500};
    double unmasked[16] = {0};
    if (chem_negative_jacobian_with_constants(&mechanism,states[cell],rates,
          unmasked) || chem_negative_jacobian_add_with_algebraic_mask(
          &mechanism,states[cell],rates,mask,jacobian)) return 3;
    for (size_t row = 0; row < 4; ++row)
      for (size_t column = 0; column < 4; ++column) {
        const double expected_value = mask[row]
            ? 500.0 : 500.0 + unmasked[row*4+column];
        if (jacobian[row*4+column] != expected_value) return 4;
      }
  }
  return 0;
}

static int check_external_constraint_equivalents(void) {
  const chem_species species[]={{.name="A_GAS"},{.name="A_AQ"}};
  const chem_term gas[]={{0,1}},aq[]={{1,1}};
  const chem_reaction reaction={.kind=CHEM_ARRHENIUS,.reactants=gas,
      .reactant_count=1,.products=aq,.product_count=1,.p={0.1,0,0,300,0}};
  const chem_mechanism mechanism={.species=species,.species_count=2,
      .reactions=&reaction,.reaction_count=1};
  const chem_constraint_term mass_terms[]={{0,1},{1,1}};
  const chem_constraint mass={.kind=CHEM_CONSTRAINT_LINEAR,.name="mass",
      .algebraic_species=1,.terms=mass_terms,.term_count=2,.constant=1.0};
  const chem_constraint_set set={.constraints=&mass,.constraint_count=1,
      .species_count=2};
  unsigned char algebraic[2]={0};
  if (chem_constraint_set_validate(&set) ||
      chem_constraint_set_algebraic_mask(&set,algebraic) ||
      algebraic[0] || !algebraic[1]) return 1;

  /* A process-only model is the same mechanism without its constraint set. */
  chem_solver_workspace ode_workspace={0};
  if (chem_solver_workspace_init(&ode_workspace,&mechanism)) return 2;
  double ode_state[]={1.0,0.0};
  const chem_solver_options ode_options=chem_solver_default_options();
  const int ode_status=chem_rosenbrock3_solve_with_workspace(&mechanism,
      &ode_workspace,ode_state,298,101325,40,NULL,10,&ode_options,NULL);
  chem_solver_workspace_destroy(&ode_workspace);
  if (ode_status || !(ode_state[0]<1.0) || !(ode_state[1]>0.0) ||
      fabs(ode_state[0]+ode_state[1]-1.0)>1.0e-10) return 3;

  /* The algebraic row contains constraint forcing/Jacobian only. */
  const unsigned char mask[]={0,1};
  const double y[]={0.7,0.2},rates[]={10.0};
  double forcing[]={0,0},negative_jacobian[4]={0};
  if (chem_forcing_add_with_algebraic_mask(&mechanism,y,rates,40,mask,forcing) ||
      chem_constraint_set_add_forcing(&set,y,NULL,1,forcing) ||
      chem_negative_jacobian_add_with_algebraic_mask(&mechanism,y,rates,mask,
                                                       negative_jacobian) ||
      chem_constraint_set_add_negative_jacobian(&set,y,NULL,1,negative_jacobian))
    return 4;
  if (fabs(forcing[0]+7.0)>1e-12 || fabs(forcing[1]+0.1)>1e-12 ||
      fabs(negative_jacobian[0]-10.0)>1e-12 ||
      fabs(negative_jacobian[2]+1.0)>1e-12 ||
      fabs(negative_jacobian[3]+1.0)>1e-12) return 5;

  /* Constraint-only and process-plus-constraint solves both preserve mass. */
  chem_solver_workspace dae_workspace={0};
  if (chem_solver_workspace_init(&dae_workspace,&mechanism)) return 6;
  const chem_rosenbrock_parameters dae=chem_rosenbrock_builtin_parameters(
      CHEM_ROSENBROCK_FOUR_STAGE_DAE);
  double y_dae[]={0.9,0.1};
  for (size_t step=0;step<20;++step) {
    if (chem_rosenbrock_solve_constrained_with_parameters(&mechanism,
          &dae_workspace,y_dae,298,101325,40,NULL,10,&dae,&set,NULL,NULL) ||
        fabs(y_dae[0]+y_dae[1]-1.0)>1e-10) {
      chem_solver_workspace_destroy(&dae_workspace);
      return 7;
    }
  }
  chem_solver_workspace_destroy(&dae_workspace);
  chem_mechanism constraint_only={.species=species,.species_count=2};
  if (chem_solver_workspace_init(&dae_workspace,&constraint_only)) return 8;
  double only_state[]={0.8,0.2};
  const int only_status=chem_rosenbrock_solve_constrained_with_parameters(
      &constraint_only,&dae_workspace,only_state,298,101325,40,NULL,50,&dae,
      &set,NULL,NULL);
  chem_solver_workspace_destroy(&dae_workspace);
  if (only_status || fabs(only_state[0]+only_state[1]-1.0)>1e-12) return 9;

  /* Three independent cells follow the same conservation constraint. */
  chem_solver_workspace group_workspace={0};
  if (chem_solver_workspace_init_group(&group_workspace,&mechanism,3)) return 10;
  double group_state[]={0.9,0.1,0.7,0.3,0.5,0.5};
  const double temperatures[]={298,298,298},pressures[]={101325,101325,101325};
  const double densities[]={40,40,40};
  const int group_status=chem_rosenbrock_solve_constrained_group_with_parameters(
      &mechanism,&group_workspace,group_state,temperatures,pressures,densities,
      NULL,3,10,&dae,&set,NULL,NULL);
  chem_solver_workspace_destroy(&group_workspace);
  if (group_status) return 11;
  for (size_t cell=0;cell<3;++cell)
    if (fabs(group_state[2*cell]+group_state[2*cell+1]-1.0)>1e-10) return 12;

  /* A built-in equilibrium and a second model constraint merge into one set. */
  const chem_species combined_species[]={{.name="A_GAS"},{.name="A_AQ"},
      {.name="B"},{.name="C"}};
  const chem_constraint_term combined_mass[]={{0,1},{1,1}};
  const chem_constraint_term eq_reactant[]={{2,1}},eq_product[]={{3,1}};
  const chem_constraint combined[]={
    {.kind=CHEM_CONSTRAINT_LINEAR,.algebraic_species=1,
      .terms=combined_mass,.term_count=2,.constant=1.0},
    {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,.algebraic_species=3,
      .reactants=eq_reactant,.reactant_count=1,.products=eq_product,
      .product_count=1,.equilibrium_constant_ref=5.0,
      .reference_temperature=298}};
  const chem_constraint_set combined_set={.constraints=combined,
      .constraint_count=2,.species_count=4};
  unsigned char combined_mask[4]={0};
  if (sizeof(combined_species)/sizeof(combined_species[0]) != combined_set.species_count ||
      chem_constraint_set_validate(&combined_set) ||
      chem_constraint_set_algebraic_mask(&combined_set,combined_mask) ||
      !combined_mask[1] || !combined_mask[3] || combined_mask[0] ||
      combined_mask[2]) return 13;
  return 0;
}

static int check_external_equilibrium_convergence(void) {
  const chem_species species[]={{.name="A"},{.name="B"},{.name="C"}};
  const chem_term a[]={{0,1}},b[]={{1,1}},c[]={{2,1}};
  const chem_term to_b[]={{1,1}},to_c[]={{2,1}},to_b_product[]={{1,1}};
  const chem_reaction kinetic_reactions[]={
    {.kind=CHEM_ARRHENIUS,.reactants=a,.reactant_count=1,.products=to_b,
      .product_count=1,.p={0.1,0,0,300,0}},
    {.kind=CHEM_ARRHENIUS,.reactants=b,.reactant_count=1,.products=to_c,
      .product_count=1,.p={10,0,0,300,0}},
    {.kind=CHEM_ARRHENIUS,.reactants=c,.reactant_count=1,.products=to_b_product,
      .product_count=1,.p={2,0,0,300,0}}};
  const chem_mechanism kinetic={.species=species,.species_count=3,
      .reactions=kinetic_reactions,.reaction_count=3};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&kinetic)) return 1;
  double kinetic_y[]={1,0,0};
  for (size_t step=0;step<200;++step)
    if (chem_rosenbrock3_solve_with_workspace(&kinetic,&workspace,kinetic_y,
          298,101325,40,NULL,1,NULL,NULL)) {
      chem_solver_workspace_destroy(&workspace); return 2;
    }
  chem_solver_workspace_destroy(&workspace);

  const chem_reaction drive={.kind=CHEM_ARRHENIUS,.reactants=a,
      .reactant_count=1,.products=to_b,.product_count=1,.p={0.1,0,0,300,0}};
  const chem_mechanism dae_mechanism={.species=species,.species_count=3,
      .reactions=&drive,.reaction_count=1};
  const chem_constraint_term eq_r[]={{1,1}},eq_p[]={{2,1}},
      mass_terms[]={{0,1},{1,1},{2,1}};
  const chem_constraint conservative_constraints[]={
    {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,.algebraic_species=2,
      .reactants=eq_r,.reactant_count=1,.products=eq_p,.product_count=1,
      .equilibrium_constant_ref=5,.reference_temperature=298},
    {.kind=CHEM_CONSTRAINT_LINEAR,.algebraic_species=1,.terms=mass_terms,
      .term_count=3,.constant=1}};
  const chem_constraint_set conservative_set={.constraints=conservative_constraints,
      .constraint_count=2,.species_count=3};
  const chem_rosenbrock_parameters dae=chem_rosenbrock_builtin_parameters(
      CHEM_ROSENBROCK_FOUR_STAGE_DAE);
  if (chem_solver_workspace_init(&workspace,&dae_mechanism)) return 3;
  double conservative_y[]={1,0,0};
  for (size_t step=0;step<200;++step)
    if (chem_rosenbrock_solve_constrained_with_parameters(&dae_mechanism,
          &workspace,conservative_y,298,101325,40,NULL,1,&dae,
          &conservative_set,NULL,NULL)) {
      chem_solver_workspace_destroy(&workspace); return 4;
    }
  chem_solver_workspace_destroy(&workspace);

  const chem_constraint simple_constraint={.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
      .algebraic_species=2,.reactants=eq_r,.reactant_count=1,
      .products=eq_p,.product_count=1,.equilibrium_constant_ref=5,
      .reference_temperature=298};
  const chem_constraint_set simple_set={.constraints=&simple_constraint,
      .constraint_count=1,.species_count=3};
  if (chem_solver_workspace_init(&workspace,&dae_mechanism)) return 5;
  double simple_y[]={1,0,0};
  for (size_t step=0;step<200;++step)
    if (chem_rosenbrock_solve_constrained_with_parameters(&dae_mechanism,
          &workspace,simple_y,298,101325,40,NULL,1,&dae,&simple_set,NULL,NULL)) {
      chem_solver_workspace_destroy(&workspace); return 6;
    }
  chem_solver_workspace_destroy(&workspace);

  const double expected_b=1.0/6.0,expected_c=5.0/6.0;
  if (fabs(kinetic_y[0])>1e-6 || fabs(conservative_y[0])>1e-6 ||
      fabs(kinetic_y[1]-expected_b)>1e-3 ||
      fabs(conservative_y[1]-expected_b)>1e-3 ||
      fabs(kinetic_y[2]-expected_c)>1e-3 ||
      fabs(conservative_y[2]-expected_c)>1e-3 ||
      fabs(kinetic_y[0]+kinetic_y[1]+kinetic_y[2]-1.0)>1e-6 ||
      fabs(conservative_y[0]+conservative_y[1]+conservative_y[2]-1.0)>1e-3 ||
      fabs(simple_y[2]/simple_y[1]-5.0)>1e-3 ||
      fabs(simple_y[0]+simple_y[1]+simple_y[2]-6.0)>0.1) return 7;
  return 0;
}

static int check_composed_equilibrium_dae(void) {
  const chem_species species[]={{.name="A"},{.name="B"},{.name="C"},{.name="D"}};
  const chem_term a[]={{0,1}},b[]={{1,1}},c[]={{2,1}},d[]={{3,1}};
  const chem_term to_b[]={{1,1}},to_c[]={{2,1}},to_d[]={{3,1}};
  const chem_reaction kinetic_reactions[]={
    {.kind=CHEM_ARRHENIUS,.reactants=a,.reactant_count=1,.products=to_b,
      .product_count=1,.p={0.1,0,0,300,0}},
    {.kind=CHEM_ARRHENIUS,.reactants=b,.reactant_count=1,.products=to_c,
      .product_count=1,.p={10,0,0,300,0}},
    {.kind=CHEM_ARRHENIUS,.reactants=c,.reactant_count=1,.products=to_b,
      .product_count=1,.p={2,0,0,300,0}},
    {.kind=CHEM_ARRHENIUS,.reactants=b,.reactant_count=1,.products=to_d,
      .product_count=1,.p={8,0,0,300,0}},
    {.kind=CHEM_ARRHENIUS,.reactants=d,.reactant_count=1,.products=to_b,
      .product_count=1,.p={8.0/3.0,0,0,300,0}}};
  const chem_mechanism kinetic={.species=species,.species_count=4,
      .reactions=kinetic_reactions,.reaction_count=5};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&kinetic)) return 1;
  double kinetic_y[]={1,0,0,0};
  for (size_t step=0;step<200;++step)
    if (chem_rosenbrock3_solve_with_workspace(&kinetic,&workspace,kinetic_y,
          298,101325,40,NULL,1,NULL,NULL)) {
      chem_solver_workspace_destroy(&workspace); return 2;
    }
  chem_solver_workspace_destroy(&workspace);

  const chem_reaction drive={.kind=CHEM_ARRHENIUS,.reactants=a,
      .reactant_count=1,.products=to_b,.product_count=1,.p={0.1,0,0,300,0}};
  const chem_mechanism dae_mechanism={.species=species,.species_count=4,
      .reactions=&drive,.reaction_count=1};
  const chem_constraint_term br[]={{1,1}},cp[]={{2,1}},dp[]={{3,1}},
      mass[]={{0,1},{1,1},{2,1},{3,1}};
  const chem_constraint constraints[]={
    {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,.algebraic_species=2,
      .reactants=br,.reactant_count=1,.products=cp,.product_count=1,
      .equilibrium_constant_ref=5,.reference_temperature=298},
    {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,.algebraic_species=3,
      .reactants=br,.reactant_count=1,.products=dp,.product_count=1,
      .equilibrium_constant_ref=3,.reference_temperature=298},
    {.kind=CHEM_CONSTRAINT_LINEAR,.algebraic_species=1,
      .terms=mass,.term_count=4,.constant=1}};
  const chem_constraint_set set={.constraints=constraints,.constraint_count=3,
      .species_count=4};
  if (chem_solver_workspace_init(&workspace,&dae_mechanism)) return 3;
  const chem_rosenbrock_parameters dae=chem_rosenbrock_builtin_parameters(
      CHEM_ROSENBROCK_FOUR_STAGE_DAE);
  double dae_y[]={1,0,0,0};
  for (size_t step=0;step<200;++step)
    if (chem_rosenbrock_solve_constrained_with_parameters(&dae_mechanism,
          &workspace,dae_y,298,101325,40,NULL,1,&dae,&set,NULL,NULL)) {
      chem_solver_workspace_destroy(&workspace); return 4;
    }
  chem_solver_workspace_destroy(&workspace);
  const double expected_b=1.0/9.0;
  if (fabs(kinetic_y[0])>1e-6 || fabs(dae_y[0])>1e-6 ||
      fabs(kinetic_y[1]-expected_b)>1e-3 || fabs(dae_y[1]-expected_b)>1e-3 ||
      fabs(kinetic_y[2]-5*expected_b)>1e-3 || fabs(dae_y[2]-5*expected_b)>1e-3 ||
      fabs(kinetic_y[3]-3*expected_b)>1e-3 || fabs(dae_y[3]-3*expected_b)>1e-3 ||
      fabs(dae_y[0]+dae_y[1]+dae_y[2]+dae_y[3]-1)>1e-3) return 5;
  return 0;
}

static int check_temperature_dependent_equilibrium(void) {
  const chem_species species[]={{.name="A"},{.name="B"},{.name="C"}};
  const chem_term a[]={{0,1}},b[]={{1,1}};
  const chem_reaction reaction={.kind=CHEM_ARRHENIUS,.reactants=a,
      .reactant_count=1,.products=b,.product_count=1,.p={0.1,0,0,300,0}};
  const chem_mechanism mechanism={.species=species,.species_count=3,
      .reactions=&reaction,.reaction_count=1};
  const chem_constraint_term eq_r[]={{1,1}},eq_p[]={{2,1}};
  const chem_constraint eq={.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
      .algebraic_species=2,.reactants=eq_r,.reactant_count=1,
      .products=eq_p,.product_count=1,.equilibrium_constant_ref=2.0,
      .enthalpy_change=-3000.0*8.31446261815324,.reference_temperature=298.15};
  const chem_constraint_set set={.constraints=&eq,.constraint_count=1,
      .species_count=3};
  const chem_rosenbrock_parameters dae=chem_rosenbrock_builtin_parameters(
      CHEM_ROSENBROCK_FOUR_STAGE_DAE);
  const double temperatures[]={298.15,350.0};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  for (size_t i=0;i<2;++i) {
    double y[]={1,0,0};
    for (size_t step=0;step<200;++step)
      if (chem_rosenbrock_solve_constrained_with_parameters(&mechanism,
            &workspace,y,temperatures[i],101325,40,NULL,1,&dae,&set,NULL,NULL)) {
        chem_solver_workspace_destroy(&workspace); return 2;
      }
    const double expected=2.0*exp(3000.0*(1.0/298.15-1.0/temperatures[i]));
    if (!(y[1]>0) || fabs(y[2]/y[1]-expected)>1e-4) {
      chem_solver_workspace_destroy(&workspace); return 3;
    }
  }
  chem_solver_workspace_destroy(&workspace);
  return 0;
}

static int check_process_jacobian_in_constraint_rows(void) {
  const chem_species species[]={
      {.name="A_GAS"},{.name="S"},{.name="B"},{.name="A_AQ"}};
  const chem_term gas[]={{0,1}},solvent[]={{1,1}};
  const chem_term to_b[]={{2,1}},to_aq[]={{3,1}};
  const chem_reaction reactions[]={
    {.kind=CHEM_ARRHENIUS,.reactants=gas,.reactant_count=1,.products=to_b,
      .product_count=1,.p={0.1,0,0,300,0}},
    {.kind=CHEM_ARRHENIUS,.reactants=solvent,.reactant_count=1,.products=to_aq,
      .product_count=1,.p={0.2,0,0,300,0}}};
  const chem_mechanism mechanism={.species=species,.species_count=4,
      .reactions=reactions,.reaction_count=2};
  const chem_constraint_term eq_r[]={{2,1}},eq_p[]={{3,1}};
  const chem_constraint constraint={.kind=CHEM_CONSTRAINT_EQUILIBRIUM,
      .algebraic_species=3,.reactants=eq_r,.reactant_count=1,
      .products=eq_p,.product_count=1,.equilibrium_constant_ref=5,
      .reference_temperature=298};
  const chem_constraint_set set={.constraints=&constraint,.constraint_count=1,
      .species_count=4};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  const chem_rosenbrock_parameters dae=chem_rosenbrock_builtin_parameters(
      CHEM_ROSENBROCK_FOUR_STAGE_DAE);
  double y[]={1,1,0,0};
  for (size_t step=0;step<20;++step)
    if (chem_rosenbrock_solve_constrained_with_parameters(&mechanism,
          &workspace,y,298,101325,40,NULL,0.5,&dae,&set,NULL,NULL)) {
      chem_solver_workspace_destroy(&workspace); return 2;
    }
  chem_solver_workspace_destroy(&workspace);
  return y[1]>=1.0 || y[2]<=0.0 || fabs(5.0*y[2]-y[3])>1e-6 ? 3 : 0;
}

static int check_constraint_workspace_reuse(void) {
  const chem_species species[]={{.name="A"},{.name="B"},{.name="C"},{.name="D"}};
  const chem_term a[]={{0,1}};
  const chem_term to_b[]={{1,1}};
  const chem_reaction reaction={.kind=CHEM_ARRHENIUS,.reactants=a,
      .reactant_count=1,.products=to_b,.product_count=1,.p={0.1,0,0,300,0}};
  const chem_mechanism mechanism={.species=species,.species_count=4,
      .reactions=&reaction,.reaction_count=1};
  const chem_constraint_term eq1_r[]={{0,1}},eq1_p[]={{3,1}},
      eq2_r[]={{1,1}},eq2_p[]={{2,1}};
  const chem_constraint constraints[]={
    {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,.algebraic_species=3,
      .reactants=eq1_r,.reactant_count=1,.products=eq1_p,.product_count=1,
      .equilibrium_constant_ref=2,.reference_temperature=298},
    {.kind=CHEM_CONSTRAINT_EQUILIBRIUM,.algebraic_species=2,
      .reactants=eq2_r,.reactant_count=1,.products=eq2_p,.product_count=1,
      .equilibrium_constant_ref=5,.reference_temperature=298}};
  const chem_constraint_set set={.constraints=constraints,.constraint_count=2,
      .species_count=4};
  const chem_rosenbrock_parameters dae=chem_rosenbrock_builtin_parameters(
      CHEM_ROSENBROCK_FOUR_STAGE_DAE);
  chem_solver_workspace workspaces[2]={{0}};
  if (chem_solver_workspace_init(&workspaces[0],&mechanism) ||
      chem_solver_workspace_init(&workspaces[1],&mechanism)) {
    chem_solver_workspace_destroy(&workspaces[0]);
    chem_solver_workspace_destroy(&workspaces[1]);
    return 1;
  }
  double states[2][4]={{1,0,0,2},{1,0,0,2}};
  int status=0;
  for (size_t step=0;step<20 && !status;++step)
    for (size_t solver=0;solver<2 && !status;++solver)
      status=chem_rosenbrock_solve_constrained_with_parameters(&mechanism,
          &workspaces[solver],states[solver],298,101325,40,NULL,1,&dae,&set,
          NULL,NULL);
  chem_solver_workspace_destroy(&workspaces[0]);
  chem_solver_workspace_destroy(&workspaces[1]);
  if (status || memcmp(states[0],states[1],sizeof(states[0]))!=0 ||
      fabs(states[0][3]-2.0*states[0][0])>1e-6 ||
      fabs(states[0][2]-5.0*states[0][1])>1e-6) return 2;
  return 0;
}

static int check_phase_system_metadata(void) {
  const chem_species species[] = {
      {.name="foo"},{.name="bar",.density_parameterized=1,.density_factor=42.0},
      {.name="baz"}};
  chem_phase_species members[3];
  if (chem_phase_species_init(&members[0],&species[0],0,0.0) ||
      chem_phase_species_init(&members[1],&species[1],0,0.0) ||
      chem_phase_species_init(&members[2],&species[2],0,0.0)) return 1;
  if (members[0].has_diffusion_coefficient ||
      chem_phase_species_set_diffusion_coefficient(&members[0],1.23e-5) ||
      !members[0].has_diffusion_coefficient ||
      fabs(members[0].diffusion_coefficient-1.23e-5)>1.0e-18 ||
      chem_phase_species_init(&members[0],&species[0],1,2.5e-6) ||
      !members[0].has_diffusion_coefficient ||
      fabs(members[0].diffusion_coefficient-2.5e-6)>1.0e-18) return 2;
  const chem_phase phase = {.name="gas",.species=members,.species_count=3};
  if (chem_phase_state_size(&phase) != 2) return 3;
  char name[64];
  if (chem_phase_unique_name(&phase,0,name,sizeof(name)) || strcmp(name,"gas.foo") ||
      chem_phase_unique_name(&phase,1,name,sizeof(name)) != 1 ||
      chem_phase_unique_name(&phase,2,name,sizeof(name)) || strcmp(name,"gas.baz"))
    return 4;
  chem_system system;
  if (chem_system_init(&system,&phase) || chem_system_state_size(&system) != 2 ||
      chem_system_unique_name(&system,0,name,sizeof(name)) || strcmp(name,"foo") ||
      chem_system_unique_name(&system,1,name,sizeof(name)) || strcmp(name,"baz") ||
      chem_system_unique_name(&system,2,name,sizeof(name)) == 0) return 5;
  return 0;
}

static int check_species_properties(void) {
  const chem_property properties[] = {
      {.name="name [units]",.kind=CHEM_PROPERTY_DOUBLE,
       .value.real_value=1.0},
      {.name="name2 [units2]",.kind=CHEM_PROPERTY_DOUBLE,
       .value.real_value=2.0},
      {.name="label",.kind=CHEM_PROPERTY_STRING,
       .value.string_value="bar"},
      {.name="count",.kind=CHEM_PROPERTY_INTEGER,
       .value.integer_value=42},
      {.name="enabled",.kind=CHEM_PROPERTY_BOOLEAN,
       .value.boolean_value=1}};
  chem_species species = {.name="thing",.properties=properties,
      .property_count=5,.temperature_parameterization_factor=100.0};
  if (strcmp(species.name,"thing") || !chem_species_is_parameterized(&species) ||
      chem_species_parameterize_conditions(&species,12.5,0.0) != 1250.0)
    return 1;
  double real_value = 0.0;
  int integer_value = 0, boolean_value = 0;
  const char *string_value = NULL;
  if (chem_species_get_double_property(&species,"name [units]",&real_value) ||
      real_value != 1.0 ||
      chem_species_get_double_property(&species,"name2 [units2]",&real_value) ||
      real_value != 2.0 ||
      chem_species_get_string_property(&species,"label",&string_value) ||
      strcmp(string_value,"bar") ||
      chem_species_get_integer_property(&species,"count",&integer_value) ||
      integer_value != 42 ||
      chem_species_get_boolean_property(&species,"enabled",&boolean_value) ||
      boolean_value != 1) return 2;
  const chem_species copied = species;
  chem_species assigned = {.name="placeholder"};
  assigned = species;
  if (strcmp(copied.name,"thing") || strcmp(assigned.name,"thing") ||
      copied.properties != properties || assigned.properties != properties ||
      copied.property_count != 5 || assigned.property_count != 5 ||
      !chem_species_is_parameterized(&copied) ||
      !chem_species_is_parameterized(&assigned) ||
      chem_species_get_integer_property(&copied,"count",&integer_value) ||
      integer_value != 42 ||
      chem_species_get_string_property(&assigned,"label",&string_value) ||
      strcmp(string_value,"bar")) return 3;
  if (chem_species_get_double_property(&species,"missing",&real_value) !=
          CHEM_PROPERTY_STATUS_NOT_FOUND ||
      chem_species_get_integer_property(&species,"name [units]",&integer_value) !=
          CHEM_PROPERTY_STATUS_TYPE_MISMATCH) return 4;
  chem_species third_body = {.name="M"};
  if (chem_species_set_third_body(&third_body) ||
      !chem_species_is_parameterized(&third_body) ||
      chem_species_parameterize(&third_body,42.4) != 42.4) return 5;
  return 0;
}

static void matrix_add_index(size_t row, size_t column, double *value,
                             void *context) {
  const double scale = *(const double *)context;
  *value += scale * (double)(row * 10 + column);
}

static int matrix_greater_than(double value, void *context) {
  return value > *(const double *)context;
}

static int matrix_is_finite(double value, void *context) {
  (void)context;
  return isfinite(value);
}

static double matrix_count(double accumulator, double value, void *context) {
  (void)value;
  (void)context;
  return accumulator + 1.0;
}

static double matrix_sum_squares(double accumulator, double value,
                                 void *context) {
  (void)context;
  return accumulator + value*value;
}

static void matrix_add_vector_rows(size_t row, double *destination,
    size_t destination_columns, size_t destination_stride,
    const double *first, size_t first_columns, size_t first_stride,
    const double *second, size_t second_columns, size_t second_stride,
    void *context) {
  (void)row;
  const double scale = context ? *(const double *)context : 1.0;
  if (destination_columns && first_columns && second_columns)
    destination[0] = first[0] + scale*second[0];
  (void)destination_stride; (void)first_stride; (void)second_stride;
}

typedef struct {
  const double *values;
  size_t count;
} matrix_vector_context;

static void matrix_add_captured_vector(size_t row, size_t column, double *value,
                                       void *context) {
  (void)row;
  const matrix_vector_context *vector = context;
  if (column < vector->count) *value += vector->values[column];
}

static void matrix_combine_rows(size_t row, double *destination,
    size_t destination_columns, size_t destination_stride, const double *first,
    size_t first_columns, size_t first_stride, const double *second,
    size_t second_columns, size_t second_stride, void *context) {
  (void)row; (void)first_stride; (void)context;
  if (destination_columns > 1 && first_columns > 0 && second_columns > 2)
    destination[destination_stride] = first[0] + second[2*second_stride];
}

static int matrix_row_sum_three(double *row, size_t columns,
                                size_t stride, void *context) {
  if (!row || columns != 3 || !context) return 1;
  const size_t output_column = *(const size_t *)context;
  if (output_column >= columns) return 1;
  row[output_column*stride] = row[0] + row[stride] + row[2*stride];
  return 0;
}

static int check_matrix_view_api(void) {
  {
    const size_t output_column = 0;
    chem_matrix_row_plan plan;
    double dense3[] = {1,2,3, 4,5,6};
    double tall_data[21], wrong4[28] = {0}, wrong2[14] = {0};
    for (size_t i = 0; i < 21; ++i) tall_data[i] = (double)i+1;
    chem_matrix_view dense_view, bad4, bad2, tall;
    if (chem_matrix_row_plan_init(&plan,3,matrix_row_sum_three,
                                  (void *)&output_column) ||
        chem_matrix_view_init(&dense_view,dense3,2,3,3,1) ||
        chem_matrix_view_init(&bad4,wrong4,2,4,4,1) ||
        chem_matrix_view_init(&bad2,wrong2,2,2,2,1) ||
        chem_matrix_view_init(&tall,tall_data,7,3,3,1) ||
        chem_matrix_row_plan_apply(&plan,&bad4) == 0 ||
        chem_matrix_row_plan_apply(&plan,&bad2) == 0 ||
        chem_matrix_row_plan_apply(&plan,&dense_view) ||
        chem_matrix_row_plan_apply(&plan,&tall) || dense3[0] != 6.0 ||
        dense3[3] != 15.0 || tall_data[18] != 60.0) return 94;
    for (size_t width = 1; width <= 4; ++width) {
      double packed[28] = {0};
      chem_vector_matrix_view packed_view;
      if (chem_vector_matrix_view_init(&packed_view,packed,7,3,width))
        return 95;
      for (size_t row = 0; row < 7; ++row)
        for (size_t column = 0; column < 3; ++column)
          if (chem_vector_matrix_view_set(&packed_view,row,column,
              (double)(row*3+column+1))) return 96;
      if (chem_vector_matrix_row_plan_apply(&plan,&packed_view)) return 97;
      double result = 0.0;
      if (chem_vector_matrix_view_get(&packed_view,6,0,&result) || result != 60.0)
        return 98;
    }
  }
  {
    const double old_values[] = {2.0, -1.0, 0.5, 4.0, 9.0, -2.0};
    const double new_values[] = {1.0, 3.0, 0.0, 2.0, 8.0, -4.0};
    const double error_values[] = {0.2, 0.4, 0.1, 0.3, 0.5, 0.6};
    const double tolerances[] = {0.01, 0.02};
    chem_const_matrix_view old_view, new_view, error_view;
    double value = 0.0, sum = 0.0;
    if (chem_const_matrix_view_init(&old_view,old_values,3,2,2,1) ||
        chem_const_matrix_view_init(&new_view,new_values,3,2,2,1) ||
        chem_const_matrix_view_init(&error_view,error_values,3,2,2,1))
      return 87;
    for (size_t row = 0; row < 3; ++row)
      for (size_t column = 0; column < 2; ++column) {
        const double scale = tolerances[column] + 1e-3*fmax(
            fabs(old_values[row*2+column]),fabs(new_values[row*2+column]));
        const double normalized = error_values[row*2+column]/scale;
        sum += normalized*normalized;
      }
    const double expected = fmax(sqrt(sum/6.0),1e-10);
    if (chem_normalized_error(&old_view,&new_view,&error_view,tolerances,2,
                              1e-3,1e-10,&value) ||
        fabs(value-expected) > 1e-13*expected ||
        chem_normalized_error(&old_view,&new_view,&error_view,tolerances,1,
                              1e-3,1e-10,&value) == 0) return 88;
    const size_t widths[] = {1,2,4};
    for (size_t w = 0; w < sizeof(widths)/sizeof(widths[0]); ++w) {
      const size_t width = widths[w];
      double packed_old[8], packed_new[8], packed_error[8];
      chem_vector_matrix_view mutable_old, mutable_new, mutable_error;
      chem_const_vector_matrix_view old_packed, new_packed, error_packed;
      if (chem_vector_matrix_view_init(&mutable_old,packed_old,3,2,width) ||
          chem_vector_matrix_view_init(&mutable_new,packed_new,3,2,width) ||
          chem_vector_matrix_view_init(&mutable_error,packed_error,3,2,width))
        return 89;
      for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 2; ++j)
          if (chem_vector_matrix_view_set(&mutable_old,i,j,old_values[i*2+j]) ||
              chem_vector_matrix_view_set(&mutable_new,i,j,new_values[i*2+j]) ||
              chem_vector_matrix_view_set(&mutable_error,i,j,error_values[i*2+j]))
            return 90;
      if (chem_const_vector_matrix_view_init(&old_packed,packed_old,3,2,width) ||
          chem_const_vector_matrix_view_init(&new_packed,packed_new,3,2,width) ||
          chem_const_vector_matrix_view_init(&error_packed,packed_error,3,2,width) ||
          chem_vector_normalized_error(&old_packed,&new_packed,&error_packed,
                                       tolerances,2,1e-3,1e-10,&value) ||
          fabs(value-expected) > 1e-13*expected) return 91;
    }
  }
  {
    double dense_data[] = {-12.2,-24.3,-0.43,-23.4};
    chem_matrix_view dense;
    if (chem_matrix_view_init(&dense,dense_data,2,2,2,1) ||
        chem_matrix_view_add_diagonal(&dense,42.042) ||
        fabs(dense_data[0]-(42.042-12.2)) > 1e-14 ||
        fabs(dense_data[1]+24.3) > 1e-14 ||
        fabs(dense_data[2]+0.43) > 1e-14 ||
        fabs(dense_data[3]-(42.042-23.4)) > 1e-14 ||
        chem_matrix_view_add_diagonal(&dense,NAN) == 0) return 92;
    double packed_data[] = {-12.2,-0.43,-24.3,-23.4};
    chem_vector_matrix_view packed;
    if (chem_vector_matrix_view_init(&packed,packed_data,2,2,2) ||
        chem_vector_matrix_view_add_diagonal(&packed,42.042) ||
        fabs(packed_data[0]-(42.042-12.2)) > 1e-14 ||
        fabs(packed_data[1]+0.43) > 1e-14 ||
        fabs(packed_data[2]+24.3) > 1e-14 ||
        fabs(packed_data[3]-(42.042-23.4)) > 1e-14) return 93;
  }
  double storage[] = {1.0, 99.0, 2.0, 99.0, 3.0, 99.0,
                      4.0, 99.0, 5.0, 99.0, 6.0, 99.0};
  chem_matrix_view view;
  if (chem_matrix_view_init(&view, storage, 2, 3, 6, 2)) return 1;
  double value = 0.0;
  if (chem_matrix_view_get(&view, 1, 2, &value) || value != 6.0 ||
      chem_matrix_view_set(&view, 0, 1, -2.0) || storage[2] != -2.0) return 2;
  const double expected[] = {1.0, -2.0, 3.0, 4.0, 5.0, 6.0};
  double flat[6];
  if (chem_matrix_view_copy_to_vector(&view, flat, 6)) return 3;
  for (size_t i = 0; i < 6; ++i) if (flat[i] != expected[i]) return 4;
  const double replacement[] = {6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
  if (chem_matrix_view_copy_from_vector(&view, replacement, 6) ||
      chem_matrix_view_copy_from_vector(&view, replacement, 5) == 0) return 5;
  double other_data[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  chem_matrix_view other;
  if (chem_matrix_view_init(&other, other_data, 2, 3, 3, 1) ||
      chem_matrix_view_axpy(&view, &other, 2.0) ||
      chem_matrix_view_axpy(&view, &other, NAN) == 0) return 6;
  if (chem_matrix_view_min(&view, &value) || value != 3.0) return 7;
  if (chem_matrix_view_max(&view, &value) || value != 8.0) return 8;
  if (chem_matrix_view_sum(&view, &value) || value != 33.0) return 9;
  if (chem_matrix_view_fill(&view, 2.5)) return 10;
  if (chem_matrix_view_sum(&view, &value) || value != 15.0) return 11;
  const double scale = 1.0;
  if (chem_matrix_view_for_each(&view, matrix_add_index, (void *)&scale) ||
      chem_matrix_view_get(&view, 1, 2, &value) || value != 14.5) return 12;
  if (chem_matrix_view_get(&view, 2, 0, &value) == 0 ||
      chem_matrix_view_init(&view, storage, (size_t)-1, 2, (size_t)-1, 1) == 0)
    return 13;
  chem_matrix_view wrong_rows, wrong_columns;
  double mismatch_storage[8] = {0};
  if (chem_matrix_view_init(&wrong_rows, mismatch_storage, 3, 3, 3, 1) ||
      chem_matrix_view_init(&wrong_columns, mismatch_storage, 2, 2, 2, 1) ||
      chem_matrix_view_axpy(&view, &wrong_rows, 1.0) == 0 ||
      chem_matrix_view_axpy(&view, &wrong_columns, 1.0) == 0) return 14;
  chem_matrix_view empty;
  if (chem_matrix_view_init(&empty, NULL, 0, 0, 0, 0) ||
      chem_matrix_view_fill(&empty, 0.0) ||
      chem_matrix_view_sum(&empty, &value) == 0 ||
      chem_matrix_view_for_each(&empty, matrix_add_index, (void *)&scale)) return 15;
  double row_storage[12] = {0};
  chem_matrix_view rows;
  const double row_values[] = {12.3,15.1,24.3};
  double row_copy[3];
  if (chem_matrix_view_init(&rows,row_storage,4,3,3,1) ||
      chem_matrix_view_copy_row_from_vector(&rows,2,row_values,3) ||
      chem_matrix_view_copy_row_to_vector(&rows,2,row_copy,3) ||
      chem_matrix_view_copy_row_from_vector(&rows,2,row_values,2) == 0)
    return 16;
  for (size_t j = 0; j < 3; ++j)
    if (row_copy[j] != row_values[j] || row_storage[6+j] != row_values[j]) return 17;
  if (row_storage[0] != 0.0 || chem_matrix_view_copy_row_to_vector(&rows,4,row_copy,3)==0)
    return 18;
  double small_data[15] = {0};
  chem_matrix_view small;
  if (chem_matrix_view_init(&small,small_data,3,5,5,1) ||
      chem_matrix_view_set(&small,1,3,64.7) ||
      chem_matrix_view_set(&small,0,0,41.2) ||
      chem_matrix_view_set(&small,2,4,102.3) ||
      small_data[0] != 41.2 || small_data[8] != 64.7 || small_data[14] != 102.3)
    return 19;
  chem_const_matrix_view const_small;
  if (chem_const_matrix_view_init(&const_small,small_data,3,5,5,1) ||
      chem_const_matrix_view_get(&const_small,1,3,&value) || value != 64.7 ||
      chem_const_matrix_view_get(&const_small,2,4,&value) || value != 102.3)
    return 20;
  double initialized_data[6] = {0};
  chem_matrix_view initialized;
  if (chem_matrix_view_init(&initialized,initialized_data,2,3,3,1) ||
      chem_matrix_view_fill(&initialized,12.4) ||
      initialized_data[0] != 12.4 || initialized_data[3] != 12.4 ||
      initialized_data[5] != 12.4) return 21;
  chem_const_matrix_view const_initialized;
  if (chem_const_matrix_view_init(&const_initialized,initialized_data,2,3,3,1) ||
      chem_const_matrix_view_get(&const_initialized,1,2,&value) || value != 12.4)
    return 22;
  double loop_data[12] = {0};
  chem_matrix_view loop;
  if (chem_matrix_view_init(&loop,loop_data,3,4,4,1)) return 23;
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 4; ++j)
      if (chem_matrix_view_set(&loop,i,j,(double)(i*100+j))) return 24;
  if (chem_matrix_view_get(&loop,0,0,&value) || value != 0.0 ||
      chem_matrix_view_get(&loop,1,2,&value) || value != 102.0 ||
      chem_matrix_view_get(&loop,2,3,&value) || value != 203.0) return 25;
  chem_const_matrix_view const_loop;
  if (chem_const_matrix_view_init(&const_loop,loop_data,3,4,4,1) ||
      chem_const_matrix_view_get(&const_loop,1,2,&value) || value != 102.0 ||
      chem_const_matrix_view_get(&const_loop,2,3,&value) || value != 203.0)
    return 26;
  double const_row[4];
  if (chem_const_matrix_view_copy_row_to_vector(&const_loop,1,const_row,4) ||
      const_row[0] != 100.0 || const_row[3] != 103.0) return 27;
  if (chem_matrix_view_copy_row_to_vector(&rows,2,row_copy,2) == 0 ||
      chem_matrix_view_copy_row_to_vector(&rows,2,row_copy,4) == 0) return 28;
  double column_data[] = {1,2,3,4,5,6};
  chem_matrix_view columns;
  double column_copy[2];
  const double column_values[] = {8.0,9.0};
  if (chem_matrix_view_init(&columns,column_data,2,3,3,1) ||
      chem_matrix_view_copy_column_to_vector(&columns,1,column_copy,2) ||
      column_copy[0] != 2.0 || column_copy[1] != 5.0 ||
      chem_matrix_view_copy_column_from_vector(&columns,1,column_values,2) ||
      column_data[1] != 8.0 || column_data[4] != 9.0 ||
      chem_matrix_view_copy_column_to_vector(&columns,1,column_copy,1) == 0)
    return 29;
  if (chem_matrix_view_copy_to_vector(&columns,NULL,0) == 0 ||
      chem_matrix_view_copy_from_vector(&columns,column_values,2) == 0 ||
      chem_matrix_view_copy_to_vector(&empty,NULL,0) ||
      chem_matrix_view_copy_from_vector(&empty,NULL,0)) return 30;
  if (chem_matrix_view_max(&columns,&value) || value != 9.0) return 31;
  double reused_a[15] = {0}, reused_b[9] = {0};
  chem_matrix_view reuse_a, reuse_b;
  if (chem_matrix_view_init(&reuse_a,reused_a,5,3,3,1) ||
      chem_matrix_view_init(&reuse_b,reused_b,3,3,3,1) ||
      chem_matrix_view_for_each(&reuse_a,matrix_add_index,(void *)&scale) ||
      chem_matrix_view_for_each(&reuse_b,matrix_add_index,(void *)&scale) ||
      chem_matrix_view_get(&reuse_a,4,2,&value) || value != 42.0 ||
      chem_matrix_view_get(&reuse_b,2,1,&value) || value != 21.0) return 32;
  double iter_data[12] = {0};
  chem_matrix_view iter;
  if (chem_matrix_view_init(&iter,iter_data,3,4,4,1)) return 33;
  double row_four[4];
  if (chem_matrix_view_copy_row_to_vector(&loop,1,row_four,4)) return 34;
  for (size_t j = 0; j < 4; ++j) row_four[j] = 42.0+(double)j;
  if (chem_matrix_view_copy_row_from_vector(&iter,1,row_four,4) ||
      iter_data[4] != 42.0 || iter_data[7] != 45.0) return 35;
  const double captured_values[] = {1.0,2.0,3.0,4.0};
  const matrix_vector_context captured = {captured_values,4};
  if (chem_matrix_view_for_each(&iter,matrix_add_captured_vector,
                                (void *)&captured) ||
      iter_data[4] != 43.0 || iter_data[7] != 49.0) return 36;
  double vectorized_data[20] = {0};
  chem_vector_matrix_view vectorized;
  size_t vectorized_size = 0;
  if (chem_vector_matrix_view_storage_size(3,5,2,&vectorized_size) ||
      vectorized_size != 20 ||
      chem_vector_matrix_view_init(&vectorized,vectorized_data,3,5,2) ||
      chem_vector_matrix_view_set(&vectorized,1,3,64.7) ||
      chem_vector_matrix_view_set(&vectorized,0,0,41.2) ||
      chem_vector_matrix_view_set(&vectorized,2,4,102.3) ||
      vectorized_data[0] != 41.2 || vectorized_data[7] != 64.7 ||
      vectorized_data[18] != 102.3 ||
      chem_vector_matrix_view_group_size(&vectorized) != 10 ||
      chem_vector_matrix_view_number_of_groups(&vectorized) != 2 ||
      chem_vector_matrix_view_row_stride(&vectorized) != 1 ||
      chem_vector_matrix_view_column_stride(&vectorized) != 2) return 37;
  chem_const_vector_matrix_view const_vectorized;
  if (chem_const_vector_matrix_view_init(&const_vectorized,vectorized_data,3,5,2) ||
      chem_const_vector_matrix_view_get(&const_vectorized,1,3,&value) ||
      value != 64.7) return 38;
  double initialize_vectorized_data[6] = {0};
  chem_vector_matrix_view initialize_vectorized;
  if (chem_vector_matrix_view_init(&initialize_vectorized,
          initialize_vectorized_data,2,3,1) ||
      chem_vector_matrix_view_fill(&initialize_vectorized,12.4) ||
      initialize_vectorized_data[0] != 12.4 ||
      initialize_vectorized_data[3] != 12.4 ||
      initialize_vectorized_data[5] != 12.4) return 39;
  double vectorized_loop_data[20] = {0};
  chem_vector_matrix_view vectorized_loop;
  if (chem_vector_matrix_view_init(&vectorized_loop,vectorized_loop_data,3,4,4))
    return 40;
  for (size_t row = 0; row < 3; ++row)
    for (size_t column = 0; column < 4; ++column)
      if (chem_vector_matrix_view_set(&vectorized_loop,row,column,
                                      (double)(row*100+column))) return 41;
  if (chem_vector_matrix_view_get(&vectorized_loop,2,3,&value) || value != 203.0 ||
      vectorized_loop_data[14] != 203.0 ||
      chem_vector_matrix_view_number_of_groups(&vectorized_loop) != 1) return 42;
  chem_const_vector_matrix_view const_vectorized_loop;
  if (chem_const_vector_matrix_view_init(&const_vectorized_loop,
          vectorized_loop_data,3,4,4) ||
      chem_const_vector_matrix_view_get(&const_vectorized_loop,1,2,&value) ||
      value != 102.0) return 43;
  double operations_data[12] = {0};
  chem_vector_matrix_view operations;
  const double flat_values[] = {1,2,3,4,5,6};
  double flat_copy[6];
  if (chem_vector_matrix_view_init(&operations,operations_data,2,3,2) ||
      chem_vector_matrix_view_copy_from_vector(&operations,flat_values,6) ||
      chem_vector_matrix_view_copy_to_vector(&operations,flat_copy,6)) return 44;
  for (size_t i = 0; i < 6; ++i) if (flat_copy[i] != flat_values[i]) return 45;
  chem_const_vector_matrix_view const_operations;
  if (chem_const_vector_matrix_view_init(&const_operations,operations_data,2,3,2) ||
      chem_const_vector_matrix_view_copy_to_vector(&const_operations,flat_copy,6))
    return 46;
  double access_data[12] = {0};
  chem_vector_matrix_view access;
  const double access_flat[] = {1,2,3,4,5,6};
  const double access_row[] = {7,8,9};
  const double access_column[] = {10,11};
  double row_result[3], column_result[2];
  if (chem_vector_matrix_view_init(&access,access_data,2,3,2) ||
      chem_vector_matrix_view_copy_from_vector(&access,access_flat,6) ||
      chem_vector_matrix_view_copy_row_to_vector(&access,1,row_result,3) ||
      row_result[0] != 4.0 || row_result[2] != 6.0 ||
      chem_vector_matrix_view_copy_row_from_vector(&access,1,access_row,3) ||
      chem_vector_matrix_view_copy_column_to_vector(&access,2,column_result,2) ||
      column_result[0] != 3.0 || column_result[1] != 9.0 ||
      chem_vector_matrix_view_copy_column_from_vector(&access,2,access_column,2) ||
      access_data[4] != 10.0 || access_data[5] != 11.0 ||
      chem_vector_matrix_view_copy_row_to_vector(&access,0,row_result,2) == 0 ||
      chem_vector_matrix_view_copy_column_to_vector(&access,0,column_result,3) == 0)
    return 53;
  if (chem_vector_matrix_view_copy_column_to_vector(&access,0,column_result,2) ||
      chem_vector_matrix_view_copy_column_to_vector(&access,1,column_result,2) ||
      column_result[0] != 2.0 || column_result[1] != 8.0) return 54;
  chem_vector_matrix_view empty_vectorized;
  if (chem_vector_matrix_view_init(&empty_vectorized,NULL,0,0,2) ||
      chem_vector_matrix_view_copy_to_vector(&empty_vectorized,NULL,0) ||
      chem_vector_matrix_view_copy_from_vector(&empty_vectorized,NULL,0) ||
      chem_vector_matrix_view_copy_to_vector(&access,NULL,0) == 0 ||
      chem_vector_matrix_view_copy_from_vector(&access,NULL,0) == 0) return 55;
  double other_operations_data[6];
  chem_vector_matrix_view other_operations;
  if (chem_vector_matrix_view_init(&other_operations,other_operations_data,2,3,1) ||
      chem_vector_matrix_view_fill(&other_operations,10.0) ||
      chem_vector_matrix_view_axpy(&operations,&other_operations,2.0) ||
      chem_vector_matrix_view_axpy(&operations,&other_operations,NAN) == 0 ||
      chem_vector_matrix_view_min(&operations,&value) || value != 21.0 ||
      chem_vector_matrix_view_max(&operations,&value) || value != 26.0) return 47;
  if (chem_vector_matrix_view_for_each(&operations,matrix_add_index,
                                       (void *)&scale) ||
      chem_vector_matrix_view_get(&operations,1,2,&value) || value != 38.0)
    return 48;
  const double captured_vector_values[] = {1.0,2.0,3.0};
  const matrix_vector_context captured_vector = {captured_vector_values,3};
  if (chem_vector_matrix_view_for_each(&operations,matrix_add_captured_vector,
                                       (void *)&captured_vector) ||
      chem_vector_matrix_view_get(&operations,1,2,&value) || value != 41.0 ||
      chem_vector_matrix_view_max(&operations,&value) || value != 41.0) return 50;
  double copied_data[12] = {0};
  chem_vector_matrix_view copied;
  if (chem_vector_matrix_view_init(&copied,copied_data,2,3,4) ||
      chem_vector_matrix_view_copy(&copied,&operations) ||
      chem_vector_matrix_view_get(&copied,1,2,&value) || value != 41.0 ||
      chem_vector_matrix_view_copy(&copied,&vectorized_loop) == 0) return 56;
  double wrong_column_data[4] = {0};
  chem_vector_matrix_view wrong_column_count;
  if (chem_vector_matrix_view_init(&wrong_column_count,wrong_column_data,2,2,2) ||
      chem_vector_matrix_view_axpy(&operations,&wrong_column_count,1.0) == 0)
    return 81;
  if (chem_vector_matrix_view_for_each(&vectorized_loop,matrix_add_index,
                                       (void *)&scale) ||
      chem_vector_matrix_view_get(&vectorized_loop,2,3,&value) || value != 226.0)
    return 51;
  double reused_group_data[12] = {0};
  chem_vector_matrix_view reused_group;
  if (chem_vector_matrix_view_init(&reused_group,reused_group_data,5,2,2) ||
      chem_vector_matrix_view_for_each(&reused_group,matrix_add_index,
                                       (void *)&scale) ||
      chem_vector_matrix_view_get(&reused_group,4,1,&value) || value != 41.0)
    return 52;
  if (chem_vector_matrix_view_copy_to_vector(&operations,flat_copy,5) == 0 ||
      chem_vector_matrix_view_axpy(&operations,&vectorized_loop,1.0) == 0)
    return 49;
  double left_data[] = {0,10,1,11,2,12};
  double right_data[] = {0,20,0,2,22,4,4,24,8};
  chem_matrix_view left, right;
  if (chem_matrix_view_init(&left,left_data,3,2,2,1) ||
      chem_matrix_view_init(&right,right_data,3,3,3,1) ||
      chem_matrix_view_for_each_row(&left,&left,&right,matrix_combine_rows,NULL) ||
      left_data[1] != 0.0 || left_data[3] != 5.0 || left_data[5] != 10.0)
    return 57;
  chem_matrix_view short_matrix;
  if (chem_matrix_view_init(&short_matrix,left_data,2,3,3,1) ||
      chem_matrix_view_for_each_row(&left,&left,&short_matrix,
                                    matrix_combine_rows,NULL) == 0) return 58;
  double vector_left_data[8] = {0}, vector_right_data[12] = {0};
  chem_vector_matrix_view vector_left, vector_right;
  if (chem_vector_matrix_view_init(&vector_left,vector_left_data,3,2,2) ||
      chem_vector_matrix_view_init(&vector_right,vector_right_data,3,3,4))
    return 59;
  for (size_t row = 0; row < 3; ++row) {
    if (chem_vector_matrix_view_set(&vector_left,row,0,(double)row) ||
        chem_vector_matrix_view_set(&vector_right,row,2,(double)(row*4))) return 60;
  }
  if (chem_vector_matrix_view_for_each_row(&vector_left,&vector_left,
        &vector_right,matrix_combine_rows,NULL) ||
      chem_vector_matrix_view_get(&vector_left,0,1,&value) || value != 0.0 ||
      chem_vector_matrix_view_get(&vector_left,1,1,&value) || value != 5.0 ||
      chem_vector_matrix_view_get(&vector_left,2,1,&value) || value != 10.0)
    return 61;
  double vector_short_data[4] = {0};
  chem_vector_matrix_view vector_short;
  if (chem_vector_matrix_view_init(&vector_short,vector_short_data,2,3,2) ||
      chem_vector_matrix_view_for_each_row(&vector_left,&vector_left,
        &vector_short,matrix_combine_rows,NULL) == 0) return 62;
  chem_const_matrix_view const_loop_view;
  int logical = 0;
  const double threshold = 150.0;
  if (chem_const_matrix_view_init(&const_loop_view,loop_data,3,4,4,1) ||
      chem_matrix_view_reduce_or(&loop,matrix_greater_than,
                                 (void *)&threshold,0,&logical) || logical != 1 ||
      chem_matrix_view_reduce_or(&loop,matrix_greater_than,
                                 (void *)&threshold,1,&logical) || logical != 1 ||
      chem_matrix_view_reduce_and(&loop,matrix_is_finite,NULL,1,&logical) ||
      logical != 1 ||
      chem_matrix_view_reduce_strict(&const_loop_view,matrix_count,NULL,0.0,
                                     &value) || value != 12.0)
    return 63;
  int dense_integer_values[12];
  for (size_t i = 0; i < 12; ++i) dense_integer_values[i] = (int)i+1;
  if (chem_matrix_view_scale_from_int_vector(&loop,dense_integer_values,12,
                                              0.5) ||
      chem_matrix_view_get(&loop,0,0,&value) || value != 0.5 ||
      chem_matrix_view_get(&loop,2,3,&value) || value != 6.0 ||
      chem_matrix_view_scale_from_int_vector(&loop,dense_integer_values,11,
                                              0.5) == 0) return 82;
  double print_matrix_data[] = {0,0,0,0,3,0};
  chem_matrix_view print_matrix_view;
  if (chem_matrix_view_init(&print_matrix_view,print_matrix_data,2,3,3,1) ||
      chem_matrix_view_set(&print_matrix_view,1,1,3.0)) return 83;
  FILE *dense_print_file = tmpfile();
  char dense_print_output[64];
  if (!dense_print_file ||
      chem_matrix_view_write(dense_print_file,&print_matrix_view) ||
      fflush(dense_print_file) || fseek(dense_print_file,0,SEEK_SET) ||
      fread(dense_print_output,1,sizeof(dense_print_output)-1,
            dense_print_file) == 0) return 84;
  const size_t dense_print_length = ftell(dense_print_file);
  if (dense_print_length >= sizeof(dense_print_output)) return 85;
  dense_print_output[dense_print_length] = '\0';
  fclose(dense_print_file);
  if (strcmp(dense_print_output,"0,0,0\n0,3,0\n")) return 86;
  const size_t lane_widths[] = {1, 2, 4, 8};
  for (size_t w = 0; w < sizeof(lane_widths)/sizeof(lane_widths[0]); ++w) {
    const size_t width = lane_widths[w], logical_rows = 3*width, columns = 4;
    double packed[96] = {0}, flat_packed[96] = {0};
    size_t storage_size = 0;
    chem_vector_matrix_view grouped;
    if (chem_vector_matrix_view_storage_size(logical_rows,columns,width,
                                               &storage_size) ||
        storage_size != logical_rows*columns ||
        chem_vector_matrix_view_init(&grouped,packed,logical_rows,columns,width))
      return 64;
    for (size_t row = 0; row < logical_rows; ++row)
      for (size_t column = 0; column < columns; ++column)
        if (chem_vector_matrix_view_set(&grouped,row,column,
                                        (double)(row*10+column))) return 65;
    if (chem_vector_matrix_view_copy_to_vector(&grouped,flat_packed,
                                                storage_size)) return 66;
    for (size_t row = 0; row < logical_rows; ++row)
      for (size_t column = 0; column < columns; ++column) {
        double got = -1.0;
        if (chem_vector_matrix_view_get(&grouped,row,column,&got) ||
            got != (double)(row*10+column) ||
            flat_packed[row*columns+column] != got) return 67;
      }
  }
  double packed_reduce_data[8] = {0};
  chem_vector_matrix_view packed_reduce;
  chem_const_vector_matrix_view packed_reduce_const;
  int any_large = 0, all_finite = 0;
  const double reduce_threshold = 5.0;
  if (chem_vector_matrix_view_init(&packed_reduce,packed_reduce_data,3,2,2))
    return 70;
  for (size_t row = 0; row < 3; ++row)
    for (size_t column = 0; column < 2; ++column)
      if (chem_vector_matrix_view_set(&packed_reduce,row,column,
                                      (double)(row*2+column+1))) return 71;
  if (chem_const_vector_matrix_view_init(&packed_reduce_const,
          packed_reduce_data,3,2,2) ||
      chem_vector_matrix_view_reduce_strict(&packed_reduce_const,
          matrix_sum_squares,NULL,0.0,&value) || value != 91.0 ||
      chem_vector_matrix_view_reduce_strict(&packed_reduce_const,
          matrix_count,NULL,0.0,&value) || value != 6.0 ||
      chem_vector_matrix_view_reduce_or(&packed_reduce_const,
          matrix_greater_than,(void *)&reduce_threshold,0,&any_large) || !any_large ||
      chem_vector_matrix_view_reduce_and(&packed_reduce_const,
          matrix_is_finite,NULL,1,&all_finite) || !all_finite)
    return 72;
  const int integer_matrix_values[] = {1,2,3,4,5,6};
  if (chem_vector_matrix_view_scale_from_int_vector(&packed_reduce,
          integer_matrix_values,6,1.5)) return 76;
  FILE *matrix_output_file = tmpfile();
  char matrix_output[128];
  const char *matrix_expected_output =
      "1.5,3\n4.5,6\n7.5,9\n";
  if (!matrix_output_file ||
      chem_vector_matrix_view_write(matrix_output_file,&packed_reduce) ||
      fflush(matrix_output_file) || fseek(matrix_output_file,0,SEEK_SET) ||
      fread(matrix_output,1,sizeof(matrix_output)-1,matrix_output_file) == 0)
    return 77;
  const size_t matrix_output_length = ftell(matrix_output_file);
  if (matrix_output_length >= sizeof(matrix_output)) return 78;
  matrix_output[matrix_output_length] = '\0';
  fclose(matrix_output_file);
  if (strcmp(matrix_output,matrix_expected_output) ||
      chem_vector_matrix_view_scale_from_int_vector(&packed_reduce,
          integer_matrix_values,5,1.5) == 0) return 79;
  chem_vector_matrix_view empty_packed_view;
  if (chem_vector_matrix_view_init(&empty_packed_view,NULL,0,0,2) ||
      chem_vector_matrix_view_for_each(&empty_packed_view,matrix_add_index,
                                       (void *)&scale)) return 80;
  double vector_a[] = {1,2,3}, vector_b[] = {10,20,30};
  double vector_out[3] = {0}, vector_mutable[] = {0,0,0};
  chem_const_matrix_view vector_a_view, vector_b_view;
  chem_const_matrix_view vector_short_view;
  chem_matrix_view vector_out_view, vector_mutable_view;
  double short_source[] = {1,2};
  const double vector_scale = 2.0;
  if (chem_const_matrix_view_init(&vector_a_view,vector_a,3,1,1,1) ||
      chem_const_matrix_view_init(&vector_b_view,vector_b,3,1,1,1) ||
      chem_const_matrix_view_init(&vector_short_view,short_source,2,1,1,1) ||
      chem_matrix_view_init(&vector_out_view,vector_out,3,1,1,1) ||
      chem_matrix_view_for_each_row_const_sources(&vector_out_view,
          &vector_a_view,&vector_b_view,matrix_add_vector_rows,
          (void *)&vector_scale) ||
      vector_out[0] != 21 || vector_out[1] != 42 || vector_out[2] != 63 ||
      chem_matrix_view_init(&vector_mutable_view,vector_mutable,3,1,1,1) ||
      chem_matrix_view_for_each_row_const_sources(&vector_mutable_view,
          &vector_a_view,&vector_b_view,matrix_add_vector_rows,NULL) ||
      vector_mutable[0] != 11 || vector_mutable[2] != 33 ||
      chem_matrix_view_for_each_row_const_sources(&vector_out_view,
          &vector_a_view,&vector_short_view,matrix_add_vector_rows,NULL) == 0)
    return 68;
  const double vector_a_again[] = {5,6,7};
  chem_const_matrix_view vector_a_again_view;
  if (chem_const_matrix_view_init(&vector_a_again_view,vector_a_again,3,1,1,1) ||
      chem_matrix_view_for_each_row_const_sources(&vector_out_view,
          &vector_a_again_view,&vector_b_view,matrix_add_vector_rows,NULL) ||
      vector_out[0] != 15 || vector_out[2] != 37) return 69;
  double grouped_a_data[8] = {0}, grouped_b_data[8] = {0},
         grouped_out_data[8] = {0};
  chem_vector_matrix_view grouped_a_mut, grouped_b_mut;
  chem_const_vector_matrix_view grouped_a_view, grouped_b_view;
  chem_vector_matrix_view grouped_out_view;
  if (chem_vector_matrix_view_init(&grouped_out_view,grouped_out_data,3,1,2) ||
      chem_vector_matrix_view_init(&grouped_a_mut,grouped_a_data,3,1,2) ||
      chem_vector_matrix_view_init(&grouped_b_mut,grouped_b_data,3,1,2))
    return 73;
  for (size_t row = 0; row < 3; ++row)
    if (chem_vector_matrix_view_set(&grouped_a_mut,row,0,(double)(row+1)) ||
        chem_vector_matrix_view_set(&grouped_b_mut,row,0,(double)(10*(row+1))))
      return 74;
  if (chem_const_vector_matrix_view_init(&grouped_a_view,grouped_a_data,3,1,2) ||
      chem_const_vector_matrix_view_init(&grouped_b_view,grouped_b_data,3,1,2) ||
      chem_vector_matrix_view_for_each_row_const_sources(&grouped_out_view,
          &grouped_a_view,&grouped_b_view,matrix_add_vector_rows,NULL) ||
      chem_vector_matrix_view_get(&grouped_out_view,0,0,&value) || value != 11.0 ||
      chem_vector_matrix_view_get(&grouped_out_view,2,0,&value) || value != 33.0)
    return 75;
  return 0;
}

static void csr_sum_sources(size_t block, double *destination,
                            const chem_csr_const_value_view *sources,
                            size_t source_count, void *context) {
  (void)context;
  double sum = 0.0;
  for (size_t i = 0; i < source_count; ++i) {
    double value = 0.0;
    (void)chem_csr_const_value_view_get(&sources[i],block,&value);
    sum += value;
  }
  *destination = sum;
}

static void csr_scale_source(size_t block, double *destination,
                             const chem_csr_const_value_view *sources,
                             size_t source_count, void *context) {
  const double scale = context ? *(const double *)context : 1.0;
  double value = 0.0;
  if (source_count)
    (void)chem_csr_const_value_view_get(&sources[0],block,&value);
  *destination = source_count ? scale*value : 0.0;
}

static void csr_count_blocks(size_t block, double *destination,
                             const chem_csr_const_value_view *sources,
                             size_t source_count, void *context) {
  (void)block; (void)destination; (void)sources; (void)source_count;
  ++*(size_t *)context;
}

static int check_vector_sparse_layout(size_t width, int compressed_by_column) {
  const size_t empty_offsets[] = {0,0,0,0,0};
  chem_csr_matrix empty;
  int is_zero = 0;
  size_t invalid_index = 0;
  if (chem_csr_matrix_init_vector(&empty,4,width+1,0,empty_offsets,NULL,NULL,
                                 width,compressed_by_column) ||
      chem_csr_matrix_fill(&empty,3.0) ||
      chem_csr_matrix_is_structural_zero(&empty,0,0,&is_zero) || !is_zero ||
      chem_csr_matrix_vector_index(&empty,0,0,0,&invalid_index) == 0)
    return 10;
  const size_t rows[] = {0,2,3,4,4};
  const size_t row_columns[] = {0,1,1,3};
  const size_t columns[] = {0,1,3,3,4};
  const size_t column_rows[] = {0,0,1,2};
  const size_t blocks = width + 1, nonzeros = 4;
  const size_t groups = (blocks + width - 1)/width;
  const size_t storage_size = groups*width*nonzeros;
  double values[64] = {0};
  chem_csr_matrix matrix;
  if (storage_size > sizeof(values)/sizeof(values[0]) ||
      chem_csr_matrix_init_vector(&matrix,4,blocks,nonzeros,
          compressed_by_column ? columns : rows,
          compressed_by_column ? column_rows : row_columns,values,width,
          compressed_by_column)) return 1;
  const size_t entry_rows[] = {0,0,1,2};
  const size_t entry_columns[] = {0,1,1,3};
  for (size_t block = 0; block < blocks; ++block)
    for (size_t entry = 0; entry < nonzeros; ++entry) {
      const size_t expected_index = (block/width)*nonzeros*width +
                                    entry*width + block%width;
      size_t actual_index = 0;
      const double value = (double)(block*10 + entry + 1);
      if (chem_csr_matrix_vector_index(&matrix,block,entry_rows[entry],
              entry_columns[entry],&actual_index) ||
          actual_index != expected_index ||
          chem_csr_matrix_set(&matrix,block,entry_rows[entry],
                              entry_columns[entry],value) ||
          values[expected_index] != value) return 2;
    }
  for (size_t block = 0; block < blocks; ++block)
    for (size_t entry = 0; entry < nonzeros; ++entry) {
      double value = 0.0;
      if (chem_csr_matrix_get(&matrix,block,entry_rows[entry],
                              entry_columns[entry],&value) ||
          value != (double)(block*10 + entry + 1)) return 3;
    }
  double diagonal_before[blocks];
  for (size_t block = 0; block < blocks; ++block)
    diagonal_before[block] = values[(block/width)*nonzeros*width + block%width];
  if (chem_csr_matrix_add_to_diagonal(&matrix,1.0)) return 4;
  for (size_t block = 0; block < blocks; ++block)
    if (values[(block/width)*nonzeros*width + block%width] !=
            diagonal_before[block] + 1.0) return 5;
  chem_csr_const_value_view view;
  double copied[10];
  if (chem_csr_matrix_get_const_value_view(&matrix,0,0,&view) ||
      chem_csr_value_view_copy_to_vector(&view,copied,blocks)) return 6;
  for (size_t block = 0; block < blocks; ++block)
    if (copied[block] != diagonal_before[block] + 1.0) return 7;
  chem_csr_value_view copy_target;
  if (chem_csr_matrix_get_value_view(&matrix,0,1,&copy_target) ||
      chem_csr_value_view_copy_from_vector(&copy_target,copied,blocks))
    return 13;
  for (size_t block = 0; block < blocks; ++block) {
    double value = 0.0;
    if (chem_csr_matrix_get(&matrix,block,0,1,&value) ||
        value != diagonal_before[block] + 1.0) return 14;
  }
  int integers[10];
  for (size_t block = 0; block < blocks; ++block)
    integers[block] = (int)(5*(block+1));
  if (chem_csr_value_view_scale_from_int_vector(&copy_target,integers,blocks,
                                                 1.5)) return 16;
  for (size_t block = 0; block < blocks; ++block) {
    double value = 0.0;
    if (chem_csr_matrix_get(&matrix,block,0,1,&value) ||
        value != 7.5*(double)(block+1)) return 17;
  }
  chem_csr_value_view callback_destination;
  chem_csr_const_value_view callback_source;
  if (chem_csr_matrix_get_value_view(&matrix,2,3,&callback_destination) ||
      chem_csr_matrix_get_const_value_view(&matrix,0,1,&callback_source) ||
      chem_csr_matrix_apply_blocks(&callback_destination,&callback_source,1,
                                   csr_scale_source,(void *)&(const double){2.0}))
    return 11;
  for (size_t block = 0; block < blocks; ++block) {
    double result = 0.0;
    if (chem_csr_matrix_get(&matrix,block,2,3,&result) ||
        result != 15.0*(double)(block+1)) return 12;
  }
  const size_t other_width = width == 4 ? 2 : 4;
  const size_t other_groups = (blocks + other_width - 1)/other_width;
  double other_values[64] = {0};
  chem_csr_matrix other_matrix;
  if (chem_csr_matrix_init_vector(&other_matrix,4,blocks,nonzeros,
          compressed_by_column ? columns : rows,
          compressed_by_column ? column_rows : row_columns,other_values,
          other_width,compressed_by_column)) return 18;
  for (size_t block = 0; block < blocks; ++block)
    if (chem_csr_matrix_set(&other_matrix,block,0,1,(double)(block+1)))
      return 19;
  chem_csr_value_view cross_destination;
  chem_csr_const_value_view cross_source;
  double dense_group_values[10];
  for (size_t block = 0; block < blocks; ++block)
    dense_group_values[block] = 100.0 + (double)block;
  chem_csr_const_value_view dense_group_source = {
      dense_group_values,blocks,1,1};
  if (other_groups*other_width*nonzeros > sizeof(other_values)/sizeof(double) ||
      chem_csr_matrix_get_value_view(&matrix,1,1,&cross_destination) ||
      chem_csr_matrix_get_const_value_view(&other_matrix,0,1,&cross_source))
    return 20;
  const chem_csr_const_value_view cross_sources[] = {
      cross_source,dense_group_source};
  if (chem_csr_matrix_apply_blocks(&cross_destination,cross_sources,2,
                                   csr_sum_sources,NULL)) return 22;
  for (size_t block = 0; block < blocks; ++block) {
    double result = 0.0;
    if (chem_csr_matrix_get(&matrix,block,1,1,&result) ||
        result != 101.0 + 2.0*(double)block) return 21;
  }
  if (chem_csr_matrix_fill(&matrix,2.0)) return 8;
  for (size_t i = 0; i < storage_size; ++i)
    if (values[i] != 2.0) return 9;
  FILE *stream = tmpfile();
  if (!stream || chem_csr_matrix_write(stream,&matrix,0) ||
      chem_csr_matrix_write(stream,&matrix,1) || ftell(stream) <= 0) {
    if (stream) fclose(stream);
    return 15;
  }
  fclose(stream);
  return 0;
}

static int check_csr_matrix_api(void) {
  const size_t empty_offsets[] = {0,0,0,0};
  chem_csr_matrix empty;
  if (chem_csr_matrix_init(&empty,3,1,0,empty_offsets,NULL,NULL) ||
      chem_csr_matrix_fill(&empty,2.0) ||
      chem_csr_matrix_is_structural_zero(&empty,2,1,(int[1]){0})) return 1;
  int is_zero = 0;
  if (chem_csr_matrix_is_structural_zero(&empty,2,1,&is_zero) || !is_zero ||
      chem_csr_matrix_vector_index(&empty,0,0,0,(size_t[1]){0}) == 0)
    return 2;
  const size_t offsets[] = {0,1,2,4,5};
  const size_t columns[] = {1,1,1,3,2};
  double values[5*53];
  for (size_t i = 0; i < 5*53; ++i) values[i] = 24.0;
  chem_csr_matrix matrix;
  if (chem_csr_matrix_init(&matrix,4,53,5,offsets,columns,values)) return 3;
  size_t index = 0;
  if (chem_csr_matrix_vector_index(&matrix,0,3,2,&index) || index != 4 ||
      chem_csr_matrix_vector_index(&matrix,0,2,3,&index) || index != 3 ||
      chem_csr_matrix_vector_index(&matrix,2,2,1,&index) || index != 12 ||
      chem_csr_matrix_vector_index(&matrix,0,2,2,&index) == 0) return 4;
  if (chem_csr_matrix_is_structural_zero(&matrix,0,1,&is_zero) || is_zero ||
      chem_csr_matrix_is_structural_zero(&matrix,2,2,&is_zero) || !is_zero ||
      chem_csr_matrix_set(&matrix,0,2,1,45.0) || values[2] != 45.0)
    return 5;
  double value = 0.0;
  if (chem_csr_matrix_get(&matrix,0,2,1,&value) || value != 45.0 ||
      chem_csr_matrix_get(&matrix,0,2,2,&value) == 0) return 6;
  if (chem_csr_matrix_set(&matrix,0,1,1,1.0) ||
      chem_csr_matrix_set(&matrix,1,1,1,3.0) ||
      chem_csr_matrix_set(&matrix,2,1,1,5.0) ||
      chem_csr_matrix_add_to_diagonal(&matrix,7.0) ||
      values[1] != 8.0 || values[5+1] != 10.0 || values[2*5+1] != 12.0 ||
      values[0] != 24.0) return 7;
  if (chem_csr_matrix_fill(&matrix,2.0)) return 8;
  for (size_t i = 0; i < 5*53; ++i) if (values[i] != 2.0) return 9;
  const size_t invalid_offsets[] = {0,2,1,1,1};
  const size_t invalid_column[] = {3};
  const size_t invalid_column_offsets[] = {0,1,1,1,1};
  if (chem_csr_matrix_init(&matrix,4,1,1,invalid_offsets,columns,values) == 0 ||
      chem_csr_matrix_init(&matrix,3,1,1,invalid_column_offsets,
                           invalid_column,values) == 0)
    return 10;
  const size_t array_offsets[] = {0,1,2,4,5};
  const size_t array_columns[] = {0,1,2,3,3};
  double array_values[15];
  for (size_t i = 0; i < 15; ++i) array_values[i] = 1.0;
  if (chem_csr_matrix_init(&matrix,4,3,5,array_offsets,array_columns,
                           array_values)) return 11;
  if (chem_csr_matrix_set(&matrix,0,0,0,1.0) ||
      chem_csr_matrix_set(&matrix,1,1,1,2.0) ||
      chem_csr_matrix_set(&matrix,2,2,2,3.0) ||
      chem_csr_matrix_set(&matrix,2,2,3,4.0) ||
      chem_csr_matrix_set(&matrix,2,3,3,5.0)) return 12;
  chem_csr_const_value_view sources[4];
  chem_csr_value_view temporary_view, output_view;
  double temporary[3] = {0};
  if (chem_csr_matrix_get_const_value_view(&matrix,0,0,&sources[0]) ||
      chem_csr_matrix_get_const_value_view(&matrix,1,1,&sources[1]) ||
      chem_csr_matrix_get_const_value_view(&matrix,2,2,&sources[2]) ||
      chem_csr_matrix_get_const_value_view(&matrix,2,3,&sources[3])) return 13;
  temporary_view = (chem_csr_value_view){temporary,3,1,1};
  if (chem_csr_matrix_apply_blocks(&temporary_view,sources,4,
                                   csr_sum_sources,NULL)) return 14;
  if (chem_csr_matrix_get_value_view(&matrix,2,3,&output_view) ||
      chem_csr_matrix_apply_blocks(&output_view,
          (chem_csr_const_value_view[1]){{temporary,3,1,1}},1,
          csr_scale_source,(void *)&(const double){2.0}) ||
      array_values[3] != 8.0 || array_values[8] != 10.0 ||
      array_values[13] != 18.0) return 15;
  chem_csr_const_value_view wrong_source = {temporary,2,1,1};
  if (chem_csr_matrix_apply_blocks(&output_view,&wrong_source,1,
                                   csr_scale_source,NULL) == 0) return 16;
  const size_t a_offsets[] = {0,1,2,2};
  const size_t a_columns[] = {1,2};
  const size_t b_offsets[] = {0,1,2,3};
  const size_t b_columns[] = {0,1,2};
  double a_values[6] = {10,20,11,22,12,24};
  double b_values[9] = {0,0,0,4,3,5,8,6,10};
  chem_csr_matrix matrix_a, matrix_b;
  chem_csr_value_view a_destination;
  chem_csr_const_value_view a_source, b_source;
  if (chem_csr_matrix_init(&matrix_a,3,3,2,a_offsets,a_columns,a_values) ||
      chem_csr_matrix_init(&matrix_b,3,3,3,b_offsets,b_columns,b_values) ||
      chem_csr_matrix_get_value_view(&matrix_a,1,2,&a_destination) ||
      chem_csr_matrix_get_const_value_view(&matrix_a,0,1,&a_source) ||
      chem_csr_matrix_get_const_value_view(&matrix_b,2,2,&b_source) ||
      chem_csr_matrix_apply_blocks(&a_destination,
          (chem_csr_const_value_view[2]){a_source,b_source},2,
          csr_sum_sources,NULL) ||
      a_values[1] != 10.0 || a_values[3] != 16.0 || a_values[5] != 22.0)
    return 17;
  chem_csr_value_view block_destination;
  chem_csr_const_value_view block_source;
  double block_copy[3];
  if (chem_csr_matrix_get_value_view(&matrix_a,1,2,&block_destination) ||
      chem_csr_matrix_get_const_value_view(&matrix_a,0,1,&block_source) ||
      chem_csr_value_view_copy(&block_destination,&block_source) ||
      a_values[1] != 10.0 || a_values[3] != 11.0 || a_values[5] != 12.0 ||
      chem_csr_value_view_copy_to_vector(&block_source,block_copy,3) ||
      block_copy[0] != 10.0 || block_copy[1] != 11.0 || block_copy[2] != 12.0 ||
      chem_csr_value_view_fill(&block_destination,3.2) ||
      a_values[1] != 3.2 || a_values[3] != 3.2 || a_values[5] != 3.2 ||
      chem_csr_value_view_copy_from_vector(&block_destination,block_copy,3) ||
      a_values[1] != 10.0 || a_values[3] != 11.0 || a_values[5] != 12.0 ||
      chem_csr_value_view_copy_from_vector(&block_destination,block_copy,2) == 0)
    return 18;
  if (chem_csr_value_view_copy_to_vector(&block_source,block_copy,2) == 0 ||
      chem_csr_value_view_copy_to_vector(&block_source,NULL,0) == 0)
    return 19;
  const size_t zero_block_offsets[] = {0,1,1,1};
  const size_t zero_block_columns[] = {1};
  chem_csr_matrix zero_block_matrix;
  chem_csr_value_view zero_destination;
  chem_csr_const_value_view zero_source;
  size_t callback_count = 0;
  if (chem_csr_matrix_init(&zero_block_matrix,3,0,1,zero_block_offsets,
                           zero_block_columns,NULL) ||
      chem_csr_matrix_get_value_view(&zero_block_matrix,0,1,
                                     &zero_destination) ||
      chem_csr_matrix_get_const_value_view(&zero_block_matrix,0,1,
                                           &zero_source) ||
      chem_csr_matrix_apply_blocks(&zero_destination,&zero_source,1,
                                   csr_count_blocks,&callback_count) ||
      callback_count != 0 ||
      chem_csr_value_view_copy_to_vector(&zero_source,NULL,0) ||
      chem_csr_value_view_copy_from_vector(&zero_destination,NULL,0)) return 20;
  const size_t csc_offsets[] = {0,0,3,4,5};
  const size_t csc_rows[] = {0,1,2,3,2};
  double csc_values[10];
  for (size_t i = 0; i < 10; ++i) csc_values[i] = 24.0;
  chem_csr_matrix csc;
  if (chem_csr_matrix_init_csc(&csc,4,2,5,csc_offsets,csc_rows,csc_values) ||
      chem_csr_matrix_vector_index(&csc,0,3,2,&index) || index != 3 ||
      chem_csr_matrix_vector_index(&csc,1,2,3,&index) || index != 9 ||
      chem_csr_matrix_vector_index(&csc,0,2,2,&index) == 0 ||
      chem_csr_matrix_set(&csc,1,2,3,31.0) || csc_values[9] != 31.0 ||
      chem_csr_matrix_add_to_diagonal(&csc,1.0) ||
      csc_values[1] != 25.0 || csc_values[6] != 25.0 ||
      chem_csr_matrix_is_structural_zero(&csc,0,0,&is_zero) || !is_zero)
    return 20;
  chem_csr_value_view csc_destination;
  chem_csr_const_value_view csc_sources[2];
  double csc_copy[2];
  if (chem_csr_matrix_get_value_view(&csc,2,1,&csc_destination) ||
      chem_csr_matrix_get_const_value_view(&csc,0,1,&csc_sources[0]) ||
      chem_csr_matrix_get_const_value_view(&csc,1,1,&csc_sources[1]) ||
      chem_csr_matrix_apply_blocks(&csc_destination,csc_sources,2,
                                   csr_sum_sources,NULL) ||
      csc_values[2] != 49.0 || csc_values[7] != 49.0 ||
      chem_csr_value_view_copy_to_vector(&csc_sources[0],csc_copy,2) ||
      csc_copy[0] != 24.0 || csc_copy[1] != 24.0 ||
      chem_csr_value_view_fill(&csc_destination,3.0) ||
      csc_values[2] != 3.0 || csc_values[7] != 3.0 ||
      chem_csr_value_view_copy_from_vector(&csc_destination,csc_copy,2) ||
      csc_values[2] != 24.0 || csc_values[7] != 24.0 ||
      chem_csr_value_view_copy_to_vector(&csc_sources[0],csc_copy,1) == 0)
    return 21;
  double print_values[15];
  for (size_t i = 0; i < 15; ++i) print_values[i] = 24.0;
  chem_csr_matrix print_matrix;
  if (chem_csr_matrix_init(&print_matrix,4,3,5,offsets,columns,print_values) ||
      chem_csr_matrix_set(&print_matrix,0,1,1,1.0) ||
      chem_csr_matrix_set(&print_matrix,1,1,1,3.0) ||
      chem_csr_matrix_set(&print_matrix,2,1,1,5.0) ||
      chem_csr_matrix_add_to_diagonal(&print_matrix,7.0)) return 22;
  const char *dense_expected =
      "Block 0\n0,24,0,0\n0,8,0,0\n0,24,0,24\n0,0,24,0\n"
      "Block 1\n0,24,0,0\n0,10,0,0\n0,24,0,24\n0,0,24,0\n"
      "Block 2\n0,24,0,0\n0,12,0,0\n0,24,0,24\n0,0,24,0\n";
  FILE *stream = tmpfile();
  char output[1024];
  if (!stream || chem_csr_matrix_write(stream,&print_matrix,0) ||
      fflush(stream) || fseek(stream,0,SEEK_SET) ||
      fread(output,1,sizeof(output)-1,stream) == 0) return 23;
  const size_t dense_length = ftell(stream);
  if (dense_length >= sizeof(output)) return 24;
  output[dense_length] = '\0';
  fclose(stream);
  if (strcmp(output,dense_expected)) return 25;
  const char *nonzero_expected =
      "Block 0\n0, 1, 24\n1, 1, 8\n2, 1, 24\n2, 3, 24\n3, 2, 24\n"
      "Block 1\n0, 1, 24\n1, 1, 10\n2, 1, 24\n2, 3, 24\n3, 2, 24\n"
      "Block 2\n0, 1, 24\n1, 1, 12\n2, 1, 24\n2, 3, 24\n3, 2, 24\n";
  stream = tmpfile();
  if (!stream || chem_csr_matrix_write(stream,&print_matrix,1) ||
      fflush(stream) || fseek(stream,0,SEEK_SET) ||
      fread(output,1,sizeof(output)-1,stream) == 0) return 26;
  const size_t nonzero_length = ftell(stream);
  if (nonzero_length >= sizeof(output)) return 27;
  output[nonzero_length] = '\0';
  fclose(stream);
  if (strcmp(output,nonzero_expected)) return 28;
  const size_t print_csc_offsets[] = {0,0,3,4,5};
  const size_t print_csc_rows[] = {0,1,2,3,2};
  double print_csc_values[15];
  for (size_t i = 0; i < 15; ++i) print_csc_values[i] = 24.0;
  chem_csr_matrix print_csc;
  if (chem_csr_matrix_init_csc(&print_csc,4,3,5,print_csc_offsets,
                                print_csc_rows,print_csc_values) ||
      chem_csr_matrix_set(&print_csc,0,1,1,8.0) ||
      chem_csr_matrix_set(&print_csc,1,1,1,10.0) ||
      chem_csr_matrix_set(&print_csc,2,1,1,12.0)) return 29;
  stream = tmpfile();
  if (!stream || chem_csr_matrix_write(stream,&print_csc,0) || fflush(stream) ||
      fseek(stream,0,SEEK_SET) || fread(output,1,sizeof(output)-1,stream) == 0)
    return 30;
  const size_t csc_length = ftell(stream);
  if (csc_length >= sizeof(output)) return 31;
  output[csc_length] = '\0';
  fclose(stream);
  if (strcmp(output,dense_expected)) return 32;
  const size_t bad_csc_offsets[] = {0,0,0,0};
  const size_t bad_csc_rows[] = {3};
  if (chem_csr_matrix_init_csc(&csc,3,1,1,bad_csc_offsets,bad_csc_rows,
                               print_values) == 0) return 33;
  const size_t widths[] = {1,2,4,8};
  for (size_t i = 0; i < sizeof(widths)/sizeof(widths[0]); ++i) {
    const int row_status = check_vector_sparse_layout(widths[i],0);
    const int column_status = check_vector_sparse_layout(widths[i],1);
    if (row_status) return 34 + (int)i*30 + row_status;
    if (column_status) return 34 + (int)i*30 + 15 + column_status;
  }
  return 0;
}

static int lu_reconstruction_matches(const double *matrix, const double *lower,
                                     const double *upper, size_t n) {
  for (size_t row = 0; row < n; ++row)
    for (size_t column = 0; column < n; ++column) {
      if (row < column && lower[row*n+column] != 0.0) return 0;
      if (column < row && upper[row*n+column] != 0.0) return 0;
      double value = 0.0;
      for (size_t k = 0; k < n; ++k)
        value += lower[row*n+k]*upper[k*n+column];
      const double error = fabs(value-matrix[row*n+column]);
      const double tolerance = 1e-9*(1.0+fabs(matrix[row*n+column]));
      if (error > tolerance) return 0;
    }
  return 1;
}

static int check_lu_input(size_t n, size_t blocks, int kind) {
  double original[5*100], lower[100], upper[100], work[100];
  const size_t matrix_size = n*n;
  for (size_t block = 0; block < blocks; ++block) {
    double *matrix = original + block*matrix_size;
    for (size_t row = 0; row < n; ++row)
      for (size_t column = 0; column < n; ++column) {
        double value = 0.0;
        if (kind == 0 && n == 3) {
          const double example[9] = {2,-1,-2,-4,6,3,-4,-2,8};
          value = example[row*n+column];
        } else if (kind == 2) {
          value = row == column ? 2.0+(double)row : 0.0;
        } else {
          const int pattern = (int)((row*13+column*7+block*3)%17)-8;
          value = row == column ? 20.0+(double)(row+block) :
                  0.1*(double)pattern;
        }
        matrix[row*n+column] = value;
      }
    for (size_t i = 0; i < matrix_size; ++i) lower[i] = upper[i] = INFINITY;
    if (chem_lu_decompose_doolittle(n,matrix,lower,upper) ||
        !lu_reconstruction_matches(matrix,lower,upper,n)) return 1;
    for (size_t i = 0; i < matrix_size; ++i) lower[i] = upper[i] = NAN;
    if (chem_lu_decompose_mozart(n,matrix,lower,upper) ||
        !lu_reconstruction_matches(matrix,lower,upper,n)) return 2;
    for (size_t i = 0; i < matrix_size; ++i) work[i] = matrix[i];
    if (chem_lu_decompose_doolittle_in_place(n,work)) return 3;
    for (size_t row = 0; row < n; ++row)
      for (size_t column = 0; column < n; ++column) {
        double value = 0.0;
        for (size_t k = 0; k <= column; ++k) {
          if (row == k) value += work[k*n+column];
          else if (k < row) value += work[row*n+k]*work[k*n+column];
        }
        if (fabs(value-matrix[row*n+column]) >
            1e-9*(1.0+fabs(matrix[row*n+column]))) return 4;
      }
    for (size_t i = 0; i < matrix_size; ++i) work[i] = matrix[i];
    if (chem_lu_decompose_mozart_in_place(n,work)) return 5;
    for (size_t row = 0; row < n; ++row)
      for (size_t column = 0; column < n; ++column) {
        double value = 0.0;
        for (size_t k = 0; k <= column; ++k) {
          if (row == k) value += work[k*n+column];
          else if (k < row) value += work[row*n+k]*work[k*n+column];
        }
        if (fabs(value-matrix[row*n+column]) >
            1e-9*(1.0+fabs(matrix[row*n+column]))) return 6;
      }
  }
  return 0;
}

static int check_lu_decomposition(void) {
  if (check_lu_input(3,1,0) || check_lu_input(10,1,1) ||
      check_lu_input(10,5,1) || check_lu_input(5,5,2)) return 1;
  double singular[] = {1.0,2.0,2.0,4.0}, lower[4], upper[4];
  if (chem_lu_decompose_doolittle(2,singular,lower,upper) == 0 ||
      chem_lu_decompose_mozart_in_place(2,singular) == 0 ||
      chem_lu_decompose_doolittle(2,NULL,lower,upper) == 0) return 2;
  return 0;
}

static int check_linear_solve_case(size_t n, size_t blocks, int diagonal) {
  double matrix[5*100], rhs[5*10], solution[10], matrix_work[100], rhs_work[10];
  for (size_t block = 0; block < blocks; ++block) {
    double *a = matrix + block*n*n;
    double *b = rhs + block*n;
    for (size_t row = 0; row < n; ++row) {
      b[row] = (double)(row+1) - 0.25*(double)block;
      for (size_t column = 0; column < n; ++column) {
        if (diagonal) a[row*n+column] = row == column ? 2.0+(double)row : 0.0;
        else {
          const int pattern = (int)((row*11+column*5+block*7)%19)-9;
          a[row*n+column] = row == column ? 25.0+(double)(row+block) :
                            0.1*(double)pattern;
        }
      }
    }
    for (size_t i = 0; i < n; ++i) solution[i] = NAN;
    if (chem_linear_solve(n,a,b,solution)) return 1;
    for (size_t row = 0; row < n; ++row) {
      double value = 0.0;
      for (size_t column = 0; column < n; ++column)
        value += a[row*n+column]*solution[column];
      if (fabs(value-b[row]) > 1e-9*(1.0+fabs(b[row]))) return 2;
    }
    for (size_t i = 0; i < n*n; ++i) matrix_work[i] = a[i];
    for (size_t i = 0; i < n; ++i) rhs_work[i] = b[i];
    if (chem_linear_solve_in_place(n,matrix_work,rhs_work)) return 3;
    for (size_t row = 0; row < n; ++row) {
      double value = 0.0;
      for (size_t column = 0; column < n; ++column)
        value += a[row*n+column]*rhs_work[column];
      if (fabs(value-b[row]) > 1e-9*(1.0+fabs(b[row]))) return 4;
    }
  }
  return 0;
}

static int check_chapman_sparse_solve(void) {
  double lower[9*9] = {0}, upper[9*9] = {0}, matrix[9*9] = {0};
  const double expected[] = {0.5,0.5,0.5,0.5,0.5,-0.5,-1.0,0.5,0.0};
  double rhs[9];
  for (size_t i = 0; i < 9; ++i) {
    lower[i*9+i] = 1.0;
    upper[i*9+i] = 1.0;
  }
  const size_t lower_positions[][2] = {
      {6,0}, {7,0}, {8,0}, {5,4}, {6,4}, {6,5}, {7,6}, {8,6}, {8,7}};
  const size_t upper_positions[][2] = {
      {5,7}, {5,8}, {6,7}, {6,8}, {7,8}};
  for (size_t i = 0; i < 9; ++i)
    lower[lower_positions[i][0]*9+lower_positions[i][1]] = 1.0;
  for (size_t i = 0; i < 5; ++i)
    upper[upper_positions[i][0]*9+upper_positions[i][1]] = 1.0;
  for (size_t row = 0; row < 9; ++row)
    for (size_t column = 0; column < 9; ++column)
      for (size_t k = 0; k < 9; ++k)
        matrix[row*9+column] += lower[row*9+k]*upper[k*9+column];
  for (size_t i = 0; i < 9; ++i) rhs[i] = 0.5;
  if (chem_linear_solve_in_place(9,matrix,rhs)) return 1;
  for (size_t i = 0; i < 9; ++i)
    if (fabs(rhs[i]-expected[i]) > 1e-14) return 2;
  return 0;
}

static int check_linear_solvers(void) {
  const double dense[] = {2,-1,-2,-4,6,3,-4,-2,8};
  const double rhs[] = {23,42,9};
  double solution[3];
  if (chem_linear_solve(3,dense,rhs,solution)) return 1;
  for (size_t i = 0; i < 3; ++i) {
    double result = 0.0;
    for (size_t j = 0; j < 3; ++j) result += dense[i*3+j]*solution[j];
    if (fabs(result-rhs[i]) > 1e-10) return 2;
  }
  if (check_linear_solve_case(10,1,0) ||
      check_linear_solve_case(10,5,0) ||
      check_linear_solve_case(5,5,1) || check_chapman_sparse_solve()) return 3;
  double singular[] = {1,2,2,4}, singular_rhs[] = {3,6};
  if (chem_linear_solve(2,singular,singular_rhs,solution) == 0 ||
      chem_linear_solve_in_place(2,singular,singular_rhs) == 0 ||
      chem_linear_solve(2,NULL,singular_rhs,solution) == 0) return 4;
  const size_t order = 50;
  unsigned char pattern[order*order], reordered[order*order];
  size_t permutation[order];
  for (size_t row = 0; row < order; ++row)
    for (size_t column = 0; column < order; ++column)
      pattern[row*order+column] = (unsigned char)(row == column ||
          ((row*17+column*31+row*column)%13 < 3));
  if (chem_markowitz_reorder(order,pattern,permutation)) return 5;
  unsigned char seen[50] = {0};
  size_t sum = 0;
  for (size_t i = 0; i < order; ++i) {
    if (permutation[i] >= order || seen[permutation[i]]) return 6;
    seen[permutation[i]] = 1;
    sum += permutation[i];
  }
  if (sum != order*(order-1)/2) return 7;
  for (size_t row = 0; row < order; ++row)
    for (size_t column = 0; column < order; ++column)
      reordered[row*order+column] =
          pattern[permutation[row]*order+permutation[column]];
  unsigned char original_factor[order*order], reordered_factor[order*order];
  memcpy(original_factor,pattern,sizeof(pattern));
  memcpy(reordered_factor,reordered,sizeof(reordered));
  for (size_t pivot = 0; pivot < order; ++pivot)
    for (size_t row = pivot+1; row < order; ++row)
      if (original_factor[row*order+pivot])
        for (size_t column = pivot+1; column < order; ++column)
          if (original_factor[pivot*order+column])
            original_factor[row*order+column] = 1;
  for (size_t pivot = 0; pivot < order; ++pivot)
    for (size_t row = pivot+1; row < order; ++row)
      if (reordered_factor[row*order+pivot])
        for (size_t column = pivot+1; column < order; ++column)
          if (reordered_factor[pivot*order+column])
            reordered_factor[row*order+column] = 1;
  size_t original_fill = 0, reordered_fill = 0;
  for (size_t row = 0; row < order; ++row)
    for (size_t column = 0; column < order; ++column) {
      original_fill += original_factor[row*order+column] != 0;
      reordered_fill += reordered_factor[row*order+column] != 0;
    }
  if (original_fill <= reordered_fill ||
      chem_markowitz_reorder(order,NULL,permutation) == 0) return 8;
  return 0;
}

static int check_rosenbrock_species_tolerances(void) {
  const chem_property foo_properties[] = {{
      .name="absolute tolerance", .kind=CHEM_PROPERTY_DOUBLE,
      .value.real_value=1e-7}};
  const chem_property bar_properties[] = {{
      .name="absolute tolerance", .kind=CHEM_PROPERTY_DOUBLE,
      .value.real_value=1e-8}};
  const chem_species species[] = {
      {.name="foo",.properties=foo_properties,.property_count=1},
      {.name="bar",.properties=bar_properties,.property_count=1}};
  const chem_term reactant[] = {{.species=0,.coefficient=1.0}};
  const chem_term product[] = {{.species=1,.coefficient=1.0}};
  const chem_reaction reaction = {.kind=CHEM_ARRHENIUS,
      .reactants=reactant,.reactant_count=1,.products=product,
      .product_count=1,.p={2e-11,0,110,300,0}};
  const chem_mechanism mechanism = {.species=species,.species_count=2,
      .reactions=&reaction,.reaction_count=1};
  double tolerances[2];
  if (chem_mechanism_absolute_tolerances(&mechanism,tolerances,2) ||
      tolerances[0]!=1e-7 || tolerances[1]!=1e-8 ||
      chem_mechanism_absolute_tolerances(&mechanism,tolerances,1)==0)
    return 1;
  chem_solver_options options = chem_solver_default_options();
  options.absolute_tolerances=tolerances;
  options.absolute_tolerance_count=2;
  double state[] = {1.0,0.0};
  chem_solver_workspace workspace = {0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 2;
  const int status = chem_rosenbrock3_solve_with_workspace(
      &mechanism,&workspace,state,298.15,101325.0,40.0,NULL,1.0,
      &options,NULL);
  chem_solver_workspace_destroy(&workspace);
  if (status || !(state[0]<1.0) || !(state[1]>0.0)) return 3;
  return 0;
}

static int check_rejected_step_determinism(void) {
  const chem_species species[] = {{.name="foo"},{.name="bar"},
      {.name="baz"},{.name="quz"},{.name="quuz"}};
  const chem_term foo_baz[] = {{.species=0,.coefficient=1},
                                {.species=2,.coefficient=1}};
  const chem_term bar_quuz[] = {{.species=1,.coefficient=1},
                                 {.species=4,.coefficient=2.4}};
  const chem_term bar[] = {{.species=1,.coefficient=1}};
  const chem_term foo_quz[] = {{.species=0,.coefficient=1},
                                {.species=3,.coefficient=1.4}};
  const chem_term baz[] = {{.species=2,.coefficient=1}};
  const chem_reaction reactions[] = {
      {.kind=CHEM_ARRHENIUS,.reactants=foo_baz,.reactant_count=2,
       .products=bar_quuz,.product_count=2,.p={2e-11,0,110,300,0}},
      {.kind=CHEM_ARRHENIUS,.reactants=bar,.reactant_count=1,
       .products=foo_quz,.product_count=2,.p={1e-6,0,0,300,0}},
      {.kind=CHEM_ARRHENIUS,.reactants=baz,.reactant_count=1,
       .p={3.5e-6,0,0,300,0}}};
  const chem_mechanism mechanism = {.species=species,.species_count=5,
      .reactions=reactions,.reaction_count=3};
  chem_solver_options options = chem_solver_default_options();
  options.initial_step=1000.0;
  double first[5]={1e12,1e12,1e12,1e12,1e12};
  double second[5];
  memcpy(second,first,sizeof(first));
  chem_solver_workspace a={0},b={0};
  if (chem_solver_workspace_init(&a,&mechanism) ||
      chem_solver_workspace_init(&b,&mechanism)) {
    chem_solver_workspace_destroy(&a);
    chem_solver_workspace_destroy(&b);
    return 1;
  }
  chem_solver_stats stats_a={0},stats_b={0};
  const int status_a=chem_rosenbrock3_solve_with_workspace(
      &mechanism,&a,first,298.15,101325.0,40.0,NULL,1000.0,&options,&stats_a);
  const int status_b=chem_rosenbrock3_solve_with_workspace(
      &mechanism,&b,second,298.15,101325.0,40.0,NULL,1000.0,&options,&stats_b);
  chem_solver_workspace_destroy(&a);
  chem_solver_workspace_destroy(&b);
  if (status_a || status_b) return 10;
  if (stats_a.rejected<2) return 20;
  if (stats_a.internal_steps!=stats_b.internal_steps ||
      stats_a.accepted!=stats_b.accepted || stats_a.rejected!=stats_b.rejected)
    return 30;
  for (size_t i=0;i<5;++i)
    if (first[i]!=second[i]) return 40;
  return 0;
}

static int check_stiff_photolysis_chain(int backward_euler) {
  const chem_species species[] = {{.name="A1"},{.name="A2"},
      {.name="B"},{.name="C"}};
  const chem_term a1[] = {{.species=0,.coefficient=1}};
  const chem_term a2[] = {{.species=1,.coefficient=1}};
  const chem_term b[] = {{.species=2,.coefficient=1}};
  const chem_term products_b[] = {{.species=2,.coefficient=1}};
  const chem_term products_c[] = {{.species=3,.coefficient=1}};
  const chem_term products_a2[] = {{.species=1,.coefficient=1}};
  const chem_term products_a1[] = {{.species=0,.coefficient=1}};
  const chem_reaction reactions[] = {
      {.kind=CHEM_USER_DEFINED,.label="photoA1B",.reactants=a1,
       .reactant_count=1,.products=products_b,.product_count=1,
       .p={[0]=1,[7]=0}},
      {.kind=CHEM_USER_DEFINED,.label="photoA2B",.reactants=a2,
       .reactant_count=1,.products=products_b,.product_count=1,
       .p={[0]=1,[7]=1}},
      {.kind=CHEM_USER_DEFINED,.label="photoB",.reactants=b,
       .reactant_count=1,.products=products_c,.product_count=1,
       .p={[0]=1,[7]=2}},
      {.kind=CHEM_ARRHENIUS,.reactants=a1,.reactant_count=1,
       .products=products_a2,.product_count=1,.p={4e10,0,0,300,0}},
      {.kind=CHEM_ARRHENIUS,.reactants=a2,.reactant_count=1,
       .products=products_a1,.product_count=1,.p={3.6e10,0,0,300,0}}};
  const char *labels[] = {"photoA1B","photoA2B","photoB"};
  const chem_mechanism mechanism = {.species=species,.species_count=4,
      .reactions=reactions,.reaction_count=5,.custom_labels=labels,
      .custom_parameter_count=3};
  const double custom[] = {2e-3,2e-3,3e-3};
  double check_rates[5];
  const double temperature=272.5, pressure=101253.3;
  const double density=pressure/(8.31446261815324*temperature);
  if (chem_rate_constants(&mechanism,temperature,pressure,density,custom,check_rates)) return 3;
  for (size_t i=0;i<3;++i)
    if (fabs(check_rates[i]-custom[i])>1e-15) return 4;
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  chem_backward_euler_options be_options=chem_backward_euler_default_options();
  be_options.initial_step=1.0;
  be_options.relative_tolerance=1e-6;
  be_options.absolute_tolerance=1e-3;
  chem_solver_options ros_options=chem_solver_default_options();
  ros_options.relative_tolerance=1e-6;
  ros_options.absolute_tolerance=1e-3;
  ros_options.initial_step=1.0;
  double states[3][4];
  const double initial[] = {1.0,2.0/3.0,1.0/3.0};
  for (size_t cell=0;cell<3;++cell) {
    states[cell][0]=0.5*initial[cell]; states[cell][1]=0.5*initial[cell];
    states[cell][2]=0.0; states[cell][3]=0.0;
  }
  int status=0;
  for (size_t step=1;step<1000 && !status;++step) {
    for (size_t cell=0;cell<3 && !status;++cell) {
      status=backward_euler
          ? chem_backward_euler_solve_with_workspace(&mechanism,&workspace,
                states[cell],temperature,pressure,density,custom,1.0,
                &be_options,NULL)
          : chem_rosenbrock3_solve_with_workspace(&mechanism,&workspace,
                states[cell],temperature,pressure,density,custom,1.0,
                &ros_options,NULL);
      if (status) break;
      const double t=(double)step, k1=2e-3, k2=3e-3;
      const double e1=exp(-k1*t), e2=exp(-k2*t);
      const double expected_a=initial[cell]*e1;
      const double expected_b=initial[cell]*k1/(k2-k1)*(e1-e2);
      const double expected_c=initial[cell]*(1.0+(k1*e2-k2*e1)/(k2-k1));
      const double actual[]={states[cell][0]+states[cell][1],
                              states[cell][2],states[cell][3]};
      const double expected[]={expected_a,expected_b,expected_c};
      for (size_t i=0;i<3;++i)
        if (fabs(actual[i]-expected[i])>1e-3) status=2;
    }
  }
  chem_solver_workspace_destroy(&workspace);
  return status;
}

static int check_analytical_robertson(chem_rosenbrock_parameter_set scheme) {
  const chem_species species[] = {{.name="A"},{.name="B"},{.name="C"}};
  const chem_term a[] = {{.species=0,.coefficient=1}};
  const chem_term b[] = {{.species=1,.coefficient=1}};
  const chem_term bb[] = {{.species=1,.coefficient=2}};
  const chem_term bc[] = {{.species=1,.coefficient=1},
                          {.species=2,.coefficient=1}};
  const chem_term bc_product[] = {{.species=1,.coefficient=1},
                                  {.species=2,.coefficient=1}};
  const chem_term ac[] = {{.species=0,.coefficient=1},
                          {.species=2,.coefficient=1}};
  const chem_reaction reactions[] = {
      {.kind=CHEM_USER_DEFINED,.reactants=a,.reactant_count=1,
       .products=b,.product_count=1,.p={[0]=1,[7]=0}},
      {.kind=CHEM_USER_DEFINED,.reactants=bb,.reactant_count=1,
       .products=bc_product,.product_count=2,.p={[0]=1,[7]=1}},
      {.kind=CHEM_USER_DEFINED,.reactants=bc,.reactant_count=2,
       .products=ac,.product_count=2,.p={[0]=1,[7]=2}}};
  const char *labels[]={"r1","r2","r3"};
  const chem_mechanism mechanism={.species=species,.species_count=3,
      .reactions=reactions,.reaction_count=3,.custom_labels=labels,
      .custom_parameter_count=3};
  const double rates[]={0.04,3e7,1e4};
  const double reference[13][3]={
      {1,0,0},
      {0.9664597373330035,0.00003074626578578675,0.03350951640121071},
      {0.8413699238414729,0.00001623390937990473,0.1586138422491472},
      {0.6172348823960878,0.000006153591274639123,0.3827589640126376},
      {0.3368745306607069,0.000002013702318261393,0.6631234556369748},
      {0.1073004285378040,0.0000004800166972571660,0.8926990914454987},
      {0.01786592114209946,0.00000007274751468436319,0.9821340061103859},
      {0.002031483924973415,0.000000008142277783356159,0.9979685079327488},
      {0.0002076093439016395,0.0000000008306077485067610,0.9997923898254906},
      {0.00002082417512179466,0.00000000008329841429908955,0.9999791757415798},
      {0.000002083229471647004,0.000000000008332935037760723,0.9999979167621954},
      {0.0000002083328471883087,0.0000000000008333315602809495,0.9999997916663195},
      {0.00000002083340149701284,0.0000000000000833336077033474,0.9999999791665152}};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  chem_solver_options options=chem_solver_default_options();
  options.relative_tolerance=1e-10;
  options.absolute_tolerance=1e-12;
  options.max_steps=1000000;
  const chem_rosenbrock_parameters parameters=
      chem_rosenbrock_builtin_parameters(scheme);
  double state[]={1,0,0};
  int status=0;
  double target=1.0, current=0.0;
  for (size_t point=1;point<13 && !status;++point) {
    status=chem_rosenbrock_solve_with_parameters(&mechanism,&workspace,state,
        272.5,101253.3,1e6,rates,target-current,&parameters,&options,NULL);
    if (status) break;
    for (size_t i=0;i<3;++i) {
      const double error=fabs(state[i]-reference[point][i]);
      const double relative=reference[point][i]>0.0
          ? error/reference[point][i] : error;
      if (error>=3e-4 && relative>=1e-6) status=2;
    }
    current=target;
    target*=10.0;
  }
  chem_solver_workspace_destroy(&workspace);
  return status;
}

static int check_analytical_hires(int backward_euler) {
  const chem_species species[]={{.name="Y0"},{.name="Y1"},{.name="Y2"},
      {.name="Y3"},{.name="Y4"},{.name="Y5"},{.name="Y6"},{.name="Y7"}};
  const chem_term y0[]={{.species=0,.coefficient=1}};
  const chem_term y1[]={{.species=1,.coefficient=1}};
  const chem_term y2[]={{.species=2,.coefficient=1}};
  const chem_term y3[]={{.species=3,.coefficient=1}};
  const chem_term y4[]={{.species=4,.coefficient=1}};
  const chem_term y5[]={{.species=5,.coefficient=1}};
  const chem_term y6[]={{.species=6,.coefficient=1}};
  const chem_term y1_products[]={{.species=1,.coefficient=1}};
  const chem_term y2_y5_products[]={{.species=2,.coefficient=0.43/1.12},
                                      {.species=5,.coefficient=0.69/1.12}};
  const chem_term y2_y3_products[]={{.species=0,.coefficient=8.32/10.03},
                                      {.species=3,.coefficient=1.71/10.03}};
  const chem_term y0_y3_products[]={{.species=0,.coefficient=0.43/8.75},
                                      {.species=3,.coefficient=8.32/8.75}};
  const chem_term y0_product[]={{.species=0,.coefficient=1}};
  const chem_term y4_product[]={{.species=4,.coefficient=1}};
  const chem_term y4_y5_products[]={{.species=4,.coefficient=0.43/1.81},
                                      {.species=5,.coefficient=0.69/1.81},
                                      {.species=7,.coefficient=1}};
  const chem_term y5_y7[]={{.species=5,.coefficient=1},
                            {.species=7,.coefficient=1}};
  const chem_term y6_product[]={{.species=6,.coefficient=1}};
  const chem_reaction reactions[]={
      {.kind=CHEM_USER_DEFINED,.reactants=y0,.reactant_count=1,
       .products=y1_products,.product_count=1,.p={[0]=1,[7]=0}},
      {.kind=CHEM_USER_DEFINED,.reactants=y1,.reactant_count=1,
       .products=y0_y3_products,.product_count=2,.p={[0]=1,[7]=1}},
      {.kind=CHEM_USER_DEFINED,.reactants=y2,.reactant_count=1,
       .products=y2_y3_products,.product_count=2,.p={[0]=1,[7]=2}},
      {.kind=CHEM_USER_DEFINED,.products=y0_product,.product_count=1,
       .p={[0]=1,[7]=3}},
      {.kind=CHEM_USER_DEFINED,.reactants=y3,.reactant_count=1,
       .products=y2_y5_products,.product_count=2,.p={[0]=1,[7]=4}},
      {.kind=CHEM_USER_DEFINED,.reactants=y4,.reactant_count=1,
       .products=(const chem_term[]){{.species=2,.coefficient=0.035/1.745},
                                     {.species=5,.coefficient=1.71/1.745}},
       .product_count=2,.p={[0]=1,[7]=5}},
      {.kind=CHEM_USER_DEFINED,.reactants=y5,.reactant_count=1,
       .products=y4_product,.product_count=1,.p={[0]=1,[7]=6}},
      {.kind=CHEM_USER_DEFINED,.reactants=y6,.reactant_count=1,
       .products=y4_y5_products,.product_count=3,.p={[0]=1,[7]=7}},
      {.kind=CHEM_USER_DEFINED,.reactants=y5_y7,.reactant_count=2,
       .products=y6_product,.product_count=1,.p={[0]=1,[7]=8}}};
  const char *labels[]={"r1","r2","r3","r4","r5","r6","r7","r8","r9"};
  const chem_mechanism mechanism={.species=species,.species_count=8,
      .reactions=reactions,.reaction_count=9,.custom_labels=labels,
      .custom_parameter_count=9};
  const double rates[]={1.71,8.75,10.03,0.0007,1.12,1.745,0.43,1.81,280.0};
  const double reference[2][8]={
      {0.000737131257332567,0.000144248572631618,0.000058887297409676,
       0.001175651343283149,0.002386356198831330,0.006238968252742796,
       0.002849998395185769,0.002850001604814231},
      {0.000670305503581864,0.000130996846986347,0.000046862231597733,
       0.001044668020551705,0.000594883830951485,0.001399628833942774,
       0.001014492757718480,0.004685507242281520}};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  double state[]={1,0,0,0,0,0,0,0.0057};
  chem_solver_options ros_options=chem_solver_default_options();
  ros_options.relative_tolerance=1e-6;
  ros_options.absolute_tolerance=1e-8;
  ros_options.max_steps=100000;
  const chem_rosenbrock_parameters ros_parameters=
      chem_rosenbrock_builtin_parameters(CHEM_ROSENBROCK_THREE_STAGE);
  chem_backward_euler_options be_options=chem_backward_euler_default_options();
  be_options.initial_step=321.8122;
  int status=0;
  for (size_t point=0;point<2 && !status;++point) {
    const double duration=321.8122+100.0*(double)point;
    status=backward_euler
        ? chem_backward_euler_solve_with_workspace(&mechanism,&workspace,state,
            298.15,101325.0,40.0,rates,duration,&be_options,NULL)
        : chem_rosenbrock_solve_with_parameters(&mechanism,&workspace,state,
            298.15,101325.0,40.0,rates,duration,
            &ros_parameters,&ros_options,NULL);
    if (status) break;
    const double tolerance=backward_euler?0.1:1e-7;
    for (size_t i=0;i<8;++i)
      if (fabs(state[i]-reference[point][i])>tolerance) status=2;
    if (point==0) be_options.initial_step=duration+100.0;
  }
  chem_solver_workspace_destroy(&workspace);
  return status;
}

static int check_analytical_oregonator(int backward_euler) {
  const chem_species species[]={{.name="X"},{.name="Y"},{.name="Z"},
                                {.name="P"},{.name="Q"}};
  const chem_term x[]={{.species=0,.coefficient=1}};
  const chem_term y[]={{.species=1,.coefficient=1}};
  const chem_term z[]={{.species=2,.coefficient=1}};
  const chem_term xx[]={{.species=0,.coefficient=2}};
  const chem_term xy[]={{.species=0,.coefficient=1},{.species=1,.coefficient=1}};
  const chem_term x_product[]={{.species=0,.coefficient=1}};
  const chem_term p_product[]={{.species=3,.coefficient=1}};
  const chem_term xz_product[]={{.species=2,.coefficient=1},
                                {.species=0,.coefficient=2}};
  const chem_term q_product[]={{.species=4,.coefficient=1}};
  const chem_term y_product[]={{.species=1,.coefficient=1}};
  const chem_reaction reactions[]={
      {.kind=CHEM_USER_DEFINED,.reactants=y,.reactant_count=1,
       .products=x_product,.product_count=1,.p={[0]=1,[7]=0}},
      {.kind=CHEM_USER_DEFINED,.reactants=xy,.reactant_count=2,
       .products=p_product,.product_count=1,.p={[0]=1,[7]=1}},
      {.kind=CHEM_USER_DEFINED,.reactants=x,.reactant_count=1,
       .products=xz_product,.product_count=2,.p={[0]=1,[7]=2}},
      {.kind=CHEM_USER_DEFINED,.reactants=xx,.reactant_count=1,
       .products=q_product,.product_count=1,.p={[0]=1,[7]=3}},
      {.kind=CHEM_USER_DEFINED,.reactants=z,.reactant_count=1,
       .products=y_product,.product_count=1,.p={[0]=1,[7]=4}}};
  const char *labels[]={"r1","r2","r3","r4","r5"};
  const chem_mechanism mechanism={.species=species,.species_count=5,
      .reactions=reactions,.reaction_count=5,.custom_labels=labels,
      .custom_parameter_count=5};
  const double s=77.27,q=8.375e-6,w=0.161,k2=1.6e9,k5=1.0;
  const double k1=k5/(s*w), k3=s*k5/w, k4=s*s*q*k2/2.0;
  const double alpha=k5/(s*w*k2), eta=s*k5/(w*k2);
  const double rho=k5*k5/(w*w*k2), tau=w/k5;
  const double rates[]={k1,k2,k3,k4,k5};
  const double references[13][3]={
      {1,2,3},
      {1.000661467180497,1512.778937348249,10358.54312767229},
      {1.000874625199626,1144.336972384497,83.72149966624639},
      {1.001890368438751,529.9926232295553,1.662279579042420},
      {1.004118022612645,243.8326079910346,1.008822224048647},
      {1.008995416634061,112.1664388662539,1.007783229065319},
      {1.019763472537298,51.59761322947535,1.016985778956376},
      {1.043985088527474,23.73442027531524,1.037691843544522},
      {1.100849071667922,10.91533805469022,1.085831969810860},
      {1.249102130020572,5.013945178605446,1.208326626237875},
      {1.779724751937019,2.281852385542403,1.613754023671725},
      {1.000889326903503,1125.438585746596,16410.49483777168},
      {1.000814870318523,1228.178521549889,132.0554942846513}};
  const double absolute_tolerances[]={alpha*1e-11,eta*1e-11,
                                      rho*1e-11,eta,eta};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  double state[]={alpha,2.0*eta,3.0*rho,0,0};
  chem_solver_options ros_options=chem_solver_default_options();
  ros_options.relative_tolerance=1e-9;
  ros_options.absolute_tolerance=1e-3;
  ros_options.absolute_tolerances=absolute_tolerances;
  ros_options.absolute_tolerance_count=5;
  ros_options.max_steps=100000;
  const chem_rosenbrock_parameters ros_parameters=
      chem_rosenbrock_builtin_parameters(CHEM_ROSENBROCK_THREE_STAGE);
  chem_backward_euler_options be_options=chem_backward_euler_default_options();
  be_options.relative_tolerance=1e-9;
  be_options.absolute_tolerance=1e-3;
  be_options.absolute_tolerances=absolute_tolerances;
  be_options.absolute_tolerance_count=5;
  const double interval=30.0*tau;
  const size_t substeps=backward_euler?18000:1;
  const double dt=interval/(double)substeps;
  int status=0;
  for (size_t point=1;point<=12 && !status;++point) {
    for (size_t substep=0;substep<substeps && !status;++substep)
      status=backward_euler
          ? chem_backward_euler_solve_with_workspace(&mechanism,&workspace,state,
              298.15,101325.0,40.0,rates,dt,&be_options,NULL)
          : chem_rosenbrock_solve_with_parameters(&mechanism,&workspace,state,
              298.15,101325.0,40.0,rates,interval,&ros_parameters,
              &ros_options,NULL);
    if (status) break;
    const double expected[]={alpha*references[point][0],
                             eta*references[point][1],
                             rho*references[point][2]};
    const double tolerance=backward_euler?5e-3:1e-6;
    for (size_t i=0;i<3;++i)
      if (fabs(state[i]-expected[i])>tolerance*fabs(expected[i])) status=2;
  }
  chem_solver_workspace_destroy(&workspace);
  return status;
}

static int check_analytical_e5(void) {
  const chem_species species[]={{.name="A1"},{.name="A2"},{.name="A3"},
      {.name="A4"},{.name="A5"},{.name="A6"}};
  const chem_term a1[]={{.species=0,.coefficient=1}};
  const chem_term a2a3[]={{.species=1,.coefficient=1},
                           {.species=2,.coefficient=1}};
  const chem_term a1a3[]={{.species=0,.coefficient=1},
                           {.species=2,.coefficient=1}};
  const chem_term a2a3_products[]={{.species=1,.coefficient=1},
                                   {.species=2,.coefficient=1}};
  const chem_term a5[]={{.species=4,.coefficient=1}};
  const chem_term a4[]={{.species=3,.coefficient=1}};
  const chem_term a3a6_products[]={{.species=2,.coefficient=1},
                                   {.species=5,.coefficient=1}};
  const chem_reaction reactions[]={
      {.kind=CHEM_USER_DEFINED,.reactants=a1,.reactant_count=1,
       .products=a2a3_products,.product_count=2,.p={[0]=1,[7]=0}},
      {.kind=CHEM_USER_DEFINED,.reactants=a2a3,.reactant_count=2,
       .products=a5,.product_count=1,.p={[0]=1,[7]=1}},
      {.kind=CHEM_USER_DEFINED,.reactants=a1a3,.reactant_count=2,
       .products=a4,.product_count=1,.p={[0]=1,[7]=2}},
      {.kind=CHEM_USER_DEFINED,.reactants=a4,.reactant_count=1,
       .products=a3a6_products,.product_count=2,.p={[0]=1,[7]=3}}};
  const char *labels[]={"r1","r2","r3","r4"};
  const chem_mechanism mechanism={.species=species,.species_count=6,
      .reactions=reactions,.reaction_count=4,.custom_labels=labels,
      .custom_parameter_count=4};
  const double rates[]={7.89e-10,1.13e9,1.1e7,1.13e3};
  const double reference[8][4]={
      {1.76e-3,0,0,0},
      {1.7599259497677897058e-3,1.3846281519376516449e-11,
       7.6370038530073911180e-13,1.3082581134075777338e-11},
      {1.6180769999072942552e-3,1.3822370304983735443e-10,
       8.2515735006838336088e-12,1.2997212954915352082e-10},
      {7.4813208224292220114e-6,2.3734781561205975019e-12,
       2.2123586689581663654e-12,1.6111948716243113653e-13},
      {4.7150333630401632232e-10,1.8188895860807021729e-14,
       1.8188812376786725407e-14,8.3484020296321693074e-20},
      {3.1317148329356996037e-14,1.4840957952870064294e-16,
       1.4840957948345691466e-16,4.5243728279782625194e-26},
      {3.8139035189787091771e-49,1.0192582567660293322e-20,
       1.0192582567660293322e-20,3.7844935507486221171e-65},
      {0,8.8612334976263783420e-23,8.8612334976263783421e-23,0}};
  const double atols[]={1e-7,1e-17,1e-17,1e-17,1e-7,1e-7};
  chem_solver_workspace workspace={0};
  if (chem_solver_workspace_init(&workspace,&mechanism)) return 1;
  double state[]={1.76e-3,0,0,0,0,0};
  chem_solver_options options=chem_solver_default_options();
  options.relative_tolerance=1e-13;
  options.absolute_tolerances=atols;
  options.absolute_tolerance_count=6;
  options.max_steps=1000000;
  const chem_rosenbrock_parameters parameters=
      chem_rosenbrock_builtin_parameters(CHEM_ROSENBROCK_THREE_STAGE);
  double target=10.0,current=0.0;
  int status=0;
  for (size_t point=1;point<8 && !status;++point) {
    status=chem_rosenbrock_solve_with_parameters(&mechanism,&workspace,state,
        298.15,101325.0,40.0,rates,target-current,&parameters,&options,NULL);
    if (status) break;
    for (size_t i=0;i<4;++i) {
      const double error=fabs(state[i]-reference[point][i]);
      const double abs_tol=i==0?1e-6:1e-13;
      const double relative=reference[point][i]>0.0
          ? error/reference[point][i] : error;
      if (error>=abs_tol && relative>=1e-10) status=2;
    }
    current=target;
    target*=100.0;
  }
  chem_solver_workspace_destroy(&workspace);
  return status;
}

int main(void) {
  if (check_rosenbrock_species_tolerances()) return 52;
  { const int status=check_rejected_step_determinism();
    if (status) return 53+status; }
  { const int status=check_stiff_photolysis_chain(0);
    if (status) return 60+status; }
  { const int status=check_stiff_photolysis_chain(1);
    if (status) return 70+status; }
  const chem_rosenbrock_parameter_set robertson_schemes[]={
      CHEM_ROSENBROCK_TWO_STAGE,CHEM_ROSENBROCK_THREE_STAGE,
      CHEM_ROSENBROCK_FOUR_STAGE,CHEM_ROSENBROCK_FOUR_STAGE_DAE,
      CHEM_ROSENBROCK_SIX_STAGE_DAE};
  for (size_t i=0;i<sizeof(robertson_schemes)/sizeof(robertson_schemes[0]);++i) {
    const int status=check_analytical_robertson(robertson_schemes[i]);
    if (status) return 80+(int)i*10+status;
  }
  if (check_analytical_hires(0)) return 90;
  if (check_analytical_hires(1)) return 91;
  if (check_analytical_oregonator(0)) return 92;
  if (check_analytical_oregonator(1)) return 93;
  if (check_analytical_e5()) return 94;
  const double t = 301.24;
  const double a = chem_arrhenius(1.0, 0.0, 0.0, 300.0, 0.0, t, 0.0);
  const double b = chem_arrhenius(3.0e-11, 0.0, -200.0, 300.0, 0.0, t, 0.0);
  const double e1 = 1.0;
  const double e2 = 3.0e-11 * exp(-200.0 / t);
  if (fabs(a - e1) > 1e-14 * e1) return 1;
  if (fabs(b - e2) > 1e-14 * e2) return 2;

  const double temperature = 301.24;
  const double pressure = 98765.0;
  const double air_density = pressure / (8.31446261815324 * temperature);
  chem_reaction formula = {.p = {1.0, 0.0, 0.0, 300.0, 0.0}};
  formula.kind = CHEM_REVERSIBLE;
  formula.p[0] = 1.0; formula.p[1] = 0.0; formula.p[2] = 0.0;
  if (chem_reaction_rate(&formula, temperature, air_density, NULL) != 0.0) return 16;
  formula.p[0] = 1.14e-2; formula.p[1] = 2300.0; formula.p[2] = 0.32;
  if (!near(chem_reaction_rate(&formula, temperature, air_density, NULL),
            1.14e-2 * exp(2300.0 / temperature) * 0.32)) return 9;
  formula.p[0] = 1.5; formula.p[1] = 2000.0; formula.p[2] = 0.5;
  for (size_t i = 0; i < 3; ++i) {
    const double t = 273.15 + 25.0 * (double)i;
    if (!near(chem_reaction_rate(&formula, t, air_density, NULL),
              1.5 * exp(2000.0 / t) * 0.5)) return 17;
  }
  formula.p[0] = 3.5; formula.p[1] = 1800.0; formula.p[2] = 1.0;
  if (!near(chem_reaction_rate(&formula, 298.15, air_density, NULL),
            3.5 * exp(1800.0 / 298.15))) return 18;
  formula.p[0] = 5.0; formula.p[1] = 1000.0; formula.p[2] = 0.0;
  if (chem_reaction_rate(&formula, 298.15, air_density, NULL) != 0.0) return 19;

  formula.kind = CHEM_TAYLOR_SERIES;
  formula.p[0] = 3.0e-11; formula.p[1] = 0.0; formula.p[2] = -200.0;
  formula.p[3] = 300.0; formula.p[4] = 0.0; formula.p[5] = 3.0;
  formula.p[6] = 12.5; formula.p[7] = 1.3e-2; formula.p[8] = 5.2e-4;
  double polynomial = 12.5 + 1.3e-2 * temperature + 5.2e-4 * temperature * temperature;
  double expected_taylor = 3.0e-11 * exp(-200.0 / temperature) * polynomial;
  if (!near(chem_reaction_rate(&formula, temperature, air_density, NULL),
            expected_taylor)) return 10;
  double constants[1];
  chem_mechanism formula_mechanism = {.reactions = &formula, .reaction_count = 1};
  formula.p[4] = 1e-8;
  expected_taylor = polynomial * formula.p[0] * exp(formula.p[2] / temperature) *
      pow(temperature / formula.p[3], formula.p[1]) * (1.0 + formula.p[4] * pressure);
  if (chem_rate_constants(&formula_mechanism, temperature, pressure, air_density,
                          NULL, constants) != 0 || !near(constants[0], expected_taylor)) return 11;

  formula.kind = CHEM_TERNARY_ACTIVATION_JPL19;
  formula.p[0] = 1.2e-30; formula.p[1] = -2.3; formula.p[2] = -302.3; formula.p[3] = 298.0;
  formula.p[4] = 2.6e-12; formula.p[5] = -3.1; formula.p[6] = -402.1; formula.p[7] = 310.0;
  formula.p[8] = 3.0e-13; formula.p[9] = 0.0; formula.p[10] = -100.0; formula.p[11] = 300.0;
  formula.p[12] = 0.9; formula.p[13] = 1.2;
  double k0 = formula.p[0] * exp(formula.p[2] / temperature) * pow(temperature / formula.p[3], formula.p[1]);
  double kinf = formula.p[4] * exp(formula.p[6] / temperature) * pow(temperature / formula.p[7], formula.p[5]);
  double kint = formula.p[8] * exp(formula.p[10] / temperature) * pow(temperature / formula.p[11], formula.p[9]);
  double ratio = k0 * air_density / kinf;
  double kf = k0 * air_density / (1.0 + ratio) * pow(formula.p[12], formula.p[13] / (formula.p[13] + pow(log10(ratio), 2.0)));
  if (!near(chem_reaction_rate(&formula, temperature, air_density, NULL),
            kf + kint * (1.0 - kf / kinf))) return 12;
  formula.p[8] = 0.0;
  const double jpl_zero_kint = chem_reaction_rate(&formula, temperature, air_density, NULL);
  formula.kind = CHEM_TROE;
  formula.p[0] = 1.2e-30; formula.p[1] = -2.3; formula.p[2] = -302.3;
  formula.p[3] = 2.6e-12; formula.p[4] = -3.1; formula.p[5] = -402.1;
  formula.p[6] = 0.9; formula.p[7] = 1.2;
  if (!near(jpl_zero_kint, chem_reaction_rate(&formula, temperature, air_density, NULL))) return 20;

  formula.kind = CHEM_LAMBDA;
  formula.rate_function = test_lambda;
  if (!near(chem_reaction_rate(&formula, temperature, air_density, NULL),
            2.0 * temperature)) return 13;

  formula.kind = CHEM_BRANCHED_ALKOXY;
  formula.p[0] = 1.2; formula.p[1] = 204.3; formula.p[2] = 1.0e-3;
  formula.p[3] = 2.0; formula.p[4] = 0.25;
  double branched_a = 2.0e-22 * 6.02214076e23 * 1.0e-6 * exp(2.0) * air_density;
  double branched_b = 0.43 * pow(temperature / 298.0, -8.0);
  double branched_av = branched_a / (1.0 + branched_a / branched_b) *
      pow(0.41, 1.0 / (1.0 + pow(log10(branched_a / branched_b), 2.0)));
  if (!near(chem_reaction_rate(&formula, temperature, air_density, NULL),
            1.2 * exp(-204.3 / temperature) * 0.25 / (0.25 + branched_av))) return 14;
  formula.kind = CHEM_BRANCHED_NITRATE;
  if (!near(chem_reaction_rate(&formula, temperature, air_density, NULL),
            1.2 * exp(-204.3 / temperature) * branched_av / (branched_av + 0.25))) return 15;

  formula.kind = CHEM_TROE;
  formula.p[0] = 1.2; formula.p[1] = 2.3; formula.p[2] = 302.3;
  formula.p[3] = 2.6; formula.p[4] = -3.1; formula.p[5] = 402.1;
  formula.p[6] = 0.9; formula.p[7] = 1.2;
  double troe_k0 = 1.2 * exp(302.3 / temperature) * pow(temperature / 300.0, 2.3);
  double troe_kinf = 2.6 * exp(402.1 / temperature) * pow(temperature / 300.0, -3.1);
  double troe_pr = air_density * troe_k0 / troe_kinf;
  double expected_troe = troe_k0 * air_density / (1.0 + troe_pr) *
      pow(0.9, 1.2 / (1.2 + pow(log10(troe_pr), 2.0)));
  if (!near(chem_reaction_rate(&formula, temperature, air_density, NULL), expected_troe)) return 21;
  formula.kind = CHEM_TERNARY_ACTIVATION;
  double expected_ternary = troe_k0 / (1.0 + troe_pr) *
      pow(0.9, 1.2 / (1.2 + pow(log10(troe_pr), 2.0)));
  if (!near(chem_reaction_rate(&formula, temperature, air_density, NULL), expected_ternary)) return 22;

  formula.kind = CHEM_TUNNELING;
  formula.p[0] = 1.2; formula.p[1] = 2.3; formula.p[2] = 302.3;
  if (!near(chem_reaction_rate(&formula, temperature, air_density, NULL),
            1.2 * exp(-2.3 / temperature + 302.3 / pow(temperature, 3.0)))) return 23;
  formula.p[0] = 1.0; formula.p[1] = 0.0; formula.p[2] = 0.0;
  if (chem_reaction_rate(&formula, temperature, air_density, NULL) != 1.0) return 24;

  formula.kind = CHEM_USER_DEFINED;
  formula.p[0] = 2.0; formula.p[7] = 0.0;
  const double user_parameter[] = {1.2};
  if (!near(chem_reaction_rate(&formula, temperature, air_density, user_parameter), 2.4)) return 25;

  const chem_species surface_species[] = {
      {.name = "aerosol", .molecular_weight = 0.025, .diffusion = 230.0}};
  formula.kind = CHEM_SURFACE;
  formula.surface_species = 0;
  formula.p[0] = 1.0;
  formula.p[7] = 0.0;
  formula_mechanism.species = surface_species;
  formula_mechanism.species_count = 1;
  const double surface_params[] = {1.0e-7, 2.5e6};
  const double surface_t = 273.65;
  const double surface_mean_speed = sqrt(8.0 * 8.31446261815324 * surface_t /
      (3.14159265358979323846 * 0.025));
  double expected_surface = 4.0 * 2.5e6 * 3.14159265358979323846 * 1e-14 /
      (1e-7 / 230.0 + 4.0 / surface_mean_speed);
  if (chem_rate_constants(&formula_mechanism, surface_t, pressure, air_density,
                          surface_params, constants) != 0 ||
      !near(constants[0], expected_surface)) return 26;
  formula.p[0] = 0.74;
  expected_surface = 4.0 * 2.5e6 * 3.14159265358979323846 * 1e-14 /
      (1e-7 / 230.0 + 4.0 / (surface_mean_speed * 0.74));
  if (chem_rate_constants(&formula_mechanism, surface_t, pressure, air_density,
                          surface_params, constants) != 0 ||
      !near(constants[0], expected_surface)) return 27;
  const chem_species missing_diffusion[] = {{.name = "aerosol", .molecular_weight = 0.025}};
  formula_mechanism.species = missing_diffusion;
  if (chem_rate_constants(&formula_mechanism, surface_t, pressure, air_density,
                          surface_params, constants) == 0) return 28;
  const chem_species missing_weight[] = {{.name = "aerosol", .diffusion = 230.0}};
  formula_mechanism.species = missing_weight;
  if (chem_rate_constants(&formula_mechanism, surface_t, pressure, air_density,
                          surface_params, constants) == 0) return 29;

  const chem_species species[] = {{.name = "A"}, {.name = "B"}};
  const chem_term reactants[] = {{.species = 0, .coefficient = 1.0}};
  const chem_term products[] = {{.species = 1, .coefficient = 1.0}};
  const chem_reaction reaction = {
      .kind = CHEM_USER_DEFINED,
      .reactants = reactants,
      .reactant_count = 1,
      .products = products,
      .product_count = 1,
      .p = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  };
  const chem_mechanism mechanism = {
      .species = species,
      .species_count = 2,
      .reactions = &reaction,
      .reaction_count = 1,
  };
  const double concentration[] = {3.0, 0.0};
  const double custom_rates[] = {2.0};
  double forcing[2];
  if (chem_forcing(&mechanism, concentration, 300.0, 1.0, custom_rates,
                   forcing) != 0) return 3;
  if (forcing[0] != -6.0 || forcing[1] != 6.0) return 4;

  const chem_species reaction_species[] = {
      {.name = "O3"}, {.name = "NO"}, {.name = "NO2"}, {.name = "O2"}};
  const chem_term built_reactants[] = {
      {.species = 0, .coefficient = 1.0}, {.species = 1, .coefficient = 1.0}};
  const chem_term built_products[] = {
      {.species = 2, .coefficient = 1.0}, {.species = 3, .coefficient = 1.0}};
  const chem_reaction built_reaction = {
      .kind = CHEM_ARRHENIUS, .label = "O3 + NO -> NO2 + O2",
      .reactants = built_reactants, .reactant_count = 2,
      .products = built_products, .product_count = 2,
      .p = {1.0, 0.0, 0.0, 300.0, 0.0}};
  const chem_mechanism built_mechanism = {
      .species = reaction_species, .species_count = 4,
      .reactions = &built_reaction, .reaction_count = 1};
  if (built_mechanism.reaction_count != 1 ||
      built_mechanism.reactions[0].reactant_count != 2 ||
      built_mechanism.reactions[0].product_count != 2 ||
      built_mechanism.species[built_mechanism.reactions[0].products[0].species].name[2] != '2') return 30;
  const chem_reaction copied_reaction = built_reaction;
  if (copied_reaction.reactants != built_reaction.reactants ||
      copied_reaction.products != built_reaction.products ||
      copied_reaction.kind != built_reaction.kind ||
      copied_reaction.p[0] != built_reaction.p[0] ||
      copied_reaction.reactants[1].species != built_reaction.reactants[1].species) return 31;

  const chem_species parameterized_species[] = {
      {.name = "a"},
      {.name = "b", .density_parameterized = 1, .density_factor = 2.0}};
  const chem_term react_a[] = {{.species = 0, .coefficient = 1.0}};
  const chem_term react_b[] = {{.species = 1, .coefficient = 1.0}};
  const chem_term product_a[] = {{.species = 0, .coefficient = 1.0}};
  const chem_term product_b[] = {{.species = 1, .coefficient = 1.0}};
  const chem_reaction parameterized_reactions[] = {
      {.kind = CHEM_ARRHENIUS, .reactants = react_a, .reactant_count = 1,
       .products = product_b, .product_count = 1,
       .p = {1.0, 0.0, 0.0, 300.0, 0.0}},
      {.kind = CHEM_ARRHENIUS, .reactants = react_b, .reactant_count = 1,
       .products = product_a, .product_count = 1,
       .p = {1.0, 0.0, 0.0, 300.0, 0.0}}};
  const chem_mechanism parameterized_mechanism = {
      .species = parameterized_species, .species_count = 2,
      .reactions = parameterized_reactions, .reaction_count = 2};
  double parameterized_rates[2], parameterized_forcing[2], parameterized_jacobian[4];
  if (chem_rate_constants(&parameterized_mechanism, 300.0, 101325.0, 5.0,
                          NULL, parameterized_rates) != 0 ||
      parameterized_rates[0] != 1.0 || parameterized_rates[1] != 10.0) return 32;
  const double parameterized_state[] = {1.0, 0.0};
  if (chem_forcing_with_constants(&parameterized_mechanism, parameterized_state,
                                  parameterized_rates, 5.0,
                                  parameterized_forcing) != 0 ||
      parameterized_forcing[0] != 9.0 || parameterized_forcing[1] != -9.0) return 33;
  if (chem_negative_jacobian_with_constants(&parameterized_mechanism,
          parameterized_state, parameterized_rates, parameterized_jacobian) != 0 ||
      parameterized_jacobian[0] != 1.0 || parameterized_jacobian[2] != -1.0 ||
      parameterized_jacobian[1] != 0.0 || parameterized_jacobian[3] != 0.0) return 34;

  const chem_term repeated_reactants[] = {
      {.species = 0, .coefficient = 1.0},
      {.species = 0, .coefficient = 1.0}};
  const chem_term repeated_reactant_coefficient[] = {
      {.species = 0, .coefficient = 2.0}};
  const chem_term repeated_products[] = {
      {.species = 0, .coefficient = 1.0},
      {.species = 1, .coefficient = 1.0}};
  const chem_reaction repeated_form_reactions[] = {
      {.kind = CHEM_ARRHENIUS, .reactants = repeated_reactants,
       .reactant_count = 2, .products = repeated_products, .product_count = 2,
       .p = {3.0, 0.0, 0.0, 300.0, 0.0}},
      {.kind = CHEM_ARRHENIUS, .reactants = repeated_reactant_coefficient,
       .reactant_count = 1, .products = repeated_products, .product_count = 2,
       .p = {3.0, 0.0, 0.0, 300.0, 0.0}}};
  const chem_mechanism repeated_mechanisms[] = {
      {.species = species, .species_count = 2,
       .reactions = &repeated_form_reactions[0], .reaction_count = 1},
      {.species = species, .species_count = 2,
       .reactions = &repeated_form_reactions[1], .reaction_count = 1}};
  const double repeated_state[] = {2.0, 0.0};
  double repeated_rates[2], repeated_jacobians[8];
  for (size_t form = 0; form < 2; ++form) {
    if (chem_rate_constants(&repeated_mechanisms[form], 300.0, 101325.0,
                            0.0, NULL, &repeated_rates[form]) != 0 ||
        chem_negative_jacobian_with_constants(&repeated_mechanisms[form],
            repeated_state, &repeated_rates[form],
            &repeated_jacobians[4*form]) != 0) return 37;
  }
  for (size_t i = 0; i < 4; ++i)
    if (repeated_jacobians[i] != repeated_jacobians[4+i]) return 38;
  if (repeated_jacobians[0] != 12.0 || repeated_jacobians[2] != -12.0 ||
      repeated_jacobians[1] != 0.0 || repeated_jacobians[3] != 0.0) return 39;

  const chem_reaction decay = {
      .kind = CHEM_ARRHENIUS,
      .reactants = reactants,
      .reactant_count = 1,
      .products = products,
      .product_count = 1,
      .p = {0.01, 0.0, 0.0, 300.0, 0.0},
  };
  const chem_mechanism decay_mechanism = {
      .species = species,
      .species_count = 2,
      .reactions = &decay,
      .reaction_count = 1,
  };
  double y[] = {1.0, 0.0};
  chem_solver_options options = {
      .relative_tolerance=1.0e-10, .absolute_tolerance=1.0e-12,
      .max_steps=1000};
  chem_solver_stats stats;
  if (chem_rosenbrock3_solve(&decay_mechanism, y, 300.0, 0.0, NULL, 1.0,
                             &options, &stats) != 0) return 5;
  if (fabs(y[0] - exp(-0.01)) > 1.0e-9 || fabs(y[0] + y[1] - 1.0) > 1.0e-10)
    return 6;

  const chem_species species3[] = {{.name = "A"}, {.name = "B"}, {.name = "C"}};
  const chem_term rb[] = {{.species = 1, .coefficient = 1.0}};
  const chem_term pc[] = {{.species = 2, .coefficient = 1.0}};
  const chem_reaction chain[] = {
      {.kind = CHEM_ARRHENIUS, .reactants = reactants, .reactant_count = 1,
       .products = products, .product_count = 1,
       .p = {2.0, 0.0, 0.0, 300.0, 0.0}},
      {.kind = CHEM_ARRHENIUS, .reactants = rb, .reactant_count = 1,
       .products = pc, .product_count = 1,
       .p = {3.0, 0.0, 0.0, 300.0, 0.0}},
  };
  const chem_mechanism chain_mechanism = {
      .species = species3, .species_count = 3,
      .reactions = chain, .reaction_count = 2,
  };
  double chain_y[] = {1.0, 0.0, 0.0};
  const chem_solver_options chain_options = {
      .relative_tolerance=1.0e-10, .absolute_tolerance=1.0e-12,
      .max_steps=10000};
  int chain_status = chem_rosenbrock3_solve(&chain_mechanism, chain_y, 300.0,
                                             0.0, NULL, 1.0, &chain_options, &stats);
  if (chain_status != 0) return 7;
  const double exact_a = exp(-2.0);
  const double exact_b = 2.0*(exp(-2.0) - exp(-3.0));
  if (fabs(chain_y[0] - exact_a) > 1.0e-8 ||
      fabs(chain_y[1] - exact_b) > 1.0e-8 ||
      fabs(chain_y[0] + chain_y[1] + chain_y[2] - 1.0) > 1.0e-9) return 8;

  /* Port of test_README_example: compare the values at the precision printed
     by the tutorial (three significant digits). */
  const chem_species foo_species[] = {
      {.name = "Foo"}, {.name = "Bar"}, {.name = "Baz"}};
  const chem_term foo_reactant[] = {{.species = 0, .coefficient = 1.0}};
  const chem_term foo_bar_reactants[] = {
      {.species = 0, .coefficient = 1.0}, {.species = 1, .coefficient = 1.0}};
  const chem_term bar_baz_products[] = {
      {.species = 1, .coefficient = 0.8}, {.species = 2, .coefficient = 0.2}};
  const chem_term baz_product[] = {{.species = 2, .coefficient = 1.0}};
  const chem_reaction foo_reactions[] = {
      {.kind = CHEM_ARRHENIUS, .reactants = foo_reactant, .reactant_count = 1,
       .products = bar_baz_products, .product_count = 2,
       .p = {1.0e-3, 0.0, 0.0, 300.0, 0.0}},
      {.kind = CHEM_ARRHENIUS, .reactants = foo_bar_reactants, .reactant_count = 2,
       .products = baz_product, .product_count = 1,
       .p = {1.0e-5, 0.0, 110.0, 300.0, 0.0}},
  };
  const chem_mechanism foo_mechanism = {
      .species = foo_species, .species_count = 3,
      .reactions = foo_reactions, .reaction_count = 2};
  const double tutorial_expected[10][3] = {
      {11.8, 5.90, 1.91}, {6.79, 9.05, 3.32}, {3.83, 10.7, 4.21},
      {2.14, 11.7, 4.74}, {1.19, 12.2, 5.04}, {0.658, 12.4, 5.21},
      {0.364, 12.6, 5.31}, {0.201, 12.7, 5.36}, {0.111, 12.7, 5.39},
      {0.0613, 12.8, 5.41}};
  double tutorial_y[] = {20.0, 0.0, 0.0};
  const double tutorial_t = 287.45, tutorial_p = 101319.9;
  const double tutorial_density = tutorial_p / (8.31446261815324 * tutorial_t);
  chem_solver_workspace tutorial_workspace;
  if (chem_solver_workspace_init(&tutorial_workspace, &foo_mechanism) != 0) return 35;
  int tutorial_error = 0;
  for (size_t step = 0; step < 10 && !tutorial_error; ++step) {
    tutorial_error = chem_rosenbrock3_solve_with_workspace(
        &foo_mechanism, &tutorial_workspace, tutorial_y, tutorial_t, tutorial_p,
        tutorial_density, NULL, 500.0, NULL, NULL);
    for (size_t species_index = 0; species_index < 3 && !tutorial_error;
         ++species_index) {
      const double expected = tutorial_expected[step][species_index];
      const double print_tolerance = 0.5 * pow(10.0, floor(log10(expected)) - 2.0);
      if (fabs(tutorial_y[species_index] - expected) > print_tolerance)
        tutorial_error = 1;
    }
  }
  chem_solver_workspace_destroy(&tutorial_workspace);
  if (tutorial_error) return 36;

  /* Port of test_solver_results: grouped cells share one adaptive controller,
     and the printed solver counters are part of this plain-program case. */
  const chem_term solve_a[] = {{.species = 0, .coefficient = 1.0}};
  const chem_term solve_b[] = {{.species = 1, .coefficient = 1.0}};
  const chem_term solve_bb[] = {{.species = 1, .coefficient = 2.0}};
  const chem_term solve_bc[] = {
      {.species = 1, .coefficient = 1.0},
      {.species = 2, .coefficient = 1.0}};
  const chem_term solve_ac[] = {
      {.species = 0, .coefficient = 1.0},
      {.species = 2, .coefficient = 1.0}};
  const chem_reaction solve_reactions[] = {
      {.kind = CHEM_USER_DEFINED, .label = "r1", .reactants = solve_a,
       .reactant_count = 1, .products = solve_b, .product_count = 1,
       .p = {1.0, 0, 0, 0, 0, 0, 0, 0}},
      {.kind = CHEM_USER_DEFINED, .label = "r2", .reactants = solve_bb,
       .reactant_count = 1, .products = solve_bc, .product_count = 2,
       .p = {1.0, 0, 0, 0, 0, 0, 0, 1}},
      {.kind = CHEM_USER_DEFINED, .label = "r3", .reactants = solve_bc,
       .reactant_count = 2, .products = solve_ac, .product_count = 2,
       .p = {1.0, 0, 0, 0, 0, 0, 0, 2}},
  };
  const char *solve_labels[] = {"r1", "r2", "r3"};
  const chem_mechanism solve_mechanism = {
      .species = species3, .species_count = 3,
      .reactions = solve_reactions, .reaction_count = 3,
      .custom_labels = solve_labels, .custom_parameter_count = 3};
  double solve_states[] = {1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.5, 0.0, 0.0};
  const double solve_temperatures[] = {272.5, 272.5, 272.5};
  const double solve_pressures[] = {101253.3, 101253.3, 101253.3};
  const double solve_densities[] = {1e6, 1e6, 1e6};
  const double solve_rates[] = {
      0.04, 3e7, 1e4, 0.04, 3e7, 1e4, 0.04, 3e7, 1e4};
  chem_solver_workspace solve_workspace = {0};
  chem_solver_stats solve_stats;
  int solve_status = chem_solver_workspace_init_group(
      &solve_workspace, &solve_mechanism, 3);
  if (!solve_status)
    solve_status = chem_rosenbrock3_solve_group_with_workspace(
        &solve_mechanism, &solve_workspace, solve_states, solve_temperatures,
        solve_pressures, solve_densities, solve_rates, 3, 200.0, NULL,
        &solve_stats);
  if (solve_status != 0 || solve_stats.internal_steps != 18 ||
      solve_stats.accepted != 18 || solve_stats.rejected != 0 ||
      solve_stats.function_calls != 36 || solve_stats.jacobian_updates != 18 ||
      solve_stats.decompositions != 18 || solve_stats.solves != 54 ||
      solve_stats.constraint_init_iterations != 0) return 40;

  /* The multiple-grid-cells tutorial prints ten states for each cell. */
  const double grid_expected[10][3][3] = {
      {{5.35e-1, 4.48e-6, 4.65e-1},
       {1.21, 6.00e-6, 7.89e-1}, {2.30e-1, 3.30e-6, 2.70e-1}},
      {{4.50e-1, 3.23e-6, 5.50e-1},
       {1.05, 4.42e-6, 9.45e-1}, {1.85e-1, 2.32e-6, 3.15e-1}},
      {{4.00e-1, 2.63e-6, 6.00e-1},
       {9.59e-1, 3.65e-6, 1.04}, {1.60e-1, 1.85e-6, 3.40e-1}},
      {{3.64e-1, 2.27e-6, 6.36e-1},
       {8.90e-1, 3.18e-6, 1.11}, {1.42e-1, 1.57e-6, 3.58e-1}},
      {{3.37e-1, 2.01e-6, 6.63e-1},
       {8.35e-1, 2.85e-6, 1.16}, {1.29e-1, 1.38e-6, 3.71e-1}},
      {{3.15e-1, 1.82e-6, 6.85e-1},
       {7.91e-1, 2.60e-6, 1.21}, {1.19e-1, 1.24e-6, 3.81e-1}},
      {{2.96e-1, 1.67e-6, 7.04e-1},
       {7.54e-1, 2.40e-6, 1.25}, {1.11e-1, 1.13e-6, 3.89e-1}},
      {{2.81e-1, 1.55e-6, 7.19e-1},
       {7.21e-1, 2.24e-6, 1.28}, {1.04e-1, 1.04e-6, 3.96e-1}},
      {{2.67e-1, 1.45e-6, 7.33e-1},
       {6.93e-1, 2.11e-6, 1.31}, {9.77e-2, 9.65e-7, 4.02e-1}},
      {{2.55e-1, 1.37e-6, 7.45e-1},
       {6.68e-1, 2.00e-6, 1.33}, {9.25e-2, 9.02e-7, 4.07e-1}}};
  int grid_error = 0;
  for (size_t step = 0; step < 10 && !grid_error; ++step) {
    if (step != 0)
      grid_error = chem_rosenbrock3_solve_group_with_workspace(
          &solve_mechanism, &solve_workspace, solve_states, solve_temperatures,
          solve_pressures, solve_densities, solve_rates, 3, 200.0, NULL, NULL);
    for (size_t cell = 0; cell < 3 && !grid_error; ++cell) {
      for (size_t variable = 0; variable < 3 && !grid_error; ++variable) {
        const double expected = grid_expected[step][cell][variable];
        const double print_tolerance =
            0.5 * pow(10.0, floor(log10(expected)) - 2.0);
        if (fabs(solve_states[cell*3 + variable] - expected) > print_tolerance)
          grid_error = 1;
      }
    }
  }
  chem_solver_workspace_destroy(&solve_workspace);
  if (grid_error) return 41;

  /* Port of test_openmp: every worker owns its solver workspace and state. */
  double thread_states[] = {
      1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
  int thread_status[] = {0, 0, 0};
  #pragma omp parallel num_threads(3)
  {
    const int thread = omp_get_thread_num();
    chem_solver_workspace thread_workspace = {0};
    thread_status[thread] = chem_solver_workspace_init(
        &thread_workspace, &solve_mechanism);
    const double thread_rates[] = {0.04, 3e7, 1e4};
    for (size_t step = 0; step < 10 && !thread_status[thread]; ++step)
      thread_status[thread] = chem_rosenbrock3_solve_with_workspace(
          &solve_mechanism, &thread_workspace, thread_states + 3*thread,
          272.5, 101253.3, 1e6, thread_rates, 200.0, NULL, NULL);
    chem_solver_workspace_destroy(&thread_workspace);
  }
  for (size_t thread = 0; thread < 3; ++thread) {
    if (thread_status[thread]) return 42;
    const double thread_expected[] = {2.55e-1, 1.37e-6, 7.45e-1};
    for (size_t variable = 0; variable < 3; ++variable) {
      const double expected = thread_expected[variable];
      const double print_tolerance =
          0.5 * pow(10.0, floor(log10(expected)) - 2.0);
      if (fabs(thread_states[3*thread + variable] - expected) > print_tolerance)
        return 43;
    }
  }

  /* Port the vectorized-matrix tutorial's observable cell and sparse layouts.
     Its Solve helper is never entered because its state starts converged. */
  const double regular_layout[] = {
      1.1, 1.2, 1.3, 2.1, 2.2, 2.3, 3.1, 3.2, 3.3};
  const double vectorized_layout[] = {
      1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.3, 2.3, 3.3};
  for (size_t cell = 0; cell < 3; ++cell)
    for (size_t variable = 0; variable < 3; ++variable)
      if (regular_layout[cell*3 + variable] !=
          vectorized_layout[variable*3 + cell]) return 44;
  const size_t sparse_rows[] = {0, 2, 2, 3};
  const size_t sparse_columns[] = {1, 1, 3, 2};
  int sparse_standard[12], sparse_vectorized[12];
  for (size_t block = 0; block < 3; ++block)
    for (size_t element = 0; element < 4; ++element) {
      const int value = (int)(100*block + 10*sparse_rows[element] +
                              sparse_columns[element]);
      sparse_standard[block*4 + element] = value;
      sparse_vectorized[element*3 + block] = value;
    }
  const int expected_standard[] = {
      1, 21, 23, 32, 101, 121, 123, 132, 201, 221, 223, 232};
  const int expected_vectorized[] = {
      1, 101, 201, 21, 121, 221, 23, 123, 223, 32, 132, 232};
  for (size_t i = 0; i < 12; ++i)
    if (sparse_standard[i] != expected_standard[i] ||
        sparse_vectorized[i] != expected_vectorized[i]) return 45;

  /* Port test_solver_configuration across all five MICM Rosenbrock sets. */
  const chem_rosenbrock_parameter_set tutorial_schemes[] = {
      CHEM_ROSENBROCK_TWO_STAGE, CHEM_ROSENBROCK_THREE_STAGE,
      CHEM_ROSENBROCK_FOUR_STAGE, CHEM_ROSENBROCK_FOUR_STAGE_DAE,
      CHEM_ROSENBROCK_SIX_STAGE_DAE};
  const double scheme_expected[5][10][3] = {
      {{5.37e-1, 4.50e-6, 4.63e-1}, {4.51e-1, 3.23e-6, 5.49e-1},
       {4.01e-1, 2.64e-6, 5.99e-1}, {3.65e-1, 2.28e-6, 6.35e-1},
       {3.38e-1, 2.02e-6, 6.62e-1}, {3.16e-1, 1.83e-6, 6.84e-1},
       {2.97e-1, 1.68e-6, 7.03e-1}, {2.82e-1, 1.56e-6, 7.18e-1},
       {2.68e-1, 1.46e-6, 7.32e-1}, {2.56e-1, 1.37e-6, 7.44e-1}},
      {{5.35e-1, 4.48e-6, 4.65e-1}, {4.50e-1, 3.23e-6, 5.50e-1},
       {4.00e-1, 2.63e-6, 6.00e-1}, {3.64e-1, 2.27e-6, 6.36e-1},
       {3.37e-1, 2.01e-6, 6.63e-1}, {3.15e-1, 1.82e-6, 6.85e-1},
       {2.96e-1, 1.67e-6, 7.04e-1}, {2.81e-1, 1.55e-6, 7.19e-1},
       {2.67e-1, 1.45e-6, 7.33e-1}, {2.55e-1, 1.37e-6, 7.45e-1}},
      {{5.36e-1, 4.49e-6, 4.64e-1}, {4.50e-1, 3.20e-6, 5.50e-1},
       {4.00e-1, 2.62e-6, 6.00e-1}, {3.64e-1, 2.26e-6, 6.36e-1},
       {3.37e-1, 2.01e-6, 6.63e-1}, {3.15e-1, 1.82e-6, 6.85e-1},
       {2.97e-1, 1.67e-6, 7.03e-1}, {2.81e-1, 1.55e-6, 7.19e-1},
       {2.67e-1, 1.45e-6, 7.33e-1}, {2.56e-1, 1.37e-6, 7.44e-1}},
      {{5.36e-1, 4.49e-6, 4.64e-1}, {4.51e-1, 3.23e-6, 5.49e-1},
       {4.00e-1, 2.63e-6, 6.00e-1}, {3.64e-1, 2.27e-6, 6.36e-1},
       {3.37e-1, 2.01e-6, 6.63e-1}, {3.15e-1, 1.82e-6, 6.85e-1},
       {2.97e-1, 1.68e-6, 7.03e-1}, {2.81e-1, 1.55e-6, 7.19e-1},
       {2.68e-1, 1.45e-6, 7.32e-1}, {2.56e-1, 1.37e-6, 7.44e-1}},
      {{5.36e-1, 4.49e-6, 4.64e-1}, {4.51e-1, 3.22e-6, 5.49e-1},
       {4.00e-1, 2.63e-6, 6.00e-1}, {3.64e-1, 2.27e-6, 6.36e-1},
       {3.37e-1, 2.01e-6, 6.63e-1}, {3.15e-1, 1.82e-6, 6.85e-1},
       {2.97e-1, 1.67e-6, 7.03e-1}, {2.81e-1, 1.55e-6, 7.19e-1},
       {2.67e-1, 1.45e-6, 7.33e-1}, {2.56e-1, 1.37e-6, 7.44e-1}}};
  const size_t scheme_counts[5][6] = {
      {149, 298, 149, 149, 298, 0}, {98, 196, 98, 98, 294, 0},
      {97, 291, 97, 97, 388, 0}, {99, 297, 99, 99, 396, 0},
      {97, 582, 97, 97, 582, 0}};
  const chem_solver_options scheme_options = chem_solver_default_options();
  for (size_t scheme = 0; scheme < 5; ++scheme) {
    const chem_rosenbrock_parameters parameters =
        chem_rosenbrock_builtin_parameters(tutorial_schemes[scheme]);
    double scheme_y[] = {1.0, 0.0, 0.0};
    chem_solver_workspace scheme_workspace = {0};
    chem_solver_stats scheme_total = {0};
    if (chem_solver_workspace_init_group(&scheme_workspace, &solve_mechanism, 4) != 0)
      return 46;
    int scheme_status = 0;
    for (size_t step = 0; step < 10 && !scheme_status; ++step) {
      chem_solver_stats one = {0};
      scheme_status = chem_rosenbrock_solve_with_parameters(
          &solve_mechanism, &scheme_workspace, scheme_y, 272.5, 101253.3,
          1e6, solve_rates, 200.0, &parameters, &scheme_options, &one);
      scheme_total.internal_steps += one.internal_steps;
      scheme_total.accepted += one.accepted;
      scheme_total.rejected += one.rejected;
      scheme_total.function_calls += one.function_calls;
      scheme_total.jacobian_updates += one.jacobian_updates;
      scheme_total.decompositions += one.decompositions;
      scheme_total.solves += one.solves;
      for (size_t variable = 0; variable < 3 && !scheme_status; ++variable) {
        const double expected = scheme_expected[scheme][step][variable];
        const double print_tolerance =
            0.5 * pow(10.0, floor(log10(expected)) - 2.0);
        if (fabs(scheme_y[variable] - expected) > print_tolerance)
          scheme_status = 1;
      }
    }
    chem_solver_workspace_destroy(&scheme_workspace);
    if (scheme_status || scheme_total.internal_steps != scheme_counts[scheme][0] ||
        scheme_total.function_calls != scheme_counts[scheme][1] ||
        scheme_total.jacobian_updates != scheme_counts[scheme][2] ||
        scheme_total.decompositions != scheme_counts[scheme][3] ||
        scheme_total.solves != scheme_counts[scheme][4] ||
        scheme_total.rejected != scheme_counts[scheme][5]) return 47;
  }
  if (check_solver_workspace_argument_errors() != 0) return 72;
  if (check_mechanism_species_usage() != 0) return 73;
  if (check_rate_constants_no_user_tutorial() != 0) return 48;
  if (check_lambda_rate_constant_position() != 0) return 49;
  if (check_jacobian_verification() != 0) return 50;
  if (check_version_api() != 0) return 51;
  if (check_one_solver_many_thread_states() != 0) return 52;
  if (check_state_view() != 0) return 53;
  if (check_process_rate_constant_layouts() != 0) return 54;
  if (check_process_set_finite_difference_jacobian() != 0) return 55;
  if (check_random_process_sets() != 0) return 56;
  if (check_process_set_forcing_and_jacobian() != 0) return 57;
  const int backward_euler_status = check_backward_euler();
  if (backward_euler_status != 0) return 70 + backward_euler_status;
  for (int variant=0;variant<5;++variant)
    if (check_superstiff_mechanism(variant,0) ||
        check_superstiff_mechanism(variant,1)) return 134+variant;
  if (check_terminator_cases() != 0) return 58;
  { const int status=check_chapman_integration();
    if (status) return 59+status; }
  if (check_integrated_reaction_rates() != 0) return 60;
  if (check_flow_tube_reference(0) != 0) return 74;
  if (check_flow_tube_reference(1) != 0) return 75;
  { const int status=check_stub_aerosol_mechanism(0); if (status) return 75+status; }
  { const int status=check_stub_aerosol_mechanism(1); if (status) return 80+status; }
  if (check_lambda_callback_override() != 0) return 61;
  if (check_phase_species_metadata() != 0) return 62;
  if (check_constraint_api() != 0) return 63;
  if (check_constrained_solver() != 0) return 65;
  for (int mode=0;mode<4;++mode)
    if (check_equilibrium_dae_case(mode)!=0) return 120+mode;
  if (check_robertson_terminator_dae()!=0) return 130;
  if (check_all_constrained_overshoot_orders()!=0) return 131;
  if (check_equilibrium_conservation_overshoot()!=0) return 132;
  if (check_process_algebraic_mask() != 0) return 66;
  { const int status=check_external_constraint_equivalents(); if (status) return 135+status; }
  { const int status=check_external_equilibrium_convergence(); if (status) return 150+status; }
  { const int status=check_composed_equilibrium_dae(); if (status) return 160+status; }
  { const int status=check_temperature_dependent_equilibrium(); if (status) return 170+status; }
  { const int status=check_process_jacobian_in_constraint_rows(); if (status) return 180+status; }
  { const int status=check_constraint_workspace_reuse(); if (status) return 190+status; }
  if (check_phase_system_metadata() != 0) return 67;
  if (check_species_properties() != 0) return 68;
  if (check_matrix_view_api() != 0) return 64;
  if (check_csr_matrix_api() != 0) return 69;
  if (check_lu_decomposition() != 0) return 70;
  if (check_linear_solvers() != 0) return 71;
  return 0;
}
