#include "chemistry.h"

#include <math.h>
#include <omp.h>
#include <string.h>

static double test_lambda(double temperature, double pressure,
                          double air_density, const double *parameters) {
  (void)parameters;
  (void)pressure;
  (void)air_density;
  return 2.0 * temperature;
}

static int near(double actual, double expected) {
  return fabs(actual - expected) <= 1e-13 * fabs(expected) + 1e-30;
}

static int tutorial_print_match(double actual, double expected) {
  if (expected == 0.0) return actual == 0.0;
  const double tolerance = 0.5 * pow(10.0, floor(log10(expected)) - 2.0);
  return fabs(actual - expected) <= tolerance;
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

int main(void) {
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
  chem_solver_options options = {1.0e-10, 1.0e-12, 1000};
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
  const chem_solver_options chain_options = {1.0e-10, 1.0e-12, 10000};
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
      solve_stats.decompositions != 18 || solve_stats.solves != 54) return 40;

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
    if (chem_solver_workspace_init(&scheme_workspace, &solve_mechanism) != 0)
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
  if (check_rate_constants_no_user_tutorial() != 0) return 48;
  if (check_lambda_rate_constant_position() != 0) return 49;
  if (check_jacobian_verification() != 0) return 50;
  if (check_version_api() != 0) return 51;
  if (check_one_solver_many_thread_states() != 0) return 52;
  if (check_state_view() != 0) return 53;
  return 0;
}
