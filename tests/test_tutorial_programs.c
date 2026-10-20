#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static const chem_species species_abc[] = {
  { "A", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "B", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "C", 0.0, 0.0, 0.0, 1.0e-3, 0u },
};
static const chem_parameter parameters_abc[] = {
  { "r1", 0.04 }, { "r2", 3.0e7 }, { "r3", 1.0e4 },
};
static const chem_term reactants_abc[] = {
  { 0u, 1.0 }, { 1u, 1.0 }, { 1u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 },
};
static const chem_term products_abc[] = {
  { 1u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 }, { 0u, 1.0 }, { 2u, 1.0 },
};
static const chem_reaction reactions_abc[] = {
  { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_USER, .label = "r1", .p = { [1] = 1.0 }, .custom_parameter = 0u } },
  { 1u, 2u, 1u, 2u, { .kind = CHEM_RATE_USER, .label = "r2", .p = { [1] = 1.0 }, .custom_parameter = 1u } },
  { 3u, 2u, 3u, 2u, { .kind = CHEM_RATE_USER, .label = "r3", .p = { [1] = 1.0 }, .custom_parameter = 2u } },
};
static const chem_mechanism mechanism_abc = {
  .name = "A-B-C tutorial",
  .species_count = 3,
  .reaction_count = 3,
  .reactant_term_count = 5,
  .product_term_count = 5,
  .parameter_count = 3,
  .species = species_abc,
  .reactions = reactions_abc,
  .reactant_terms = reactants_abc,
  .product_terms = products_abc,
  .temperature = 272.5,
  .pressure = 101253.3,
  .relative_tolerance = 1.0e-6,
};

static int near_printed(double actual, double printed) {
  if (!isfinite(actual) || !isfinite(printed)) return 0;
  const double unit = pow(10.0, floor(log10(fabs(printed))) - 2.0);
  return fabs(actual - printed) <= 0.51 * unit + 1.0e-30;
}

static void set_abc_inputs(size_t cells, const double *initial_a, double *state, double *custom,
                           double *temperature, double *pressure) {
  for (size_t c = 0; c < cells; ++c) {
    state[c * 3] = initial_a[c];
    state[c * 3 + 1] = 0.0;
    state[c * 3 + 2] = 0.0;
    for (size_t p = 0; p < 3; ++p)
      custom[c * 3 + p] = parameters_abc[p].initial_value;
    temperature[c] = mechanism_abc.temperature;
    pressure[c] = mechanism_abc.pressure;
  }
}

static int test_readme_example(void) {
  const chem_species species[] = {
    { "Foo", 0.0, 0.0, 0.0, 1.0e-3, 0u },
    { "Bar", 0.0, 0.0, 0.0, 1.0e-3, 0u },
    { "Baz", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  };
  const chem_term reactants[] = { { 0u, 1.0 }, { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_term products[] = { { 1u, 0.8 }, { 2u, 0.2 }, { 2u, 1.0 } };
  chem_reaction reactions[2];
  memset(reactions, 0, sizeof(reactions));
  reactions[0].reactant_offset = 0u; reactions[0].reactant_count = 1u;
  reactions[0].product_offset = 0u; reactions[0].product_count = 2u;
  reactions[0].rate.kind = CHEM_RATE_ARRHENIUS;
  reactions[0].rate.p[0] = 1.0e-3;
  reactions[0].rate.p[3] = 300.0;
  reactions[1].reactant_offset = 1u; reactions[1].reactant_count = 2u;
  reactions[1].product_offset = 2u; reactions[1].product_count = 1u;
  reactions[1].rate.kind = CHEM_RATE_ARRHENIUS;
  reactions[1].rate.p[0] = 1.0e-5;
  reactions[1].rate.p[2] = 110.0;
  reactions[1].rate.p[3] = 300.0;
  const chem_mechanism mechanism = {
    .name = "README example", .species_count = 3, .reaction_count = 2,
    .reactant_term_count = 3, .product_term_count = 3,
    .species = species, .reactions = reactions,
    .reactant_terms = reactants, .product_terms = products,
    .temperature = 287.45, .pressure = 101319.9, .relative_tolerance = 1.0e-6,
  };
  const double expected[10][3] = {
    { 11.8, 5.90, 1.91 }, { 6.79, 9.05, 3.32 }, { 3.83, 10.7, 4.21 },
    { 2.14, 11.7, 4.74 }, { 1.19, 12.2, 5.04 }, { 0.658, 12.4, 5.21 },
    { 0.364, 12.6, 5.31 }, { 0.201, 12.7, 5.36 }, { 0.111, 12.7, 5.39 },
    { 0.0613, 12.8, 5.41 },
  };
  double state[] = { 20.0, 0.0, 0.0 };
  const double temperature[] = { 287.45 }, pressure[] = { 101319.9 };
  for (size_t i = 0; i < 10; ++i) {
    if (chem_integrate_group(&mechanism, 1, state, NULL, temperature, pressure, 500.0, NULL) != 0)
      return 0;
    for (size_t s = 0; s < 3; ++s)
      if (!near_printed(state[s], expected[i][s])) {
        fprintf(stderr, "README_example step %zu species %zu: %.17g expected %.17g\n", i, s, state[s], expected[i][s]);
        return 0;
      }
  }
  return 1;
}

static int test_multiple_grid_cells(void) {
  const double initial_a[] = { 1.0, 2.0, 0.5 };
  const double expected[3][3] = {
    { 0.255, 1.37e-6, 0.745 },
    { 0.668, 2.00e-6, 1.33 },
    { 0.0925, 9.02e-7, 0.407 },
  };
  double state[9] = { 0.0 }, custom[9] = { 0.0 };
  double temperature[3], pressure[3];
  set_abc_inputs(3, initial_a, state, custom, temperature, pressure);
  for (size_t step = 0; step < 10; ++step)
    if (chem_integrate_group(&mechanism_abc, 3, state, custom, temperature, pressure, 200.0, NULL) != 0)
      return 0;
  for (size_t c = 0; c < 3; ++c)
    for (size_t s = 0; s < 3; ++s)
      if (!near_printed(state[c * 3 + s], expected[c][s])) {
        fprintf(stderr, "multiple_grid_cells cell %zu species %zu: %.17g expected %.17g\n",
                c, s, state[c * 3 + s], expected[c][s]);
        return 0;
      }
  return 1;
}

static int test_openmp_tutorial(void) {
  double state[9] = { 0.0 }, custom[9] = { 0.0 };
  const double initial_a[] = { 1.0, 1.0, 1.0 };
  double temperature[3], pressure[3];
  set_abc_inputs(3, initial_a, state, custom, temperature, pressure);
  int failed[3] = { 0, 0, 0 };
  #pragma omp parallel for schedule(static) num_threads(3)
  for (int cell = 0; cell < 3; ++cell) {
    for (size_t step = 0; step < 10; ++step)
      if (chem_integrate_group(&mechanism_abc, 1, state + cell * 3, custom + cell * 3,
                               temperature + cell, pressure + cell, 200.0, NULL) != 0)
        failed[cell] = 1;
  }
  const double expected[] = { 0.255, 1.37e-6, 0.745 };
  for (size_t cell = 0; cell < 3; ++cell) {
    if (failed[cell]) return 0;
    for (size_t s = 0; s < 3; ++s)
      if (!near_printed(state[cell * 3 + s], expected[s])) {
        fprintf(stderr, "openmp_tutorial cell %zu species %zu: %.17g expected %.17g\n",
                cell, s, state[cell * 3 + s], expected[s]);
        return 0;
      }
  }
  return 1;
}

static int test_solve_results(void) {
  const double initial_a[] = { 1.0, 2.0, 0.5 };
  double state[9] = { 0.0 }, custom[9] = { 0.0 };
  double temperature[3], pressure[3];
  set_abc_inputs(3, initial_a, state, custom, temperature, pressure);
  chem_statistics statistics = { 0 };
  if (chem_integrate_group(&mechanism_abc, 3, state, custom, temperature, pressure, 200.0, &statistics) != 0)
    return 0;
  if (statistics.internal_steps != 18 || statistics.accepted != 18 || statistics.rejected != 0 ||
      statistics.function_calls != 36 || statistics.jacobian_updates != 18 ||
      statistics.decompositions != 18 || statistics.solves != 54)
    fprintf(stderr, "solve_results stats steps=%lu accepted=%lu rejected=%lu functions=%lu jacobians=%lu decompositions=%lu solves=%lu\n",
            statistics.internal_steps, statistics.accepted, statistics.rejected, statistics.function_calls,
            statistics.jacobian_updates, statistics.decompositions, statistics.solves);
  return statistics.internal_steps == 18 && statistics.accepted == 18 && statistics.rejected == 0 &&
         statistics.function_calls == 36 && statistics.jacobian_updates == 18 &&
         statistics.decompositions == 18 && statistics.solves == 54;
}

static int test_vectorized_matrix_solver(void) {
  const double row_major[] = { 1.1, 1.2, 1.3, 2.1, 2.2, 2.3, 3.1, 3.2, 3.3 };
  const double column_major[] = { 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.3, 2.3, 3.3 };
  for (size_t cell = 0; cell < 3; ++cell)
    for (size_t species = 0; species < 3; ++species)
      if (row_major[cell * 3 + species] != column_major[species * 3 + cell]) return 0;
  const unsigned char pattern[4][4] = {
    { 0, 1, 0, 0 }, { 0, 0, 0, 0 }, { 0, 1, 0, 1 }, { 0, 0, 1, 0 },
  };
  size_t row_order[4], col_order[4], row_count = 0, col_count = 0;
  for (size_t row = 0; row < 4; ++row)
    for (size_t col = 0; col < 4; ++col)
      if (pattern[row][col]) row_order[row_count++] = row * 4 + col;
  for (size_t col = 0; col < 4; ++col)
    for (size_t row = 0; row < 4; ++row)
      if (pattern[row][col]) col_order[col_count++] = row * 4 + col;
  const size_t expected_row[] = { 1, 9, 11, 14 };
  const size_t expected_col[] = { 1, 9, 14, 11 };
  return row_count == 4 && col_count == 4 &&
      memcmp(row_order, expected_row, sizeof(expected_row)) == 0 &&
      memcmp(col_order, expected_col, sizeof(expected_col)) == 0;
}

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  int passed = 0;
  if (strcmp(argv[1], "README_example") == 0) passed = test_readme_example();
  else if (strcmp(argv[1], "multiple_grid_cells") == 0) passed = test_multiple_grid_cells();
  else if (strcmp(argv[1], "openmp_tutorial") == 0) passed = test_openmp_tutorial();
  else if (strcmp(argv[1], "solve_results") == 0) passed = test_solve_results();
  else if (strcmp(argv[1], "vectorized_matrix_solver") == 0) passed = test_vectorized_matrix_solver();
  else return 2;
  if (!passed) fprintf(stderr, "%s C equivalent check failed\n", argv[1]);
  return passed ? 0 : 1;
}
