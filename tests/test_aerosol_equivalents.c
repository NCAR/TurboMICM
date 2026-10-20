#include "chem.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define COUNT(items) (sizeof(items) / sizeof((items)[0]))
#define SPECIES_COUNT 22u
#define REACTION_COUNT 4u
#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "%s:%d check failed: %s\n", __FILE__, __LINE__, #condition); return 1; \
} } while (0)

enum {
  FO2 = 0, BAR = 1, STUB1_MODE1_QUUX_BAZ = 2, STUB1_MODE1_QUUX_QUX = 3,
  STUB1_MODE2_QUUX_BAZ = 4, STUB1_MODE2_QUUX_QUX = 5,
  STUB1_MODE2_CORGE_FO2 = 6, STUB1_MODE2_CORGE_BAZ = 7, STUB1_MODE2_CORGE_QUX = 8,
  STUB2_MODE1_NUMBER = 9, STUB2_MODE1_QUUX_BAZ = 10, STUB2_MODE1_QUUX_QUX = 11,
  STUB2_MODE2_NUMBER = 12, STUB2_MODE2_CORGE_FO2 = 13, STUB2_MODE2_CORGE_BAZ = 14,
  STUB2_MODE2_CORGE_QUX = 15, STUB2_MODE3_NUMBER = 16, STUB2_MODE3_QUUX_BAZ = 17,
  STUB2_MODE3_QUUX_QUX = 18, STUB2_MODE3_CORGE_FO2 = 19, STUB2_MODE3_CORGE_BAZ = 20,
  STUB2_MODE3_CORGE_QUX = 21,
};

static const char *const state_names[SPECIES_COUNT] = {
  "FO2", "BAR", "STUB1.MODE1.QUUX.BAZ", "STUB1.MODE1.QUUX.QUX",
  "STUB1.MODE2.QUUX.BAZ", "STUB1.MODE2.QUUX.QUX", "STUB1.MODE2.CORGE.FO2",
  "STUB1.MODE2.CORGE.BAZ", "STUB1.MODE2.CORGE.QUX", "STUB2.MODE1.NUMBER",
  "STUB2.MODE1.QUUX.BAZ", "STUB2.MODE1.QUUX.QUX", "STUB2.MODE2.NUMBER",
  "STUB2.MODE2.CORGE.FO2", "STUB2.MODE2.CORGE.BAZ", "STUB2.MODE2.CORGE.QUX",
  "STUB2.MODE3.NUMBER", "STUB2.MODE3.QUUX.BAZ", "STUB2.MODE3.QUUX.QUX",
  "STUB2.MODE3.CORGE.FO2", "STUB2.MODE3.CORGE.BAZ", "STUB2.MODE3.CORGE.QUX",
};
static const chem_species species[SPECIES_COUNT] = {
  { "FO2", 0.0, 0.0, 0.0, 1e-10, 0u }, { "BAR", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB1.MODE1.QUUX.BAZ", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB1.MODE1.QUUX.QUX", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB1.MODE2.QUUX.BAZ", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB1.MODE2.QUUX.QUX", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB1.MODE2.CORGE.FO2", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB1.MODE2.CORGE.BAZ", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB1.MODE2.CORGE.QUX", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE1.NUMBER", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE1.QUUX.BAZ", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE1.QUUX.QUX", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE2.NUMBER", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE2.CORGE.FO2", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE2.CORGE.BAZ", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE2.CORGE.QUX", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE3.NUMBER", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE3.QUUX.BAZ", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE3.QUUX.QUX", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE3.CORGE.FO2", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE3.CORGE.BAZ", 0.0, 0.0, 0.0, 1e-10, 0u },
  { "STUB2.MODE3.CORGE.QUX", 0.0, 0.0, 0.0, 1e-10, 0u },
};
static const chem_parameter parameters[] = {
  { "STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT", 0.01 },
};
static const chem_term reactants[] = {
  { FO2, 1.0 }, { STUB1_MODE1_QUUX_BAZ, 1.0 },
  { STUB2_MODE2_CORGE_FO2, 1.0 }, { STUB2_MODE3_QUUX_BAZ, 1.0 },
};
static const chem_term products[] = {
  { STUB1_MODE2_CORGE_FO2, 1.0 }, { STUB1_MODE2_QUUX_BAZ, 1.0 },
  { STUB2_MODE2_CORGE_BAZ, 1.0 }, { STUB2_MODE3_QUUX_QUX, 1.0 },
};

static double temperature_rate(double temperature, double pressure, double air_density, void *context) {
  (void)pressure;
  (void)air_density;
  (void)context;
  return 0.005 * temperature;
}

static const chem_reaction reactions[REACTION_COUNT] = {
  { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 1.0e-3, [3] = 298.15 } } },
  { 1u, 1u, 1u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 2.0e-3, [3] = 298.15 } } },
  { 2u, 1u, 2u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 0u } },
  { 3u, 1u, 3u, 1u, { .kind = CHEM_RATE_CALLBACK, .callback = temperature_rate } },
};
static const chem_mechanism mechanism = {
  .name = "external aerosol process equivalents", .species_count = SPECIES_COUNT,
  .reaction_count = REACTION_COUNT, .reactant_term_count = COUNT(reactants),
  .product_term_count = COUNT(products), .parameter_count = COUNT(parameters),
  .species = species, .reactions = reactions, .reactant_terms = reactants,
  .product_terms = products, .parameters = parameters, .temperature = 298.15,
  .pressure = 101325.0, .relative_tolerance = 1e-9,
};

static size_t find_species(const char *name) {
  for (size_t i = 0; i < SPECIES_COUNT; ++i)
    if (strcmp(state_names[i], name) == 0) return i;
  return SIZE_MAX;
}

static int state_includes(void) {
  CHECK(SPECIES_COUNT == 2u + 7u + 13u);
  CHECK(COUNT(parameters) == 1u);
  for (size_t i = 0; i < SPECIES_COUNT; ++i) {
    CHECK(strcmp(species[i].name, state_names[i]) == 0);
    for (size_t j = i + 1u; j < SPECIES_COUNT; ++j) CHECK(strcmp(species[i].name, species[j].name) != 0);
  }
  CHECK(find_species("FO2") == FO2);
  CHECK(find_species("STUB1.MODE1.QUUX.BAZ") == STUB1_MODE1_QUUX_BAZ);
  CHECK(find_species("STUB1.MODE2.CORGE.FO2") == STUB1_MODE2_CORGE_FO2);
  CHECK(find_species("STUB2.MODE1.NUMBER") == STUB2_MODE1_NUMBER);
  CHECK(find_species("STUB2.MODE2.CORGE.BAZ") == STUB2_MODE2_CORGE_BAZ);
  CHECK(find_species("STUB2.MODE3.CORGE.QUX") == STUB2_MODE3_CORGE_QUX);
  CHECK(strcmp(parameters[0].name, "STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT") == 0);
  return 0;
}

static int state_updates(int cells) {
  double values[3u * SPECIES_COUNT] = { 0.0 };
  const double foo[] = { 1.34, 1.35, 1.36 }, bar[] = { 2.53, 2.54, 2.55 };
  const double mode1_baz[] = { 0.75, 0.76, 0.77 }, mode2_fo2[] = { 1.23, 1.24, 1.25 };
  const double mode3_qux[] = { 0.42, 0.43, 0.44 }, mode3_baz[] = { 0.33, 0.34, 0.35 };
  const double mode1_number[] = { 1000.0, 1001.0, 1002.0 }, mode2_number[] = { 500.0, 501.0, 502.0 };
  for (size_t cell = 0; cell < (size_t)cells; ++cell) {
    double *row = values + cell * SPECIES_COUNT;
    row[FO2] = foo[cell]; row[BAR] = bar[cell]; row[STUB1_MODE1_QUUX_BAZ] = mode1_baz[cell];
    row[STUB1_MODE2_CORGE_FO2] = mode2_fo2[cell]; row[STUB2_MODE3_CORGE_QUX] = mode3_qux[cell];
    row[STUB2_MODE3_CORGE_BAZ] = mode3_baz[cell]; row[STUB2_MODE1_NUMBER] = mode1_number[cell];
    row[STUB2_MODE2_NUMBER] = mode2_number[cell];
  }
  for (size_t cell = 0; cell < (size_t)cells; ++cell) {
    const double *row = values + cell * SPECIES_COUNT;
    CHECK(row[FO2] == foo[cell] && row[BAR] == bar[cell]);
    CHECK(row[STUB1_MODE1_QUUX_BAZ] == mode1_baz[cell]);
    CHECK(row[STUB1_MODE2_CORGE_FO2] == mode2_fo2[cell]);
    CHECK(row[STUB2_MODE3_CORGE_QUX] == mode3_qux[cell]);
    CHECK(row[STUB2_MODE3_CORGE_BAZ] == mode3_baz[cell]);
    CHECK(row[STUB2_MODE1_NUMBER] == mode1_number[cell]);
    CHECK(row[STUB2_MODE2_NUMBER] == mode2_number[cell]);
  }
  return 0;
}

static int direct_forcing(void) {
  const double state[SPECIES_COUNT] = {
    [FO2] = 1.0, [STUB1_MODE1_QUUX_BAZ] = 0.5,
  };
  const double rates[] = { 1.0e-3, 2.0e-3, 1.0e-2, 1.0 };
  double forcing[SPECIES_COUNT];
  CHECK(chem_evaluate_forcing(&mechanism, 1u, state, rates, forcing) == 0);
  CHECK(forcing[FO2] == -1.0e-3 && forcing[STUB1_MODE2_CORGE_FO2] == 1.0e-3);
  CHECK(forcing[STUB1_MODE1_QUUX_BAZ] == -1.0e-3 && forcing[STUB1_MODE2_QUUX_BAZ] == 1.0e-3);
  for (size_t species_index = 0; species_index < SPECIES_COUNT; ++species_index)
    if (species_index != FO2 && species_index != STUB1_MODE2_CORGE_FO2 &&
        species_index != STUB1_MODE1_QUUX_BAZ && species_index != STUB1_MODE2_QUUX_BAZ)
      CHECK(forcing[species_index] == 0.0);
  return 0;
}

static int direct_jacobian(void) {
  const double state[SPECIES_COUNT] = { [FO2] = 1.0, [STUB1_MODE1_QUUX_BAZ] = 0.5 };
  const double rates[] = { 1.0e-3, 2.0e-3, 1.0e-2, 1.0 };
  double negative_jacobian[SPECIES_COUNT * SPECIES_COUNT] = { 0.0 };
  CHECK(chem_accumulate_negative_jacobian(&mechanism, 1u, state, rates, NULL, negative_jacobian) == 0);
  CHECK(negative_jacobian[FO2 * SPECIES_COUNT + FO2] == 1.0e-3);
  CHECK(negative_jacobian[STUB1_MODE2_CORGE_FO2 * SPECIES_COUNT + FO2] == -1.0e-3);
  CHECK(negative_jacobian[STUB1_MODE1_QUUX_BAZ * SPECIES_COUNT + STUB1_MODE1_QUUX_BAZ] == 2.0e-3);
  CHECK(negative_jacobian[STUB1_MODE2_QUUX_BAZ * SPECIES_COUNT + STUB1_MODE1_QUUX_BAZ] == -2.0e-3);
  return 0;
}

static chem_rosenbrock_parameters make_solver(size_t method, chem_lu_method lu_method) {
  chem_rosenbrock_parameters solver;
  switch (method) {
    case 0u: solver = chem_rosenbrock_two_stage_parameters(); break;
    case 1u: solver = chem_rosenbrock_default_parameters(); break;
    case 2u: solver = chem_rosenbrock_four_stage_parameters(); break;
    case 3u: solver = chem_rosenbrock_four_stage_da_parameters(); break;
    default: solver = chem_rosenbrock_six_stage_da_parameters(); break;
  }
  solver.lu_method = lu_method;
  solver.max_internal_steps = 10000u;
  return solver;
}

static int integrate(const char *method, size_t method_index, chem_lu_method lu_method,
    size_t cells, double *state, double *rate_parameters, const double *temperature) {
  const double pressure[3] = { 101325.0, 101325.0, 101325.0 };
  if (strcmp(method, "backward_euler") == 0)
    return chem_integrate_group_backward_euler(&mechanism, cells, state, rate_parameters,
        temperature, pressure, 10.0, NULL);
  chem_rosenbrock_parameters solver = make_solver(method_index, lu_method);
  return chem_integrate_group_with_parameters(&mechanism, cells, state, rate_parameters,
      temperature, pressure, 10.0, &solver, NULL);
}

static int check_value(const char *name, size_t cell, size_t species_index,
    double actual, double expected, double tolerance) {
  if (isfinite(actual) && fabs(actual - expected) <= tolerance) return 0;
  fprintf(stderr, "%s cell %zu %s: %.17g vs %.17g tolerance %.3g\n",
      name, cell, species[species_index].name, actual, expected, tolerance);
  return 1;
}

static void initialize_first_model(double *state, size_t cells) {
  memset(state, 0, cells * SPECIES_COUNT * sizeof(*state));
  for (size_t cell = 0; cell < cells; ++cell) {
    double *row = state + cell * SPECIES_COUNT;
    row[FO2] = cells == 1u ? 1.0 : 1.0 + 0.5 * (double)cell;
    row[BAR] = 2.0 + 0.5 * (double)cell;
    row[STUB1_MODE1_QUUX_BAZ] = cells == 1u ? 0.5 : 0.5 + 0.2 * (double)cell;
    row[STUB1_MODE2_CORGE_FO2] = cells == 1u ? 0.8 : 0.8 + 0.2 * (double)cell;
    row[STUB2_MODE3_CORGE_QUX] = 0.3 + 0.1 * (double)cell;
    row[STUB2_MODE3_CORGE_BAZ] = 0.2 + 0.1 * (double)cell;
    row[STUB2_MODE1_NUMBER] = 1000.0 + 100.0 * (double)cell;
    row[STUB2_MODE2_NUMBER] = 500.0 + 50.0 * (double)cell;
  }
}

static void initialize_two_models(double *state, size_t cells) {
  initialize_first_model(state, cells);
  for (size_t cell = 0; cell < cells; ++cell) {
    double *row = state + cell * SPECIES_COUNT;
    row[STUB2_MODE2_CORGE_FO2] = 0.4 + 0.2 * (double)cell;
    row[STUB2_MODE2_CORGE_BAZ] = 0.1 + 0.05 * (double)cell;
    row[STUB2_MODE3_QUUX_BAZ] = 0.2 + 0.1 * (double)cell;
    row[STUB2_MODE3_QUUX_QUX] = 0.3 + 0.1 * (double)cell;
  }
}

static int solve_first_model(const char *method, size_t cells, size_t method_index, chem_lu_method lu_method) {
  double state[3u * SPECIES_COUNT];
  double rate_parameters[] = { 0.01, 0.01, 0.01 };
  const double temperature[] = { 275.0, 285.0, 295.0 };
  initialize_first_model(state, cells);
  double initial_fo2[3], initial_baz[3];
  for (size_t cell = 0; cell < cells; ++cell) {
    initial_fo2[cell] = state[cell * SPECIES_COUNT + FO2];
    initial_baz[cell] = state[cell * SPECIES_COUNT + STUB1_MODE1_QUUX_BAZ];
  }
  if (integrate(method, method_index, lu_method, cells, state, rate_parameters, temperature) != 0) return 1;
  const double delta_fo2_fraction = 1.0 - exp(-1.0e-3 * 10.0);
  const double delta_baz_fraction = 1.0 - exp(-2.0e-3 * 10.0);
  const double tolerance = 2.0e-4;
  for (size_t cell = 0; cell < cells; ++cell) {
    const double *row = state + cell * SPECIES_COUNT;
    const double fo2_delta = initial_fo2[cell] * delta_fo2_fraction;
    const double baz_delta = initial_baz[cell] * delta_baz_fraction;
    if (check_value("aerosol1 FO2", cell, FO2, row[FO2], initial_fo2[cell] - fo2_delta, tolerance * initial_fo2[cell]) ||
        check_value("aerosol1 FO2 transfer", cell, STUB1_MODE2_CORGE_FO2, row[STUB1_MODE2_CORGE_FO2],
          (cells == 1u ? 0.8 : 0.8 + 0.2 * (double)cell) + fo2_delta, tolerance * initial_fo2[cell]) ||
        check_value("aerosol1 BAZ loss", cell, STUB1_MODE1_QUUX_BAZ, row[STUB1_MODE1_QUUX_BAZ],
          initial_baz[cell] - baz_delta, tolerance * initial_baz[cell]) ||
        check_value("aerosol1 BAZ gain", cell, STUB1_MODE2_QUUX_BAZ, row[STUB1_MODE2_QUUX_BAZ],
          baz_delta, tolerance * initial_baz[cell])) return 1;
  }
  return 0;
}

static int solve_two_models(const char *method, size_t cells, size_t method_index, chem_lu_method lu_method) {
  double state[3u * SPECIES_COUNT];
  double parameters_by_cell[3u] = { 0.01, 0.015, 0.02 };
  const double temperature[] = { 275.0, 285.0, 295.0 };
  initialize_two_models(state, cells);
  double initial[3u][4u];
  for (size_t cell = 0; cell < cells; ++cell) {
    const double *row = state + cell * SPECIES_COUNT;
    initial[cell][0] = row[FO2];
    initial[cell][1] = row[STUB1_MODE1_QUUX_BAZ];
    initial[cell][2] = row[STUB2_MODE2_CORGE_FO2];
    initial[cell][3] = row[STUB2_MODE3_QUUX_BAZ];
  }
  if (integrate(method, method_index, lu_method, cells, state, parameters_by_cell, temperature) != 0) return 1;
  const double tolerance = strcmp(method, "backward_euler") == 0 ? 1.0e-1 :
      (method_index == 0u ? 5.0e-4 : 5.0e-5);
  for (size_t cell = 0; cell < cells; ++cell) {
    const double *row = state + cell * SPECIES_COUNT;
    const double delta_fo2 = initial[cell][0] * (1.0 - exp(-1.0e-3 * 10.0));
    const double delta_baz = initial[cell][1] * (1.0 - exp(-2.0e-3 * 10.0));
    const double delta_custom = initial[cell][2] * (1.0 - exp(-parameters_by_cell[cell] * 10.0));
    const double delta_temperature = initial[cell][3] * (1.0 - exp(-temperature[cell] * 0.005 * 10.0));
    if (check_value("two-model FO2", cell, FO2, row[FO2], initial[cell][0] - delta_fo2, tolerance * initial[cell][0]) ||
        check_value("two-model BAZ mode1", cell, STUB1_MODE1_QUUX_BAZ, row[STUB1_MODE1_QUUX_BAZ],
          initial[cell][1] - delta_baz, tolerance * initial[cell][1]) ||
        check_value("two-model BAZ mode2", cell, STUB1_MODE2_QUUX_BAZ, row[STUB1_MODE2_QUUX_BAZ],
          delta_baz, tolerance * initial[cell][1]) ||
        check_value("two-model FO2 transfer", cell, STUB1_MODE2_CORGE_FO2, row[STUB1_MODE2_CORGE_FO2],
          0.8 + 0.2 * (double)cell + delta_fo2, tolerance * initial[cell][0]) ||
        check_value("two-model custom FO2", cell, STUB2_MODE2_CORGE_FO2, row[STUB2_MODE2_CORGE_FO2],
          initial[cell][2] - delta_custom, tolerance * initial[cell][2]) ||
        check_value("two-model custom BAZ", cell, STUB2_MODE2_CORGE_BAZ, row[STUB2_MODE2_CORGE_BAZ],
          0.1 + 0.05 * (double)cell + delta_custom, tolerance * initial[cell][2]) ||
        check_value("two-model thermal BAZ", cell, STUB2_MODE3_QUUX_BAZ, row[STUB2_MODE3_QUUX_BAZ],
          initial[cell][3] - delta_temperature, tolerance * initial[cell][3]) ||
        check_value("two-model thermal QUX", cell, STUB2_MODE3_QUUX_QUX, row[STUB2_MODE3_QUUX_QUX],
          0.3 + 0.1 * (double)cell + delta_temperature, tolerance * initial[cell][3])) return 1;
  }
  return 0;
}

static int run_suite_case(const char *method, const char *test_name) {
  if (strcmp(test_name, "StateIncludesStubAerosolModel") == 0) return state_includes();
  if (strcmp(test_name, "CanUpdateStateWithStubAerosolModel") == 0) return state_updates(1);
  if (strcmp(test_name, "CanUpdateMultiCellStateWithStubAerosolModel") == 0) return state_updates(3);
  if (strcmp(test_name, "CanCalculateSingleGridCellForcingWithStubAerosolModel") == 0) return direct_forcing();
  if (strcmp(test_name, "CanCalculateSingleGridCellJacobianWithStubAerosolModel") == 0) return direct_jacobian();
  if (strcmp(test_name, "CanSolveSingleGridCellWithStubAerosolModel1") == 0)
    return solve_first_model(method, 1u, 1u, CHEM_LU_DOOLITTLE);
  if (strcmp(test_name, "CanSolveSingleGridCellWithTwoStubAerosolModels") == 0)
    return solve_two_models(method, 1u, 1u, CHEM_LU_DOOLITTLE);
  if (strcmp(test_name, "CanSolveMultiGridCellWithStubAerosolModel1") == 0)
    return solve_first_model(method, 3u, 1u, CHEM_LU_DOOLITTLE);
  if (strcmp(test_name, "CanSolveMultiGridCellWithTwoStubAerosolModels") == 0)
    return solve_two_models(method, 3u, 1u, CHEM_LU_DOOLITTLE);
  return 2;
}

int main(int argc, char **argv) {
  if (argc != 3) return 2;
  return run_suite_case(argv[1], argv[2]);
}
