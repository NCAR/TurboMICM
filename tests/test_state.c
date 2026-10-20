#include "chem.h"

#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
  return 1; \
} } while (0)

static const char *const species[] = { "foo", "bar", "baz", "quz" };
static const char *const parameters[] = { "quux", "corge" };
static const char *const five_parameters[] = { "O1", "O2", "O3", "AAA", "BBB" };

static int init_state(chem_state *state, size_t cells, const char *const *parameter_names,
                      size_t parameter_count) {
  const double tolerance[] = { 1e-10, 1e-10, 1e-10, 1e-10 };
  return chem_state_init(state, species, 4u, 10u, parameter_names, parameter_count,
                         cells, 1e-5, tolerance);
}

static int default_constructor(void) {
  chem_state state = { 0 };
  chem_state states[10] = { { 0 } };
  chem_state_destroy(&state);
  for (size_t i = 0; i < 10u; ++i) chem_state_destroy(&states[i]);
  return 0;
}

static int constructor(void) {
  chem_state state = { 0 };
  CHECK(init_state(&state, 3u, parameters, 2u) == 0);
  CHECK(state.cell_count == 3u && state.species_count == 4u);
  CHECK(state.parameter_count == 2u && state.rate_constant_count == 10u);
  CHECK(state.concentrations != NULL && state.custom_parameters != NULL && state.rate_constants != NULL);
  CHECK(state.relative_tolerance == 1e-5 && state.absolute_tolerances[3] == 1e-10);
  CHECK(strcmp(state.species_names[0], "foo") == 0 && strcmp(state.species_names[3], "quz") == 0);
  CHECK(strcmp(state.parameter_names[0], "quux") == 0 && strcmp(state.parameter_names[1], "corge") == 0);
  chem_state_destroy(&state);
  return 0;
}

static int copy_state(void) {
  chem_state original = { 0 }, copy = { 0 };
  CHECK(init_state(&original, 3u, parameters, 2u) == 0);
  for (size_t i = 0; i < 12u; ++i) original.concentrations[i] = (double)i + 0.25;
  for (size_t i = 0; i < 6u; ++i) original.custom_parameters[i] = (double)i + 10.0;
  for (size_t i = 0; i < 30u; ++i) original.rate_constants[i] = (double)i + 20.0;
  CHECK(chem_state_copy(&copy, &original) == 0);
  CHECK(copy.cell_count == original.cell_count && copy.species_count == original.species_count);
  CHECK(copy.parameter_count == original.parameter_count && copy.rate_constant_count == original.rate_constant_count);
  CHECK(copy.relative_tolerance == original.relative_tolerance);
  CHECK(memcmp(copy.absolute_tolerances, original.absolute_tolerances, 4u * sizeof(double)) == 0);
  CHECK(memcmp(copy.concentrations, original.concentrations, 12u * sizeof(double)) == 0);
  CHECK(memcmp(copy.custom_parameters, original.custom_parameters, 6u * sizeof(double)) == 0);
  CHECK(memcmp(copy.rate_constants, original.rate_constants, 30u * sizeof(double)) == 0);
  CHECK(copy.concentrations != original.concentrations && copy.species_names[0] != original.species_names[0]);
  chem_state_destroy(&original);
  chem_state_destroy(&copy);
  return 0;
}

static int move_state(void) {
  chem_state original = { 0 }, moved = { 0 };
  CHECK(init_state(&original, 3u, parameters, 2u) == 0);
  original.concentrations[4] = 42.0;
  chem_state_move(&moved, &original);
  CHECK(moved.cell_count == 3u && moved.species_count == 4u && moved.parameter_count == 2u);
  CHECK(moved.concentrations[4] == 42.0 && strcmp(moved.species_names[0], "foo") == 0);
  CHECK(original.cell_count == 0u && original.concentrations == NULL && original.species_names == NULL);
  chem_state_destroy(&original);
  chem_state_destroy(&moved);
  return 0;
}

static int invalid_concentration_arguments(void) {
  chem_state state = { 0 };
  const double wrong[] = { 1.0, 2.0 };
  const double right[] = { 1.0, 2.0, 3.0, 4.0 };
  CHECK(init_state(&state, 3u, parameters, 2u) == 0);
  CHECK(chem_state_set_concentration(&state, "foo", 1.0) != 0);
  CHECK(chem_state_set_concentrations(&state, "foo", wrong, 2u) != 0);
  CHECK(chem_state_set_concentration(&state, "not foo", 1.0) != 0);
  CHECK(chem_state_set_concentrations(&state, "not foo", right, 4u) != 0);
  chem_state_destroy(&state);
  return 0;
}

static int set_single_concentration(void) {
  chem_state state = { 0 };
  const double concentrations[] = { 12.0, 42.0, 35.2 };
  double value = 0.0;
  CHECK(init_state(&state, 3u, parameters, 2u) == 0);
  CHECK(chem_state_set_concentrations(&state, "bar", concentrations, 3u) == 0);
  for (size_t i = 0; i < 3u; ++i) {
    CHECK(chem_state_get_concentration(&state, "bar", i, &value) == 0);
    CHECK(value == concentrations[i]);
  }
  chem_state_destroy(&state);
  memset(&state, 0, sizeof(state));
  CHECK(init_state(&state, 1u, parameters, 2u) == 0);
  CHECK(chem_state_set_concentration(&state, "bar", 324.2) == 0);
  CHECK(chem_state_get_concentration(&state, "bar", 0u, &value) == 0 && value == 324.2);
  chem_state_destroy(&state);
  return 0;
}

static int set_concentration_element_single(void) {
  chem_state state = { 0 };
  double baz = 0.0, quz = 0.0;
  CHECK(init_state(&state, 1u, parameters, 2u) == 0);
  CHECK(chem_state_set_concentration(&state, "baz", 42.0) == 0);
  CHECK(chem_state_set_concentration(&state, "quz", 12.0) == 0);
  CHECK(chem_state_get_concentration(&state, "baz", 0u, &baz) == 0 && baz == 42.0);
  CHECK(chem_state_get_concentration(&state, "quz", 0u, &quz) == 0 && quz == 12.0);
  chem_state_destroy(&state);
  return 0;
}

static int set_concentration_element_vector(void) {
  chem_state state = { 0 };
  const double values[] = { 12.0, 42.0, 35.2 };
  double actual = 0.0;
  CHECK(init_state(&state, 3u, parameters, 2u) == 0);
  CHECK(chem_state_set_concentrations(&state, "foo", values, 3u) == 0);
  for (size_t cell = 0; cell < 3u; ++cell) {
    CHECK(chem_state_get_concentration(&state, "foo", cell, &actual) == 0);
    CHECK(actual == values[cell]);
  }
  chem_state_destroy(&state);
  return 0;
}

static int invalid_bulk_concentrations(void) {
  chem_state state = { 0 };
  const double one[] = { 0.1 }, bar[] = { 0.2, 0.2, 0.2 };
  chem_named_values values[] = { { "FUU", one, 1u }, { "bar", bar, 3u } };
  CHECK(init_state(&state, 3u, parameters, 2u) == 0);
  CHECK(chem_state_set_all_concentrations(&state, values, 2u) != 0);
  CHECK(state.concentrations[1] == 0.0);
  chem_state_destroy(&state);
  return 0;
}

static int set_bulk_concentrations(void) {
  chem_state state = { 0 };
  const double bar[] = { 0.2, 0.22, 0.222 }, baz[] = { 0.3, 0.33, 0.333 };
  const double foo[] = { 0.9, 0.99, 0.999 }, quz[] = { 0.4, 0.44, 0.444 };
  const chem_named_values values[] = {
    { "bar", bar, 3u }, { "baz", baz, 3u }, { "foo", foo, 3u }, { "quz", quz, 3u }
  };
  const double expected[] = { 0.9, 0.2, 0.3, 0.4, 0.99, 0.22, 0.33, 0.44, 0.999, 0.222, 0.333, 0.444 };
  CHECK(init_state(&state, 3u, parameters, 2u) == 0);
  CHECK(chem_state_set_all_concentrations(&state, values, 4u) == 0);
  CHECK(memcmp(state.concentrations, expected, sizeof(expected)) == 0);
  chem_state_destroy(&state);
  return 0;
}

static int bracket_single_cell(void) {
  chem_state state = { 0 };
  double value = 0.0;
  CHECK(init_state(&state, 1u, parameters, 2u) == 0);
  CHECK(chem_state_set_concentration(&state, "foo", 42.0) == 0);
  CHECK(chem_state_set_concentration(&state, "bar", 12.0) == 0);
  CHECK(chem_state_set_concentration(&state, "baz", 35.2) == 0);
  CHECK(chem_state_set_concentration(&state, "quz", 24.2) == 0);
  CHECK(chem_state_get_concentration(&state, "foo", 0u, &value) == 0 && value == 42.0);
  CHECK(chem_state_get_concentration(&state, "bar", 0u, &value) == 0 && value == 12.0);
  CHECK(chem_state_get_concentration(&state, "baz", 0u, &value) == 0 && value == 35.2);
  CHECK(chem_state_get_concentration(&state, "quz", 0u, &value) == 0 && value == 24.2);
  CHECK(chem_state_set_concentration(&state, "foo", 42.0 + 8.0) == 0);
  CHECK(chem_state_get_concentration(&state, "foo", 0u, &value) == 0 && value == 50.0);
  CHECK(chem_state_set_concentration(&state, "foo", value - 20.0) == 0);
  CHECK(chem_state_get_concentration(&state, "foo", 0u, &value) == 0 && value == 30.0);
  CHECK(chem_state_set_concentration(&state, "foo", value * 2.0) == 0);
  CHECK(chem_state_get_concentration(&state, "foo", 0u, &value) == 0 && value == 60.0);
  CHECK(chem_state_set_concentration(&state, "foo", value / 3.0) == 0);
  CHECK(chem_state_get_concentration(&state, "foo", 0u, &value) == 0 && value == 20.0);
  chem_state_destroy(&state);
  return 0;
}

static int bracket_multi_cell(void) {
  chem_state state = { 0 };
  const double foo[] = { 1.0, 2.0, 3.0 }, bar[] = { 4.0, 5.0, 6.0 };
  double value = 0.0;
  CHECK(init_state(&state, 3u, parameters, 2u) == 0);
  CHECK(chem_state_set_concentrations(&state, "foo", foo, 3u) == 0);
  CHECK(chem_state_set_concentrations(&state, "bar", bar, 3u) == 0);
  for (size_t cell = 0; cell < 3u; ++cell) {
    CHECK(chem_state_get_concentration(&state, "foo", cell, &value) == 0 && value == foo[cell]);
    CHECK(chem_state_get_concentration(&state, "bar", cell, &value) == 0 && value == bar[cell]);
  }
  chem_state_destroy(&state);
  return 0;
}

static int invalid_custom_size(void) {
  chem_state state = { 0 };
  const double bad[] = { 0.3 };
  const chem_named_values values[] = {
    { "O3", bad, 1u }, { "O1", bad, 1u }, { "O2", bad, 1u }, { "AAA", bad, 1u }, { "BBB", bad, 1u }
  };
  CHECK(init_state(&state, 3u, five_parameters, 5u) == 0);
  CHECK(chem_state_set_all_custom_parameters(&state, values, 5u) != 0);
  chem_state_destroy(&state);
  return 0;
}

static int invalid_custom_label(void) {
  chem_state state = { 0 };
  const double val[] = { 0.7 };
  const chem_named_values values[] = {
    { "O3", val, 1u }, { "O1", val, 1u }, { "O2", val, 1u }, { "CCC", val, 1u }, { "AAA", val, 1u }
  };
  CHECK(init_state(&state, 1u, five_parameters, 5u) == 0);
  CHECK(chem_state_set_all_custom_parameters(&state, values, 5u) != 0);
  chem_state_destroy(&state);
  return 0;
}

static int set_custom_single(void) {
  chem_state state = { 0 };
  CHECK(init_state(&state, 1u, five_parameters, 5u) == 0);
  CHECK(chem_state_set_custom_parameter(&state, "O2", 42.3) == 0);
  CHECK(state.custom_parameters[1] == 42.3);
  chem_state_destroy(&state);
  return 0;
}

static int set_custom_bulk(void) {
  chem_state state = { 0 };
  const double o3[] = { 0.3, 0.33, 0.333 }, o1[] = { 0.1, 0.11, 0.111 };
  const double o2[] = { 0.5, 0.55, 0.555 }, aaa[] = { 0.5, 0.55, 0.555 };
  const double bbb[] = { 0.7, 0.77, 0.777 };
  const chem_named_values values[] = {
    { "O3", o3, 3u }, { "O1", o1, 3u }, { "O2", o2, 3u }, { "AAA", aaa, 3u }, { "BBB", bbb, 3u }
  };
  const double expected[] = { 0.1, 0.5, 0.3, 0.5, 0.7, 0.11, 0.55, 0.33, 0.55, 0.77,
                              0.111, 0.555, 0.333, 0.555, 0.777 };
  CHECK(init_state(&state, 3u, five_parameters, 5u) == 0);
  CHECK(chem_state_set_all_custom_parameters(&state, values, 5u) == 0);
  CHECK(memcmp(state.custom_parameters, expected, sizeof(expected)) == 0);
  chem_state_destroy(&state);
  return 0;
}

static int raw_custom_parameters(void) {
  chem_state state = { 0 };
  const double one_cell[] = { 0.1, 0.2, 0.3, 0.4, 0.5 };
  const double multi_cell[] = { 0.1, 0.2, 0.3, 0.4, 0.5,
                                0.1, 0.2, 0.3, 0.4, 0.5,
                                0.1, 0.2, 0.3, 0.4, 0.5 };
  CHECK(init_state(&state, 1u, five_parameters, 5u) == 0);
  CHECK(chem_state_set_custom_parameters_raw(&state, one_cell, 1u, 5u) == 0);
  CHECK(memcmp(state.custom_parameters, one_cell, sizeof(one_cell)) == 0);
  chem_state_destroy(&state);
  memset(&state, 0, sizeof(state));
  CHECK(init_state(&state, 3u, five_parameters, 5u) == 0);
  CHECK(chem_state_set_custom_parameters_raw(&state, multi_cell, 3u, 5u) == 0);
  CHECK(memcmp(state.custom_parameters, multi_cell, sizeof(multi_cell)) == 0);
  chem_state_destroy(&state);
  return 0;
}

static int raw_custom_invalid_rows(void) {
  chem_state state = { 0 };
  const double value[] = { 0.1, 0.2, 0.3, 0.4, 0.5 };
  CHECK(init_state(&state, 2u, five_parameters, 5u) == 0);
  CHECK(chem_state_set_custom_parameters_raw(&state, value, 1u, 5u) != 0);
  chem_state_destroy(&state);
  return 0;
}

static int raw_custom_invalid_columns(void) {
  chem_state state = { 0 };
  const double value[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };
  CHECK(init_state(&state, 2u, five_parameters, 5u) == 0);
  CHECK(chem_state_set_custom_parameters_raw(&state, value, 2u, 3u) != 0);
  chem_state_destroy(&state);
  return 0;
}

typedef int (*state_test)(void);
typedef struct { const char *name; state_test test; } named_test;

int main(int argc, char **argv) {
  static const named_test tests[] = {
    { "DefaultConstructor", default_constructor }, { "Constructor", constructor },
    { "CopyConstructor", copy_state }, { "CopyAssignmentOperator", copy_state },
    { "MoveConstructor", move_state }, { "MoveAssignmentOperator", move_state },
    { "SettingSingleConcentrationWithInvalidArgumentsThowsException", invalid_concentration_arguments },
    { "SetSingleConcentration", set_single_concentration },
    { "SetConcentrationByElementSingleValue", set_concentration_element_single },
    { "SetConcentrationByElementVector", set_concentration_element_vector },
    { "SettingConcentrationsWithInvalidArguementsThrowsException", invalid_bulk_concentrations },
    { "SetConcentrations", set_bulk_concentrations },
    { "SetStateWithSquareBracketOperator", bracket_single_cell },
    { "SetMultiCellStateWithSquareBrackets", bracket_multi_cell },
    { "SettingCustomRateParameterWithInvalidVectorSizeThrowsException", invalid_custom_size },
    { "SettingCustomRateParameterWithInvalidLabelThrowsException", invalid_custom_label },
    { "SetCustomRateParameter", set_custom_single }, { "SetCustomRateParameters", set_custom_bulk },
    { "UnsafelySetCustomRateParameterOneCell", raw_custom_parameters },
    { "UnsafelySetCustomRateParameterMultiCell", raw_custom_parameters },
    { "UnsafelySetCustomRateParameterCatchesTooFewGridCells", raw_custom_invalid_rows },
    { "UnsafelySetCustomRateParameterCatchesTooParameters", raw_custom_invalid_columns },
  };
  if (argc != 2) return 2;
  const char *requested = strrchr(argv[1], '.');
  requested = requested == NULL ? argv[1] : requested + 1;
  for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    if (strcmp(requested, tests[i].name) == 0) return tests[i].test();
  fprintf(stderr, "unknown state case: %s\n", argv[1]);
  return 2;
}
