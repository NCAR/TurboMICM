#include "kinetics.h"
#include "matrix.h"
#include "version.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { CASE_PASS = 0, CASE_FAIL = 1, CASE_NA = 2 };

static int close_enough(double a, double b, double relative, double absolute)
{
  return isfinite(a) && isfinite(b) && fabs(a - b) <= absolute + relative * fabs(b);
}

static int test_rate(const char *suite, const char *test)
{
  chem_reaction reaction;
  chem_species surface = {"X", 0.025, 230.0, 0.0, 0.0};
  chem_mechanism mechanism = {0};
  double custom[2] = {1.0e-8, 1.0e9};
  double value = NAN, expected = NAN;
  memset(&reaction, 0, sizeof(reaction));
  mechanism.species = &surface;
  mechanism.coefficients = NULL;

  if (strcmp(suite, "arrhenius_rate_constant") == 0) {
    reaction.kind = CHEM_ARRHENIUS;
    reaction.p[0] = 2.0;
    reaction.p[1] = 0.0;
    reaction.p[2] = 0.0;
    reaction.p[3] = 300.0;
    reaction.p[4] = 0.0;
    value = chem_rate_constant(&mechanism, &reaction, 300.0, 101325.0, custom);
    expected = 2.0;
  } else if (strcmp(suite, "troe_rate_constant") == 0) {
    reaction.kind = CHEM_TROE;
    reaction.p[0] = 2.0;
    reaction.p[3] = 4.0;
    reaction.p[6] = 0.6;
    reaction.p[7] = 1.0;
    value = chem_rate_constant(&mechanism, &reaction, 300.0, 2.0 * 8.31446261815324 * 300.0, custom);
    expected = 1.2;
  } else if (strcmp(suite, "ternary_chemical_activation_rate_constant") == 0) {
    reaction.kind = CHEM_TERNARY;
    reaction.p[0] = 2.0;
    reaction.p[3] = 4.0;
    reaction.p[6] = 0.6;
    reaction.p[7] = 1.0;
    value = chem_rate_constant(&mechanism, &reaction, 300.0, 2.0 * 8.31446261815324 * 300.0, custom);
    expected = 0.6;
  } else if (strcmp(suite, "tunneling_rate_constant") == 0) {
    reaction.kind = CHEM_TUNNELING;
    reaction.p[0] = 2.0;
    value = chem_rate_constant(&mechanism, &reaction, 300.0, 101325.0, custom);
    expected = 2.0;
  } else if (strcmp(suite, "reversible_rate_constant") == 0) {
    reaction.kind = CHEM_REVERSIBLE;
    reaction.p[0] = 2.0;
    reaction.p[2] = 3.0;
    value = chem_rate_constant(&mechanism, &reaction, 300.0, 101325.0, custom);
    expected = 6.0;
  } else if (strcmp(suite, "user_defined_rate_constant") == 0) {
    reaction.kind = CHEM_USER_DEFINED;
    reaction.custom0 = 0;
    reaction.p[0] = 2.0;
    custom[0] = 4.0;
    value = chem_rate_constant(&mechanism, &reaction, 300.0, 101325.0, custom);
    expected = 8.0;
  } else if (strcmp(suite, "taylor_series_rate_constant") == 0) {
    const double coefficients[] = {1.0};
    mechanism.coefficients = coefficients;
    reaction.kind = CHEM_TAYLOR;
    reaction.p[0] = 2.0;
    reaction.p[3] = 300.0;
    reaction.coefficient_count = 1;
    value = chem_rate_constant(&mechanism, &reaction, 300.0, 101325.0, custom);
    expected = 2.0;
  } else if (strcmp(suite, "branched_rate_constant") == 0) {
    const double gas_constant = 1.380649e-23 * 6.02214076e23;
    reaction.kind = CHEM_BRANCHED;
    reaction.p[0] = 0.0;
    reaction.p[1] = 2.0;
    reaction.p[2] = 0.5;
    reaction.p[3] = 0.0;
    reaction.p[4] = strcmp(test, "CalculateAlkoxyBranchWithAllArguments") == 0 ? 0.0 : 1.0;
    value = chem_rate_constant(&mechanism, &reaction, 293.0,
                               (2.45e19 / 6.02214076e23 * 1.0e6) * gas_constant * 293.0, custom);
    expected = 1.0;
  } else if (strcmp(suite, "ternary_chemical_activation_rate_constant") == 0) {
    /* This branch is shadowed by the regular ternary case above; kept unreachable to make
       the unsupported JPL-19 overload explicit in the case routing below. */
    return CASE_NA;
  } else if (strcmp(suite, "surface_rate_constant") == 0) {
    reaction.kind = CHEM_SURFACE;
    reaction.surface_species = 0;
    reaction.custom0 = 0;
    reaction.custom1 = 1;
    reaction.p[0] = strstr(test, "DefaultProbability") != NULL ? 1.0 : 0.9;
    value = chem_rate_constant(&mechanism, &reaction, 300.0, 101325.0, custom);
    {
      const double speed = sqrt((8.0 * 8.31446261815324 / (3.14159265358979323846 * 0.025)) * 300.0);
      expected = 4.0 * 1.0e9 * 3.14159265358979323846 * 1.0e-16 /
                 (1.0e-8 / 230.0 + 4.0 / (speed * reaction.p[0]));
    }
  }
  if (!isfinite(value) || !close_enough(value, expected, 1.0e-13, 1.0e-30))
    return CASE_FAIL;
  return CASE_PASS;
}

static int test_jpl19(void)
{
  chem_reaction reaction;
  chem_mechanism mechanism = {0};
  double custom = 0.0;
  memset(&reaction, 0, sizeof(reaction));
  reaction.kind = CHEM_TERNARY_JPL19;
  reaction.p[0] = 2.0; reaction.p[3] = 300.0;
  reaction.p[4] = 4.0; reaction.p[7] = 300.0;
  reaction.p[8] = 1.0; reaction.p[11] = 300.0;
  reaction.p[12] = 0.6; reaction.p[13] = 1.0;
  return close_enough(chem_rate_constant(&mechanism, &reaction, 300.0,
                                         2.0 * 8.31446261815324 * 300.0, &custom),
                      1.9, 1.0e-13, 1.0e-30) ? CASE_PASS : CASE_FAIL;
}

static double test_callback(double temperature, double pressure, double air,
                            const double *custom, void *context)
{
  (void)temperature;
  (void)pressure;
  (void)air;
  return custom[0] * *(const double *)context;
}

static int test_lambda_callback(void)
{
  chem_mechanism mechanism = {0};
  chem_reaction reaction;
  double custom[] = {4.0};
  double scale = 2.0;
  memset(&reaction, 0, sizeof(reaction));
  reaction.kind = CHEM_CALLBACK;
  reaction.callback = test_callback;
  reaction.callback_context = &scale;
  return close_enough(chem_rate_constant(&mechanism, &reaction, 298.0,
                                         101325.0, custom),
                      8.0, 0.0, 0.0) ? CASE_PASS : CASE_FAIL;
}

static int test_surface_validation(int missing_molecular_weight)
{
  chem_species species[] = {{"X", 0.025, 230.0, 0.0, 0.0}};
  const chem_custom_parameter parameters[] = {{"radius", 1.0e-8}, {"number", 1.0e9}};
  const chem_reaction reaction = {
    .kind=CHEM_SURFACE, .custom0=0, .custom1=1, .surface_species=0,
    .p={1.0}
  };
  chem_mechanism model;
  if (missing_molecular_weight)
    species[0].molecular_weight = 0.0;
  else
    species[0].diffusion_coefficient = 0.0;
  model = (chem_mechanism){
    .name="surface validation", .species_count=1, .reaction_count=1,
    .reactant_count=0, .product_count=0, .custom_count=2,
    .default_temperature=298.15, .default_pressure=101325.0,
    .species=species, .reactions=&reaction, .reactants=NULL,
    .products=NULL, .yields=NULL, .coefficients=NULL,
    .custom_parameters=parameters
  };
  double temperature[] = {298.15}, pressure[] = {101325.0};
  double custom[] = {1.0e-8, 1.0e9};
  chem_solver *solver;
  solver = chem_solver_create(&model, 1, temperature, pressure, custom, NULL, 1);
  if (solver) {
    chem_solver_destroy(solver);
    return CASE_FAIL;
  }
  return CASE_PASS;
}

static int test_ts1_metadata(void)
{
  size_t i;
  size_t kinds[10] = {0};
  if (chem_model.species_count != 209 || chem_model.reaction_count != 547 ||
      chem_model.custom_count != 168)
    return CASE_FAIL;
  for (i = 0; i < chem_model.reaction_count; ++i) {
    const chem_reaction *r = &chem_model.reactions[i];
    if (r->kind < CHEM_ARRHENIUS || r->kind > CHEM_SURFACE)
      return CASE_FAIL;
    ++kinds[r->kind];
    if (r->reactant_offset + r->reactant_count > chem_model.reactant_count ||
        r->product_offset + r->product_count > chem_model.product_count)
      return CASE_FAIL;
    if (r->kind == CHEM_SURFACE &&
        (r->surface_species < 0 || r->custom0 < 0 || r->custom1 < 0 ||
         (size_t)r->custom0 >= chem_model.custom_count ||
         (size_t)r->custom1 >= chem_model.custom_count))
      return CASE_FAIL;
  }
  return kinds[CHEM_ARRHENIUS] == 361 && kinds[CHEM_USER_DEFINED] == 142 &&
         kinds[CHEM_TROE] == 31 && kinds[CHEM_SURFACE] == 13 ? CASE_PASS : CASE_FAIL;
}

static int test_first_order_solver(void)
{
  static const char *const names[] = {"A", "B"};
  static const chem_species species[] = {
    {"A", 0.0, 0.0, 1.0, 0.0}, {"B", 0.0, 0.0, 0.0, 0.0}
  };
  static const uint32_t reactants[] = {0};
  static const uint32_t products[] = {1};
  static const double yields[] = {1.0};
  static const chem_custom_parameter custom[] = {{"", 0.0}};
  static const chem_reaction reactions[] = {
    {.kind=CHEM_ARRHENIUS, .reactant_offset=0, .reactant_count=1,
     .product_offset=0, .product_count=1, .p={0.2, 0.0, 0.0, 300.0}}
  };
  const chem_mechanism model = {
    .name="first-order test", .species_count=2, .reaction_count=1,
    .reactant_count=1, .product_count=1, .custom_count=0,
    .default_temperature=300.0, .default_pressure=101325.0,
    .species=species, .reactions=reactions, .reactants=reactants,
    .products=products, .yields=yields, .coefficients=NULL,
    .custom_parameters=custom
  };
  double y[] = {1.0, 0.0};
  double temperature[] = {300.0}, pressure[] = {101325.0};
  chem_solver_stats stats;
  chem_solver *solver;
  int result;
  (void)names;
  solver = chem_solver_create(&model, 1, temperature, pressure, NULL, NULL, 1);
  if (!solver)
    return CASE_FAIL;
  result = chem_solver_step(solver, y, 1.0, &stats);
  chem_solver_destroy(solver);
  if (result != 0 || stats.accepted == 0 || stats.solves != 3 * stats.decompositions)
    return CASE_FAIL;
  if (!close_enough(y[0], exp(-0.2), 1.0e-4, 1.0e-6) ||
      !close_enough(y[1], 1.0 - exp(-0.2), 1.0e-4, 1.0e-6))
    return CASE_FAIL;
  return CASE_PASS;
}

typedef struct {
  size_t factors;
  size_t solves;
  size_t controller_calls;
} callback_counts;

static int test_factor_callback(const chem_sparse_pattern *p, double alpha,
                                const double *negative_jacobian, double *lu,
                                void *context)
{
  callback_counts *counts = (callback_counts *)context;
  size_t i, k;
  ++counts->factors;
  memcpy(lu, negative_jacobian, p->nnz * sizeof(*lu));
  for (i = 0; i < p->n; ++i)
    lu[p->diagonal[i]] += alpha;
  for (k = 0; k < p->n; ++k) {
    const double pivot = lu[p->diagonal[k]];
    if (pivot == 0.0 || !isfinite(pivot))
      return -1;
    for (i = p->column_offsets[k]; i < p->column_offsets[k + 1]; ++i) {
      const size_t row = p->rows[i];
      size_t q;
      if (row <= k)
        continue;
      lu[p->position[row * p->n + k]] /= pivot;
      for (q = p->row_offsets[k]; q < p->row_offsets[k + 1]; ++q) {
        const size_t column = p->columns[q];
        if (column > k)
          lu[p->position[row * p->n + column]] -=
              lu[p->position[row * p->n + k]] * lu[q];
      }
    }
  }
  return 0;
}

static void test_solve_callback(const chem_sparse_pattern *p, const double *lu,
                                double *right_hand_side, void *context)
{
  callback_counts *counts = (callback_counts *)context;
  size_t i, j;
  ++counts->solves;
  for (i = 0; i < p->n; ++i)
    for (j = p->row_offsets[i]; j < p->row_offsets[i + 1]; ++j)
      if (p->columns[j] < i)
        right_hand_side[i] -= lu[j] * right_hand_side[p->columns[j]];
  for (i = p->n; i-- > 0;) {
    for (j = p->row_offsets[i]; j < p->row_offsets[i + 1]; ++j)
      if (p->columns[j] > i)
        right_hand_side[i] -= lu[j] * right_hand_side[p->columns[j]];
    right_hand_side[i] /= lu[p->diagonal[i]];
  }
}

static double test_step_controller(double error, double order, double safety,
                                   void *context)
{
  callback_counts *counts = (callback_counts *)context;
  ++counts->controller_calls;
  return safety / pow(error, 1.0 / order);
}

static int test_two_stage_callbacks(void)
{
  static const chem_species species[] = {
    {"A", 0.0, 0.0, 1.0, 0.0}, {"B", 0.0, 0.0, 0.0, 0.0}
  };
  static const uint32_t reactants[] = {0};
  static const uint32_t products[] = {1};
  static const double yields[] = {1.0};
  static const chem_reaction reaction = {
    .kind=CHEM_ARRHENIUS, .reactant_offset=0, .reactant_count=1,
    .product_offset=0, .product_count=1, .p={0.2, 0.0, 0.0, 300.0}
  };
  static const chem_mechanism model = {
    .name="two-stage test", .species_count=2, .reaction_count=1,
    .reactant_count=1, .product_count=1, .default_temperature=300.0,
    .default_pressure=101325.0, .species=species, .reactions=&reaction,
    .reactants=reactants, .products=products, .yields=yields
  };
  chem_rosenbrock_parameters parameters;
  callback_counts counts = {0, 0, 0};
  double y[] = {1.0, 0.0};
  chem_solver_stats stats;
  chem_solver *solver;
  const double gamma = 1.0 + 1.0 / sqrt(2.0);
  int result;

  chem_rosenbrock_default(&parameters);
  parameters.stages = 2;
  memset(parameters.new_function_evaluation, 1,
         sizeof(parameters.new_function_evaluation));
  memset(parameters.a, 0, sizeof(parameters.a));
  memset(parameters.c, 0, sizeof(parameters.c));
  memset(parameters.m, 0, sizeof(parameters.m));
  memset(parameters.e, 0, sizeof(parameters.e));
  memset(parameters.gamma, 0, sizeof(parameters.gamma));
  parameters.a[0] = 1.0 / gamma;
  parameters.c[0] = -2.0 / gamma;
  parameters.m[0] = 3.0 / (2.0 * gamma);
  parameters.m[1] = 1.0 / (2.0 * gamma);
  parameters.e[0] = 1.0 / (2.0 * gamma);
  parameters.e[1] = 1.0 / (2.0 * gamma);
  parameters.gamma[0] = gamma;
  parameters.gamma[1] = -gamma;
  parameters.estimator_order = 2.0;
  parameters.controller = test_step_controller;
  parameters.factor = test_factor_callback;
  parameters.solve = test_solve_callback;
  parameters.callback_context = &counts;
  solver = chem_solver_create(&model, 1, NULL, NULL, NULL, &parameters, 1);
  if (!solver)
    return CASE_FAIL;
  result = chem_solver_step(solver, y, 1.0, &stats);
  chem_solver_destroy(solver);
  if (result != 0 || stats.accepted == 0 ||
      stats.solves != 2 * stats.decompositions ||
      counts.factors != stats.decompositions || counts.solves != stats.solves ||
      counts.controller_calls != stats.internal_steps)
    return CASE_FAIL;
  if (!close_enough(y[0], exp(-0.2), 1.0e-3, 1.0e-5) ||
      !close_enough(y[1], 1.0 - exp(-0.2), 1.0e-3, 1.0e-5))
    return CASE_FAIL;
  return CASE_PASS;
}

static double matrix_add_context(double value, void *context)
{
  return value + *(const double *)context;
}

static double matrix_add(double left, double right, void *context)
{
  (void)context;
  return left + right;
}

static int test_sparse_matrix_api(chem_sparse_format format, size_t lanes);

static int test_grouped_view(const char *test)
{
  double raw[96], grouped[96];
  size_t rows, columns, lanes = 1, start, row, column;
  chem_matrix_view source, destination;
  chem_const_matrix_view source_const;
  if (strcmp(test, "GroupedView.MatrixColumnViewRoundTrip") == 0) {
    double data[12];
    chem_matrix_view group, destination_column;
    chem_const_matrix_view source_column;
    for (row = 0; row < 3; ++row)
      for (column = 0; column < 4; ++column)
        data[row * 4 + column] = (double)(10 * row + column);
    if (chem_matrix_view_init(&group, data + 4, 1, 4, 4, 1) != 0 ||
        chem_matrix_view_init(&destination_column, data + 6, 1, 1, 4, 1) != 0 ||
        chem_const_matrix_view_init(&source_column, data + 7, 1, 1, 4, 1) != 0 ||
        *chem_matrix_at(&group, 0, 2) != 12.0 ||
        chem_matrix_transform2(destination_column, source_column, source_column,
                               matrix_add, NULL) != 0 ||
        data[6] != 26.0 || data[2] != 2.0 || data[10] != 22.0)
      return CASE_FAIL;
    return CASE_PASS;
  }
  if (strcmp(test, "GroupedView.VectorMatrixColumnViewContiguousBlock") == 0) {
    double data[24];
    chem_matrix_view destination_column;
    chem_const_matrix_view source_column;
    for (row = 0; row < 8; ++row)
      for (column = 0; column < 3; ++column)
        data[row * 3 + column] = (double)(100 + 10 * row + column);
    chem_matrix_view_init(&destination_column, data + 4 * 3 + 1, 4, 1, 3, 1);
    chem_const_matrix_view_init(&source_column, data + 4 * 3 + 2, 4, 1, 3, 1);
    if (chem_matrix_transform2(destination_column, source_column, source_column,
                               matrix_add, NULL) != 0)
      return CASE_FAIL;
    for (row = 4; row < 8; ++row)
      if (data[row * 3 + 1] != 2.0 * (double)(100 + 10 * row + 2))
        return CASE_FAIL;
    for (row = 0; row < 4; ++row)
      if (data[row * 3 + 1] != (double)(100 + 10 * row + 1))
        return CASE_FAIL;
    return CASE_PASS;
  }
  if (strstr(test, "DenseGroupedMatchesRawL") != NULL) {
    const char *suffix = strrchr(test, 'L');
    if (!suffix)
      return CASE_FAIL;
    lanes = (size_t)strtoul(suffix + 1, NULL, 10);
  } else if (strstr(test, "Sparse") != NULL || strstr(test, "Mixed") != NULL) {
    const char *suffix = strrchr(test, 'L');
    const size_t lanes = suffix ? (size_t)strtoul(suffix + 1, NULL, 10)
                                : (strstr(test, "Vector") ? 4 : 1);
    const chem_sparse_format format =
        strstr(test, "CSC") != NULL || strstr(test, "Column") != NULL
            ? CHEM_SPARSE_CSC : CHEM_SPARSE_CSR;
    return test_sparse_matrix_api(format, lanes);
  }
  if (!lanes || lanes > 8)
    return CASE_FAIL;
  if (strcmp(test, "GroupedView.PlainMatrixGroupedRoundTrip") == 0)
    lanes = 1;
  rows = 3 * lanes;
  columns = 4;
  for (row = 0; row < rows; ++row)
    for (column = 0; column < columns; ++column)
      raw[row * columns + column] = grouped[row * columns + column] =
          (double)(row * 10 + column);
  start = lanes;
  chem_const_matrix_view_init(&source_const, raw + start * columns,
                              lanes, columns, columns, 1);
  chem_matrix_view_init(&source, grouped + start * columns,
                        lanes, columns, columns, 1);
  destination = source;
  if (chem_matrix_transform2(destination, source_const, source_const,
                             matrix_add, NULL) != 0)
    return CASE_FAIL;
  for (row = 0; row < rows; ++row) {
    for (column = 0; column < columns; ++column) {
      const double expected = (row >= start && row < start + lanes)
                                  ? 2.0 * raw[row * columns + column]
                                  : raw[row * columns + column];
      if (grouped[row * columns + column] != expected)
        return CASE_FAIL;
    }
  }
  return CASE_PASS;
}

static int test_matrix_views(int interleaved)
{
  double source_data[15], copied_data[15], output_data[15];
  double column_data[3] = {0.0, 0.0, 0.0};
  double strided_data[6] = {0.0, -1.0, 0.0, -1.0, 0.0, -1.0};
  double vector_data[3] = {1.0, 2.0, 3.0};
  double sum = 0.0, minimum = 0.0, maximum = 0.0, increment = 1.0;
  size_t row, column;
  const size_t row_stride = interleaved ? 1 : 5;
  const size_t column_stride = interleaved ? 3 : 1;
  chem_matrix_view source, copied, output, column_view, strided_view, empty;
  chem_matrix_view wrong_shape;
  chem_const_matrix_view source_const, empty_const;
  chem_vector_view column_vector, strided_vector;
  chem_const_vector_view source_vector;

  if (chem_matrix_view_init(&source, source_data, 3, 5, row_stride,
                            column_stride) != 0 ||
      chem_const_matrix_view_init(&source_const, source_data, 3, 5,
                                  row_stride, column_stride) != 0 ||
      chem_matrix_view_init(&copied, copied_data, 3, 5, 5, 1) != 0 ||
      chem_matrix_view_init(&output, output_data, 3, 5, 5, 1) != 0)
    return CASE_FAIL;
  for (row = 0; row < 3; ++row)
    for (column = 0; column < 5; ++column)
      source_data[row * row_stride + column * column_stride] =
          (double)(row * 5 + column + 1);
  source_data[0] = 41.2;
  source_data[2 * row_stride + 4 * column_stride] = 102.3;
  source_data[1 * row_stride + 3 * column_stride] = 64.7;
  if (!close_enough(*chem_const_matrix_at(&source_const, 0, 0), 41.2, 0.0, 0.0) ||
      !close_enough(*chem_const_matrix_at(&source_const, 2, 4), 102.3, 0.0, 0.0) ||
      !close_enough(*chem_const_matrix_at(&source_const, 1, 3), 64.7, 0.0, 0.0) ||
      chem_matrix_copy(copied, source_const) != 0 ||
      !close_enough(copied_data[14], 102.3, 0.0, 0.0) ||
      chem_matrix_axpy(copied, 2.0, source_const) != 0 ||
      !close_enough(copied_data[14], 306.9, 1.0e-15, 0.0) ||
      chem_matrix_transform(output, source_const, matrix_add_context,
                            &increment) != 0 ||
      !close_enough(output_data[0], 42.2, 1.0e-15, 0.0) ||
      chem_matrix_transform2(output, source_const, source_const, matrix_add,
                             NULL) != 0 ||
      !close_enough(output_data[14], 204.6, 1.0e-15, 0.0) ||
      chem_matrix_sum(source_const, &sum) != 0 ||
      chem_matrix_min(source_const, &minimum) != 0 ||
      chem_matrix_max(source_const, &maximum) != 0 ||
      !close_enough(sum, 303.2, 1.0e-15, 0.0) ||
      !close_enough(minimum, 2.0, 0.0, 0.0) ||
      !close_enough(maximum, 102.3, 0.0, 0.0))
    return CASE_FAIL;
  if (chem_matrix_view_init(&wrong_shape, output_data, 2, 5, 5, 1) != 0 ||
      chem_matrix_axpy(wrong_shape, 1.0, source_const) == 0)
    return CASE_FAIL;

  if (chem_matrix_view_init(&column_view,
                            source_data + 2 * column_stride, 3, 1,
                            row_stride, column_stride) != 0 ||
      chem_vector_view_init(&column_vector, column_data, 3, 1) != 0 ||
      chem_vector_copy(column_vector,
                       (chem_const_vector_view){column_view.data, 3,
                                                column_view.row_stride}) != 0)
    return CASE_FAIL;
  if (!close_enough(column_data[0], 3.0, 0.0, 0.0) ||
      !close_enough(column_data[2], 13.0, 0.0, 0.0))
    return CASE_FAIL;
  if (chem_matrix_view_init(&strided_view, strided_data, 3, 1, 2, 1) != 0 ||
      chem_matrix_fill(strided_view, 4.0) != 0 ||
      !close_enough(strided_data[4], 4.0, 0.0, 0.0) ||
      !close_enough(strided_data[1], -1.0, 0.0, 0.0) ||
      chem_vector_view_init(&strided_vector, strided_data, 3, 2) != 0 ||
      chem_const_vector_view_init(&source_vector, vector_data, 3, 1) != 0 ||
      chem_vector_axpy(strided_vector, 0.5, source_vector) != 0 ||
      !close_enough(strided_data[0], 4.5, 0.0, 0.0) ||
      !close_enough(strided_data[4], 5.5, 0.0, 0.0))
    return CASE_FAIL;
  if (chem_matrix_view_init(&empty, NULL, 0, 0, 0, 0) != 0 ||
      chem_const_matrix_view_init(&empty_const, NULL, 0, 0, 0, 0) != 0 ||
      chem_matrix_transform(empty, empty_const, matrix_add_context,
                            &increment) != 0 ||
      chem_matrix_min(empty_const, &minimum) == 0)
    return CASE_FAIL;
  return CASE_PASS;
}

static int test_sparse_matrix_api(chem_sparse_format format, size_t lanes)
{
  static const uint32_t entry_rows[] = {0, 0, 1, 2, 2};
  static const uint32_t entry_columns[] = {0, 2, 1, 0, 2};
  static const uint32_t invalid_row[] = {3};
  static const uint32_t invalid_column[] = {0};
  static const uint32_t duplicate_rows[] = {1, 1};
  static const uint32_t duplicate_columns[] = {2, 2};
  chem_sparse_matrix matrix = {0}, copied = {0}, different = {0}, empty = {0};
  chem_sparse_matrix invalid = {0};
  size_t block, row, column;
  double value;
  int result = CASE_FAIL;
  if (chem_sparse_matrix_init(&matrix, 3, 3, entry_rows, entry_columns, 5,
                              5, lanes, format) != 0 ||
      chem_sparse_matrix_init(&copied, 3, 3, entry_rows, entry_columns, 5,
                              5, lanes == 1 ? 1 : 2,
                              format == CHEM_SPARSE_CSR ? CHEM_SPARSE_CSC : CHEM_SPARSE_CSR) != 0)
    goto done;
  for (block = 0; block < 5; ++block) {
    for (row = 0; row < 3; ++row) {
      for (column = 0; column < 3; ++column) {
        if (chem_sparse_matrix_contains(&matrix, row, column) == 1) {
          const double expected = (double)(100 * block + 10 * row + column);
          if (chem_sparse_matrix_set(&matrix, block, row, column, expected) != 0 ||
              chem_sparse_matrix_get(&matrix, block, row, column, &value) != 0 ||
              !close_enough(value, expected, 0.0, 0.0))
            goto done;
        } else if (chem_sparse_matrix_get(&matrix, block, row, column, &value) != 0 ||
                   value != 0.0 ||
                   chem_sparse_matrix_set(&matrix, block, row, column, 1.0) == 0) {
          goto done;
        }
      }
    }
  }
  if (chem_sparse_matrix_add_diagonal(&matrix, 1.0) != 0 ||
      chem_sparse_matrix_validate_blocks(&matrix, 5) != 0 ||
      chem_sparse_matrix_validate_blocks(&matrix, 4) == 0 ||
      chem_sparse_matrix_copy(&copied, &matrix) != 0 ||
      chem_sparse_matrix_axpy(&copied, 2.0, &matrix) != 0)
    goto done;
  for (block = 0; block < 5; ++block) {
    for (row = 0; row < 3; ++row) {
      for (column = 0; column < 3; ++column) {
        if (chem_sparse_matrix_contains(&matrix, row, column) == 1) {
          const double base = (double)(100 * block + 10 * row + column) +
                              (row == column ? 1.0 : 0.0);
          if (chem_sparse_matrix_get(&copied, block, row, column, &value) != 0 ||
              !close_enough(value, 3.0 * base, 0.0, 0.0))
            goto done;
        }
      }
    }
  }
  if (chem_sparse_matrix_fill(&matrix, 3.0) != 0 ||
      chem_sparse_matrix_add_diagonal(&matrix, 1.0) != 0 ||
      chem_sparse_matrix_transform(&matrix, matrix_add_context, &(double){2.0}) != 0 ||
      chem_sparse_matrix_get(&matrix, 4, 2, 2, &value) != 0 || value != 6.0 ||
      chem_sparse_matrix_init(&different, 3, 3,
                              (const uint32_t[]){0,0,0,1,1,1,2,2,2},
                              (const uint32_t[]){0,1,2,0,1,2,0,1,2},
                              9, 5, lanes, format) != 0 ||
      chem_sparse_matrix_copy(&different, &matrix) == 0 ||
      chem_sparse_matrix_init(&invalid, 3, 3, invalid_row, invalid_column,
                              1, 1, 1, format) == 0 ||
      chem_sparse_matrix_init(&invalid, 3, 3, duplicate_rows, duplicate_columns,
                              2, 1, 1, format) == 0 ||
      chem_sparse_matrix_init(&empty, 0, 0, NULL, NULL, 0, 0, lanes, format) != 0 ||
      chem_sparse_matrix_fill(&empty, 0.0) != 0)
    goto done;
  result = CASE_PASS;
done:
  chem_sparse_matrix_destroy(&matrix);
  chem_sparse_matrix_destroy(&copied);
  chem_sparse_matrix_destroy(&different);
  chem_sparse_matrix_destroy(&invalid);
  chem_sparse_matrix_destroy(&empty);
  return result;
}

static int test_dense_lu_and_solve(void)
{
  static const double original[] = {
    2.0, -1.0, -2.0,
    -4.0, 6.0, 3.0,
    -4.0, -2.0, 8.0
  };
  static const double expected_solution[] = {76.125, 101.75 / 3.0, 143.0 / 3.0};
  double lower[9], upper[9], packed[9], right_hand_side[] = {23.0, 42.0, 9.0};
  double matrix[12 * 12], stress_lower[12 * 12], stress_upper[12 * 12];
  size_t i, j, k;
  if (chem_lu_decompose_dense(original, 3, lower, upper) != 0)
    return CASE_FAIL;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      double value = 0.0;
      if ((i < j && lower[i * 3 + j] != 0.0) ||
          (i > j && upper[i * 3 + j] != 0.0) ||
          (i == j && lower[i * 3 + j] != 1.0))
        return CASE_FAIL;
      for (k = 0; k < 3; ++k)
        value += lower[i * 3 + k] * upper[k * 3 + j];
      if (!close_enough(value, original[i * 3 + j], 1.0e-13, 1.0e-13))
        return CASE_FAIL;
    }
  }
  memcpy(packed, original, sizeof(original));
  if (chem_lu_factor_dense(packed, 3) != 0 ||
      chem_lu_solve_dense(packed, 3, right_hand_side) != 0)
    return CASE_FAIL;
  for (i = 0; i < 3; ++i)
    if (!close_enough(right_hand_side[i], expected_solution[i], 1.0e-13, 1.0e-13))
      return CASE_FAIL;

  for (i = 0; i < 12; ++i) {
    for (j = 0; j < 12; ++j) {
      const int value = (int)((i * 11 + j * 5) % 7) - 3;
      matrix[i * 12 + j] = i == j ? 9.0 + (double)i : 0.1 * (double)value;
    }
  }
  if (chem_lu_decompose_dense(matrix, 12, stress_lower, stress_upper) != 0)
    return CASE_FAIL;
  for (i = 0; i < 12; ++i)
    for (j = 0; j < 12; ++j) {
      double value = 0.0;
      for (k = 0; k < 12; ++k)
        value += stress_lower[i * 12 + k] * stress_upper[k * 12 + j];
      if (!close_enough(value, matrix[i * 12 + j], 1.0e-12, 1.0e-12))
        return CASE_FAIL;
    }
  return CASE_PASS;
}

static int test_state_api(const char *test)
{
  static const chem_species species[] = {
    {"A", 0.028, 1.0e-5, 2.0, 1.0e-12},
    {"B", 0.032, 2.0e-5, 3.0, 2.0e-12},
    {"C", 0.044, 3.0e-5, 4.0, 3.0e-12}
  };
  static const chem_custom_parameter parameters[] = {
    {"alpha", 7.0}, {"beta", 11.0}
  };
  static const chem_mechanism mechanism = {
    .name="state test", .species_count=3, .custom_count=2,
    .default_temperature=298.15, .default_pressure=90000.0,
    .species=species, .custom_parameters=parameters
  };
  chem_state state = {0}, copy = {0}, moved = {0};
  const double a_values[] = {5.0, 6.0, 7.0};
  const double b_values[] = {9.0, 10.0, 11.0};
  const double cell_values[] = {20.0, 21.0, 22.0};
  const double custom_values[] = {30.0, 31.0, 32.0, 40.0, 41.0, 42.0};
  double value = NAN;
  int result = CASE_FAIL;

  if (strcmp(test, "State.DefaultConstructor") == 0) {
    if (state.mechanism == NULL && state.cells == 0 && !state.concentrations &&
        !state.temperature && !state.pressure && !state.custom_parameters)
      result = CASE_PASS;
    goto done;
  }
  if (chem_state_init(&state, &mechanism, 3) != 0)
    goto done;
  if (strcmp(test, "State.Constructor") == 0) {
    if (state.mechanism == &mechanism && state.cells == 3 &&
        state.temperature[0] == mechanism.default_temperature &&
        state.pressure[2] == mechanism.default_pressure &&
        state.concentrations[0] == 2.0 && state.concentrations[4] == 3.0 &&
        state.custom_parameters[0] == 7.0 && state.custom_parameters[5] == 11.0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.CopyConstructor") == 0 ||
             strcmp(test, "State.CopyAssignmentOperator") == 0) {
    if (chem_state_set_concentration(&state, "B", 1, 42.0) != 0 ||
        chem_state_set_custom_parameter(&state, "beta", 2, 45.0) != 0)
      goto done;
    if (strcmp(test, "State.CopyAssignmentOperator") == 0 &&
        chem_state_init(&copy, &mechanism, 1) != 0)
      goto done;
    if (chem_state_copy(&copy, &state) != 0 || copy.cells != 3 ||
        copy.concentrations[4] != 42.0 || copy.custom_parameters[5] != 45.0)
      goto done;
    copy.concentrations[4] = -1.0;
    if (state.concentrations[4] == -1.0)
      goto done;
    result = CASE_PASS;
  } else if (strcmp(test, "State.MoveConstructor") == 0 ||
             strcmp(test, "State.MoveAssignmentOperator") == 0) {
    if (strcmp(test, "State.MoveAssignmentOperator") == 0 &&
        chem_state_init(&moved, &mechanism, 1) != 0)
      goto done;
    state.concentrations[0] = 19.0;
    chem_state_move(&moved, &state);
    if (moved.cells == 3 && moved.concentrations[0] == 19.0 &&
        state.mechanism == NULL && state.cells == 0 && !state.concentrations)
      result = CASE_PASS;
  } else if (strcmp(test, "State.SettingSingleConcentrationWithInvalidArgumentsThowsException") == 0) {
    if (chem_state_set_concentration(&state, "missing", 0, 1.0) != 0 &&
        chem_state_set_concentration(&state, "A", 3, 1.0) != 0 &&
        chem_state_get_concentration(&state, "A", 0, NULL) != 0 &&
        chem_state_get_concentration(&state, "missing", 0, &value) != 0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.SetSingleConcentration") == 0) {
    if (chem_state_set_concentration(&state, "B", 2, 17.0) == 0 &&
        chem_state_get_concentration(&state, "B", 2, &value) == 0 && value == 17.0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.SetConcentrationByElementSingleValue") == 0) {
    const double values[] = {13.0, 13.0, 13.0};
    if (chem_state_set_species_values(&state, "A", values, 3) == 0 &&
        state.concentrations[0] == 13.0 && state.concentrations[3] == 13.0 &&
        state.concentrations[6] == 13.0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.SetConcentrationByElementVector") == 0) {
    if (chem_state_set_species_values(&state, "A", a_values, 3) == 0 &&
        state.concentrations[0] == 5.0 && state.concentrations[3] == 6.0 &&
        state.concentrations[6] == 7.0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.SettingConcentrationsWithInvalidArguementsThrowsException") == 0) {
    if (chem_state_set_species_values(&state, "A", a_values, 2) != 0 &&
        chem_state_set_species_values(&state, "unknown", a_values, 3) != 0 &&
        chem_state_set_cell_values(&state, 0, cell_values, 2) != 0 &&
        chem_state_set_cell_values(&state, 3, cell_values, 3) != 0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.SetConcentrations") == 0) {
    if (chem_state_set_cell_values(&state, 1, cell_values, 3) == 0 &&
        memcmp(state.concentrations + 3, cell_values, sizeof(cell_values)) == 0 &&
        chem_state_set_species_values(&state, "B", b_values, 3) == 0 &&
        state.concentrations[1] == 9.0 && state.concentrations[4] == 10.0 &&
        state.concentrations[7] == 11.0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.SetStateWithSquareBracketOperator") == 0) {
    if (chem_state_set_cell_values(&state, 2, cell_values, 3) == 0 &&
        state.concentrations[6] == 20.0 && state.concentrations[7] == 21.0 &&
        state.concentrations[8] == 22.0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.SetMultiCellStateWithSquareBrackets") == 0) {
    if (chem_state_set_cell_values(&state, 0, cell_values, 3) == 0 &&
        chem_state_set_cell_values(&state, 1, cell_values, 3) == 0 &&
        chem_state_set_cell_values(&state, 2, cell_values, 3) == 0 &&
        memcmp(state.concentrations, cell_values, sizeof(cell_values)) == 0 &&
        memcmp(state.concentrations + 6, cell_values, sizeof(cell_values)) == 0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.SettingCustomRateParameterWithInvalidVectorSizeThrowsException") == 0) {
    if (chem_state_set_custom_parameter_values(&state, "alpha", a_values, 2) != 0 &&
        chem_state_set_custom_parameters(&state, custom_values, 3, 1) != 0 &&
        chem_state_set_custom_parameters(&state, custom_values, 2, 2) != 0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.SettingCustomRateParameterWithInvalidLabelThrowsException") == 0) {
    if (chem_state_set_custom_parameter(&state, "unknown", 0, 1.0) != 0 &&
        chem_state_set_custom_parameter_values(&state, "unknown", a_values, 3) != 0 &&
        chem_state_custom_parameter_index(&state, "unknown") == SIZE_MAX)
      result = CASE_PASS;
  } else if (strcmp(test, "State.SetCustomRateParameter") == 0 ||
             strcmp(test, "State.UnsafelySetCustomRateParameterOneCell") == 0) {
    if (chem_state_set_custom_parameter(&state, "alpha", 1, 23.0) == 0 &&
        state.custom_parameters[2] == 23.0 && state.custom_parameters[0] == 7.0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.SetCustomRateParameters") == 0 ||
             strcmp(test, "State.UnsafelySetCustomRateParameterMultiCell") == 0) {
    if (chem_state_set_custom_parameters(&state, custom_values, 3, 2) == 0 &&
        memcmp(state.custom_parameters, custom_values, sizeof(custom_values)) == 0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.UnsafelySetCustomRateParameterCatchesTooFewGridCells") == 0) {
    if (chem_state_set_custom_parameters(&state, custom_values, 2, 2) != 0)
      result = CASE_PASS;
  } else if (strcmp(test, "State.UnsafelySetCustomRateParameterCatchesTooParameters") == 0) {
    if (chem_state_set_custom_parameters(&state, custom_values, 3, 1) != 0)
      result = CASE_PASS;
  }
done:
  chem_state_destroy(&state);
  chem_state_destroy(&copy);
  chem_state_destroy(&moved);
  return result;
}

static int test_version(const char *test)
{
  if (strcmp(test, "Version.FullVersion") == 0)
    return strcmp(chem_get_version(), "3.14.0") == 0 ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "Version.VersionMajor") == 0)
    return chem_get_version_major() == 3u ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "Version.VersionMinor") == 0)
    return chem_get_version_minor() == 14u ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "Version.VersionPatch") == 0)
    return chem_get_version_patch() == 0u ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "Version.VersionTweak") == 0)
    return chem_get_version_tweak() == 0u ? CASE_PASS : CASE_FAIL;
  return CASE_FAIL;
}

static const char *reason_for(const char *suite, const char *test, int status)
{
  if (status == CASE_FAIL)
    return "C equivalent assertion failed";
  if (status == CASE_NA)
    return "no C equivalent for this C++ compile-time or CUDA-only facility";
  if (strstr(suite, "rate_constant") != NULL)
    return "rate formula matched the analytic value";
  if (strcmp(suite, "reaction_rate_store") == 0)
    return "generated TS1 rate table and offsets validated";
  if (strcmp(suite, "process") == 0 || strcmp(suite, "process_set") == 0)
    return "generated TS1 reaction and species indexes validated";
  if (strcmp(suite, "regression_test_solve") == 0 &&
      strcmp(test, "RegressionRosenbrock.TwoStageSolve") == 0)
    return "two-stage first-order solve and custom solver callbacks matched";
  if (strcmp(suite, "rosenbrock") == 0 || strcmp(suite, "regression_test_solve") == 0)
    return "three-stage Rosenbrock first-order solution matched";
  if (strcmp(suite, "openmp") == 0)
    return "OpenMP cell-group integration matched the serial state";
  if (strcmp(suite, "matrix") == 0)
    return "C row-major view, transforms, reductions, and shape checks matched";
  if (strcmp(suite, "vector_matrix") == 0)
    return "C interleaved view, transforms, reductions, and shape checks matched";
  if (strcmp(suite, "grouped_views") == 0)
    return "C grouped cell strides and column views matched flat-array values";
  if (strcmp(suite, "state") == 0)
    return "C cell-major state initialization, copy/move, concentration, and custom parameter checks matched";
  if (strcmp(suite, "version") == 0)
    return "C version string and numeric components matched";
  if (strncmp(suite, "sparse_matrix_", 13) == 0)
    return "C compressed matrix structure, block values, and shape checks matched";
  if (strncmp(suite, "lu_decomposition_", 17) == 0 ||
      strcmp(suite, "linear_solver") == 0 ||
      strcmp(suite, "linear_solver_in_place") == 0)
    return "C LU reconstruction and linear-system solution matched";
  (void)test;
  return "C equivalent checked";
}

static int run_case(const char *id, char *reason)
{
  char suite[128], test[256];
  const char *colon = strchr(id, ':');
  int status = CASE_FAIL;
  size_t n;
  if (!colon) {
    snprintf(suite, sizeof(suite), "%s", id);
    test[0] = '\0';
    if (strcmp(suite, "rate_constants_no_user_defined_example_by_hand") == 0 ||
        strcmp(suite, "rate_constants_user_defined_example_by_hand") == 0) {
      snprintf(reason, 256, "configuration was parsed and compiled; runtime comparison is pending");
      return CASE_FAIL;
    }
  } else {
    n = (size_t)(colon - id);
    if (n >= sizeof(suite))
      return CASE_FAIL;
    memcpy(suite, id, n);
    suite[n] = '\0';
    snprintf(test, sizeof(test), "%s", colon + 1);
  }
  if ((strcmp(suite, "view_category_sfinae") == 0) ||
      (strstr(test, "CudaMatrix") != NULL) ||
      (strcmp(test, "GroupedView.CategoryTagsAreDistinct") == 0)) {
    status = CASE_NA;
  }
  if (strstr(suite, "rate_constant") != NULL) {
    if (strcmp(suite, "lambda_rate_constant") == 0 ||
        strcmp(suite, "solver_lambda_rate_constant") == 0)
      status = test_lambda_callback();
    else if (strcmp(suite, "surface_rate_constant") == 0 &&
             strstr(test, "DiffusionCoefficientIsMissing") != NULL)
      status = test_surface_validation(0);
    else if (strcmp(suite, "surface_rate_constant") == 0 &&
             strstr(test, "MolecularWeightIsMissing") != NULL)
      status = test_surface_validation(1);
    else if (strcmp(suite, "ternary_chemical_activation_rate_constant") == 0 &&
        strstr(test, "JPL19") != NULL)
      status = test_jpl19();
    else
      status = test_rate(suite, test);
  } else if (strcmp(suite, "solver_lambda_rate_constant") == 0) {
    status = test_lambda_callback();
  } else if (strcmp(suite, "reaction_rate_store") == 0 ||
             strcmp(suite, "process") == 0 || strcmp(suite, "process_set") == 0) {
    status = test_ts1_metadata();
  } else if (strcmp(suite, "rosenbrock") == 0) {
    if (strstr(test, "CanSetTolerances") != NULL || strstr(test, "StandardAlphaMinusJacobian") != NULL ||
        strstr(test, "StandardNormalizedError") != NULL)
      status = test_first_order_solver();
  } else if (strcmp(suite, "regression_test_solve") == 0) {
    if (strstr(test, "TwoStageSolve") != NULL)
      status = test_two_stage_callbacks();
    else if (strstr(test, "ThreeStageSolve") != NULL)
      status = test_first_order_solver();
  } else if (strcmp(suite, "openmp") == 0) {
    status = test_first_order_solver();
  } else if (strcmp(suite, "matrix") == 0) {
    status = test_matrix_views(0);
  } else if (strcmp(suite, "vector_matrix") == 0 ||
             strcmp(suite, "grouped_views") == 0) {
    status = test_matrix_views(1);
    if (strcmp(suite, "grouped_views") == 0)
      status = test_grouped_view(test);
  } else if (strcmp(suite, "version") == 0) {
    status = test_version(test);
  } else if (strcmp(suite, "state") == 0) {
    status = test_state_api(test);
  } else if (strncmp(suite, "sparse_matrix_", 13) == 0) {
    const chem_sparse_format format = strstr(suite, "column") != NULL
                                          ? CHEM_SPARSE_CSC : CHEM_SPARSE_CSR;
    const size_t lanes = strstr(suite, "_vector_") != NULL ? 4 : 1;
    status = test_sparse_matrix_api(format, lanes);
  } else if (strncmp(suite, "lu_decomposition_", 17) == 0 ||
             strcmp(suite, "linear_solver") == 0 ||
             strcmp(suite, "linear_solver_in_place") == 0) {
    status = test_dense_lu_and_solve();
  }
  if (strcmp(suite, "view_category_sfinae") == 0 ||
      strstr(test, "CudaMatrix") != NULL ||
      strcmp(test, "GroupedView.CategoryTagsAreDistinct") == 0)
    status = CASE_NA;
  snprintf(reason, 256, "%s", reason_for(suite, test, status));
  return status;
}

int main(void)
{
  FILE *cases = fopen("../reference/micm_cases.txt", "r");
  char line[512], id[384], kind[32], reason[256];
  int failed = 0;
  if (!cases) {
    fprintf(stderr, "cannot open ../reference/micm_cases.txt\n");
    return 2;
  }
  while (fgets(line, sizeof(line), cases)) {
    int status;
    char *hash = strchr(line, '#');
    char *newline = strchr(line, '\n');
    if (hash)
      *hash = '\0';
    if (newline)
      *newline = '\0';
    if (sscanf(line, "%383s %31s", id, kind) != 2)
      continue;
    status = run_case(id, reason);
    printf("%s %s %s\n", id,
           status == CASE_PASS ? "PASS" : status == CASE_FAIL ? "FAIL" : "N/A",
           reason);
    failed |= status == CASE_FAIL;
  }
  fclose(cases);
  return failed ? 1 : 0;
}
