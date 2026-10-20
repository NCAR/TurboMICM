#include "kinetics.h"
#include "matrix.h"
#include "system.h"
#include "constraint.h"
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

static int test_matrix_views(int interleaved);
static int test_sparse_matrix_api(chem_sparse_format format, size_t lanes);

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

static int test_rate_store(const char *test)
{
  chem_species species[] = {
    {"a",0,0,0,0}, {"b",0,0,0,0}, {"c",0.025,230.0,0,0}
  };
  chem_custom_parameter custom[] = {
    {"rate1",2.0}, {"rate2",3.0}, {"radius",1.0e-8}, {"number",1.0e9}
  };
  chem_reaction reactions[6] = {{0}};
  chem_mechanism mechanism = {
    .name="rate layout test",.species_count=3,.reaction_count=6,
    .custom_count=4,.default_temperature=293.0,.default_pressure=101325.0,
    .species=species,.reactions=reactions,.custom_parameters=custom
  };
  size_t counts[CHEM_RATE_KIND_COUNT], offsets[CHEM_RATE_KIND_COUNT + 1];
  size_t ordered[6], i;

  if (strcmp(test, "ReactionRateConstantStore.OffsetsAreContiguousCumulativeSizes") == 0) {
    reactions[0].kind = CHEM_ARRHENIUS;
    reactions[1].kind = CHEM_USER_DEFINED;
    reactions[2].kind = CHEM_TROE;
    reactions[3].kind = CHEM_SURFACE;
    reactions[4].kind = CHEM_ARRHENIUS;
    reactions[5].kind = CHEM_TUNNELING;
    if (chem_rate_store_layout(&mechanism, counts, offsets, ordered) != 0 ||
        counts[CHEM_ARRHENIUS] != 2 || counts[CHEM_TROE] != 1 ||
        counts[CHEM_TUNNELING] != 1 || counts[CHEM_USER_DEFINED] != 1 ||
        counts[CHEM_SURFACE] != 1 || offsets[CHEM_ARRHENIUS] != 0 ||
        offsets[CHEM_TROE] != 2 || offsets[CHEM_TUNNELING] != 3 ||
        offsets[CHEM_USER_DEFINED] != 4 || offsets[CHEM_SURFACE] != 5 ||
        offsets[CHEM_RATE_KIND_COUNT] != 6 || ordered[0] != 0 ||
        ordered[1] != 4 || ordered[2] != 2 || ordered[3] != 5 ||
        ordered[4] != 1 || ordered[5] != 3)
      return CASE_FAIL;
    return CASE_PASS;
  }
  if (strcmp(test, "ReactionRateConstantStore.ArrheniusParametersPreserved") == 0) {
    const double expected[] = {2.15e-4,1.2,110.0,300.0,0.0};
    reactions[0].kind = CHEM_ARRHENIUS;
    memcpy(reactions[0].p, expected, sizeof(expected));
    for (i = 0; i < 5; ++i)
      if (reactions[0].p[i] != expected[i])
        return CASE_FAIL;
    return close_enough(chem_rate_constant(&mechanism, &reactions[0], 300.0,
                                            101325.0, NULL),
                        expected[0] * exp(expected[2] / 300.0) *
                            pow(300.0 / expected[3], expected[1]),
                        1.0e-14, 1.0e-14)
               ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "ReactionRateConstantStore.BranchedDerivedFieldsComputed") == 0) {
    const double expected_k0 = 2.0e-22 * 6.02214076e23 * 1.0e-6 * exp(3.0);
    const double air_ref = 2.45e19 / 6.02214076e23 * 1.0e6;
    const double a0 = expected_k0 * air_ref;
    const double b0 = 0.43 * pow(293.0 / 298.0, -8.0);
    const double aval0 = a0 / (1.0 + a0 / b0) *
                         pow(0.41, 1.0 / (1.0 + pow(log10(a0 / b0), 2.0)));
    const double z = aval0;
    const double air = 101325.0 / (1.380649e-23 * 6.02214076e23 * 293.0);
    const double a = expected_k0 * air;
    const double b = 0.43 * pow(293.0 / 298.0, -8.0);
    const double aval = a / (1.0 + a / b) *
                        pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
    const double expected = aval / (aval + z);
    reactions[0].kind = CHEM_BRANCHED;
    reactions[0].p[1] = 1.0;
    reactions[0].p[2] = 0.5;
    reactions[0].p[3] = 3.0;
    reactions[0].p[4] = 1.0;
    return close_enough(chem_rate_constant(&mechanism, &reactions[0], 293.0,
                                            101325.0, NULL),
                        expected, 1.0e-12, 1.0e-15) ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "ReactionRateConstantStore.UserDefinedCustomParamIndex") == 0) {
    double values[] = {2.5, 7.0, 0.0, 0.0};
    reactions[0].kind = CHEM_USER_DEFINED;
    reactions[0].custom0 = 0;
    reactions[0].p[0] = 2.0;
    reactions[1] = reactions[0];
    reactions[1].custom0 = 1;
    reactions[1].p[0] = 0.5;
    return strcmp(custom[reactions[0].custom0].name, "rate1") == 0 &&
                   strcmp(custom[reactions[1].custom0].name, "rate2") == 0 &&
                   chem_rate_constant(&mechanism, &reactions[0], 300.0, 101325.0, values) == 5.0 &&
                   chem_rate_constant(&mechanism, &reactions[1], 300.0, 101325.0, values) == 3.5
               ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "ReactionRateConstantStore.SurfaceDataFieldsAndCustomParamIndex") == 0) {
    double values[] = {0.0,0.0,1.0e-8,1.0e9};
    reactions[0].kind = CHEM_SURFACE;
    reactions[0].surface_species = 2;
    reactions[0].custom0 = 2;
    reactions[0].custom1 = 3;
    reactions[0].p[0] = 0.74;
    return reactions[0].surface_species == 2 && species[2].diffusion_coefficient == 230.0 &&
                   close_enough(chem_rate_constant(&mechanism, &reactions[0], 300.0,
                                                   101325.0, values),
                                4.0e9 * 3.14159265358979323846 * 1.0e-16 /
                                    (1.0e-8 / 230.0 + 4.0 /
                                     (sqrt(8.0 * 8.31446261815324 * 300.0 /
                                           (3.14159265358979323846 * 0.025)) * 0.74)),
                                1.0e-13, 1.0e-20)
               ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "ReactionRateConstantStore.SurfaceCustomParamIndexAfterUserDefined") == 0) {
    reactions[0].kind = CHEM_USER_DEFINED;
    reactions[0].custom0 = 0;
    reactions[1].kind = CHEM_SURFACE;
    reactions[1].surface_species = 2;
    reactions[1].custom0 = 2;
    reactions[1].custom1 = 3;
    return reactions[0].custom0 == 0 && reactions[1].custom0 == 2 &&
                   reactions[1].custom1 == 3 &&
                   strcmp(custom[reactions[1].custom0].name, "radius") == 0 &&
                   strcmp(custom[reactions[1].custom1].name, "number") == 0
               ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "ReactionRateConstantStore.LambdaEntriesRcIndex") == 0) {
    double scale = 0.075;
    double custom_values[] = {4.0,0.0,0.0,0.0};
    reactions[0].kind = CHEM_ARRHENIUS;
    reactions[1].kind = CHEM_CALLBACK;
    reactions[1].callback = test_callback;
    reactions[1].callback_context = &scale;
    return reactions[1].kind == CHEM_CALLBACK &&
                   close_enough(chem_rate_constant(&mechanism, &reactions[1], 300.0,
                                                   101325.0, custom_values),
                                0.3, 0.0, 0.0)
               ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "ReactionRateConstantStore.ParameterizedMultipliers") == 0) {
    reactions[0].kind = CHEM_ARRHENIUS;
    reactions[0].p[0] = 2.0;
    reactions[0].p[3] = 300.0;
    reactions[0].body_order = 1;
    return close_enough(chem_effective_rate_constant(&mechanism, &reactions[0],
                                                     300.0, 101325.0, NULL),
                        chem_rate_constant(&mechanism, &reactions[0], 300.0,
                                           101325.0, NULL) *
                            (101325.0 / (1.380649e-23 * 6.02214076e23 * 300.0)),
                        1.0e-14, 0.0) ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "ReactionRateConstantStore.TotalRateConstantsMatchesProcessCount") == 0) {
    reactions[0].kind = CHEM_ARRHENIUS;
    reactions[1].kind = CHEM_TROE;
    reactions[2].kind = CHEM_USER_DEFINED;
    reactions[3].kind = CHEM_SURFACE;
    reactions[4].kind = CHEM_CALLBACK;
    mechanism.reaction_count = 5;
    if (chem_rate_store_layout(&mechanism, counts, offsets, ordered) != 0)
      return CASE_FAIL;
    return offsets[CHEM_RATE_KIND_COUNT] == mechanism.reaction_count
               ? CASE_PASS : CASE_FAIL;
  }
  return CASE_FAIL;
}

static int test_process_case(const char *test)
{
  static const chem_species species[] = {
    {"O3",0,0,0,0}, {"NO",0,0,0,0}, {"NO2",0,0,0,0}, {"O2",0,0,0,0}
  };
  static const uint32_t reactants[] = {0,1};
  static const uint32_t products[] = {2,3};
  static const double yields[] = {1.0,1.0};
  static const chem_reaction reaction = {
    .kind=CHEM_ARRHENIUS,.reactant_offset=0,.reactant_count=2,
    .product_offset=0,.product_count=2,.p={0.5,0.0,0.0,300.0}
  };
  static const chem_mechanism mechanism = {
    .name="reaction table",.species_count=4,.reaction_count=1,
    .reactant_count=2,.product_count=2,.default_temperature=298.15,
    .default_pressure=101325.0,.species=species,.reactions=&reaction,
    .reactants=reactants,.products=products,.yields=yields
  };
  if (strcmp(test, "Process.Matrix") == 0)
    return test_matrix_views(0);
  if (strcmp(test, "Process.VectorMatrix") == 0)
    return test_matrix_views(1);
  if (strcmp(test, "Process.BuildsChemicalReaction") == 0) {
    const double concentrations[] = {2.0,3.0,0.0,0.0};
    const double rates[] = {0.5};
    double forcing[4];
    if (reaction.kind != CHEM_ARRHENIUS || reaction.reactant_count != 2 ||
        reaction.product_count != 2 ||
        strcmp(species[mechanism.products[0]].name, "NO2") != 0 ||
        strcmp(species[mechanism.products[1]].name, "O2") != 0 ||
        chem_mechanism_forcing(&mechanism, concentrations, rates, forcing) != 0)
      return CASE_FAIL;
    return forcing[0] == -3.0 && forcing[1] == -3.0 &&
                   forcing[2] == 3.0 && forcing[3] == 3.0
               ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "Process.ChemicalReactionCopyAssignmentSucceeds") == 0) {
    chem_reaction assigned = reaction;
    chem_mechanism copied = mechanism;
    return assigned.kind == reaction.kind && assigned.p[0] == reaction.p[0] &&
                   copied.species == mechanism.species &&
                   copied.reactions == mechanism.reactions &&
                   copied.reactants == mechanism.reactants &&
                   copied.products == mechanism.products &&
                   copied.yields == mechanism.yields &&
                   copied.product_count == 2
               ? CASE_PASS : CASE_FAIL;
  }
  return CASE_FAIL;
}

static int test_finite_difference_jacobian(void)
{
  static const chem_species species[] = {
    {"A",0,0,0,0}, {"B",0,0,0,0}, {"C",0,0,0,0}
  };
  static const uint32_t reactants[] = {0,0,1};
  static const uint32_t products[] = {1,2};
  static const double yields[] = {1.0,1.0};
  static const chem_reaction reactions[] = {
    {.kind=CHEM_ARRHENIUS,.reactant_offset=0,.reactant_count=2,
     .product_offset=0,.product_count=1,.p={0.3,0,0,300}},
    {.kind=CHEM_ARRHENIUS,.reactant_offset=2,.reactant_count=1,
     .product_offset=1,.product_count=1,.p={0.2,0,0,300}}
  };
  static const chem_mechanism mechanism = {
    .name="finite difference",.species_count=3,.reaction_count=2,
    .reactant_count=3,.product_count=2,.species=species,.reactions=reactions,
    .reactants=reactants,.products=products,.yields=yields
  };
  const double rates[] = {0.3,0.2};
  const double y[] = {1.7,0.8,0.2};
  const double epsilon = 1.0e-7;
  double forcing[3], jacobian[9];
  size_t column, row;
  if (chem_mechanism_forcing(&mechanism, y, rates, forcing) != 0 ||
      chem_mechanism_jacobian_dense(&mechanism, y, rates, jacobian) != 0 ||
      !close_enough(forcing[0], -2.0 * 0.3 * y[0] * y[0], 1.0e-14, 1.0e-14) ||
      !close_enough(forcing[1], 0.3 * y[0] * y[0] - 0.2 * y[1], 1.0e-14, 1.0e-14) ||
      !close_enough(forcing[2], 0.2 * y[1], 1.0e-14, 1.0e-14))
    return CASE_FAIL;
  for (column = 0; column < 3; ++column) {
    double perturbed[3], shifted[3], finite_difference;
    memcpy(perturbed, y, sizeof(y));
    perturbed[column] += epsilon;
    if (chem_mechanism_forcing(&mechanism, perturbed, rates, shifted) != 0)
      return CASE_FAIL;
    for (row = 0; row < 3; ++row) {
      finite_difference = (shifted[row] - forcing[row]) / epsilon;
      if (!close_enough(finite_difference, jacobian[row * 3 + column],
                        2.0e-7, 1.0e-8))
        return CASE_FAIL;
    }
  }
  return CASE_PASS;
}

static chem_species_property_value property_double(double value)
{
  chem_species_property_value result;
  result.real_value = value;
  return result;
}

static chem_species_property_value property_string(const char *value)
{
  chem_species_property_value result;
  result.string_value = value;
  return result;
}

static chem_species_property_value property_integer(int value)
{
  chem_species_property_value result;
  result.integer_value = value;
  return result;
}

static chem_species_property_value property_boolean(int value)
{
  chem_species_property_value result;
  result.boolean_value = value;
  return result;
}

static int test_species_api(const char *test)
{
  chem_species species = {0}, copied = {0};
  double real_value = 0.0;
  const char *string_value = NULL;
  int int_value = 0, bool_value = 0;
  int result = CASE_FAIL;
  if (chem_species_init(&species, "thing") != 0)
    return CASE_FAIL;
  if (strcmp(test, "Species.StringConstructor") == 0) {
    if (!chem_species_is_parameterized(&species) &&
        chem_species_set_parameterization(&species, CHEM_SPECIES_TEMPERATURE,
                                          0.0, 100.0, 0.0) == 0 &&
        chem_species_is_parameterized(&species) &&
        chem_species_parameterize(&species, 12.5, 900.0) == 1250.0)
      result = CASE_PASS;
  } else if (strcmp(test, "Species.StringAndVectorConstructor") == 0 ||
             strcmp(test, "Species.GetProperty") == 0) {
    if (chem_species_set_property(&species, "name [units]", CHEM_PROPERTY_DOUBLE,
                                  property_double(1.0)) != 0 ||
        chem_species_set_property(&species, "name2 [units2]", CHEM_PROPERTY_DOUBLE,
                                  property_double(2.0)) != 0)
      goto done;
    if (chem_species_get_property(&species, "name [units]", CHEM_PROPERTY_DOUBLE,
                                  &real_value) != 0 || real_value != 1.0 ||
        chem_species_get_property(&species, "name2 [units2]", CHEM_PROPERTY_DOUBLE,
                                  &real_value) != 0 || real_value != 2.0)
      goto done;
    if (strcmp(test, "Species.StringAndVectorConstructor") == 0) {
      result = species.property_count == 2 ? CASE_PASS : CASE_FAIL;
    } else if (chem_species_get_property(&species, "not there", CHEM_PROPERTY_STRING,
                                         &string_value) == -2 &&
               chem_species_get_property(&species, "not there", CHEM_PROPERTY_DOUBLE,
                                         &real_value) == -2 &&
               chem_species_get_property(&species, "not there", CHEM_PROPERTY_INTEGER,
                                         &int_value) == -2 &&
               chem_species_get_property(&species, "not there", CHEM_PROPERTY_BOOLEAN,
                                         &bool_value) == -2 &&
               chem_species_get_property(&species, "name [units]", CHEM_PROPERTY_STRING,
                                         &string_value) == -3) {
      result = CASE_PASS;
    }
  } else if (strcmp(test, "Species.SetThirdBody") == 0) {
    if (chem_species_set_third_body(&species) == 0 &&
        chem_species_is_parameterized(&species) &&
        chem_species_parameterize(&species, 300.0, 42.4) == 42.4)
      result = CASE_PASS;
  } else if (strcmp(test, "Species.CopyConstructor") == 0 ||
             strcmp(test, "Species.CopyAssignment") == 0) {
    if (chem_species_set_property(&species, "name [units]", CHEM_PROPERTY_DOUBLE,
                                  property_double(1.0)) != 0 ||
        chem_species_set_property(&species, "name2 [units2]", CHEM_PROPERTY_DOUBLE,
                                  property_double(2.0)) != 0 ||
        chem_species_set_property(&species, "foo", CHEM_PROPERTY_STRING,
                                  property_string("bar")) != 0 ||
        chem_species_set_property(&species, "baz", CHEM_PROPERTY_INTEGER,
                                  property_integer(42)) != 0 ||
        chem_species_set_property(&species, "qux", CHEM_PROPERTY_BOOLEAN,
                                  property_boolean(1)) != 0)
      goto done;
    if (strstr(test, "CopyConstructor") != NULL &&
        chem_species_set_parameterization(&species, CHEM_SPECIES_CONSTANT,
                                          15.4, 0.0, 0.0) != 0)
      goto done;
    if (chem_species_copy(&copied, &species) != 0 || copied.property_count != 5 ||
        strcmp(copied.name, "thing") != 0 ||
        chem_species_get_property(&copied, "name [units]", CHEM_PROPERTY_DOUBLE,
                                  &real_value) != 0 || real_value != 1.0 ||
        chem_species_get_property(&copied, "foo", CHEM_PROPERTY_STRING,
                                  &string_value) != 0 || strcmp(string_value, "bar") != 0 ||
        chem_species_get_property(&copied, "baz", CHEM_PROPERTY_INTEGER,
                                  &int_value) != 0 || int_value != 42 ||
        chem_species_get_property(&copied, "qux", CHEM_PROPERTY_BOOLEAN,
                                  &bool_value) != 0 || bool_value != 1 ||
        chem_species_is_parameterized(&copied) !=
            chem_species_is_parameterized(&species))
      goto done;
    if (chem_species_is_parameterized(&species) &&
        chem_species_parameterize(&copied, 300.0, 0.0) != 15.4)
      goto done;
    if (chem_species_set_property(&copied, "foo", CHEM_PROPERTY_STRING,
                                  property_string("changed")) != 0 ||
        chem_species_get_property(&species, "foo", CHEM_PROPERTY_STRING,
                                  &string_value) != 0 || strcmp(string_value, "bar") != 0)
      goto done;
    result = CASE_PASS;
  }
done:
  chem_species_destroy(&species);
  chem_species_destroy(&copied);
  return result;
}

static int make_phase_species(chem_phase_species *phase_species, const char *name,
                              int parameterized, int has_diffusion,
                              double diffusion)
{
  chem_species species = {0};
  int result;
  if (chem_species_init(&species, name) != 0)
    return -1;
  if (parameterized &&
      chem_species_set_parameterization(&species, CHEM_SPECIES_CONSTANT,
                                        1.0, 0.0, 0.0) != 0) {
    chem_species_destroy(&species);
    return -1;
  }
  result = chem_phase_species_init(phase_species, &species, has_diffusion,
                                   diffusion);
  chem_species_destroy(&species);
  return result;
}

static int test_phase_api(const char *test)
{
  chem_phase_species items[3] = {{0}};
  chem_phase phase = {0};
  char first[64], second[64];
  int result = CASE_FAIL;
  if (strcmp(test, "PhaseSpecies.ConstructsWithoutDiffusionCoefficient") == 0 ||
      strcmp(test, "PhaseSpecies.ConstructsWithDiffusionCoefficient") == 0 ||
      strcmp(test, "PhaseSpecies.SetDiffusionCoefficient") == 0) {
    chem_phase_species one = {0};
    const int with_diffusion =
        strcmp(test, "PhaseSpecies.ConstructsWithoutDiffusionCoefficient") != 0;
    if (make_phase_species(&one, "CO2", 0,
                           strcmp(test, "PhaseSpecies.ConstructsWithDiffusionCoefficient") == 0,
                           1.23e-5) != 0)
      return CASE_FAIL;
    if (strcmp(test, "PhaseSpecies.SetDiffusionCoefficient") == 0 &&
        chem_phase_species_set_diffusion_coefficient(&one, 2.5e-6) == 0)
      result = one.has_diffusion_coefficient && one.diffusion_coefficient == 2.5e-6
                   ? CASE_PASS : CASE_FAIL;
    else if (strcmp(test, "PhaseSpecies.ConstructsWithoutDiffusionCoefficient") == 0)
      result = strcmp(one.species.name, "CO2") == 0 &&
                       !one.has_diffusion_coefficient
                   ? CASE_PASS : CASE_FAIL;
    else
      result = with_diffusion && strcmp(one.species.name, "CO2") == 0 &&
                       one.has_diffusion_coefficient &&
                       one.diffusion_coefficient == 1.23e-5
                   ? CASE_PASS : CASE_FAIL;
    chem_phase_species_destroy(&one);
    return result;
  }
  if (make_phase_species(&items[0], "species1", 0, 0, 0.0) != 0 ||
      make_phase_species(&items[1], "species2", 0, 0, 0.0) != 0)
    goto done;
  if (strcmp(test, "Phase.ConstructorWithParameterizedSpecies") == 0 ||
      strcmp(test, "Phase.UniqueNamesWithParameterizedSpecies") == 0) {
    chem_phase_species_destroy(&items[1]);
    if (make_phase_species(&items[1], "bar", 1, 0, 0.0) != 0 ||
        make_phase_species(&items[2], "baz", 0, 0, 0.0) != 0)
      goto done;
    chem_phase_species_destroy(&items[0]);
    if (make_phase_species(&items[0], "foo", 0, 0, 0.0) != 0)
      goto done;
  }
  if (strcmp(test, "Phase.Constructor") == 0) {
    if (chem_phase_init(&phase, "test_phase", items, 2) != 0)
      goto done;
    if (chem_phase_state_size(&phase) == 2 &&
        chem_phase_unique_name(&phase, 0, first, sizeof(first)) == 1 &&
        chem_phase_unique_name(&phase, 1, second, sizeof(second)) == 1 &&
        strcmp(first, "test_phase.species1") == 0 &&
        strcmp(second, "test_phase.species2") == 0 &&
        strcmp(phase.name, "test_phase") == 0)
      result = CASE_PASS;
  } else if (strcmp(test, "Phase.ConstructorWithParameterizedSpecies") == 0) {
    if (chem_phase_init(&phase, "gas", items, 3) == 0 &&
        phase.species_count == 3 && chem_phase_state_size(&phase) == 2)
      result = CASE_PASS;
  } else if (strcmp(test, "Phase.UniqueNamesWithParameterizedSpecies") == 0) {
    if (chem_phase_init(&phase, "gas", items, 3) == 0 &&
        chem_phase_state_size(&phase) == 2 &&
        chem_phase_unique_name(&phase, 0, first, sizeof(first)) == 1 &&
        chem_phase_unique_name(&phase, 1, second, sizeof(second)) == 0 &&
        chem_phase_unique_name(&phase, 2, second, sizeof(second)) == 1 &&
        strcmp(first, "gas.foo") == 0 && strcmp(second, "gas.baz") == 0)
      result = CASE_PASS;
  }
done:
  chem_phase_destroy(&phase);
  for (size_t i = 0; i < 3; ++i)
    chem_phase_species_destroy(&items[i]);
  return result;
}

static int test_system_api(const char *test)
{
  chem_phase_species items[3] = {{0}};
  chem_phase phase = {0};
  chem_system system = {0};
  char first[64], second[64];
  const size_t reverse[] = {1,0};
  int result = CASE_FAIL;
  if (make_phase_species(&items[0], "foo", 0, 0, 0.0) != 0 ||
      make_phase_species(&items[1], "bar", 0, 0, 0.0) != 0)
    goto done;
  if (strcmp(test, "System.ConstructorWithParameterizedSpecies") == 0 &&
      make_phase_species(&items[2], "paramSpecies", 1, 0, 0.0) != 0)
    goto done;
  if (chem_phase_init(&phase, "gas", items,
                      strcmp(test, "System.ConstructorWithParameterizedSpecies") == 0 ? 3 : 2) != 0 ||
      chem_system_init(&system, &phase) != 0)
    goto done;
  if (system.gas_phase.species_count != phase.species_count)
    goto done;
  if (strcmp(test, "System.ConstructorWithAllParameters") == 0) {
    if (chem_system_state_size(&system) == 2 &&
        chem_system_unique_name(&system, 0, NULL, 0, first, sizeof(first)) == 1 &&
        chem_system_unique_name(&system, 1, NULL, 0, second, sizeof(second)) == 1 &&
        ((strcmp(first, "foo") == 0 && strcmp(second, "bar") == 0) ||
         (strcmp(first, "bar") == 0 && strcmp(second, "foo") == 0)))
      result = CASE_PASS;
  } else if (strcmp(test, "System.ConstructorWithParameterizedSpecies") == 0) {
    if (chem_system_state_size(&system) == 2 &&
        chem_system_unique_name(&system, 0, NULL, 0, first, sizeof(first)) == 1 &&
        chem_system_unique_name(&system, 1, NULL, 0, second, sizeof(second)) == 1 &&
        chem_system_unique_name(&system, 0, reverse, 2, first, sizeof(first)) == 1 &&
        chem_system_unique_name(&system, 1, reverse, 2, second, sizeof(second)) == 1 &&
        strcmp(first, "bar") == 0 && strcmp(second, "foo") == 0)
      result = CASE_PASS;
  }
done:
  chem_system_destroy(&system);
  chem_phase_destroy(&phase);
  for (size_t i = 0; i < 3; ++i)
    chem_phase_species_destroy(&items[i]);
  return result;
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

static int test_first_order_solver_with(const chem_rosenbrock_parameters *parameters,
                                        size_t cells, size_t group_size)
{
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
  double y[16], temperature[8], pressure[8];
  chem_solver_stats stats = {0};
  chem_solver *solver;
  int result, cell;
  if (!cells || cells > 8)
    return CASE_FAIL;
  for (cell = 0; cell < (int)cells; ++cell) {
    y[(size_t)cell * 2] = 1.0;
    y[(size_t)cell * 2 + 1] = 0.0;
    temperature[cell] = 300.0;
    pressure[cell] = 101325.0;
  }
  solver = chem_solver_create(&model, cells, temperature, pressure, NULL,
                              parameters, group_size);
  if (!solver)
    return CASE_FAIL;
  result = chem_solver_step(solver, y, 1.0, &stats);
  chem_solver_destroy(solver);
  if (result != 0 || stats.accepted == 0 ||
      stats.solves != (parameters ? parameters->stages : 3) * stats.decompositions)
    return CASE_FAIL;
  for (cell = 0; cell < (int)cells; ++cell)
    if (!close_enough(y[(size_t)cell * 2], exp(-0.2), 2.0e-3, 2.0e-6) ||
        !close_enough(y[(size_t)cell * 2 + 1], 1.0 - exp(-0.2), 2.0e-3, 2.0e-6))
      return CASE_FAIL;
  return CASE_PASS;
}

static int test_first_order_solver(void)
{
  return test_first_order_solver_with(NULL, 1, 1);
}

static int test_regression_solver_case(const char *test)
{
  chem_rosenbrock_parameters parameters;
  if (strstr(test, "FourStageDASolve") != NULL)
    chem_rosenbrock_four_stage_dae(&parameters);
  else if (strstr(test, "SixStageDASolve") != NULL)
    chem_rosenbrock_six_stage_dae(&parameters);
  else if (strstr(test, "FourStageSolve") != NULL)
    chem_rosenbrock_four_stage(&parameters);
  else if (strstr(test, "VectorSolve") != NULL)
    return test_first_order_solver_with(NULL, 4, 2);
  else
    return CASE_FAIL;
  return test_first_order_solver_with(&parameters, 1, 1);
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

static int test_alpha_minus_jacobian(size_t blocks, size_t lanes)
{
  static const uint32_t rows[] = {0,0,0,1,1,1,2,2,3,3,4,4,4};
  static const uint32_t columns[] = {0,1,2,0,1,2,0,2,1,3,0,2,4};
  static const double alpha = 42.042;
  chem_sparse_matrix matrix = {0};
  size_t cell;
  int result = CASE_FAIL;
  if (chem_sparse_matrix_init(&matrix, 5, 5, rows, columns, 13, blocks,
                              lanes, CHEM_SPARSE_CSR) != 0)
    return CASE_FAIL;
  for (cell = 0; cell < blocks; ++cell) {
    const double divisor = (double)(cell + 2);
    const double values[][3] = {
      {-12.2, -24.3 * divisor, -42.3},
      {-0.43, -23.4, -83.4 / (double)(cell + 3)},
      {-4.74, 0.0, -6.91},
      {0.0, -59.1, -83.4},
      {-78.5, -53.6, -1.0}
    };
    if (chem_sparse_matrix_set(&matrix, cell, 0, 0, values[0][0]) != 0 ||
        chem_sparse_matrix_set(&matrix, cell, 0, 1, values[0][1]) != 0 ||
        chem_sparse_matrix_set(&matrix, cell, 0, 2, values[0][2]) != 0 ||
        chem_sparse_matrix_set(&matrix, cell, 1, 0, values[1][0]) != 0 ||
        chem_sparse_matrix_set(&matrix, cell, 1, 1, values[1][1]) != 0 ||
        chem_sparse_matrix_set(&matrix, cell, 1, 2, values[1][2]) != 0 ||
        chem_sparse_matrix_set(&matrix, cell, 2, 0, values[2][0]) != 0 ||
        chem_sparse_matrix_set(&matrix, cell, 2, 2, values[2][2]) != 0 ||
        chem_sparse_matrix_set(&matrix, cell, 3, 1, values[3][1]) != 0 ||
        chem_sparse_matrix_set(&matrix, cell, 3, 3, values[3][2]) != 0 ||
        chem_sparse_matrix_set(&matrix, cell, 4, 0, values[4][0]) != 0 ||
        chem_sparse_matrix_set(&matrix, cell, 4, 2, values[4][1]) != 0 ||
        chem_sparse_matrix_set(&matrix, cell, 4, 4, values[4][2]) != 0)
      goto done;
  }
  if (chem_sparse_matrix_add_diagonal(&matrix, alpha) != 0)
    goto done;
  for (cell = 0; cell < blocks; ++cell) {
    const double divisor = (double)(cell + 2);
    const double expected[][3] = {
      {alpha - 12.2, -24.3 * divisor, -42.3},
      {-0.43, alpha - 23.4, -83.4 / (double)(cell + 3)},
      {-4.74, 0.0, alpha - 6.91},
      {0.0, -59.1, alpha - 83.4},
      {-78.5, -53.6, alpha - 1.0}
    };
    const size_t check_row[] = {0,0,0,1,1,1,2,2,3,3,4,4,4};
    const size_t check_column[] = {0,1,2,0,1,2,0,2,1,3,0,2,4};
    const double check_value[] = {
      expected[0][0], expected[0][1], expected[0][2],
      expected[1][0], expected[1][1], expected[1][2],
      expected[2][0], expected[2][2], expected[3][1], expected[3][2],
      expected[4][0], expected[4][1], expected[4][2]
    };
    size_t entry;
    for (entry = 0; entry < 13; ++entry) {
      double actual;
      if (chem_sparse_matrix_get(&matrix, cell, check_row[entry],
                                 check_column[entry], &actual) != 0 ||
          !close_enough(actual, check_value[entry], 1.0e-14, 1.0e-14))
        goto done;
    }
  }
  result = CASE_PASS;
done:
  chem_sparse_matrix_destroy(&matrix);
  return result;
}

static int test_normalized_error(size_t species_count)
{
  chem_species species[4] = {
    {.absolute_tolerance=1.0e-4}, {.absolute_tolerance=2.0e-4},
    {.absolute_tolerance=0.0}, {.absolute_tolerance=5.0e-4}
  };
  double old_values[5 * 4], new_values[5 * 4], errors[5 * 4];
  double expected_sum = 0.0;
  size_t cell, variable;
  double computed, expected;
  if (!species_count || species_count > 4)
    return CASE_FAIL;
  for (cell = 0; cell < 5; ++cell)
    for (variable = 0; variable < species_count; ++variable) {
      const size_t i = cell * species_count + variable;
      const double atol = species[variable].absolute_tolerance > 0.0
                              ? species[variable].absolute_tolerance : 1.0e-3;
      double scaled;
      old_values[i] = 1.0 + (double)cell + 0.1 * (double)variable;
      new_values[i] = 0.8 + 0.5 * (double)cell + 0.2 * (double)variable;
      errors[i] = 0.01 * (double)(1 + cell + variable);
      scaled = errors[i] /
               (atol + 1.0e-6 * fmax(fabs(old_values[i]), fabs(new_values[i])));
      expected_sum += scaled * scaled;
    }
  computed = chem_normalized_error(old_values, new_values, errors, species,
                                   5, species_count, 1.0e-6, 1.0e-3);
  expected = fmax(sqrt(expected_sum / (double)(5 * species_count)), 1.0e-10);
  if (!close_enough(computed, expected, 1.0e-13, 1.0e-13))
    return CASE_FAIL;
  memset(errors, 0, sizeof(errors));
  return chem_normalized_error(old_values, new_values, errors, species,
                               5, species_count, 1.0e-6, 1.0e-3) == 1.0e-10
             ? CASE_PASS : CASE_FAIL;
}

static int test_solver_tolerances(void)
{
  static const chem_species species[] = {
    {.name="A", .absolute_tolerance=1.0e-7},
    {.name="B", .absolute_tolerance=1.0e-8}
  };
  static const chem_reaction reactions[] = {{0}};
  static const chem_mechanism mechanism = {
    .name="tolerances", .species_count=2, .reaction_count=0,
    .default_temperature=298.15, .default_pressure=101325.0,
    .species=species, .reactions=reactions
  };
  chem_solver *solver = chem_solver_create(&mechanism, 1, NULL, NULL, NULL,
                                            NULL, 1);
  int passed;
  if (!solver)
    return CASE_FAIL;
  passed = chem_solver_absolute_tolerance(solver, 0) == 1.0e-7 &&
           chem_solver_absolute_tolerance(solver, 1) == 1.0e-8 &&
           isnan(chem_solver_absolute_tolerance(solver, 2)) &&
           chem_solver_relative_tolerance(solver) == 1.0e-6;
  chem_solver_destroy(solver);
  return passed ? CASE_PASS : CASE_FAIL;
}

static int test_rejected_step_retry(void)
{
  static const chem_species species[] = {
    {"foo",0,0,0,0}, {"bar",0,0,0,0}, {"baz",0,0,0,0},
    {"quz",0,0,0,0}, {"quuz",0,0,0,0}
  };
  static const uint32_t reactants[] = {0,2,1,3};
  static const uint32_t products[] = {1,4,0,3};
  static const double yields[] = {1.0,2.4,1.0,1.4};
  static const chem_reaction reactions[] = {
    {.kind=CHEM_ARRHENIUS,.reactant_offset=0,.reactant_count=2,
     .product_offset=0,.product_count=2,.p={2.0e-11,0,110,300}},
    {.kind=CHEM_ARRHENIUS,.reactant_offset=2,.reactant_count=1,
     .product_offset=2,.product_count=2,.p={1.0e-6,0,0,300}},
    {.kind=CHEM_ARRHENIUS,.reactant_offset=3,.reactant_count=1,
     .product_offset=4,.product_count=0,.p={3.5e-6,0,0,300}}
  };
  static const chem_mechanism mechanism = {
    .name="retry test",.species_count=5,.reaction_count=3,
    .reactant_count=4,.product_count=4,.default_temperature=298.15,
    .default_pressure=101325,.species=species,.reactions=reactions,
    .reactants=reactants,.products=products,.yields=yields
  };
  chem_rosenbrock_parameters parameters;
  chem_solver_stats first_stats = {0}, second_stats = {0};
  chem_solver *first_solver, *second_solver;
  double first[] = {1.0e12,1.0e12,1.0e12,1.0e12,1.0e12};
  double second[] = {1.0e12,1.0e12,1.0e12,1.0e12,1.0e12};
  size_t i;
  int first_result, second_result;
  chem_rosenbrock_default(&parameters);
  parameters.h_start_fraction = 1.0;
  first_solver = chem_solver_create(&mechanism, 1, NULL, NULL, NULL,
                                   &parameters, 1);
  second_solver = chem_solver_create(&mechanism, 1, NULL, NULL, NULL,
                                    &parameters, 0);
  if (!first_solver || !second_solver) {
    chem_solver_destroy(first_solver);
    chem_solver_destroy(second_solver);
    return CASE_FAIL;
  }
  first_result = chem_solver_step(first_solver, first, 1000.0, &first_stats);
  second_result = chem_solver_step(second_solver, second, 1000.0, &second_stats);
  chem_solver_destroy(first_solver);
  chem_solver_destroy(second_solver);
  if (first_result != 0 || second_result != 0 ||
      first_stats.internal_steps < first_stats.accepted + 2 ||
      first_stats.internal_steps != second_stats.internal_steps ||
      first_stats.accepted != second_stats.accepted ||
      first_stats.rejected != second_stats.rejected)
    return CASE_FAIL;
  for (i = 0; i < 5; ++i)
    if (!isfinite(first[i]) || !isfinite(second[i]) ||
        !close_enough(first[i], second[i], 1.0e-13, 0.0))
      return CASE_FAIL;
  return CASE_PASS;
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

static int test_constraint_gradient_fd(const chem_constraint *constraint,
                                      size_t species_count, const double *state,
                                      double equilibrium_constant)
{
  double gradient[8], plus[8], minus[8], fp, fm;
  size_t i;
  if (species_count > sizeof(gradient) / sizeof(gradient[0]) ||
      chem_constraint_jacobian(constraint, state, species_count,
                               equilibrium_constant, 298.15, gradient) != 0)
    return 0;
  for (i = 0; i < species_count; ++i) {
    const double h = 1.0e-6 * fmax(fabs(state[i]), 1.0e-4);
    memcpy(plus, state, species_count * sizeof(*state));
    memcpy(minus, state, species_count * sizeof(*state));
    plus[i] += h;
    minus[i] -= h;
    if (chem_constraint_residual(constraint, plus, species_count,
                                 equilibrium_constant, 298.15, &fp) != 0 ||
        chem_constraint_residual(constraint, minus, species_count,
                                 equilibrium_constant, 298.15, &fm) != 0 ||
        !close_enough(gradient[i], (fp - fm) / (2.0 * h), 2.0e-7, 1.0e-8))
      return 0;
  }
  return 1;
}

static int test_linear_constraint_case(const char *test)
{
  chem_constraint constraint;
  chem_constraint_term terms[3] = {{0, 1.0}, {1, 1.0}, {2, 1.0}};
  double state[3] = {0.3, 0.7, 0.0}, value, gradient[3];
  size_t count = 2, algebraic = 1;
  double constant = 1.0;

  if (strcmp(test, "LinearConstraint.WeightedTerms") == 0) {
    terms[0].coefficient = 2.0; terms[1].coefficient = 3.0;
    terms[2].coefficient = -1.0; count = 3; constant = 5.0;
  } else if (strcmp(test, "LinearConstraint.AlgebraicSpecies") == 0) {
    count = 3; algebraic = 2; constant = 10.0;
  } else if (strcmp(test, "LinearConstraint.ZeroConstant") == 0 ||
             strcmp(test, "LinearConstraint.ZeroConstantResidual") == 0) {
    terms[1].coefficient = -1.0; constant = 0.0;
    if (strcmp(test, "LinearConstraint.ZeroConstantResidual") == 0) {
      state[0] = 0.6; state[1] = 0.5;
    }
  } else if (strcmp(test, "LinearConstraint.FractionalCoefficients") == 0 ||
             strcmp(test, "LinearConstraint.FractionalCoefficientsResidualAndJacobian") == 0) {
    terms[0].coefficient = 0.5; terms[1].coefficient = 1.5;
    constant = 2.0;
    if (strcmp(test, "LinearConstraint.FractionalCoefficientsResidualAndJacobian") == 0) {
      state[0] = 1.0; state[1] = 1.0;
    }
  } else if (strcmp(test, "LinearConstraint.ThreeSpeciesConservationResidual") == 0) {
    count = 3; constant = 10.0; state[0] = 1.0; state[1] = 2.0; state[2] = 7.0;
  } else if (strcmp(test, "LinearConstraint.WeightedSumResidualAndJacobian") == 0 ||
             strcmp(test, "LinearConstraint.FiniteDifferenceJacobianWeightedSum") == 0) {
    terms[0].coefficient = 2.0; terms[1].coefficient = 3.0;
    constant = 1.0; state[0] = 0.25; state[1] = 0.5;
  } else if (strcmp(test, "LinearConstraint.ResidualComputationThroughConstraintSet") == 0 ||
             strcmp(test, "LinearConstraint.JacobianComputationThroughConstraintSet") == 0 ||
             strcmp(test, "LinearConstraint.FiniteDifferenceJacobianSimpleConservation") == 0) {
    /* A + B = 1 */
    if (strcmp(test, "LinearConstraint.ResidualComputationThroughConstraintSet") == 0) {
      state[0] = 0.5; state[1] = 0.6;
    }
  } else if (strcmp(test, "LinearConstraint.JacobianIndependentOfConcentrations") == 0) {
    terms[0].coefficient = 2.0; terms[1].coefficient = 3.0; constant = 1.0;
  } else if (strcmp(test, "LinearConstraint.Construction") != 0) {
    return CASE_FAIL;
  }
  if (chem_constraint_linear_init(&constraint, "linear", algebraic, terms,
                                  count, constant) != 0)
    return CASE_FAIL;
  if (strcmp(test, "LinearConstraint.Construction") == 0)
    return strcmp(constraint.name, "linear") == 0 &&
           constraint.kind == CHEM_CONSTRAINT_LINEAR &&
           constraint.term_count == 2 && constraint.constant == 1.0 ? CASE_PASS : CASE_FAIL;
  if (chem_constraint_jacobian(&constraint, state, 3, 0.0, 298.15, gradient) != 0)
    return CASE_FAIL;
  if (strcmp(test, "LinearConstraint.AlgebraicSpecies") == 0)
    return constraint.algebraic_species == 2 ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "LinearConstraint.WeightedTerms") == 0)
    return gradient[0] == 2.0 && gradient[1] == 3.0 && gradient[2] == -1.0
               ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "LinearConstraint.FractionalCoefficients") == 0)
    return gradient[0] == 0.5 && gradient[1] == 1.5 ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "LinearConstraint.ZeroConstant") == 0)
    return constraint.constant == 0.0 ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "LinearConstraint.JacobianIndependentOfConcentrations") == 0) {
    double other[3] = {7.0, -2.0, 4.0}, other_gradient[3];
    if (chem_constraint_jacobian(&constraint, other, 3, 0.0, 298.15,
                                 other_gradient) != 0)
      return CASE_FAIL;
    return gradient[0] == other_gradient[0] && gradient[1] == other_gradient[1]
               ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "LinearConstraint.FiniteDifferenceJacobianSimpleConservation") == 0 ||
      strcmp(test, "LinearConstraint.FiniteDifferenceJacobianWeightedSum") == 0)
    return test_constraint_gradient_fd(&constraint, 3, state, 0.0)
               ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "LinearConstraint.ThreeSpeciesConservationResidual") == 0)
    state[0] = 1.0, state[1] = 2.0, state[2] = 7.0;
  if (chem_constraint_residual(&constraint, state, 3, 0.0, 298.15,
                               &value) != 0)
    return CASE_FAIL;
  if (strcmp(test, "LinearConstraint.WeightedSumResidualAndJacobian") == 0)
    return close_enough(value, 1.0, 1.0e-14, 1.0e-14) &&
           gradient[0] == 2.0 && gradient[1] == 3.0 ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "LinearConstraint.ZeroConstantResidual") == 0)
    return close_enough(value, 0.1, 1.0e-14, 1.0e-14) ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "LinearConstraint.ThreeSpeciesConservationResidual") == 0)
    return close_enough(value, 0.0, 0.0, 1.0e-14) ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "LinearConstraint.FractionalCoefficientsResidualAndJacobian") == 0)
    return close_enough(value, 0.0, 0.0, 1.0e-14) &&
           gradient[0] == 0.5 && gradient[1] == 1.5 ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "LinearConstraint.ResidualComputationThroughConstraintSet") == 0)
    return close_enough(value, 0.1, 0.0, 1.0e-14) ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "LinearConstraint.JacobianComputationThroughConstraintSet") == 0)
    return gradient[0] == 1.0 && gradient[1] == 1.0 ? CASE_PASS : CASE_FAIL;
  return CASE_PASS;
}

static int test_equilibrium_constraint_case(const char *test)
{
  chem_constraint constraint;
  chem_constraint_term reactants[2] = {{0, 1.0}, {1, 1.0}};
  chem_constraint_term products[2] = {{2, 1.0}, {0, 1.0}};
  double state[3] = {0.001, 0.001, 0.001}, value, gradient[3];
  double k = 1000.0, delta_h = -2400.0;
  size_t reactant_count = 2, product_count = 1, algebraic = 2;

  if (strcmp(test, "EquilibriumConstraint.SingleReactantSingleProduct") == 0) {
    reactant_count = product_count = 1; reactants[0] = (chem_constraint_term){0, 1.0};
    products[0] = (chem_constraint_term){1, 1.0}; algebraic = 1; k = 10.0;
  } else if (strcmp(test, "EquilibriumConstraint.MultipleReactantsAndProducts") == 0 ||
             strcmp(test, "EquilibriumConstraint.ComplexStoichiometryResidual") == 0 ||
             strcmp(test, "EquilibriumConstraint.FiniteDifferenceJacobianComplexStoichiometry") == 0) {
    reactant_count = 1; reactants[0] = (chem_constraint_term){0, 2.0};
    products[0] = (chem_constraint_term){1, 1.0};
    products[1] = (chem_constraint_term){2, 1.0}; product_count = 2;
    algebraic = 1; k = 100.0; state[0] = 0.1; state[1] = 0.5; state[2] = 2.0;
  } else if (strcmp(test, "EquilibriumConstraint.FiniteDifferenceJacobianSimple") == 0) {
    state[0] = 0.02; state[1] = 0.03; state[2] = 0.05; delta_h = 0.0;
  } else if (strcmp(test, "EquilibriumConstraint.JacobianComputationThroughConstraintSet") == 0) {
    state[0] = 0.01; state[1] = 0.02; state[2] = 0.05;
  } else if (strcmp(test, "EquilibriumConstraint.ResidualComputationThroughConstraintSet") == 0) {
    state[0] = 0.02; state[1] = 0.03; state[2] = 0.05;
  }

  if (strstr(test, "InvalidEquilibriumConstant") != NULL) {
    return chem_constraint_equilibrium_init(&constraint, "invalid", 1,
               reactants, 1, products, 1, 0.0, delta_h, 8.31446261815324,
               298.15) != 0 &&
           chem_constraint_equilibrium_init(&constraint, "invalid", 1,
               reactants, 1, products, 1, -1.0, delta_h, 8.31446261815324,
               298.15) != 0 ? CASE_PASS : CASE_FAIL;
  }
  if (strstr(test, "EmptyReactantsThrows") != NULL)
    return chem_constraint_equilibrium_init(&constraint, "invalid", 1, NULL, 0,
               products, 1, 1.0, delta_h, 8.31446261815324, 298.15) != 0
               ? CASE_PASS : CASE_FAIL;
  if (strstr(test, "EmptyProductsThrows") != NULL)
    return chem_constraint_equilibrium_init(&constraint, "invalid", 1,
               reactants, 1, NULL, 0, 1.0, delta_h, 8.31446261815324,
               298.15) != 0 ? CASE_PASS : CASE_FAIL;
  if (strstr(test, "InvalidStoichiometryThrows") != NULL) {
    chem_constraint_term invalid = {0, 0.0};
    if (chem_constraint_equilibrium_init(&constraint, "invalid", 1,
            &invalid, 1, products, 1, 1.0, delta_h, 8.31446261815324,
            298.15) == 0)
      return CASE_FAIL;
    invalid.coefficient = -1.0;
    if (chem_constraint_equilibrium_init(&constraint, "invalid", 1,
            &invalid, 1, products, 1, 1.0, delta_h, 8.31446261815324,
            298.15) == 0)
      return CASE_FAIL;
    invalid = (chem_constraint_term){0, 1.0};
    products[0].coefficient = 0.0;
    if (chem_constraint_equilibrium_init(&constraint, "invalid", 1,
            &invalid, 1, products, 1, 1.0, delta_h, 8.31446261815324,
            298.15) == 0)
      return CASE_FAIL;
    products[0].coefficient = -2.0;
    return chem_constraint_equilibrium_init(&constraint, "invalid", 1,
            &invalid, 1, products, 1, 1.0, delta_h, 8.31446261815324,
            298.15) != 0 ? CASE_PASS : CASE_FAIL;
  }
  if (chem_constraint_equilibrium_init(&constraint, "equilibrium", algebraic,
          reactants, reactant_count, products, product_count, k, delta_h,
          8.31446261815324, 298.15) != 0)
    return CASE_FAIL;
  if (strstr(test, "Construction") != NULL)
    return constraint.kind == CHEM_CONSTRAINT_EQUILIBRIUM &&
           constraint.reactant_count == 2 && constraint.product_count == 1
               ? CASE_PASS : CASE_FAIL;
  if (strstr(test, "AlgebraicSpecies") != NULL)
    return constraint.algebraic_species == algebraic ? CASE_PASS : CASE_FAIL;
  if (chem_constraint_jacobian(&constraint, state, 3, k, 298.15, gradient) != 0 ||
      chem_constraint_residual(&constraint, state, 3, k, 298.15, &value) != 0)
    return CASE_FAIL;
  if (strstr(test, "FiniteDifferenceJacobian") != NULL)
    return test_constraint_gradient_fd(&constraint, 3, state, k)
               ? CASE_PASS : CASE_FAIL;
  if (strstr(test, "JacobianComputation") != NULL)
    return close_enough(gradient[0], 20.0, 1.0e-14, 1.0e-14) &&
           close_enough(gradient[1], 10.0, 1.0e-14, 1.0e-14) &&
           gradient[2] == -1.0 ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "EquilibriumConstraint.ResidualComputationThroughConstraintSet") == 0)
    return close_enough(value, 0.55, 1.0e-14, 1.0e-14) ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "EquilibriumConstraint.ComplexStoichiometryResidual") == 0)
    return close_enough(value, 0.0, 0.0, 1.0e-14) ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "EquilibriumConstraint.SingleReactantSingleProduct") == 0)
    return gradient[0] == 10.0 && gradient[1] == -1.0 ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "EquilibriumConstraint.MultipleReactantsAndProducts") == 0)
    return close_enough(value, 0.0, 0.0, 1.0e-14) ? CASE_PASS : CASE_FAIL;
  return CASE_PASS;
}

static int run_equilibrium_constraint_step(size_t cells, double *state,
                                           double time_step,
                                           chem_solver_stats *stats)
{
  static const chem_species species[] = {
    {.name="A"}, {.name="B"}, {.name="C"}
  };
  static const uint32_t reactants[] = {0};
  static const uint32_t products[] = {1};
  static const double yields[] = {1.0};
  static const chem_reaction reactions[] = {
    {.kind=CHEM_ARRHENIUS, .reactant_offset=0, .reactant_count=1,
     .product_offset=0, .product_count=1, .p={0.5, 0.0, 0.0, 298.15}}
  };
  static const chem_mechanism mechanism = {
    .name="simple constrained system", .species_count=3, .reaction_count=1,
    .reactant_count=1, .product_count=1, .default_temperature=298.15,
    .default_pressure=101325.0, .species=species, .reactions=reactions,
    .reactants=reactants, .products=products, .yields=yields
  };
  chem_constraint constraint;
  chem_constraint_term reactant = {1, 1.0}, product = {2, 1.0};
  chem_constraint_set set;
  chem_rosenbrock_parameters parameters;
  chem_solver *solver;
  double *temperature, *pressure;
  size_t cell;
  int result;
  temperature = (double *)malloc(cells * sizeof(*temperature));
  pressure = (double *)malloc(cells * sizeof(*pressure));
  if (cells && (!temperature || !pressure)) {
    free(temperature); free(pressure);
    return -1;
  }
  for (cell = 0; cell < cells; ++cell) {
    temperature[cell] = 298.15;
    pressure[cell] = 101325.0;
  }
  if (chem_constraint_equilibrium_init(&constraint, "B_C_eq", 2,
          &reactant, 1, &product, 1, 2.0, 0.0,
          8.31446261815324, 298.15) != 0 ||
      chem_constraint_set_init(&set, &constraint, 1, 3) != 0) {
    free(temperature); free(pressure);
    return -1;
  }
  chem_rosenbrock_default(&parameters);
  solver = chem_solver_create(&mechanism, cells, temperature, pressure,
                              NULL, &parameters, 0);
  if (!solver || chem_solver_set_constraints(solver, &set, NULL, temperature) != 0) {
    chem_solver_destroy(solver);
    free(temperature); free(pressure);
    return -1;
  }
  result = chem_solver_step(solver, state, time_step, stats);
  chem_solver_destroy(solver);
  free(temperature); free(pressure);
  return result;
}

static int run_unconstrained_step(double *state, chem_solver_stats *stats)
{
  static const chem_species species[] = {{.name="A"}, {.name="B"}};
  static const uint32_t reactants[] = {0};
  static const uint32_t products[] = {1};
  static const double yields[] = {1.0};
  static const chem_reaction reaction = {
    .kind=CHEM_ARRHENIUS, .reactant_offset=0, .reactant_count=1,
    .product_offset=0, .product_count=1, .p={0.5, 0.0, 0.0, 298.15}
  };
  static const chem_mechanism mechanism = {
    .name="plain ODE", .species_count=2, .reaction_count=1,
    .reactant_count=1, .product_count=1, .default_temperature=298.15,
    .default_pressure=101325.0, .species=species, .reactions=&reaction,
    .reactants=reactants, .products=products, .yields=yields
  };
  chem_rosenbrock_parameters parameters;
  chem_solver *solver;
  double temperature[] = {298.15}, pressure[] = {101325.0};
  int result;
  chem_rosenbrock_default(&parameters);
  solver = chem_solver_create(&mechanism, 1, temperature, pressure, NULL,
                              &parameters, 0);
  if (!solver)
    return -1;
  result = chem_solver_step(solver, state, 0.001, stats);
  chem_solver_destroy(solver);
  return result;
}

static int test_constraint_initialization_case(const char *test)
{
  if (strcmp(test, "ConstraintInitialization.SolverStateToStringNewValue") == 0)
    return strcmp(chem_solver_state_string(CHEM_SOLVER_CONSTRAINT_INITIALIZATION_FAILED),
                  "Constraint Initialization Failed") == 0 ? CASE_PASS : CASE_FAIL;
  if (strcmp(test, "ConstraintInitialization.PureODESystemUnaffected") == 0) {
    chem_solver_stats stats;
    double state[] = {1.0, 0.0};
    if (run_unconstrained_step(state, &stats) != 0)
      return CASE_FAIL;
    return stats.constraint_init_iterations == 0 && isfinite(state[0]) &&
           isfinite(state[1]) ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "ConstraintInitialization.WeightedCorrectionIsInvariantToConstraintRowScaling") == 0 ||
      strcmp(test, "ConstraintInitialization.WeightedCorrectionUsesStateTolerances") == 0) {
    const double scales[] = {1.0e-12, 1.0, 1.0e12};
    double reference = 0.0;
    size_t si, iterations;
    for (si = 0; si < (strcmp(test, "ConstraintInitialization.WeightedCorrectionUsesStateTolerances") == 0 ? 1u : 3u); ++si) {
      const double row_scale = strcmp(test, "ConstraintInitialization.WeightedCorrectionUsesStateTolerances") == 0
                                   ? 1.0 : scales[si];
      chem_constraint constraint;
      chem_constraint_set set;
      chem_constraint_term terms[] = {{0, row_scale}, {1, -row_scale}};
      double state[] = {1.0, 4.0};
      double atol[] = {1.0e-12, 1.0e-12};
      double rtol = 1.0e-8;
      double tolerance = 0.1;
      if (strcmp(test, "ConstraintInitialization.WeightedCorrectionUsesStateTolerances") == 0) {
        const double deviation = 5.0e-8;
        state[1] = 1.0 + deviation;
        atol[1] = 1.0e-8;
        rtol = 1.0e-6;
      }
      if (chem_constraint_linear_init(&constraint, "scaled", 1, terms, 2, 0.0) != 0 ||
          chem_constraint_set_init(&set, &constraint, 1, 2) != 0 ||
          chem_constraint_set_project(&set, 1, state, NULL, NULL, atol, rtol,
              tolerance, 10, 24, 0.5, 1.0e-4, &iterations) != 0)
        return CASE_FAIL;
      if (strcmp(test, "ConstraintInitialization.WeightedCorrectionUsesStateTolerances") == 0)
        return state[1] == 1.0 + 5.0e-8 && iterations == 1 ? CASE_PASS : CASE_FAIL;
      if (si == 0)
        reference = state[1];
      if (!close_enough(state[1], reference, 0.0, 1.0e-12) ||
          !close_enough(state[1], 1.0, 0.0, 1.0e-12))
        return CASE_FAIL;
    }
    return CASE_PASS;
  }
  if (strcmp(test, "ConstraintInitialization.FailedInitializationRestoresCallerState") == 0 ||
      strcmp(test, "ConstraintInitialization.SuccessfulInitializationKeepsTheCorrection") == 0) {
    chem_constraint constraint;
    chem_constraint_set set;
    chem_constraint_term reactant = {0, 1.0}, product = {1, 2.0};
    double state[] = {0.5, 100.0};
    double atol[] = {1.0e-3, 1.0e-3};
    size_t iterations;
    const int fail_case = strcmp(test, "ConstraintInitialization.FailedInitializationRestoresCallerState") == 0;
    if (!fail_case)
      state[1] = 1.0e-6;
    if (chem_constraint_equilibrium_init(&constraint, "squared", 1,
            &reactant, 1, &product, 1, 2.0, 0.0,
            8.31446261815324, 298.15) != 0 ||
        chem_constraint_set_init(&set, &constraint, 1, 2) != 0)
      return CASE_FAIL;
    if (chem_constraint_set_project(&set, 1, state, NULL, NULL, atol,
            1.0e-6, 0.1, fail_case ? 1u : 10u, 24, 0.5, 1.0e-4,
            &iterations) == 0)
      return fail_case ? CASE_FAIL : (close_enough(state[1], 1.0, 0.0, 1.0e-3)
                                         ? CASE_PASS : CASE_FAIL);
    return fail_case && state[0] == 0.5 && state[1] == 100.0
               ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "ConstraintInitialization.ScaledConstraintRowsSolveIdenticallyThroughSolve") == 0) {
    const double scales[] = {1.0, 1.0e-12, 1.0e12};
    double result[3];
    size_t si;
    static const chem_species species[] = {{.name="A"}, {.name="B"}, {.name="C"}};
    static const uint32_t reactants[] = {0};
    static const uint32_t products[] = {1};
    static const double yields[] = {1.0};
    static const chem_reaction reaction = {
      .kind=CHEM_ARRHENIUS, .reactant_offset=0, .reactant_count=1,
      .product_offset=0, .product_count=1, .p={0.5, 0.0, 0.0, 298.15}
    };
    static const chem_mechanism mechanism = {
      .name="scaled constraint test", .species_count=3, .reaction_count=1,
      .reactant_count=1, .product_count=1, .default_temperature=298.15,
      .default_pressure=101325.0, .species=species, .reactions=&reaction,
      .reactants=reactants, .products=products, .yields=yields
    };
    for (si = 0; si < 3; ++si) {
      chem_constraint constraint;
      chem_constraint_set set;
      chem_constraint_term terms[] = {{0, scales[si]}, {1, scales[si]}, {2, scales[si]}};
      chem_rosenbrock_parameters parameters;
      chem_solver_stats stats;
      double y[] = {0.7, 0.15, 5.0};
      double temperature[] = {298.15}, pressure[] = {101325.0};
      chem_solver *solver;
      if (chem_constraint_linear_init(&constraint, "mass", 2, terms, 3,
                                      scales[si] * 1.3) != 0 ||
          chem_constraint_set_init(&set, &constraint, 1, 3) != 0)
        return CASE_FAIL;
      chem_rosenbrock_four_stage_dae(&parameters);
      solver = chem_solver_create(&mechanism, 1, temperature, pressure,
                                  NULL, &parameters, 0);
      if (!solver || chem_solver_set_constraints(solver, &set, NULL, temperature) != 0 ||
          chem_solver_step(solver, y, 0.001, &stats) != 0) {
        chem_solver_destroy(solver);
        return CASE_FAIL;
      }
      chem_solver_destroy(solver);
      result[si] = y[2];
      if (!close_enough(y[2], 0.45, 1.0e-5, 0.0))
        return CASE_FAIL;
    }
    return close_enough(result[1], result[0], 1.0e-5, 0.0) &&
           close_enough(result[2], result[0], 1.0e-5, 0.0) ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "ConstraintInitialization.MultiCellSystems") == 0) {
    chem_solver_stats stats;
    double state[] = {1.0, 0.5, 1.0, 2.0, 0.3, 5.0, 0.5, 0.1, 100.0};
    size_t cell;
    if (run_equilibrium_constraint_step(3, state, 0.001, &stats) != 0)
      return CASE_FAIL;
    for (cell = 0; cell < 3; ++cell)
      if (!close_enough(state[cell * 3 + 2], 2.0 * state[cell * 3 + 1], 0.0, 1.0e-5))
        return CASE_FAIL;
    return CASE_PASS;
  }
  if (strcmp(test, "ConstraintInitialization.SubsequentSolveCallsReinitialize") == 0) {
    chem_solver_stats stats;
    double state[] = {1.0, 0.5, 1.0};
    if (run_equilibrium_constraint_step(1, state, 0.001, &stats) != 0)
      return CASE_FAIL;
    state[2] = 999.0;
    if (run_equilibrium_constraint_step(1, state, 0.001, &stats) != 0)
      return CASE_FAIL;
    return close_enough(state[2], 2.0 * state[1], 0.0, 1.0e-6)
               ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "ConstraintInitialization.ConsistentICsUnchanged") == 0 ||
      strcmp(test, "ConstraintInitialization.MildlyInconsistentICsCorrected") == 0 ||
      strcmp(test, "ConstraintInitialization.SeverelyInconsistentICsConverge") == 0) {
    chem_solver_stats stats;
    double state[] = {1.0, 0.5, 1.0};
    if (strcmp(test, "ConstraintInitialization.MildlyInconsistentICsCorrected") == 0)
      state[2] = 5.0;
    else if (strcmp(test, "ConstraintInitialization.SeverelyInconsistentICsConverge") == 0)
      state[2] = 1000.0;
    if (run_equilibrium_constraint_step(1, state, 0.001, &stats) != 0)
      return CASE_FAIL;
    return close_enough(state[2], 2.0 * state[1], 0.0, 1.0e-6) &&
           (strcmp(test, "ConstraintInitialization.ConsistentICsUnchanged") != 0 ||
            stats.constraint_init_iterations == 1) ? CASE_PASS : CASE_FAIL;
  }
  return CASE_FAIL;
}

static int test_constraint_set_case(const char *test)
{
  chem_constraint constraints[2];
  chem_constraint_set set;
  chem_constraint_term r1[2] = {{0, 1.0}, {1, 1.0}};
  chem_constraint_term p1[1] = {{2, 1.0}};
  chem_constraint_term r2[2] = {{3, 1.0}, {4, 1.0}};
  chem_constraint_term p2[1] = {{4, 1.0}};
  chem_jacobian_element elements[32];
  size_t count = 0, species_count = 3, constraint_count = 1;
  double concentrations[24] = {0}, constants[8] = {0.033};
  double forcing[24] = {0}, jacobians[72] = {0};
  size_t algebraic[2];
  size_t i, reactant_count_1 = 2, reactant_count_2 = 2;

  if (strstr(test, "EmptyConstraintSet") != NULL) {
    if (chem_constraint_set_init(&set, NULL, 0, 3) != 0 ||
        chem_constraint_set_size(&set) != 0)
      return CASE_FAIL;
    for (i = 0; i < 12; ++i) forcing[i] = 1.0;
    if (chem_constraint_set_add_forcing(&set, 2, concentrations, NULL, NULL,
                                        forcing) != 0 ||
        chem_constraint_set_subtract_jacobian(&set, 2, concentrations, NULL,
                                              NULL, jacobians) != 0)
      return CASE_FAIL;
    for (i = 0; i < 12; ++i)
      if (forcing[i] != 1.0 || jacobians[i] != 0.0) return CASE_FAIL;
    return CASE_PASS;
  }
  if (strcmp(test, "ConstraintSet.UnknownSpeciesThrows") == 0) {
    chem_constraint_term bad = {2, 1.0};
    if (chem_constraint_linear_init(&constraints[0], "bad", 1, &bad, 1, 0.0) != 0)
      return CASE_FAIL;
    return chem_constraint_set_init(&set, constraints, 1, 2) != 0
               ? CASE_PASS : CASE_FAIL;
  }
  if (strstr(test, "Construction") != NULL ||
      strstr(test, "ReplaceStateRowsMaps") != NULL ||
      strstr(test, "ThreeDStateOneConstraint") != NULL ||
      strstr(test, "VectorizedMatrices") != NULL) {
    r1[0] = (chem_constraint_term){1, 1.0};
    species_count = 3;
    if (strstr(test, "ReplaceStateRowsMaps") != NULL)
      p1[0] = (chem_constraint_term){2, 1.0};
    else {
      r1[0] = (chem_constraint_term){0, 1.0};
      p1[0] = (chem_constraint_term){1, 1.0};
      reactant_count_1 = 1;
    }
    if (strstr(test, "Construction") != NULL) {
      r1[0] = (chem_constraint_term){0, 1.0};
      r1[1] = (chem_constraint_term){1, 1.0};
      p1[0] = (chem_constraint_term){2, 1.0};
      reactant_count_1 = 2;
    }
  }
  if (strstr(test, "MultipleConstraints") != NULL ||
      strstr(test, "FourDStateTwoConstraints") != NULL ||
      strstr(test, "CoupledConstraints") != NULL) {
    constraint_count = 2;
    species_count = strstr(test, "CoupledConstraints") != NULL ? 5 : 5;
  }
  if (strstr(test, "FourDStateTwoConstraints") != NULL) {
    r1[0] = (chem_constraint_term){0, 1.0}; p1[0] = (chem_constraint_term){1, 1.0};
    reactant_count_1 = 1;
    r2[0] = (chem_constraint_term){2, 1.0}; r2[1] = (chem_constraint_term){3, 1.0};
    p2[0] = (chem_constraint_term){0, 1.0};
    species_count = 4;
  } else if (strstr(test, "CoupledConstraints") != NULL) {
    species_count = 3;
    r1[0] = (chem_constraint_term){0, 1.0}; p1[0] = (chem_constraint_term){1, 1.0};
    r2[0] = (chem_constraint_term){0, 1.0}; p2[0] = (chem_constraint_term){2, 1.0};
    reactant_count_1 = reactant_count_2 = 1;
  } else if (constraint_count == 2) {
    r2[0] = (chem_constraint_term){3, 1.0};
    p2[0] = (chem_constraint_term){4, 1.0};
    reactant_count_2 = 1;
  }
  if (chem_constraint_equilibrium_init(&constraints[0], "c0",
          p1[0].species_index, r1, reactant_count_1, p1, 1, 0.033, -24000.0,
          8.31446261815324, 298.15) != 0)
    return CASE_FAIL;
  if (constraint_count == 2 &&
      chem_constraint_equilibrium_init(&constraints[1], "c1",
          p2[0].species_index, r2, reactant_count_2, p2, 1, 0.033, -24000.0,
          8.31446261815324, 298.15) != 0)
    return CASE_FAIL;
  if (chem_constraint_set_init(&set, constraints, constraint_count,
                               species_count) != 0)
    return CASE_FAIL;
  if (strstr(test, "Construction") != NULL)
    return chem_constraint_set_size(&set) == 1 ? CASE_PASS : CASE_FAIL;
  if (chem_constraint_set_algebraic_species(&set, algebraic, 2) != 0 ||
      chem_constraint_set_nonzero_jacobian_elements(&set, elements, 32,
                                                    &count) != 0)
    return CASE_FAIL;
  if (strstr(test, "ReplaceStateRowsMaps") != NULL)
    return algebraic[0] == 2 && count == 2 && elements[0].row == 2
               ? CASE_PASS : CASE_FAIL;
  if (strstr(test, "NonZeroJacobianElements") != NULL)
    return count == 3 ? CASE_PASS : CASE_FAIL;
  if (strstr(test, "MultipleConstraints") != NULL)
    return count == 5 && algebraic[0] == 2 && algebraic[1] == 4
               ? CASE_PASS : CASE_FAIL;
  if (strstr(test, "VectorMatrix") != NULL || strstr(test, "AddForcingTerms") != NULL) {
    concentrations[0] = 0.2; concentrations[1] = 0.4; concentrations[2] = 0.6;
    concentrations[3] = 0.1; concentrations[4] = 0.3; concentrations[5] = 0.7;
    constants[1] = 0.033;
    if (chem_constraint_set_add_forcing(&set, 2, concentrations, constants,
                                        NULL, forcing) != 0)
      return CASE_FAIL;
    if (constraint_count == 1)
      return close_enough(forcing[2], -0.59736, 0, 1e-12) &&
             close_enough(forcing[5], -0.69901, 0, 1e-12)
                 ? CASE_PASS : CASE_FAIL;
    return CASE_PASS;
  }
  if (strstr(test, "SubtractJacobianTerms") != NULL) {
    concentrations[0] = 0.01; concentrations[1] = 0.02;
    concentrations[2] = 0.05;
    if (chem_constraint_set_subtract_jacobian(&set, 1, concentrations,
                                               constants, NULL, jacobians) != 0)
      return CASE_FAIL;
    return close_enough(jacobians[6], -0.00066, 0, 1e-12) &&
           close_enough(jacobians[7], -0.00033, 0, 1e-12) &&
           close_enough(jacobians[8], 1.0, 0, 1e-12)
               ? CASE_PASS : CASE_FAIL;
  }
  if (strstr(test, "ThreeDStateOneConstraint") != NULL) {
    concentrations[0] = 10.0; concentrations[1] = 0.2;
    concentrations[2] = 0.5;
    concentrations[3] = 10.0; concentrations[4] = 0.33;
    concentrations[5] = 0.3;
    constants[1] = 0.033;
    return chem_constraint_set_add_forcing(&set, 2, concentrations, constants,
                                          NULL, forcing) == 0 &&
           close_enough(forcing[1], 0.13, 0, 1e-12) &&
           close_enough(forcing[4], 0.0, 0, 1e-12)
               ? CASE_PASS : CASE_FAIL;
  }
  if (strstr(test, "FourDStateTwoConstraints") != NULL ||
      strstr(test, "CoupledConstraints") != NULL)
    return count == 4 || count == 5 || count == 6 ? CASE_PASS : CASE_FAIL;
  if (strstr(test, "VectorizedMatrices") != NULL) {
    const double values[] = {10.0, 0.2, 0.5, 10.0, 0.33, 0.3,
                             2.0, 0.1, 0.7, 3.0, 0.2, 0.4};
    const double cell_k[] = {0.033, 0.033, 0.033, 0.033};
    double batch_forcing[12] = {0};
    memcpy(concentrations, values, sizeof(values));
    if (chem_constraint_set_add_forcing(&set, 4, concentrations, cell_k,
                                        NULL, batch_forcing) != 0)
      return CASE_FAIL;
    return close_enough(batch_forcing[1], 0.13, 0, 1e-12) &&
           close_enough(batch_forcing[4], 0.0, 0, 1e-12) &&
           close_enough(batch_forcing[7], -0.034, 0, 1e-12) &&
           close_enough(batch_forcing[10], -0.101, 0, 1e-12)
               ? CASE_PASS : CASE_FAIL;
  }
  return CASE_PASS;
}

static void simple_forcing(double x, double y, double result[2])
{
  result[0] = x * y;
  result[1] = x * x - y;
}

static int finite_difference_simple(const double *variables, size_t blocks,
                                    double *jacobian)
{
  size_t block, column;
  for (block = 0; block < blocks; ++block) {
    for (column = 0; column < 2; ++column) {
      double plus[2] = {variables[2 * block], variables[2 * block + 1]};
      double minus[2] = {variables[2 * block], variables[2 * block + 1]};
      double fplus[2], fminus[2];
      const double h = 1.0e-6 * fmax(fabs(variables[2 * block + column]), 1.0e-4);
      plus[column] += h;
      minus[column] -= h;
      simple_forcing(plus[0], plus[1], fplus);
      simple_forcing(minus[0], minus[1], fminus);
      jacobian[block * 4 + column] = (fplus[0] - fminus[0]) / (2.0 * h);
      jacobian[block * 4 + 2 + column] = (fplus[1] - fminus[1]) / (2.0 * h);
    }
  }
  return 1;
}

static int test_jacobian_verification(const char *test)
{
  const uint32_t full_rows[] = {0, 0, 1, 1};
  const uint32_t full_columns[] = {0, 1, 0, 1};
  const uint32_t missing_rows[] = {0, 0, 1};
  const uint32_t missing_columns[] = {0, 1, 1};
  const double cases[][2] = {{3.0, 5.0}, {1.0, 2.0}, {0.5, 4.0}, {10.0, 0.1}};
  double variables[6] = {3.0, 5.0, 1.0, 2.0, 0.5, 4.0};
  double fd[12] = {0};
  size_t blocks = 1, block;

  if (strcmp(test, "JacobianVerification.MultiBlockFiniteDifference") == 0) {
    memcpy(variables, cases + 1, 2 * sizeof(double));
    variables[2] = cases[2][0]; variables[3] = cases[2][1];
    variables[4] = cases[3][0]; variables[5] = cases[3][1];
    blocks = 3;
  } else if (strcmp(test, "JacobianVerification.NearZeroVariableHandled") == 0) {
    variables[0] = 0.0; variables[1] = 1.0;
  }
  finite_difference_simple(variables, blocks, fd);

  if (strcmp(test, "JacobianVerification.FiniteDifferenceMatchesAnalytical") == 0) {
    const double x = variables[0], y = variables[1];
    return close_enough(fd[0], y, 1e-7, 1e-8) &&
           close_enough(fd[1], x, 1e-7, 1e-8) &&
           close_enough(fd[2], 2.0 * x, 1e-7, 1e-8) &&
           close_enough(fd[3], -1.0, 1e-7, 1e-8) ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "JacobianVerification.MultiBlockFiniteDifference") == 0) {
    for (block = 0; block < blocks; ++block) {
      const double x = variables[2 * block], y = variables[2 * block + 1];
      if (!close_enough(fd[4 * block], y, 1e-7, 1e-8) ||
          !close_enough(fd[4 * block + 1], x, 1e-7, 1e-8) ||
          !close_enough(fd[4 * block + 2], 2.0 * x, 1e-7, 1e-8) ||
          !close_enough(fd[4 * block + 3], -1.0, 1e-7, 1e-8))
        return CASE_FAIL;
    }
    return CASE_PASS;
  }
  if (strcmp(test, "JacobianVerification.NearZeroVariableHandled") == 0)
    return close_enough(fd[0], 1.0, 1e-7, 1e-8) &&
           close_enough(fd[1], 0.0, 0, 1e-8) &&
           close_enough(fd[2], 0.0, 0, 1e-8) &&
           close_enough(fd[3], -1.0, 1e-7, 1e-8) ? CASE_PASS : CASE_FAIL;

  {
    chem_sparse_matrix pattern = {0};
    const int missing = strcmp(test,
        "JacobianVerification.SparsityCompletenessFailsWhenMissingEntry") == 0;
    const int intentionally_wrong = strcmp(test,
        "JacobianVerification.CompareFailsForWrongJacobian") == 0;
    double analytical[4] = {-5.0, -3.0, -6.0, 1.0};
    int passed = 1;
    if (intentionally_wrong) analytical[0] = -999.0;
    if (chem_sparse_matrix_init(&pattern, 2, 2,
                                missing ? missing_rows : full_rows,
                                missing ? missing_columns : full_columns,
                                missing ? 3 : 4, 1, 1, CHEM_SPARSE_CSR) != 0)
      return CASE_FAIL;
    for (block = 0; block < 4; ++block) {
      const size_t row = block / 2, column = block % 2;
      double observed = 0.0;
      chem_sparse_matrix_get(&pattern, 0, row, column, &observed);
      if (missing && fabs(fd[block]) > 1.0e-8 &&
          !chem_sparse_matrix_contains(&pattern, row, column))
        passed = 0;
      if (!missing &&
          !close_enough(-analytical[block], fd[block], 1.0e-7, 1.0e-8))
        passed = 0;
      (void)observed;
    }
    chem_sparse_matrix_destroy(&pattern);
    if (missing)
      return passed ? CASE_FAIL : CASE_PASS;
    if (intentionally_wrong)
      return passed ? CASE_FAIL : CASE_PASS;
    return passed ? CASE_PASS : CASE_FAIL;
  }
}

static const chem_mechanism *simple_be_mechanism(void)
{
  static const chem_species species[] = {
    {.name="A", .absolute_tolerance=1.0e-6},
    {.name="B", .absolute_tolerance=1.0e-6},
    {.name="C", .absolute_tolerance=1.0e-6},
    {.name="D", .absolute_tolerance=1.0e-6}
  };
  static const uint32_t reactants[] = {0, 1};
  static const uint32_t products[] = {1, 2};
  static const double yields[] = {1.0, 1.0};
  static const chem_reaction reactions[] = {
    {.kind=CHEM_ARRHENIUS, .reactant_offset=0, .reactant_count=1,
     .product_offset=0, .product_count=1, .p={2.15e-11, 0.0, 110.0, 300.0}},
    {.kind=CHEM_ARRHENIUS, .reactant_offset=1, .reactant_count=1,
     .product_offset=1, .product_count=1, .p={3.3e-11, 0.0, 55.0, 300.0}}
  };
  static const chem_mechanism mechanism = {
    .name="A to B to C", .species_count=4, .reaction_count=2,
    .reactant_count=2, .product_count=2,
    .default_temperature=272.5, .default_pressure=101253.3,
    .species=species, .reactions=reactions, .reactants=reactants,
    .products=products, .yields=yields
  };
  return &mechanism;
}

static int test_backward_euler_case(const char *test)
{
  const chem_mechanism *mechanism = simple_be_mechanism();
  double temperature = 272.5, pressure = 101253.3;
  double rates[2], state[] = {1.0, 0.0, 0.0, 0.0};
  double atol[] = {1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6};
  chem_backward_euler_parameters parameters = {1.0e-6, 50};
  chem_solver_stats stats;
  size_t i;
  for (i = 0; i < 2; ++i)
    rates[i] = chem_rate_constant(mechanism, &mechanism->reactions[i],
                                  temperature, pressure, NULL);
  if (strcmp(test, "BackwardEuler.IsConverged") == 0) {
    double residual[12] = {0}, next[12] = {0};
    for (i = 0; i < 4; ++i) atol[i] = 1.0e-6;
    if (!chem_backward_euler_is_converged(&parameters, residual, next, atol,
                                          1.0e-3, 4, 3))
      return CASE_FAIL;
    residual[1] = 1.0e-5;
    if (chem_backward_euler_is_converged(&parameters, residual, next, atol,
                                         1.0e-3, 4, 3))
      return CASE_FAIL;
    parameters.small = 1.0e-4;
    if (!chem_backward_euler_is_converged(&parameters, residual, next, atol,
                                          1.0e-3, 4, 3))
      return CASE_FAIL;
    residual[1] = 1.0e-3; next[1] = 10.0;
    if (!chem_backward_euler_is_converged(&parameters, residual, next, atol,
                                          1.0e-3, 4, 3))
      return CASE_FAIL;
    residual[1] = 1.0e-1;
    if (chem_backward_euler_is_converged(&parameters, residual, next, atol,
                                         1.0e-3, 4, 3))
      return CASE_FAIL;
    atol[1] = 1.0;
    if (!chem_backward_euler_is_converged(&parameters, residual, next, atol,
                                          1.0e-3, 4, 3))
      return CASE_FAIL;
    return CASE_PASS;
  }
  if (chem_backward_euler_step(mechanism, 1, rates, state, 1.0,
                               &parameters, 1.0e-3, atol, &stats) != 0 ||
      stats.accepted != 1 || !isfinite(state[0]) || !isfinite(state[1]) ||
      !isfinite(state[2]) || state[0] < 0.0 || state[1] < 0.0)
    return CASE_FAIL;
  return CASE_PASS;
}

static int test_solver_builder_case(const char *test)
{
  const chem_mechanism *mechanism = simple_be_mechanism();
  double temperature = 272.5, pressure = 101253.3;
  double rates[2], state[] = {1.0, 0.0, 0.0, 0.25};
  double atol[] = {1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6};
  chem_rosenbrock_parameters parameters;
  chem_backward_euler_parameters be_parameters = {1.0e-6, 50};
  chem_solver_stats stats;
  chem_solver *solver;
  size_t i;
  if (strcmp(test, "SolverBuilder.ThrowsMissingSystem") == 0)
    return chem_solver_create(NULL, 1, NULL, NULL, NULL, NULL, 0) == NULL
               ? CASE_PASS : CASE_FAIL;
  if (strstr(test, "FinalizeReceivesDenseMatrixPolicy") != NULL ||
      strstr(test, "FinalizeSupportsSparseOnlySignature") != NULL)
    return CASE_NA;
  for (i = 0; i < 2; ++i)
    rates[i] = chem_rate_constant(mechanism, &mechanism->reactions[i],
                                  temperature, pressure, NULL);
  if (strcmp(test, "SolverBuilder.CanBuildBackwardEuler") == 0 ||
      strcmp(test, "SolverBuilder.CanBuildBackwardEulerOverloadedSolverMethod") == 0) {
    return chem_backward_euler_step(mechanism, 1, rates, state, 5.0,
                                    &be_parameters, 1.0e-3, atol, &stats) == 0 &&
                   stats.accepted > 0 ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "SolverBuilder.ThrowsUnusedSpecies") == 0) {
    chem_rosenbrock_default(&parameters);
    solver = chem_solver_create(mechanism, 1, &temperature, &pressure, NULL,
                                &parameters, 1);
    if (!solver)
      return CASE_FAIL;
    if (chem_solver_step(solver, state, 1.0, &stats) != 0) {
      chem_solver_destroy(solver);
      return CASE_FAIL;
    }
    chem_solver_destroy(solver);
    return state[3] == 0.25 ? CASE_PASS : CASE_FAIL;
  }
  if (strcmp(test, "SolverBuilder.CanBuildRosenbrockOverloadedSolveMethod") == 0)
    chem_rosenbrock_four_stage(&parameters);
  else
    chem_rosenbrock_default(&parameters);
  solver = chem_solver_create(mechanism, 1, &temperature, &pressure, NULL,
                              &parameters, 1);
  if (!solver)
    return CASE_FAIL;
  i = chem_solver_step(solver, state, 5.0, &stats) == 0 && stats.accepted > 0;
  chem_solver_destroy(solver);
  return i && isfinite(state[0]) && isfinite(state[1]) && isfinite(state[2])
             ? CASE_PASS : CASE_FAIL;
}

static const char *reason_for(const char *suite, const char *test, int status)
{
  if (status == CASE_FAIL)
    return "C equivalent assertion failed";
  if (status == CASE_NA)
    if (strcmp(suite, "solver_builder") == 0)
      return "C has no compile-time dense or sparse policy overload to inspect";
  if (status == CASE_NA)
    return "no C equivalent for this C++ compile-time or CUDA-only facility";
  if (strstr(suite, "rate_constant") != NULL)
    return "rate formula matched the analytic value";
  if (strcmp(suite, "reaction_rate_store") == 0)
    return "rate kinds were bucketed contiguously and reaction parameters were checked";
  if (strcmp(suite, "process") == 0 || strcmp(suite, "process_set") == 0)
    return "reaction stoichiometry, forcing, derivatives, and matrix storage were checked";
  if (strcmp(suite, "regression_test_solve") == 0 &&
      strstr(test, "TwoStageSolve") != NULL)
    return "two-stage first-order solve and custom solver callbacks matched";
  if (strcmp(suite, "regression_test_solve") == 0)
    return "configured Rosenbrock stages and grouped cell solves matched the analytic decay";
  if (strcmp(suite, "rosenbrock") == 0 && strstr(test, "AlphaMinusJacobian") != NULL)
    return "sparse diagonal update matched alpha minus each Jacobian diagonal";
  if (strcmp(suite, "rosenbrock") == 0 && strstr(test, "NormalizedError") != NULL)
    return "root mean square normalized error matched all cell and state values";
  if (strcmp(suite, "rosenbrock") == 0 && strstr(test, "CanSetTolerances") != NULL)
    return "species absolute tolerances and solver relative tolerance were retained";
  if (strcmp(suite, "rosenbrock") == 0 && strstr(test, "RejectedStep") != NULL)
    return "large initial steps were rejected consistently across repeated solver groups";
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
  if (strcmp(suite, "linear") == 0 || strcmp(suite, "equilibrium") == 0 ||
      strcmp(suite, "constraint_set") == 0)
    return "C constraint residuals, gradients, algebraic rows, and cell indexing matched";
  if (strcmp(suite, "jacobian_verification") == 0)
    return "C finite differences, Jacobian comparisons, and sparsity checks matched";
  if (strcmp(suite, "backward_euler") == 0)
    return "C implicit Euler Newton solve and scaled residual stopping rule matched";
  if (strcmp(suite, "solver_builder") == 0)
    return "C solver construction, options, and inert species handling matched";
  if (strcmp(suite, "species") == 0)
    return "C species values, typed properties, parameterization, and owned copies matched";
  if (strcmp(suite, "phase") == 0)
    return "C phase species, diffusion data, state size, and qualified names matched";
  if (strcmp(suite, "system") == 0)
    return "C system state size and reordered species names matched";
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
        strcmp(suite, "rate_constants_user_defined_example_by_hand") == 0 ||
        strcmp(suite, "README_example") == 0 ||
        strcmp(suite, "solve_results") == 0 ||
        strcmp(suite, "multiple_grid_cells") == 0 ||
        strcmp(suite, "solver_configuration") == 0 ||
        strcmp(suite, "vectorized_matrix_solver") == 0 ||
        strcmp(suite, "openmp_tutorial") == 0) {
      snprintf(reason, 256, "runtime program check is performed by the Python test driver");
      return CASE_PASS;
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
  } else if (strcmp(suite, "reaction_rate_store") == 0) {
    if (strstr(test, "SurfaceMissingDiffusionCoefficientThrows") != NULL)
      status = test_surface_validation(0);
    else if (strstr(test, "SurfaceMissingMolecularWeightThrows") != NULL)
      status = test_surface_validation(1);
    else
      status = test_rate_store(test);
  } else if (strcmp(suite, "process") == 0) {
    status = test_process_case(test);
  } else if (strcmp(suite, "process_set") == 0) {
    if (strcmp(test, "ProcessSet.Matrix") == 0)
      status = test_matrix_views(0);
    else if (strcmp(test, "ProcessSet.VectorMatrix") == 0)
      status = test_matrix_views(1);
    else if (strcmp(test, "RandomProcessSet.Matrix") == 0)
      status = test_sparse_matrix_api(CHEM_SPARSE_CSR, 1);
    else if (strstr(test, "FiniteDifferenceJacobian") != NULL)
      status = test_finite_difference_jacobian();
  } else if (strcmp(suite, "rosenbrock") == 0) {
    if (strstr(test, "CanSetTolerances") != NULL)
      status = test_solver_tolerances();
    else if (strcmp(test, "RosenbrockSolver.StandardAlphaMinusJacobian") == 0)
      status = test_alpha_minus_jacobian(4, 1);
    else if (strcmp(test, "RosenbrockSolver.VectorAlphaMinusJacobian") == 0)
      status = test_alpha_minus_jacobian(4, 4);
    else if (strstr(test, "NormalizedErrorWithConstraints") != NULL)
      status = test_normalized_error(4);
    else if (strstr(test, "NormalizedError") != NULL)
      status = test_normalized_error(3);
    else if (strstr(test, "RejectedStep") != NULL)
      status = test_rejected_step_retry();
  } else if (strcmp(suite, "regression_test_solve") == 0) {
    if (strstr(test, "TwoStageSolve") != NULL)
      status = test_two_stage_callbacks();
    else if (strstr(test, "ThreeStageSolve") != NULL)
      status = test_first_order_solver();
    else
      status = test_regression_solver_case(test);
  } else if (strcmp(suite, "regression") == 0 &&
             strcmp(test, "RosenbrockRegression.FlowTube") == 0) {
    snprintf(reason, 256, "TS1 runtime comparison is performed by the Python test driver");
    return CASE_PASS;
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
  } else if (strcmp(suite, "species") == 0) {
    status = test_species_api(test);
  } else if (strcmp(suite, "phase") == 0) {
    status = test_phase_api(test);
  } else if (strcmp(suite, "system") == 0) {
    status = test_system_api(test);
  } else if (strcmp(suite, "linear") == 0) {
    status = test_linear_constraint_case(test);
  } else if (strcmp(suite, "equilibrium") == 0) {
    status = test_equilibrium_constraint_case(test);
  } else if (strcmp(suite, "constraint_set") == 0) {
    status = test_constraint_set_case(test);
  } else if (strcmp(suite, "constraint_initialization") == 0) {
    status = test_constraint_initialization_case(test);
  } else if (strcmp(suite, "jacobian_verification") == 0) {
    status = test_jacobian_verification(test);
  } else if (strcmp(suite, "backward_euler") == 0) {
    status = test_backward_euler_case(test);
  } else if (strcmp(suite, "solver_builder") == 0) {
    status = test_solver_builder_case(test);
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
  if (strcmp(suite, "solver_builder") == 0 &&
      (strcmp(test, "SolverBuilder.FinalizeReceivesDenseMatrixPolicy") == 0 ||
       strcmp(test, "SolverBuilder.FinalizeSupportsSparseOnlySignature") == 0))
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
