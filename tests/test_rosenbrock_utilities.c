#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static const chem_species error_species[] = {
  { "A", 0.0, 0.0, 0.0, 1.0e-8, 0u },
  { "B", 0.0, 0.0, 0.0, 2.0e-8, 0u },
  { "C", 0.0, 0.0, 0.0, 3.0e-8, 0u },
};

static const chem_mechanism error_mechanism = {
  .name = "Rosenbrock normalized error",
  .species_count = 3u,
  .species = error_species,
  .relative_tolerance = 1.0e-6,
};

static int alpha_minus_jacobian(size_t block_count) {
  double matrix[4u * 25u];
  double expected[4u * 25u];
  for (size_t block = 0; block < block_count; ++block) {
    for (size_t row = 0; row < 5u; ++row) {
      for (size_t column = 0; column < 5u; ++column) {
        const size_t i = block * 25u + row * 5u + column;
        matrix[i] = (double)(i + 1u) / 7.0;
        expected[i] = matrix[i] + (row == column ? 42.042 : 0.0);
      }
    }
  }
  if (chem_alpha_minus_jacobian(matrix, block_count, 5u, 42.042) != 0) return 1;
  for (size_t i = 0; i < block_count * 25u; ++i)
    if (matrix[i] != expected[i]) return 1;
  return 0;
}

static int normalized_error(size_t cell_count) {
  double old_state[15], new_state[15], errors[15];
  double sum = 0.0;
  const size_t species_count = error_mechanism.species_count;
  for (size_t cell = 0; cell < cell_count; ++cell) {
    for (size_t species = 0; species < species_count; ++species) {
      const size_t i = cell * species_count + species;
      old_state[i] = 7.7 * (double)cell + (double)species;
      new_state[i] = -13.9 / (double)(species + 1u) - (double)cell;
      errors[i] = 81.57 / (double)(cell + 7u) / (double)(species + 3u);
      const double maximum = fmax(fabs(old_state[i]), fabs(new_state[i]));
      const double scale = error_species[species].absolute_tolerance +
          error_mechanism.relative_tolerance * maximum;
      const double scaled = errors[i] / scale;
      sum += scaled * scaled;
    }
  }
  const double expected = fmax(sqrt(sum / (double)(cell_count * species_count)), 1.0e-10);
  double actual = 0.0;
  if (chem_calculate_normalized_error(
          &error_mechanism, cell_count, old_state, new_state, errors, &actual) != 0) return 1;
  const double relative_error = fabs(actual - expected) / fmax(fabs(actual), fabs(expected));
  if (relative_error > 1.0e-11) {
    fprintf(stderr, "normalized error %.17g != %.17g\n", actual, expected);
    return 1;
  }

  memset(errors, 0, cell_count * species_count * sizeof(double));
  if (chem_calculate_normalized_error(
          &error_mechanism, cell_count, old_state, new_state, errors, &actual) != 0 || actual != 1.0e-10)
    return 1;
  return 0;
}

static int can_set_tolerances(void) {
  const chem_species species[] = {
    { "foo", 0.0, 0.0, 0.0, 1.0e-7, 0u },
    { "bar", 0.0, 0.0, 0.0, 1.0e-8, 0u },
  };
  return species[0].absolute_tolerance == 1.0e-7 && species[1].absolute_tolerance == 1.0e-8 ? 0 : 1;
}

static const chem_species stiff_species[] = {
  { "foo", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "bar", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "baz", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "quz", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "quuz", 0.0, 0.0, 0.0, 1.0e-3, 0u },
};
static const chem_term stiff_reactants[] = {
  { 0u, 1.0 }, { 2u, 1.0 }, { 1u, 1.0 }, { 3u, 1.0 },
};
static const chem_term stiff_products[] = {
  { 1u, 1.0 }, { 4u, 2.4 }, { 0u, 1.0 }, { 3u, 1.4 },
};
static const chem_reaction stiff_reactions[] = {
  { 0u, 2u, 0u, 2u, { .kind = CHEM_RATE_ARRHENIUS,
      .p = { [0] = 2.0e-11, [2] = 110.0, [3] = 300.0 } } },
  { 2u, 1u, 2u, 2u, { .kind = CHEM_RATE_ARRHENIUS,
      .p = { [0] = 1.0e-6, [3] = 300.0 } } },
  { 3u, 1u, 4u, 0u, { .kind = CHEM_RATE_ARRHENIUS,
      .p = { [0] = 3.5e-6, [3] = 300.0 } } },
};
static const chem_mechanism stiff_mechanism = {
  .name = "Rosenbrock rejected-step alpha",
  .species_count = 5u,
  .reaction_count = 3u,
  .reactant_term_count = 4u,
  .product_term_count = 4u,
  .species = stiff_species,
  .reactions = stiff_reactions,
  .reactant_terms = stiff_reactants,
  .product_terms = stiff_products,
  .temperature = 298.15,
  .pressure = 101325.0,
  .relative_tolerance = 1.0e-6,
};

static int rejected_step_alpha_matches(void) {
  chem_rosenbrock_parameters standard = chem_rosenbrock_default_parameters();
  chem_rosenbrock_parameters in_place = standard;
  standard.h_start_fraction = 1.0;
  in_place.h_start_fraction = 1.0;
  standard.lu_method = CHEM_LU_DOOLITTLE;
  in_place.lu_method = CHEM_LU_MOZART_IN_PLACE;
  double standard_state[] = { 1.0e12, 1.0e12, 1.0e12, 1.0e12, 1.0e12 };
  double in_place_state[5];
  memcpy(in_place_state, standard_state, sizeof(standard_state));
  const double temperature[] = { 298.15 }, pressure[] = { 101325.0 };
  chem_statistics standard_stats = { 0 }, in_place_stats = { 0 };
  if (chem_integrate_group_with_parameters(
          &stiff_mechanism, 1u, standard_state, NULL, temperature, pressure, 1000.0,
          &standard, &standard_stats) != 0 ||
      chem_integrate_group_with_parameters(
          &stiff_mechanism, 1u, in_place_state, NULL, temperature, pressure, 1000.0,
          &in_place, &in_place_stats) != 0) return 1;
  if (standard_stats.internal_steps < standard_stats.accepted + 2u ||
      standard_stats.internal_steps != in_place_stats.internal_steps ||
      standard_stats.accepted != in_place_stats.accepted ||
      standard_stats.rejected != in_place_stats.rejected) {
    fprintf(stderr, "rejected-step counts differ: standard %lu/%lu/%lu, in-place %lu/%lu/%lu\n",
        standard_stats.internal_steps, standard_stats.accepted, standard_stats.rejected,
        in_place_stats.internal_steps, in_place_stats.accepted, in_place_stats.rejected);
    return 1;
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  const char *test = strrchr(argv[1], '.');
  test = test ? test + 1 : argv[1];
  if (strcmp(test, "StandardAlphaMinusJacobian") == 0)
    return alpha_minus_jacobian(1u) || alpha_minus_jacobian(2u) || alpha_minus_jacobian(3u) || alpha_minus_jacobian(4u);
  if (strcmp(test, "VectorAlphaMinusJacobian") == 0)
    return alpha_minus_jacobian(1u) || alpha_minus_jacobian(4u) || alpha_minus_jacobian(3u) || alpha_minus_jacobian(2u);
  if (strcmp(test, "CanSetTolerances") == 0) return can_set_tolerances();
  if (strcmp(test, "StandardNormalizedError") == 0)
    return normalized_error(1u) || normalized_error(2u) || normalized_error(3u) || normalized_error(4u);
  if (strcmp(test, "VectorNormalizedError") == 0)
    return normalized_error(1u) || normalized_error(2u) || normalized_error(3u) || normalized_error(4u) ||
        normalized_error(5u);
  if (strcmp(test, "StandardNormalizedErrorWithConstraints") == 0)
    return normalized_error(2u) || normalized_error(5u);
  if (strcmp(test, "VectorNormalizedErrorWithConstraints") == 0)
    return normalized_error(3u) || normalized_error(5u);
  if (strcmp(test, "RejectedStepAlphaMatchesInPlaceSolver") == 0) return rejected_step_alpha_matches();
  return 2;
}
