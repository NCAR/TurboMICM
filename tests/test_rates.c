#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define GAS_CONSTANT (1.380649e-23 * 6.02214076e23)
#define AVOGADRO 6.02214076e23
#define PI 3.14159265358979323846264338327950288

static int close_to(double actual, double expected) {
  return isfinite(actual) && fabs(actual - expected) <= 1.0e-13 * fabs(expected) + 1.0e-30;
}

static int evaluate(chem_rate rate, const chem_species *species, size_t species_count,
                    const double *parameters, size_t parameter_count,
                    double temperature, double air_density, double *value) {
  chem_reaction reaction;
  memset(&reaction, 0, sizeof(reaction));
  reaction.rate = rate;
  chem_mechanism mechanism;
  memset(&mechanism, 0, sizeof(mechanism));
  mechanism.species_count = species_count;
  mechanism.reaction_count = 1;
  mechanism.parameter_count = parameter_count;
  mechanism.species = species;
  mechanism.reactions = &reaction;
  return chem_evaluate_rate(&mechanism, 0, parameters, temperature,
                            air_density * GAS_CONSTANT * temperature, value);
}

static chem_rate basic_rate(chem_rate_kind kind) {
  chem_rate rate;
  memset(&rate, 0, sizeof(rate));
  rate.kind = kind;
  rate.p[0] = 1.0;
  rate.p[3] = 300.0;
  rate.p[5] = 1.0;
  rate.p[8] = 1.0;
  rate.p[11] = 0.6;
  rate.p[12] = 1.0;
  rate.p[17] = 0.0;
  rate.p[18] = 298.0;
  rate.p[19] = 298.0;
  rate.p[20] = 1.0;
  rate.p[23] = 298.0;
  return rate;
}

static int arrhenius_cases(void) {
  const double t = 301.24;
  chem_rate p = basic_rate(CHEM_RATE_ARRHENIUS);
  const double cases[][2] = {
    { 1.0, 1.0 }, { 2.2e-10, 2.2e-10 },
    { 3.0e-11, 3.0e-11 * exp(-200.0 / t) },
    { 2.6e-12, 2.6e-12 * exp(-350.0 / t) }
  };
  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
    p.p[0] = cases[i][0];
    p.p[2] = i == 2 ? -200.0 : (i == 3 ? -350.0 : 0.0);
    double value;
    if (evaluate(p, NULL, 0, NULL, 0, t, 42.2, &value) || !close_to(value, cases[i][1])) return 0;
  }
  return 1;
}

static int branched_cases(void) {
  const double t = 301.24, density = 42.2, a0 = 1.0e-3;
  chem_rate p = basic_rate(CHEM_RATE_BRANCHED);
  p.p[13] = 1.2; p.p[14] = 204.3; p.p[15] = a0; p.p[16] = 2.0;
  const double k0 = 2.0e-22 * AVOGADRO * 1.0e-6 * exp(2.0);
  const double aref = k0 * (2.45e19 / AVOGADRO * 1.0e6);
  const double bref = 0.43 * pow(293.0 / 298.0, -8.0);
  const double avref = aref / (1.0 + aref / bref) * pow(0.41, 1.0 / (1.0 + pow(log10(aref / bref), 2.0)));
  const double z = avref * (1.0 - a0) / a0;
  const double a = k0 * density;
  const double b = 0.43 * pow(t / 298.0, -8.0);
  const double aval = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
  const double pre = 1.2 * exp(-204.3 / t);
  double value;
  p.branch = 1;
  if (evaluate(p, NULL, 0, NULL, 0, t, density, &value) || !close_to(value, pre * z / (z + aval))) return 0;
  p.branch = 0;
  return !evaluate(p, NULL, 0, NULL, 0, t, density, &value) && close_to(value, pre * aval / (z + aval));
}

static int falloff_cases(void) {
  const double t = 301.24, density = 42.2;
  chem_rate p = basic_rate(CHEM_RATE_TROE);
  double value;
  p.p[5] = p.p[8] = 1.0;
  if (evaluate(p, NULL, 0, NULL, 0, t, density, &value)) return 0;
  double ratio = density;
  double expected = density / (1.0 + ratio) * pow(0.6, 1.0 / (1.0 + pow(log10(ratio), 2.0)));
  if (!close_to(value, expected)) return 0;
  p.p[5] = 1.2; p.p[6] = 2.3; p.p[7] = 302.3;
  p.p[8] = 2.6; p.p[9] = -3.1; p.p[10] = 402.1; p.p[11] = 0.9; p.p[12] = 1.2;
  const double k0 = 1.2 * exp(302.3 / t) * pow(t / 300.0, 2.3);
  const double kinf = 2.6 * exp(402.1 / t) * pow(t / 300.0, -3.1);
  ratio = k0 * density / kinf;
  expected = density * k0 / (1.0 + ratio) * pow(0.9, 1.0 / (1.0 + 1.0 / 1.2 * pow(log10(ratio), 2.0)));
  if (evaluate(p, NULL, 0, NULL, 0, t, density, &value) || !close_to(value, expected)) return 0;
  p = basic_rate(CHEM_RATE_TROE);
  p.p[5] = 4.0e-18;
  if (evaluate(p, NULL, 0, NULL, 0, t, density, &value)) return 0;
  ratio = p.p[5] * density;
  expected = p.p[5] * density / (1.0 + ratio) * pow(0.6, 1.0 / (1.0 + pow(log10(ratio), 2.0)));
  if (!close_to(value, expected)) return 0;
  p.p[5] = 1.2e-12; p.p[6] = 167.0; p.p[7] = 3.0;
  p.p[8] = 136.0; p.p[9] = 5.0; p.p[10] = 24.0; p.p[11] = 0.9; p.p[12] = 0.8;
  const double k0b = p.p[5] * exp(3.0 / t) * pow(t / 300.0, 167.0);
  const double kinfb = 136.0 * exp(24.0 / t) * pow(t / 300.0, 5.0);
  ratio = k0b * density / kinfb;
  expected = k0b * density / (1.0 + ratio) * pow(0.9, 1.0 / (1.0 + 1.0 / 0.8 * pow(log10(ratio), 2.0)));
  return !evaluate(p, NULL, 0, NULL, 0, t, density, &value) && close_to(value, expected);
}

static int ternary_cases(void) {
  const double t = 301.24, density = 42.2;
  chem_rate p = basic_rate(CHEM_RATE_TERNARY);
  double value;
  if (evaluate(p, NULL, 0, NULL, 0, t, density, &value)) return 0;
  double ratio = density;
  double expected = 1.0 / (1.0 + ratio) * pow(0.6, 1.0 / (1.0 + pow(log10(ratio), 2.0)));
  if (!close_to(value, expected)) return 0;
  p.p[5] = 1.2; p.p[6] = 2.3; p.p[7] = 302.3;
  p.p[8] = 2.6; p.p[9] = -3.1; p.p[10] = 402.1; p.p[11] = 0.9; p.p[12] = 1.2;
  const double k0 = 1.2 * exp(302.3 / t) * pow(t / 300.0, 2.3);
  const double kinf = 2.6 * exp(402.1 / t) * pow(t / 300.0, -3.1);
  ratio = k0 * density / kinf;
  expected = k0 / (1.0 + ratio) * pow(0.9, 1.0 / (1.0 + 1.0 / 1.2 * pow(log10(ratio), 2.0)));
  if (evaluate(p, NULL, 0, NULL, 0, t, density, &value) || !close_to(value, expected)) return 0;

  p = basic_rate(CHEM_RATE_TERNARY_JPL19);
  p.p[5] = 1.2; p.p[6] = 2.3; p.p[7] = 302.3; p.p[18] = 298.0;
  p.p[8] = 2.6; p.p[9] = -3.1; p.p[10] = 402.1; p.p[19] = 310.0;
  p.p[20] = 0.7; p.p[21] = 0.5; p.p[22] = 120.0; p.p[23] = 290.0; p.p[11] = 0.9; p.p[12] = 1.2;
  const double k0j = 1.2 * exp(302.3 / t) * pow(t / 298.0, 2.3);
  const double kinfj = 2.6 * exp(402.1 / t) * pow(t / 310.0, -3.1);
  const double kint = 0.7 * exp(120.0 / t) * pow(t / 290.0, 0.5);
  ratio = k0j * density / kinfj;
  const double kf = k0j * density / (1.0 + ratio) * pow(0.9, 1.2 / (1.2 + pow(log10(ratio), 2.0)));
  expected = kf + kint * (1.0 - kf / kinfj);
  if (evaluate(p, NULL, 0, NULL, 0, t, density, &value) || !close_to(value, expected)) return 0;
  p = basic_rate(CHEM_RATE_TERNARY_JPL19);
  p.p[5] = 3.0e-3; p.p[6] = -2.0; p.p[8] = 1.5; p.p[9] = -1.0;
  p.p[18] = p.p[19] = 298.0; p.p[20] = 0.0; p.p[21] = 0.0; p.p[11] = 0.6; p.p[12] = 1.0;
  const double tj = 250.0, dj = 30.0;
  const double k0z = 3.0e-3 * pow(tj / 298.0, -2.0);
  const double kinfz = 1.5 * pow(tj / 298.0, -1.0);
  ratio = k0z * dj / kinfz;
  expected = k0z * dj / (1.0 + ratio) * pow(0.6, 1.0 / (1.0 + pow(log10(ratio), 2.0)));
  return !evaluate(p, NULL, 0, NULL, 0, tj, dj, &value) && close_to(value, expected);
}

static int remaining_rate_cases(void) {
  double value;
  chem_rate p = basic_rate(CHEM_RATE_TUNNELING);
  if (evaluate(p, NULL, 0, NULL, 0, 301.24, 1.0, &value) || !close_to(value, 1.0)) return 0;
  p.p[0] = 1.2; p.p[1] = 2.3; p.p[2] = 302.3;
  if (evaluate(p, NULL, 0, NULL, 0, 301.24, 1.0, &value) ||
      !close_to(value, 1.2 * exp(-2.3 / 301.24 + 302.3 / pow(301.24, 3.0)))) return 0;

  p = basic_rate(CHEM_RATE_TAYLOR);
  const double coeffs[] = { 12.5, 1.3e-2, 5.2e-4 };
  p.p[0] = 3.0e-11; p.p[2] = -200.0; p.coefficients = coeffs; p.coefficient_count = 3;
  double polynomial = coeffs[0] + coeffs[1] * 301.24 + coeffs[2] * pow(301.24, 2.0);
  if (evaluate(p, NULL, 0, NULL, 0, 301.24, 1.0, &value) ||
      !close_to(value, 3.0e-11 * exp(-200.0 / 301.24) * polynomial)) return 0;

  p = basic_rate(CHEM_RATE_REVERSIBLE);
  p.p[0] = 1.14e-2; p.p[2] = 2300.0; p.p[17] = 0.32;
  if (evaluate(p, NULL, 0, NULL, 0, 301.24, 1.0, &value) ||
      !close_to(value, 1.14e-2 * exp(2300.0 / 301.24) * 0.32)) return 0;
  if (evaluate(p, NULL, 0, NULL, 0, 298.15, 1.0, &value) ||
      !close_to(value, 1.14e-2 * exp(2300.0 / 298.15) * 0.32)) return 0;
  p.p[0] = 1.5; p.p[2] = 2000.0; p.p[17] = 0.5;
  const double temperatures[] = { 273.15, 298.15, 323.15 };
  for (size_t i = 0; i < sizeof(temperatures) / sizeof(temperatures[0]); ++i) {
    if (evaluate(p, NULL, 0, NULL, 0, temperatures[i], 1.0, &value) ||
        !close_to(value, 1.5 * exp(2000.0 / temperatures[i]) * 0.5)) return 0;
  }
  p.p[0] = 3.5; p.p[2] = 1800.0; p.p[17] = 1.0;
  if (evaluate(p, NULL, 0, NULL, 0, 298.15, 1.0, &value) ||
      !close_to(value, 3.5 * exp(1800.0 / 298.15))) return 0;
  p.p[17] = 0.0;
  if (evaluate(p, NULL, 0, NULL, 0, 298.15, 1.0, &value) || value != 0.0) return 0;

  p = basic_rate(CHEM_RATE_USER);
  p.custom_parameter = 0; p.p[1] = 2.0;
  const double user_params[] = { 0.6 };
  if (evaluate(p, NULL, 0, user_params, 1, 301.24, 1.0, &value) || value != 1.2) return 0;
  if (evaluate(p, NULL, 0, NULL, 0, 301.24, 1.0, &value) == 0) return 0;

  const chem_species species[] = { { "foo", 0.025, 230.0, 0.0, 1.0e-3, 0u } };
  p = basic_rate(CHEM_RATE_SURFACE);
  p.surface_species = 0; p.custom_parameter = 0;
  const double surface_params[] = { 1.0e-7, 2.5e6 };
  const double speed = sqrt(8.0 * GAS_CONSTANT * 273.65 / (PI * 0.025));
  double expected_surface = 4.0 * 2.5e6 * PI * pow(1.0e-7, 2.0) /
      (1.0e-7 / 230.0 + 4.0 / speed);
  if (evaluate(p, species, 1, surface_params, 2, 273.65, 1.0, &value) ||
      !close_to(value, expected_surface)) return 0;
  p.p[0] = 0.74;
  expected_surface = 4.0 * 2.5e6 * PI * pow(1.0e-7, 2.0) /
      (1.0e-7 / 230.0 + 4.0 / (speed * 0.74));
  if (evaluate(p, species, 1, surface_params, 2, 273.65, 1.0, &value) ||
      !close_to(value, expected_surface)) return 0;
  const chem_species missing_diffusion[] = { { "foo", 0.025, 0.0, 0.0, 1.0e-3, 0u } };
  const chem_species missing_weight[] = { { "foo", 0.0, 230.0, 0.0, 1.0e-3, 0u } };
  if (evaluate(p, missing_diffusion, 1, surface_params, 2, 273.65, 1.0, &value) != -2) return 0;
  return evaluate(p, missing_weight, 1, surface_params, 2, 273.65, 1.0, &value) == -2;
}

static int descriptor_cases(void) {
  const chem_mechanism *mechanism = chem_generated_mechanism();
  if (!mechanism || mechanism->species_count != 209 || mechanism->reaction_count != 547 ||
      mechanism->parameter_count != 168)
    return 0;
  size_t reactant_offset = 0, product_offset = 0;
  size_t user_rates = 0, surface_rates = 0, branched_rates = 0, density_terms = 0;
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    const chem_reaction *reaction = &mechanism->reactions[r];
    if (reaction->reactant_offset != reactant_offset || reaction->product_offset != product_offset)
      return 0;
    reactant_offset += reaction->reactant_count;
    product_offset += reaction->product_count;
    if (reaction->rate.kind == CHEM_RATE_USER) {
      if (reaction->rate.custom_parameter >= mechanism->parameter_count) return 0;
      ++user_rates;
    } else if (reaction->rate.kind == CHEM_RATE_SURFACE) {
      if (reaction->rate.custom_parameter + 1u >= mechanism->parameter_count ||
          reaction->rate.surface_species >= mechanism->species_count)
        return 0;
      ++surface_rates;
    } else if (reaction->rate.kind == CHEM_RATE_BRANCHED) {
      ++branched_rates;
    }
    density_terms += reaction->rate.air_density_power;
  }
  return reactant_offset == mechanism->reactant_term_count &&
      product_offset == mechanism->product_term_count && user_rates == 142 &&
      surface_rates == 13 && branched_rates == 0 && density_terms > 0;
}

int main(void) {
  const int results[] = { arrhenius_cases(), branched_cases(), falloff_cases(), ternary_cases(),
                          remaining_rate_cases(), descriptor_cases() };
  const char *names[] = { "arrhenius", "branched", "troe", "ternary", "other rate constants", "mechanism descriptor" };
  int failed = 0;
  for (size_t i = 0; i < sizeof(results) / sizeof(results[0]); ++i) {
    if (!results[i]) {
      fprintf(stderr, "rate case group failed: %s\n", names[i]);
      failed = 1;
    }
  }
  return failed;
}
