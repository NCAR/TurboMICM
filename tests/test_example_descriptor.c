#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef EXPECT_USER_DEFINED
#define EXPECT_USER_DEFINED 0
#endif

static int parameter_index(const chem_mechanism *mechanism, const char *name) {
  for (size_t i = 0; i < mechanism->parameter_count; ++i)
    if (strcmp(mechanism->parameters[i].name, name) == 0) return (int)i;
  return -1;
}

static int close_to(double a, double b) {
  return isfinite(a) && fabs(a - b) <= 1.0e-13 * fabs(b) + 1.0e-30;
}

static chem_rate make_rate(chem_rate_kind kind) {
  chem_rate rate;
  memset(&rate, 0, sizeof(rate));
  rate.kind = kind;
  rate.p[3] = 300.0;
  rate.p[5] = 1.0;
  rate.p[8] = 1.0;
  rate.p[11] = 0.6;
  rate.p[12] = 1.0;
  rate.p[18] = rate.p[19] = rate.p[23] = 298.0;
  rate.p[20] = 1.0;
  return rate;
}

static void set_reaction(chem_reaction *reaction, size_t roff, size_t rc, size_t poff, size_t pc,
                         chem_rate rate) {
  reaction->reactant_offset = (unsigned)roff;
  reaction->reactant_count = (unsigned)rc;
  reaction->product_offset = (unsigned)poff;
  reaction->product_count = (unsigned)pc;
  reaction->rate = rate;
}

static int matches_hand_built_result(const chem_mechanism *generated) {
  chem_species species[7] = {
    { "A", 0.0, 0.0, 0.0, 1.0e-3, 0u }, { "B", 0.0, 0.0, 0.0, 1.0e-3, 0u },
    { "C", 0.025, 230.0, 0.0, 1.0e-3, 0u }, { "D", 0.0, 0.0, 0.0, 1.0e-3, 0u },
    { "E", 0.0, 0.0, 0.0, 1.0e-3, 0u }, { "F", 0.0, 0.0, 0.0, 1.0e-3, 0u },
    { "G", 0.0, 0.0, 0.0, 1.0e-3, 0u }
  };
  chem_term reactants[12];
  chem_term products[10];
  chem_reaction reactions[10];
  size_t roff = 0, poff = 0, r = 0;
#define REACTANT(id) (reactants[roff++] = (chem_term){ (id), 1.0 })
#define PRODUCT(id, coeff) (products[poff++] = (chem_term){ (id), (coeff) })
#define START_REACTION(rcount, pcount, kind) \
  set_reaction(&reactions[r++], roff, (rcount), poff, (pcount), make_rate(kind))

  START_REACTION(1, 1, CHEM_RATE_ARRHENIUS); REACTANT(0u); PRODUCT(1u, 1.0);
  reactions[0].rate.p[0] = 2.15e-4; reactions[0].rate.p[2] = 110.0;
  START_REACTION(1, 1, CHEM_RATE_BRANCHED); REACTANT(1u); PRODUCT(2u, 1.0);
  reactions[1].rate.p[13] = 1.2; reactions[1].rate.p[14] = 204.3; reactions[1].rate.p[15] = 1.0e-3;
  reactions[1].rate.p[16] = 2.0; reactions[1].rate.branch = 1u;
  START_REACTION(1, 1, CHEM_RATE_BRANCHED); REACTANT(1u); PRODUCT(3u, 1.0);
  reactions[2].rate.p[13] = 1.2; reactions[2].rate.p[14] = 204.3; reactions[2].rate.p[15] = 1.0e-3;
  reactions[2].rate.p[16] = 2.0; reactions[2].rate.branch = 0u;
  START_REACTION(1, 1, CHEM_RATE_SURFACE); REACTANT(2u); PRODUCT(4u, 1.0);
  reactions[3].rate.p[0] = 0.90; reactions[3].rate.surface_species = 2u;
  reactions[3].rate.custom_parameter = 0u;
  START_REACTION(1, 1, CHEM_RATE_TERNARY); REACTANT(3u); PRODUCT(5u, 2.0);
  reactions[4].rate.p[5] = 1.2; reactions[4].rate.p[6] = 2.3; reactions[4].rate.p[7] = 302.3;
  reactions[4].rate.p[8] = 2.6 / (6.02214076e23 * 1.0e-6); reactions[4].rate.p[9] = -3.1;
  reactions[4].rate.p[10] = 402.1; reactions[4].rate.p[11] = 0.9; reactions[4].rate.p[12] = 1.2;
  START_REACTION(2, 1, CHEM_RATE_TROE); REACTANT(4u); REACTANT(4u); PRODUCT(6u, 1.0);
  reactions[5].rate.p[5] = 1.8e-30 * pow(6.02214076e17, 2.0); reactions[5].rate.p[6] = -3.0;
  reactions[5].rate.p[8] = 2.8e-11 * 6.02214076e17; reactions[5].rate.p[11] = 0.6;
  reactions[5].rate.p[12] = 1.0;
  START_REACTION(1, 1, CHEM_RATE_TUNNELING); REACTANT(5u); PRODUCT(6u, 1.0);
  reactions[6].rate.p[0] = 1.2; reactions[6].rate.p[1] = 2.3; reactions[6].rate.p[2] = 302.3;
#if EXPECT_USER_DEFINED
  START_REACTION(1, 1, CHEM_RATE_USER); REACTANT(2u); PRODUCT(6u, 1.0);
  reactions[7].rate.label = "my photolysis rate"; reactions[7].rate.custom_parameter = 2u; reactions[7].rate.p[1] = 1.0;
  START_REACTION(1, 0, CHEM_RATE_USER); REACTANT(1u);
  reactions[8].rate.label = "my loss rate"; reactions[8].rate.custom_parameter = 3u; reactions[8].rate.p[1] = 1.0;
  START_REACTION(0, 1, CHEM_RATE_USER); PRODUCT(0u, 1.0);
  reactions[9].rate.label = "my emission rate"; reactions[9].rate.custom_parameter = 4u; reactions[9].rate.p[1] = 1.0;
#endif

  chem_mechanism hand = *generated;
  hand.species = species;
  hand.reactions = reactions;
  hand.reactant_terms = reactants;
  hand.product_terms = products;
  hand.reactant_term_count = roff;
  hand.product_term_count = poff;
  hand.temperature = 287.45;
  hand.pressure = 101319.9;
  hand.relative_tolerance = 1.0e-6;

  double generated_state[7] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double hand_state[7];
  memcpy(hand_state, generated_state, sizeof(hand_state));
  double generated_parameters[8] = { 0.0 }, hand_parameters[8] = { 0.0 };
  for (size_t p = 0; p < generated->parameter_count; ++p) {
    const char *name = generated->parameters[p].name;
    const double value = strstr(name, "effective radius") ? 1.0e-7 :
        (strstr(name, "particle number") ? 2.5e6 :
         (strcmp(name, "my photolysis rate") == 0 ? 1.0e-10 :
          (strcmp(name, "my emission rate") == 0 ? 1.0e-20 :
           (strcmp(name, "my loss rate") == 0 ? 1.0e-23 : 0.0))));
    generated_parameters[p] = value;
    hand_parameters[p] = value;
  }
  const double temperature = 287.45, pressure = 101319.9;
  for (size_t i = 0; i < generated->reaction_count; ++i) {
    double configured_rate = 0.0, hand_rate = 0.0;
    if (chem_evaluate_rate(generated, i, generated_parameters, temperature, pressure, &configured_rate) ||
        chem_evaluate_rate(&hand, i, hand_parameters, temperature, pressure, &hand_rate) ||
        !close_to(configured_rate, hand_rate)) {
      fprintf(stderr, "example rate %zu mismatch: %.17g vs %.17g\n", i, configured_rate, hand_rate);
      return 0;
    }
  }
  for (size_t step = 0; step < 10; ++step) {
#if EXPECT_USER_DEFINED
    const int photo = parameter_index(generated, "my photolysis rate");
    generated_parameters[photo] = hand_parameters[photo] = 1.0e-10 * pow(1.5, (double)step);
#endif
    if (chem_integrate_group(generated, 1, generated_state, generated_parameters,
                             &temperature, &pressure, 500.0, NULL) != 0 ||
        chem_integrate_group(&hand, 1, hand_state, hand_parameters,
                             &temperature, &pressure, 500.0, NULL) != 0)
      return 0;
  }
  for (size_t i = 0; i < 7; ++i) {
    if (!close_to(generated_state[i], hand_state[i])) {
      fprintf(stderr, "example state mismatch %s: %.17g vs %.17g\n",
              generated->species[i].name, generated_state[i], hand_state[i]);
      return 0;
    }
  }
  return 1;
#undef REACTANT
#undef PRODUCT
#undef START_REACTION
}

int main(void) {
  const chem_mechanism *mechanism = chem_generated_mechanism();
  const size_t expected_reactions = EXPECT_USER_DEFINED ? 10u : 7u;
  const size_t expected_parameters = EXPECT_USER_DEFINED ? 5u : 2u;
  if (!mechanism || mechanism->species_count != 7 || mechanism->reaction_count != expected_reactions ||
      mechanism->parameter_count != expected_parameters) {
    fputs("configuration descriptor dimensions differ from the hand-built mechanism\n", stderr);
    return 1;
  }
  size_t counts[CHEM_RATE_SURFACE + 1] = { 0 };
  size_t reactant_offset = 0, product_offset = 0;
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    const chem_reaction *reaction = &mechanism->reactions[r];
    if (reaction->reactant_offset != reactant_offset || reaction->product_offset != product_offset ||
        reaction->rate.kind > CHEM_RATE_SURFACE) return 1;
    ++counts[reaction->rate.kind];
    reactant_offset += reaction->reactant_count;
    product_offset += reaction->product_count;
  }
  if (counts[CHEM_RATE_ARRHENIUS] != 1 || counts[CHEM_RATE_BRANCHED] != 2 ||
      counts[CHEM_RATE_SURFACE] != 1 || counts[CHEM_RATE_TERNARY] != 1 ||
      counts[CHEM_RATE_TROE] != 1 || counts[CHEM_RATE_TUNNELING] != 1 ||
      counts[CHEM_RATE_USER] != (EXPECT_USER_DEFINED ? 3u : 0u)) return 1;

  double parameters[8] = { 0.0 };
  const int radius = parameter_index(mechanism, "C surface.effective radius [m]");
  const int number = parameter_index(mechanism, "C surface.particle number concentration [# m-3]");
  if (radius < 0 || number < 0) return 1;
  parameters[radius] = 1.0e-7;
  parameters[number] = 2.5e6;
#if EXPECT_USER_DEFINED
  const int photo = parameter_index(mechanism, "my photolysis rate");
  const int loss = parameter_index(mechanism, "my loss rate");
  const int emission = parameter_index(mechanism, "my emission rate");
  if (photo < 0 || loss < 0 || emission < 0) return 1;
  parameters[photo] = 1.0e-10;
  parameters[loss] = 1.0e-23;
  parameters[emission] = 1.0e-20;
#endif
  const double temperature = 287.45;
  const double pressure = 101319.9;
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    double value = 0.0;
    if (chem_evaluate_rate(mechanism, r, parameters, temperature, pressure, &value) ||
        !isfinite(value) || value < 0.0) return 1;
    const chem_reaction *reaction = &mechanism->reactions[r];
    if (reaction->rate.kind == CHEM_RATE_USER) {
      const double expected = parameters[reaction->rate.custom_parameter] * reaction->rate.p[1];
      if (!close_to(value, expected)) return 1;
    }
    if (reaction->rate.kind == CHEM_RATE_SURFACE && value <= 0.0) return 1;
  }
  if (!matches_hand_built_result(mechanism)) {
    fputs("configuration and hand-built mechanism integrations differ\n", stderr);
    return 1;
  }
  puts("configuration descriptor matches expected by-hand structure and evaluates all rate types");
  return 0;
}
