#include "chemistry.h"

#include <math.h>

int main(void) {
  const double t = 301.24;
  const double a = chem_arrhenius(1.0, 0.0, 0.0, 300.0, 0.0, t, 0.0);
  const double b = chem_arrhenius(3.0e-11, 0.0, -200.0, 300.0, 0.0, t, 0.0);
  const double e1 = 1.0;
  const double e2 = 3.0e-11 * exp(-200.0 / t);
  if (fabs(a - e1) > 1e-14 * e1) return 1;
  if (fabs(b - e2) > 1e-14 * e2) return 2;

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
  return 0;
}
