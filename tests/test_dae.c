#include "test.h"
#include <float.h>

static int equilibrium_integration(const char *test)
{
  int api = IS("SetConstraintsAPIWorks"), multiple = IS("SetConstraintsAPIMultipleConstraints"),
      coupled = IS("DAESolveWithTwoCoupledConstraints"), nonunit = IS("DAESolveWithNonUnitStoichiometry"),
      reorder = IS("DAESolveWithConstraintAndReorderState");
  size_t n = multiple ? 6 : coupled ? 4 : 3, nc = multiple || coupled ? 2 : 1;
  ck_species sp[] = {{.name = "A"}, {.name = "B"}, {.name = "C"},
                     {.name = "D"}, {.name = "E"}, {.name = "F"}};
  size_t reactants[] = {nonunit ? 2 : 0, 3};
  ck_product products[] = {{nonunit ? 0 : 1, 1}, {4, 1}};
  ck_reaction rx[2];
  for (size_t r = 0; r < 2; ++r)
  {
    rx[r] = (ck_reaction){.nreactants = 1,
                          .reactants = reactants + r,
                          .nproducts = 1,
                          .products = products + r,
                          .rate = ck_rate_defaults(CK_ARRHENIUS)};
    rx[r].rate.A = r ? .2 : api ? .1 : multiple || nonunit ? .5 : 1;
  }
  ck_product cr[] = {{nonunit ? 0 : 1, nonunit ? 2 : 1}, {multiple ? 4 : 1, 1}},
             cp[] = {{nonunit ? 1 : 2, 1}, {multiple ? 5 : 3, 1}};
  ck_constraint constraints[2];
  double eq[] = {api || multiple ? .034 : coupled ? 3 : nonunit ? 10 : 2, multiple ? .012 : 5};
  for (size_t k = 0; k < 2; ++k)
    constraints[k] = (ck_constraint){.type = CK_EQUILIBRIUM_CONSTRAINT,
                                     .name = k         ? (multiple ? "E_F_eq" : "B_D_eq")
                                             : nonunit ? "A2_B_eq"
                                                       : "B_C_eq",
                                     .algebraic_species = cp[k].species,
                                     .nreactants = 1,
                                     .reactants = cr + k,
                                     .nproducts = 1,
                                     .products = cp + k,
                                     .equilibrium_ref = eq[k],
                                     .enthalpy = k && multiple ? -2000 : -2400};
  ck_mechanism d = {.name = "DAE",
                    .nspecies = n,
                    .species = sp,
                    .nreactions = multiple ? 2 : 1,
                    .reactions = rx,
                    .nconstraints = nc,
                    .constraints = constraints,
                    .temperature = api        ? 300
                                   : multiple ? 310
                                   : coupled  ? 300
                                   : nonunit  ? 298
                                   : reorder  ? 400
                                              : 270,
                    .pressure = 101325};
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create_with_order(&d, reorder, &m) == CK_OK);
  CHECK(ck_state_create(m, 1, &s) == CK_OK);
  CHECK(ck_model_nspecies(m) == n && ck_model_nconstraints(m) == nc);
  size_t index[6], parameter[2];
  for (size_t i = 0; i < n; ++i)
  {
    CHECK(ck_species_index(m, sp[i].name, index + i) == CK_OK);
    CHECK(ck_is_algebraic(m, index[i]) == (i == (nonunit ? 1u : 2u) || (nc == 2 && i == cp[1].species)));
  }
  for (size_t k = 0; k < nc; ++k)
    CHECK(ck_custom_index(m, constraints[k].name, parameter + k) == CK_OK);
  if (nonunit)
  {
    s->y[index[0]] = .1;
    s->y[index[1]] = .1;
    s->y[index[2]] = 1;
  }
  else
  {
    s->y[index[0]] = 1;
    if (api || multiple)
    {
      s->y[index[1]] = .1;
      s->y[index[2]] = eq[0] * .1;
      if (multiple)
      {
        s->y[index[3]] = .5;
        s->y[index[4]] = .05;
        s->y[index[5]] = eq[1] * .05;
      }
    }
  }
  CHECK(ck_update_rates(s) == CK_OK);
  for (size_t k = 0; k < nc; ++k)
  {
    double expected =
        eq[k] * exp(constraints[k].enthalpy / CK_GAS_CONSTANT * (1 / d.temperature - 1 / 298.15));
    NEAR(s->custom[parameter[k]], expected, fabs(expected) * 100 * DBL_EPSILON);
  }
  if (!api && !multiple)
  {
    ck_rosenbrock_parameters p;
    CHECK(ck_rosenbrock_parameters_init(3, &p) == CK_OK);
    for (double time = 0; time < (nonunit ? .05 : .1); time += .001)
    {
      CHECK(ck_update_rates(s) == CK_OK);
      ck_stats stats;
      CHECK(ck_solve(s, .001, &p, &stats) == CK_OK);
      for (size_t k = 0; k < nc; ++k)
      {
        double residual = s->custom[parameter[k]] * pow(s->y[index[cr[k].species]], cr[k].yield) -
                          s->y[index[cp[k].species]];
        NEAR(residual, 0, nonunit ? 1e-5 : 1e-6);
      }
    }
    if (!nonunit)
      NEAR(s->y[index[0]] + s->y[index[1]], 1, .01);
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
static int budget_system(int method, int equilibrium, int cascade, double balance_atol, double rate,
                         double K1, double K2, ck_stats *total, double *result)
{
  size_t n = cascade ? 4 : 3;
  ck_species sp[] = {
      {.name = "gas"}, {.name = "aqueous"}, {.name = cascade ? "aqueous2" : "P"}, {.name = "P"}};
  size_t reactant = cascade ? 2 : equilibrium ? 1 : 0;
  ck_product prod = {cascade ? 3 : equilibrium ? 2 : 1, 1};
  ck_reaction rx = {.nreactants = 1,
                    .reactants = &reactant,
                    .nproducts = 1,
                    .products = &prod,
                    .rate = ck_rate_defaults(CK_ARRHENIUS)};
  rx.rate.A = rate;
  ck_product eqr[] = {{0, 1}, {1, 1}}, eqp[] = {{1, 1}, {2, 1}},
             mass[] = {{1, 1}, {cascade ? 2 : 2, 1}, {cascade ? 3 : 0, 1}, {0, 1}};
  if (!equilibrium && !cascade)
  {
    mass[0].species = 0;
    mass[1].species = 1;
    mass[2].species = 2;
  }
  ck_constraint c[3];
  size_t nc = 0;
  if (equilibrium || cascade)
    c[nc++] = (ck_constraint){.type = CK_EQUILIBRIUM_CONSTRAINT,
                              .name = "eq1",
                              .algebraic_species = 1,
                              .nreactants = 1,
                              .reactants = eqr,
                              .nproducts = 1,
                              .products = eqp,
                              .equilibrium_ref = K1};
  if (cascade)
    c[nc++] = (ck_constraint){.type = CK_EQUILIBRIUM_CONSTRAINT,
                              .name = "eq2",
                              .algebraic_species = 2,
                              .nreactants = 1,
                              .reactants = eqr + 1,
                              .nproducts = 1,
                              .products = eqp + 1,
                              .equilibrium_ref = K2};
  c[nc++] = (ck_constraint){.type = CK_LINEAR_CONSTRAINT,
                            .name = "mass",
                            .algebraic_species = equilibrium || cascade ? 0 : 2,
                            .nreactants = n,
                            .reactants = mass,
                            .constant = 1e-6};
  ck_mechanism d = {.name = "conservation",
                    .nspecies = n,
                    .species = sp,
                    .nreactions = 1,
                    .reactions = &rx,
                    .nconstraints = nc,
                    .constraints = c,
                    .temperature = 298,
                    .pressure = 101325};
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  CHECK(ck_state_create(m, 1, &s) == CK_OK);
  s->relative_tolerance = 1e-6;
  for (size_t j = 0; j < n; ++j)
    s->absolute_tolerance[j] = 1e-12;
  if (cascade)
  {
    s->absolute_tolerance[0] = balance_atol;
    s->absolute_tolerance[1] = s->absolute_tolerance[2] = 1e-8;
    s->y[0] = 1e-6 / (1 + K1 + K1 * K2);
    s->y[1] = K1 * s->y[0];
    s->y[2] = K2 * s->y[1];
  }
  else if (equilibrium)
  {
    s->absolute_tolerance[0] = s->absolute_tolerance[1] = 1e-8;
    s->y[0] = 1e-6 / (1 + K1);
    s->y[1] = K1 * s->y[0];
  }
  else
  {
    s->y[0] = .9e-6;
    s->y[2] = .1e-6;
  }
  CHECK(ck_update_rates(s) == CK_OK);
  ck_rosenbrock_parameters p;
  CHECK(ck_rosenbrock_parameters_init(method, &p) == CK_OK);
  memset(total, 0, sizeof *total);
  double advanced = 0, min_gas = 1e99;
  while (advanced < 30)
  {
    ck_stats stats;
    CHECK(ck_solve(s, 30 - advanced, &p, &stats) == CK_OK);
    CHECK(stats.final_time > 0);
    advanced += stats.final_time;
    ck_stats_add(total, &stats);
    min_gas = fmin(min_gas, s->y[0]);
    double sum = 0;
    for (size_t j = 0; j < n; ++j)
      sum += s->y[j];
    if (!cascade)
      NEAR(sum, 1e-6, 1e-12);
    if (!cascade && equilibrium)
    {
      CHECK(s->y[0] >= -1e-18);
      CHECK(s->y[1] >= -1e-18);
    }
    if (!cascade && !equilibrium)
      CHECK(s->y[2] >= -1e-18);
  }
  if (!cascade)
  {
    if (equilibrium)
      NEAR(s->y[2], 1e-6, 1e-8);
    else
    {
      CHECK(s->y[0] < 1e-12);
      NEAR(s->y[1], .9e-6, 1e-10);
      NEAR(s->y[2], .1e-6, 1e-10);
    }
  }
  if (result)
  {
    result[0] = s->y[0];
    result[1] = s->y[n - 1];
    result[2] = min_gas;
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
int test_dae(const char *suite, const char *test)
{
  if (!strcmp(suite, "EquilibriumIntegration"))
    return equilibrium_integration(test);
  ck_stats stats;
  if (!strcmp(suite, "DAEConstraintOvershoot"))
  {
    if (IS("AllRosenbrockOrdersConstrained"))
    {
      int methods[] = {3, 4, 43, 6};
      for (size_t k = 0; k < 4; ++k)
        CHECK(budget_system(methods[k], 0, 0, 0, 1e4, 0, 0, &stats, NULL) == 0);
      return 0;
    }
    return budget_system(43, IS("EquilibriumPlusConservation"), 0, 0,
                         IS("EquilibriumPlusConservation") ? 1e3 : 1e4, 10, 0, &stats, NULL);
  }
  if (!strcmp(suite, "DAEAlgebraicError"))
  {
    if (IS("ErrorSensitiveToBalanceAtol"))
    {
      ck_stats tight;
      CHECK(budget_system(43, 0, 1, 1e-3, 100, 2, 2, &stats, NULL) == 0);
      CHECK(budget_system(43, 0, 1, 1e-8, 100, 2, 2, &tight, NULL) == 0);
      CHECK(tight.accepted > stats.accepted);
      return 0;
    }
    double result[3];
    CHECK(budget_system(43, 0, 1, 1e-8, 1e4, 100, 50, &stats, result) == 0);
    CHECK(result[2] >= -1e-8);
    NEAR(result[1], 1e-6, 1e-8);
    NEAR(result[0], 0, 1e-10);
    return 0;
  }
  return 2;
}
