#include "test.h"
int test_initialization(const char *suite, const char *test)
{
  if (strcmp(suite, "ConstraintInitialization"))
    return 2;
  if (IS("SolverStateToStringNewValue"))
  {
    CHECK(!strcmp(ck_status_string(CK_CONSTRAINT_FAILED), "Constraint Initialization Failed"));
    return 0;
  }
  const ck_species species[] = {{.name = "A"}, {.name = "B"}, {.name = "C"}};
  const size_t r[] = {0};
  const ck_product product[] = {{1, 1}};
  ck_reaction rx = {.nreactants = 1,
                    .reactants = r,
                    .nproducts = 1,
                    .products = product,
                    .rate = ck_rate_defaults(CK_ARRHENIUS)};
  rx.rate.A = 0.5;
  ck_product reactants[] = {{1, 1}, {0, 1}, {2, 1}}, products[] = {{2, 1}};
  ck_constraint constraint = {.type = CK_EQUILIBRIUM_CONSTRAINT,
                              .name = "B_C_eq",
                              .algebraic_species = 2,
                              .nreactants = 1,
                              .reactants = reactants,
                              .nproducts = 1,
                              .products = products,
                              .equilibrium_ref = 2};
  ck_mechanism mechanism = {.name = "initialization",
                            .nspecies = 3,
                            .species = species,
                            .nreactions = 1,
                            .reactions = &rx,
                            .nconstraints = 1,
                            .constraints = &constraint,
                            .temperature = 298.15,
                            .pressure = 101325};
  ck_rosenbrock_parameters p;
  CHECK(ck_rosenbrock_parameters_init(3, &p) == CK_OK);
  int squared =
      IS("FailedInitializationRestoresCallerState") || IS("SuccessfulInitializationKeepsTheCorrection");
  int scale = IS("WeightedCorrectionIsInvariantToConstraintRowScaling") ||
              IS("WeightedCorrectionUsesStateTolerances") ||
              IS("ScaledConstraintRowsSolveIdenticallyThroughSolve");
  int through = IS("ScaledConstraintRowsSolveIdenticallyThroughSolve");
  if (squared)
  {
    products[0].yield = 2;
    mechanism.nreactions = 0;
  }
  if (scale)
  {
    constraint.type = CK_LINEAR_CONSTRAINT;
    constraint.nproducts = 0;
    constraint.nreactants = 2;
    reactants[0] = (ck_product){0, 1};
    reactants[1] = (ck_product){2, -1};
    mechanism.nreactions = 0;
  }
  if (through)
  {
    constraint.nreactants = 3;
    reactants[0] = (ck_product){0, 1};
    reactants[1] = (ck_product){1, 1};
    reactants[2] = (ck_product){2, 1};
    mechanism.nreactions = 1;
    CHECK(ck_rosenbrock_parameters_init(43, &p) == CK_OK);
  }
  if (IS("PureODESystemUnaffected"))
    mechanism.nconstraints = 0;
  size_t runs = IS("WeightedCorrectionUsesStateTolerances") ? 2 : scale ? 3 : 1;
  double scales[] = {1, 1e-12, 1e12}, reference = 0;
  for (size_t run = 0; run < runs; ++run)
  {
    if (scale && !IS("WeightedCorrectionUsesStateTolerances"))
    {
      reactants[0].yield = scales[run];
      reactants[1].yield = through ? scales[run] : -scales[run];
      if (through)
      {
        reactants[2].yield = scales[run];
        constraint.constant = scales[run] * 1.3;
      }
    }
    ck_model *m;
    ck_state *s;
    CHECK(ck_model_create_with_order(&mechanism, 0, &m) == CK_OK);
    size_t cells = IS("MultiCellSystems") ? 3 : 1;
    CHECK(ck_state_create(m, cells, &s) == CK_OK);
    s->y[0] = 1;
    s->y[cells] = 0.5;
    s->y[2 * cells] = 1;
    if (IS("MildlyInconsistentICsCorrected"))
      s->y[2 * cells] = 5;
    if (IS("SeverelyInconsistentICsConverge"))
      s->y[2 * cells] = 1000;
    if (cells == 3)
    {
      s->y[1] = 2;
      s->y[2] = 0.5;
      s->y[4] = 0.3;
      s->y[5] = 0.1;
      s->y[7] = 5;
      s->y[8] = 100;
    }
    if (squared)
    {
      s->y[2] = IS("FailedInitializationRestoresCallerState") ? 100 : 1e-6;
      if (IS("FailedInitializationRestoresCallerState"))
        p.constraint_init_max_iterations = 1;
    }
    if (scale)
    {
      s->y[0] = 1;
      s->y[2] = 4;
      double atol[] = {1e-12, 1e-12, 1e-12};
      CHECK(ck_state_tolerances(s, 1e-8, atol, 3) == CK_OK);
    }
    if (through)
    {
      s->y[0] = 0.7;
      s->y[1] = 0.15;
      s->y[2] = 5;
    }
    if (IS("WeightedCorrectionUsesStateTolerances"))
    {
      double atol[] = {run ? 1e-12 : 1e-8, run ? 1e-12 : 1e-8, run ? 1e-12 : 1e-8};
      CHECK(ck_state_tolerances(s, run ? 1e-8 : 1e-6, atol, 3) == CK_OK);
      s->y[2] = 1 + 5e-8;
      NEAR(p.constraint_init_tolerance, 0.1, 0);
    }
    CHECK(ck_update_rates(s) == CK_OK);
    ck_stats stats = {0};
    ck_status status = squared || (scale && !through) ? ck_initialize_constraints(s, &p, &stats)
                                                      : ck_solve(s, 0.001, &p, &stats);
    if (IS("FailedInitializationRestoresCallerState"))
    {
      CHECK(status == CK_CONSTRAINT_FAILED);
      NEAR(s->y[1], 0.5, 0);
      NEAR(s->y[2], 100, 0);
    }
    else
    {
      CHECK(status == CK_OK);
      if (IS("ConsistentICsUnchanged"))
        CHECK(stats.constraint_init_iterations <= 1);
      if (IS("MildlyInconsistentICsCorrected"))
        CHECK(stats.constraint_init_iterations > 0);
      if (IS("PureODESystemUnaffected"))
        CHECK(stats.constraint_init_iterations == 0);
      else if (scale)
      {
        double expected = through ? 1.3 - 0.7 - 0.15 : 1;
        if (IS("WeightedCorrectionUsesStateTolerances"))
        {
          if (run == 0)
          {
            NEAR(s->y[2], 1 + 5e-8, 0);
            CHECK(stats.constraint_init_iterations == 1);
          }
          else
          {
            NEAR(s->y[2], 1, 1e-12);
            CHECK(stats.constraint_init_iterations == 2);
          }
        }
        else
        {
          NEAR(s->y[2], expected, through ? 1e-5 * expected : 1e-12);
          if (run == 0)
            reference = s->y[2];
          else if (through)
            NEAR(s->y[2], reference, 1e-5 * expected);
        }
      }
      else if (squared)
      {
        NEAR(s->y[2], 1, 1e-3);
        const size_t backtracks[] = {0, 1};
        for (size_t b = 0; b < 2; ++b)
        {
          p.constraint_init_max_backtracks = backtracks[b];
          s->y[2] = 1e-6;
          CHECK(ck_initialize_constraints(s, &p, &stats) == CK_CONSTRAINT_FAILED);
          NEAR(s->y[1], 0.5, 0);
          NEAR(s->y[2], 1e-6, 0);
        }
      }
      else
        for (size_t c = 0; c < cells; ++c)
          NEAR(2 * s->y[cells + c] - s->y[2 * cells + c], 0, cells == 3 ? 1e-5 : 1e-6);
      if (IS("SubsequentSolveCallsReinitialize"))
      {
        s->y[2] = 999;
        CHECK(ck_solve(s, 0.001, &p, &stats) == CK_OK);
        NEAR(2 * s->y[1] - s->y[2], 0, 1e-6);
      }
    }
    ck_state_free(s);
    ck_model_free(m);
  }
  return 0;
}
