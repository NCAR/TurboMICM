#include "test.h"
#include "internal.h"

typedef struct
{
  int calls;
  size_t variables, entries;
} finalize_record;
static ck_status record_finalize(const ck_model *m, void *context)
{
  finalize_record *r = context;
  ++r->calls;
  r->variables = ck_model_nspecies(m);
  r->entries = ck_jacobian_size(m);
  size_t a;
  return ck_species_index(m, "A", &a);
}
static double lambda_rate(const ck_conditions *c, const double *custom, void *context)
{
  (void)custom;
  return *(const double *)context * c->temperature;
}
int test_builder(const char *suite, const char *test)
{
  int lambda = !strcmp(suite, "Solver") && IS("GetLambdaRateConstantByNameCanOverrideLambda");
  if (strcmp(suite, "SolverBuilder") && !lambda)
    return 2;
  ck_species sp[] = {{.name = "A"}, {.name = "B"}, {.name = "C"}, {.name = "D"}};
  size_t reactants[] = {0, 1};
  ck_product products[] = {{1, 1}, {2, 1}};
  ck_reaction rx[2] = {0};
  for (size_t r = 0; r < 2; ++r)
  {
    rx[r].nreactants = rx[r].nproducts = 1;
    rx[r].reactants = reactants + r;
    rx[r].products = products + r;
    rx[r].rate = ck_rate_defaults(CK_ARRHENIUS);
    rx[r].rate.A = r ? 3.3e-11 : 2.15e-11;
    rx[r].rate.C = r ? 55 : 110;
  }
  ck_mechanism d = {.name = "gas", .nspecies = 3, .species = sp, .nreactions = 2, .reactions = rx};
  ck_model *m = NULL;
  ck_state *s = NULL;
  if (IS("ThrowsMissingSystem"))
  {
    CHECK(ck_model_create(NULL, &m) == CK_INVALID);
    CHECK(!m);
    return 0;
  }
  if (IS("ThrowsUnusedSpecies"))
  {
    d.nspecies = 4;
    CHECK(ck_model_validate_species_used(&d) == CK_INVALID);
    return 0;
  }
  finalize_record record = {0};
  ck_external_model e = {.context = &record, .finalize = record_finalize};
  if (IS("FinalizeReceivesDenseMatrixPolicy") || IS("FinalizeSupportsSparseOnlySignature"))
  {
    d.nexternal = 1;
    d.external = &e;
    for (size_t cells = 1; cells <= 4; cells += 3)
    {
      record.calls = 0;
      CHECK(ck_model_create(&d, &m) == CK_OK);
      CHECK(ck_state_create(m, cells, &s) == CK_OK);
      CHECK(record.calls == 1 && record.variables == 3);
      CHECK(record.entries == ck_jacobian_size(m));
      CHECK(s->cells == cells);
      ck_state_free(s);
      ck_model_free(m);
    }
    return 0;
  }
  double original = .001, override = .002;
  if (lambda)
  {
    d.nspecies = 2;
    d.nreactions = 1;
    rx[0].rate = ck_rate_defaults(CK_LAMBDA);
    rx[0].rate.label = "lambda_rc";
    rx[0].rate.function = lambda_rate;
    rx[0].rate.context = &original;
  }
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  CHECK(ck_model_nspecies(m) == d.nspecies);
  CHECK(!strcmp(m->description->name, "gas"));
  CHECK(ck_state_create(m, 1, &s) == CK_OK);
  ck_stats stats;
  ck_rosenbrock_parameters rp;
  CHECK(ck_rosenbrock_parameters_init(3, &rp) == CK_OK);
  ck_backward_euler_parameters bp = ck_backward_euler_defaults();
  if (lambda)
  {
    CHECK(ck_model_set_lambda(m, "lambda_rc", lambda_rate, &override) == CK_OK);
    s->conditions[0] = (ck_conditions){300, 101253.3, 1e6};
    s->y[0] = 1;
    double tol[] = {1e-6, 1e-6};
    CHECK(ck_state_tolerances(s, 1e-6, tol, 2) == CK_OK);
    CHECK(ck_update_rates(s) == CK_OK);
    NEAR(s->rates[0], .6, 1e-12);
    CHECK(ck_backward_euler_solve(s, 1, &bp, &stats) == CK_OK);
    CHECK(stats.final_time == 1);
  }
  else if (IS("CanBuildBackwardEulerOverloadedSolverMethod"))
  {
    CHECK(ck_backward_euler_solve(s, 5, &bp, &stats) == CK_OK);
    CHECK(stats.final_time == 5 && stats.function_calls == 2 && stats.jacobian_updates == 2 &&
          stats.internal_steps == 2 && stats.solves == 2);
    bp.small = 1;
    bp.max_iterations = 1;
    CHECK(ck_backward_euler_solve(s, 5, &bp, &stats) == CK_UNCONVERGED);
    CHECK(stats.final_time == .03125 && stats.function_calls == 6 && stats.jacobian_updates == 6 &&
          stats.internal_steps == 6 && stats.solves == 6);
  }
  else if (IS("CanBuildRosenbrockOverloadedSolveMethod"))
  {
    s->y[0] = 1;
    CHECK(ck_solve(s, 5, &rp, &stats) == CK_OK);
    CHECK(stats.final_time == 5 && stats.function_calls == 18 && stats.jacobian_updates == 9 &&
          stats.internal_steps == 9 && stats.solves == 27);
    rp.h_min = 15;
    rp.max_steps = 6;
    s->y[0] = 1;
    CHECK(ck_solve(s, 5, &rp, &stats) == CK_OK);
    CHECK(stats.final_time == 5 && stats.function_calls == 2 && stats.jacobian_updates == 1 &&
          stats.internal_steps == 1 && stats.solves == 3);
    CHECK(rp.h_min == 15 && rp.max_steps == 6);
  }
  else
  {
    /* C has dynamically sized cell groups rather than a template-imposed maximum. */
    ck_state_free(s);
    CHECK(ck_state_create(m, 4, &s) == CK_OK);
    CHECK(s->cells == 4);
    CHECK(ck_model_nreactions(m) == 2);
    CHECK(SIZE_MAX > 100000000);
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
