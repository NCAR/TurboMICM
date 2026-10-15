#include "test.h"
#include "internal.h"
#include <float.h>

static double density_scale(const ck_conditions *c, const double *v, void *context)
{
  (void)v;
  return c->air_density * (*(double *)context);
}
static double lambda(const ck_conditions *c, const double *v, void *context)
{
  (void)v;
  (void)context;
  return c->temperature * .001;
}
static int process_update(void)
{
  double scale = .82;
  ck_species sp[] = {{.name = "foo", .molecular_weight = .025, .diffusion = 2.3e2},
                     {.name = "bar", .parameterized = 2, .function = density_scale, .context = &scale}};
  size_t r[] = {0, 1, 0, 1};
  ck_reaction rx[3] = {0};
  for (size_t k = 0; k < 3; ++k)
  {
    rx[k].nreactants = k ? 1 : 2;
    rx[k].reactants = r + (k == 0 ? 0 : k == 1 ? 2 : 3);
  }
  rx[0].rate = ck_rate_defaults(CK_ARRHENIUS);
  rx[0].rate.A = 12.2;
  rx[0].rate.C = 300;
  rx[1].rate = ck_rate_defaults(CK_SURFACE);
  rx[1].rate.label = "foo_surf";
  rx[1].rate.diffusion = 2.3e2;
  rx[1].rate.molecular_weight = .025;
  rx[2].rate = ck_rate_defaults(CK_USER);
  rx[2].rate.label = "bar_user";
  ck_mechanism d = {.name = "gas", .nspecies = 2, .species = sp, .nreactions = 3, .reactions = rx};
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create(&d, &m) == CK_OK);
  CHECK(ck_state_create(m, 5, &s) == CK_OK);
  for (size_t c = 0; c < 5; ++c)
  {
    s->conditions[c] = (ck_conditions){285 * (1 + .001 * c), 101100 * (1 + .003 * c), 10 * (1 + .002 * c)};
    s->custom[c] = 1e-2 * (1 + .004 * c);
    s->custom[5 + c] = 1e-8 * (1 + .005 * c);
    s->custom[10 + c] = 1e5 * (1 + .006 * c);
  }
  CHECK(ck_update_rates(s) == CK_OK);
  for (size_t c = 0; c < 5; ++c)
  {
    double T = s->conditions[c].temperature, rho = s->conditions[c].air_density;
    double expected = 12.2 * exp(300 / T) * rho * .82;
    NEAR(s->rates[c], expected, 1e-12 * expected);
    expected = s->custom[c] * rho * .82;
    NEAR(s->rates[5 + c], expected, 1e-12 * expected);
    double speed = sqrt(8 * CK_GAS_CONSTANT / (CK_PI * .025) * T), radius = s->custom[5 + c],
           number = s->custom[10 + c];
    expected = 4 * number * CK_PI * radius * radius / (radius / 2.3e2 + 4 / speed);
    NEAR(s->rates[10 + c], expected, 1e-10 * expected);
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
static int store_case(const char *test)
{
  ck_species sp[] = {
      {.name = "a"}, {.name = "b"}, {.name = "c", .molecular_weight = .025, .diffusion = 1e-5}};
  size_t reactants[] = {0, 0, 0, 0, 0, 2};
  ck_product product = {1, 1};
  ck_reaction rx[6] = {0};
  size_t nr = 1;
  ck_rate_type types[] = {CK_ARRHENIUS, CK_ARRHENIUS, CK_TROE, CK_TUNNELING, CK_USER, CK_SURFACE};
  for (size_t k = 0; k < 6; ++k)
  {
    rx[k] = (ck_reaction){.nreactants = 1,
                          .reactants = reactants + k,
                          .nproducts = 1,
                          .products = &product,
                          .rate = ck_rate_defaults(types[k])};
    rx[k].rate.label = k == 4 ? "ud1" : "surf1";
    rx[k].rate.diffusion = 1e-5;
    rx[k].rate.molecular_weight = .025;
  }
  if (IS("OffsetsAreContiguousCumulativeSizes"))
    nr = 6;
  if (IS("ArrheniusParametersPreserved"))
  {
    rx[0].rate.A = 2.15e-4;
    rx[0].rate.B = 1.2;
    rx[0].rate.C = 110;
  }
  if (IS("BranchedDerivedFieldsComputed"))
  {
    rx[0].rate = ck_rate_defaults(CK_BRANCH_ALKOXY);
    rx[0].rate.X = 1;
    rx[0].rate.a0 = .5;
    rx[0].rate.n = 3;
  }
  if (IS("UserDefinedCustomParamIndex"))
  {
    nr = 2;
    for (size_t k = 0; k < 2; ++k)
    {
      rx[k].rate = ck_rate_defaults(CK_USER);
      rx[k].rate.label = k ? "rate2" : "rate1";
      rx[k].rate.scaling = k ? .5 : 2;
    }
  }
  if (strstr(test, "Surface") && !IS("OffsetsAreContiguousCumulativeSizes"))
  {
    rx[0].rate = ck_rate_defaults(CK_SURFACE);
    rx[0].rate.label = "surf";
    rx[0].rate.diffusion = 2.3e2;
    rx[0].rate.molecular_weight = .025;
    rx[0].rate.probability = .74;
    if (IS("SurfaceMissingDiffusionCoefficientThrows"))
      rx[0].rate.diffusion = 0;
    if (IS("SurfaceMissingMolecularWeightThrows"))
      rx[0].rate.molecular_weight = 0;
    if (IS("SurfaceCustomParamIndexAfterUserDefined"))
    {
      nr = 2;
      rx[1].rate = rx[0].rate;
      rx[0].rate = ck_rate_defaults(CK_USER);
      rx[0].rate.label = "ud";
    }
  }
  if (IS("LambdaEntriesRcIndex"))
  {
    nr = 2;
    rx[1].rate = ck_rate_defaults(CK_LAMBDA);
    rx[1].rate.label = "lam";
    rx[1].rate.function = lambda;
  }
  double multiplier = 2;
  if (IS("ParameterizedMultipliers"))
  {
    nr = 2;
    sp[1].parameterized = 2;
    sp[1].function = density_scale;
    sp[1].context = &multiplier;
    reactants[1] = 1;
  }
  if (IS("TotalRateConstantsMatchesProcessCount"))
  {
    nr = 5;
    rx[1].rate = rx[2].rate;
    rx[2].rate = rx[4].rate;
    rx[3].rate = rx[5].rate;
    rx[4].rate = ck_rate_defaults(CK_LAMBDA);
    rx[4].rate.function = lambda;
    rx[4].rate.label = "lam";
  }
  ck_mechanism d = {.name = "store", .nspecies = 3, .species = sp, .nreactions = nr, .reactions = rx};
  ck_model *m = NULL;
  if (strstr(test, "Missing"))
  {
    CHECK(ck_model_create(&d, &m) == CK_INVALID);
    CHECK(!m);
    return 0;
  }
  CHECK(ck_model_create(&d, &m) == CK_OK);
  CHECK(m->nr == nr);
  if (IS("OffsetsAreContiguousCumulativeSizes"))
  {
    for (size_t i = 0; i < 6; ++i)
      CHECK(m->rates[i].type == types[i]);
    const size_t offsets[] = {0, 2, 3, 3, 3, 3, 4, 4, 4, 5, 6};
    const ck_rate_type order[] = {CK_ARRHENIUS,     CK_TROE,      CK_TERNARY, CK_TERNARY_JPL19,
                                  CK_BRANCH_ALKOXY, CK_TUNNELING, CK_TAYLOR,  CK_REVERSIBLE,
                                  CK_USER,          CK_SURFACE,   CK_LAMBDA};
    size_t cursor = 0;
    for (size_t t = 0; t < 11; ++t)
    {
      CHECK(cursor == offsets[t]);
      while (cursor < nr && m->rates[cursor].type == order[t])
        ++cursor;
    }
    CHECK(cursor == nr);
  }
  if (IS("ArrheniusParametersPreserved"))
  {
    CHECK(m->rates[0].A == 2.15e-4 && m->rates[0].B == 1.2 && m->rates[0].C == 110 && m->rates[0].D == 300);
  }
  if (IS("BranchedDerivedFieldsComputed"))
  {
    double k = 2e-22 * CK_AVOGADRO * 1e-6 * exp(3), a = k * 2.45e19 / CK_AVOGADRO * 1e6,
           b = .43 * pow(293. / 298, -8);
    double z = a / (1 + a / b) * pow(.41, 1 / (1 + pow(log10(a / b), 2)));
    NEAR(m->rates[0].k0, k, 1e-10 * k);
    NEAR(m->rates[0].z, z, 1e-10 * fabs(z));
  }
  if (IS("UserDefinedCustomParamIndex"))
  {
    CHECK(m->custom_offset[0] == 0 && m->custom_offset[1] == 1);
    CHECK(m->rates[0].scaling == 2 && m->rates[1].scaling == .5);
  }
  if (IS("SurfaceDataFieldsAndCustomParamIndex"))
  {
    CHECK(m->rates[0].diffusion == 2.3e2 && m->rates[0].molecular_weight == .025 &&
          m->rates[0].probability == .74 && m->custom_offset[0] == 0);
    NEAR(8 * CK_GAS_CONSTANT / (CK_PI * m->rates[0].molecular_weight), 8 * CK_GAS_CONSTANT / (CK_PI * .025),
         1e-14);
  }
  if (IS("SurfaceCustomParamIndexAfterUserDefined"))
    CHECK(m->custom_offset[0] == 0 && m->custom_offset[1] == 1 && m->nc == 3);
  if (IS("LambdaEntriesRcIndex"))
  {
    CHECK(m->rates[1].type == CK_LAMBDA && m->rates[1].function);
    ck_conditions c = {.temperature = 300};
    NEAR(ck_rate_calculate(m->rates + 1, &c, NULL), .3, 1e-14);
  }
  if (IS("ParameterizedMultipliers"))
  {
    ck_state *s;
    CHECK(ck_state_create(m, 1, &s) == CK_OK);
    s->conditions[0] = (ck_conditions){300, 1, 5};
    CHECK(ck_update_rates(s) == CK_OK);
    CHECK(m->nreactants[0] == 1 && m->nreactants[1] == 0);
    NEAR(s->rates[1], 10, 1e-14);
    ck_state_free(s);
  }
  ck_model_free(m);
  return 0;
}
int test_process(const char *suite, const char *test)
{
  if (!strcmp(suite, "ReactionRateConstantStore"))
    return store_case(test);
  if (strcmp(suite, "Process"))
    return 2;
  if (IS("Matrix") || IS("VectorMatrix"))
    return process_update();
  if (IS("ChemicalReactionCopyAssignmentSucceeds"))
  {
    double scale = .82;
    ck_species species[] = {
        {.name = "foo", .molecular_weight = .025, .diffusion = 2.3e2},
        {.name = "bar", .parameterized = 2, .function = density_scale, .context = &scale}};
    size_t reactants[] = {0, 1};
    ck_reaction original = {.nreactants = 2, .reactants = reactants, .rate = ck_rate_defaults(CK_ARRHENIUS)};
    original.rate.A = 12.2;
    original.rate.C = 300;
    ck_reaction copy = original;
    CHECK(!strcmp(species[copy.reactants[0]].name, species[original.reactants[0]].name));
    CHECK(copy.nproducts == original.nproducts && copy.nproducts == 0);
    CHECK(copy.rate.type == CK_ARRHENIUS && copy.rate.A == 12.2 && copy.rate.C == 300);
    return 0;
  }
  ck_species sp[] = {{.name = "O3"}, {.name = "NO"}, {.name = "NO2"}, {.name = "O2"}};
  size_t r[] = {0, 1};
  ck_product p[] = {{2, 1}, {3, 1}};
  ck_reaction rx = {
      .nreactants = 2, .reactants = r, .nproducts = 2, .products = p, .rate = ck_rate_defaults(CK_ARRHENIUS)};
  ck_mechanism d = {.name = "gas", .nspecies = 4, .species = sp, .nreactions = 1, .reactions = &rx};
  ck_model *m;
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  CHECK(m->nreactants[0] == 2);
  CHECK(!strcmp(ck_species_name(m, m->products[0][0].species), "NO2"));
  CHECK(!strcmp(m->description->name, "gas"));
  if (IS("ChemicalReactionCopyAssignmentSucceeds"))
  {
    ck_reaction copy = rx;
    CHECK(copy.reactants[0] == rx.reactants[0] && copy.nproducts == rx.nproducts && copy.rate.A == rx.rate.A);
  }
  ck_model_free(m);
  return 0;
}
