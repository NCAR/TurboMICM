/* External aerosol models ported from stub_aerosol_{1,2}.hpp. */
#include "test.h"
typedef struct
{
  const char *reactants[2], *products[2], *parameters[2];
  double rates[2];
  size_t r[2], p[2], parameter[2], diagonal[2], offdiagonal[2];
  int temperature_rate;
} transfer_model;
static void *transfer_clone(const void *source)
{
  transfer_model *copy = malloc(sizeof *copy);
  if (copy)
    *copy = *(const transfer_model *)source;
  return copy;
}
static ck_status transfer_finalize(const ck_model *m, void *context)
{
  transfer_model *t = context;
  for (size_t k = 0; k < 2; ++k)
  {
    if (ck_species_index(m, t->reactants[k], t->r + k) || ck_species_index(m, t->products[k], t->p + k))
      return CK_UNKNOWN;
    t->parameter[k] = SIZE_MAX;
    if (t->parameters[k] && ck_custom_index(m, t->parameters[k], t->parameter + k))
      return CK_UNKNOWN;
    t->diagonal[k] = ck_jacobian_index(m, t->r[k], t->r[k]);
    t->offdiagonal[k] = ck_jacobian_index(m, t->p[k], t->r[k]);
    if (t->diagonal[k] == SIZE_MAX || t->offdiagonal[k] == SIZE_MAX)
      return CK_INVALID;
  }
  return CK_OK;
}
static ck_status transfer_update(ck_state *s, void *context)
{
  transfer_model *t = context;
  if (t->temperature_rate)
    for (size_t c = 0; c < s->cells; ++c)
      s->custom[t->parameter[1] * s->cells + c] = .005 * s->conditions[c].temperature;
  return CK_OK;
}
static void transfer_forcing(const ck_state *s, const double *y, double *f, void *context)
{
  const transfer_model *t = context;
  for (size_t k = 0; k < 2; ++k)
    for (size_t c = 0; c < s->cells; ++c)
    {
      double rate = (t->parameter[k] == SIZE_MAX ? t->rates[k] : s->custom[t->parameter[k] * s->cells + c]) *
                    y[t->r[k] * s->cells + c];
      f[t->r[k] * s->cells + c] -= rate;
      f[t->p[k] * s->cells + c] += rate;
    }
}
static void transfer_jacobian(const ck_state *s, const double *y, double *a, void *context)
{
  (void)y;
  const transfer_model *t = context;
  for (size_t k = 0; k < 2; ++k)
    for (size_t c = 0; c < s->cells; ++c)
    {
      double rate = t->parameter[k] == SIZE_MAX ? t->rates[k] : s->custom[t->parameter[k] * s->cells + c];
      a[t->diagonal[k] * s->cells + c] += rate;
      a[t->offdiagonal[k] * s->cells + c] -= rate;
    }
}
static const ck_species gas[] = {{.name = "FO2"}, {.name = "BAR"}};
static const ck_species stub1[] = {{.name = "STUB1.MODE1.QUUX.BAZ"},  {.name = "STUB1.MODE1.QUUX.QUX"},
                                   {.name = "STUB1.MODE2.CORGE.BAZ"}, {.name = "STUB1.MODE2.CORGE.FO2"},
                                   {.name = "STUB1.MODE2.CORGE.QUX"}, {.name = "STUB1.MODE2.QUUX.BAZ"},
                                   {.name = "STUB1.MODE2.QUUX.QUX"}};
static const ck_species stub2[] = {
    {.name = "STUB2.MODE1.NUMBER"},    {.name = "STUB2.MODE1.QUUX.BAZ"},  {.name = "STUB2.MODE1.QUUX.QUX"},
    {.name = "STUB2.MODE2.CORGE.BAZ"}, {.name = "STUB2.MODE2.CORGE.FO2"}, {.name = "STUB2.MODE2.CORGE.QUX"},
    {.name = "STUB2.MODE2.NUMBER"},    {.name = "STUB2.MODE3.CORGE.BAZ"}, {.name = "STUB2.MODE3.CORGE.FO2"},
    {.name = "STUB2.MODE3.CORGE.QUX"}, {.name = "STUB2.MODE3.NUMBER"},    {.name = "STUB2.MODE3.QUUX.BAZ"},
    {.name = "STUB2.MODE3.QUUX.QUX"}};
static const char *const parameters[] = {"STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT",
                                         "STUB2.PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT"};
static const ck_named_jacobian_entry pattern1[] = {{"FO2", "FO2", 0},
                                                   {"STUB1.MODE2.CORGE.FO2", "FO2", 0},
                                                   {"STUB1.MODE1.QUUX.BAZ", "STUB1.MODE1.QUUX.BAZ", 0},
                                                   {"STUB1.MODE2.QUUX.BAZ", "STUB1.MODE1.QUUX.BAZ", 0}};
static const ck_named_jacobian_entry pattern2[] = {{"STUB2.MODE2.CORGE.FO2", "STUB2.MODE2.CORGE.FO2", 0},
                                                   {"STUB2.MODE2.CORGE.BAZ", "STUB2.MODE2.CORGE.FO2", 0},
                                                   {"STUB2.MODE3.QUUX.BAZ", "STUB2.MODE3.QUUX.BAZ", 0},
                                                   {"STUB2.MODE3.QUUX.QUX", "STUB2.MODE3.QUUX.BAZ", 0}};
static int aerosol_case(const char *test, int method, ck_lu_algorithm lu)
{
  transfer_model t1 = {.reactants = {"FO2", "STUB1.MODE1.QUUX.BAZ"},
                       .products = {"STUB1.MODE2.CORGE.FO2", "STUB1.MODE2.QUUX.BAZ"},
                       .rates = {.001, .002}};
  transfer_model t2 = {.reactants = {"STUB2.MODE2.CORGE.FO2", "STUB2.MODE3.QUUX.BAZ"},
                       .products = {"STUB2.MODE2.CORGE.BAZ", "STUB2.MODE3.QUUX.QUX"},
                       .parameters = {parameters[0], parameters[1]},
                       .temperature_rate = 1};
  ck_external_model external[2] = {{.name = "STUB1",
                                    .nspecies = 7,
                                    .species = stub1,
                                    .njacobian = 4,
                                    .jacobian = pattern1,
                                    .context = &t1},
                                   {.name = "STUB2",
                                    .nspecies = 13,
                                    .species = stub2,
                                    .ncustom = 2,
                                    .custom_names = parameters,
                                    .njacobian = 4,
                                    .jacobian = pattern2,
                                    .context = &t2}};
  for (size_t k = 0; k < 2; ++k)
  {
    external[k].clone_context = transfer_clone;
    external[k].destroy_context = free;
    external[k].finalize = transfer_finalize;
    external[k].update = transfer_update;
    external[k].forcing = transfer_forcing;
    external[k].jacobian_terms = transfer_jacobian;
  }
  ck_mechanism d = {.name = "aerosol", .nspecies = 2, .species = gas, .nexternal = 2, .external = external};
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create(&d, &m) == CK_OK);
  size_t cells = strstr(test, "Multi") ? 3 : 1;
  CHECK(ck_state_create(m, cells, &s) == CK_OK);
  CHECK(ck_model_nspecies(m) == 22 && ck_model_ncustom(m) == 2);
  for (size_t k = 0; k < 22; ++k)
  {
    const char *name = k < 2 ? gas[k].name : k < 9 ? stub1[k - 2].name : stub2[k - 9].name;
    size_t index;
    CHECK(ck_species_index(m, name, &index) == CK_OK);
    CHECK(!strcmp(ck_species_name(m, index), name));
  }
  const char *names[] = {"FO2",
                         "BAR",
                         "STUB1.MODE1.QUUX.BAZ",
                         "STUB1.MODE2.CORGE.FO2",
                         "STUB2.MODE3.CORGE.QUX",
                         "STUB2.MODE3.CORGE.BAZ",
                         "STUB2.MODE1.NUMBER",
                         "STUB2.MODE2.NUMBER",
                         "STUB1.MODE2.QUUX.BAZ",
                         "STUB2.MODE2.CORGE.FO2",
                         "STUB2.MODE2.CORGE.BAZ",
                         "STUB2.MODE3.QUUX.BAZ",
                         "STUB2.MODE3.QUUX.QUX"};
  double initial[13][3] = {{1, 1.5, 2},   {2, 2.5, 3},        {.5, .7, .9},    {.8, 1, 1.2}, {.3, .4, .5},
                           {.2, .3, .4},  {1000, 1100, 1200}, {500, 550, 600}, {0, 0, 0},    {.4, .6, .8},
                           {.1, .15, .2}, {.2, .3, .4},       {.3, .4, .5}};
  int two = strstr(test, "TwoStub") != NULL, update = strstr(test, "Update") != NULL;
  if (!two)
    for (size_t k = 9; k < 13; ++k)
      for (size_t c = 0; c < cells; ++c)
        initial[k][c] = 0;
  if (two)
  {
    initial[5][1] = .25;
    initial[5][2] = .3;
  }
  if (update)
  {
    double base[] = {1.34, 2.53, .75, 1.23, .42, .33, 1000, 500};
    for (size_t k = 0; k < 8; ++k)
      for (size_t c = 0; c < cells; ++c)
        initial[k][c] = base[k] + (k < 6 ? .01 : 1) * c;
  }
  size_t indices[13];
  for (size_t k = 0; k < 13; ++k)
  {
    CHECK(ck_species_index(m, names[k], indices + k) == CK_OK);
    CHECK(ck_state_set(s, names[k], initial[k], cells) == CK_OK);
  }
  if (two)
  {
    double k[] = {.01, .015, .02};
    CHECK(ck_state_set_custom(s, parameters[0], k, cells) == CK_OK);
    for (size_t c = 0; c < cells; ++c)
      s->conditions[c].temperature = 275 + 10 * c;
  }
  CHECK(ck_update_rates(s) == CK_OK);
  if (strstr(test, "Forcing") || strstr(test, "Jacobian"))
  {
    double f[22];
    size_t nnz = ck_jacobian_size(m);
    double *j = calloc(nnz, sizeof(double));
    CHECK(j);
    ck_forcing(s, s->y, f);
    ck_jacobian(s, s->y, j);
    NEAR(f[indices[0]], -.001, 1e-18);
    NEAR(f[indices[3]], .001, 1e-18);
    NEAR(f[indices[2]], -.001, 1e-18);
    NEAR(f[indices[8]], .001, 1e-18);
    NEAR(j[ck_jacobian_index(m, indices[0], indices[0])], .001, 1e-18);
    NEAR(j[ck_jacobian_index(m, indices[3], indices[0])], -.001, 1e-18);
    NEAR(j[ck_jacobian_index(m, indices[2], indices[2])], .002, 1e-18);
    NEAR(j[ck_jacobian_index(m, indices[8], indices[2])], -.002, 1e-18);
    free(j);
  }
  else if (strstr(test, "CanSolve"))
  {
    ck_stats stats;
    ck_status status;
    if (method)
    {
      ck_rosenbrock_parameters p;
      CHECK(ck_rosenbrock_parameters_init(method, &p) == CK_OK);
      p.lu = lu;
      status = ck_solve(s, 10, &p, &stats);
    }
    else
    {
      ck_backward_euler_parameters p = ck_backward_euler_defaults();
      p.lu = lu;
      status = ck_backward_euler_solve(s, 10, &p, &stats);
    }
    CHECK(status == CK_OK);
    CHECK(stats.final_time == 10);
    double tolerance = method ? (two && method == 2 ? 5e-4 : 5e-5) : (two ? .1 : 2e-4);
    for (size_t c = 0; c < cells; ++c)
    {
      double delta1 = initial[0][c] * (1 - exp(-.001 * 10)), delta2 = initial[2][c] * (1 - exp(-.002 * 10));
      double expected[13];
      for (size_t k = 0; k < 13; ++k)
        expected[k] = initial[k][c];
      expected[0] -= delta1;
      expected[3] += delta1;
      expected[2] -= delta2;
      expected[8] += delta2;
      if (two)
      {
        double delta3 = initial[9][c] * (1 - exp(-(.01 + .005 * c) * 10)),
               delta4 = initial[11][c] * (1 - exp(-(275.0 + 10 * c) * .005 * 10));
        expected[9] -= delta3;
        expected[10] += delta3;
        expected[11] -= delta4;
        expected[12] += delta4;
      }
      size_t ref[] = {0, 1, 2, 0, 4, 5, 6, 7, 2, 9, 9, 11, 11};
      for (size_t k = 0; k < 13; ++k)
      {
        double tol =
            (k == 0 || k == 2 || k == 3 || k == 8 || (two && k >= 9)) ? tolerance * initial[ref[k]][c] : 0;
        NEAR(s->y[indices[k] * cells + c], expected[k], tol);
      }
    }
  }
  else
    for (size_t k = 0; k < 13; ++k)
      for (size_t c = 0; c < cells; ++c)
        CHECK(s->y[indices[k] * cells + c] == initial[k][c]);
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
int test_aerosol(const char *id, const char *test)
{
  if (strncmp(id, "aerosol_model_", 14))
    return 2;
  int backward = strstr(id, "backward") != NULL, methods[] = {2, 3, 4, 43, 6};
  for (int lu = 0; lu < 2; ++lu)
    for (size_t k = 0; k < (backward ? 1u : 5u); ++k)
      CHECK(aerosol_case(test, backward ? 0 : methods[k], (ck_lu_algorithm)lu) == 0);
  return 0;
}
