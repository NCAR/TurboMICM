#include "test.h"
#include "analytical_expected.h"
#include "tutorials.h"

static int integrate_interval(ck_state *s, double dt, int method, ck_lu_algorithm lu, size_t substeps)
{
  ck_rosenbrock_parameters rp;
  ck_backward_euler_parameters bp = ck_backward_euler_defaults();
  bp.lu = lu;
  if (method)
  {
    CHECK(ck_rosenbrock_parameters_init(method, &rp) == CK_OK);
    rp.lu = lu;
  }
  double advanced = 0, max_step = dt / substeps;
  while (advanced < dt)
  {
    ck_stats stats;
    double step = fmin(max_step, dt - advanced);
    ck_status status =
        method ? ck_solve(s, step, &rp, &stats) : ck_backward_euler_solve(s, step, &bp, &stats);
    CHECK(status == CK_OK || status == CK_MAX_STEPS);
    CHECK(advanced + stats.final_time > advanced);
    advanced += stats.final_time;
  }
  return 0;
}
static int stiff_case(const char *test, int method, ck_lu_algorithm lu)
{
  int robertson = IS("Robertson"), e5 = IS("E5"), hires = IS("HIRES"), oregon = IS("Oregonator");
  size_t n = robertson ? 3 : e5 ? 6 : hires ? 8 : 5, nr = robertson ? 3 : e5 ? 4 : hires ? 9 : 5;
  ck_species sp[8];
  char names[8][8];
  for (size_t j = 0; j < n; ++j)
  {
    snprintf(names[j], sizeof names[j], "Y%zu", j);
    sp[j] = (ck_species){.name = names[j]};
  }
  size_t re[9][2] = {{0, 0}, {1, 2}, {0, 2}, {3, 0}};
  size_t nre[9] = {1, 2, 2, 1};
  ck_product pr[9][3] = {{{1, 1}, {2, 1}}, {{4, 1}}, {{3, 1}}, {{2, 1}, {5, 1}}};
  size_t npr[9] = {2, 1, 1, 2};
  double rates[9] = {7.89e-10, 1.13e9, 1.1e7, 1.13e3};
  double initial[8] = {e5 ? 1.76e-3 : 1, 0, 0, 0, 0, 0, 0, hires ? .0057 : 0}, atol[8],
         rtol = robertson ? 1e-10
                : e5      ? 1e-13
                : hires   ? 1e-6
                          : 1e-9;
  double scales[3] = {1, 1, 1};
  if (robertson)
  {
    size_t r[][2] = {{0, 0}, {1, 1}, {1, 2}};
    memcpy(re, r, sizeof r);
    size_t counts[] = {1, 2, 2};
    memcpy(nre, counts, sizeof counts);
    ck_product p[][3] = {{{1, 1}}, {{1, 1}, {2, 1}}, {{0, 1}, {2, 1}}};
    memcpy(pr, p, sizeof p);
    size_t pc[] = {1, 2, 2};
    memcpy(npr, pc, sizeof pc);
    rates[0] = .04;
    rates[1] = 3e7;
    rates[2] = 1e4;
  }
  if (hires)
  {
    size_t r[][2] = {{0, 0}, {1, 0}, {2, 0}, {0, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}, {5, 7}};
    memcpy(re, r, sizeof r);
    size_t rc[] = {1, 1, 1, 0, 1, 1, 1, 1, 2};
    memcpy(nre, rc, sizeof rc);
    ck_product p[][3] = {
        {{1, 1}}, {{0, .43 / 8.75}, {3, 8.32 / 8.75}},        {{0, 8.32 / 10.03}, {3, 1.71 / 10.03}},
        {{0, 1}}, {{2, .43 / 1.12}, {5, .69 / 1.12}},         {{2, .035 / 1.745}, {5, 1.71 / 1.745}},
        {{4, 1}}, {{4, .43 / 1.81}, {5, .69 / 1.81}, {7, 1}}, {{6, 1}}};
    memcpy(pr, p, sizeof p);
    size_t pc[] = {1, 2, 2, 1, 2, 2, 1, 3, 1};
    memcpy(npr, pc, sizeof pc);
    double k[] = {1.71, 8.75, 10.03, .0007, 1.12, 1.745, .43, 1.81, 280};
    memcpy(rates, k, sizeof k);
  }
  if (oregon)
  {
    size_t r[][2] = {{1, 0}, {0, 1}, {0, 0}, {0, 0}, {2, 0}};
    memcpy(re, r, sizeof r);
    size_t rc[] = {1, 2, 1, 2, 1};
    memcpy(nre, rc, sizeof rc);
    ck_product p[][3] = {{{0, 1}}, {{3, 1}}, {{2, 1}, {0, 2}}, {{4, 1}}, {{1, 1}}};
    memcpy(pr, p, sizeof p);
    size_t pc[] = {1, 1, 2, 1, 1};
    memcpy(npr, pc, sizeof pc);
    rates[0] = 1 / (77.27 * .161);
    rates[1] = 1.6e9;
    rates[2] = 77.27 / .161;
    rates[3] = 77.27 * 77.27 * 8.375e-6 * 1.6e9 / 2;
    rates[4] = 1;
    scales[0] = 1 / (77.27 * .161 * 1.6e9);
    scales[1] = 77.27 / (.161 * 1.6e9);
    scales[2] = 1 / (.161 * .161 * 1.6e9);
    initial[0] = scales[0];
    initial[1] = 2 * scales[1];
    initial[2] = 3 * scales[2];
  }
  ck_reaction rx[9];
  char labels[9][8];
  for (size_t r = 0; r < nr; ++r)
  {
    snprintf(labels[r], sizeof labels[r], "r%zu", r + 1);
    rx[r] = (ck_reaction){.nreactants = nre[r],
                          .reactants = re[r],
                          .nproducts = npr[r],
                          .products = pr[r],
                          .rate = ck_rate_defaults(CK_USER)};
    rx[r].rate.label = labels[r];
  }
  ck_mechanism d = {
      .name = "stiff reference", .nspecies = n, .species = sp, .nreactions = nr, .reactions = rx};
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  CHECK(ck_state_create(m, 1, &s) == CK_OK);
  memcpy(s->y, initial, n * sizeof(double));
  for (size_t j = 0; j < n; ++j)
    atol[j] = rtol * 1e-2;
  if (e5)
  {
    for (size_t j = 0; j < n; ++j)
      atol[j] = j == 0 || j >= 4 ? 1e-7 : 1e-17;
  }
  if (oregon)
  {
    for (size_t j = 0; j < 3; ++j)
      atol[j] = scales[j] * rtol * 1e-2;
    atol[3] = atol[4] = scales[1];
  }
  CHECK(ck_state_tolerances(s, rtol, atol, n) == CK_OK);
  for (size_t r = 0; r < nr; ++r)
    CHECK(ck_state_set_custom(s, labels[r], rates + r, 1) == CK_OK);
  CHECK(ck_update_rates(s) == CK_OK);
  size_t steps = robertson || oregon ? 12 : e5 ? 7 : 2;
  double current = 0, target = e5 ? 10 : robertson ? 1 : hires ? 321.8122 : 30 * .161;
  for (size_t t = 0; t < steps; ++t)
  {
    double dt = robertson || e5 ? target - current : target;
    CHECK(integrate_interval(s, dt, method, lu, oregon && !method ? 18000 : 1) == 0);
    for (size_t j = 0; j < (e5 ? 4 : oregon ? 3 : n); ++j)
    {
      double expected = robertson ? expected_robertson[t + 1][j]
                        : e5      ? expected_e5[t + 1][j]
                        : hires   ? expected_hires[t + 1][j]
                                  : expected_oregonator[t + 1][j] * scales[j];
      double absolute = robertson ? .3e-3
                        : e5      ? (j ? 1e-13 : 1e-6)
                        : hires   ? (method == 0                                  ? .1
                                     : method == 2 || method == 43 || method == 6 ? 1e-6
                                                                                  : 1e-7)
                                  : 0;
      double relative = oregon ? (method ? 1e-6 : .005) : robertson ? 1e-6 : e5 ? 1e-10 : 0;
      if (hires)
        NEAR(s->y[j], expected, absolute);
      else
        CHECK(fabs(s->y[j] - expected) < absolute || fabs(s->y[j] - expected) < relative * fabs(expected));
    }
    if (robertson || e5)
    {
      current = target;
      target *= e5 ? 100 : 10;
    }
    if (hires)
      target += 100;
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
static int flow_tube(int backward)
{
  ck_species sp[] = {{.name = "SOA1"}, {.name = "SOA2"}, {.name = "a-pinene"}, {.name = "O3"}};
  size_t re[][2] = {{0, 0}, {1, 0}, {2, 3}};
  ck_product products[] = {{0, .18}, {1, .09}};
  ck_reaction rx[3];
  for (size_t r = 0; r < 3; ++r)
  {
    rx[r] = (ck_reaction){.nreactants = r == 2 ? 2 : 1,
                          .reactants = re[r],
                          .nproducts = r == 2 ? 2 : 0,
                          .products = products,
                          .rate = ck_rate_defaults(r == 2 ? CK_ARRHENIUS : CK_USER)};
    rx[r].rate.label = r ? "r2" : "r1";
  }
  rx[2].rate.A = 8.8e-17 * 1e-6 * CK_AVOGADRO;
  ck_mechanism d = {.name = "flow tube",
                    .nspecies = 4,
                    .species = sp,
                    .nreactions = 3,
                    .reactions = rx,
                    .temperature = 298.15,
                    .pressure = 101325};
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  CHECK(ck_state_create(m, 1, &s) == CK_OK);
  s->y[2] = 8e-8;
  s->y[3] = .00002;
  double k1 = .01, k2 = .05;
  CHECK(ck_state_set_custom(s, "r1", &k1, 1) == CK_OK);
  CHECK(ck_state_set_custom(s, "r2", &k2, 1) == CK_OK);
  CHECK(ck_update_rates(s) == CK_OK);
  FILE *f = fopen(backward ? "tests/expected/backward_euler_wall_loss.csv"
                           : "tests/expected/rosenbrock_wall_loss.csv",
                  "r");
  CHECK(f);
  char line[512];
  CHECK(fgets(line, sizeof line, f));
  CHECK(!strcmp(line, "time,SOA1,SOA2,a-pinene,O3\n"));
  for (size_t t = 0; t <= 3600; ++t)
  {
    double time, expected[4];
    CHECK(fscanf(f, "%lf,%lf,%lf,%lf,%lf", &time, expected, expected + 1, expected + 2, expected + 3) == 5);
    CHECK(time == (double)t);
    if (t < 10)
      for (size_t j = 0; j < 4; ++j)
        NEAR(s->y[j], expected[j], 1e-10);
    if (t < 3600)
      CHECK(integrate_interval(s, 1, backward ? 0 : 2, CK_LU_MOZART, 1) == 0);
  }
  CHECK(fclose(f) == 0);
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
int test_stiff(const char *id, const char *suite, const char *test)
{
  if (!strcmp(suite, "RosenbrockRegression") && IS("FlowTube"))
    return flow_tube(0);
  if (!strcmp(suite, "BackwardEulerRegression") && IS("FlowTube"))
    return flow_tube(1);
  if (strcmp(suite, "AnalyticalExamples") ||
      !(IS("Robertson") || IS("E5") || IS("Oregonator") || IS("HIRES")))
    return 2;
  int backward = strstr(id, "backward") != NULL, methods[] = {2, 3, 4, 43, 6};
  for (int lu = 0; lu < 2; ++lu)
    for (size_t k = 0; k < (backward ? 1u : 5u); ++k)
      CHECK(stiff_case(test, backward ? 0 : methods[k], (ck_lu_algorithm)lu) == 0);
  return 0;
}
