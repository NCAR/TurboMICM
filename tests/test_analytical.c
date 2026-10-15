/* Analytical integration cases translated from analytical_policy.hpp. */
#include "test.h"

static double branch_reference(double t, double rho, int nitrate)
{
  double a = 2e-22 * exp(2) * 2.45e19, b = .43 * pow(293. / 298., -8);
  double z = a / (1 + a / b) * pow(.41, 1 / (1 + pow(log10(a / b), 2))) * .999 / .001;
  a = 2e-22 * exp(2) * rho * CK_AVOGADRO * 1e-6;
  b = .43 * pow(t / 298., -8);
  double A = a / (1 + a / b) * pow(.41, 1 / (1 + pow(log10(a / b), 2)));
  return 1e-4 * exp(-204.3 / t) * (nitrate ? A : z) / (z + A);
}
static double reference_rate(int kind, int second, int stiff, const ck_conditions *c)
{
  double t = c->temperature, rho = c->air_density;
  if (kind == 0 || kind == 2)
  {
    double k0 = second ? 1.2e-3 * exp(3 / t) * pow(t / 300, 1.6) : (kind == 0 ? 4e-11 : 4e-5);
    double ki = second ? 136 * exp(24 / t) * pow(t / 300, 5) : 1;
    double fc = second ? .9 : .6, N = second ? .8 : 1;
    return k0 * (kind == 0 ? rho : 1) / (1 + k0 * rho / ki) * pow(fc, N / (N + pow(log10(k0 * rho / ki), 2)));
  }
  if (kind == 1)
    return second ? 3e-3 : 2e-3;
  if (kind == 3)
    return second ? 1.2e-4 * exp(-167 / t + 1e8 / pow(t, 3)) : 4e-3;
  if (kind == 4)
    return second ? 1.2e-4 * exp(75 / t) * pow(t / 50, stiff ? 1.6 : 7) * (1 + .5 * c->pressure)
                  : 4e-3 * exp(50 / t);
  return branch_reference(t, rho, second);
}
static ck_rate analytical_rate(int kind, int second, int stiff, const char *label)
{
  ck_rate_type types[] = {CK_TROE, CK_USER, CK_TERNARY, CK_TUNNELING, CK_ARRHENIUS, CK_BRANCH_ALKOXY};
  ck_rate p = ck_rate_defaults(types[kind]);
  if (kind == 0 || kind == 2)
  {
    p.k0_A = second ? 1.2e-3 : (kind == 0 ? 4e-11 : 4e-5);
    if (second)
    {
      p.k0_B = 1.6;
      p.k0_C = 3;
      p.kinf_A = 136;
      p.kinf_B = 5;
      p.kinf_C = 24;
      p.Fc = .9;
      p.N = .8;
    }
  }
  else if (kind == 1)
    p.label = label;
  else if (kind == 3)
  {
    p.A = second ? 1.2e-4 : 4e-3;
    if (second)
    {
      p.B = 167;
      p.C = 1e8;
    }
  }
  else if (kind == 4)
  {
    p.A = second ? 1.2e-4 : 4e-3;
    p.C = second ? 75 : 50;
    if (second)
    {
      p.B = stiff ? 1.6 : 7;
      p.D = 50;
      p.E = .5;
    }
  }
  else
  {
    p.type = second ? CK_BRANCH_NITRATE : CK_BRANCH_ALKOXY;
    p.X = 1e-4;
    p.Y = 204.3;
    p.a0 = 1e-3;
    p.n = 2;
  }
  return p;
}
static double analytical_tolerance(int kind, int stiff, int method)
{
  if (method == 0)
  {
    double t[] = {1e-6, 1e-3, 1e-5, 1e-3, 1e-3, 1e-5};
    if (stiff)
      return kind == 0 ? 1e-5 : (kind == 2 || kind == 5 ? 1e-2 : 1e-3);
    return t[kind];
  }
  double regular[] = {1e-10, 2e-6, 1e-8, 1e-8, 1e-9, 1e-13};
  double hard[] = {1e-5, 2e-5, 2e-3, 1e-4, 2e-5, 2e-3};
  if (stiff && kind == 4)
    return method == 2 ? 1e-4 : method == 6 ? 1e-5 : 2e-5;
  if (!stiff && method == 2)
    return kind == 3 ? 2e-5 : kind == 4 ? 4e-6 : kind == 5 ? 1e-10 : regular[kind];
  return stiff ? hard[kind] : regular[kind];
}
static int simple_system(int kind, int stiff, int method, ck_lu_algorithm lu)
{
  size_t n = stiff ? 4 : 3, nr = stiff ? 5 : 2, cells = 3;
  ck_species sp[] = {
      {.name = stiff ? "A1" : "A"}, {.name = stiff ? "A2" : "B"}, {.name = stiff ? "B" : "C"}, {.name = "C"}};
  size_t reactants[] = {0, stiff ? 1 : 1, 2, 0, 1};
  ck_product products[] = {{stiff ? 2 : 1, 1}, {stiff ? 2 : 2, 1}, {3, 1}, {1, 1}, {0, 1}};
  ck_reaction rx[5] = {0};
  for (size_t r = 0; r < nr; ++r)
  {
    rx[r].nreactants = rx[r].nproducts = 1;
    rx[r].reactants = reactants + r;
    rx[r].products = products + r;
    rx[r].rate = r >= 3 ? ck_rate_defaults(CK_ARRHENIUS)
                        : analytical_rate(kind, stiff ? r == 2 : r == 1, stiff,
                                          r == 0            ? "photoA1B"
                                          : r == 1 && stiff ? "photoA2B"
                                                            : "photoB");
    if (r >= 3)
      rx[r].rate.A = r == 3 ? 4e10 : 3.6e10;
  }
  ck_mechanism d = {.name = "analytical", .nspecies = n, .species = sp, .nreactions = nr, .reactions = rx};
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create(&d, &m) == CK_OK);
  CHECK(ck_state_create(m, cells, &s) == CK_OK);
  size_t a1, a2 = 0, b, c;
  CHECK(ck_species_index(m, stiff ? "A1" : "A", &a1) == CK_OK);
  if (stiff)
    CHECK(ck_species_index(m, "A2", &a2) == CK_OK);
  CHECK(ck_species_index(m, "B", &b) == CK_OK);
  CHECK(ck_species_index(m, "C", &c) == CK_OK);
  double temps[] = {272.5, 254.7, 312.6}, pressures[] = {101253.3, 100672.5, 101319.8}, k1[3], k2[3];
  for (size_t cell = 0; cell < cells; ++cell)
  {
    s->conditions[cell] = (ck_conditions){temps[cell], pressures[cell], 0};
    ck_ideal_air_density(s->conditions + cell);
    k1[cell] = reference_rate(kind, 0, stiff, s->conditions + cell);
    k2[cell] = reference_rate(kind, 1, stiff, s->conditions + cell);
    s->y[a1 * cells + cell] = (stiff ? .5 : 1) * (1 - (double)cell / cells);
    if (stiff)
      s->y[a2 * cells + cell] = .5 * (1 - (double)cell / cells);
  }
  if (kind == 1)
  {
    double rates[3] = {.002, .002, .002};
    CHECK(ck_state_set_custom(s, "photoA1B", rates, 3) == CK_OK);
    if (stiff)
      CHECK(ck_state_set_custom(s, "photoA2B", rates, 3) == CK_OK);
    for (size_t cell = 0; cell < cells; ++cell)
      rates[cell] = .003;
    CHECK(ck_state_set_custom(s, "photoB", rates, 3) == CK_OK);
  }
  ck_rosenbrock_parameters p;
  ck_backward_euler_parameters bp = ck_backward_euler_defaults();
  bp.lu = lu;
  if (method)
  {
    CHECK(ck_rosenbrock_parameters_init(method, &p) == CK_OK);
    p.lu = lu;
  }
  for (size_t time = 1; time < 1000; ++time)
  {
    CHECK(ck_update_rates(s) == CK_OK);
    ck_stats stats;
    ck_status status = method ? ck_solve(s, 1, &p, &stats) : ck_backward_euler_solve(s, 1, &bp, &stats);
    CHECK(status == CK_OK);
    for (size_t cell = 0; cell < cells; ++cell)
    {
      double initial = 1 - (double)cell / cells, t = (double)time;
      double e1 = exp(-k1[cell] * t), e2 = exp(-k2[cell] * t);
      double expected[] = {initial * e1, initial * k1[cell] / (k2[cell] - k1[cell]) * (e1 - e2),
                           initial * (1 + (k1[cell] * e2 - k2[cell] * e1) / (k2[cell] - k1[cell]))};
      double actual[] = {s->y[a1 * cells + cell] + (stiff ? s->y[a2 * cells + cell] : 0),
                         s->y[b * cells + cell], s->y[c * cells + cell]};
      for (size_t j = 0; j < 3; ++j)
        NEAR(actual[j], expected[j], analytical_tolerance(kind, stiff, method));
      if (!stiff)
      {
        NEAR(s->rates[cell], k1[cell], 1e-9 * fabs(k1[cell]));
        NEAR(s->rates[cells + cell], k2[cell], 1e-9 * fabs(k2[cell]));
      }
    }
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
static int surface_system(int method, ck_lu_algorithm lu)
{
  ck_species sp[] = {{.name = "foo"}, {.name = "bar"}, {.name = "baz"}};
  size_t r = 0;
  ck_product products[] = {{1, 1}, {2, .4}};
  ck_reaction rx = {.nreactants = 1, .reactants = &r, .nproducts = 2, .products = products};
  rx.rate = ck_rate_defaults(CK_SURFACE);
  rx.rate.label = "foo";
  rx.rate.molecular_weight = .04607;
  rx.rate.diffusion = .95e-5;
  rx.rate.probability = .02;
  ck_mechanism d = {.name = "surface",
                    .nspecies = 3,
                    .species = sp,
                    .nreactions = 1,
                    .reactions = &rx,
                    .temperature = 272.5,
                    .pressure = 101253.3};
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  CHECK(ck_state_create(m, 1, &s) == CK_OK);
  double radius = 1e-6 / 2 * exp(5 * log(.1) * log(.1) / 2),
         number = 6 / (CK_PI * pow(1e-6, 3) * exp(4.5 * log(.1) * log(.1))) * (.002 / 1000 + .003 / 1000);
  CHECK(ck_state_set_custom(s, "foo.effective radius [m]", &radius, 1) == CK_OK);
  CHECK(ck_state_set_custom(s, "foo.particle number concentration [# m-3]", &number, 1) == CK_OK);
  s->y[0] = 1;
  double speed = sqrt(8 * CK_GAS_CONSTANT / (CK_PI * .04607) * 272.5);
  double k = 4 * number * CK_PI * radius * radius / (radius / .95e-5 + 4 / (speed * .02)), dt = .1 / k;
  double tolerance = method == 0 ? .05 : method == 2 ? .01 : method == 4 ? 1e-6 : method == 6 ? 1e-7 : 1e-5;
  ck_rosenbrock_parameters p;
  ck_backward_euler_parameters bp = ck_backward_euler_defaults();
  bp.lu = lu;
  if (method)
  {
    CHECK(ck_rosenbrock_parameters_init(method, &p) == CK_OK);
    p.lu = lu;
  }
  for (size_t t = 1; t <= 10; ++t)
  {
    CHECK(ck_update_rates(s) == CK_OK);
    ck_stats stats;
    CHECK((method ? ck_solve(s, dt, &p, &stats) : ck_backward_euler_solve(s, dt, &bp, &stats)) == CK_OK);
    NEAR(s->rates[0], k, 1e-8);
    double expected[] = {exp(-k * t * dt), 1 - exp(-k * t * dt), .4 * (1 - exp(-k * t * dt))};
    for (size_t i = 0; i < 3; ++i)
      NEAR(s->y[i], expected[i], tolerance * fabs(expected[i]));
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
int test_analytical(const char *id, const char *test)
{
  if (strncmp(id, "analytical_", 11))
    return 2;
  int backward = strstr(id, "backward") != NULL, stiff = strstr(test, "SuperStiff") != NULL;
  const char *names[] = {"Troe",      "Photolysis", "TernaryChemicalActivation",
                         "Tunneling", "Arrhenius",  "Branched"};
  int kind = -1;
  for (int j = 0; j < 6; ++j)
    if (!strncmp(test, names[j], strlen(names[j])))
      kind = j;
  if (kind < 0 && !IS("SurfaceRxn"))
    return 2;
  int methods[] = {2, 3, 4, 43, 6};
  for (int lu = 0; lu < 2; ++lu)
    for (size_t k = 0; k < (backward ? 1u : 5u); ++k)
    {
      int method = backward ? 0 : methods[k];
      CHECK((kind < 0 ? surface_system(method, (ck_lu_algorithm)lu)
                      : simple_system(kind, stiff, method, (ck_lu_algorithm)lu)) == 0);
    }
  return 0;
}
