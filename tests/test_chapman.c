#include "test.h"
#include "internal.h"
static const ck_species chapman_species[] = {{.name = "M"},   {.name = "Ar"}, {.name = "CO2"},
                                             {.name = "H2O"}, {.name = "N2"}, {.name = "O1D"},
                                             {.name = "O"},   {.name = "O2"}, {.name = "O3"}};
static const size_t reactants[][3] = {{5, 4, 0}, {5, 7, 0}, {6, 8, 0}, {6, 7, 0},
                                      {7, 0, 0}, {8, 0, 0}, {8, 0, 0}};
static const ck_product products[][2] = {{{6, 1}, {4, 1}}, {{6, 1}, {7, 1}}, {{7, 2}, {0, 0}},
                                         {{8, 1}, {0, 1}}, {{6, 2}, {0, 0}}, {{5, 1}, {7, 1}},
                                         {{6, 1}, {7, 1}}};
static const char *const photo[] = {"jO2", "jO3a", "jO3b"};
static const size_t positions[][2] = {{0, 0}, {6, 0}, {7, 0}, {8, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4},
                                      {5, 4}, {6, 4}, {5, 5}, {6, 5}, {6, 6}, {7, 6}, {8, 6}, {5, 7},
                                      {6, 7}, {7, 7}, {8, 7}, {5, 8}, {6, 8}, {7, 8}, {8, 8}};
static void chapman_description(ck_mechanism *d, ck_reaction *rx, int integration)
{
  double A[] = {2.15e-11, 3.3e-11, 8e-12, 6e-34}, C[] = {110, 55, -2060, 0};
  for (size_t r = 0; r < 7; ++r)
  {
    rx[r] = (ck_reaction){.nreactants = r == 3  ? 3
                                        : r < 3 ? 2
                                                : 1,
                          .reactants = reactants[r],
                          .nproducts = r == 2 || r == 4 ? 1 : 2,
                          .products = products[r],
                          .rate = ck_rate_defaults(r < 4 ? CK_ARRHENIUS : CK_USER)};
    if (r < 4)
    {
      rx[r].rate.A = A[r];
      rx[r].rate.C = C[r];
      if (r == 3)
      {
        rx[r].rate.B = integration ? 0 : 2.4;
        rx[r].rate.C = integration ? 2.4 : 0;
      }
    }
    else
      rx[r].rate.label = photo[r - 4];
  }
  *d = (ck_mechanism){
      .name = "Chapman", .nspecies = 9, .species = chapman_species, .nreactions = 7, .reactions = rx};
}
/* Independent equations for the four reactive species; spectators have zero net stoichiometry. */
static void chapman_expected(const double *k, const double *y, double *f, double *jac)
{
  memset(f, 0, 9 * sizeof(double));
  memset(jac, 0, 81 * sizeof(double));
  f[5] = k[5] * y[8] - k[0] * y[4] * y[5] - k[1] * y[5] * y[7];
  f[6] = 2 * k[4] * y[7] + k[6] * y[8] + k[0] * y[4] * y[5] + k[1] * y[5] * y[7] - k[2] * y[6] * y[8] -
         k[3] * y[0] * y[6] * y[7];
  f[7] = -k[4] * y[7] + k[5] * y[8] + k[6] * y[8] + 2 * k[2] * y[6] * y[8] - k[3] * y[0] * y[6] * y[7];
  f[8] = -k[5] * y[8] - k[6] * y[8] - k[2] * y[6] * y[8] + k[3] * y[0] * y[6] * y[7];
  jac[6 * 9 + 0] = -k[3] * y[6] * y[7];
  jac[7 * 9 + 0] = jac[6 * 9 + 0];
  jac[8 * 9 + 0] = -jac[6 * 9 + 0];
  jac[5 * 9 + 4] = -k[0] * y[5];
  jac[6 * 9 + 4] = -jac[5 * 9 + 4];
  jac[5 * 9 + 5] = -k[0] * y[4] - k[1] * y[7];
  jac[6 * 9 + 5] = -jac[5 * 9 + 5];
  jac[6 * 9 + 6] = -k[2] * y[8] - k[3] * y[0] * y[7];
  jac[7 * 9 + 6] = 2 * k[2] * y[8] - k[3] * y[0] * y[7];
  jac[8 * 9 + 6] = -k[2] * y[8] + k[3] * y[0] * y[7];
  jac[5 * 9 + 7] = -k[1] * y[5];
  jac[6 * 9 + 7] = 2 * k[4] + k[1] * y[5] - k[3] * y[0] * y[6];
  jac[7 * 9 + 7] = -k[4] - k[3] * y[0] * y[6];
  jac[8 * 9 + 7] = k[3] * y[0] * y[6];
  jac[5 * 9 + 8] = k[5];
  jac[6 * 9 + 8] = k[6] - k[2] * y[6];
  jac[7 * 9 + 8] = k[5] + k[6] + 2 * k[2] * y[6];
  jac[8 * 9 + 8] = -k[5] - k[6] - k[2] * y[6];
}
static int matrix_cases(const char *test)
{
  unsigned char pattern[81] = {0};
  for (size_t i = 0; i < 23; ++i)
    pattern[positions[i][0] * 9 + positions[i][1]] = 1;
  ck_model *m;
  CHECK(ck_linear_model_create(9, pattern, &m) == CK_OK);
  size_t nnz = ck_jacobian_size(m);
  double *a = calloc(nnz, sizeof(double));
  CHECK(a);
  double x[23], b[23] = {0};
  for (size_t i = 0; i < 23; ++i)
    a[ck_jacobian_index(m, positions[i][0], positions[i][1])] = 1;
  for (size_t j = 0; j < 23; ++j)
    x[j] = .5;
  if (IS("lin_solve"))
  {
    memcpy(b, x, 9 * sizeof(double));
    ck_linear_solve(m, 1, a, b);
    double expected[] = {.5, .5, .5, .5, .5, -.5, -1, .5, 0};
    for (size_t j = 0; j < 9; ++j)
      CHECK(b[j] == expected[j]);
  }
  if (IS("factored_alpha_minus_jac"))
  {
    for (size_t i = 0; i < nnz; ++i)
      a[i] = -a[i];
    ck_alpha_minus_jacobian(m, 1, a, 2);
    ck_factor(m, 1, a, CK_LU_MOZART);
    double expected[] = {1, -1, -1, -1, 1,  1,  1, 1,  -1, -1, 1,   -1,
                         1, -1, -1, -1, -2, -1, 3, -1, -2, -3, .125};
    for (size_t i = 0; i < 23; ++i)
    {
      double value = a[ck_jacobian_index(m, positions[i][0], positions[i][1])];
      if (positions[i][0] == positions[i][1])
        value = 1 / value;
      NEAR(value, expected[i], .01);
    }
  }
  if (IS("dforce_dy_times_vector"))
  {
    for (size_t j = 0; j < 5; ++j)
      a[ck_jacobian_index(m, j, j)] = 0;
    ck_matrix_vector(m, 1, a, x, b);
    double expected[23] = {0, 0, 0, 0, 0, 2, 3, 2, 2};
    for (size_t j = 0; j < 23; ++j)
      NEAR(b[j], expected[j], .01);
  }
  free(a);
  ck_model_free(m);
  return 0;
}
static int chapman_unit(const char *test)
{
  if (IS("lin_solve") || IS("factored_alpha_minus_jac") || IS("dforce_dy_times_vector"))
    return matrix_cases(test);
  ck_mechanism d;
  ck_reaction rx[7];
  chapman_description(&d, rx, 0);
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  CHECK(ck_state_create(m, 1, &s) == CK_OK);
  ck_rosenbrock_parameters p;
  CHECK(ck_rosenbrock_parameters_init(3, &p) == CK_OK);
  if (IS("DefaultConstructor"))
    CHECK(p.stages == 3);
  else if (IS("ReactionNames"))
    CHECK(m->nr == 7);
  else if (IS("PhotolysisNames"))
    CHECK(m->nc == 3);
  else if (IS("SpeciesNames"))
    CHECK(m->n == 9);
  else if (IS("simple_force") || IS("smaller_force"))
  {
    for (size_t r = 0; r < 7; ++r)
      s->rates[r] = IS("simple_force") ? 1 : 3e-8;
    for (size_t j = 0; j < 9; ++j)
      s->y[j] = IS("simple_force") ? 1 : 5e-6;
    double f[9], expected[] = {0, 0, 0, 0, 0, -1, 3, 2, -2},
                 small[] = {0, 0, 0, 0, 0, 1.49e-13, 4.55e-13, 1.5e-13, -3e-13};
    ck_forcing(s, s->y, f);
    for (size_t j = 0; j < 9; ++j)
    {
      if (IS("simple_force"))
        CHECK(f[j] == expected[j]);
      else if (j < 5)
        CHECK(f[j] == 0);
      else
        NEAR(f[j], small[j], .01);
    }
  }
  else
  {
    double initial[] = {1, .392, .0169, 0, 32.9, 0, 0, 8.84, 0}, scale = IS("solve_10_times_larger")    ? 10
                                                                         : IS("solve_10_times_smaller") ? .1
                                                                                                        : 1;
    for (size_t j = 0; j < 9; ++j)
      s->y[j] = IS("solve_with_random_number_densities") ? .01 + .1 * j : initial[j] * scale;
    s->conditions[0] = (ck_conditions){273.15, 100000, 2.7e19};
    double photos[] = {1e-4, 1e-5, 1e-6};
    for (size_t j = 0; j < 3; ++j)
      CHECK(ck_state_set_custom(s, photo[j], photos + j, 1) == CK_OK);
    CHECK(ck_update_rates(s) == CK_OK);
    s->relative_tolerance = 1e-4;
    p.h_start = 1e-6;
    p.max_steps = 100;
    ck_stats stats;
    CHECK(ck_solve(s, 1, &p, &stats) == CK_OK);
    if (!IS("solve_with_random_number_densities"))
    {
      double expected[] = {1, .392, .0169, 0, 32.9, 1.8039e-41, .00176789, 8.83912, 4.5031e-36};
      for (size_t j = 0; j < 9; ++j)
        NEAR(s->y[j], expected[j] * scale, j == 5 || j == 8 ? 1e-4 : 1e-4);
    }
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
static int chapman_integration(void)
{
  ck_mechanism d;
  ck_reaction rx[7];
  chapman_description(&d, rx, 1);
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create(&d, &m) == CK_OK);
  CHECK(ck_state_create(m, 1, &s) == CK_OK);
  double initial[] = {.1, .1, .1, .2, .2, .2, .3, .3, .3};
  memcpy(s->y, initial, sizeof initial);
  s->conditions[0] = (ck_conditions){2, 3, 0};
  double photos[] = {.1, .2, .3};
  for (size_t j = 0; j < 3; ++j)
    CHECK(ck_state_set_custom(s, photo[j], photos + j, 1) == CK_OK);
  size_t indices[9];
  double original[9];
  for (size_t j = 0; j < 9; ++j)
  {
    CHECK(ck_species_index(m, chapman_species[j].name, indices + j) == CK_OK);
    original[j] = s->y[indices[j]];
  }
  double oxygen = original[6] + original[5] + 2 * original[7] + 3 * original[8];
  ck_rosenbrock_parameters p;
  CHECK(ck_rosenbrock_parameters_init(3, &p) == CK_OK);
  for (size_t t = 0; t < 100; ++t)
  {
    CHECK(ck_update_rates(s) == CK_OK);
    ck_stats stats;
    CHECK(ck_solve(s, 30, &p, &stats) == CK_OK);
    for (size_t j = 0; j < 9; ++j)
    {
      CHECK(isfinite(s->y[j]) && s->y[j] >= 0);
      if (j < 5)
        NEAR(s->y[indices[j]], original[j], 1e-8 * fabs(original[j]));
    }
    double atoms = s->y[indices[6]] + s->y[indices[5]] + 2 * s->y[indices[7]] + 3 * s->y[indices[8]];
    NEAR(atoms, oxygen, 1e-8 * fabs(oxygen));
  }
  CHECK(s->y[indices[7]] < .5 * original[7] && s->y[indices[6]] > original[6] &&
        s->y[indices[6]] > .5 * oxygen);
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
static int chapman_regression(const char *test)
{
  ck_mechanism d;
  ck_reaction rx[7];
  chapman_description(&d, rx, 0);
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  CHECK(ck_state_create(m, 3, &s) == CK_OK);
  double temperature[] = {284.19, 215.02, 299.31}, pressure[] = {101245, 100789.2, 101398},
         photos[][3] = {{1e-4, 1e-5, 1e-6}, {3.2e-4, 7.3e-5, 3.2e-6}, {5.2e-4, 8.2e-5, 4.6e-6}};
  for (size_t c = 0; c < 3; ++c)
  {
    s->conditions[c] = (ck_conditions){temperature[c], pressure[c], 0};
    for (size_t j = 0; j < 3; ++j)
      s->custom[j * 3 + c] = photos[c][j];
  }
  CHECK(ck_update_rates(s) == CK_OK);
  if (strstr(test, "RateConstants"))
  {
    for (size_t c = 0; c < 3; ++c)
    {
      double expected[] = {2.15e-11 * exp(110 / temperature[c]),
                           3.3e-11 * exp(55 / temperature[c]),
                           8e-12 * exp(-2060 / temperature[c]),
                           6e-34 * pow(temperature[c] / 300, 2.4),
                           photos[c][0],
                           photos[c][1],
                           photos[c][2]};
      for (size_t r = 0; r < 7; ++r)
        CHECK(s->rates[r * 3 + c] == expected[r]);
    }
  }
  else if (strstr(test, "Forcing") || strstr(test, "Jacobian"))
  {
    for (size_t j = 0; j < 9; ++j)
      for (size_t c = 0; c < 3; ++c)
        s->y[j * 3 + c] = exp(-2 + .43 * j - .2 * c);
    for (size_t r = 0; r < 7; ++r)
      for (size_t c = 0; c < 3; ++c)
        s->rates[r * 3 + c] = exp(-2 + .61 * r - .3 * c);
    double f[27];
    double *j = calloc(ck_jacobian_size(m) * 3, sizeof(double));
    CHECK(j);
    ck_forcing(s, s->y, f);
    ck_jacobian(s, s->y, j);
    for (size_t c = 0; c < 3; ++c)
    {
      double y[9], k[7], ef[9], ej[81];
      for (size_t v = 0; v < 9; ++v)
        y[v] = s->y[v * 3 + c];
      for (size_t r = 0; r < 7; ++r)
        k[r] = s->rates[r * 3 + c];
      chapman_expected(k, y, ef, ej);
      for (size_t v = 0; v < 9; ++v)
        NEAR(f[v * 3 + c], ef[v], (fabs(f[v * 3 + c]) + fabs(ef[v])) * 1e-8 + 1e-12);
      for (size_t p = 0; p < 23; ++p)
      {
        size_t row = positions[p][0], col = positions[p][1], off = ck_jacobian_index(m, row, col);
        NEAR(j[off * 3 + c], -ej[row * 9 + col],
             (fabs(j[off * 3 + c]) + fabs(ej[row * 9 + col])) * 1e-8 + 1e-12);
      }
    }
    free(j);
  }
  else
  {
    int method = IS("TwoStageSolve")      ? 2
                 : IS("FourStageSolve")   ? 4
                 : IS("FourStageDASolve") ? 43
                 : IS("SixStageDASolve")  ? 6
                                          : 3;
    ck_rosenbrock_parameters p, baseline;
    CHECK(ck_rosenbrock_parameters_init(method, &p) == CK_OK);
    CHECK(ck_rosenbrock_parameters_init(3, &baseline) == CK_OK);
    baseline.h_start = 1e-6;
    baseline.max_steps = 100;
    double initial[27], atol = 0;
    for (size_t v = 0; v < 9; ++v)
      for (size_t c = 0; c < 3; ++c)
      {
        double value = exp(-2 + .27 * v + .31 * c);
        initial[v * 3 + c] = s->y[v * 3 + c] = value;
        atol = fmax(atol, value);
      }
    atol *= 1e-12;
    ck_stats stats;
    CHECK(ck_solve(s, 500, &p, &stats) == CK_OK);
    for (size_t c = 0; c < 3; ++c)
    {
      ck_state *reference;
      CHECK(ck_state_create(m, 1, &reference) == CK_OK);
      for (size_t v = 0; v < 9; ++v)
        reference->y[v] = initial[v * 3 + c];
      reference->conditions[0] = s->conditions[c];
      reference->relative_tolerance = 1e-4;
      for (size_t j = 0; j < 3; ++j)
        reference->custom[j] = photos[c][j];
      CHECK(ck_update_rates(reference) == CK_OK);
      CHECK(ck_solve(reference, 500, &baseline, &stats) == CK_OK);
      for (size_t v = 0; v < 9; ++v)
        NEAR(s->y[v * 3 + c], reference->y[v],
             (fabs(s->y[v * 3 + c]) + fabs(reference->y[v])) * (method == 2 ? 1e-2 : 1e-4) + atol);
      ck_state_free(reference);
    }
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
int test_chapman(const char *id, const char *suite, const char *test)
{
  if (!strcmp(suite, "ChapmanMechanismHardCodedAndGeneral") || !strcmp(suite, "RegressionChapmanODESolver"))
    return chapman_unit(test);
  if (!strncmp(id, "chapman_integration:", 20))
    return chapman_integration();
  if (!strncmp(id, "regression_test_", 16))
    return chapman_regression(test);
  return 2;
}
