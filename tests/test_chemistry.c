#include "test.h"
#include "tutorials.h"
#include <omp.h>
#include "random.h"

static int shared_solver(void)
{
  const ck_mechanism *d = tutorial_hand_robertson();
  ck_model *m;
  CHECK(ck_model_create_with_order(d, 0, &m) == CK_OK);
  ck_state *states[8] = {0};
  double results[8][3] = {{0}};
  int failures[8] = {0};
  ck_rosenbrock_parameters p;
  CHECK(ck_rosenbrock_parameters_init(3, &p) == CK_OK);
  for (size_t i = 0; i < 8; ++i)
    CHECK(ck_state_create(m, 1, states + i) == CK_OK);
#pragma omp parallel for num_threads(8) schedule(static)
  for (size_t thread = 0; thread < 8; ++thread)
  {
    ck_state *s = states[thread];
    s->y[0] = 1;
    s->conditions[0] = (ck_conditions){272.5, 101253.3, 1e6};
    const double rates[] = {.04, 3e7, 1e4};
    const char *names[] = {"r1", "r2", "r3"};
    for (size_t k = 0; k < 3; ++k)
      if (ck_state_set_custom(s, names[k], rates + k, 1))
        failures[thread] = 1;
    for (size_t step = 0; step < 10 && !failures[thread]; ++step)
    {
      if (ck_update_rates(s))
      {
        failures[thread] = 1;
        break;
      }
      double elapsed = 0;
      while (elapsed < 200)
      {
        ck_stats stats;
        ck_status status = ck_solve(s, 200 - elapsed, &p, &stats);
        if ((status != CK_OK && status != CK_MAX_STEPS) || stats.final_time <= 0)
        {
          failures[thread] = 1;
          break;
        }
        elapsed += stats.final_time;
      }
    }
    memcpy(results[thread], s->y, sizeof results[thread]);
  }
  for (size_t thread = 0; thread < 8; ++thread)
  {
    CHECK(!failures[thread]);
    for (size_t j = 0; j < 3; ++j)
      CHECK(results[thread][j] == results[0][j]);
    ck_state_free(states[thread]);
  }
  ck_model_free(m);
  return 0;
}
static int tracers(void)
{
  ck_species sp[] = {{.name = "a"}, {.name = "b"}, {.name = "c"}, {.name = "irr_1"}, {.name = "irr_2"}};
  size_t reactants[] = {0, 1};
  ck_product products[][2] = {{{1, 1}, {3, 1}}, {{2, 1}, {4, 1}}};
  ck_reaction rx[2];
  for (size_t r = 0; r < 2; ++r)
    rx[r] = (ck_reaction){.nreactants = 1,
                          .reactants = reactants + r,
                          .nproducts = 2,
                          .products = products[r],
                          .rate = ck_rate_defaults(r ? CK_USER : CK_ARRHENIUS)};
  rx[0].rate.A = 2.15e-11;
  rx[0].rate.C = 110;
  rx[1].rate.label = "r2";
  ck_mechanism d = {.name = "integrated reaction rates",
                    .nspecies = 5,
                    .species = sp,
                    .nreactions = 2,
                    .reactions = rx,
                    .temperature = 273,
                    .pressure = 1000};
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create(&d, &m) == CK_OK);
  CHECK(ck_state_create(m, 1, &s) == CK_OK);
  size_t a, r1, r2;
  CHECK(ck_species_index(m, "a", &a) == CK_OK);
  CHECK(ck_species_index(m, "irr_1", &r1) == CK_OK);
  CHECK(ck_species_index(m, "irr_2", &r2) == CK_OK);
  s->y[a] = 1;
  ck_rosenbrock_parameters p;
  CHECK(ck_rosenbrock_parameters_init(3, &p) == CK_OK);
  for (size_t step = 0; step < 100; ++step)
  {
    double rate = step > 50 ? 0 : 1;
    CHECK(ck_state_set_custom(s, "r2", &rate, 1) == CK_OK);
    CHECK(ck_update_rates(s) == CK_OK);
    double first = s->y[r1], second = s->y[r2];
    ck_stats stats;
    CHECK(ck_solve(s, 30, &p, &stats) == CK_OK);
    CHECK(s->y[r1] >= first && s->y[r2] >= second);
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
static int terminator(size_t cells, int combined)
{
  ck_species sp[] = {{.name = "Cl2", .absolute_tolerance = 1e-20},
                     {.name = "Cl", .absolute_tolerance = 1e-20},
                     {.name = "A"},
                     {.name = "B"},
                     {.name = "C"}};
  size_t reactants[][2] = {{0, 0}, {1, 1}, {2, 0}, {3, 3}, {3, 4}};
  ck_product products[][2] = {
      {{1, 2}, {0, 0}}, {{0, 1}, {0, 0}}, {{3, 1}, {0, 0}}, {{3, 1}, {4, 1}}, {{2, 1}, {4, 1}}};
  const char *labels[] = {"terminator_k1", NULL, "robertson_r1", "robertson_r2", "robertson_r3"};
  ck_reaction rx[5];
  for (size_t r = 0; r < 5; ++r)
  {
    rx[r] = (ck_reaction){.nreactants = r == 1 || r >= 3 ? 2 : 1,
                          .reactants = reactants[r],
                          .nproducts = r >= 3 ? 2 : 1,
                          .products = products[r],
                          .rate = ck_rate_defaults(r == 1 ? CK_ARRHENIUS : CK_USER)};
    rx[r].rate.label = labels[r];
  }
  ck_product mass[] = {{2, 1}, {3, 1}, {4, 1}};
  ck_constraint constraint = {.type = CK_LINEAR_CONSTRAINT,
                              .name = "mass",
                              .algebraic_species = 4,
                              .nreactants = 3,
                              .reactants = mass,
                              .constant = 1};
  ck_mechanism d = {.name = "terminator",
                    .nspecies = combined ? 5 : 2,
                    .species = sp,
                    .nreactions = combined ? 5 : 2,
                    .reactions = rx,
                    .nconstraints = combined ? 1 : 0,
                    .constraints = &constraint,
                    .temperature = 298,
                    .pressure = 101300};
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  CHECK(ck_state_create(m, cells, &s) == CK_OK);
  s->relative_tolerance = 1e-8;
  ck_rosenbrock_parameters p;
  CHECK(ck_rosenbrock_parameters_init(combined ? 43 : 3, &p) == CK_OK);
  p.max_steps = combined ? p.max_steps : 100000;
  if (combined)
  {
    s->y[0] = 1.8e-10;
    s->y[1] = 1.2e-6;
    s->y[2] = 1;
    s->conditions[0].air_density = 42;
    const double rates[] = {.04, 3e7, 1e4};
    for (size_t r = 0; r < 3; ++r)
      CHECK(ck_state_set_custom(s, labels[r + 2], rates + r, 1) == CK_OK);
    CHECK(ck_update_rates(s) == CK_OK);
    double dt = 1, total = 2 * 1.8e-10 + 1.2e-6;
    for (size_t step = 0; step < 12; ++step, dt *= 10)
    {
      double advanced = 0;
      while (advanced < dt)
      {
        ck_stats stats;
        ck_status status = ck_solve(s, dt - advanced, &p, &stats);
        CHECK(status == CK_OK || status == CK_MAX_STEPS);
        CHECK(advanced + stats.final_time > advanced);
        advanced += stats.final_time;
      }
      NEAR(s->y[2] + s->y[3] + s->y[4], 1, 1e-10);
      NEAR(2 * s->y[0] + s->y[1], total, 1e-6 * fabs(total));
      for (size_t j = 0; j < 5; ++j)
        CHECK(isfinite(s->y[j]) && s->y[j] >= 0);
      if (step == 0)
      {
        CHECK(s->y[2] < .99 && s->y[3] > 0 && s->y[4] > 1e-3);
      }
    }
    CHECK(s->y[4] > .95 && s->y[2] < .05 && s->y[0] > .4 * total && s->y[1] < .1 * total);
  }
  else
  {
    uint32_t seed = 739;
    double initial[8] = {0};
    for (size_t c = 0; c < cells; ++c)
    {
      initial[c] = test_lognormal(&seed, -2, 2) * 1e-6;
      initial[cells + c] = test_lognormal(&seed, -2, 2) * 1e-10;
    }
    for (size_t step = 0; step < (size_t)floor(2 * CK_PI / .3); ++step)
    {
      memcpy(s->y, initial, 2 * cells * sizeof(double));
      double k1[4];
      for (size_t c = 0; c < cells; ++c)
      {
        double lat = CK_PI / 180 * (c * (90. / cells)), lon = step * .3;
        k1[c] = fmax(0, sin(lat) * sin(CK_PI * 20 / 180) +
                            cos(lat) * cos(CK_PI * 20 / 180) * cos(lon - CK_PI * 300 / 180));
        s->conditions[c].air_density = 42;
      }
      CHECK(ck_state_set_custom(s, "terminator_k1", k1, cells) == CK_OK);
      CHECK(ck_update_rates(s) == CK_OK);
      ck_stats stats;
      CHECK(ck_solve(s, 30, &p, &stats) == CK_OK);
      for (size_t c = 0; c < cells; ++c)
      {
        double r = k1[c] / 4, cl = initial[cells + c], cl2 = initial[c], cly = cl + 2 * cl2,
               det = sqrt(r * r + 2 * r * cly), e = exp(-4 * det * 30),
               l = det * 30 > 1e-16 ? (1 - e) / det / 30 : 4;
        double change = -l * (cl - det + r) * (cl + det + r) / (1 + e + 30 * l * (cl + r));
        double clref = cl + 30 * change, cl2ref = cl2 - 15 * change;
        NEAR(s->y[cells + c], clref, fabs(clref) * 1e-6 + 1e-14);
        NEAR(s->y[c], cl2ref, fabs(cl2ref) * 1e-6 + 1e-14);
      }
    }
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
int test_chemistry(const char *id, const char *test)
{
  if (!strncmp(id, "openmp:", 7))
    return shared_solver();
  if (!strncmp(id, "integrated_reaction_rates:", 26))
    return tracers();
  if (!strncmp(id, "linear_constraint:", 18))
    return terminator(1, 1);
  if (!strncmp(id, "terminator:", 11))
  {
    for (size_t cells = 1; cells <= 4; ++cells)
      CHECK(terminator(cells, 0) == 0);
    return 0;
  }
  (void)test;
  return 2;
}
