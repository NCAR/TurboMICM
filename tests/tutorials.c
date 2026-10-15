/* C ports of the eight tutorial programs and independent configuration checks. */
#include "test.h"
#include "tutorials.h"
#include <omp.h>
extern const ck_mechanism tutorial_no, tutorial_user, tutorial_robertson;
const ck_mechanism *tutorial_hand_rates(int user)
{
  static const ck_species species[] = {
      {.name = "A"}, {.name = "B"}, {.name = "C", .molecular_weight = 0.025, .diffusion = 2.3e2},
      {.name = "D"}, {.name = "E"}, {.name = "F"},
      {.name = "G"}};
  static const size_t reactants[][2] = {{0, 0}, {1, 0}, {1, 0}, {2, 0}, {3, 0},
                                        {4, 4}, {5, 0}, {2, 0}, {0, 0}, {1, 0}};
  static const ck_product products[][2] = {{{1, 1}}, {{2, 1}}, {{3, 1}}, {{4, 1}}, {{5, 2}},
                                           {{6, 1}}, {{6, 1}}, {{6, 1}}, {{0, 1}}};
  static ck_reaction reactions[10];
  static ck_mechanism mechanism;
  for (size_t j = 0; j < 10; ++j)
    reactions[j] = (ck_reaction){.nreactants = j == 8   ? 0
                                               : j == 5 ? 2
                                                        : 1,
                                 .reactants = reactants[j],
                                 .nproducts = j == 9 ? 0 : 1,
                                 .products = products[j],
                                 .rate = ck_rate_defaults(CK_ARRHENIUS)};
  reactions[0].rate.A = 2.15e-4;
  reactions[0].rate.C = 110;
  reactions[1].rate = ck_rate_defaults(CK_BRANCH_ALKOXY);
  reactions[1].rate.X = 1.2;
  reactions[1].rate.Y = 204.3;
  reactions[1].rate.a0 = 1e-3;
  reactions[1].rate.n = 2;
  reactions[2].rate = reactions[1].rate;
  reactions[2].rate.type = CK_BRANCH_NITRATE;
  reactions[3].rate = ck_rate_defaults(CK_SURFACE);
  reactions[3].rate.label = "C";
  reactions[3].rate.molecular_weight = 0.025;
  reactions[3].rate.diffusion = 2.3e2;
  reactions[3].rate.probability = 0.9;
  double factor = 1e-6 * CK_AVOGADRO;
  reactions[4].rate = ck_rate_defaults(CK_TERNARY);
  reactions[4].rate.k0_A = 1.2;
  reactions[4].rate.k0_B = 2.3;
  reactions[4].rate.k0_C = 302.3;
  reactions[4].rate.kinf_A = 2.6 / factor;
  reactions[4].rate.kinf_B = -3.1;
  reactions[4].rate.kinf_C = 402.1;
  reactions[4].rate.Fc = 0.9;
  reactions[4].rate.N = 1.2;
  reactions[5].rate = ck_rate_defaults(CK_TROE);
  reactions[5].rate.k0_A = 1.8e-30 * factor * factor;
  reactions[5].rate.k0_B = -3;
  reactions[5].rate.kinf_A = 2.8e-11 * factor;
  reactions[6].rate = ck_rate_defaults(CK_TUNNELING);
  reactions[6].rate.A = 1.2;
  reactions[6].rate.B = 2.3;
  reactions[6].rate.C = 302.3;
  const char *labels[] = {"my photolysis rate", "my emission rate", "my loss rate"};
  for (size_t j = 7; j < 10; ++j)
  {
    reactions[j].rate = ck_rate_defaults(CK_USER);
    reactions[j].rate.label = labels[j - 7];
  }
  mechanism = (ck_mechanism){.name = "rates by hand",
                             .nspecies = 7,
                             .species = species,
                             .nreactions = user ? 10 : 7,
                             .reactions = reactions};
  return &mechanism;
}
const ck_mechanism *tutorial_hand_robertson(void)
{
  static const ck_species species[] = {{.name = "A"}, {.name = "B"}, {.name = "C"}};
  static const size_t reactants[][2] = {{0, 0}, {1, 1}, {1, 2}};
  static const ck_product products[][2] = {{{1, 1}, {0, 0}}, {{1, 1}, {2, 1}}, {{0, 1}, {2, 1}}};
  static ck_reaction reactions[3];
  static const char *labels[] = {"r1", "r2", "r3"};
  static ck_mechanism mechanism;
  for (size_t j = 0; j < 3; ++j)
  {
    reactions[j] = (ck_reaction){.nreactants = j ? 2 : 1,
                                 .reactants = reactants[j],
                                 .nproducts = j ? 2 : 1,
                                 .products = products[j],
                                 .rate = ck_rate_defaults(CK_USER)};
    reactions[j].rate.label = labels[j];
  }
  mechanism = (ck_mechanism){.name = "Robertson by hand",
                             .nspecies = 3,
                             .species = species,
                             .nreactions = 3,
                             .reactions = reactions};
  return &mechanism;
}
static void header(FILE *out, const ck_mechanism *desc, size_t cells)
{
  fprintf(out, "%6s", "time");
  if (cells > 1)
    fprintf(out, ",%6s", "grid");
  for (size_t i = 0; i < desc->nspecies; ++i)
    if (!desc->species[i].parameterized)
      fprintf(out, ",%11s", desc->species[i].name);
  fputc('\n', out);
}
static void print_state(FILE *out, const ck_state *s, const ck_mechanism *desc, double time)
{
  for (size_t c = 0; c < s->cells; ++c)
  {
    fprintf(out, "%6g,", time);
    if (s->cells > 1)
      fprintf(out, "%6zu,", c);
    int first = 1;
    for (size_t j = 0; j < desc->nspecies; ++j)
      if (!desc->species[j].parameterized)
      {
        size_t index;
        if (ck_species_index(s->model, desc->species[j].name, &index))
          continue;
        fprintf(out, "%s%11.2e", first ? "" : ",", s->y[index * s->cells + c]);
        first = 0;
      }
    fputc('\n', out);
  }
}
static int initialize_robertson(ck_state *s, int configuration, int multi)
{
  double init[] = {1, 2, 0.5}, zero[3] = {0}, k1[] = {0.04, 0.04, 0.04}, k2[] = {3e7, 3e7, 3e7},
         k3[] = {1e4, 1e4, 1e4};
  if (!multi)
    init[1] = init[2] = 1;
  if (configuration)
    for (size_t c = 0; c < s->cells; ++c)
    {
      k2[c] /= 1e-6 * CK_AVOGADRO;
      k3[c] /= 1e-6 * CK_AVOGADRO;
    }
  CHECK(ck_state_set(s, "A", init, s->cells) == CK_OK);
  CHECK(ck_state_set(s, "B", zero, s->cells) == CK_OK);
  CHECK(ck_state_set(s, "C", zero, s->cells) == CK_OK);
  CHECK(ck_state_set_custom(s, "r1", k1, s->cells) == CK_OK);
  CHECK(ck_state_set_custom(s, "r2", k2, s->cells) == CK_OK);
  CHECK(ck_state_set_custom(s, "r3", k3, s->cells) == CK_OK);
  for (size_t c = 0; c < s->cells; ++c)
    s->conditions[c] = (ck_conditions){272.5, 101253.3, 1e6};
  return ck_update_rates(s);
}
static int solve_interval(ck_state *s, double dt, const ck_rosenbrock_parameters *p, ck_stats *sum)
{
  double elapsed = 0;
  size_t attempts = 0;
  while (elapsed < dt)
  {
    ck_stats stats;
    ck_status status = ck_solve(s, dt - elapsed, p, &stats);
    if (sum)
      ck_stats_add(sum, &stats);
    CHECK(status == CK_OK || status == CK_MAX_STEPS);
    CHECK(stats.final_time > 0 && ++attempts < 1000);
    elapsed += stats.final_time;
  }
  return 0;
}
static void print_stats(FILE *out, const ck_stats *s)
{
  fprintf(out,
          "accepted: %llu\nfunction_calls: %llu\njacobian_updates: %llu\nnumber_of_steps: %llu\naccepted: "
          "%llu\nrejected: %llu\ndecompositions: %llu\nsolves: %llu\n",
          (unsigned long long)s->accepted, (unsigned long long)s->function_calls,
          (unsigned long long)s->jacobian_updates, (unsigned long long)s->internal_steps,
          (unsigned long long)s->accepted, (unsigned long long)s->rejected,
          (unsigned long long)s->decompositions, (unsigned long long)s->solves);
}
static int rate_tutorial(FILE *out, int user, int configuration)
{
  const ck_mechanism *desc =
      configuration ? (user ? &tutorial_user : &tutorial_no) : tutorial_hand_rates(user);
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create(desc, &m) == CK_OK);
  CHECK(ck_state_create(m, 1, &s) == CK_OK);
  s->conditions[0] = (ck_conditions){287.45, 101319.9, 0};
  ck_ideal_air_density(&s->conditions[0]);
  double one = 1;
  CHECK(ck_state_set(s, "A", &one, 1) == CK_OK);
  double radius = 1e-7, number = 2.5e6;
  CHECK(ck_state_set_custom(s, configuration ? "C surface.effective radius [m]" : "C.effective radius [m]",
                            &radius, 1) == CK_OK);
  CHECK(ck_state_set_custom(s,
                            configuration ? "C surface.particle number concentration [# m-3]"
                                          : "C.particle number concentration [# m-3]",
                            &number, 1) == CK_OK);
  if (user)
  {
    double emission = 1e-20, loss = 1e-23;
    if (configuration)
      emission *= 1e-6 * CK_AVOGADRO;
    CHECK(ck_state_set_custom(s, "my emission rate", &emission, 1) == CK_OK);
    CHECK(ck_state_set_custom(s, "my loss rate", &loss, 1) == CK_OK);
  }
  header(out, desc, 1);
  print_state(out, s, desc, 0);
  ck_rosenbrock_parameters p;
  CHECK(ck_rosenbrock_parameters_init(3, &p) == CK_OK);
  double photo = 1e-10;
  for (size_t step = 0; step < 10; ++step)
  {
    if (user)
      CHECK(ck_state_set_custom(s, "my photolysis rate", &photo, 1) == CK_OK);
    CHECK(ck_update_rates(s) == CK_OK);
    CHECK(solve_interval(s, 500, &p, NULL) == 0);
    print_state(out, s, desc, 500 * (step + 1));
    photo *= 1.5;
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
int run_tutorial(const char *name, FILE *out, int configuration)
{
  if (!strcmp(name, "rate_constants_no_user_defined_example_by_hand"))
    return rate_tutorial(out, 0, configuration);
  if (!strcmp(name, "rate_constants_user_defined_example_by_hand"))
    return rate_tutorial(out, 1, configuration);
  int readme = !strcmp(name, "README_example"), results = !strcmp(name, "solve_results"),
      multi = !strcmp(name, "multiple_grid_cells"), config = !strcmp(name, "solver_configuration");
  int vector = !strcmp(name, "vectorized_matrix_solver"), threads = !strcmp(name, "openmp_tutorial");
  if (!readme && !results && !multi && !config && !vector && !threads)
    return 2;
  const ck_mechanism *desc = configuration ? &tutorial_robertson : tutorial_hand_robertson();
  static const ck_species species[] = {{.name = "Foo"}, {.name = "Bar"}, {.name = "Baz"}};
  static const size_t r0[] = {0}, r1[] = {0, 1};
  static const ck_product p0[] = {{1, 0.8}, {2, 0.2}}, p1[] = {{2, 1}};
  ck_reaction reactions[2];
  ck_mechanism readme_desc = {
      .name = "README", .nspecies = 3, .species = species, .nreactions = 2, .reactions = reactions};
  if (readme)
  {
    reactions[0] = (ck_reaction){.nreactants = 1,
                                 .reactants = r0,
                                 .nproducts = 2,
                                 .products = p0,
                                 .rate = ck_rate_defaults(CK_ARRHENIUS)};
    reactions[0].rate.A = 1e-3;
    reactions[1] = (ck_reaction){.nreactants = 2,
                                 .reactants = r1,
                                 .nproducts = 1,
                                 .products = p1,
                                 .rate = ck_rate_defaults(CK_ARRHENIUS)};
    reactions[1].rate.A = 1e-5;
    reactions[1].rate.C = 110;
    desc = &readme_desc;
  }
  ck_model *m;
  CHECK(ck_model_create(desc, &m) == CK_OK);
  if (threads)
  {
    ck_state *states[3] = {0};
    int statuses[3] = {0};
    for (size_t i = 0; i < 3; ++i)
      CHECK(ck_state_create(m, 1, &states[i]) == CK_OK);
#pragma omp parallel for num_threads(3) schedule(static)
    for (int i = 0; i < 3; ++i)
    {
      ck_rosenbrock_parameters p;
      ck_rosenbrock_parameters_init(3, &p);
      statuses[i] = initialize_robertson(states[i], configuration, 0);
      for (size_t j = 0; j < 10 && !statuses[i]; ++j)
        statuses[i] = solve_interval(states[i], 200, &p, NULL);
    }
    for (size_t i = 0; i < 3; ++i)
      fprintf(out, "Running solver on thread %zu\n", i);
    for (size_t i = 0; i < 3; ++i)
    {
      CHECK(!statuses[i]);
      fprintf(out, "Thread %zu\n%10s,%10s,%10s\n", i + 1, "A", "B", "C");
      size_t a, b, c;
      CHECK(ck_species_index(m, "A", &a) == CK_OK);
      CHECK(ck_species_index(m, "B", &b) == CK_OK);
      CHECK(ck_species_index(m, "C", &c) == CK_OK);
      fprintf(out, "%10.2e,%10.2e,%10.2e\n", states[i]->y[a], states[i]->y[b], states[i]->y[c]);
      ck_state_free(states[i]);
    }
    ck_model_free(m);
    return 0;
  }
  if (vector)
  {
    /* The original program's while condition skips integration; preserve its
       printed behavior while separate integration tests exercise both layouts. */
    double values[3][3] = {{1.1, 1.2, 1.3}, {2.1, 2.2, 2.3}, {3.1, 3.2, 3.3}};
    for (size_t c = 0; c < 3; ++c)
      for (size_t j = 0; j < 3; ++j)
        fprintf(out, "%g\n", values[c][j]);
    fputc('\n', out);
    for (size_t j = 0; j < 3; ++j)
      for (size_t c = 0; c < 3; ++c)
        fprintf(out, "%g\n", values[c][j]);
    fputc('\n', out);
    for (size_t c = 0; c < 3; ++c)
    {
      fprintf(out, "Cell %zu\n%10s%20s%20s\n", c, "Species", "Regular", "Vectorized");
      for (size_t j = 0; j < 3; ++j)
        fprintf(out, "%10s%20g%20g\n", desc->species[j].name, values[c][j], values[c][j]);
      fputc('\n', out);
    }
    fputs("\nSparse matrix standard ordering elements\n", out);
    const size_t rows[] = {0, 2, 2, 3}, cols[] = {1, 1, 3, 2};
    for (size_t c = 0; c < 3; ++c)
      for (size_t j = 0; j < 4; ++j)
        fprintf(out, "%zu.%zu.%zu\n", c, rows[j], cols[j]);
    fputs("\nSparse matrix vector ordering elements\n", out);
    for (size_t j = 0; j < 4; ++j)
      for (size_t c = 0; c < 3; ++c)
        fprintf(out, "%zu.%zu.%zu\n", c, rows[j], cols[j]);
    ck_model_free(m);
    return 0;
  }
  const int methods[] = {2, 3, 4, 43, 6};
  const char *method_names[] = {
      "Two stages: ", "Three stages: ", "Four stages: ", "Four stages differential algebraic: ",
      "Six stages differential algebraic: "};
  for (size_t method = 0; method < (config ? 5u : 1u); ++method)
  {
    ck_state *s;
    CHECK(ck_state_create(m, multi || results ? 3 : 1, &s) == CK_OK);
    if (readme)
    {
      s->conditions[0] = (ck_conditions){287.45, 101319.9, 0};
      ck_ideal_air_density(&s->conditions[0]);
      double init = 20;
      CHECK(ck_state_set(s, "Foo", &init, 1) == CK_OK);
    }
    else
      CHECK(initialize_robertson(s, configuration, multi || results) == CK_OK);
    ck_rosenbrock_parameters p;
    CHECK(ck_rosenbrock_parameters_init(config ? methods[method] : 3, &p) == CK_OK);
    if (config)
    {
      if (method)
        fputc('\n', out);
      fprintf(out, "%s\n", method_names[method]);
    }
    if (!results)
    {
      header(out, desc, s->cells);
      if (!readme)
        print_state(out, s, desc, 0);
    }
    ck_stats sum = {0};
    size_t nsteps = results ? 1 : 10;
    for (size_t step = 0; step < nsteps; ++step)
    {
      CHECK(ck_update_rates(s) == CK_OK);
      CHECK(solve_interval(s, readme ? 500 : 200, &p, &sum) == 0);
      if (!results)
        print_state(out, s, desc, readme ? step * 500 : (step + 1) * 200);
    }
    if (results)
    {
      fprintf(out, "Solver state: Converged\n");
      print_stats(out, &sum);
      fputs("final simulation time: 200\n", out);
    }
    if (config)
    {
      fputs("Total solve time: 0 nanoseconds\n", out);
      print_stats(out, &sum);
    }
    ck_state_free(s);
  }
  ck_model_free(m);
  return 0;
}
static int compare_output(FILE *actual, FILE *expected, int threaded)
{
  rewind(actual);
  rewind(expected);
  char a[1024], b[1024];
  if (threaded)
  {
    while (fgets(a, sizeof a, actual) && strncmp(a, "Thread 1\n", 9))
    {
    }
    while (fgets(b, sizeof b, expected) && strncmp(b, "Thread 1\n", 9))
    {
    }
    CHECK(!strcmp(a, b));
  }
  for (;;)
  {
    char *aa = fgets(a, sizeof a, actual), *bb = fgets(b, sizeof b, expected);
    CHECK((aa == NULL) == (bb == NULL));
    if (!aa)
      break;
    if (!strncmp(a, "Total solve time:", 17) && !strncmp(b, "Total solve time:", 17))
      continue;
    if (strcmp(a, b))
    {
      fprintf(stderr, "tutorial output mismatch\nactual: %sexpected: %s", a, b);
      return 1;
    }
  }
  return 0;
}
int test_tutorial(const char *name)
{
  if (!strcmp(name, "rate_constants_no_user_defined_example_by_hand"))
    CHECK(test_configuration_equivalence(0) == 0);
  if (!strcmp(name, "rate_constants_user_defined_example_by_hand"))
    CHECK(test_configuration_equivalence(1) == 0);
  if (!strcmp(name, "solver_configuration"))
    CHECK(test_configuration_equivalence(2) == 0);
  char path[512];
  snprintf(path, sizeof path, "tests/expected/%s.txt", name);
  FILE *expected = fopen(path, "r");
  if (!expected)
    return 2;
  FILE *actual = tmpfile();
  CHECK(actual);
  int status = run_tutorial(name, actual, 0);
  if (!status)
    status = compare_output(actual, expected, !strcmp(name, "openmp_tutorial"));
  if (!status && (!strncmp(name, "rate_constants_", 15) || !strcmp(name, "solver_configuration")))
  {
    FILE *configured = tmpfile();
    CHECK(configured);
    status = run_tutorial(name, configured, 1);
    if (!status)
      status = compare_output(configured, actual, 0);
    fclose(configured);
  }
  fclose(actual);
  fclose(expected);
  return status;
}
