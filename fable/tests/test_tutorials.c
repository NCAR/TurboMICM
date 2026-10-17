/* MICM's eight tutorial programs (test/tutorial/test_*.cpp), ported one function per
 * program. Each port writes what the program prints to std::cout into a
 * buffer and compares it, line by line, with the program's output recorded
 * from MICM's build in tests/data/tutorials/<program>.txt. Only the fields no
 * two runs agree on are left out: the measured solve times of
 * solver_configuration (the line is compared up to its number) and the
 * thread start-up messages of openmp_tutorial, which interleave (the
 * comparison starts at "Thread 1").
 *
 * The two *_by_hand programs build the example configurations' mechanisms in
 * code; their cases also build the mechanisms generated from the
 * configuration files (configs/examples) and check that both give the same
 * rate constants and concentrations to 1e-10 relative. */
#include "harness.h"
#include "rate_constants_no_user_defined.h"
#include "rate_constants_user_defined.h"
#include "tmc_tutorial_ext.h"

#include <stdarg.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* ---- a buffer standing for std::cout ---------------------------------------- */

typedef struct outbuf
{
  char *data;
  size_t length, capacity;
} outbuf;

static void out_reserve(outbuf *o, size_t extra)
{
  if (o->length + extra + 1 <= o->capacity)
    return;
  size_t capacity = o->capacity ? o->capacity : 4096;
  while (capacity < o->length + extra + 1)
    capacity *= 2;
  char *data = realloc(o->data, capacity);
  if (!data)
    abort();
  o->data = data;
  o->capacity = capacity;
  o->data[o->length] = '\0';
}

static void out_printf(outbuf *o, const char *format, ...)
{
  va_list args, copy;
  va_start(args, format);
  va_copy(copy, args);
  const int n = vsnprintf(NULL, 0, format, copy);
  va_end(copy);
  if (n > 0)
  {
    out_reserve(o, (size_t)n);
    vsnprintf(o->data + o->length, (size_t)n + 1, format, args);
    o->length += (size_t)n;
  }
  va_end(args);
}

/* state.PrintHeader() and state.PrintState(time) */
static void out_header(outbuf *o, const tmc_state *state)
{
  const size_t n = tmc_state_format_header(state, NULL, 0);
  out_reserve(o, n);
  tmc_state_format_header(state, o->data + o->length, n + 1);
  o->length += n;
}

static void out_state(outbuf *o, const tmc_state *state, double time)
{
  const size_t n = tmc_state_format_state(state, time, NULL, 0);
  out_reserve(o, n);
  tmc_state_format_state(state, time, o->data + o->length, n + 1);
  o->length += n;
}

static void out_free(outbuf *o)
{
  free(o->data);
  o->data = NULL;
  o->length = o->capacity = 0;
}

/* ---- comparison with the recorded output ------------------------------------ */

static char *read_file(const char *path)
{
  FILE *file = fopen(path, "rb");
  if (!file)
    return NULL;
  size_t capacity = 8192, length = 0;
  char *text = malloc(capacity);
  if (!text)
    abort();
  size_t n;
  while ((n = fread(text + length, 1, capacity - length - 1, file)) > 0)
  {
    length += n;
    if (length + 1 >= capacity)
    {
      capacity *= 2;
      char *grown = realloc(text, capacity);
      if (!grown)
        abort();
      text = grown;
    }
  }
  fclose(file);
  text[length] = '\0';
  return text;
}

/* The next line of *text into line (without its newline); false at the end. */
static bool next_line(const char **text, char *line, size_t size)
{
  if (**text == '\0')
    return false;
  const size_t n = strcspn(*text, "\n");
  snprintf(line, size, "%.*s", (int)(n < size - 1 ? n : size - 1), *text);
  *text += n + ((*text)[n] == '\n' ? 1 : 0);
  return true;
}

/* A line that is compared up to this prefix when it starts with it (its
 * remainder is a measured time). */
static void apply_prefix_filter(char *line, const char *prefix)
{
  if (prefix && strncmp(line, prefix, strlen(prefix)) == 0)
    line[strlen(prefix)] = '\0';
}

/* Compares the program's output with tests/data/tutorials/<name>.txt. When
 * start_marker is given, both texts are compared from the first line equal to
 * it; lines starting with time_prefix are compared up to the prefix. */
static void compare_output(const char *name, const outbuf *actual, const char *start_marker, const char *time_prefix)
{
  char relative[128];
  snprintf(relative, sizeof relative, "tutorials/%s.txt", name);
  char *expected_text = read_file(tmc_test_data_path(relative));
  if (!expected_text)
  {
    tmc_test_fail(__FILE__, __LINE__, "%s: cannot read %s", name, tmc_test_data_path(relative));
    return;
  }
  const char *expected = expected_text, *got = actual->data ? actual->data : "";
  char expected_line[1024], got_line[1024];
  if (start_marker)
  {
    bool found_expected = false, found_got = false;
    while (!found_expected && next_line(&expected, expected_line, sizeof expected_line))
      found_expected = strcmp(expected_line, start_marker) == 0;
    while (!found_got && next_line(&got, got_line, sizeof got_line))
      found_got = strcmp(got_line, start_marker) == 0;
    if (!found_expected || !found_got)
    {
      tmc_test_fail(__FILE__, __LINE__, "%s: \"%s\" not found in the %s output", name, start_marker, found_expected ? "program's" : "recorded");
      free(expected_text);
      return;
    }
  }
  int line_number = start_marker ? 1 : 0;
  for (;;)
  {
    const bool more_expected = next_line(&expected, expected_line, sizeof expected_line);
    const bool more_got = next_line(&got, got_line, sizeof got_line);
    ++line_number;
    if (!more_expected && !more_got)
      break;
    if (more_expected != more_got)
    {
      tmc_test_fail(__FILE__, __LINE__, "%s: line %d: %s", name, line_number, more_expected ? "program output ends early" : "program output has extra lines");
      break;
    }
    apply_prefix_filter(expected_line, time_prefix);
    apply_prefix_filter(got_line, time_prefix);
    if (strcmp(expected_line, got_line) != 0)
    {
      tmc_test_fail(__FILE__, __LINE__, "%s: line %d: expected \"%s\", got \"%s\"", name, line_number, expected_line, got_line);
      break;
    }
  }
  free(expected_text);
}

/* ---- helpers ----------------------------------------------------------------- */

static tmc_rate_constant arrhenius(double A, double B, double C)
{
  tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
  p.A = A;
  p.B = B;
  p.C = C;
  return tmc_rate_arrhenius(p);
}

/* Builds the solver, frees the builder, records a failure when the build fails */
static tmc_solver *build(tmc_builder *builder, const tmc_system *system, const tmc_reaction *reactions, int n)
{
  tmc_builder_set_system(builder, system);
  tmc_builder_set_reactions(builder, reactions, n);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_solver *solver = tmc_builder_build(builder, &err);
  if (!solver)
    tmc_test_fail(__FILE__, __LINE__, "build failed: %s", err.message);
  tmc_builder_free(builder);
  return solver;
}

/* MICM's CpuSolverBuilder<RosenbrockSolverParameters>(parameters): standard layout */
static tmc_builder *standard_builder(tmc_rosenbrock_parameters parameters)
{
  tmc_builder *builder = tmc_builder_create_rosenbrock(parameters);
  tmc_builder_set_width(builder, 1);
  return builder;
}

static void set_concentration(tmc_state *state, const char *name, double value)
{
  tmc_error err;
  tmc_error_clear(&err);
  TMC_EXPECT_EQ(tmc_state_set_concentration(state, name, value, &err), TMC_OK);
}

static void set_concentrations(tmc_state *state, const char *name, const double *values, int n)
{
  tmc_error err;
  tmc_error_clear(&err);
  TMC_EXPECT_EQ(tmc_state_set_concentrations(state, name, values, n, &err), TMC_OK);
}

static void set_custom(tmc_state *state, const char *label, double value)
{
  tmc_error err;
  tmc_error_clear(&err);
  if (tmc_state_set_custom_rate_parameter(state, label, value, &err) != TMC_OK)
    tmc_test_fail(__FILE__, __LINE__, "custom rate parameter %s: %s", label, err.message);
}

static void set_customs(tmc_state *state, const char *label, const double *values, int n)
{
  tmc_error err;
  tmc_error_clear(&err);
  if (tmc_state_set_custom_rate_parameters(state, label, values, n, &err) != TMC_OK)
    tmc_test_fail(__FILE__, __LINE__, "custom rate parameters %s: %s", label, err.message);
}

/* The A -> B, B + B -> B + C, B + C -> A + C mechanism with user-defined
 * rates r1, r2, r3 that five of the programs build */
typedef struct abc_mechanism
{
  tmc_phase gas;
  tmc_system system;
  tmc_reaction reactions[3];
} abc_mechanism;

static abc_mechanism abc_create(void)
{
  abc_mechanism m;
  tmc_species a = tmc_species_create("A"), b = tmc_species_create("B"), c = tmc_species_create("C");
  m.gas = tmc_phase_create("gas");
  tmc_phase_add_species(&m.gas, &a);
  tmc_phase_add_species(&m.gas, &b);
  tmc_phase_add_species(&m.gas, &c);
  m.reactions[0] = tmc_reaction_create();
  tmc_reaction_add_reactant(&m.reactions[0], &a);
  tmc_reaction_add_product(&m.reactions[0], &b, 1);
  tmc_reaction_set_rate_constant(&m.reactions[0], tmc_rate_user_defined("r1", 1.0));
  tmc_reaction_set_phase(&m.reactions[0], &m.gas);
  m.reactions[1] = tmc_reaction_create();
  tmc_reaction_add_reactant(&m.reactions[1], &b);
  tmc_reaction_add_reactant(&m.reactions[1], &b);
  tmc_reaction_add_product(&m.reactions[1], &b, 1);
  tmc_reaction_add_product(&m.reactions[1], &c, 1);
  tmc_reaction_set_rate_constant(&m.reactions[1], tmc_rate_user_defined("r2", 1.0));
  tmc_reaction_set_phase(&m.reactions[1], &m.gas);
  m.reactions[2] = tmc_reaction_create();
  tmc_reaction_add_reactant(&m.reactions[2], &b);
  tmc_reaction_add_reactant(&m.reactions[2], &c);
  tmc_reaction_add_product(&m.reactions[2], &a, 1);
  tmc_reaction_add_product(&m.reactions[2], &c, 1);
  tmc_reaction_set_rate_constant(&m.reactions[2], tmc_rate_user_defined("r3", 1.0));
  tmc_reaction_set_phase(&m.reactions[2], &m.gas);
  m.system = tmc_system_create(&m.gas);
  return m;
}

static void abc_free(abc_mechanism *m)
{
  for (int i = 0; i < 3; ++i)
    tmc_reaction_free(&m->reactions[i]);
  tmc_system_free(&m->system);
  tmc_phase_free(&m->gas);
}

static tmc_solver *abc_solver(abc_mechanism *m, tmc_rosenbrock_parameters parameters, int L)
{
  tmc_builder *builder = tmc_builder_create_rosenbrock(parameters);
  tmc_builder_set_width(builder, L);
  return build(builder, &m->system, m->reactions, 3);
}

/* MICM's SolverStateToString */
static const char *state_name(tmc_solver_state s)
{
  return tmc_solver_state_name(s);
}

/* ---- README_example ----------------------------------------------------------- */

TMC_TEST(readme_example)
{
  outbuf out = { NULL, 0, 0 };
  tmc_species foo = tmc_species_create("Foo"), bar = tmc_species_create("Bar"), baz = tmc_species_create("Baz");
  tmc_phase gas_phase = tmc_phase_create("gas");
  tmc_phase_add_species(&gas_phase, &foo);
  tmc_phase_add_species(&gas_phase, &bar);
  tmc_phase_add_species(&gas_phase, &baz);
  tmc_system chemical_system = tmc_system_create(&gas_phase);

  tmc_reaction r[2];
  r[0] = tmc_reaction_create();
  tmc_reaction_add_reactant(&r[0], &foo);
  tmc_reaction_add_product(&r[0], &bar, 0.8);
  tmc_reaction_add_product(&r[0], &baz, 0.2);
  tmc_reaction_set_rate_constant(&r[0], arrhenius(1.0e-3, 0, 0));
  tmc_reaction_set_phase(&r[0], &gas_phase);
  r[1] = tmc_reaction_create();
  tmc_reaction_add_reactant(&r[1], &foo);
  tmc_reaction_add_reactant(&r[1], &bar);
  tmc_reaction_add_product(&r[1], &baz, 1);
  tmc_reaction_set_rate_constant(&r[1], arrhenius(1.0e-5, 0, 110.0));
  tmc_reaction_set_phase(&r[1], &gas_phase);

  tmc_solver *solver = build(standard_builder(tmc_rosenbrock_three_stage()), &chemical_system, r, 2);
  tmc_reaction_free(&r[0]);
  tmc_reaction_free(&r[1]);
  tmc_system_free(&chemical_system);
  tmc_phase_free(&gas_phase);
  TMC_ASSERT_TRUE(solver != NULL);

  tmc_state *state = tmc_solver_create_state(solver, 1);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = 287.45; /* K */
  conditions[0].pressure = 101319.9;  /* Pa */
  tmc_conditions_calculate_ideal_air_density(&conditions[0]);
  set_concentration(state, "Foo", 20.0); /* mol m-3 */

  out_header(&out, state);
  for (int i = 0; i < 10; ++i)
  {
    tmc_solver_update_state_parameters(solver, state);
    tmc_result result = tmc_solver_solve(solver, 500.0, state);
    (void)result;
    out_state(&out, state, i * 500);
  }
  compare_output("README_example", &out, NULL, NULL);
  out_free(&out);
  tmc_state_free(state);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

/* ---- solve_results (test_solver_results.cpp) ---------------------------------- */

TMC_TEST(solve_results)
{
  outbuf out = { NULL, 0, 0 };
  abc_mechanism m = abc_create();
  const int number_of_grid_cells = 3;
  tmc_solver *solver = abc_solver(&m, tmc_rosenbrock_three_stage(), 1);
  abc_free(&m);
  TMC_ASSERT_TRUE(solver != NULL);
  tmc_state *state = tmc_solver_create_state(solver, number_of_grid_cells);

  /* mol m-3 */
  const double a0[3] = { 1, 2, 0.5 }, zeros[3] = { 0, 0, 0 };
  set_concentrations(state, "A", a0, 3);
  set_concentrations(state, "B", zeros, 3);
  set_concentrations(state, "C", zeros, 3);
  const double k1[3] = { 0.04, 0.04, 0.04 }, k2[3] = { 3e7, 3e7, 3e7 }, k3[3] = { 1e4, 1e4, 1e4 };
  set_customs(state, "r1", k1, 3);
  set_customs(state, "r2", k2, 3);
  set_customs(state, "r3", k3, 3);
  tmc_conditions *conditions = tmc_state_conditions(state);
  for (int cell = 0; cell < number_of_grid_cells; ++cell)
  {
    conditions[cell].temperature = 272.5; /* [K] */
    conditions[cell].pressure = 101253.3; /* [Pa] */
    conditions[cell].air_density = 1e6;   /* [mol m-3] */
  }
  tmc_solver_update_state_parameters(solver, state);

  const double time_step = 200; /* s */
  tmc_result result = tmc_solver_solve(solver, time_step, state);
  out_printf(&out, "Solver state: %s\n", state_name(result.state));
  out_printf(&out, "accepted: %ld\n", result.stats.accepted);
  out_printf(&out, "function_calls: %ld\n", result.stats.function_calls);
  out_printf(&out, "jacobian_updates: %ld\n", result.stats.jacobian_updates);
  out_printf(&out, "number_of_steps: %ld\n", result.stats.number_of_steps);
  out_printf(&out, "accepted: %ld\n", result.stats.accepted);
  out_printf(&out, "rejected: %ld\n", result.stats.rejected);
  out_printf(&out, "decompositions: %ld\n", result.stats.decompositions);
  out_printf(&out, "solves: %ld\n", result.stats.solves);
  out_printf(&out, "final simulation time: %g\n", result.stats.final_time);
  compare_output("solve_results", &out, NULL, NULL);
  out_free(&out);
  tmc_state_free(state);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

/* ---- multiple_grid_cells ------------------------------------------------------- */

TMC_TEST(multiple_grid_cells)
{
  outbuf out = { NULL, 0, 0 };
  abc_mechanism m = abc_create();
  const int number_of_grid_cells = 3;
  tmc_solver *solver = abc_solver(&m, tmc_rosenbrock_three_stage(), 1);
  abc_free(&m);
  TMC_ASSERT_TRUE(solver != NULL);
  tmc_state *state = tmc_solver_create_state(solver, number_of_grid_cells);

  const double a0[3] = { 1, 2, 0.5 }, zeros[3] = { 0, 0, 0 };
  set_concentrations(state, "A", a0, 3);
  set_concentrations(state, "B", zeros, 3);
  set_concentrations(state, "C", zeros, 3);
  const double k1[3] = { 0.04, 0.04, 0.04 }, k2[3] = { 3e7, 3e7, 3e7 }, k3[3] = { 1e4, 1e4, 1e4 };
  set_customs(state, "r1", k1, 3);
  set_customs(state, "r2", k2, 3);
  set_customs(state, "r3", k3, 3);
  tmc_conditions *conditions = tmc_state_conditions(state);
  for (int cell = 0; cell < number_of_grid_cells; ++cell)
  {
    conditions[cell].temperature = 272.5;
    conditions[cell].pressure = 101253.3;
    conditions[cell].air_density = 1e6;
  }

  const double time_step = 200; /* s */
  out_header(&out, state);
  out_state(&out, state, 0);
  for (int i = 0; i < 10; ++i)
  {
    double elapsed_solve_time = 0;
    while (elapsed_solve_time < time_step)
    {
      tmc_solver_update_state_parameters(solver, state);
      tmc_result result = tmc_solver_solve(solver, time_step - elapsed_solve_time, state);
      if (result.stats.final_time <= 0)
      {
        tmc_test_fail(__FILE__, __LINE__, "Solver made no progress: %s", state_name(result.state));
        break;
      }
      elapsed_solve_time += result.stats.final_time;
    }
    out_state(&out, state, time_step * (i + 1));
  }
  compare_output("multiple_grid_cells", &out, NULL, NULL);
  out_free(&out);
  tmc_state_free(state);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

/* ---- openmp_tutorial (test_openmp.cpp) ------------------------------------------ */

/* run_solver_on_thread_with_own_state: the "Running solver on thread" message
 * goes into the shared buffer, which the comparison skips */
static void run_solver_on_thread_with_own_state(outbuf *out, tmc_solver *solver, tmc_state *state, int thread, double *results)
{
#ifdef _OPENMP
#pragma omp critical(tutorial_output)
#endif
  out_printf(out, "Running solver on thread %d\n", thread);

  /* mol m-3 */
  static const double initial[3] = { 1, 0, 0 };
  tmc_dense_set_row(tmc_state_variables(state), 0, initial);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_state_set_custom_rate_parameter(state, "r1", 0.04, &err);
  tmc_state_set_custom_rate_parameter(state, "r2", 3e7, &err);
  tmc_state_set_custom_rate_parameter(state, "r3", 1e4, &err);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = 272.5; /* [K] */
  conditions[0].pressure = 101253.3; /* [Pa] */
  conditions[0].air_density = 1e6;   /* [mol m-3] */

  const double time_step = 200; /* s */
  for (int i = 0; i < 10; ++i)
  {
    double elapsed_solve_time = 0;
    tmc_solver_update_state_parameters(solver, state);
    while (elapsed_solve_time < time_step)
    {
      tmc_result result = tmc_solver_solve(solver, time_step - elapsed_solve_time, state);
      if (result.stats.final_time <= 0)
        break;
      /* MICM assigns here rather than accumulating */
      elapsed_solve_time = result.stats.final_time;
    }
  }
  tmc_dense_get_row(tmc_state_variables(state), 0, results);
}

TMC_TEST(openmp_tutorial)
{
  outbuf out = { NULL, 0, 0 };
  const int n_threads = 3;
  abc_mechanism m = abc_create();
  tmc_solver *solver = abc_solver(&m, tmc_rosenbrock_three_stage(), 1);
  abc_free(&m);
  TMC_ASSERT_TRUE(solver != NULL);

  double results[3][3];
  bool ran[3] = { false, false, false };
#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
  {
    /* each thread uses its own state */
    const int thread = omp_get_thread_num();
    if (thread < n_threads)
    {
      tmc_state *state = tmc_solver_create_state(solver, 1);
      run_solver_on_thread_with_own_state(&out, solver, state, thread, results[thread]);
      ran[thread] = true;
      tmc_state_free(state);
    }
  }
#endif
  /* without OpenMP, or when fewer threads were given, the remaining threads' work runs here */
  for (int thread = 0; thread < n_threads; ++thread)
    if (!ran[thread])
    {
      tmc_state *state = tmc_solver_create_state(solver, 1);
      run_solver_on_thread_with_own_state(&out, solver, state, thread, results[thread]);
      tmc_state_free(state);
    }

  for (int thread = 0; thread < n_threads; ++thread)
  {
    out_printf(&out, "Thread %d\n", thread + 1);
    out_printf(&out, "%10s,%10s,%10s\n", "A", "B", "C");
    out_printf(&out, "%10.2e,%10.2e,%10.2e\n", results[thread][0], results[thread][1], results[thread][2]);
  }
  compare_output("openmp_tutorial", &out, "Thread 1", NULL);
  out_free(&out);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

/* ---- rate_constants_*_example_by_hand ----------------------------------------- */

/* Conversion factor from moles m-3 to molecules cm-3 for consistency with the configuration file */
#define MOLES_M3_TO_MOLECULES_CM3 (1.0e-6 * 6.02214076e23)

/* The by-hand mechanism of the two programs: seven reactions, plus the three
 * user-defined ones of rate_constants_user_defined_example_by_hand */
static tmc_solver *by_hand_solver(bool user_defined)
{
  tmc_species a = tmc_species_create("A"), b = tmc_species_create("B");
  tmc_species c = tmc_species_create_with_property("C", TMC_PROPERTY_MOLECULAR_WEIGHT, 0.025);
  tmc_species d = tmc_species_create("D"), e = tmc_species_create("E"), f = tmc_species_create("F"), g = tmc_species_create("G");
  tmc_phase gas_phase = tmc_phase_create("gas");
  tmc_phase_add_species(&gas_phase, &a);
  tmc_phase_add_species(&gas_phase, &b);
  const double c_diffusion_coefficient = 2.3e2;
  tmc_phase_add_species_with_diffusion(&gas_phase, &c, c_diffusion_coefficient);
  tmc_phase_add_species(&gas_phase, &d);
  tmc_phase_add_species(&gas_phase, &e);
  tmc_phase_add_species(&gas_phase, &f);
  tmc_phase_add_species(&gas_phase, &g);

  tmc_reaction r[10];
  int n = 0;
  for (int i = 0; i < 10; ++i)
  {
    r[i] = tmc_reaction_create();
    tmc_reaction_set_phase(&r[i], &gas_phase);
  }
  /* r1: A -> B */
  tmc_reaction_add_reactant(&r[n], &a);
  tmc_reaction_add_product(&r[n], &b, 1);
  tmc_reaction_set_rate_constant(&r[n++], arrhenius(2.15e-4, 0, 110));
  /* a branched reaction has two output pathways, represented as two reactions */
  tmc_branched_parameters branched_params = { TMC_BRANCH_ALKOXY, 1.2, 204.3, 1.0e-3, 2 };
  tmc_reaction_add_reactant(&r[n], &b);
  tmc_reaction_add_product(&r[n], &c, 1);
  tmc_reaction_set_rate_constant(&r[n++], tmc_rate_branched(branched_params));
  branched_params.branch = TMC_BRANCH_NITRATE;
  tmc_reaction_add_reactant(&r[n], &b);
  tmc_reaction_add_product(&r[n], &d, 1);
  tmc_reaction_set_rate_constant(&r[n++], tmc_rate_branched(branched_params));
  /* r4: C -> E on a surface; the effective radius and particle number concentration are set later */
  tmc_reaction_add_reactant(&r[n], &c);
  tmc_reaction_add_product(&r[n], &e, 1);
  tmc_reaction_set_rate_constant(&r[n++], tmc_rate_surface("C", tmc_phase_find(&gas_phase, "C"), 0.90));
  /* r5: D -> 2 F */
  tmc_falloff_parameters ternary = tmc_falloff_defaults();
  ternary.k0_A = 1.2;
  ternary.k0_B = 2.3;
  ternary.k0_C = 302.3;
  ternary.kinf_A = 2.6 / MOLES_M3_TO_MOLECULES_CM3;
  ternary.kinf_B = -3.1;
  ternary.kinf_C = 402.1;
  ternary.Fc = 0.9;
  ternary.N = 1.2;
  tmc_reaction_add_reactant(&r[n], &d);
  tmc_reaction_add_product(&r[n], &f, 2);
  tmc_reaction_set_rate_constant(&r[n++], tmc_rate_ternary_chemical_activation(ternary));
  /* r6: E + E -> G (JPL 19-5's OH + NO2 + M falloff parameters) */
  tmc_falloff_parameters troe = tmc_falloff_defaults();
  troe.k0_A = 1.8e-30 * MOLES_M3_TO_MOLECULES_CM3 * MOLES_M3_TO_MOLECULES_CM3;
  troe.k0_B = -3.0;
  troe.k0_C = 0.0;
  troe.kinf_A = 2.8e-11 * MOLES_M3_TO_MOLECULES_CM3;
  troe.kinf_B = 0.0;
  troe.kinf_C = 0.0;
  troe.Fc = 0.6;
  troe.N = 1.0;
  tmc_reaction_add_reactant(&r[n], &e);
  tmc_reaction_add_reactant(&r[n], &e);
  tmc_reaction_add_product(&r[n], &g, 1);
  tmc_reaction_set_rate_constant(&r[n++], tmc_rate_troe(troe));
  /* r7: F -> G */
  tmc_tunneling_parameters tunneling = tmc_tunneling_defaults();
  tunneling.A = 1.2;
  tunneling.B = 2.3;
  tunneling.C = 302.3;
  tmc_reaction_add_reactant(&r[n], &f);
  tmc_reaction_add_product(&r[n], &g, 1);
  tmc_reaction_set_rate_constant(&r[n++], tmc_rate_tunneling(tunneling));
  if (user_defined)
  {
    /* r8: C -> G by photolysis */
    tmc_reaction_add_reactant(&r[n], &c);
    tmc_reaction_add_product(&r[n], &g, 1);
    tmc_reaction_set_rate_constant(&r[n++], tmc_rate_user_defined("my photolysis rate", 1.0));
    /* r9: emission of A */
    tmc_reaction_add_product(&r[n], &a, 1);
    tmc_reaction_set_rate_constant(&r[n++], tmc_rate_user_defined("my emission rate", 1.0));
    /* r10: first-order loss of B */
    tmc_reaction_add_reactant(&r[n], &b);
    tmc_reaction_set_rate_constant(&r[n++], tmc_rate_user_defined("my loss rate", 1.0));
  }

  tmc_system chemical_system = tmc_system_create(&gas_phase);
  tmc_solver *solver = build(standard_builder(tmc_rosenbrock_three_stage()), &chemical_system, r, n);
  for (int i = 0; i < 10; ++i)
    tmc_reaction_free(&r[i]);
  tmc_system_free(&chemical_system);
  tmc_phase_free(&gas_phase);
  return solver;
}

/* The same mechanism generated from configs/examples/<name> (generated/<name>.c) */
static tmc_solver *generated_solver(bool user_defined)
{
  tmc_system system = user_defined ? rate_constants_user_defined_system() : rate_constants_no_user_defined_system();
  tmc_reaction *reactions = NULL;
  const int n = user_defined ? rate_constants_user_defined_reactions(&reactions) : rate_constants_no_user_defined_reactions(&reactions);
  tmc_solver *solver = build(standard_builder(tmc_rosenbrock_three_stage()), &system, reactions, n);
  tmc_reactions_free(reactions, n);
  tmc_system_free(&system);
  return solver;
}

enum
{
  BH_STEPS = 10,
  BH_SPECIES = 7,
  BH_MAX_REACTIONS = 10
};
static const char *const bh_species[BH_SPECIES] = { "A", "B", "C", "D", "E", "F", "G" };

/* What a run of the program computes, for comparing two mechanisms */
typedef struct bh_record
{
  int n_reactions;
  char keys[BH_MAX_REACTIONS][TMC_NAME_MAX];          /* which reaction each rate constant column is */
  double rate_constants[BH_STEPS][BH_MAX_REACTIONS]; /* after each step's update */
  double concentrations[BH_STEPS + 1][BH_SPECIES];   /* by species A .. G, at each printed time */
} bh_record;

/* A name for a reaction of these mechanisms, the same in both builds */
static void reaction_key(const tmc_reaction *r, char *key, size_t size)
{
  const tmc_rate_constant *rc = &r->rate_constant;
  switch (rc->type)
  {
    case TMC_RATE_ARRHENIUS: snprintf(key, size, "arrhenius"); break;
    case TMC_RATE_TROE: snprintf(key, size, "troe"); break;
    case TMC_RATE_TERNARY_CHEMICAL_ACTIVATION: snprintf(key, size, "ternary chemical activation"); break;
    case TMC_RATE_BRANCHED: snprintf(key, size, "%s", rc->p.branched.branch == TMC_BRANCH_ALKOXY ? "branched alkoxy" : "branched nitrate"); break;
    case TMC_RATE_TUNNELING: snprintf(key, size, "tunneling"); break;
    case TMC_RATE_SURFACE: snprintf(key, size, "surface"); break;
    case TMC_RATE_USER_DEFINED: snprintf(key, size, "user-defined %s", rc->p.user_defined.label); break;
    default: snprintf(key, size, "type %d", (int)rc->type); break;
  }
}

/* The program's body after the solver is built. surface_label is the surface
 * reaction's label ("C" by hand, "C surface" in the configuration), which
 * names its two custom parameters. */
static void run_rate_constants_program(outbuf *out, tmc_solver *solver, bool user_defined, const char *surface_label, bh_record *record)
{
  memset(record, 0, sizeof *record);
  int n_reactions = 0;
  const tmc_reaction *reactions = tmc_solver_reactions(solver, &n_reactions);
  record->n_reactions = n_reactions <= BH_MAX_REACTIONS ? n_reactions : BH_MAX_REACTIONS;
  for (int k = 0; k < record->n_reactions; ++k)
    reaction_key(&reactions[k], record->keys[k], sizeof record->keys[k]);

  tmc_state *state = tmc_solver_create_state(solver, 1);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = 287.45; /* K */
  conditions[0].pressure = 101319.9;  /* Pa */
  tmc_conditions_calculate_ideal_air_density(&conditions[0]);
  for (int s = 0; s < BH_SPECIES; ++s)
    set_concentration(state, bh_species[s], s == 0 ? 1.0 : 0.0); /* mol m-3 */
  char label[TMC_NAME_MAX];
  snprintf(label, sizeof label, "%s.effective radius [m]", surface_label);
  set_custom(state, label, 1e-7);
  snprintf(label, sizeof label, "%s.particle number concentration [# m-3]", surface_label);
  set_custom(state, label, 2.5e6);

  const double time_step = 500; /* s */
  out_header(out, state);
  out_state(out, state, 0);
  for (int s = 0; s < BH_SPECIES; ++s)
    record->concentrations[0][s] = tmc_state_get(state, 0, tmc_state_variable_index(state, bh_species[s]));

  double photo_rate = 1e-10;
  if (user_defined)
  {
    const double emission_rate = 1e-20;
    const double loss = emission_rate * 1e-3;
    /* these rates are constant through the simulation */
    set_custom(state, "my emission rate", emission_rate);
    set_custom(state, "my loss rate", loss);
  }

  for (int i = 0; i < BH_STEPS; ++i)
  {
    double elapsed_solve_time = 0;
    if (user_defined)
      set_custom(state, "my photolysis rate", photo_rate); /* updated at each time step */
    tmc_solver_update_state_parameters(solver, state);
    tmc_dense_get_row(tmc_state_rate_constants(state), 0, record->rate_constants[i]);
    while (elapsed_solve_time < time_step)
    {
      tmc_result result = tmc_solver_solve(solver, time_step - elapsed_solve_time, state);
      if (result.stats.final_time <= 0)
      {
        tmc_test_fail(__FILE__, __LINE__, "Solver made no progress: %s", state_name(result.state));
        break;
      }
      elapsed_solve_time += result.stats.final_time;
    }
    out_state(out, state, time_step * (i + 1));
    for (int s = 0; s < BH_SPECIES; ++s)
      record->concentrations[i + 1][s] = tmc_state_get(state, 0, tmc_state_variable_index(state, bh_species[s]));
    if (user_defined)
      photo_rate *= 1.5;
  }
  tmc_state_free(state);
}

/* The generated mechanism gives the same rate constants and concentrations as
 * the one built in code, to 1e-10 relative */
static void compare_records(const bh_record *by_hand, const bh_record *generated)
{
  TMC_EXPECT_EQ(by_hand->n_reactions, generated->n_reactions);
  for (int k = 0; k < by_hand->n_reactions; ++k)
  {
    int match = -1;
    for (int j = 0; j < generated->n_reactions && match < 0; ++j)
      if (strcmp(by_hand->keys[k], generated->keys[j]) == 0)
        match = j;
    if (match < 0)
    {
      tmc_test_fail(__FILE__, __LINE__, "the generated mechanism has no reaction \"%s\"", by_hand->keys[k]);
      continue;
    }
    for (int i = 0; i < BH_STEPS; ++i)
    {
      const double a = by_hand->rate_constants[i][k], b = generated->rate_constants[i][match];
      if (!(fabs(a - b) <= 1e-10 * fmax(fabs(a), fabs(b))))
        tmc_test_fail(__FILE__, __LINE__, "step %d, %s: rate constant %.17g by hand, %.17g generated", i, by_hand->keys[k], a, b);
    }
  }
  for (int i = 0; i <= BH_STEPS; ++i)
    for (int s = 0; s < BH_SPECIES; ++s)
    {
      const double a = by_hand->concentrations[i][s], b = generated->concentrations[i][s];
      if (!(fabs(a - b) <= 1e-10 * fmax(fabs(a), fabs(b))))
        tmc_test_fail(__FILE__, __LINE__, "time %g, %s: %.17g by hand, %.17g generated", 500.0 * i, bh_species[s], a, b);
    }
}

static int rate_constants_example(bool user_defined, const char *name)
{
  outbuf out = { NULL, 0, 0 }, generated_out = { NULL, 0, 0 };
  static bh_record by_hand_record, generated_record;
  tmc_solver *solver = by_hand_solver(user_defined);
  TMC_ASSERT_TRUE(solver != NULL);
  run_rate_constants_program(&out, solver, user_defined, "C", &by_hand_record);
  compare_output(name, &out, NULL, NULL);
  tmc_solver_free(solver);

  tmc_solver *generated = generated_solver(user_defined);
  TMC_ASSERT_TRUE(generated != NULL);
  run_rate_constants_program(&generated_out, generated, user_defined, "C surface", &generated_record);
  compare_output(name, &generated_out, NULL, NULL);
  compare_records(&by_hand_record, &generated_record);
  tmc_solver_free(generated);
  out_free(&out);
  out_free(&generated_out);
  return TMC_TEST_PASS;
}

TMC_TEST(rate_constants_no_user_defined_example_by_hand)
{
  return rate_constants_example(false, "rate_constants_no_user_defined_example_by_hand");
}

TMC_TEST(rate_constants_user_defined_example_by_hand)
{
  return rate_constants_example(true, "rate_constants_user_defined_example_by_hand");
}

/* ---- solver_configuration ------------------------------------------------------ */

static double nanoseconds_now(void)
{
  struct timespec ts;
  if (timespec_get(&ts, TIME_UTC) != TIME_UTC)
    return 0.0;
  return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}

static void test_solver_type(outbuf *out, tmc_solver *solver)
{
  tmc_state *state = tmc_solver_create_state(solver, 1);
  static const double initial[3] = { 1, 0, 0 }; /* mol m-3 */
  tmc_dense_set_row(tmc_state_variables(state), 0, initial);
  set_custom(state, "r1", 0.04);
  set_custom(state, "r2", 3e7);
  set_custom(state, "r3", 1e4);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = 272.5; /* [K] */
  conditions[0].pressure = 101253.3; /* [Pa] */
  conditions[0].air_density = 1e6;   /* [mol m-3] */

  const double time_step = 200; /* s */
  out_header(out, state);
  out_state(out, state, 0);
  tmc_stats total_stats;
  memset(&total_stats, 0, sizeof total_stats);
  double total_solve_time = 0;
  for (int i = 0; i < 10; ++i)
  {
    double elapsed_solve_time = 0;
    while (elapsed_solve_time < time_step)
    {
      const double start = nanoseconds_now();
      tmc_solver_update_state_parameters(solver, state);
      tmc_result result = tmc_solver_solve(solver, time_step - elapsed_solve_time, state);
      total_solve_time += nanoseconds_now() - start;
      total_stats.function_calls += result.stats.function_calls;
      total_stats.jacobian_updates += result.stats.jacobian_updates;
      total_stats.number_of_steps += result.stats.number_of_steps;
      total_stats.accepted += result.stats.accepted;
      total_stats.rejected += result.stats.rejected;
      total_stats.decompositions += result.stats.decompositions;
      total_stats.solves += result.stats.solves;
      if (result.stats.final_time <= 0)
      {
        tmc_test_fail(__FILE__, __LINE__, "Solver made no progress: %s", state_name(result.state));
        break;
      }
      elapsed_solve_time += result.stats.final_time;
    }
    out_state(out, state, time_step * (i + 1));
  }
  out_printf(out, "Total solve time: %g nanoseconds\n", total_solve_time);
  out_printf(out, "accepted: %ld\n", total_stats.accepted);
  out_printf(out, "function_calls: %ld\n", total_stats.function_calls);
  out_printf(out, "jacobian_updates: %ld\n", total_stats.jacobian_updates);
  out_printf(out, "number_of_steps: %ld\n", total_stats.number_of_steps);
  out_printf(out, "accepted: %ld\n", total_stats.accepted);
  out_printf(out, "rejected: %ld\n", total_stats.rejected);
  out_printf(out, "decompositions: %ld\n", total_stats.decompositions);
  out_printf(out, "solves: %ld\n", total_stats.solves);
  tmc_state_free(state);
}

TMC_TEST(solver_configuration)
{
  outbuf out = { NULL, 0, 0 };
  abc_mechanism m = abc_create();
  tmc_solver *two_stage = abc_solver(&m, tmc_rosenbrock_two_stage(), 1);
  tmc_solver *three_stage = abc_solver(&m, tmc_rosenbrock_three_stage(), 1);
  tmc_solver *four_stage = abc_solver(&m, tmc_rosenbrock_four_stage(), 1);
  tmc_solver *four_stage_da = abc_solver(&m, tmc_rosenbrock_four_stage_dae(), 1);
  tmc_solver *six_stage_da = abc_solver(&m, tmc_rosenbrock_six_stage_dae(), 1);
  abc_free(&m);
  TMC_ASSERT_TRUE(two_stage && three_stage && four_stage && four_stage_da && six_stage_da);

  out_printf(&out, "Two stages: \n");
  test_solver_type(&out, two_stage);
  out_printf(&out, "\nThree stages: \n");
  test_solver_type(&out, three_stage);
  out_printf(&out, "\nFour stages: \n");
  test_solver_type(&out, four_stage);
  out_printf(&out, "\nFour stages differential algebraic: \n");
  test_solver_type(&out, four_stage_da);
  out_printf(&out, "\nSix stages differential algebraic: \n");
  test_solver_type(&out, six_stage_da);
  compare_output("solver_configuration", &out, NULL, "Total solve time:");
  out_free(&out);
  tmc_solver_free(two_stage);
  tmc_solver_free(three_stage);
  tmc_solver_free(four_stage);
  tmc_solver_free(four_stage_da);
  tmc_solver_free(six_stage_da);
  return TMC_TEST_PASS;
}

/* ---- vectorized_matrix_solver -------------------------------------------------- */

/* The program's Solve: its inner loop never runs, because solver_state starts
 * as Converged, so the states keep their initial values; kept as written. */
static void vms_solve(tmc_solver *solver, tmc_state *state, int number_of_grid_cells)
{
  const double k1[3] = { 0.04, 0.04, 0.04 }, k2[3] = { 3e7, 3e7, 3e7 }, k3[3] = { 1e4, 1e4, 1e4 };
  set_customs(state, "r1", k1, 3);
  set_customs(state, "r2", k2, 3);
  set_customs(state, "r3", k3, 3);
  tmc_conditions *conditions = tmc_state_conditions(state);
  for (int cell = 0; cell < number_of_grid_cells; ++cell)
  {
    conditions[cell].temperature = 272.5; /* [K] */
    conditions[cell].pressure = 101253.3; /* [Pa] */
    conditions[cell].air_density = 1e6;   /* [mol m-3] */
  }
  const double time_step = 100; /* s */
  tmc_solver_state solver_state = TMC_CONVERGED;
  for (int i = 0; i < 10 && solver_state == TMC_CONVERGED; ++i)
  {
    double elapsed_solve_time = 0;
    tmc_solver_update_state_parameters(solver, state);
    while (elapsed_solve_time < time_step && solver_state != TMC_CONVERGED)
    {
      tmc_result result = tmc_solver_solve(solver, time_step - elapsed_solve_time, state);
      elapsed_solve_time += result.stats.final_time;
      solver_state = result.state;
    }
  }
  if (solver_state != TMC_CONVERGED)
    tmc_test_fail(__FILE__, __LINE__, "Solver did not converge");
}

/* The strings "block.row.column" of a 4 x 4 sparse matrix with elements (0, 1),
 * (2, 1), (2, 3), (3, 2) and 3 blocks, in the order of its storage (MICM's
 * AsVector of a SparseMatrix<std::string>) */
static void print_sparse_storage(outbuf *out, int L)
{
  tmc_sparse_builder *b = tmc_sparse_builder_create(4);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_sparse_builder_with_element(b, 0, 1, &err);
  tmc_sparse_builder_with_element(b, 2, 1, &err);
  tmc_sparse_builder_with_element(b, 2, 3, &err);
  tmc_sparse_builder_with_element(b, 3, 2, &err);
  tmc_sparse_builder_set_number_of_blocks(b, 3);
  tmc_sparse *m = tmc_sparse_build(b, TMC_SPARSE_CSR, L, false);
  tmc_sparse_builder_free(b);
  char (*labels)[16] = calloc(m->size, sizeof *labels);
  if (!labels)
    abort();
  for (int block = 0; block < m->n_blocks; ++block)
    for (int e = 0; e < m->nnz; ++e)
    {
      const size_t index = tmc_sparse_data_index(m, block, e);
      snprintf(labels[index], sizeof labels[index], "%d.%d.%d", block, m->row_of[e], m->col_of[e]);
    }
  for (size_t k = 0; k < m->size; ++k)
    out_printf(out, "%s\n", labels[k]);
  free(labels);
  tmc_sparse_free(m);
}

TMC_TEST(vectorized_matrix_solver)
{
  outbuf out = { NULL, 0, 0 };
  abc_mechanism m = abc_create();
  const int number_of_grid_cells = 3;
  tmc_solver *solver = abc_solver(&m, tmc_rosenbrock_three_stage(), 1);
  TMC_ASSERT_TRUE(solver != NULL);
  tmc_state *state = tmc_solver_create_state(solver, number_of_grid_cells);
  const double a0[3] = { 1.1, 2.1, 3.1 }, b0[3] = { 1.2, 2.2, 3.2 }, c0[3] = { 1.3, 2.3, 3.3 };
  set_concentrations(state, "A", a0, 3);
  set_concentrations(state, "B", b0, 3);
  set_concentrations(state, "C", c0, 3);
  const tmc_dense *variables = tmc_state_variables(state);
  for (size_t k = 0; k < variables->size; ++k)
    out_printf(&out, "%g\n", variables->data[k]);
  out_printf(&out, "\n");

  /* VectorMatrix<double, 3> and SparseMatrixVectorOrdering<3> */
  tmc_solver *vectorized_solver = abc_solver(&m, tmc_rosenbrock_three_stage(), 3);
  TMC_ASSERT_TRUE(vectorized_solver != NULL);
  tmc_state *vectorized_state = tmc_solver_create_state(vectorized_solver, number_of_grid_cells);
  set_concentrations(vectorized_state, "A", a0, 3);
  set_concentrations(vectorized_state, "B", b0, 3);
  set_concentrations(vectorized_state, "C", c0, 3);
  const tmc_dense *vectorized_variables = tmc_state_variables(vectorized_state);
  for (size_t k = 0; k < vectorized_variables->size; ++k)
    out_printf(&out, "%g\n", vectorized_variables->data[k]);
  out_printf(&out, "\n");

  vms_solve(solver, state, number_of_grid_cells);
  vms_solve(vectorized_solver, vectorized_state, number_of_grid_cells);

  char names[3][TMC_NAME_MAX];
  const int n_names = tmc_system_unique_names(&m.system, names, 3);
  for (int cell = 0; cell < number_of_grid_cells; ++cell)
  {
    out_printf(&out, "Cell %d\n", cell);
    out_printf(&out, "%10s%20s%20s\n", "Species", "Regular", "Vectorized");
    for (int s = 0; s < n_names; ++s)
    {
      const int regular_idx = tmc_state_variable_index(state, names[s]);
      const int vectorized_idx = tmc_state_variable_index(vectorized_state, names[s]);
      out_printf(&out, "%10s%20g%20g\n", names[s], tmc_state_get(state, cell, regular_idx), tmc_state_get(vectorized_state, cell, vectorized_idx));
    }
    out_printf(&out, "\n");
  }

  out_printf(&out, "\nSparse matrix standard ordering elements\n");
  print_sparse_storage(&out, 1);
  out_printf(&out, "\nSparse matrix vector ordering elements\n");
  print_sparse_storage(&out, 3);

  compare_output("vectorized_matrix_solver", &out, NULL, NULL);
  out_free(&out);
  tmc_state_free(state);
  tmc_state_free(vectorized_state);
  tmc_solver_free(solver);
  tmc_solver_free(vectorized_solver);
  abc_free(&m);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("README_example", readme_example),
  TMC_CASE("solve_results", solve_results),
  TMC_CASE("multiple_grid_cells", multiple_grid_cells),
  TMC_CASE("rate_constants_no_user_defined_example_by_hand", rate_constants_no_user_defined_example_by_hand),
  TMC_CASE("rate_constants_user_defined_example_by_hand", rate_constants_user_defined_example_by_hand),
  TMC_CASE("solver_configuration", solver_configuration),
  TMC_CASE("vectorized_matrix_solver", vectorized_matrix_solver),
  TMC_CASE("openmp_tutorial", openmp_tutorial),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_tutorials(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
