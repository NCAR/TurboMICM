/* MICM test/integration/analytical_policy.hpp and analytical_surface_rxn_policy.hpp,
 * see analytical_policy.h. The mechanisms, initial states, tolerances, and
 * closed-form solutions are MICM's; the solver runs through the C API. */
#include "analytical_policy.h"

enum
{
  NSTEPS = 1000,
  NUM_CELLS = 3
};

/* MICM's constants.hpp */
#define BOLTZMANN_CONSTANT 1.380649e-23
#define AVOGADRO_CONSTANT 6.02214076e23
#define GAS_CONSTANT (BOLTZMANN_CONSTANT * AVOGADRO_CONSTANT)
#define PI 3.14159265358979323846

static const char *current; /* the configuration under test, for messages */

#define FAIL(...) tmc_test_fail(__FILE__, __LINE__, __VA_ARGS__)

/* EXPECT_NEAR(actual, reference, tol) */
static void expect_near(double actual, double reference, double tol, const char *what, int step, int cell, int species)
{
  if (!(fabs(actual - reference) <= tol))
    FAIL("%s [%s]: step %d cell %d species %d: %.17g vs %.17g, tolerance %g", what, current, step, cell, species, actual, reference, tol);
}

/* MICM's EXPECT_REAL_REL in a double build: |actual - reference| <= rel |reference| */
static void expect_rel(double actual, double reference, double rel, const char *what, int step, int cell, int species)
{
  expect_near(actual, reference, rel * fabs(reference), what, step, cell, species);
}

/* MICM's RelativeError: |a - b| / |b| */
static double relative_error(double a, double b)
{
  return fabs(a - b) / fabs(b);
}

static double air_density_mol_m3(double pressure, double temperature)
{
  return pressure / (GAS_CONSTANT * temperature);
}

/* ---- configurations ----------------------------------------------------- */

static tmc_rosenbrock_parameters rosenbrock_parameters(int stages)
{
  switch (stages)
  {
    case 2: return tmc_rosenbrock_two_stage();
    case 4: return tmc_rosenbrock_four_stage();
    case 44: return tmc_rosenbrock_four_stage_dae();
    case 6: return tmc_rosenbrock_six_stage_dae();
    default: return tmc_rosenbrock_three_stage();
  }
}

tmc_builder *ap_create_builder(const ap_config *c)
{
  tmc_builder *b = c->method == TMC_METHOD_BACKWARD_EULER ? tmc_builder_create_backward_euler(tmc_backward_euler_defaults())
                                                          : tmc_builder_create_rosenbrock(rosenbrock_parameters(c->stages));
  tmc_builder_set_width(b, c->L);
  tmc_builder_set_lu_algorithm(b, c->lu);
  return b;
}

/* ---- mechanism helpers -------------------------------------------------- */

static tmc_phase make_gas_phase(const char *const *names, int n)
{
  tmc_phase gas = tmc_phase_create("gas");
  for (int i = 0; i < n; ++i)
  {
    tmc_species s = tmc_species_create(names[i]);
    tmc_phase_add_species(&gas, &s);
  }
  return gas;
}

static const tmc_species *species_of(const tmc_phase *gas, const char *name)
{
  return &tmc_phase_find(gas, name)->species;
}

/* A reaction with no reactants or products yet, in the gas phase */
static tmc_reaction reaction(const tmc_phase *gas, tmc_rate_constant rc)
{
  tmc_reaction r = tmc_reaction_create();
  tmc_reaction_set_phase(&r, gas);
  tmc_reaction_set_rate_constant(&r, rc);
  return r;
}

/* reactant -> product, yield 1 */
static tmc_reaction unimolecular(const tmc_phase *gas, const char *reactant, const char *product, tmc_rate_constant rc)
{
  tmc_reaction r = reaction(gas, rc);
  tmc_reaction_add_reactant(&r, species_of(gas, reactant));
  tmc_reaction_add_product(&r, species_of(gas, product), 1.0);
  return r;
}

static tmc_rate_constant arrhenius(double A, double B, double C, double D, double E)
{
  tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
  p.A = A;
  p.B = B;
  p.C = C;
  p.D = D;
  p.E = E;
  return tmc_rate_arrhenius(p);
}

static tmc_rate_constant arrhenius_A(double A)
{
  tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
  p.A = A;
  return tmc_rate_arrhenius(p);
}

static tmc_falloff_parameters falloff(double k0_A, double k0_B, double k0_C, double kinf_A, double kinf_B, double kinf_C, double Fc, double N)
{
  tmc_falloff_parameters p = tmc_falloff_defaults();
  p.k0_A = k0_A;
  p.k0_B = k0_B;
  p.k0_C = k0_C;
  p.kinf_A = kinf_A;
  p.kinf_B = kinf_B;
  p.kinf_C = kinf_C;
  p.Fc = Fc;
  p.N = N;
  return p;
}

/* The second reaction's falloff parameters of the Troe and ternary tests */
static tmc_falloff_parameters falloff_k2(void)
{
  return falloff(1.2e-3, 1.6, 3, 136, 5, 24, 0.9, 0.8);
}

static tmc_rate_constant tunneling(double A, double B, double C)
{
  tmc_tunneling_parameters p = tmc_tunneling_defaults();
  p.A = A;
  p.B = B;
  p.C = C;
  return tmc_rate_tunneling(p);
}

static tmc_rate_constant branched(tmc_branch branch)
{
  tmc_branched_parameters p;
  p.branch = branch;
  p.X = 1e-4;
  p.Y = 204.3;
  p.a0 = 1.0e-3;
  p.n = 2;
  return tmc_rate_branched(p);
}

static tmc_rate_constant user_defined(const char *label)
{
  return tmc_rate_user_defined(label, 1.0);
}

/* builder.SetSystem(System(gas)).SetReactions(reactions).Build(): frees the
 * builder, the reactions, and the phase; returns NULL after recording a failure. */
static tmc_solver *build(tmc_builder *builder, tmc_phase *gas, tmc_reaction *reactions, int n_reactions)
{
  tmc_system system = tmc_system_create(gas);
  tmc_builder_set_system(builder, &system);
  tmc_builder_set_reactions(builder, reactions, n_reactions);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_solver *solver = tmc_builder_build(builder, &err);
  if (!solver)
    FAIL("[%s] build failed: %s", current, err.message);
  tmc_builder_free(builder);
  for (int i = 0; i < n_reactions; ++i)
    tmc_reaction_free(&reactions[i]);
  tmc_system_free(&system);
  tmc_phase_free(gas);
  return solver;
}

static int variable(const tmc_state *state, const char *name)
{
  int i = tmc_state_variable_index(state, name);
  if (i < 0)
    FAIL("[%s] no state variable %s", current, name);
  return i;
}

static void set_custom(tmc_state *state, const char *label, double value)
{
  tmc_error err;
  tmc_error_clear(&err);
  if (tmc_state_set_custom_rate_parameter(state, label, value, &err) != TMC_OK)
    FAIL("[%s] custom rate parameter %s: %s", current, label, err.message);
}

static void set_concentration(tmc_state *state, const char *name, double value)
{
  tmc_error err;
  tmc_error_clear(&err);
  if (tmc_state_set_concentration(state, name, value, &err) != TMC_OK)
    FAIL("[%s] concentration %s: %s", current, name, err.message);
}

static void set_absolute_tolerances(tmc_state *state, const double *atol, int n)
{
  tmc_error err;
  tmc_error_clear(&err);
  if (tmc_state_set_absolute_tolerances(state, atol, n, &err) != TMC_OK)
    FAIL("[%s] absolute tolerances: %s", current, err.message);
}

/* Solves delta_t in as many calls as the solver needs (MICM's
 * `while (actual_solve < delta_t)` loops); a call that advances no time is a
 * failure here, where MICM would spin. */
static void solve_through(tmc_solver *solver, tmc_state *state, double delta_t, const char *what)
{
  double actual_solve = 0;
  while (actual_solve < delta_t)
  {
    tmc_result result = tmc_solver_solve(solver, delta_t - actual_solve, state);
    if (!(result.stats.final_time > 0))
    {
      FAIL("%s [%s]: solve advanced no time (%s)", what, current, tmc_solver_state_name(result.state));
      return;
    }
    actual_solve += result.stats.final_time;
  }
}

/* ---- A -> B -> C against its closed form (TestSimpleSystem and TestSimpleStiffSystem) ---- */

typedef double (*rate_fn)(double temperature, double pressure, double air_density);

typedef struct custom_parameter
{
  const char *label;
  double value;
} custom_parameter;

/* The stiff variant splits A into A1 and A2, kept in fast equilibrium, and
 * compares A1 + A2 with A; its rate constants are not checked. */
static void run_simple_system(const char *label, tmc_solver *solver, double abs_tol, rate_fn calculate_k1, rate_fn calculate_k2, const custom_parameter *custom, int n_custom, bool stiff)
{
  if (!solver)
    return;
  static const double temperatures[3] = { 272.5, 254.7, 312.6 };
  static const double pressures[3] = { 101253.3, 100672.5, 101319.8 };
  double k1[NUM_CELLS], k2[NUM_CELLS];
  for (int i = 0; i < NUM_CELLS; ++i)
  {
    const double temperature = temperatures[i % 3], pressure = pressures[i % 3];
    const double air_density = air_density_mol_m3(pressure, temperature);
    k1[i] = calculate_k1(temperature, pressure, air_density);
    k2[i] = calculate_k2(temperature, pressure, air_density);
  }

  const double time_step = 1.0;
  tmc_state *state = tmc_solver_create_state(solver, NUM_CELLS);
  tmc_error err;
  tmc_error_clear(&err);
  for (int k = 0; k < n_custom; ++k)
  {
    double values[NUM_CELLS];
    for (int i = 0; i < NUM_CELLS; ++i)
      values[i] = custom[k].value;
    if (tmc_state_set_custom_rate_parameters(state, custom[k].label, values, NUM_CELLS, &err) != TMC_OK)
      FAIL("%s [%s]: custom rate parameter %s: %s", label, current, custom[k].label, err.message);
  }

  const int a1 = variable(state, stiff ? "A1" : "A"), a2 = stiff ? variable(state, "A2") : -1;
  const int b = variable(state, "B"), c = variable(state, "C");
  double initial_A[NUM_CELLS];
  tmc_conditions *conditions = tmc_state_conditions(state);
  for (int i = 0; i < NUM_CELLS; ++i)
  {
    initial_A[i] = 1.0 - (double)i / (double)NUM_CELLS;
    if (stiff)
    {
      tmc_state_set(state, i, a1, 0.5 * initial_A[i]);
      tmc_state_set(state, i, a2, 0.5 * initial_A[i]);
    }
    else
      tmc_state_set(state, i, a1, initial_A[i]);
    tmc_state_set(state, i, b, 0.0);
    tmc_state_set(state, i, c, 0.0);
    conditions[i].temperature = temperatures[i % 3];
    conditions[i].pressure = pressures[i % 3];
    conditions[i].air_density = air_density_mol_m3(pressures[i % 3], temperatures[i % 3]);
  }

  const tmc_dense *rate_constants = tmc_state_rate_constants(state);
  for (int i_time = 1; i_time < NSTEPS; ++i_time)
  {
    tmc_solver_update_state_parameters(solver, state);
    tmc_result result = tmc_solver_solve(solver, time_step, state);
    if (result.state != TMC_CONVERGED)
      FAIL("%s [%s]: step %d: solver state %s", label, current, i_time, tmc_solver_state_name(result.state));
    const double t = i_time * time_step;
    for (int i = 0; i < NUM_CELLS; ++i)
    {
      if (!stiff)
      {
        /* MICM compares the rate constants relatively, to 1e-9 */
        expect_rel(k1[i], tmc_dense_get(rate_constants, i, 0), 1e-9, label, i_time, i, 0);
        expect_rel(k2[i], tmc_dense_get(rate_constants, i, 1), 1e-9, label, i_time, i, 1);
      }
      const double model_A = stiff ? tmc_state_get(state, i, a1) + tmc_state_get(state, i, a2) : tmc_state_get(state, i, a1);
      const double model_B = tmc_state_get(state, i, b), model_C = tmc_state_get(state, i, c);
      const double kk1 = k1[i], kk2 = k2[i];
      const double analytical_A = initial_A[i] * exp(-kk1 * t);
      const double analytical_B = initial_A[i] * (kk1 / (kk2 - kk1)) * (exp(-kk1 * t) - exp(-kk2 * t));
      const double analytical_C = initial_A[i] * (1.0 + (kk1 * exp(-kk2 * t) - kk2 * exp(-kk1 * t)) / (kk2 - kk1));
      expect_near(model_A, analytical_A, abs_tol, label, i_time, i, 0);
      expect_near(model_B, analytical_B, abs_tol, label, i_time, i, 1);
      expect_near(model_C, analytical_C, abs_tol, label, i_time, i, 2);
    }
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
}

/* A -> B (k1), B -> C (k2) */
static tmc_solver *simple_solver(const ap_config *c, tmc_rate_constant k1, tmc_rate_constant k2)
{
  static const char *const names[3] = { "A", "B", "C" };
  tmc_phase gas = make_gas_phase(names, 3);
  tmc_reaction r[2];
  r[0] = unimolecular(&gas, "A", "B", k1);
  r[1] = unimolecular(&gas, "B", "C", k2);
  return build(ap_create_builder(c), &gas, r, 2);
}

/* A1 -> B (k1), A2 -> B (k1), B -> C (k2), A1 <-> A2 fast */
static tmc_solver *stiff_solver(const ap_config *c, tmc_rate_constant k1a, tmc_rate_constant k1b, tmc_rate_constant k2)
{
  static const char *const names[4] = { "A1", "A2", "B", "C" };
  tmc_phase gas = make_gas_phase(names, 4);
  tmc_reaction r[5];
  r[0] = unimolecular(&gas, "A1", "B", k1a);
  r[1] = unimolecular(&gas, "A2", "B", k1b);
  r[2] = unimolecular(&gas, "B", "C", k2);
  r[3] = unimolecular(&gas, "A1", "A2", arrhenius_A(4.0e10));
  r[4] = unimolecular(&gas, "A2", "A1", arrhenius_A(0.9 * 4.0e10));
  return build(ap_create_builder(c), &gas, r, 5);
}

/* ---- Troe ----------------------------------------------------------------- */

static double troe_k1(double temperature, double pressure, double air_density)
{
  (void)temperature;
  (void)pressure;
  const double k_0 = 4.0e-11, k_inf = 1;
  return k_0 * air_density / (1.0 + k_0 * air_density / k_inf) * pow(0.6, 1.0 / (1.0 + pow(log10(k_0 * air_density / k_inf), 2)));
}

static double troe_k2(double temperature, double pressure, double air_density)
{
  (void)pressure;
  const double k_0 = 1.2e-3 * exp(3.0 / temperature) * pow(temperature / 300.0, 1.6);
  const double k_inf = 136.0 * exp(24.0 / temperature) * pow(temperature / 300.0, 5.0);
  return k_0 * air_density / (1.0 + k_0 * air_density / k_inf) * pow(0.9, 0.8 / (0.8 + pow(log10(k_0 * air_density / k_inf), 2)));
}

/* the stiff test writes the exponent as 1 / (1 + (1 / 0.8) log10^2) */
static double stiff_troe_k2(double temperature, double pressure, double air_density)
{
  (void)pressure;
  const double k_0 = 1.2e-3 * exp(3.0 / temperature) * pow(temperature / 300.0, 1.6);
  const double k_inf = 136.0 * exp(24.0 / temperature) * pow(temperature / 300.0, 5.0);
  return k_0 * air_density / (1.0 + k_0 * air_density / k_inf) * pow(0.9, 1.0 / (1.0 + (1.0 / 0.8) * pow(log10(k_0 * air_density / k_inf), 2)));
}

void ap_troe(const ap_config *c, double abs_tol)
{
  current = c->name;
  tmc_falloff_parameters k1 = tmc_falloff_defaults();
  k1.k0_A = 4.0e-11;
  tmc_solver *solver = simple_solver(c, tmc_rate_troe(k1), tmc_rate_troe(falloff_k2()));
  run_simple_system("troe", solver, abs_tol, troe_k1, troe_k2, NULL, 0, false);
}

void ap_stiff_troe(const ap_config *c, double abs_tol)
{
  current = c->name;
  tmc_falloff_parameters k1 = tmc_falloff_defaults();
  k1.k0_A = 4.0e-11;
  tmc_solver *solver = stiff_solver(c, tmc_rate_troe(k1), tmc_rate_troe(k1), tmc_rate_troe(falloff_k2()));
  run_simple_system("troe", solver, abs_tol, troe_k1, stiff_troe_k2, NULL, 0, true);
}

/* ---- Photolysis ----------------------------------------------------------- */

static double photolysis_k1(double temperature, double pressure, double air_density)
{
  (void)temperature;
  (void)pressure;
  (void)air_density;
  return 2e-3;
}

static double photolysis_k2(double temperature, double pressure, double air_density)
{
  (void)temperature;
  (void)pressure;
  (void)air_density;
  return 3e-3;
}

void ap_photolysis(const ap_config *c, double abs_tol)
{
  current = c->name;
  static const custom_parameter custom[2] = { { "photoA", 2e-3 }, { "photoB", 3e-3 } };
  tmc_solver *solver = simple_solver(c, user_defined("photoA"), user_defined("photoB"));
  run_simple_system("photolysis", solver, abs_tol, photolysis_k1, photolysis_k2, custom, 2, false);
}

void ap_stiff_photolysis(const ap_config *c, double abs_tol)
{
  current = c->name;
  static const custom_parameter custom[3] = { { "photoA1B", 2e-3 }, { "photoA2B", 2e-3 }, { "photoB", 3e-3 } };
  tmc_solver *solver = stiff_solver(c, user_defined("photoA1B"), user_defined("photoA2B"), user_defined("photoB"));
  run_simple_system("photolysis", solver, abs_tol, photolysis_k1, photolysis_k2, custom, 3, true);
}

/* ---- Ternary chemical activation ------------------------------------------ */

static double ternary_k1(double temperature, double pressure, double air_density)
{
  (void)temperature;
  (void)pressure;
  const double k_0 = 4.0e-5, k_inf = 1;
  return k_0 / (1.0 + k_0 * air_density / k_inf) * pow(0.6, 1.0 / (1.0 + (1.0 / 1.0) * pow(log10(k_0 * air_density / k_inf), 2)));
}

static double ternary_k2(double temperature, double pressure, double air_density)
{
  (void)pressure;
  const double k_0 = 1.2e-3 * exp(3.0 / temperature) * pow(temperature / 300.0, 1.6);
  const double k_inf = 136.0 * exp(24.0 / temperature) * pow(temperature / 300.0, 5.0);
  return k_0 / (1.0 + k_0 * air_density / k_inf) * pow(0.9, 1.0 / (1.0 + (1.0 / 0.8) * pow(log10(k_0 * air_density / k_inf), 2)));
}

void ap_ternary_chemical_activation(const ap_config *c, double abs_tol)
{
  current = c->name;
  tmc_falloff_parameters k1 = tmc_falloff_defaults();
  k1.k0_A = 4.0e-5;
  k1.kinf_A = 1;
  tmc_solver *solver = simple_solver(c, tmc_rate_ternary_chemical_activation(k1), tmc_rate_ternary_chemical_activation(falloff_k2()));
  run_simple_system("ternary_chemical_activation", solver, abs_tol, ternary_k1, ternary_k2, NULL, 0, false);
}

void ap_stiff_ternary_chemical_activation(const ap_config *c, double abs_tol)
{
  current = c->name;
  tmc_falloff_parameters k1 = tmc_falloff_defaults();
  k1.k0_A = 4.0e-5;
  k1.kinf_A = 1;
  tmc_solver *solver = stiff_solver(c, tmc_rate_ternary_chemical_activation(k1), tmc_rate_ternary_chemical_activation(k1), tmc_rate_ternary_chemical_activation(falloff_k2()));
  run_simple_system("ternary_chemical_activation", solver, abs_tol, ternary_k1, ternary_k2, NULL, 0, true);
}

/* ---- Tunneling ------------------------------------------------------------ */

static double tunneling_k1(double temperature, double pressure, double air_density)
{
  (void)temperature;
  (void)pressure;
  (void)air_density;
  return 4.0e-3;
}

static double tunneling_k2(double temperature, double pressure, double air_density)
{
  (void)pressure;
  (void)air_density;
  return 1.2e-4 * exp(-167 / temperature + 1.0e8 / pow(temperature, 3));
}

void ap_tunneling(const ap_config *c, double abs_tol)
{
  current = c->name;
  tmc_solver *solver = simple_solver(c, tunneling(4.0e-3, 0, 0), tunneling(1.2e-4, 167, 1.0e8));
  run_simple_system("tunneling", solver, abs_tol, tunneling_k1, tunneling_k2, NULL, 0, false);
}

void ap_stiff_tunneling(const ap_config *c, double abs_tol)
{
  current = c->name;
  tmc_solver *solver = stiff_solver(c, tunneling(4.0e-3, 0, 0), tunneling(4.0e-3, 0, 0), tunneling(1.2e-4, 167, 1.0e8));
  run_simple_system("tunneling", solver, abs_tol, tunneling_k1, tunneling_k2, NULL, 0, true);
}

/* ---- Arrhenius ------------------------------------------------------------ */

static double arrhenius_k1(double temperature, double pressure, double air_density)
{
  (void)pressure;
  (void)air_density;
  return 4.0e-3 * exp(50 / temperature);
}

static double arrhenius_k2(double temperature, double pressure, double air_density)
{
  (void)air_density;
  return 1.2e-4 * exp(75 / temperature) * pow(temperature / 50, 7) * (1.0 + 0.5 * pressure);
}

/* the stiff test's B -> C has B = 1.6 */
static double stiff_arrhenius_k2(double temperature, double pressure, double air_density)
{
  (void)air_density;
  return 1.2e-4 * exp(75 / temperature) * pow(temperature / 50, 1.6) * (1.0 + 0.5 * pressure);
}

void ap_arrhenius(const ap_config *c, double abs_tol)
{
  current = c->name;
  tmc_solver *solver = simple_solver(c, arrhenius(4.0e-3, 0, 50, 300, 0), arrhenius(1.2e-4, 7, 75, 50, 0.5));
  run_simple_system("arrhenius", solver, abs_tol, arrhenius_k1, arrhenius_k2, NULL, 0, false);
}

void ap_stiff_arrhenius(const ap_config *c, double abs_tol)
{
  current = c->name;
  tmc_solver *solver = stiff_solver(c, arrhenius(4.0e-3, 0, 50, 300, 0), arrhenius(4.0e-3, 0, 50, 300, 0), arrhenius(1.2e-4, 1.6, 75, 50, 0.5));
  run_simple_system("arrhenius", solver, abs_tol, arrhenius_k1, stiff_arrhenius_k2, NULL, 0, true);
}

/* ---- Branched ------------------------------------------------------------- */

/* MICM's Wennberg branching: z from the reference conditions, A from the cell's */
static void branched_z_A(double temperature, double air_density, double *z_out, double *A_out)
{
  const double air_dens_n_cm3 = air_density * AVOGADRO_CONSTANT * 1.0e-6;
  double a_ = 2.0e-22 * exp(2) * 2.45e19;
  double b_ = 0.43 * pow(293.0 / 298.0, -8.0);
  const double z = a_ / (1.0 + a_ / b_) * pow(0.41, 1.0 / (1.0 + pow(log10(a_ / b_), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
  a_ = 2.0e-22 * exp(2) * air_dens_n_cm3;
  b_ = 0.43 * pow(temperature / 298.0, -8.0);
  const double A = a_ / (1.0 + a_ / b_) * pow(0.41, 1.0 / (1.0 + pow(log10(a_ / b_), 2)));
  *z_out = z;
  *A_out = A;
}

static double branched_k1(double temperature, double pressure, double air_density)
{
  (void)pressure;
  double z, A;
  branched_z_A(temperature, air_density, &z, &A);
  return 1e-4 * exp(-204.3 / temperature) * (z / (z + A));
}

static double branched_k2(double temperature, double pressure, double air_density)
{
  (void)pressure;
  double z, A;
  branched_z_A(temperature, air_density, &z, &A);
  return 1e-4 * exp(-204.3 / temperature) * (A / (z + A));
}

void ap_branched(const ap_config *c, double abs_tol)
{
  current = c->name;
  tmc_solver *solver = simple_solver(c, branched(TMC_BRANCH_ALKOXY), branched(TMC_BRANCH_NITRATE));
  run_simple_system("branched", solver, abs_tol, branched_k1, branched_k2, NULL, 0, false);
}

void ap_stiff_branched(const ap_config *c, double abs_tol)
{
  current = c->name;
  tmc_solver *solver = stiff_solver(c, branched(TMC_BRANCH_ALKOXY), branched(TMC_BRANCH_ALKOXY), branched(TMC_BRANCH_NITRATE));
  run_simple_system("branched", solver, abs_tol, branched_k1, branched_k2, NULL, 0, true);
}

/* ---- Surface reaction (analytical_surface_rxn_policy.hpp) ----------------- */

void ap_surface_rxn(const ap_config *c, double tolerance)
{
  current = c->name;
  /* parameters, from CAMP/test/unit_rxn_data/test_rxn_surface.F90 */
  const double mode_GMD = 1.0e-6;           /* mode geometric mean diameter [m] */
  const double mode_GSD = 0.1;              /* mode geometric standard deviation [unitless] */
  const double DENSITY_stuff = 1000.0;      /* [kg m-3] */
  const double DENSITY_more_stuff = 1000.0; /* [kg m-3] */
  const double MW_foo = 0.04607;            /* [kg mol-1] */
  const double Dg_foo = 0.95e-5;            /* diffusion coefficient [m2 s-1] */
  const double rxn_gamma = 2.0e-2;          /* [unitless] */
  const double bar_yield = 1.0;             /* [unitless] */
  const double baz_yield = 0.4;             /* [unitless] */
  const double temperature = 272.5;         /* K */
  const double pressure = 101253.3;         /* Pa */
  const double conc_foo = 1.0;
  const double conc_stuff = 2.0e-3;
  const double conc_more_stuff = 3.0e-3;

  const double radius = mode_GMD / 2.0 * exp(5.0 * log(mode_GSD) * log(mode_GSD) / 2.0);
  const double number_conc = 6.0 / (PI * pow(mode_GMD, 3.0) * exp(9.0 / 2.0 * log(mode_GSD) * log(mode_GSD))) * (conc_stuff / DENSITY_stuff + conc_more_stuff / DENSITY_more_stuff);

  tmc_species foo = tmc_species_create_with_property("foo", TMC_PROPERTY_MOLECULAR_WEIGHT, MW_foo);
  tmc_species bar = tmc_species_create("bar"), baz = tmc_species_create("baz");
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species_with_diffusion(&gas, &foo, Dg_foo);
  tmc_phase_add_species(&gas, &bar);
  tmc_phase_add_species(&gas, &baz);

  tmc_reaction r = reaction(&gas, tmc_rate_surface("foo", tmc_phase_find(&gas, "foo"), rxn_gamma));
  tmc_reaction_add_reactant(&r, &foo);
  tmc_reaction_add_product(&r, &bar, bar_yield);
  tmc_reaction_add_product(&r, &baz, baz_yield);
  tmc_solver *solver = build(ap_create_builder(c), &gas, &r, 1);
  if (!solver)
    return;

  tmc_state *state = tmc_solver_create_state(solver, 1);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = temperature;
  conditions[0].pressure = pressure;
  set_custom(state, "foo.effective radius [m]", radius);
  set_custom(state, "foo.particle number concentration [# m-3]", number_conc);
  set_concentration(state, "foo", conc_foo);

  const double mean_free_speed = sqrt(8.0 * GAS_CONSTANT / (PI * MW_foo) * temperature);
  const double k1 = 4.0 * number_conc * PI * radius * radius / (radius / Dg_foo + 4.0 / (mean_free_speed * rxn_gamma));
  const double time_step = 0.1 / k1;
  const int nstep = 10;
  const int idx_foo = variable(state, "foo"), idx_bar = variable(state, "bar"), idx_baz = variable(state, "baz");

  for (int i = 1; i <= nstep; ++i)
  {
    tmc_solver_update_state_parameters(solver, state);
    tmc_result result = tmc_solver_solve(solver, time_step, state);
    if (result.state != TMC_CONVERGED)
      FAIL("surface [%s]: step %d: solver state %s", current, i, tmc_solver_state_name(result.state));
    expect_near(k1, tmc_dense_get(tmc_state_rate_constants(state), 0, 0), 1e-8, "surface rate constant", i, 0, 0);
    const double t = i * time_step;
    const double curr_foo = conc_foo * exp(-k1 * t);
    const double curr_bar = bar_yield * (1.0 - curr_foo);
    const double curr_baz = baz_yield * (1.0 - curr_foo);
    expect_rel(tmc_state_get(state, 0, idx_foo), curr_foo, tolerance, "surface", i, 0, 0);
    expect_rel(tmc_state_get(state, 0, idx_bar), curr_bar, tolerance, "surface", i, 0, 1);
    expect_rel(tmc_state_get(state, 0, idx_baz), curr_baz, tolerance, "surface", i, 0, 2);
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
}

/* ---- Robertson ------------------------------------------------------------ */

void ap_robertson(const ap_config *c, double relative_tolerance)
{
  current = c->name;
  static const char *const names[3] = { "A", "B", "C" };
  tmc_phase gas = make_gas_phase(names, 3);
  tmc_reaction r[3];
  r[0] = unimolecular(&gas, "A", "B", user_defined("r1"));
  r[1] = reaction(&gas, user_defined("r2")); /* B + B -> B + C */
  tmc_reaction_add_reactant(&r[1], species_of(&gas, "B"));
  tmc_reaction_add_reactant(&r[1], species_of(&gas, "B"));
  tmc_reaction_add_product(&r[1], species_of(&gas, "B"), 1);
  tmc_reaction_add_product(&r[1], species_of(&gas, "C"), 1);
  r[2] = reaction(&gas, user_defined("r3")); /* B + C -> A + C */
  tmc_reaction_add_reactant(&r[2], species_of(&gas, "B"));
  tmc_reaction_add_reactant(&r[2], species_of(&gas, "C"));
  tmc_reaction_add_product(&r[2], species_of(&gas, "A"), 1);
  tmc_reaction_add_product(&r[2], species_of(&gas, "C"), 1);
  tmc_builder *builder = ap_create_builder(c);
  tmc_builder_set_reorder_state(builder, false);
  tmc_solver *solver = build(builder, &gas, r, 3);
  if (!solver)
    return;

  tmc_state *state = tmc_solver_create_state(solver, 1);
  tmc_state_set_relative_tolerance(state, 1e-10);
  const double atol[3] = { 1e-10 * 1e-2, 1e-10 * 1e-2, 1e-10 * 1e-2 };
  set_absolute_tolerances(state, atol, 3);
  set_custom(state, "r1", 0.04);
  set_custom(state, "r2", 3e7);
  set_custom(state, "r3", 1e4);

  enum
  {
    N = 12
  };
  static const double analytical[N + 1][3] = {
    { 1, 0, 0 },
    { 0.9664597373330035E+00, 0.3074626578578675E-04, 0.3350951640121071E-01 },
    { 0.8413699238414729E+00, 0.1623390937990473E-04, 0.1586138422491472E+00 },
    { 0.6172348823960878E+00, 0.6153591274639123E-05, 0.3827589640126376E+00 },
    { 0.3368745306607069E+00, 0.2013702318261393E-05, 0.6631234556369748E+00 },
    { 0.1073004285378040E+00, 0.4800166972571660E-06, 0.8926990914454987E+00 },
    { 0.1786592114209946E-01, 0.7274751468436319E-07, 0.9821340061103859E+00 },
    { 0.2031483924973415E-02, 0.8142277783356159E-08, 0.9979685079327488E+00 },
    { 0.2076093439016395E-03, 0.8306077485067610E-09, 0.9997923898254906E+00 },
    { 0.2082417512179460E-04, 0.8329841429908955E-10, 0.9999791757415798E+00 },
    { 0.2083229471647004E-05, 0.8332935037760723E-11, 0.9999979167621954E+00 },
    { 0.2083328471883087E-06, 0.8333315602809495E-12, 0.9999997916663195E+00 },
    { 0.2083340149701284E-07, 0.8333360770334744E-13, 0.9999999791665152E+00 }
  };
  static const double initial[3] = { 1, 0, 0 };
  tmc_dense_set_row(tmc_state_variables(state), 0, initial);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = 272.5;
  conditions[0].pressure = 101253.3;
  conditions[0].air_density = 1e6;

  double model[N + 1][3];
  memcpy(model[0], initial, sizeof initial);
  double target_time = 1.0, current_time = 0.0;
  tmc_solver_update_state_parameters(solver, state);
  for (int i_time = 0; i_time < N; ++i_time)
  {
    solve_through(solver, state, target_time - current_time, "robertson");
    tmc_dense_get_row(tmc_state_variables(state), 0, model[i_time + 1]);
    current_time = target_time;
    target_time *= 10;
  }

  const int idx[3] = { variable(state, "A"), variable(state, "B"), variable(state, "C") };
  const double absolute_tolerance = 0.3e-3; /* average of the starting concentration * 1e-3 */
  for (int i = 1; i <= N; ++i)
    for (int j = 0; j < 3; ++j)
    {
      const double rel_error = relative_error(model[i][idx[j]], analytical[i][j]);
      const double abs_error = fabs(model[i][idx[j]] - analytical[i][j]);
      if (!(abs_error < absolute_tolerance || rel_error < relative_tolerance))
        FAIL("robertson [%s]: output %d species %d: relative error %g, absolute error %g", current, i, j, rel_error, abs_error);
    }
  tmc_state_free(state);
  tmc_solver_free(solver);
}

/* ---- Oregonator ----------------------------------------------------------- */

void ap_oregonator(const ap_config *c, double relative_tolerance, int substeps_per_output, double solver_relative_tolerance)
{
  current = c->name;
  static const char *const names[5] = { "X", "Y", "Z", "P", "Q" };
  tmc_phase gas = make_gas_phase(names, 5);
  tmc_reaction r[5];
  r[0] = unimolecular(&gas, "Y", "X", user_defined("r1"));
  r[1] = reaction(&gas, user_defined("r2")); /* X + Y -> P */
  tmc_reaction_add_reactant(&r[1], species_of(&gas, "X"));
  tmc_reaction_add_reactant(&r[1], species_of(&gas, "Y"));
  tmc_reaction_add_product(&r[1], species_of(&gas, "P"), 1);
  r[2] = reaction(&gas, user_defined("r3")); /* X -> Z + 2 X */
  tmc_reaction_add_reactant(&r[2], species_of(&gas, "X"));
  tmc_reaction_add_product(&r[2], species_of(&gas, "Z"), 1);
  tmc_reaction_add_product(&r[2], species_of(&gas, "X"), 2);
  r[3] = reaction(&gas, user_defined("r4")); /* 2 X -> Q */
  tmc_reaction_add_reactant(&r[3], species_of(&gas, "X"));
  tmc_reaction_add_reactant(&r[3], species_of(&gas, "X"));
  tmc_reaction_add_product(&r[3], species_of(&gas, "Q"), 1);
  r[4] = unimolecular(&gas, "Z", "Y", user_defined("r5"));
  tmc_builder *builder = ap_create_builder(c);
  tmc_builder_set_reorder_state(builder, false);
  tmc_solver *solver = build(builder, &gas, r, 5);
  if (!solver)
    return;

  /* The parameters of equation.f in the Hairer test set */
  const double s_const = 77.27, q_const = 8.375e-6, w_const = 0.161;
  /* The two rate constants that s, q, and w leave free, at the Field and Noyes values */
  const double k2_const = 1.6e9, k5_const = 1.0;
  /* The three rate constants that s, q, and w fix */
  const double k1a_const = k5_const / (s_const * w_const);
  const double k3b_const = s_const * k5_const / w_const;
  const double k4_const = s_const * s_const * q_const * k2_const / 2.0;
  /* The scale of each dimensionless variable, and of the dimensionless time */
  const double alpha_const = k5_const / (s_const * w_const * k2_const);
  const double eta_const = s_const * k5_const / (w_const * k2_const);
  const double rho_const = k5_const * k5_const / (w_const * w_const * k2_const);
  const double tau = w_const / k5_const;

  const double time_step = 30 * tau;
  enum
  {
    N = 12
  };
  double analytical[N + 1][5] = {
    { 1, 2, 3, 0, 0 },
    { 0.1000661467180497E+01, 0.1512778937348249E+04, 0.1035854312767229E+05, 0, 0 },
    { 0.1000874625199626E+01, 0.1144336972384497E+04, 0.8372149966624639E+02, 0, 0 },
    { 0.1001890368438751E+01, 0.5299926232295553E+03, 0.1662279579042420E+01, 0, 0 },
    { 0.1004118022612645E+01, 0.2438326079910346E+03, 0.1008822224048647E+01, 0, 0 },
    { 0.1008995416634061E+01, 0.1121664388662539E+03, 0.1007783229065319E+01, 0, 0 },
    { 0.1019763472537298E+01, 0.5159761322947535E+02, 0.1016985778956374E+01, 0, 0 },
    { 0.1043985088527474E+01, 0.2373442027531524E+02, 0.1037691843544522E+01, 0, 0 },
    { 0.1100849071667922E+01, 0.1091533805469020E+02, 0.1085831969810860E+01, 0, 0 },
    { 0.1249102130020572E+01, 0.5013945178605446E+01, 0.1208326626237875E+01, 0, 0 },
    { 0.1779724751937019E+01, 0.2281852385542403E+01, 0.1613754023671725E+01, 0, 0 },
    { 0.1000889326903503E+01, 0.1125438585746596E+04, 0.1641049483777168E+05, 0, 0 },
    { 0.1000814870318523E+01, 0.1228178521549889E+04, 0.1320554942846513E+03, 0, 0 },
  };
  for (int i = 0; i <= N; ++i)
  {
    analytical[i][0] *= alpha_const;
    analytical[i][1] *= eta_const;
    analytical[i][2] *= rho_const;
  }
  const double initial[5] = { 1 * alpha_const, 2 * eta_const, 3 * rho_const, 0, 0 };

  tmc_state *state = tmc_solver_create_state(solver, 1);
  const double achievable_relative_tolerance = fmax(solver_relative_tolerance, 100.0 * 2.220446049250313e-16);
  tmc_state_set_relative_tolerance(state, achievable_relative_tolerance);
  const double tolerance_floor = achievable_relative_tolerance * 1e-2;
  const double atol[5] = { alpha_const * tolerance_floor, eta_const * tolerance_floor, rho_const * tolerance_floor, eta_const, eta_const };
  set_absolute_tolerances(state, atol, 5);
  set_custom(state, "r1", k1a_const);
  set_custom(state, "r2", k2_const);
  set_custom(state, "r3", k3b_const);
  set_custom(state, "r4", k4_const);
  set_custom(state, "r5", k5_const);
  tmc_dense_set_row(tmc_state_variables(state), 0, initial);
  tmc_solver_update_state_parameters(solver, state);

  double model[N + 1][5];
  memcpy(model[0], initial, sizeof initial);
  bool ok = true;
  for (int i_time = 0; i_time < N && ok; ++i_time)
  {
    /* A solver that controls its own step size takes the full output interval
     * in one call; a first-order solver takes sub-steps. */
    double actual_solve = 0;
    const double max_substep = time_step / substeps_per_output;
    while (actual_solve < time_step)
    {
      const double dt = fmin(max_substep, time_step - actual_solve);
      tmc_result result = tmc_solver_solve(solver, dt, state);
      if (!(result.state == TMC_CONVERGED || result.state == TMC_CONVERGENCE_EXCEEDED_MAX_STEPS))
      {
        FAIL("oregonator [%s]: solver returned %s at t = %g", current, tmc_solver_state_name(result.state), i_time * time_step + actual_solve);
        ok = false;
        break;
      }
      const double advanced_to = actual_solve + result.stats.final_time;
      if (!(advanced_to > actual_solve))
      {
        FAIL("oregonator [%s]: solve advanced only %g s at t = %g", current, result.stats.final_time, i_time * time_step + actual_solve);
        ok = false;
        break;
      }
      actual_solve = advanced_to;
    }
    tmc_dense_get_row(tmc_state_variables(state), 0, model[i_time + 1]);
  }

  if (ok)
  {
    const int idx[3] = { variable(state, "X"), variable(state, "Y"), variable(state, "Z") };
    for (int i = 1; i <= N; ++i)
      for (int j = 0; j < 3; ++j)
      {
        const double rel_error = relative_error(model[i][idx[j]], analytical[i][j]);
        if (!(rel_error < relative_tolerance))
          FAIL("oregonator [%s]: output %d species %s: relative error %g (%.17g vs %.17g)", current, i, names[j], rel_error, model[i][idx[j]], analytical[i][j]);
      }
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
}

/* ---- HIRES ---------------------------------------------------------------- */

void ap_hires(const ap_config *c, double absolute_tolerance)
{
  current = c->name;
  static const char *const names[8] = { "Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7" };
  tmc_phase gas = make_gas_phase(names, 8);
  tmc_reaction r[9];
  r[0] = unimolecular(&gas, "Y0", "Y1", user_defined("r1"));
  r[1] = reaction(&gas, user_defined("r2")); /* y1 -> (0.43 / 8.75) y0 + (8.32 / 8.75) y3 */
  tmc_reaction_add_reactant(&r[1], species_of(&gas, "Y1"));
  tmc_reaction_add_product(&r[1], species_of(&gas, "Y0"), 0.43 / 8.75);
  tmc_reaction_add_product(&r[1], species_of(&gas, "Y3"), 8.32 / 8.75);
  r[2] = reaction(&gas, user_defined("r3")); /* y2 -> (8.32 / 10.03) y0 + (1.71 / 10.03) y3 */
  tmc_reaction_add_reactant(&r[2], species_of(&gas, "Y2"));
  tmc_reaction_add_product(&r[2], species_of(&gas, "Y0"), 8.32 / 10.03);
  tmc_reaction_add_product(&r[2], species_of(&gas, "Y3"), 1.71 / 10.03);
  r[3] = reaction(&gas, user_defined("r4")); /* -> y0 */
  tmc_reaction_add_product(&r[3], species_of(&gas, "Y0"), 1);
  r[4] = reaction(&gas, user_defined("r5")); /* y3 -> (0.43 / 1.12) y2 + (0.69 / 1.12) y5 */
  tmc_reaction_add_reactant(&r[4], species_of(&gas, "Y3"));
  tmc_reaction_add_product(&r[4], species_of(&gas, "Y2"), 0.43 / 1.12);
  tmc_reaction_add_product(&r[4], species_of(&gas, "Y5"), 0.69 / 1.12);
  r[5] = reaction(&gas, user_defined("r6")); /* y4 -> (0.035 / 1.745) y2 + (1.71 / 1.745) y5 */
  tmc_reaction_add_reactant(&r[5], species_of(&gas, "Y4"));
  tmc_reaction_add_product(&r[5], species_of(&gas, "Y2"), 0.035 / 1.745);
  tmc_reaction_add_product(&r[5], species_of(&gas, "Y5"), 1.71 / 1.745);
  r[6] = unimolecular(&gas, "Y5", "Y4", user_defined("r7"));
  r[7] = reaction(&gas, user_defined("r8")); /* y6 -> (0.43 / 1.81) y4 + (0.69 / 1.81) y5 + y7 */
  tmc_reaction_add_reactant(&r[7], species_of(&gas, "Y6"));
  tmc_reaction_add_product(&r[7], species_of(&gas, "Y4"), 0.43 / 1.81);
  tmc_reaction_add_product(&r[7], species_of(&gas, "Y5"), 0.69 / 1.81);
  tmc_reaction_add_product(&r[7], species_of(&gas, "Y7"), 1);
  r[8] = reaction(&gas, user_defined("r9")); /* y5 + y7 -> y6 */
  tmc_reaction_add_reactant(&r[8], species_of(&gas, "Y5"));
  tmc_reaction_add_reactant(&r[8], species_of(&gas, "Y7"));
  tmc_reaction_add_product(&r[8], species_of(&gas, "Y6"), 1);
  tmc_builder *builder = ap_create_builder(c);
  tmc_builder_set_reorder_state(builder, false);
  tmc_solver *solver = build(builder, &gas, r, 9);
  if (!solver)
    return;

  enum
  {
    N = 2
  };
  static const double analytical[N + 1][8] = {
    { 1, 0, 0, 0, 0, 0, 0, 0.0057 },
    { 0.000737131257332567, 0.000144248572631618, 0.000058887297409676, 0.001175651343283149, 0.002386356198831330, 0.006238968252742796, 0.002849998395185769, 0.002850001604814231 },
    { 0.000670305503581864, 0.000130996846986347, 0.000046862231597733, 0.001044668020551705, 0.000594883830951485, 0.001399628833942774, 0.001014492757718480, 0.004685507242281520 },
  };
  tmc_state *state = tmc_solver_create_state(solver, 1);
  tmc_state_set_relative_tolerance(state, 1e-6);
  double atol[8];
  for (int i = 0; i < 8; ++i)
    atol[i] = 1e-6 * 1e-2;
  set_absolute_tolerances(state, atol, 8);
  static const double rates[9] = { 1.71, 8.75, 10.03, 0.0007, 1.12, 1.745, 0.43, 1.81, 280.0 };
  for (int i = 0; i < 9; ++i)
  {
    char label[8];
    snprintf(label, sizeof label, "r%d", i + 1);
    set_custom(state, label, rates[i]);
  }
  tmc_dense_set_row(tmc_state_variables(state), 0, analytical[0]);
  tmc_solver_update_state_parameters(solver, state);

  double model[N + 1][8];
  memcpy(model[0], analytical[0], sizeof model[0]);
  double time_step = 321.8122;
  for (int i_time = 0; i_time < N; ++i_time)
  {
    solve_through(solver, state, time_step, "hires");
    tmc_dense_get_row(tmc_state_variables(state), 0, model[i_time + 1]);
    time_step += 100;
  }
  for (int i = 1; i <= N; ++i)
    for (int j = 0; j < 8; ++j)
      expect_near(model[i][variable(state, names[j])], analytical[i][j], absolute_tolerance, "hires", i, 0, j);
  tmc_state_free(state);
  tmc_solver_free(solver);
}

/* ---- E5 ------------------------------------------------------------------- */

void ap_e5(const ap_config *c, double relative_tolerance)
{
  current = c->name;
  static const char *const names[6] = { "A1", "A2", "A3", "A4", "A5", "A6" };
  tmc_phase gas = make_gas_phase(names, 6);
  tmc_reaction r[4];
  r[0] = reaction(&gas, user_defined("r1")); /* A1 -> A2 + A3 */
  tmc_reaction_add_reactant(&r[0], species_of(&gas, "A1"));
  tmc_reaction_add_product(&r[0], species_of(&gas, "A2"), 1);
  tmc_reaction_add_product(&r[0], species_of(&gas, "A3"), 1);
  r[1] = reaction(&gas, user_defined("r2")); /* A2 + A3 -> A5 */
  tmc_reaction_add_reactant(&r[1], species_of(&gas, "A2"));
  tmc_reaction_add_reactant(&r[1], species_of(&gas, "A3"));
  tmc_reaction_add_product(&r[1], species_of(&gas, "A5"), 1);
  r[2] = reaction(&gas, user_defined("r3")); /* A1 + A3 -> A4 */
  tmc_reaction_add_reactant(&r[2], species_of(&gas, "A1"));
  tmc_reaction_add_reactant(&r[2], species_of(&gas, "A3"));
  tmc_reaction_add_product(&r[2], species_of(&gas, "A4"), 1);
  r[3] = reaction(&gas, user_defined("r4")); /* A4 -> A3 + A6 */
  tmc_reaction_add_reactant(&r[3], species_of(&gas, "A4"));
  tmc_reaction_add_product(&r[3], species_of(&gas, "A3"), 1);
  tmc_reaction_add_product(&r[3], species_of(&gas, "A6"), 1);
  tmc_builder *builder = ap_create_builder(c);
  tmc_builder_set_reorder_state(builder, false);
  tmc_solver *solver = build(builder, &gas, r, 4);
  if (!solver)
    return;

  enum
  {
    N = 7
  };
  /* the concentrations of A5 and A6 are ignored */
  static const double analytical[N + 1][6] = {
    { 1.76e-3, 0, 0, 0, 0, 0 },
    { 1.7599259497677897058e-003, 1.3846281519376516449e-011, 7.6370038530073911180e-013, 1.3082581134075777338e-011, 0, 0 },
    { 1.6180769999072942552e-003, 1.3822370304983735443e-010, 8.2515735006838336088e-012, 1.2997212954915352082e-010, 0, 0 },
    { 7.4813208224292220114e-006, 2.3734781561205975019e-012, 2.2123586689581663654e-012, 1.6111948716243113653e-013, 0, 0 },
    { 4.7150333630401632232e-010, 1.8188895860807021729e-014, 1.8188812376786725407e-014, 8.3484020296321693074e-020, 0, 0 },
    { 3.1317148329356996037e-014, 1.4840957952870064294e-016, 1.4840957948345691466e-016, 4.5243728279782625194e-026, 0, 0 },
    { 3.8139035189787091771e-049, 1.0192582567660293322e-020, 1.0192582567660293322e-020, 3.7844935507486221171e-065, 0, 0 },
    { 0.0000000000000000000e-000, 8.8612334976263783420e-023, 8.8612334976263783421e-023, 0.0000000000000000000e-000, 0, 0 }
  };
  tmc_state *state = tmc_solver_create_state(solver, 1);
  tmc_state_set_relative_tolerance(state, 1e-13);
  double atol[6];
  for (int i = 0; i < 6; ++i)
    atol[i] = 1e-17;
  set_absolute_tolerances(state, atol, 6);
  memcpy(atol, tmc_state_absolute_tolerances(state), sizeof atol);
  atol[0] = 1e-7;
  atol[4] = 1e-7;
  atol[5] = 1e-7;
  set_absolute_tolerances(state, atol, 6);
  set_custom(state, "r1", 7.89e-10);
  set_custom(state, "r2", 1.13e9);
  set_custom(state, "r3", 1.1e7);
  set_custom(state, "r4", 1.13e3);
  tmc_dense_set_row(tmc_state_variables(state), 0, analytical[0]);
  tmc_solver_update_state_parameters(solver, state);

  double model[N + 1][6];
  memcpy(model[0], analytical[0], sizeof model[0]);
  double target_time = 10.0, current_time = 0.0;
  for (int i_time = 0; i_time < N; ++i_time)
  {
    solve_through(solver, state, target_time - current_time, "e5");
    tmc_dense_get_row(tmc_state_variables(state), 0, model[i_time + 1]);
    current_time = target_time;
    target_time *= 100;
  }
  static const double absolute_tolerances[4] = { 1e-6, 1e-13, 1e-13, 1e-13 };
  for (int i = 0; i <= N; ++i)
    for (int j = 0; j < 4; ++j)
    {
      const double rel_error = relative_error(model[i][j], analytical[i][j]);
      const double abs_error = fabs(model[i][j] - analytical[i][j]);
      if (!(abs_error < absolute_tolerances[j] || rel_error < relative_tolerance))
        FAIL("e5 [%s]: output %d species %d: relative error %g, absolute error %g", current, i, j, rel_error, abs_error);
    }
  tmc_state_free(state);
  tmc_solver_free(solver);
}
