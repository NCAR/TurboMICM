#include "kinetics.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static const chem_species three_species[] = {
  {.name="A"}, {.name="B"}, {.name="C"}
};
static const chem_custom_parameter rates[] = {
  {.name="r1", .default_value=0.04},
  {.name="r2", .default_value=3.0e7},
  {.name="r3", .default_value=1.0e4}
};
static const uint32_t chain_reactants[] = {0, 1, 1, 1, 2};
static const uint32_t chain_products[] = {1, 2, 1, 2, 0, 2};
static const double chain_yields[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
static const chem_reaction chain_reactions[] = {
  {.kind=CHEM_USER_DEFINED, .reactant_offset=0, .reactant_count=1,
   .product_offset=0, .product_count=1, .custom0=0, .p={1.0}},
  {.kind=CHEM_USER_DEFINED, .reactant_offset=1, .reactant_count=2,
   .product_offset=1, .product_count=2, .custom0=1, .p={1.0}},
  {.kind=CHEM_USER_DEFINED, .reactant_offset=3, .reactant_count=2,
   .product_offset=3, .product_count=2, .custom0=2, .p={1.0}}
};
static const chem_mechanism stiff_chain = {
  .name="three-species chain", .species_count=3, .reaction_count=3,
  .reactant_count=5, .product_count=6, .custom_count=3,
  .default_temperature=272.5, .default_pressure=101253.3,
  .species=three_species, .reactions=chain_reactions,
  .reactants=chain_reactants, .products=chain_products,
  .yields=chain_yields, .custom_parameters=rates
};

static int finite_state(const double *values, size_t count)
{
  size_t i;
  for (i = 0; i < count; ++i)
    if (!isfinite(values[i]) || values[i] < -1.0e-10)
      return 0;
  return 1;
}

static int run_chain(const char *id)
{
  const size_t cells = strcmp(id, "multiple_grid_cells") == 0 ||
                       strcmp(id, "vectorized_matrix_solver") == 0 ? 3 : 1;
  double temperature[3] = {272.5, 272.5, 272.5};
  double pressure[3] = {101253.3, 101253.3, 101253.3};
  double custom[9], state[9] = {0};
  chem_rosenbrock_parameters parameters;
  chem_solver_stats stats;
  chem_solver *solver;
  size_t cell, step, steps = strcmp(id, "multiple_grid_cells") == 0 ? 10 : 1;
  double dt = strcmp(id, "multiple_grid_cells") == 0 ? 100.0 : 200.0;
  size_t group_size = strcmp(id, "vectorized_matrix_solver") == 0 ? 3 : 0;
  if (strcmp(id, "vectorized_matrix_solver") == 0)
    dt = 0.1;

  for (cell = 0; cell < cells; ++cell) {
    custom[cell * 3] = 0.04;
    custom[cell * 3 + 1] = 3.0e7;
    custom[cell * 3 + 2] = 1.0e4;
    if (strcmp(id, "multiple_grid_cells") == 0 ||
        strcmp(id, "vectorized_matrix_solver") == 0) {
      state[cell * 3] = 1.0 + cell;
      state[cell * 3 + 1] = 0.0;
      state[cell * 3 + 2] = 0.0;
    } else {
      state[cell * 3] = 1.0;
    }
  }
  if (strcmp(id, "solver_configuration") == 0)
    chem_rosenbrock_four_stage(&parameters);
  else
    chem_rosenbrock_default(&parameters);
  solver = chem_solver_create(&stiff_chain, cells, temperature, pressure,
                              custom, &parameters, group_size);
  if (!solver)
    return 0;
  for (step = 0; step < steps; ++step)
    if (chem_solver_step(solver, state, dt, &stats) != 0) {
      chem_solver_destroy(solver);
      return 0;
    }
  chem_solver_destroy(solver);
  if (!finite_state(state, cells * 3) || stats.accepted == 0)
    return 0;
  if (cells == 3 && !(state[0] != state[3] && state[3] != state[6]))
    return 0;
  printf("%s: integrated %zu cell(s), accepted %llu steps\n", id, cells,
         (unsigned long long)stats.accepted);
  return 1;
}

static int run_readme_example(void)
{
  static const chem_species species[] = {
    {.name="Foo"}, {.name="Bar"}, {.name="Baz"}
  };
  static const chem_reaction example_reactions[] = {
    {.kind=CHEM_ARRHENIUS, .reactant_offset=0, .reactant_count=1,
     .product_offset=0, .product_count=2, .p={1.0e-3, 0.0, 0.0, 300.0}},
    {.kind=CHEM_ARRHENIUS, .reactant_offset=1, .reactant_count=2,
     .product_offset=2, .product_count=1, .p={1.0e-5, 0.0, 110.0, 300.0}}
  };
  static const uint32_t reactants[] = {0, 0, 1};
  static const uint32_t products[] = {1, 2, 2};
  static const double yields[] = {0.8, 0.2, 1.0};
  static const chem_mechanism mechanism = {
    .name="README tutorial", .species_count=3, .reaction_count=2,
    .reactant_count=3, .product_count=3,
    .default_temperature=287.45, .default_pressure=101319.9,
    .species=species, .reactions=example_reactions,
    .reactants=reactants, .products=products, .yields=yields
  };
  double temperature = 287.45, pressure = 101319.9, state[] = {20.0, 0.0, 0.0};
  chem_rosenbrock_parameters parameters;
  chem_solver_stats stats;
  chem_solver *solver;
  int step;
  chem_rosenbrock_default(&parameters);
  solver = chem_solver_create(&mechanism, 1, &temperature, &pressure, NULL,
                              &parameters, 0);
  if (!solver)
    return 0;
  for (step = 0; step < 10; ++step)
    if (chem_solver_step(solver, state, 500.0, &stats) != 0) {
      chem_solver_destroy(solver);
      return 0;
    }
  chem_solver_destroy(solver);
  if (!finite_state(state, 3) || state[0] >= 20.0 || state[1] <= 0.0 ||
      state[2] <= 0.0)
    return 0;
  printf("README_example: final Foo %.12g Bar %.12g Baz %.12g\n",
         state[0], state[1], state[2]);
  return 1;
}

static int run_openmp_tutorial(void)
{
  double states[3][3] = {{0}};
  int thread, status[3] = {0};
#ifdef _OPENMP
#pragma omp parallel num_threads(3)
#endif
  {
#ifdef _OPENMP
    const int index = omp_get_thread_num();
#else
    const int index = 0;
#endif
    double temperature = 272.5, pressure = 101253.3;
    double custom[] = {0.04, 3.0e7, 1.0e4};
    chem_rosenbrock_parameters parameters;
    chem_solver_stats stats;
    chem_solver *solver;
    int step;
    states[index][0] = 1.0;
    chem_rosenbrock_default(&parameters);
    solver = chem_solver_create(&stiff_chain, 1, &temperature, &pressure,
                                custom, &parameters, 1);
    if (solver) {
      status[index] = 1;
      for (step = 0; step < 10; ++step)
        if (chem_solver_step(solver, states[index], 200.0, &stats) != 0) {
          status[index] = 0;
          break;
        }
      chem_solver_destroy(solver);
      status[index] = status[index] && finite_state(states[index], 3);
    }
  }
  for (thread = 0; thread < 3; ++thread)
    if (!status[thread] || states[thread][0] >= 1.0 || states[thread][1] <= 0.0)
      return 0;
  printf("openmp_tutorial: independent states converged on three worker threads\n");
  return 1;
}

int main(int argc, char **argv)
{
  const char *id;
  int ok;
  if (argc != 2)
    return 2;
  id = argv[1];
  if (strcmp(id, "README_example") == 0)
    ok = run_readme_example();
  else if (strcmp(id, "solve_results") == 0 ||
           strcmp(id, "multiple_grid_cells") == 0 ||
           strcmp(id, "solver_configuration") == 0 ||
           strcmp(id, "vectorized_matrix_solver") == 0)
    ok = run_chain(id);
  else if (strcmp(id, "openmp_tutorial") == 0)
    ok = run_openmp_tutorial();
  else
    return 2;
  if (!ok) {
    fprintf(stderr, "%s failed to produce a finite integrated state\n", id);
    return 1;
  }
  return 0;
}
