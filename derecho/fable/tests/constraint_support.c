/* Support for the ports of MICM's constraint tests; see constraint_support.h. */
#include "constraint_support.h"

#include "internal.h"

/* ---- systems, reactions, constraints ------------------------------------ */

tmc_stoich_species cst_stoich(const char *name, double coefficient)
{
  tmc_stoich_species s;
  s.species = tmc_species_create(name);
  s.coefficient = coefficient;
  return s;
}

tmc_system cst_system(const tmc_species *species, int n)
{
  tmc_phase gas = tmc_phase_create("gas");
  for (int i = 0; i < n; ++i)
    tmc_phase_add_species(&gas, &species[i]);
  tmc_system system = tmc_system_create(&gas);
  tmc_phase_free(&gas);
  return system;
}

tmc_reaction cst_reaction(const tmc_species *reactants, int n_reactants, const tmc_stoich_species *products, int n_products, tmc_rate_constant rate_constant)
{
  tmc_reaction r = tmc_reaction_create();
  for (int i = 0; i < n_reactants; ++i)
    tmc_reaction_add_reactant(&r, &reactants[i]);
  for (int i = 0; i < n_products; ++i)
    tmc_reaction_add_product(&r, &products[i].species, products[i].coefficient);
  tmc_reaction_set_rate_constant(&r, rate_constant);
  tmc_strlcpy(r.phase, "gas", sizeof r.phase);
  return r;
}

tmc_rate_constant cst_arrhenius(double A)
{
  tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
  p.A = A;
  p.B = 0.0;
  p.C = 0.0;
  return tmc_rate_arrhenius(p);
}

tmc_constraint cst_equilibrium(const char *name, const tmc_species *algebraic, const tmc_stoich_species *reactants, int n_reactants, const tmc_stoich_species *products, int n_products, double K_HLC_ref, double delta_H)
{
  tmc_constraint c;
  tmc_error err;
  if (tmc_equilibrium_constraint(&c, name, algebraic, reactants, n_reactants, products, n_products, (tmc_vant_hoff){ K_HLC_ref, delta_H, 0.0, 0.0 }, &err) != TMC_OK)
    tmc_test_fail(__FILE__, __LINE__, "equilibrium constraint '%s': %s", name, err.message);
  return c;
}

tmc_constraint cst_linear(const char *name, const tmc_species *algebraic, const tmc_stoich_species *terms, int n_terms, double constant)
{
  tmc_constraint c;
  tmc_error err;
  if (tmc_linear_constraint(&c, name, algebraic, terms, n_terms, constant, &err) != TMC_OK)
    tmc_test_fail(__FILE__, __LINE__, "linear constraint '%s': %s", name, err.message);
  return c;
}

tmc_solver *cst_build(tmc_rosenbrock_parameters parameters, const tmc_system *system, const tmc_reaction *reactions, int n_reactions, const tmc_constraint *constraints, int n_constraints, int L, bool reorder_state, tmc_error *err)
{
  tmc_builder *b = tmc_builder_create_rosenbrock(parameters);
  tmc_builder_set_system(b, system);
  tmc_builder_set_reactions(b, reactions, n_reactions);
  for (int i = 0; i < n_constraints; ++i)
    tmc_builder_add_constraint(b, &constraints[i]);
  tmc_builder_set_reorder_state(b, reorder_state);
  tmc_builder_set_width(b, L);
  tmc_builder_set_lu_algorithm(b, TMC_LU_DOOLITTLE);
  tmc_solver *s = tmc_builder_build(b, err);
  tmc_builder_free(b);
  return s;
}

double cst_equilibrium_constant(double K_HLC_ref, double delta_H, double temperature)
{
  const double gas_constant = 1.380649e-23 * 6.02214076e23;
  return K_HLC_ref * exp((delta_H / gas_constant) * (1.0 / temperature - 1.0 / 298.15));
}

/* ---- MICM's stand-alone ConstraintSet ------------------------------------ */

tmc_solver *cst_constraint_set(const char *const *species, int n_species, const tmc_constraint *constraints, int n_constraints, int L, tmc_error *err)
{
  tmc_species *list = tmc_malloc((size_t)(n_species > 0 ? n_species : 1) * sizeof *list);
  for (int i = 0; i < n_species; ++i)
    list[i] = tmc_species_create(species[i]);
  tmc_system system = cst_system(list, n_species);
  free(list);
  tmc_builder *b = tmc_builder_create_rosenbrock(tmc_rosenbrock_four_stage_dae());
  tmc_builder_set_system(b, &system);
  for (int i = 0; i < n_constraints; ++i)
    tmc_builder_add_constraint(b, &constraints[i]);
  tmc_builder_set_reorder_state(b, false);
  tmc_builder_set_ignore_unused_species(b, true);
  tmc_builder_set_width(b, L);
  tmc_solver *s = tmc_builder_build(b, err);
  tmc_builder_free(b);
  tmc_system_free(&system);
  for (int i = 0; s && i < n_species; ++i)
    if (strcmp(tmc_solver_variable_name(s, i), species[i]) != 0)
    {
      tmc_error_set(err, TMC_ERROR_INTERNAL, TMC_INTERNAL_GENERAL, "species '%s' is not at index %d", species[i], i);
      tmc_solver_free(s);
      s = NULL;
    }
  return s;
}

tmc_index_pair *cst_nonzero_elements(const tmc_solver *s, int *n)
{
  *n = tmc_solver_constraint_jacobian_elements(s, NULL, 0);
  tmc_index_pair *pairs = tmc_malloc((size_t)(*n > 0 ? *n : 1) * sizeof *pairs);
  tmc_solver_constraint_jacobian_elements(s, pairs, *n);
  return pairs;
}

bool cst_has_pair(const tmc_index_pair *pairs, int n, int row, int col)
{
  for (int i = 0; i < n; ++i)
    if (pairs[i].row == row && pairs[i].col == col)
      return true;
  return false;
}

int cst_parameter_index(const tmc_solver *s, const char *label)
{
  return tmc_name_map_lookup(tmc_solver_parameter_map(s), label);
}

tmc_dense *cst_parameters(const tmc_solver *s, int n_cells, double value)
{
  return tmc_dense_create(n_cells, tmc_solver_num_custom_parameters(s), tmc_solver_width(s), value);
}

tmc_sparse *cst_jacobian(const tmc_solver *s, int n_blocks)
{
  return tmc_sparse_create_like(tmc_solver_lu_pattern(s), n_blocks, 0.0);
}

double cst_jacobian_at(const tmc_sparse *jacobian, int block, int row, int col, const char *file, int line)
{
  if (tmc_sparse_is_zero(jacobian, row, col))
  {
    tmc_test_fail(file, line, "Jacobian element (%d, %d) is not in the sparsity pattern", row, col);
    return NAN;
  }
  return tmc_sparse_get(jacobian, block, row, col);
}

void cst_add_forcing_terms(const tmc_solver *s, const tmc_dense *variables, const tmc_dense *parameters, tmc_dense *forcing)
{
  tmc_constraints_add_forcing(s, variables, parameters, forcing);
}

void cst_subtract_jacobian_terms(const tmc_solver *s, const tmc_dense *variables, const tmc_dense *parameters, tmc_sparse *jacobian)
{
  tmc_constraints_subtract_jacobian(s, variables, parameters, jacobian);
}

/* ---- finite-difference Jacobians ---------------------------------------- */

void cst_constraint_forcing(void *context, const tmc_dense *variables, tmc_dense *forcing)
{
  const cst_forcing_context *c = context;
  cst_add_forcing_terms(c->solver, variables, c->parameters, forcing);
}

tmc_dense *cst_finite_difference_jacobian(cst_forcing_function forcing, void *context, const tmc_dense *base, int n_species, double perturbation)
{
  const int n_blocks = base->n_rows, L = base->L;
  tmc_dense *result = tmc_dense_create(n_blocks, n_species * n_species, L, 0.0);
  tmc_dense *plus = tmc_dense_create(n_blocks, n_species, L, 0.0);
  tmc_dense *minus = tmc_dense_create(n_blocks, n_species, L, 0.0);
  tmc_dense *forcing_plus = tmc_dense_create(n_blocks, n_species, L, 0.0);
  tmc_dense *forcing_minus = tmc_dense_create(n_blocks, n_species, L, 0.0);
  tmc_dense *forcing_base = tmc_dense_create(n_blocks, n_species, L, 0.0);
  for (int col = 0; col < n_species; ++col)
  {
    tmc_dense_copy_into(plus, base);
    tmc_dense_copy_into(minus, base);
    for (int block = 0; block < n_blocks; ++block)
    {
      const double x = tmc_dense_get(base, block, col);
      const double h = perturbation * fmax(1.0, fabs(x));
      tmc_dense_set(plus, block, col, x + h);
      tmc_dense_set(minus, block, col, x - h);
    }
    tmc_dense_fill(forcing_plus, 0.0);
    tmc_dense_fill(forcing_minus, 0.0);
    forcing(context, plus, forcing_plus);
    forcing(context, minus, forcing_minus);
    for (int block = 0; block < n_blocks; ++block)
    {
      const double x = tmc_dense_get(base, block, col);
      const double h = perturbation * fmax(1.0, fabs(x));
      if (x - h >= 0.0)
      {
        const double two_h = 2.0 * h;
        for (int row = 0; row < n_species; ++row)
          tmc_dense_set(result, block, row * n_species + col, (tmc_dense_get(forcing_plus, block, row) - tmc_dense_get(forcing_minus, block, row)) / two_h);
      }
      else
      {
        /* a negative perturbation would leave the domain: forward differences */
        tmc_dense_fill(forcing_base, 0.0);
        forcing(context, base, forcing_base);
        for (int row = 0; row < n_species; ++row)
          tmc_dense_set(result, block, row * n_species + col, (tmc_dense_get(forcing_plus, block, row) - tmc_dense_get(forcing_base, block, row)) / h);
      }
    }
  }
  tmc_dense_free(plus);
  tmc_dense_free(minus);
  tmc_dense_free(forcing_plus);
  tmc_dense_free(forcing_minus);
  tmc_dense_free(forcing_base);
  return result;
}

cst_jacobian_comparison cst_compare_jacobian_to_finite_difference(const tmc_sparse *analytical, const tmc_dense *fd, int n_species, double atol, double rtol)
{
  cst_jacobian_comparison r;
  memset(&r, 0, sizeof r);
  r.passed = true;
  for (int block = 0; block < analytical->n_blocks; ++block)
    for (int row = 0; row < n_species; ++row)
      for (int col = 0; col < n_species; ++col)
      {
        if (tmc_sparse_is_zero(analytical, row, col))
          continue;
        /* the analytical Jacobian holds -df/dx */
        const double a = -tmc_sparse_get(analytical, block, row, col);
        const double f = tmc_dense_get(fd, block, row * n_species + col);
        const double abs_error = fabs(a - f);
        const double scale = fmax(fabs(a), fabs(f));
        const double rel_error = scale > 0.0 ? abs_error / scale : 0.0;
        if (abs_error > r.max_abs_error)
        {
          r.max_abs_error = abs_error;
          r.worst_block = block;
          r.worst_row = row;
          r.worst_col = col;
          r.worst_analytical = a;
          r.worst_fd = f;
        }
        if (rel_error > r.max_rel_error)
          r.max_rel_error = rel_error;
        if (abs_error > atol + rtol * scale)
          r.passed = false;
      }
  return r;
}

cst_jacobian_comparison cst_check_jacobian_sparsity_completeness(const tmc_sparse *analytical, const tmc_dense *fd, int n_species, double threshold)
{
  cst_jacobian_comparison r;
  memset(&r, 0, sizeof r);
  r.passed = true;
  for (int block = 0; block < analytical->n_blocks; ++block)
    for (int row = 0; row < n_species; ++row)
      for (int col = 0; col < n_species; ++col)
      {
        if (!tmc_sparse_is_zero(analytical, row, col))
          continue;
        const double f = fabs(tmc_dense_get(fd, block, row * n_species + col));
        if (f > threshold)
        {
          r.passed = false;
          if (f > r.max_abs_error)
          {
            r.max_abs_error = f;
            r.worst_block = block;
            r.worst_row = row;
            r.worst_col = col;
            r.worst_analytical = 0.0;
            r.worst_fd = f;
          }
        }
      }
  return r;
}
