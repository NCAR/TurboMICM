/* MICM's stub external models (test/integration/stub_aerosol_1.hpp,
 * stub_aerosol_2.hpp, stub_aerosol_with_constraints.hpp); see stub_aerosol.h. */
#include "stub_aerosol.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- helpers ---------------------------------------------------------------- */

static void set_name(char dst[TMC_NAME_MAX], const char *src)
{
  size_t n = strlen(src);
  if (n >= TMC_NAME_MAX)
    n = TMC_NAME_MAX - 1;
  memcpy(dst, src, n);
  dst[n] = '\0';
}

void stub_name_set_insert(stub_name_set *set, const char *name)
{
  if (set->n < STUB_MAX_NAMES)
    set_name(set->names[set->n++], name);
}

static int compare_names(const void *a, const void *b)
{
  return strcmp((const char *)a, (const char *)b);
}

int stub_name_set_finish(stub_name_set *set, char (*names)[TMC_NAME_MAX], int max)
{
  qsort(set->names, (size_t)set->n, sizeof set->names[0], compare_names);
  int k = 0;
  for (int i = 0; i < set->n; ++i)
    if (k == 0 || strcmp(set->names[k - 1], set->names[i]) != 0)
    {
      if (k != i)
        memcpy(set->names[k], set->names[i], TMC_NAME_MAX);
      ++k;
    }
  set->n = k;
  if (names)
    for (int i = 0; i < k && i < max; ++i)
      memcpy(names[i], set->names[i], TMC_NAME_MAX);
  return k;
}

void stub_pair_list_add(stub_pair_list *list, int row, int col)
{
  if (list->n < (int)(sizeof list->pairs / sizeof list->pairs[0]))
  {
    list->pairs[list->n].row = row;
    list->pairs[list->n].col = col;
    ++list->n;
  }
}

int stub_pair_list_finish(const stub_pair_list *list, tmc_index_pair *out, int max)
{
  if (out)
    for (int i = 0; i < list->n && i < max; ++i)
      out[i] = list->pairs[i];
  return list->n;
}

int stub_dense_rows(const tmc_dense *m)
{
  return m->n_groups * m->L;
}

int stub_sparse_blocks(const tmc_sparse *m)
{
  return m->n_groups * m->L;
}

double *stub_dense_at(tmc_dense *m, int row, int col)
{
  return &m->data[tmc_dense_index(m, row, col)];
}

double stub_dense_value(const tmc_dense *m, int row, int col)
{
  return m->data[tmc_dense_index(m, row, col)];
}

double *stub_sparse_at(tmc_sparse *m, int block, int element)
{
  return &m->data[tmc_sparse_data_index(m, block, element)];
}

stub_phase stub_phase_from(const tmc_phase *phase)
{
  stub_phase p;
  memset(&p, 0, sizeof p);
  set_name(p.name, phase->name);
  int n = tmc_phase_species_names(phase, true, NULL, 0);
  if (n > STUB_MAX_PHASE_SPECIES)
    n = STUB_MAX_PHASE_SPECIES;
  p.n_unique_names = tmc_phase_species_names(phase, true, p.unique_names, n);
  if (p.n_unique_names > n)
    p.n_unique_names = n;
  return p;
}

/* "<prefix>.<suffix>" */
static void join_name(char out[TMC_NAME_MAX], const char *prefix, const char *suffix)
{
  char buffer[4 * TMC_NAME_MAX];
  snprintf(buffer, sizeof buffer, "%s.%s", prefix, suffix);
  set_name(out, buffer);
}

/* "<name>.MODE<mode + 1>" */
static void mode_name(char out[TMC_NAME_MAX], const char *name, int mode)
{
  char buffer[4 * TMC_NAME_MAX];
  snprintf(buffer, sizeof buffer, "%s.MODE%d", name, mode + 1);
  set_name(out, buffer);
}

/* name_ + ".MODE" + std::to_string(mode + 1) + "." + phase.name_ + "." + species.name_ */
static void species_name(char out[TMC_NAME_MAX], const char *name, int mode, const tmc_phase *phase, const tmc_species *species)
{
  char buffer[4 * TMC_NAME_MAX];
  snprintf(buffer, sizeof buffer, "%s.MODE%d.%s.%s", name, mode + 1, phase->name, species->name);
  set_name(out, buffer);
}

/* Insert "<name>.MODE<mode + 1>.<unique name>" for each unique name of the phase. */
static void insert_mode_names(stub_name_set *set, const char *name, int mode, const stub_phase *phase)
{
  char prefix[TMC_NAME_MAX], full[TMC_NAME_MAX];
  mode_name(prefix, name, mode);
  for (int i = 0; i < phase->n_unique_names; ++i)
  {
    join_name(full, prefix, phase->unique_names[i]);
    stub_name_set_insert(set, full);
  }
}

static int no_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  (void)context;
  (void)names;
  (void)max;
  return 0;
}

static void no_parameter_update(void *context, const tmc_conditions *conditions, int n_cells, tmc_dense *state_parameters)
{
  (void)context;
  (void)conditions;
  (void)n_cells;
  (void)state_parameters;
}

/* The first-order transfers both aerosol stubs share: rate = k reactant, the
 * reactant loses it and the product gains it. k is the fixed rate constant,
 * or the state parameter rate_param_index when that is not negative. */
static void transfer_forcing(const stub_forcing_info *info, int n_info, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_dense *forcing)
{
  const int rows = stub_dense_rows(forcing);
  for (int i = 0; i < n_info; ++i)
  {
    const stub_forcing_info *f = &info[i];
    for (int row = 0; row < rows; ++row)
    {
      const double k = f->rate_param_index >= 0 ? stub_dense_value(state_parameters, row, f->rate_param_index) : f->rate_constant;
      const double rate = k * stub_dense_value(state_variables, row, f->reactant_index);
      *stub_dense_at(forcing, row, f->reactant_index) -= rate;
      *stub_dense_at(forcing, row, f->product_index) += rate;
    }
  }
}

static void transfer_jacobian(const stub_forcing_info *info, int n_info, const tmc_dense *state_parameters, tmc_sparse *jacobian)
{
  const int blocks = stub_sparse_blocks(jacobian);
  for (int i = 0; i < n_info; ++i)
  {
    const stub_forcing_info *f = &info[i];
    if (f->reactant_reactant_flat < 0 || f->product_reactant_flat < 0)
      continue;
    for (int block = 0; block < blocks; ++block)
    {
      const double k = f->rate_param_index >= 0 ? stub_dense_value(state_parameters, block, f->rate_param_index) : f->rate_constant;
      *stub_sparse_at(jacobian, block, f->reactant_reactant_flat) -= -k;
      *stub_sparse_at(jacobian, block, f->product_reactant_flat) -= k;
    }
  }
}

/* forcing_info_.push_back({ reactant, product, k or parameter, flat (r, r), flat (p, r) })
 * when both species (and the parameter, if any) are known */
static void add_transfer(stub_forcing_info *info, int *n_info, const tmc_name_map *variables, const tmc_sparse *jacobian, const char *reactant, const char *product, double rate_constant, int rate_param_index)
{
  int r = tmc_name_map_lookup(variables, reactant);
  int p = tmc_name_map_lookup(variables, product);
  if (r < 0 || p < 0)
    return;
  stub_forcing_info *f = &info[(*n_info)++];
  f->reactant_index = r;
  f->product_index = p;
  f->rate_constant = rate_constant;
  f->rate_param_index = rate_param_index;
  f->reactant_reactant_flat = tmc_sparse_element(jacobian, r, r);
  f->product_reactant_flat = tmc_sparse_element(jacobian, p, r);
}

/* ---- StubAerosolModel ------------------------------------------------------- */

void stub_aerosol_1_init(stub_aerosol_1 *m, const char *name, const tmc_phase *phases, int n_phases, stub1_rate_constants rate_constants)
{
  memset(m, 0, sizeof *m);
  set_name(m->name, name);
  m->n_phases = n_phases < 2 ? n_phases : 2;
  for (int i = 0; i < m->n_phases; ++i)
    m->phases[i] = stub_phase_from(&phases[i]);
  m->rate_constants = rate_constants;
}

void stub_aerosol_1_state_size(const stub_aerosol_1 *m, int *n_variables, int *n_parameters)
{
  *n_variables = m->phases[0].n_unique_names + m->phases[0].n_unique_names + m->phases[1].n_unique_names;
  *n_parameters = 0;
}

int stub_aerosol_1_state_variable_names(const stub_aerosol_1 *m, char (*names)[TMC_NAME_MAX], int max)
{
  stub_name_set set;
  set.n = 0;
  insert_mode_names(&set, m->name, 0, &m->phases[0]);
  insert_mode_names(&set, m->name, 1, &m->phases[0]);
  insert_mode_names(&set, m->name, 1, &m->phases[1]);
  return stub_name_set_finish(&set, names, max);
}

int stub_aerosol_1_state_parameter_names(const stub_aerosol_1 *m, char (*names)[TMC_NAME_MAX], int max)
{
  (void)m;
  return no_names(NULL, names, max);
}

void stub_aerosol_1_species(const stub_aerosol_1 *m, int mode, const tmc_phase *phase, const tmc_species *species, char out[TMC_NAME_MAX])
{
  species_name(out, m->name, mode, phase, species);
}

static int s1_state_variable_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  return stub_aerosol_1_state_variable_names(context, names, max);
}

static int s1_state_parameter_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  return stub_aerosol_1_state_parameter_names(context, names, max);
}

/* SpeciesUsed() is StateVariableNames() */
static int s1_species_used(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  return stub_aerosol_1_state_variable_names(context, names, max);
}

static int s1_nonzero_jacobian_elements(void *context, const tmc_name_map *variables, tmc_index_pair *out, int max)
{
  (void)context;
  stub_pair_list list;
  list.n = 0;
  int fo2_gas = tmc_name_map_lookup(variables, "FO2");
  int fo2_mode2 = tmc_name_map_lookup(variables, "STUB1.MODE2.CORGE.FO2");
  if (fo2_gas >= 0 && fo2_mode2 >= 0)
  {
    stub_pair_list_add(&list, fo2_gas, fo2_gas);
    stub_pair_list_add(&list, fo2_mode2, fo2_gas);
  }
  int baz_mode1 = tmc_name_map_lookup(variables, "STUB1.MODE1.QUUX.BAZ");
  int baz_mode2 = tmc_name_map_lookup(variables, "STUB1.MODE2.QUUX.BAZ");
  if (baz_mode1 >= 0 && baz_mode2 >= 0)
  {
    stub_pair_list_add(&list, baz_mode1, baz_mode1);
    stub_pair_list_add(&list, baz_mode2, baz_mode1);
  }
  return stub_pair_list_finish(&list, out, max);
}

static void s1_finalize_process_setup(void *context, const tmc_name_map *parameters, const tmc_name_map *variables, const tmc_sparse *jacobian)
{
  stub_aerosol_1 *m = context;
  (void)parameters;
  m->n_forcing = 0;
  add_transfer(m->forcing_info, &m->n_forcing, variables, jacobian, "FO2", "STUB1.MODE2.CORGE.FO2", m->rate_constants.fo2_gas_to_mode2_corge, -1);
  add_transfer(m->forcing_info, &m->n_forcing, variables, jacobian, "STUB1.MODE1.QUUX.BAZ", "STUB1.MODE2.QUUX.BAZ", m->rate_constants.baz_mode1_to_mode2_quux, -1);
}

static void s1_add_forcing_terms(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_dense *forcing)
{
  const stub_aerosol_1 *m = context;
  transfer_forcing(m->forcing_info, m->n_forcing, state_parameters, state_variables, forcing);
}

static void s1_subtract_jacobian_terms(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_sparse *jacobian)
{
  const stub_aerosol_1 *m = context;
  (void)state_variables;
  transfer_jacobian(m->forcing_info, m->n_forcing, state_parameters, jacobian);
}

tmc_external_model stub_aerosol_1_model(stub_aerosol_1 *m)
{
  tmc_external_model model;
  memset(&model, 0, sizeof model);
  model.context = m;
  model.state_variable_names = s1_state_variable_names;
  model.state_parameter_names = s1_state_parameter_names;
  model.species_used = s1_species_used;
  model.nonzero_jacobian_elements = s1_nonzero_jacobian_elements;
  model.finalize_process_setup = s1_finalize_process_setup;
  model.update_state_parameters = no_parameter_update;
  model.add_forcing_terms = s1_add_forcing_terms;
  model.subtract_jacobian_terms = s1_subtract_jacobian_terms;
  return model;
}

/* ---- AnotherStubAerosolModel ------------------------------------------------ */

void stub_aerosol_2_init(stub_aerosol_2 *m, const char *name, const tmc_phase *phases, int n_phases)
{
  memset(m, 0, sizeof *m);
  set_name(m->name, name);
  m->n_phases = n_phases < 2 ? n_phases : 2;
  for (int i = 0; i < m->n_phases; ++i)
    m->phases[i] = stub_phase_from(&phases[i]);
  m->baz_to_qux_param_index = -1;
}

void stub_aerosol_2_state_size(const stub_aerosol_2 *m, int *n_variables, int *n_parameters)
{
  int size = 0;
  size += 1;
  size += m->phases[0].n_unique_names;
  size += 1;
  size += m->phases[1].n_unique_names;
  size += 1;
  size += m->phases[0].n_unique_names;
  size += m->phases[1].n_unique_names;
  *n_variables = size;
  *n_parameters = 2;
}

void stub_aerosol_2_number(const stub_aerosol_2 *m, int mode, char out[TMC_NAME_MAX])
{
  char prefix[TMC_NAME_MAX];
  mode_name(prefix, m->name, mode);
  join_name(out, prefix, "NUMBER");
}

int stub_aerosol_2_state_variable_names(const stub_aerosol_2 *m, char (*names)[TMC_NAME_MAX], int max)
{
  stub_name_set set;
  set.n = 0;
  char number[TMC_NAME_MAX];
  stub_aerosol_2_number(m, 0, number);
  stub_name_set_insert(&set, number);
  insert_mode_names(&set, m->name, 0, &m->phases[0]);
  stub_aerosol_2_number(m, 1, number);
  stub_name_set_insert(&set, number);
  insert_mode_names(&set, m->name, 1, &m->phases[1]);
  stub_aerosol_2_number(m, 2, number);
  stub_name_set_insert(&set, number);
  insert_mode_names(&set, m->name, 2, &m->phases[0]);
  insert_mode_names(&set, m->name, 2, &m->phases[1]);
  return stub_name_set_finish(&set, names, max);
}

static void fo2_to_baz_parameter(const stub_aerosol_2 *m, char out[TMC_NAME_MAX])
{
  join_name(out, m->name, "PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT");
}

static void baz_to_qux_parameter(const stub_aerosol_2 *m, char out[TMC_NAME_MAX])
{
  join_name(out, m->name, "PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT");
}

int stub_aerosol_2_state_parameter_names(const stub_aerosol_2 *m, char (*names)[TMC_NAME_MAX], int max)
{
  stub_name_set set;
  set.n = 0;
  char label[TMC_NAME_MAX];
  fo2_to_baz_parameter(m, label);
  stub_name_set_insert(&set, label);
  baz_to_qux_parameter(m, label);
  stub_name_set_insert(&set, label);
  return stub_name_set_finish(&set, names, max);
}

void stub_aerosol_2_species(const stub_aerosol_2 *m, int mode, const tmc_phase *phase, const tmc_species *species, char out[TMC_NAME_MAX])
{
  species_name(out, m->name, mode, phase, species);
}

static int s2_state_variable_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  return stub_aerosol_2_state_variable_names(context, names, max);
}

static int s2_state_parameter_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  return stub_aerosol_2_state_parameter_names(context, names, max);
}

static int s2_species_used(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  return stub_aerosol_2_state_variable_names(context, names, max);
}

static int s2_nonzero_jacobian_elements(void *context, const tmc_name_map *variables, tmc_index_pair *out, int max)
{
  (void)context;
  stub_pair_list list;
  list.n = 0;
  int fo2_mode2 = tmc_name_map_lookup(variables, "STUB2.MODE2.CORGE.FO2");
  int baz_mode2 = tmc_name_map_lookup(variables, "STUB2.MODE2.CORGE.BAZ");
  if (fo2_mode2 >= 0 && baz_mode2 >= 0)
  {
    stub_pair_list_add(&list, fo2_mode2, fo2_mode2);
    stub_pair_list_add(&list, baz_mode2, fo2_mode2);
  }
  int baz_mode3 = tmc_name_map_lookup(variables, "STUB2.MODE3.QUUX.BAZ");
  int qux_mode3 = tmc_name_map_lookup(variables, "STUB2.MODE3.QUUX.QUX");
  if (baz_mode3 >= 0 && qux_mode3 >= 0)
  {
    stub_pair_list_add(&list, baz_mode3, baz_mode3);
    stub_pair_list_add(&list, qux_mode3, baz_mode3);
  }
  return stub_pair_list_finish(&list, out, max);
}

static void s2_finalize_process_setup(void *context, const tmc_name_map *parameters, const tmc_name_map *variables, const tmc_sparse *jacobian)
{
  stub_aerosol_2 *m = context;
  char fo2_to_baz[TMC_NAME_MAX], baz_to_qux[TMC_NAME_MAX];
  fo2_to_baz_parameter(m, fo2_to_baz);
  baz_to_qux_parameter(m, baz_to_qux);
  m->n_forcing = 0;
  int p = tmc_name_map_lookup(parameters, fo2_to_baz);
  if (p >= 0)
    add_transfer(m->forcing_info, &m->n_forcing, variables, jacobian, "STUB2.MODE2.CORGE.FO2", "STUB2.MODE2.CORGE.BAZ", 0.0, p);
  p = tmc_name_map_lookup(parameters, baz_to_qux);
  if (p >= 0)
    add_transfer(m->forcing_info, &m->n_forcing, variables, jacobian, "STUB2.MODE3.QUUX.BAZ", "STUB2.MODE3.QUUX.QUX", 0.0, p);
  m->baz_to_qux_param_index = tmc_name_map_lookup(parameters, baz_to_qux);
}

/* the BAZ -> QUX rate constant is 0.005 T */
static void s2_update_state_parameters(void *context, const tmc_conditions *conditions, int n_cells, tmc_dense *state_parameters)
{
  const stub_aerosol_2 *m = context;
  if (m->baz_to_qux_param_index < 0)
    return;
  for (int cell = 0; cell < n_cells; ++cell)
    *stub_dense_at(state_parameters, cell, m->baz_to_qux_param_index) = 0.005 * conditions[cell].temperature;
}

static void s2_add_forcing_terms(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_dense *forcing)
{
  const stub_aerosol_2 *m = context;
  transfer_forcing(m->forcing_info, m->n_forcing, state_parameters, state_variables, forcing);
}

static void s2_subtract_jacobian_terms(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_sparse *jacobian)
{
  const stub_aerosol_2 *m = context;
  (void)state_variables;
  transfer_jacobian(m->forcing_info, m->n_forcing, state_parameters, jacobian);
}

tmc_external_model stub_aerosol_2_model(stub_aerosol_2 *m)
{
  tmc_external_model model;
  memset(&model, 0, sizeof model);
  model.context = m;
  model.state_variable_names = s2_state_variable_names;
  model.state_parameter_names = s2_state_parameter_names;
  model.species_used = s2_species_used;
  model.nonzero_jacobian_elements = s2_nonzero_jacobian_elements;
  model.finalize_process_setup = s2_finalize_process_setup;
  model.update_state_parameters = s2_update_state_parameters;
  model.add_forcing_terms = s2_add_forcing_terms;
  model.subtract_jacobian_terms = s2_subtract_jacobian_terms;
  return model;
}

/* ---- StubAerosolWithConstraints --------------------------------------------- */

void stub_aerosol_constraints_init(stub_aerosol_constraints *m, double rate_constant)
{
  memset(m, 0, sizeof *m);
  m->rate_constant = rate_constant;
  m->has_total_mass = false;
  m->i_gas = m->i_aq = -1;
  m->gas_gas_flat = m->aq_gas_flat = m->aq_aq_flat = -1;
}

void stub_aerosol_constraints_init_with_total(stub_aerosol_constraints *m, double rate_constant, double total_mass)
{
  stub_aerosol_constraints_init(m, rate_constant);
  m->has_total_mass = true;
  m->total_mass = total_mass;
}

static int sc_state_variable_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  (void)context;
  stub_name_set set;
  set.n = 0;
  stub_name_set_insert(&set, "AEROSOL.A_AQ");
  return stub_name_set_finish(&set, names, max);
}

static int sc_species_used(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  (void)context;
  stub_name_set set;
  set.n = 0;
  stub_name_set_insert(&set, "A_GAS");
  stub_name_set_insert(&set, "AEROSOL.A_AQ");
  return stub_name_set_finish(&set, names, max);
}

static int sc_nonzero_jacobian_elements(void *context, const tmc_name_map *variables, tmc_index_pair *out, int max)
{
  (void)context;
  stub_pair_list list;
  list.n = 0;
  int i_gas = tmc_name_map_lookup(variables, "A_GAS");
  int i_aq = tmc_name_map_lookup(variables, "AEROSOL.A_AQ");
  if (i_gas >= 0 && i_aq >= 0)
  {
    stub_pair_list_add(&list, i_gas, i_gas);
    stub_pair_list_add(&list, i_aq, i_gas);
  }
  return stub_pair_list_finish(&list, out, max);
}

static void sc_finalize_process_setup(void *context, const tmc_name_map *parameters, const tmc_name_map *variables, const tmc_sparse *jacobian)
{
  stub_aerosol_constraints *m = context;
  (void)parameters;
  m->i_gas = tmc_name_map_lookup(variables, "A_GAS");
  m->i_aq = tmc_name_map_lookup(variables, "AEROSOL.A_AQ");
  m->gas_gas_flat = tmc_sparse_element(jacobian, m->i_gas, m->i_gas);
  m->aq_gas_flat = tmc_sparse_element(jacobian, m->i_aq, m->i_gas);
}

static void sc_add_forcing_terms(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_dense *forcing)
{
  const stub_aerosol_constraints *m = context;
  (void)state_parameters;
  const int rows = stub_dense_rows(forcing);
  for (int row = 0; row < rows; ++row)
  {
    const double rate = m->rate_constant * stub_dense_value(state_variables, row, m->i_gas);
    *stub_dense_at(forcing, row, m->i_gas) -= rate;
    *stub_dense_at(forcing, row, m->i_aq) += rate;
  }
}

static void sc_subtract_jacobian_terms(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_sparse *jacobian)
{
  const stub_aerosol_constraints *m = context;
  (void)state_parameters;
  (void)state_variables;
  if (m->gas_gas_flat < 0 || m->aq_gas_flat < 0)
    return;
  const double k = m->rate_constant;
  const int blocks = stub_sparse_blocks(jacobian);
  for (int block = 0; block < blocks; ++block)
  {
    *stub_sparse_at(jacobian, block, m->gas_gas_flat) -= -k;
    *stub_sparse_at(jacobian, block, m->aq_gas_flat) -= k;
  }
}

static int sc_constraint_algebraic_variable_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  const stub_aerosol_constraints *m = context;
  stub_name_set set;
  set.n = 0;
  if (m->has_total_mass)
    stub_name_set_insert(&set, "AEROSOL.A_AQ");
  return stub_name_set_finish(&set, names, max);
}

static int sc_constraint_species_dependencies(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  const stub_aerosol_constraints *m = context;
  stub_name_set set;
  set.n = 0;
  if (m->has_total_mass)
  {
    stub_name_set_insert(&set, "A_GAS");
    stub_name_set_insert(&set, "AEROSOL.A_AQ");
  }
  return stub_name_set_finish(&set, names, max);
}

static int sc_nonzero_constraint_jacobian_elements(void *context, const tmc_name_map *variables, tmc_index_pair *out, int max)
{
  const stub_aerosol_constraints *m = context;
  stub_pair_list list;
  list.n = 0;
  if (m->has_total_mass)
  {
    int i_gas = tmc_name_map_lookup(variables, "A_GAS");
    int i_aq = tmc_name_map_lookup(variables, "AEROSOL.A_AQ");
    if (i_gas < 0 || i_aq < 0) /* MICM's at() throws */
      return 0;
    stub_pair_list_add(&list, i_aq, i_gas);
    stub_pair_list_add(&list, i_aq, i_aq);
  }
  return stub_pair_list_finish(&list, out, max);
}

static void sc_finalize_constraint_setup(void *context, const tmc_name_map *parameters, const tmc_name_map *variables, const tmc_sparse *jacobian)
{
  stub_aerosol_constraints *m = context;
  (void)parameters;
  if (!m->has_total_mass)
    return;
  m->i_gas = tmc_name_map_lookup(variables, "A_GAS");
  m->i_aq = tmc_name_map_lookup(variables, "AEROSOL.A_AQ");
  m->aq_gas_flat = tmc_sparse_element(jacobian, m->i_aq, m->i_gas);
  m->aq_aq_flat = tmc_sparse_element(jacobian, m->i_aq, m->i_aq);
}

/* G(y) = [A_GAS] + [A_AQ] - total = 0 */
static void sc_add_constraint_residual(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_dense *forcing)
{
  const stub_aerosol_constraints *m = context;
  (void)state_parameters;
  if (!m->has_total_mass)
    return;
  const int rows = stub_dense_rows(forcing);
  for (int row = 0; row < rows; ++row)
    *stub_dense_at(forcing, row, m->i_aq) = stub_dense_value(state_variables, row, m->i_gas) + stub_dense_value(state_variables, row, m->i_aq) - m->total_mass;
}

static void sc_subtract_constraint_jacobian(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_sparse *jacobian)
{
  const stub_aerosol_constraints *m = context;
  (void)state_parameters;
  (void)state_variables;
  if (!m->has_total_mass || m->aq_gas_flat < 0 || m->aq_aq_flat < 0)
    return;
  const int blocks = stub_sparse_blocks(jacobian);
  for (int block = 0; block < blocks; ++block)
  {
    *stub_sparse_at(jacobian, block, m->aq_gas_flat) -= 1.0;
    *stub_sparse_at(jacobian, block, m->aq_aq_flat) -= 1.0;
  }
}

tmc_external_model stub_aerosol_constraints_model(stub_aerosol_constraints *m)
{
  tmc_external_model model;
  memset(&model, 0, sizeof model);
  model.context = m;
  model.state_variable_names = sc_state_variable_names;
  model.state_parameter_names = no_names;
  model.species_used = sc_species_used;
  model.nonzero_jacobian_elements = sc_nonzero_jacobian_elements;
  model.finalize_process_setup = sc_finalize_process_setup;
  model.update_state_parameters = no_parameter_update;
  model.add_forcing_terms = sc_add_forcing_terms;
  model.subtract_jacobian_terms = sc_subtract_jacobian_terms;
  model.constraint_algebraic_variable_names = sc_constraint_algebraic_variable_names;
  model.constraint_species_dependencies = sc_constraint_species_dependencies;
  model.nonzero_constraint_jacobian_elements = sc_nonzero_constraint_jacobian_elements;
  model.constraint_state_parameter_names = no_names;
  model.finalize_constraint_setup = sc_finalize_constraint_setup;
  model.update_constraint_state_parameters = no_parameter_update;
  model.add_constraint_residual = sc_add_constraint_residual;
  model.subtract_constraint_jacobian = sc_subtract_constraint_jacobian;
  return model;
}

/* ---- StubAerosolWithSolvent ------------------------------------------------- */

void stub_aerosol_solvent_init(stub_aerosol_solvent *m, double rate_constant, double total_mass)
{
  memset(m, 0, sizeof *m);
  m->rate_constant = rate_constant;
  m->total_mass = total_mass;
  m->i_gas = m->i_aq = m->i_s = -1;
  m->gas_gas_flat = m->gas_s_flat = m->aq_gas_flat = m->aq_s_flat = m->aq_aq_flat = -1;
}

static int sv_species_used(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  (void)context;
  stub_name_set set;
  set.n = 0;
  stub_name_set_insert(&set, "A_GAS");
  stub_name_set_insert(&set, "AEROSOL.A_AQ");
  stub_name_set_insert(&set, "S");
  return stub_name_set_finish(&set, names, max);
}

static int sv_nonzero_jacobian_elements(void *context, const tmc_name_map *variables, tmc_index_pair *out, int max)
{
  (void)context;
  stub_pair_list list;
  list.n = 0;
  int i_gas = tmc_name_map_lookup(variables, "A_GAS");
  int i_aq = tmc_name_map_lookup(variables, "AEROSOL.A_AQ");
  int i_s = tmc_name_map_lookup(variables, "S");
  if (i_gas < 0 || i_aq < 0 || i_s < 0) /* MICM's at() throws */
    return 0;
  stub_pair_list_add(&list, i_gas, i_gas);
  stub_pair_list_add(&list, i_gas, i_s);
  stub_pair_list_add(&list, i_aq, i_gas);
  stub_pair_list_add(&list, i_aq, i_s);
  return stub_pair_list_finish(&list, out, max);
}

static void sv_finalize_process_setup(void *context, const tmc_name_map *parameters, const tmc_name_map *variables, const tmc_sparse *jacobian)
{
  stub_aerosol_solvent *m = context;
  (void)parameters;
  m->i_gas = tmc_name_map_lookup(variables, "A_GAS");
  m->i_aq = tmc_name_map_lookup(variables, "AEROSOL.A_AQ");
  m->i_s = tmc_name_map_lookup(variables, "S");
  m->gas_gas_flat = tmc_sparse_element(jacobian, m->i_gas, m->i_gas);
  m->gas_s_flat = tmc_sparse_element(jacobian, m->i_gas, m->i_s);
  m->aq_gas_flat = tmc_sparse_element(jacobian, m->i_aq, m->i_gas);
  m->aq_s_flat = tmc_sparse_element(jacobian, m->i_aq, m->i_s);
}

/* rate = k [A_GAS] [S] */
static void sv_add_forcing_terms(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_dense *forcing)
{
  const stub_aerosol_solvent *m = context;
  (void)state_parameters;
  const int rows = stub_dense_rows(forcing);
  for (int row = 0; row < rows; ++row)
  {
    const double rate = m->rate_constant * stub_dense_value(state_variables, row, m->i_gas) * stub_dense_value(state_variables, row, m->i_s);
    *stub_dense_at(forcing, row, m->i_gas) -= rate;
    *stub_dense_at(forcing, row, m->i_aq) += rate;
  }
}

static void sv_subtract_jacobian_terms(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_sparse *jacobian)
{
  const stub_aerosol_solvent *m = context;
  (void)state_parameters;
  if (m->gas_gas_flat < 0 || m->gas_s_flat < 0 || m->aq_gas_flat < 0 || m->aq_s_flat < 0)
    return;
  const double k = m->rate_constant;
  const int blocks = stub_sparse_blocks(jacobian);
  for (int block = 0; block < blocks; ++block)
  {
    const double gas = stub_dense_value(state_variables, block, m->i_gas);
    const double s = stub_dense_value(state_variables, block, m->i_s);
    *stub_sparse_at(jacobian, block, m->gas_gas_flat) -= -k * s;
    *stub_sparse_at(jacobian, block, m->aq_gas_flat) -= k * s;
    *stub_sparse_at(jacobian, block, m->gas_s_flat) -= -k * gas;
    *stub_sparse_at(jacobian, block, m->aq_s_flat) -= k * gas;
  }
}

static int sv_constraint_algebraic_variable_names(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  return sc_state_variable_names(context, names, max); /* { "AEROSOL.A_AQ" } */
}

static int sv_constraint_species_dependencies(void *context, char (*names)[TMC_NAME_MAX], int max)
{
  (void)context;
  stub_name_set set;
  set.n = 0;
  stub_name_set_insert(&set, "A_GAS");
  stub_name_set_insert(&set, "AEROSOL.A_AQ");
  return stub_name_set_finish(&set, names, max);
}

static int sv_nonzero_constraint_jacobian_elements(void *context, const tmc_name_map *variables, tmc_index_pair *out, int max)
{
  (void)context;
  stub_pair_list list;
  list.n = 0;
  int i_gas = tmc_name_map_lookup(variables, "A_GAS");
  int i_aq = tmc_name_map_lookup(variables, "AEROSOL.A_AQ");
  if (i_gas < 0 || i_aq < 0) /* MICM's at() throws */
    return 0;
  stub_pair_list_add(&list, i_aq, i_gas);
  stub_pair_list_add(&list, i_aq, i_aq);
  return stub_pair_list_finish(&list, out, max);
}

static void sv_finalize_constraint_setup(void *context, const tmc_name_map *parameters, const tmc_name_map *variables, const tmc_sparse *jacobian)
{
  stub_aerosol_solvent *m = context;
  (void)parameters;
  m->i_gas = tmc_name_map_lookup(variables, "A_GAS");
  m->i_aq = tmc_name_map_lookup(variables, "AEROSOL.A_AQ");
  m->aq_gas_flat = tmc_sparse_element(jacobian, m->i_aq, m->i_gas);
  m->aq_aq_flat = tmc_sparse_element(jacobian, m->i_aq, m->i_aq);
}

static void sv_add_constraint_residual(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_dense *forcing)
{
  const stub_aerosol_solvent *m = context;
  (void)state_parameters;
  const int rows = stub_dense_rows(forcing);
  for (int row = 0; row < rows; ++row)
    *stub_dense_at(forcing, row, m->i_aq) = stub_dense_value(state_variables, row, m->i_gas) + stub_dense_value(state_variables, row, m->i_aq) - m->total_mass;
}

static void sv_subtract_constraint_jacobian(void *context, const tmc_dense *state_parameters, const tmc_dense *state_variables, tmc_sparse *jacobian)
{
  const stub_aerosol_solvent *m = context;
  (void)state_parameters;
  (void)state_variables;
  if (m->aq_gas_flat < 0 || m->aq_aq_flat < 0)
    return;
  const int blocks = stub_sparse_blocks(jacobian);
  for (int block = 0; block < blocks; ++block)
  {
    *stub_sparse_at(jacobian, block, m->aq_gas_flat) -= 1.0;
    *stub_sparse_at(jacobian, block, m->aq_aq_flat) -= 1.0;
  }
}

tmc_external_model stub_aerosol_solvent_model(stub_aerosol_solvent *m)
{
  tmc_external_model model;
  memset(&model, 0, sizeof model);
  model.context = m;
  model.state_variable_names = sc_state_variable_names; /* { "AEROSOL.A_AQ" } */
  model.state_parameter_names = no_names;
  model.species_used = sv_species_used;
  model.nonzero_jacobian_elements = sv_nonzero_jacobian_elements;
  model.finalize_process_setup = sv_finalize_process_setup;
  model.update_state_parameters = no_parameter_update;
  model.add_forcing_terms = sv_add_forcing_terms;
  model.subtract_jacobian_terms = sv_subtract_jacobian_terms;
  model.constraint_algebraic_variable_names = sv_constraint_algebraic_variable_names;
  model.constraint_species_dependencies = sv_constraint_species_dependencies;
  model.nonzero_constraint_jacobian_elements = sv_nonzero_constraint_jacobian_elements;
  model.constraint_state_parameter_names = no_names;
  model.finalize_constraint_setup = sv_finalize_constraint_setup;
  model.update_constraint_state_parameters = no_parameter_update;
  model.add_constraint_residual = sv_add_constraint_residual;
  model.subtract_constraint_jacobian = sv_subtract_constraint_jacobian;
  return model;
}
