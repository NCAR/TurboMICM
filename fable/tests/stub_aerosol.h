/* MICM's stub external models as C structs with the library's callback table:
 *
 *   stub_aerosol_1               StubAerosolModel             (test/integration/stub_aerosol_1.hpp)
 *   stub_aerosol_2               AnotherStubAerosolModel      (test/integration/stub_aerosol_2.hpp)
 *   stub_aerosol_constraints     StubAerosolWithConstraints   (test/integration/stub_aerosol_with_constraints.hpp)
 *   stub_aerosol_solvent         StubAerosolWithSolvent       (same file)
 *
 * Each struct holds the model's data and the indices it binds in its finalize
 * callbacks; the table from stub_*_model() has the struct as its context. The
 * library copies the table but borrows the context, so a model must outlive
 * the solvers built with it. A plain struct copy is an independent model, as
 * a copy of the C++ object is.
 *
 * The callbacks follow MICM's DenseMatrixPolicy::Function and
 * SparseMatrixPolicy::Function: they act on every row (block) of the matrices
 * they receive, padding rows of the last vector group included. */
#ifndef TMC_TEST_STUB_AEROSOL_H
#define TMC_TEST_STUB_AEROSOL_H

#include "tmc.h"

#include <stdbool.h>

#define STUB_MAX_PHASE_SPECIES 8

/* Rate constants used by the tests (stub_aerosol_1.hpp) */
#define STUB1_RATE_CONSTANT_FO2_CORGE 1e-3
#define STUB1_RATE_CONSTANT_BAZ_QUUX 2e-3

/* What the stubs keep of a micm::Phase: its name and UniqueNames(), the
 * "PHASE.SPECIES" names of the species that hold a state variable. */
typedef struct stub_phase
{
  char name[TMC_NAME_MAX];
  int n_unique_names;
  char unique_names[STUB_MAX_PHASE_SPECIES][TMC_NAME_MAX];
} stub_phase;

stub_phase stub_phase_from(const tmc_phase *phase);

/* One first-order transfer reactant -> product, with the Jacobian elements of
 * (reactant, reactant) and (product, reactant). */
typedef struct stub_forcing_info
{
  int reactant_index;
  int product_index;
  double rate_constant;  /* StubAerosolModel */
  int rate_param_index;  /* AnotherStubAerosolModel */
  int reactant_reactant_flat;
  int product_reactant_flat;
} stub_forcing_info;

/* ---- StubAerosolModel ------------------------------------------------------
 * A single-moment two-mode aerosol model. Mode 1 holds the first phase, mode 2
 * both phases. Two processes: FO2 gas -> FO2 mode 2 CORGE and BAZ mode 1 ->
 * BAZ mode 2 (QUUX), first order with fixed rate constants. */
typedef struct stub1_rate_constants
{
  double fo2_gas_to_mode2_corge;
  double baz_mode1_to_mode2_quux;
} stub1_rate_constants;

typedef struct stub_aerosol_1
{
  char name[TMC_NAME_MAX];
  int n_phases;
  stub_phase phases[2];
  stub1_rate_constants rate_constants;
  int n_forcing;
  stub_forcing_info forcing_info[2];
} stub_aerosol_1;

/* phases: the two phases of the model (MICM passes { quux, corge }) */
void stub_aerosol_1_init(stub_aerosol_1 *m, const char *name, const tmc_phase *phases, int n_phases, stub1_rate_constants rate_constants);
tmc_external_model stub_aerosol_1_model(stub_aerosol_1 *m);
/* StateSize(): the number of state variables and of state parameters */
void stub_aerosol_1_state_size(const stub_aerosol_1 *m, int *n_variables, int *n_parameters);
/* StateVariableNames() and StateParameterNames(), sorted and unique as a
 * std::set; names may be NULL to ask for the count. */
int stub_aerosol_1_state_variable_names(const stub_aerosol_1 *m, char (*names)[TMC_NAME_MAX], int max);
int stub_aerosol_1_state_parameter_names(const stub_aerosol_1 *m, char (*names)[TMC_NAME_MAX], int max);
/* Species(mode, phase, species): "<name>.MODE<mode + 1>.<phase>.<species>" */
void stub_aerosol_1_species(const stub_aerosol_1 *m, int mode, const tmc_phase *phase, const tmc_species *species, char out[TMC_NAME_MAX]);

/* ---- AnotherStubAerosolModel -----------------------------------------------
 * A two-moment three-mode aerosol model with temperature-dependent state
 * parameters. Processes: FO2 -> BAZ in mode 2 CORGE at the rate parameter
 * "<name>.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT" (set by the user), and
 * BAZ -> QUX in mode 3 QUUX at "<name>.PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT"
 * = 0.005 T (set by update_state_parameters). */
typedef struct stub_aerosol_2
{
  char name[TMC_NAME_MAX];
  int n_phases;
  stub_phase phases[2];
  int n_forcing;
  stub_forcing_info forcing_info[2];
  int baz_to_qux_param_index;
} stub_aerosol_2;

void stub_aerosol_2_init(stub_aerosol_2 *m, const char *name, const tmc_phase *phases, int n_phases);
tmc_external_model stub_aerosol_2_model(stub_aerosol_2 *m);
void stub_aerosol_2_state_size(const stub_aerosol_2 *m, int *n_variables, int *n_parameters);
int stub_aerosol_2_state_variable_names(const stub_aerosol_2 *m, char (*names)[TMC_NAME_MAX], int max);
int stub_aerosol_2_state_parameter_names(const stub_aerosol_2 *m, char (*names)[TMC_NAME_MAX], int max);
void stub_aerosol_2_species(const stub_aerosol_2 *m, int mode, const tmc_phase *phase, const tmc_species *species, char out[TMC_NAME_MAX]);
/* Number(mode): "<name>.MODE<mode + 1>.NUMBER" */
void stub_aerosol_2_number(const stub_aerosol_2 *m, int mode, char out[TMC_NAME_MAX]);

/* ---- StubAerosolWithConstraints --------------------------------------------
 * A_GAS -> AEROSOL.A_AQ at rate k. With a total mass, the mass-conservation
 * constraint [A_GAS] + [A_AQ] - total = 0 replaces the A_AQ row (DAE mode);
 * without one, the constraint callbacks report no algebraic variable and the
 * model contributes processes only. */
typedef struct stub_aerosol_constraints
{
  double rate_constant;
  bool has_total_mass;
  double total_mass;
  int i_gas, i_aq;
  int gas_gas_flat, aq_gas_flat, aq_aq_flat;
} stub_aerosol_constraints;

void stub_aerosol_constraints_init(stub_aerosol_constraints *m, double rate_constant);
void stub_aerosol_constraints_init_with_total(stub_aerosol_constraints *m, double rate_constant, double total_mass);
tmc_external_model stub_aerosol_constraints_model(stub_aerosol_constraints *m);

/* ---- StubAerosolWithSolvent ------------------------------------------------
 * StubAerosolWithConstraints with a solvent S: A_GAS -> A_AQ at rate
 * k [A_GAS] [S]. The process Jacobian has (A_AQ, S) in the algebraic A_AQ row,
 * an element no constraint declares. */
typedef struct stub_aerosol_solvent
{
  double rate_constant;
  double total_mass;
  int i_gas, i_aq, i_s;
  int gas_gas_flat, gas_s_flat, aq_gas_flat, aq_s_flat, aq_aq_flat;
} stub_aerosol_solvent;

void stub_aerosol_solvent_init(stub_aerosol_solvent *m, double rate_constant, double total_mass);
tmc_external_model stub_aerosol_solvent_model(stub_aerosol_solvent *m);

/* ---- helpers shared with the tests' own models ------------------------------ */

/* A std::set<std::string> of names under construction. */
#define STUB_MAX_NAMES 32
typedef struct stub_name_set
{
  int n;
  char names[STUB_MAX_NAMES][TMC_NAME_MAX];
} stub_name_set;

void stub_name_set_insert(stub_name_set *set, const char *name);
/* Sort and drop duplicates, then copy up to max names out (names may be NULL);
 * returns the set's size. */
int stub_name_set_finish(stub_name_set *set, char (*names)[TMC_NAME_MAX], int max);

/* A std::set<std::pair<Index, Index>> under construction (pairs may repeat;
 * the library's pattern has set semantics). Returns the count; out may be
 * NULL to ask for it. */
typedef struct stub_pair_list
{
  int n;
  tmc_index_pair pairs[32];
} stub_pair_list;

void stub_pair_list_add(stub_pair_list *list, int row, int col);
int stub_pair_list_finish(const stub_pair_list *list, tmc_index_pair *out, int max);

/* Rows of a dense matrix and blocks of a sparse matrix, padding included. */
int stub_dense_rows(const tmc_dense *m);
int stub_sparse_blocks(const tmc_sparse *m);
/* Element (row, col) of a dense matrix, and element e of a block of a sparse
 * matrix (e from tmc_sparse_element on the pattern a model was bound to). */
double *stub_dense_at(tmc_dense *m, int row, int col);
double stub_dense_value(const tmc_dense *m, int row, int col);
double *stub_sparse_at(tmc_sparse *m, int block, int element);

#endif
