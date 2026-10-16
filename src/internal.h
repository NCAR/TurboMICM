/* SPDX-License-Identifier: Apache-2.0 */
#ifndef CH_INTERNAL_H
#define CH_INTERNAL_H
#include "chem.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

typedef struct {
  char *name;
  double atol, molecular_weight, diffusion;
  int third_body, index, has_weight, has_diffusion, used;
  ch_species_spec *spec;
} ch_species;
typedef struct { char *key; int type; union { double real; int integer; char *string; } value; } ch_property;
struct ch_species_spec {
  char *name;
  size_t nproperties;
  ch_property *properties;
  int parameterized;
  double coefficients[4];
};
struct ch_phase { char *name; size_t count; ch_species_spec **species; char **qualified; };
typedef struct { char *name; double initial; } ch_parameter;
typedef struct {
  char *name;
  int kind, row, custom_index;
  size_t nr, np;
  int *ids, *flat;
  double *coefficients, constant, enthalpy;
} ch_constraint;
typedef struct {
  size_t nr, np, nparameterized;
  int *reactants, *products, *parameterized;
  double *yields;
  ch_rate rate;
} ch_reaction;
typedef struct { int destination, left, right, kind; } ch_lu_op;
struct ch_sparse {
  size_t n, nnz, nops;
  int *index, *row_start, *columns, *diagonal;
  ch_lu_op *ops;
};
typedef struct {
  int reaction, independent;
  size_t nr, np;
  int *reactants, *destinations;
  double *yields;
} ch_derivative;
struct ch_model {
  size_t ns, nv, nr, nparam, nderiv, nconstraint, nexternal, nalgebraic;
  ch_species *species;
  ch_reaction *reactions;
  ch_parameter *parameters;
  ch_derivative *derivatives;
  ch_constraint *constraints;
  ch_external *externals;
  int *order;
  unsigned char *algebraic;
  ch_sparse *sparse;
  int finalized;
};
struct ch_state {
  const ch_model *model;
  size_t cells;
  double *y, *custom, *rates, *atol, rtol;
  ch_conditions *conditions;
  double *new_y, *forcing, *error, *stages, *jacobian, *scratch, *inverse_diagonal;
  int lu_algorithm;
  ch_factor_function factor;
  ch_linear_function linear_solve;
  void *linear_context;
};
void *ch_allocate(size_t n, size_t size);
char *ch_string_copy(const char *);
void ch_update_constraints(ch_state *);
void ch_post_solve(ch_state *);
int ch_factor_state(ch_state *);
void ch_linear_state(ch_state *, double *);
int ch_lu_factor_work(const ch_sparse *, size_t, double *, int, double *, double *inverse_diagonal);
void ch_lu_solve_work(const ch_sparse *, size_t, const double *, double *, const double *inverse_diagonal);
#endif
