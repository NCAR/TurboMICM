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
  int third_body, index, has_weight, has_diffusion;
} ch_species;
typedef struct { char *name; double initial; } ch_parameter;
typedef struct {
  size_t nr, np, nparameterized;
  int *reactants, *products;
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
  size_t ns, nv, nr, nparam, nderiv;
  ch_species *species;
  ch_reaction *reactions;
  ch_parameter *parameters;
  ch_derivative *derivatives;
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
  double *new_y, *forcing, *error, *stages, *jacobian, *scratch;
};
void *ch_allocate(size_t n, size_t size);
char *ch_string_copy(const char *);
#endif
