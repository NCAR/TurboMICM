#ifndef CK_INTERNAL_H
#define CK_INTERNAL_H
#include "chem.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

typedef struct { size_t target, left, right; } ck_operation;
typedef struct { size_t reaction, independent, nreactants, *reactants, nproducts; ck_product *products; size_t *targets; } ck_derivative;
typedef struct { ck_constraint description; size_t row, parameter; ck_product *reactants, *products; } ck_bound_constraint;
struct ck_model {
  const ck_mechanism *description;
  size_t nconstraints;
  ck_bound_constraint *constraints;
  unsigned char *algebraic;
  size_t n, nr, nc, nnz, njobs;
  size_t *species_to_variable, *variable_to_species;
  size_t *reaction_order, *custom_offset, *nreactants, **reactants, *nproducts;
  ck_product **products;
  ck_rate *rates;
  char **custom_names;
  size_t *map, *diagonal, *row_start, *columns;
  ck_derivative *jobs;
  size_t nops;
  ck_operation *ops;
  size_t *divide_start, *divides, *update_start;
};
typedef struct { double *initial, *next, *error, *stage[6], *jacobian, *scratch; } ck_workspace;
void *ck_allocate(size_t n, size_t size);
char *ck_string_copy(const char *text);
#endif
