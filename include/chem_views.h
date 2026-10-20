#ifndef CHEM_VIEWS_H
#define CHEM_VIEWS_H

#include <stddef.h>

typedef struct {
  double *base;
  size_t length;
  size_t stride;
} chem_group_view;

typedef struct {
  const double *base;
  size_t length;
  size_t stride;
} chem_const_group_view;

typedef double (*chem_group_unary_function)(double value, void *context);
typedef double (*chem_group_binary_function)(double left, double right, void *context);

/* Create a non-owning strided view after checking it fits in the backing
   array. Zero-length views may point at NULL. */
int chem_group_view_init(
    double *storage, size_t storage_length, size_t offset,
    size_t length, size_t stride, chem_group_view *view);
int chem_const_group_view_init(
    const double *storage, size_t storage_length, size_t offset,
    size_t length, size_t stride, chem_const_group_view *view);

int chem_group_view_fill(chem_group_view view, double value);
int chem_group_view_copy(chem_group_view destination, chem_const_group_view source);
int chem_group_view_apply(
    chem_group_view destination, chem_const_group_view source,
    chem_group_unary_function function, void *context);
int chem_group_view_apply_binary(
    chem_group_view destination, chem_const_group_view left, chem_const_group_view right,
    chem_group_binary_function function, void *context);

#endif
