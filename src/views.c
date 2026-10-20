#include "chem_views.h"

#include <stdint.h>
#include <stdlib.h>

static int view_extent_valid(size_t storage_length, size_t offset, size_t length, size_t stride) {
  if (offset > storage_length) return 0;
  if (length == 0u) return 1;
  if (offset >= storage_length || (stride != 0u && length - 1u > SIZE_MAX / stride)) return 0;
  const size_t span = (length - 1u) * stride;
  return span <= SIZE_MAX - offset && offset + span < storage_length;
}

static int mutable_view_valid(chem_group_view view) {
  return (view.length == 0u || view.base != NULL) &&
      (view.stride == 0u || view.length <= 1u || view.length - 1u <= SIZE_MAX / view.stride);
}

static int const_view_valid(chem_const_group_view view) {
  return (view.length == 0u || view.base != NULL) &&
      (view.stride == 0u || view.length <= 1u || view.length - 1u <= SIZE_MAX / view.stride);
}

int chem_group_view_init(
    double *storage, size_t storage_length, size_t offset,
    size_t length, size_t stride, chem_group_view *view) {
  if (!view || (storage_length && !storage) || !view_extent_valid(storage_length, offset, length, stride)) return -1;
  view->base = storage ? storage + offset : NULL;
  view->length = length;
  view->stride = stride;
  return 0;
}

int chem_const_group_view_init(
    const double *storage, size_t storage_length, size_t offset,
    size_t length, size_t stride, chem_const_group_view *view) {
  if (!view || (storage_length && !storage) || !view_extent_valid(storage_length, offset, length, stride)) return -1;
  view->base = storage ? storage + offset : NULL;
  view->length = length;
  view->stride = stride;
  return 0;
}

int chem_group_view_fill(chem_group_view view, double value) {
  if (!mutable_view_valid(view)) return -1;
  for (size_t i = 0; i < view.length; ++i) view.base[i * view.stride] = value;
  return 0;
}

int chem_group_view_copy(chem_group_view destination, chem_const_group_view source) {
  if (!mutable_view_valid(destination) || !const_view_valid(source) || destination.length != source.length ||
      destination.length > SIZE_MAX / sizeof(double)) return -1;
  double *temporary = destination.length ? malloc(destination.length * sizeof(*temporary)) : NULL;
  if (destination.length && !temporary) return -2;
  for (size_t i = 0; i < destination.length; ++i) temporary[i] = source.base[i * source.stride];
  for (size_t i = 0; i < destination.length; ++i) destination.base[i * destination.stride] = temporary[i];
  free(temporary);
  return 0;
}

int chem_group_view_apply(
    chem_group_view destination, chem_const_group_view source,
    chem_group_unary_function function, void *context) {
  if (!mutable_view_valid(destination) || !const_view_valid(source) || !function ||
      destination.length != source.length || destination.length > SIZE_MAX / sizeof(double)) return -1;
  double *temporary = destination.length ? malloc(destination.length * sizeof(*temporary)) : NULL;
  if (destination.length && !temporary) return -2;
  for (size_t i = 0; i < destination.length; ++i)
    temporary[i] = function(source.base[i * source.stride], context);
  for (size_t i = 0; i < destination.length; ++i) destination.base[i * destination.stride] = temporary[i];
  free(temporary);
  return 0;
}

int chem_group_view_apply_binary(
    chem_group_view destination, chem_const_group_view left, chem_const_group_view right,
    chem_group_binary_function function, void *context) {
  if (!mutable_view_valid(destination) || !const_view_valid(left) || !const_view_valid(right) || !function ||
      destination.length != left.length || destination.length != right.length ||
      destination.length > SIZE_MAX / sizeof(double)) return -1;
  double *temporary = destination.length ? malloc(destination.length * sizeof(*temporary)) : NULL;
  if (destination.length && !temporary) return -2;
  for (size_t i = 0; i < destination.length; ++i)
    temporary[i] = function(left.base[i * left.stride], right.base[i * right.stride], context);
  for (size_t i = 0; i < destination.length; ++i) destination.base[i * destination.stride] = temporary[i];
  free(temporary);
  return 0;
}
