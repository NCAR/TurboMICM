#include "chem.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int valid_dense_order(size_t order) {
  return order != 0u && order <= SIZE_MAX / order && order * order <= SIZE_MAX / sizeof(double);
}

static int factor_mozart(double *matrix, size_t order) {
  for (size_t pivot_index = 0; pivot_index < order; ++pivot_index) {
    const double pivot = matrix[pivot_index * order + pivot_index];
    if (pivot == 0.0 || !isfinite(pivot)) return -1;
    for (size_t row = pivot_index + 1u; row < order; ++row) {
      matrix[row * order + pivot_index] /= pivot;
      for (size_t column = pivot_index + 1u; column < order; ++column)
        matrix[row * order + column] -=
            matrix[row * order + pivot_index] * matrix[pivot_index * order + column];
    }
  }
  return 0;
}

static int factor_doolittle(double *matrix, size_t order) {
  for (size_t row = 0; row < order; ++row) {
    for (size_t column = row; column < order; ++column) {
      double value = matrix[row * order + column];
      for (size_t inner = 0; inner < row; ++inner)
        value -= matrix[row * order + inner] * matrix[inner * order + column];
      matrix[row * order + column] = value;
    }
    const double pivot = matrix[row * order + row];
    if (pivot == 0.0 || !isfinite(pivot)) return -1;
    for (size_t lower_row = row + 1u; lower_row < order; ++lower_row) {
      double value = matrix[lower_row * order + row];
      for (size_t inner = 0; inner < row; ++inner)
        value -= matrix[lower_row * order + inner] * matrix[inner * order + row];
      matrix[lower_row * order + row] = value / pivot;
    }
  }
  return 0;
}

int chem_lu_factor_dense_in_place(double *combined_lu, size_t order, chem_lu_method method) {
  if (combined_lu == NULL || !valid_dense_order(order) ||
      (method != CHEM_LU_MOZART_IN_PLACE && method != CHEM_LU_DOOLITTLE)) return -1;
  return method == CHEM_LU_DOOLITTLE
      ? factor_doolittle(combined_lu, order)
      : factor_mozart(combined_lu, order);
}

int chem_lu_factor_dense(
    const double *matrix,
    size_t order,
    double *lower,
    double *upper,
    chem_lu_method method) {
  if (matrix == NULL || lower == NULL || upper == NULL || lower == upper ||
      matrix == lower || matrix == upper || !valid_dense_order(order)) return -1;
  double *combined = malloc(order * order * sizeof(*combined));
  if (combined == NULL) return -2;
  memcpy(combined, matrix, order * order * sizeof(*combined));
  const int status = chem_lu_factor_dense_in_place(combined, order, method);
  if (status == 0) {
    for (size_t row = 0; row < order; ++row) {
      for (size_t column = 0; column < order; ++column) {
        lower[row * order + column] = row == column ? 1.0 :
            (row > column ? combined[row * order + column] : 0.0);
        upper[row * order + column] = row <= column ? combined[row * order + column] : 0.0;
      }
    }
  }
  free(combined);
  return status;
}

int chem_lu_solve_dense_in_place(const double *combined_lu, size_t order, double *right_hand_side) {
  if (combined_lu == NULL || right_hand_side == NULL || !valid_dense_order(order)) return -1;
  for (size_t row = 0; row < order; ++row) {
    for (size_t column = 0; column < row; ++column)
      right_hand_side[row] -= combined_lu[row * order + column] * right_hand_side[column];
  }
  for (size_t row = order; row-- > 0u;) {
    for (size_t column = row + 1u; column < order; ++column)
      right_hand_side[row] -= combined_lu[row * order + column] * right_hand_side[column];
    const double pivot = combined_lu[row * order + row];
    if (pivot == 0.0 || !isfinite(pivot)) return -1;
    right_hand_side[row] /= pivot;
    if (!isfinite(right_hand_side[row])) return -1;
  }
  return 0;
}
