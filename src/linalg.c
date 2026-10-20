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

static int sparse_position_map(
    size_t order, size_t nonzero_count, const size_t *rows, const size_t *columns, size_t **position_out) {
  if (!order || order > SIZE_MAX / order || order * order > SIZE_MAX / sizeof(size_t) ||
      !rows || !columns || !position_out) return -1;
  size_t *position = malloc(order * order * sizeof(*position));
  if (!position) return -2;
  for (size_t i = 0; i < order * order; ++i) position[i] = SIZE_MAX;
  for (size_t entry = 0; entry < nonzero_count; ++entry) {
    if (rows[entry] >= order || columns[entry] >= order) {
      free(position);
      return -1;
    }
    const size_t key = rows[entry] * order + columns[entry];
    if (position[key] != SIZE_MAX) {
      free(position);
      return -1;
    }
    position[key] = entry;
  }
  *position_out = position;
  return 0;
}

int chem_lu_factor_sparse_in_place(
    size_t order, size_t nonzero_count, const size_t *rows, const size_t *columns, double *values) {
  if (!values) return -1;
  size_t *position = NULL;
  const int map_status = sparse_position_map(order, nonzero_count, rows, columns, &position);
  if (map_status != 0) return map_status;

  int status = 0;
  for (size_t pivot = 0; pivot < order; ++pivot) {
    const size_t diagonal = position[pivot * order + pivot];
    if (diagonal == SIZE_MAX || values[diagonal] == 0.0 || !isfinite(values[diagonal])) {
      status = -1;
      break;
    }
    values[diagonal] = 1.0 / values[diagonal];
    if (!isfinite(values[diagonal])) {
      status = -1;
      break;
    }
    for (size_t row = pivot + 1u; row < order; ++row) {
      const size_t lower = position[row * order + pivot];
      if (lower == SIZE_MAX) continue;
      values[lower] *= values[diagonal];
      if (!isfinite(values[lower])) {
        status = -1;
        break;
      }
      for (size_t column = pivot + 1u; column < order; ++column) {
        const size_t upper = position[pivot * order + column];
        if (upper == SIZE_MAX) continue;
        const double update = values[lower] * values[upper];
        const size_t target = position[row * order + column];
        if (target == SIZE_MAX) {
          if (update != 0.0) status = -1;
          if (status != 0) break;
        } else {
          values[target] -= update;
          if (!isfinite(values[target])) {
            status = -1;
            break;
          }
        }
      }
      if (status != 0) break;
    }
    if (status != 0) break;
  }
  free(position);
  return status;
}

int chem_lu_solve_sparse_in_place(
    size_t order, size_t nonzero_count, const size_t *rows, const size_t *columns,
    const double *factor_values, double *right_hand_side) {
  if (!factor_values || !right_hand_side) return -1;
  size_t *position = NULL;
  const int map_status = sparse_position_map(order, nonzero_count, rows, columns, &position);
  if (map_status != 0) return map_status;
  for (size_t diagonal = 0; diagonal < order; ++diagonal) {
    const size_t entry = position[diagonal * order + diagonal];
    if (entry == SIZE_MAX || factor_values[entry] == 0.0 || !isfinite(factor_values[entry])) {
      free(position);
      return -1;
    }
  }
  for (size_t row = 0; row < order; ++row) {
    for (size_t column = 0; column < row; ++column) {
      const size_t entry = position[row * order + column];
      if (entry != SIZE_MAX) right_hand_side[row] -= factor_values[entry] * right_hand_side[column];
    }
  }
  for (size_t row = order; row-- > 0u;) {
    for (size_t column = row + 1u; column < order; ++column) {
      const size_t entry = position[row * order + column];
      if (entry != SIZE_MAX) right_hand_side[row] -= factor_values[entry] * right_hand_side[column];
    }
    right_hand_side[row] *= factor_values[position[row * order + row]];
    if (!isfinite(right_hand_side[row])) {
      free(position);
      return -1;
    }
  }
  free(position);
  return 0;
}

int chem_sparse_matrix_vector(
    size_t row_count, size_t column_count, size_t nonzero_count,
    const size_t *rows, const size_t *columns, const double *values,
    const double *vector, double *result) {
  if (!row_count || !column_count || !rows || !columns || !values || !vector || !result) return -1;
  memset(result, 0, row_count * sizeof(*result));
  for (size_t entry = 0; entry < nonzero_count; ++entry) {
    if (rows[entry] >= row_count || columns[entry] >= column_count) return -1;
    result[rows[entry]] += values[entry] * vector[columns[entry]];
  }
  return 0;
}
