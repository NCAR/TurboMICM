#include "chemistry.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

int chem_linear_solve_in_place(size_t dimension, double *matrix,
                                double *right_hand_side) {
  if (dimension && (dimension > (size_t)-1/dimension || !matrix ||
                    !right_hand_side)) return 1;
  for (size_t pivot = 0; pivot < dimension; ++pivot) {
    size_t best = pivot;
    double magnitude = fabs(matrix[pivot*dimension+pivot]);
    for (size_t row = pivot+1; row < dimension; ++row) {
      const double candidate = fabs(matrix[row*dimension+pivot]);
      if (candidate > magnitude) {
        best = row;
        magnitude = candidate;
      }
    }
    if (magnitude == 0.0 || !isfinite(magnitude)) return 1;
    if (best != pivot) {
      for (size_t column = 0; column < dimension; ++column) {
        const double value = matrix[pivot*dimension+column];
        matrix[pivot*dimension+column] = matrix[best*dimension+column];
        matrix[best*dimension+column] = value;
      }
      const double value = right_hand_side[pivot];
      right_hand_side[pivot] = right_hand_side[best];
      right_hand_side[best] = value;
    }
    for (size_t row = pivot+1; row < dimension; ++row) {
      const double factor = matrix[row*dimension+pivot]/
                            matrix[pivot*dimension+pivot];
      matrix[row*dimension+pivot] = 0.0;
      for (size_t column = pivot+1; column < dimension; ++column)
        matrix[row*dimension+column] -=
            factor*matrix[pivot*dimension+column];
      right_hand_side[row] -= factor*right_hand_side[pivot];
    }
  }
  for (size_t row = dimension; row-- > 0;) {
    double value = right_hand_side[row];
    for (size_t column = row+1; column < dimension; ++column)
      value -= matrix[row*dimension+column]*right_hand_side[column];
    const double diagonal = matrix[row*dimension+row];
    if (diagonal == 0.0 || !isfinite(diagonal)) return 1;
    right_hand_side[row] = value/diagonal;
  }
  return 0;
}

int chem_linear_solve(size_t dimension, const double *matrix,
                       const double *right_hand_side, double *solution) {
  if (dimension == 0) return 0;
  if (!matrix || !right_hand_side || !solution ||
      dimension > (size_t)-1/dimension ||
      dimension*dimension > (size_t)-1/sizeof(double) ||
      dimension > (size_t)-1/sizeof(double)) return 1;
  double *matrix_copy = malloc(dimension*dimension*sizeof(double));
  double *rhs_copy = malloc(dimension*sizeof(double));
  if (!matrix_copy || !rhs_copy) {
    free(matrix_copy);
    free(rhs_copy);
    return 1;
  }
  memcpy(matrix_copy,matrix,dimension*dimension*sizeof(double));
  memcpy(rhs_copy,right_hand_side,dimension*sizeof(double));
  const int status = chem_linear_solve_in_place(dimension,matrix_copy,rhs_copy);
  if (!status) memcpy(solution,rhs_copy,dimension*sizeof(double));
  free(matrix_copy);
  free(rhs_copy);
  return status;
}
