#include "chemistry.h"

#include <math.h>
#include <stdint.h>
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

int chem_markowitz_reorder(size_t dimension,
                           const unsigned char *nonzero_pattern,
                           size_t *new_to_old_permutation) {
  if (dimension == 0) return 0;
  if (!nonzero_pattern || !new_to_old_permutation ||
      dimension > (size_t)-1/dimension) return 1;
  const size_t count = dimension*dimension;
  unsigned char *graph = malloc(count);
  unsigned char *alive = malloc(dimension);
  size_t *permutation = malloc(dimension*sizeof(size_t));
  if (!graph || !alive || !permutation) {
    free(graph); free(alive); free(permutation);
    return 1;
  }
  memcpy(graph,nonzero_pattern,count);
  memset(alive,1,dimension);
  for (size_t i = 0; i < dimension; ++i)
    if (!graph[i*dimension+i]) {
      free(graph); free(alive); free(permutation);
      return 1;
    }
  for (size_t step = 0; step < dimension; ++step) {
    size_t pivot = dimension;
    size_t best_cost = SIZE_MAX;
    for (size_t candidate = 0; candidate < dimension; ++candidate) {
      if (!alive[candidate]) continue;
      size_t row_degree = 0, column_degree = 0;
      for (size_t other = 0; other < dimension; ++other) {
        if (!alive[other]) continue;
        row_degree += graph[candidate*dimension+other] != 0;
        column_degree += graph[other*dimension+candidate] != 0;
      }
      if (!row_degree || !column_degree) {
        free(graph); free(alive); free(permutation);
        return 1;
      }
      const size_t row_cost = row_degree-1, column_cost = column_degree-1;
      const size_t cost = column_cost && row_cost > SIZE_MAX/column_cost
          ? SIZE_MAX : row_cost*column_cost;
      if (pivot == dimension || cost < best_cost) {
        pivot = candidate;
        best_cost = cost;
      }
    }
    if (pivot == dimension) {
      free(graph); free(alive); free(permutation);
      return 1;
    }
    permutation[step] = pivot;
    for (size_t row = 0; row < dimension; ++row) {
      if (!alive[row] || row == pivot || !graph[row*dimension+pivot]) continue;
      for (size_t column = 0; column < dimension; ++column)
        if (alive[column] && column != pivot &&
            graph[pivot*dimension+column])
          graph[row*dimension+column] = 1;
    }
    alive[pivot] = 0;
  }
  memcpy(new_to_old_permutation,permutation,dimension*sizeof(size_t));
  free(graph); free(alive); free(permutation);
  return 0;
}
