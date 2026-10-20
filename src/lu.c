#include "chemistry.h"

#include <math.h>

static int dense_size_valid(size_t dimension) {
  return dimension == 0 || dimension <= (size_t)-1 / dimension;
}

int chem_lu_decompose_doolittle(size_t dimension, const double *matrix,
                                double *lower, double *upper) {
  if (!dense_size_valid(dimension) ||
      (dimension && (!matrix || !lower || !upper))) return 1;
  const size_t count = dimension*dimension;
  for (size_t i = 0; i < count; ++i) lower[i] = upper[i] = 0.0;
  for (size_t pivot = 0; pivot < dimension; ++pivot) {
    lower[pivot*dimension+pivot] = 1.0;
    for (size_t column = pivot; column < dimension; ++column) {
      double value = matrix[pivot*dimension+column];
      for (size_t k = 0; k < pivot; ++k)
        value -= lower[pivot*dimension+k]*upper[k*dimension+column];
      upper[pivot*dimension+column] = value;
    }
    const double diagonal = upper[pivot*dimension+pivot];
    if (diagonal == 0.0 || !isfinite(diagonal)) return 1;
    for (size_t row = pivot+1; row < dimension; ++row) {
      double value = matrix[row*dimension+pivot];
      for (size_t k = 0; k < pivot; ++k)
        value -= lower[row*dimension+k]*upper[k*dimension+pivot];
      lower[row*dimension+pivot] = value/diagonal;
    }
  }
  return 0;
}

int chem_lu_decompose_mozart(size_t dimension, const double *matrix,
                             double *lower, double *upper) {
  if (!dense_size_valid(dimension) ||
      (dimension && (!matrix || !lower || !upper))) return 1;
  const size_t count = dimension*dimension;
  for (size_t i = 0; i < count; ++i) {
    lower[i] = 0.0;
    upper[i] = matrix[i];
  }
  for (size_t pivot = 0; pivot < dimension; ++pivot) {
    lower[pivot*dimension+pivot] = 1.0;
    const double diagonal = upper[pivot*dimension+pivot];
    if (diagonal == 0.0 || !isfinite(diagonal)) return 1;
    for (size_t row = pivot+1; row < dimension; ++row)
      lower[row*dimension+pivot] = upper[row*dimension+pivot]/diagonal;
    for (size_t row = pivot+1; row < dimension; ++row)
      for (size_t column = pivot+1; column < dimension; ++column)
        upper[row*dimension+column] -=
            lower[row*dimension+pivot]*upper[pivot*dimension+column];
    for (size_t row = pivot+1; row < dimension; ++row)
      upper[row*dimension+pivot] = 0.0;
  }
  return 0;
}

int chem_lu_decompose_doolittle_in_place(size_t dimension, double *matrix) {
  if (!dense_size_valid(dimension) || (dimension && !matrix)) return 1;
  for (size_t pivot = 0; pivot < dimension; ++pivot) {
    const double diagonal = matrix[pivot*dimension+pivot];
    if (diagonal == 0.0 || !isfinite(diagonal)) return 1;
    for (size_t row = pivot+1; row < dimension; ++row)
      matrix[row*dimension+pivot] /= diagonal;
    for (size_t row = pivot+1; row < dimension; ++row)
      for (size_t column = pivot+1; column < dimension; ++column)
        matrix[row*dimension+column] -=
            matrix[row*dimension+pivot]*matrix[pivot*dimension+column];
  }
  return 0;
}

int chem_lu_decompose_mozart_in_place(size_t dimension, double *matrix) {
  if (!dense_size_valid(dimension) || (dimension && !matrix)) return 1;
  for (size_t pivot = 0; pivot < dimension; ++pivot) {
    const double diagonal = matrix[pivot*dimension+pivot];
    if (diagonal == 0.0 || !isfinite(diagonal)) return 1;
    for (size_t row = pivot+1; row < dimension; ++row)
      matrix[row*dimension+pivot] /= diagonal;
    for (size_t column = pivot+1; column < dimension; ++column)
      for (size_t row = pivot+1; row < dimension; ++row)
        matrix[row*dimension+column] -=
            matrix[row*dimension+pivot]*matrix[pivot*dimension+column];
  }
  return 0;
}
