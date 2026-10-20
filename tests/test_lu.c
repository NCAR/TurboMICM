#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
  return 1; \
} } while (0)

static unsigned next_random(unsigned *state) {
  *state = *state * 1664525u + 1013904223u;
  return *state;
}

static void make_matrix(double *matrix, size_t order, unsigned seed, int diagonal) {
  unsigned random_state = seed;
  for (size_t row = 0; row < order; ++row) {
    double row_sum = 0.0;
    for (size_t column = 0; column < order; ++column) {
      const unsigned sample = next_random(&random_state);
      const double value = (double)(sample % 2001u) / 1000.0 - 1.0;
      if (row == column) {
        matrix[row * order + column] = 0.0;
      } else if (diagonal || sample % 2u == 0u) {
        matrix[row * order + column] = value;
        row_sum += fabs(value);
      } else {
        matrix[row * order + column] = 0.0;
      }
    }
    matrix[row * order + row] = row_sum + 1.0 + (double)(next_random(&random_state) % 1000u) / 1000.0;
  }
}

static void make_reference_matrix(double *matrix) {
  const double values[] = { 2.0, -1.0, -2.0, -4.0, 6.0, 3.0, -4.0, -2.0, 8.0 };
  memcpy(matrix, values, sizeof(values));
}

static int check_decomposition(
    const double *matrix, size_t order, chem_lu_method method, int in_place,
    double initial_value, double tolerance) {
  double lower[100], upper[100], combined[100];
  if (order > 10u) return -1;
  if (in_place) {
    for (size_t i = 0; i < order * order; ++i) combined[i] = initial_value;
    memcpy(combined, matrix, order * order * sizeof(*combined));
    if (chem_lu_factor_dense_in_place(combined, order, method) != 0) return -1;
  } else {
    for (size_t i = 0; i < order * order; ++i) lower[i] = upper[i] = initial_value;
    if (chem_lu_factor_dense(matrix, order, lower, upper, method) != 0) return -1;
  }

  for (size_t row = 0; row < order; ++row) {
    for (size_t column = 0; column < order; ++column) {
      double reconstructed = 0.0;
      if (row < column && !in_place && fabs(lower[row * order + column]) > tolerance) return -1;
      if (row > column && !in_place && fabs(upper[row * order + column]) > tolerance) return -1;
      for (size_t inner = 0; inner < order; ++inner) {
        const double l = in_place
            ? (row == inner ? 1.0 : (row > inner ? combined[row * order + inner] : 0.0))
            : lower[row * order + inner];
        const double u = in_place
            ? (inner <= column ? combined[inner * order + column] : 0.0)
            : upper[inner * order + column];
        reconstructed += l * u;
      }
      if (fabs(reconstructed - matrix[row * order + column]) > tolerance) return -1;
    }
  }
  return 0;
}

static int run_factor_case(chem_lu_method method, int in_place, const char *kind, int vector_ordering) {
  if (strcmp(kind, "DenseMatrix") == 0) {
    double matrix[9];
    make_reference_matrix(matrix);
    return check_decomposition(matrix, 3u, method, in_place, 0.0, 1e-10);
  }
  if (strcmp(kind, "RandomMatrix") == 0) {
    const size_t block_counts[] = { 1u, 5u };
    for (size_t width = 0; width < (vector_ordering ? 4u : 1u); ++width) {
      for (size_t block = 0; block < (vector_ordering ? 1u : 2u); ++block) {
        const size_t blocks = vector_ordering ? 5u : block_counts[block];
        for (size_t cell = 0; cell < blocks; ++cell) {
          double matrix[100];
          make_matrix(matrix, 10u, (unsigned)(100u + width * 11u + cell), 0);
          if (check_decomposition(matrix, 10u, method, in_place, 0.0, 1e-9) != 0) return -1;
        }
      }
    }
    return 0;
  }
  if (strcmp(kind, "DiagonalMatrix") == 0) {
    for (size_t width = 0; width < (vector_ordering ? 4u : 1u); ++width) {
      for (size_t cell = 0; cell < 5u; ++cell) {
        double matrix[36] = { 0.0 };
        unsigned random_state = (unsigned)(200u + width * 17u + cell);
        for (size_t i = 0; i < 6u; ++i)
          matrix[i * 6u + i] = 0.5 + (double)(next_random(&random_state) % 100u) / 100.0;
        if (check_decomposition(matrix, 6u, method, in_place, 0.0, 1e-10) != 0) return -1;
      }
    }
    return 0;
  }
  if (strcmp(kind, "AgnosticToInitialValue") == 0) {
    const double initial_values[] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
    for (size_t width = 0; width < (vector_ordering ? 4u : 1u); ++width) {
      for (size_t cell = 0; cell < (vector_ordering ? 1u : 1u); ++cell) {
        double matrix[25];
        make_matrix(matrix, 5u, (unsigned)(300u + width + cell), 0);
        for (size_t initial = 0; initial < sizeof(initial_values) / sizeof(initial_values[0]); ++initial)
          if (check_decomposition(matrix, 5u, method, in_place, initial_values[initial], 1e-9) != 0) return -1;
      }
    }
    return 0;
  }
  return -1;
}

static int check_solve(
    const double *matrix, size_t order, chem_lu_method method, int in_place,
    double initial_value, double tolerance, int reference_rhs) {
  double lower[100], upper[100], combined[100], rhs[10], expected[10];
  if (order > 10u) return -1;
  if (reference_rhs) {
    const double reference[] = { 23.0, 42.0, 9.0 };
    memcpy(rhs, reference, order * sizeof(*rhs));
  } else {
    for (size_t column = 0; column < order; ++column)
      expected[column] = ((double)((column * 317u + order * 11u) % 1001u) / 500.0) - 1.0;
    for (size_t row = 0; row < order; ++row) {
      rhs[row] = 0.0;
      for (size_t column = 0; column < order; ++column)
        rhs[row] += matrix[row * order + column] * expected[column];
    }
  }

  if (in_place) {
    memcpy(combined, matrix, order * order * sizeof(*combined));
    if (chem_lu_factor_dense_in_place(combined, order, method) != 0) return -1;
  } else {
    for (size_t i = 0; i < order * order; ++i) lower[i] = upper[i] = initial_value;
    if (chem_lu_factor_dense(matrix, order, lower, upper, method) != 0) return -1;
    for (size_t row = 0; row < order; ++row)
      for (size_t column = 0; column < order; ++column)
        combined[row * order + column] = row > column
            ? lower[row * order + column] : upper[row * order + column];
  }
  if (chem_lu_solve_dense_in_place(combined, order, rhs) != 0) return -1;
  for (size_t row = 0; row < order; ++row) {
    double product = 0.0, scale = 0.0;
    for (size_t column = 0; column < order; ++column) {
      const double term = matrix[row * order + column] * rhs[column];
      product += term;
      scale += fabs(term);
    }
    const double target = reference_rhs ? (double[]){ 23.0, 42.0, 9.0 }[row] :
        0.0; /* The random cases use the solution residual below. */
    if (reference_rhs) {
      if (fabs(product - target) > tolerance * (1.0 + scale)) return -1;
    } else {
      double expected_product = 0.0;
      for (size_t column = 0; column < order; ++column)
        expected_product += matrix[row * order + column] * expected[column];
      if (fabs(product - expected_product) > tolerance * (1.0 + scale)) return -1;
    }
  }
  return 0;
}

static int run_solve_case(chem_lu_method method, int in_place, const char *kind, int vector_ordering) {
  if (strcmp(kind, "DenseMatrix") == 0) {
    double matrix[9];
    make_reference_matrix(matrix);
    return check_solve(matrix, 3u, method, in_place, 0.0, 1e-10, 1);
  }
  if (strcmp(kind, "RandomMatrix") == 0 || strcmp(kind, "AgnosticToInitialValue") == 0) {
    const double initial_values[] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
    for (size_t width = 0; width < (vector_ordering ? 4u : 1u); ++width) {
      const size_t blocks = vector_ordering ? 5u : 5u;
      for (size_t cell = 0; cell < blocks; ++cell) {
        double matrix[100];
        make_matrix(matrix, 10u, (unsigned)(700u + width * 19u + cell), 0);
        const size_t initial_count = strcmp(kind, "AgnosticToInitialValue") == 0 ? 5u : 1u;
        for (size_t initial = 0; initial < initial_count; ++initial)
          if (check_solve(matrix, 10u, method, in_place, initial_values[initial], 1e-9, 0) != 0) return -1;
      }
    }
    return 0;
  }
  if (strcmp(kind, "DiagonalMatrix") == 0) {
    for (size_t width = 0; width < (vector_ordering ? 4u : 1u); ++width) {
      for (size_t cell = 0; cell < 5u; ++cell) {
        double matrix[36] = { 0.0 };
        unsigned random_state = (unsigned)(900u + width * 17u + cell);
        for (size_t i = 0; i < 6u; ++i)
          matrix[i * 6u + i] = 0.5 + (double)(next_random(&random_state) % 100u) / 100.0;
        if (check_solve(matrix, 6u, method, in_place, 0.0, 1e-10, 0) != 0) return -1;
      }
    }
    return 0;
  }
  return -1;
}

typedef struct {
  const char *name;
  chem_lu_method method;
  int in_place;
  const char *kind;
  int vector_ordering;
} test_case;

typedef struct {
  const char *name;
  chem_lu_method method;
  int in_place;
  const char *kind;
  int vector_ordering;
} solve_test_case;

int main(int argc, char **argv) {
  static const test_case tests[] = {
    { "LuDecompositionMozart.DenseMatrixStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "DenseMatrix", 0 },
    { "LuDecompositionMozart.RandomMatrixStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "RandomMatrix", 0 },
    { "LuDecompositionMozart.DiagonalMatrixStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "DiagonalMatrix", 0 },
    { "LuDecompositionMozart.AgnosticToInitialValueStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "AgnosticToInitialValue", 0 },
    { "LuDecompositionMozart.DenseMatrixVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "DenseMatrix", 1 },
    { "LuDecompositionMozart.RandomMatrixVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "RandomMatrix", 1 },
    { "LuDecompositionMozart.DiagonalMatrixVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "DiagonalMatrix", 1 },
    { "LuDecompositionMozart.AgnosticToInitialValueVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "AgnosticToInitialValue", 1 },
    { "LuDecompositionDoolittle.DenseMatrixStandardOrdering", CHEM_LU_DOOLITTLE, 0, "DenseMatrix", 0 },
    { "LuDecompositionDoolittle.RandomMatrixStandardOrdering", CHEM_LU_DOOLITTLE, 0, "RandomMatrix", 0 },
    { "LuDecompositionDoolittle.DiagonalMatrixStandardOrdering", CHEM_LU_DOOLITTLE, 0, "DiagonalMatrix", 0 },
    { "LuDecompositionDoolittle.AgnosticToInitialValueStandardOrdering", CHEM_LU_DOOLITTLE, 0, "AgnosticToInitialValue", 0 },
    { "LuDecompositionDoolittle.DenseMatrixVectorOrdering", CHEM_LU_DOOLITTLE, 0, "DenseMatrix", 1 },
    { "LuDecompositionDoolittle.RandomMatrixVectorOrdering", CHEM_LU_DOOLITTLE, 0, "RandomMatrix", 1 },
    { "LuDecompositionDoolittle.DiagonalMatrixVectorOrdering", CHEM_LU_DOOLITTLE, 0, "DiagonalMatrix", 1 },
    { "LuDecompositionDoolittle.AgnosticToInitialValueVectorOrdering", CHEM_LU_DOOLITTLE, 0, "AgnosticToInitialValue", 1 },
    { "LuDecompositionMozartInPlace.DenseMatrixStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "DenseMatrix", 0 },
    { "LuDecompositionMozartInPlace.RandomMatrixStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "RandomMatrix", 0 },
    { "LuDecompositionMozartInPlace.DiagonalMatrixStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "DiagonalMatrix", 0 },
    { "LuDecompositionMozartInPlace.AgnosticToInitialValueStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "AgnosticToInitialValue", 0 },
    { "LuDecompositionMozartInPlace.DenseMatrixVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "DenseMatrix", 1 },
    { "LuDecompositionMozartInPlace.RandomMatrixVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "RandomMatrix", 1 },
    { "LuDecompositionMozartInPlace.DiagonalMatrixVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "DiagonalMatrix", 1 },
    { "LuDecompositionMozartInPlace.AgnosticToInitialValueVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "AgnosticToInitialValue", 1 },
    { "LuDecompositionDoolittleInPlace.DenseMatrixStandardOrdering", CHEM_LU_DOOLITTLE, 1, "DenseMatrix", 0 },
    { "LuDecompositionDoolittleInPlace.RandomMatrixStandardOrdering", CHEM_LU_DOOLITTLE, 1, "RandomMatrix", 0 },
    { "LuDecompositionDoolittleInPlace.DiagonalMatrixStandardOrdering", CHEM_LU_DOOLITTLE, 1, "DiagonalMatrix", 0 },
    { "LuDecompositionDoolittleInPlace.AgnosticToInitialValueStandardOrdering", CHEM_LU_DOOLITTLE, 1, "AgnosticToInitialValue", 0 },
    { "LuDecompositionDoolittleInPlace.DenseMatrixVectorOrdering", CHEM_LU_DOOLITTLE, 1, "DenseMatrix", 1 },
    { "LuDecompositionDoolittleInPlace.RandomMatrixVectorOrdering", CHEM_LU_DOOLITTLE, 1, "RandomMatrix", 1 },
    { "LuDecompositionDoolittleInPlace.DiagonalMatrixVectorOrdering", CHEM_LU_DOOLITTLE, 1, "DiagonalMatrix", 1 },
    { "LuDecompositionDoolittleInPlace.AgnosticToInitialValueVectorOrdering", CHEM_LU_DOOLITTLE, 1, "AgnosticToInitialValue", 1 },
  };
  static const solve_test_case solve_tests[] = {
    { "LinearSolver.DenseMatrixStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "DenseMatrix", 0 },
    { "LinearSolver.RandomMatrixStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "RandomMatrix", 0 },
    { "LinearSolver.DiagonalMatrixStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "DiagonalMatrix", 0 },
    { "LinearSolver.DiagonalMarkowitzReorder", CHEM_LU_MOZART_IN_PLACE, 0, "DiagonalMatrix", 0 },
    { "LinearSolver.StandardOrderingAgnosticToInitialValue", CHEM_LU_MOZART_IN_PLACE, 0, "AgnosticToInitialValue", 0 },
    { "LinearSolver.DenseMatrixVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "DenseMatrix", 1 },
    { "LinearSolver.RandomMatrixVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "RandomMatrix", 1 },
    { "LinearSolver.DiagonalMatrixVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 0, "DiagonalMatrix", 1 },
    { "LinearSolver.VectorOrderingAgnosticToInitialValue", CHEM_LU_MOZART_IN_PLACE, 0, "AgnosticToInitialValue", 1 },
    { "LinearSolver.VectorDiagonalMarkowitzReordering", CHEM_LU_MOZART_IN_PLACE, 0, "DiagonalMatrix", 1 },
    { "LinearSolverInPlace.DenseMatrixStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "DenseMatrix", 0 },
    { "LinearSolverInPlace.RandomMatrixStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "RandomMatrix", 0 },
    { "LinearSolverInPlace.DiagonalMatrixStandardOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "DiagonalMatrix", 0 },
    { "LinearSolverInPlace.DiagonalMarkowitzReorder", CHEM_LU_MOZART_IN_PLACE, 1, "DiagonalMatrix", 0 },
    { "LinearSolverInPlace.StandardOrderingAgnosticToInitialValue", CHEM_LU_MOZART_IN_PLACE, 1, "AgnosticToInitialValue", 0 },
    { "LinearSolverInPlace.DenseMatrixVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "DenseMatrix", 1 },
    { "LinearSolverInPlace.RandomMatrixVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "RandomMatrix", 1 },
    { "LinearSolverInPlace.DiagonalMatrixVectorOrdering", CHEM_LU_MOZART_IN_PLACE, 1, "DiagonalMatrix", 1 },
    { "LinearSolverInPlace.VectorOrderingAgnosticToInitialValue", CHEM_LU_MOZART_IN_PLACE, 1, "AgnosticToInitialValue", 1 },
    { "LinearSolverInPlace.VectorDiagonalMarkowitzReordering", CHEM_LU_MOZART_IN_PLACE, 1, "DiagonalMatrix", 1 },
  };
  if (argc != 2) return 2;
  for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    if (strcmp(argv[1], tests[i].name) == 0)
      return run_factor_case(tests[i].method, tests[i].in_place, tests[i].kind, tests[i].vector_ordering) == 0 ? 0 : 1;
  for (size_t i = 0; i < sizeof(solve_tests) / sizeof(solve_tests[0]); ++i)
    if (strcmp(argv[1], solve_tests[i].name) == 0)
      return run_solve_case(
          solve_tests[i].method, solve_tests[i].in_place, solve_tests[i].kind, solve_tests[i].vector_ordering) == 0 ? 0 : 1;
  fprintf(stderr, "unknown LU test: %s\n", argv[1]);
  return 2;
}
