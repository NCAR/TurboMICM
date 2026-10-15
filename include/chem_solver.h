/* SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_SOLVER_H
#define CHEM_SOLVER_H
#include "chem.h"

typedef struct ChemSparse ChemSparse;
typedef enum { CHEM_LU_MOZART, CHEM_LU_DOOLITTLE } ChemLuAlgorithm;
/* Pattern is an n*n array of bytes; nonzero means a structural entry.
 * The constructor adds the diagonal and symbolic LU fill. */
ChemStatus chem_sparse_create(size_t n, const unsigned char *pattern, ChemSparse **out);
void chem_sparse_free(ChemSparse *);
size_t chem_sparse_order(const ChemSparse *);
size_t chem_sparse_elements(const ChemSparse *);
size_t chem_sparse_index(const ChemSparse *, size_t row, size_t column);
/* Multiply unfactored sparse values by a column-major vector. Output must not
 * overlap either input; writes order*cells values, including structural zeros. */
void chem_sparse_multiply(const ChemSparse *,size_t cells,const double *values,const double *vector,double *output);
ChemStatus chem_lu_factor(const ChemSparse *, size_t cells, double *values, ChemLuAlgorithm);
/* Out-of-place factors use the same compact L/U storage; input is preserved.
 * Output must contain chem_sparse_elements(pattern)*cells doubles. */
ChemStatus chem_lu_factor_copy(const ChemSparse *,size_t cells,const double *input,double *output,ChemLuAlgorithm);
double chem_lu_lower(const ChemSparse *,size_t cells,const double *factors,size_t cell,size_t row,size_t column);
double chem_lu_upper(const ChemSparse *,size_t cells,const double *factors,size_t cell,size_t row,size_t column);
void chem_lu_solve(const ChemSparse *, size_t cells, const double *values, double *rhs);
ChemStatus chem_markowitz(size_t n, const unsigned char *pattern, size_t *permutation);
/* Reorder once before constructing any states or solvers. */
ChemStatus chem_mechanism_reorder(ChemMechanism *);

typedef struct {
    size_t stages, max_steps;
    double round_off, factor_min, factor_max, rejection_factor, safety_factor;
    double h_min, h_max, h_start, local_order;
    int new_function[6];
    double a[15], c[15], m[6], e[6], alpha[6], gamma[6];
    size_t constraint_max_iterations, constraint_max_backtracks;
    double constraint_tolerance, constraint_backtrack_factor, constraint_sufficient_decrease;
} ChemRosenbrockParameters;
typedef enum { CHEM_ROS2, CHEM_ROS3, CHEM_ROS4, CHEM_ROS4_DA, CHEM_ROS6_DA } ChemRosenbrockMethod;
ChemRosenbrockParameters chem_rosenbrock_parameters(ChemRosenbrockMethod method);
typedef struct {
    size_t internal_steps, accepted, rejected, function_calls, jacobian_updates, decompositions, solves;
    size_t constraint_init_iterations;
    double final_time;
} ChemSolverStats;

typedef double (*ChemStepController)(double error, const ChemRosenbrockParameters *, void *context);
double chem_default_step_controller(double error, const ChemRosenbrockParameters *, void *context);
typedef struct {
    ChemStatus (*factor)(const ChemSparse *, size_t, double *, void *context);
    void (*solve)(const ChemSparse *, size_t, const double *, double *, void *context);
    void *context;
} ChemLinearSolver;
typedef struct ChemSolver ChemSolver;
ChemStatus chem_solver_create(const ChemMechanism *, size_t cells, ChemSolver **out);
void chem_solver_free(ChemSolver *);
void chem_solver_set_controller(ChemSolver *, ChemStepController, void *context);
void chem_solver_set_linear_solver(ChemSolver *, ChemLinearSolver);
ChemStatus chem_solve(ChemSolver *, ChemState *, double dt,
                      const ChemRosenbrockParameters *, ChemSolverStats *);
ChemStatus chem_initialize_constraints(ChemSolver *, ChemState *,
                                       const ChemRosenbrockParameters *, ChemSolverStats *);
/* Add alpha*M to a negative sparse Jacobian; M is zero on algebraic rows. */
void chem_add_mass_diagonal(const ChemMechanism *,const ChemSparse *,size_t cells,double alpha,double *matrix);
double chem_normalized_error(const ChemState *, const double *y, const double *next, const double *error);
typedef struct {
    double small, h_start;
    size_t max_iterations;
    double time_step_reductions[5];
} ChemBackwardEulerParameters;
ChemBackwardEulerParameters chem_backward_euler_parameters(void);
int chem_backward_euler_converged(const ChemState *, const double *residual,
                                  const double *values, const ChemBackwardEulerParameters *);
ChemStatus chem_solve_backward_euler(ChemSolver *, ChemState *, double dt,
                                     const ChemBackwardEulerParameters *, ChemSolverStats *);
#endif
