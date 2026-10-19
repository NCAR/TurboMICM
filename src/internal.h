/* Internal structures shared by the library's sources. */
#ifndef TMC_INTERNAL_H
#define TMC_INTERNAL_H

#include "tmc.h"

#include <stdlib.h>
#include <string.h>

/* Allocation helpers that abort on failure: the library has no way to recover
 * from exhausted memory, and callers would have to check every call. */
void *tmc_malloc(size_t bytes);
void *tmc_calloc(size_t count, size_t bytes);
void *tmc_realloc(void *p, size_t bytes);
char *tmc_strdup(const char *s);
/* The arrays the kernels stream through; one place to change their allocation. */
double *tmc_alloc_doubles(size_t count);
void tmc_free_doubles(double *p);
void tmc_strlcpy(char *dst, const char *src, size_t size);

/* Growable int / double arrays used while building plans. */
typedef struct tmc_ivec
{
  int *v;
  int n, cap;
} tmc_ivec;
typedef struct tmc_dvec
{
  double *v;
  int n, cap;
} tmc_dvec;
void tmc_ivec_push(tmc_ivec *a, int x);
void tmc_dvec_push(tmc_dvec *a, double x);
void tmc_ivec_free(tmc_ivec *a);
void tmc_dvec_free(tmc_dvec *a);

/* FNV-1a over bytes, for fingerprints shared with tools/mechanism.py. */
unsigned long long tmc_fnv1a(unsigned long long h, const void *bytes, size_t n);
/* MICM's ParameterizedMultiplier::kMaxFuncs */
#define TMC_MAX_PARAMETERIZED_REACTANTS 4
#define TMC_FNV_OFFSET 14695981039346656037ULL

/* ---- Symbolic LU ------------------------------------------------------- */

struct tmc_lu
{
  tmc_lu_algorithm algorithm;
  int n;
  tmc_sparse *pattern; /* combined L and U pattern, 1 block, indexing only */
  /* Mozart (right-looking) plan, MICM's index lists:
   *   for each pivot i: diag[i], nji[i] column elements aji[], nki[i] row elements,
   *   for each row element: aik, njk, then njk pairs (ajk, aji). */
  int *mz_diag, *mz_nji, *mz_nki; /* [n] */
  int *mz_aji;                    /* [sum nji] */
  int *mz_aik, *mz_njk;           /* [sum nki] */
  int *mz_ajk, *mz_ajk_aji;       /* [sum njk] */
  /* Doolittle plan (row-wise), MICM's lists:
   *   row i: for j < i (L part): niLU pairs... see lu.c */
  int *dl_n_lower, *dl_n_upper;    /* [n] number of L entries (j < i) and U entries (j >= i) per row */
  int *dl_elem;                    /* [nnz] element ids, row by row, L entries then U entries */
  int *dl_nterms;                  /* [nnz] number of (L[i][k], U[k][j]) terms per entry */
  int *dl_terms;                   /* [2 * total terms] element id pairs */
  int *dl_diag_of_col;             /* [nnz] for L entries: element id of U[j][j] */
  /* Solve plan: forward rows then backward rows (reverse order). */
  int *fw_count;  /* [n] L entries per row */
  int *fw_elem;   /* [sum] element ids L[i][j] */
  int *fw_col;    /* [sum] j */
  int *bw_count;  /* [n] U entries per row (strict upper), rows in reverse order */
  int *bw_elem;   /* [sum] */
  int *bw_col;    /* [sum] */
  int *bw_diag;   /* [n] diagonal element id, rows in reverse order */
  int n_lower, n_upper; /* strict lower / upper element counts */
};

/* ---- Compiled mechanism (solver) --------------------------------------- */

typedef struct tmc_constraint_plan
{
  int row;            /* the algebraic variable */
  int n_dependencies; /* reactants then products, or the terms */
  int *dependencies;  /* state ids */
  int *flat_ids;      /* LU element of (row, dependency) */
  int parameter_index; /* the equilibrium constant's state parameter, or -1 */
} tmc_constraint_plan;

struct tmc_name_map
{
  int n;
  char (*names)[TMC_NAME_MAX];
};

typedef struct tmc_process_plan
{
  int n_reactants; /* non-parameterized reactants (state variable ids) */
  int *reactants;
  int n_products;
  int *products;
  double *yields;
  int n_parameterized; /* parameterized reactants, multiplied into the rate constant */
  tmc_species parameterized[TMC_MAX_PARAMETERIZED_REACTANTS];
  int custom_index; /* first custom rate parameter, or -1 */
} tmc_process_plan;

struct tmc_solver
{
  tmc_method method;
  tmc_rosenbrock_parameters ros;
  tmc_backward_euler_parameters be;
  tmc_step_controller controller;
  void *controller_context;

  tmc_system system;
  int n_reactions;
  tmc_reaction *reactions; /* stably sorted by rate constant type */
  int n_vars;
  char (*variable_names)[TMC_NAME_MAX];
  int n_custom;
  char (*custom_labels)[TMC_NAME_MAX];
  double *absolute_tolerance; /* [n_vars] defaults */
  int width;                  /* L of this solver's states */
  bool reorder_state;

  tmc_process_plan *processes; /* [n_reactions] */
  /* Jacobian jobs: one per (independent variable, process) pair, ordered by
   * independent variable then process, as MICM's ProcessSet. */
  int n_jobs;
  int *job_process, *job_independent, *job_n_dependent, *job_n_products;
  int *job_dependent_ids; /* concatenated */
  int *job_product_ids;
  double *job_yields;
  int *job_flat_ids; /* per job: n_dependent + 1 + n_products LU element ids */
  int n_job_dependent_total, n_job_products_total, n_job_flat_total;

  tmc_sparse *jacobian_pattern; /* pre-fill, 1 block, indexing only */
  tmc_lu *lu;
  const tmc_sparse *lu_pattern; /* lu->pattern */
  unsigned long long fingerprint;
  const tmc_kernels *kernels; /* NULL when the generic kernels are used */

  /* constraints and algebraic variables */
  int n_constraints;
  tmc_constraint *constraints;
  tmc_constraint_plan *constraint_plans;
  int n_algebraic;
  int *algebraic_ids;          /* ascending */
  unsigned char *is_algebraic; /* [n_vars] */
  double *mass_diagonal;       /* [n_vars] 1 or 0 */
  /* external models */
  int n_external;
  tmc_external_model *externals;
  unsigned char *external_constraint_active; /* [n_external] */
  bool has_external_processes;
  int n_algebraic_flat;   /* external process elements in algebraic rows, cleared after evaluation */
  int *algebraic_flat_ids;
  tmc_name_map variable_map, parameter_map;
};

struct tmc_state
{
  const tmc_solver *solver;
  int n_cells, L, n_blocks, n_vars, n_reactions, n_custom;
  tmc_dense *variables;
  tmc_dense *rate_constants;
  tmc_dense *custom_rate_parameters;
  tmc_conditions *conditions; /* [n_blocks * L] */
  tmc_sparse *jacobian;       /* LU pattern, n_cells blocks */
  double *reciprocals;        /* [n_blocks * L * n_vars] */
  double relative_tolerance;
  double *absolute_tolerance; /* [n_vars] */
  /* Rosenbrock temporaries */
  tmc_dense *ynew, *initial_forcing, *yerror;
  tmc_dense *K[TMC_MAX_STAGES];
  /* Backward Euler temporaries */
  tmc_dense *yn, *yn1, *forcing;
  /* step-size memory, factorization reuse, and controller history (row.c, extrapolation.c) */
  double h_next;             /* the step the next solve starts from when warm-started (0: none) */
  double lu_h;               /* the step the factorization in `jacobian` was built for (0: none valid) */
  int lu_age;                /* accepted steps since that factorization */
  double rc_change;          /* compounded relative change of the rate constants since it */
  double err_prev, h_prev;   /* the last accepted step's error and size (0: none) */
  tmc_sparse *jacobian_copy; /* extrapolation: -J unshifted, refactored for each substep size */
  tmc_dense *extrap_prev[TMC_MAX_EXTRAPOLATION], *extrap_cur[TMC_MAX_EXTRAPOLATION]; /* the tableau's last two rows */
  tmc_sparse *extrap_lu[TMC_MAX_EXTRAPOLATION]; /* with reuse: the factorization of n_j / H I - J for each entry */
  double *extrap_recip[TMC_MAX_EXTRAPOLATION];  /* and its pivot reciprocals */
};

/* ---- Internal kernels (process.c) --------------------------------------- */

/* The complete forcing F(y) and the complete shifted-free negative Jacobian
 * -dF/dy of a state: kinetic terms (generic or generated kernels), external
 * models, and constraints. forcing and jacobian are overwritten. */
void tmc_eval_forcing(const tmc_solver *s, tmc_state *st, const tmc_dense *y, tmc_dense *forcing);
void tmc_eval_jacobian(const tmc_solver *s, tmc_state *st, const tmc_dense *y, tmc_sparse *jacobian);

/* Rate constants of every cell from conditions and custom parameters. */
void tmc_calculate_rate_constants(const tmc_solver *s, tmc_state *state);
/* Generic forcing: forcing += F(variables) over all blocks. */
void tmc_generic_forcing(const tmc_solver *s, const tmc_dense *rate_constants, const tmc_dense *variables, tmc_dense *forcing);
/* Generic negative Jacobian: jacobian -= dF/dy (added to existing values). */
void tmc_generic_jacobian(const tmc_solver *s, const tmc_dense *rate_constants, const tmc_dense *variables, tmc_sparse *jacobian);

/* ---- Rosenbrock internals (rosenbrock.c), also reached by tmc_solver_ext.h ---- */

/* alpha M - J: adds alpha times the mass-matrix diagonal to the diagonal of
 * the state's Jacobian (which holds -J), in every lane. */
void tmc_rosenbrock_alpha_minus_jacobian(const tmc_solver *s, tmc_state *st, double alpha);
/* MICM's normalized error over the state's active lanes (see rosenbrock.c). */
double tmc_rosenbrock_normalized_error(const tmc_state *st, const tmc_dense *y, const tmc_dense *ynew, const tmc_dense *errors);
/* MICM's InitializeConstraints: the Newton projection of the algebraic
 * variables with the given parameters; on failure the variables are restored.
 * stats is required. */
tmc_solver_state tmc_rosenbrock_initialize_constraints(const tmc_solver *s, tmc_state *st, const tmc_rosenbrock_parameters *p, tmc_stats *stats);

/* ---- Constraints (constraints.c) ---------------------------------------- */
void tmc_constraints_update_parameters(const tmc_solver *s, const tmc_conditions *conditions, int n_cells, tmc_dense *params);
void tmc_constraints_add_forcing(const tmc_solver *s, const tmc_dense *variables, const tmc_dense *params, tmc_dense *forcing);
void tmc_constraints_subtract_jacobian(const tmc_solver *s, const tmc_dense *variables, const tmc_dense *params, tmc_sparse *jacobian);
void tmc_set_algebraic_errors(const tmc_solver *s, tmc_dense *yerror, const tmc_dense *y, const tmc_dense *ynew);

/* ---- External models (external.c) -------------------------------------- */
void tmc_external_update_parameters(const tmc_solver *s, const tmc_conditions *conditions, int n_cells, tmc_dense *params);
void tmc_external_add_forcing(const tmc_solver *s, const tmc_dense *params, const tmc_dense *y, tmc_dense *forcing);
void tmc_external_subtract_jacobian(const tmc_solver *s, const tmc_dense *params, const tmc_dense *y, tmc_sparse *jacobian);
void tmc_external_update_constraint_parameters(const tmc_solver *s, const tmc_conditions *conditions, int n_cells, tmc_dense *params);
void tmc_external_add_constraint_residual(const tmc_solver *s, const tmc_dense *params, const tmc_dense *y, tmc_dense *forcing);
void tmc_external_subtract_constraint_jacobian(const tmc_solver *s, const tmc_dense *params, const tmc_dense *y, tmc_sparse *jacobian);
void tmc_external_initialize_constraint_parameters(const tmc_solver *s, const tmc_dense *y, tmc_dense *params);
/* Sorted unique names from a name callback (NULL callback: none). Returns the count. */
int tmc_external_names(const tmc_external_model *m, int (*cb)(void *, char (*)[TMC_NAME_MAX], int), char (**out)[TMC_NAME_MAX]);

/* ---- LU internals (lu.c) ------------------------------------------------ */
void tmc_lu_decompose_mozart_in_place(const tmc_lu *lu, tmc_sparse *ALU, double *reciprocals);
void tmc_lu_decompose_doolittle_in_place(const tmc_lu *lu, tmc_sparse *ALU, double *reciprocals);

/* ---- Pieces of the Rosenbrock step shared by the loops (rosenbrock.c) ---- */
void tmc_rosenbrock_forcing(const tmc_solver *s, tmc_state *st, const tmc_dense *y, tmc_dense *f);
/* alpha M - J(y) into the state's LU matrix, then its factorization. */
void tmc_rosenbrock_factor(const tmc_solver *s, tmc_state *st, const tmc_dense *y, double alpha, tmc_stats *stats);
void tmc_rosenbrock_linear_solve(const tmc_solver *s, tmc_state *st, tmc_dense *x, tmc_stats *stats);
/* The stages of one attempt with step H: K_i from the state's variables and
 * initial_forcing, with the matrix already factored. */
void tmc_rosenbrock_stages(const tmc_solver *s, tmc_state *st, const tmc_rosenbrock_parameters *p, double H, tmc_stats *stats);
void tmc_rosenbrock_combine_solution_and_error(tmc_dense *ynew, tmc_dense *yerror, const tmc_dense *y, int stages, const double *m, const double *e, tmc_dense *const *k);

/* ---- The ROW loop's pieces, shared with the extrapolation loop (row.c) ---- */
double tmc_row_initial_step(const tmc_rosenbrock_parameters *p, const tmc_state *st, double time_step, double h_min, double h_max);
double tmc_row_quantize(double h, double interval);
double tmc_row_step_factor(const tmc_rosenbrock_parameters *p, double error, double H, const tmc_state *st);
double tmc_row_error_norm(const tmc_state *st, const tmc_dense *y, const tmc_dense *ynew, const tmc_dense *errors, int norm);

/* ---- Solvers ----------------------------------------------------------- */
tmc_result tmc_rosenbrock_solve(const tmc_solver *s, double time_step, tmc_state *state);
tmc_result tmc_row_solve(const tmc_solver *s, double time_step, tmc_state *state);
tmc_result tmc_extrapolation_solve(const tmc_solver *s, double time_step, tmc_state *state);
tmc_result tmc_backward_euler_solve(const tmc_solver *s, double time_step, tmc_state *state);

#endif
