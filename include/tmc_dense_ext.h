/*
 * tmc: the rest of MICM's dense container interface.
 *
 * MICM's Matrix is tmc_dense with L = 1 and its VectorMatrix<L> is tmc_dense
 * with width L (tmc.h). The solver needs only the operations in tmc.h; this
 * header adds the parts of MICM's container interface that it does not use:
 * construction from rows, checked row assignment, group geometry, element-wise
 * ForEach, printing, MICM's PaddedVector, and the Function facility, which
 * applies a kernel to dense matrices, sparse matrices, and vectors one group
 * at a time.
 *
 * Lanes and groups. A lane is one row of a dense matrix, one block of a sparse
 * matrix, or one value of a vector: one grid cell of the solver's arrays.
 * Lanes are stored in groups of L, the lane index fastest, so a group's L
 * values of one column (dense), one non-zero element (sparse), or one vector
 * are contiguous. When the lane count is not a multiple of L, the last group
 * ends in padding lanes. As in MICM, the operations that work on a whole group
 * (for_each, reduce, fill, copy) also touch its padding lanes, which exist in
 * storage but belong to no row; the _strict variants visit only the real
 * lanes. A padding lane holds the value the array was created with (0 for
 * tmc_padded_vector_from_values) until something writes it.
 */
#ifndef TMC_DENSE_EXT_H
#define TMC_DENSE_EXT_H

#include "tmc.h"

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------- */
/* Dense matrices                                                            */
/* ------------------------------------------------------------------------- */

/* A matrix from n_rows rows of values (MICM's constructor from a vector of
 * vectors): row i is rows[i][0 .. row_lengths[i]). The column count is
 * row_lengths[0], or 0 when there are no rows (rows and row_lengths may then
 * be NULL). Returns NULL and sets err (TMC_ERROR_MATRIX,
 * TMC_MATRIX_INVALID_VECTOR) when a row's length differs from the first. */
tmc_dense *tmc_dense_create_from_rows(int n_rows, const double *const *rows, const int *row_lengths, int L, tmc_error *err);

/* Sets row `row` from the first n_cols of n_values values (MICM's row
 * assignment from a vector; extra values are ignored). Fails with
 * TMC_MATRIX_ROW_SIZE_MISMATCH when n_values < n_cols and with
 * TMC_MATRIX_ELEMENT_OUT_OF_RANGE when the row does not exist; the matrix is
 * then unchanged. */
int tmc_dense_assign_row(tmc_dense *m, int row, const double *values, int n_values, tmc_error *err);

/* Group geometry, as MICM's VectorMatrix reports it (NumberOfGroups,
 * GroupSize, GroupVectorSize): the number of groups of L rows (a partial last
 * group included), the values per group (L * n_cols), and the width L. */
int tmc_dense_number_of_groups(const tmc_dense *m);
int tmc_dense_group_size(const tmc_dense *m);
int tmc_dense_group_width(const tmc_dense *m);

/* Element-wise ForEach (MICM's Matrix::ForEach with one or two other
 * matrices): calls f once per element of m, with pointers to the element and
 * to the elements at the same position of a (and b). Padding lanes are
 * skipped. All matrices must have the same rows, columns, and L. */
typedef void (*tmc_dense_unary_fn)(double *y, const double *a, void *ctx);
typedef void (*tmc_dense_binary_fn)(double *y, const double *a, const double *b, void *ctx);
void tmc_dense_for_each(tmc_dense *m, tmc_dense_unary_fn f, void *ctx, const tmc_dense *a);
void tmc_dense_for_each2(tmc_dense *m, tmc_dense_binary_fn f, void *ctx, const tmc_dense *a, const tmc_dense *b);

/* MICM's stream output: each row's values formatted as %g (the C++ stream
 * default), separated by commas, each row ending in a newline.
 * tmc_dense_format writes at most size - 1 characters and a terminating NUL
 * (nothing when size is 0) and returns the length of the whole text, as
 * snprintf does. tmc_dense_print writes to a stream and returns the number of
 * characters written, or a negative value on an output error. */
size_t tmc_dense_format(const tmc_dense *m, char *buffer, size_t size);
int tmc_dense_print(const tmc_dense *m, FILE *stream);

/* ------------------------------------------------------------------------- */
/* Padded vectors                                                            */
/* ------------------------------------------------------------------------- */

/* MICM's PaddedVector<T, L>: `size` values stored in padded_size =
 * ceil(size / L) * L doubles, so a vector lines up with the groups of a
 * width-L matrix of `size` rows. */
typedef struct tmc_padded_vector
{
  int size;        /* values */
  int L;           /* group width */
  int padded_size; /* stored values, padding included */
  double *data;
} tmc_padded_vector;

/* Every stored value, padding included, set to initial_value. */
tmc_padded_vector *tmc_padded_vector_create(int size, int L, double initial_value);
/* The values copied in; the padding set to 0. `values` may be NULL when size is 0. */
tmc_padded_vector *tmc_padded_vector_from_values(const double *values, int size, int L);
void tmc_padded_vector_free(tmc_padded_vector *v);

/* ------------------------------------------------------------------------- */
/* Functions over arrays (MICM's Matrix::Function and SparseMatrix::Function) */
/* ------------------------------------------------------------------------- */

/* One argument of an array function: a dense matrix, a sparse matrix, or a
 * padded vector. The const constructors mark the argument read-only; C cannot
 * enforce it, so a kernel must not write a read-only argument's lanes. */
typedef enum tmc_array_kind
{
  TMC_ARRAY_DENSE,  /* lanes are rows; a kernel reads columns */
  TMC_ARRAY_SPARSE, /* lanes are blocks; a kernel reads non-zero elements */
  TMC_ARRAY_VECTOR  /* lanes are values */
} tmc_array_kind;

typedef struct tmc_array_arg
{
  tmc_array_kind kind;
  bool read_only;
  tmc_dense *dense;          /* TMC_ARRAY_DENSE */
  tmc_sparse *sparse;        /* TMC_ARRAY_SPARSE */
  tmc_padded_vector *vector; /* TMC_ARRAY_VECTOR */
} tmc_array_arg;

tmc_array_arg tmc_array_dense(tmc_dense *m);
tmc_array_arg tmc_array_const_dense(const tmc_dense *m);
tmc_array_arg tmc_array_sparse(tmc_sparse *m);
tmc_array_arg tmc_array_const_sparse(const tmc_sparse *m);
tmc_array_arg tmc_array_vector(tmc_padded_vector *v);
tmc_array_arg tmc_array_const_vector(const tmc_padded_vector *v);

/* A kernel runs once per group, with the group (for reaching its arguments'
 * lanes) and the context given when the function was created. */
typedef struct tmc_array_group tmc_array_group;
typedef void (*tmc_array_kernel)(tmc_array_group *g, void *ctx);
typedef struct tmc_array_function tmc_array_function;

/* Creates a function of n_args arguments shaped like `args`. It records each
 * argument's kind and each dense argument's column count; lane counts may
 * differ here and are checked when the function is applied. Every argument
 * must have the same width L (a dense or sparse matrix's L, a vector's L);
 * otherwise this returns NULL and sets err (TMC_ERROR_MATRIX,
 * TMC_MATRIX_INVALID_VECTOR). ctx is passed to every kernel call. */
tmc_array_function *tmc_array_function_create(tmc_array_kernel kernel, void *ctx, int n_args, const tmc_array_arg *args, tmc_error *err);

/* Applies the function to arguments of the recorded kinds and width: the
 * kernel runs on every group of lanes in order, the last group possibly
 * partial. Before any kernel runs, the arguments are checked as MICM does when
 * a function is invoked: every argument must have the first argument's lane
 * count (a dense matrix's rows, a sparse matrix's blocks, a vector's size), and
 * every dense argument the column count recorded at creation; a mismatch
 * fails with TMC_ERROR_MATRIX, TMC_MATRIX_INVALID_VECTOR. An error raised
 * inside the kernel (see the accessors below) stops the application after
 * that kernel call and is returned. Returns TMC_OK on success. Applying does
 * not modify the function (row variables live for one application), so
 * threads may apply one function to different arguments at once when the
 * kernel's writes through its context are safe to share. */
int tmc_array_function_apply(const tmc_array_function *f, int n_args, const tmc_array_arg *args, tmc_error *err);

/* As tmc_array_function_apply, but runs the kernel on one group only (MICM's
 * GroupView of one group of a matrix); fails with
 * TMC_MATRIX_ELEMENT_OUT_OF_RANGE when the group does not exist. */
int tmc_array_function_apply_group(const tmc_array_function *f, int n_args, const tmc_array_arg *args, int group, tmc_error *err);

void tmc_array_function_free(tmc_array_function *f);

/* ---- Inside a kernel ------------------------------------------------------
 *
 * The accessors return a pointer to the L lanes of this group: lane i is
 * lanes[i], and the lanes from tmc_array_group_cells(g) on are padding. A
 * misuse (an argument index out of range, an argument of another kind, a
 * column or element out of range) returns NULL and records an error, which
 * tmc_array_function_apply returns once the kernel finishes; the helpers that
 * take lanes do nothing once an error is recorded, so a kernel need not check
 * each pointer. */

/* The group's index, its width L, and how many of its lanes are real cells
 * (L, or fewer in the last group). */
int tmc_array_group_index(const tmc_array_group *g);
int tmc_array_group_width(const tmc_array_group *g);
int tmc_array_group_cells(const tmc_array_group *g);

/* Column `column` of dense argument `arg` (MICM's GetColumnView and
 * GetConstColumnView on a GroupView): TMC_MATRIX_ELEMENT_OUT_OF_RANGE when the
 * column does not exist. */
double *tmc_array_group_column(tmc_array_group *g, int arg, int column);
/* Non-zero element `element` (an element id, as tmc_sparse_element returns) of
 * sparse argument `arg` (MICM's GetBlockView): TMC_MATRIX_ZERO_ELEMENT_ACCESS
 * for a negative id (a structural zero), TMC_MATRIX_ELEMENT_OUT_OF_RANGE for an
 * id past the last element. */
double *tmc_array_group_element(tmc_array_group *g, int arg, int element);
/* Vector argument `arg`. */
double *tmc_array_group_vector(tmc_array_group *g, int arg);
/* A new row variable (MICM's GetRowVariable and GetBlockVariable): L lanes of
 * scratch whose values are unspecified until written, valid until the kernel
 * returns. Each call returns a different variable. */
double *tmc_array_group_variable(tmc_array_group *g);

/* MICM's ForEachRow / ForEachBlock: for each lane i of the group, calls f with
 * x[k] = lanes[k] + i for k < n, so f reads and writes lane i of each of its n
 * operands. for_each visits all L lanes, padding included; for_each_strict
 * visits only the real cells. */
typedef void (*tmc_array_lane_fn)(double *const *x, void *ctx);
void tmc_array_group_for_each(tmc_array_group *g, tmc_array_lane_fn f, void *ctx, int n, double *const *lanes);
void tmc_array_group_for_each_strict(tmc_array_group *g, tmc_array_lane_fn f, void *ctx, int n, double *const *lanes);

/* MICM's Reduce and ReduceStrict: as for_each, with the accumulator passed to
 * every call; f combines lane i into it (a sum, a maximum, a logical and/or,
 * of whatever type the caller's accumulator has). reduce includes the padding
 * lanes, so its function must leave the result unchanged on the padding's
 * values; reduce_strict visits only the real cells. */
typedef void (*tmc_array_reduce_fn)(double *const *x, void *accumulator, void *ctx);
void tmc_array_group_reduce(tmc_array_group *g, tmc_array_reduce_fn f, void *ctx, int n, double *const *lanes, void *accumulator);
void tmc_array_group_reduce_strict(tmc_array_group *g, tmc_array_reduce_fn f, void *ctx, int n, double *const *lanes, void *accumulator);

/* MICM's GroupView Fill and Copy: fill sets all L lanes; copy copies all L
 * lanes of src into dst; copy_strict copies only the real cells (MICM's copy
 * of a vector into a dense column, which leaves the column's padding as it
 * was). */
void tmc_array_group_fill(tmc_array_group *g, double *lanes, double value);
void tmc_array_group_copy(tmc_array_group *g, double *dst, const double *src);
void tmc_array_group_copy_strict(tmc_array_group *g, double *dst, const double *src);

#ifdef __cplusplus
}
#endif

#endif
