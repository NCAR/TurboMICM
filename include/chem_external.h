/* SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_EXTERNAL_H
#define CHEM_EXTERNAL_H
#include "chem_solver.h"

typedef struct { const char *row, *column; } ChemNamedJacobianElement;
typedef struct { size_t row, column; } ChemJacobianElement;
/* Resolved indices follow the corresponding descriptor arrays, including after
 * mechanism reordering. Callbacks must not retain these borrowed pointers. */
typedef struct {
    size_t variable_count, parameter_count, algebraic_count;
    const size_t *variables, *parameters, *algebraic;
    size_t process_element_count, constraint_element_count;
    const ChemJacobianElement *process_elements, *constraint_elements;
} ChemModelBinding;
typedef struct {
    size_t order, cells;
    const ChemSparse *sparse; /* NULL denotes dense row/column blocks. */
    double *values;
    double sign;
} ChemModelJacobian;
/* Add a positive derivative; sign supplies the solver's negative-J convention.
 * Every element written by a model must be declared in its descriptor. */
ChemStatus chem_model_jacobian_add(ChemModelJacobian *,size_t row,size_t column,size_t cell,double derivative);
typedef void (*ChemModelUpdate)(const ChemModelBinding *,const void *,ChemState *);
typedef void (*ChemModelForcing)(const ChemModelBinding *,const void *,const ChemState *,const double *,double *);
typedef void (*ChemModelDerivative)(const ChemModelBinding *,const void *,const ChemState *,const double *,ChemModelJacobian *);
/* Solver storage descriptors replace C++ policy type information. Dense
 * values are variable-major contiguous cell blocks; sparse values use the
 * supplied pattern with the same contiguous cell blocks. No padding is added. */
typedef struct { size_t maximum_cells;int variable_major; } ChemDenseStorage;
typedef struct { size_t maximum_cells;const ChemSparse *pattern; } ChemSparseStorage;
typedef struct {
    size_t additional_species_count;
    const ChemSpecies *additional_species;
    size_t variable_count, parameter_count, algebraic_count;
    const char *const *variables, *const *parameters, *const *algebraic;
    size_t process_element_count, constraint_element_count;
    const ChemNamedJacobianElement *process_elements, *constraint_elements;
    /* A context is owned by each constructed mechanism. Plain data may use
     * context_size; pointer-containing data must supply clone and destroy. */
    const void *context;
    size_t context_size;
    void *(*clone)(const void *);
    void (*destroy)(void *);
    ChemModelUpdate update_process_parameters, update_constraint_parameters;
    ChemModelForcing add_forcing, constraint_residual;
    ChemModelDerivative process_jacobian, constraint_jacobian;
    /* Called once by chem_solver_build after index resolution and reordering.
     * Store indices, not borrowed descriptor pointers. Full signature takes
     * precedence when both are supplied. The context is the owned model copy. */
    ChemStatus (*finalize)(const ChemModelBinding *,void *,const ChemDenseStorage *,const ChemSparseStorage *);
    ChemStatus (*finalize_sparse)(const ChemModelBinding *,void *,const ChemSparseStorage *);
} ChemExternalModel;
/* Inputs may be reused for multiple builds and released after construction.
 * Reaction indices address the initial species followed by each model's added
 * species, in descriptor order. Additional species must be solution variables. */
ChemStatus chem_mechanism_create_with_models(size_t,const ChemSpecies *,size_t,const ChemReaction *,
                                             size_t,const ChemExternalModel *,ChemMechanism **);
size_t chem_external_model_count(const ChemMechanism *);
/* Separate process accumulation is useful for model Jacobian verification.
 * The composed forcing/Jacobian APIs replace all algebraic rows afterward. */
void chem_external_add_forcing(const ChemState *,const double *,double *);
void chem_external_subtract_jacobian(const ChemState *,const double *,double *);
int chem_external_has_element(const ChemMechanism *,size_t row,size_t column,int constraints);
#endif
