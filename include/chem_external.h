#ifndef CHEM_EXTERNAL_H
#define CHEM_EXTERNAL_H
#include "chem.h"

typedef struct { const char *row,*column; } chem_named_element;
/* External models may supply any combination of state, kinetic processes and
 * algebraic constraints. Names and declared patterns are consumed at model
 * creation. Finalize runs once, after reordering and symbolic factorization.
 * Kernels use the same SoA and negative-Jacobian conventions as core chemistry.
 * A constraint residual overwrites its algebraic rows; its Jacobian adds -dG/dy.
 * Callback contexts are immutable during updates/solves except for per-state
 * data accessed through the state argument. For independent model bindings,
 * either context_size copies a plain value or clone/destroy manage ownership.
 * With neither, context is borrowed and must outlive the compiled model. */
typedef struct chem_external_model {
    size_t n_species; const chem_species *species;
    size_t n_parameters; const chem_named_value *parameters;
    size_t n_process_elements; const chem_named_element *process_elements;
    size_t n_algebraic; const char *const *algebraic;
    size_t n_constraint_elements; const chem_named_element *constraint_elements;
    chem_external_kernel forcing,jacobian,residual,constraint_jacobian;
    chem_status (*finalize)(const chem_model *,void *);
    void (*update)(chem_state *,void *);
    const void *context;
    size_t context_size;
    void *(*clone)(const void *);
    void (*destroy)(void *);
} chem_external_model;
/* Borrowed binding context, useful for implementing model-specific queries. */
const void *chem_model_external_context(const chem_model *,size_t);
#endif
