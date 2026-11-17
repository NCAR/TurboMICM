#ifndef CHEM_MODEL_H
#define CHEM_MODEL_H
#include "composite.h"
#include "sparse_matrix.h"
#include "sparse.h"
#include "state.h"

typedef struct {
    ChemStateFragment state;
    ChemOdeComponent process;
    size_t jacobian_count;
    const ChemSparseEntry *jacobian_entries; /* indices in the combined state */
    size_t constraint_count;
    const ChemConstraint *constraints; /* indices in the combined state */
    /* Called once per build after state and sparse structure are constructed. */
    void (*finalize)(void *context,const ChemState *state,const ChemSparseMatrix *structure);
} ChemExternalModel;
typedef enum { CHEM_MODEL_OK, CHEM_MODEL_INVALID, CHEM_MODEL_NO_MEMORY, CHEM_MODEL_STATE_ERROR } ChemModelStatus;
typedef struct {
    ChemState state;
    ChemCompositeOde ode;
    ChemOdeComponent *components;
    ChemConstraintSet constraints;
    ChemConstraint *constraint_rows;
    ChemSparseMatrix structure;
    ChemSparsePlan *plan;
} ChemModelSystem;
/* Initialize destination to zero. Inputs and component contexts are borrowed
   for its lifetime, including constraint terms and callback contexts. State
   names, lists, and sparse structures are independently owned. base_pattern,
   if present, is a combined-state row-major Jacobian mask. External process
   entries are kept even on algebraic rows. Missing absolute tolerances default
   to 1e-3. Finalizers receive the destination's owned state and structure. */
ChemModelStatus chem_model_system_init(ChemModelSystem *destination,const ChemStateConfig *base,
    const ChemOdeComponent *base_process,const unsigned char *base_pattern,
    const ChemConstraintSet *built_in,size_t model_count,const ChemExternalModel *models,
    size_t group_width,ChemSparseOrdering ordering);
void chem_model_system_destroy(ChemModelSystem *system);
#endif
