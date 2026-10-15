/* SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_BUILDER_H
#define CHEM_BUILDER_H
#include "chem_system.h"
#include "chem_external.h"
#include "chem_constraints.h"

typedef enum { CHEM_SOLVER_ROSENBROCK, CHEM_SOLVER_BACKWARD_EULER } ChemSolverMethod;
typedef struct {
    ChemSolverMethod method;
    union { ChemRosenbrockParameters rosenbrock;ChemBackwardEulerParameters backward_euler; } parameters;
} ChemSolverOptions;
/* Builder pointers are borrowed until build returns. A plan owns an independent
 * mechanism, system and model snapshot; the builder may be reused afterward. */
typedef struct {
    const ChemSystem *system;
    size_t reaction_count;
    const ChemReactionInfo *const *reactions;
    size_t model_count;
    const ChemExternalModel *models;
    size_t constraint_count;
    const ChemConstraint *constraints;
    size_t maximum_cells; /* zero: no explicit cell limit */
    int ignore_unused_species, reorder;
    ChemSolverOptions options;
} ChemSolverBuilder;
typedef struct ChemSolverPlan ChemSolverPlan;
ChemSolverBuilder chem_solver_builder(ChemSolverMethod);
ChemStatus chem_solver_build(const ChemSolverBuilder *,ChemSolverPlan **);
void chem_solver_plan_free(ChemSolverPlan *);
const ChemMechanism *chem_solver_plan_mechanism(const ChemSolverPlan *);
const ChemSystem *chem_solver_plan_system(const ChemSolverPlan *);
size_t chem_solver_plan_maximum_cells(const ChemSolverPlan *);
const ChemSolverOptions *chem_solver_plan_options(const ChemSolverPlan *);
/* Free the returned state and workspace before freeing their plan. Independent
 * workspaces can share a plan; changing retained options requires serialization. */
ChemStatus chem_solver_plan_state(const ChemSolverPlan *,size_t cells,ChemState **);
ChemStatus chem_solver_plan_workspace(const ChemSolverPlan *,size_t cells,ChemSolver **);
/* A non-NULL override replaces retained options. NULL reuses the last options.
 * Call chem_update_rates explicitly when conditions or custom parameters change. */
ChemStatus chem_solver_plan_solve(ChemSolverPlan *,ChemSolver *,ChemState *,double dt,
                                  const ChemSolverOptions *override,ChemSolverStats *);
#endif
