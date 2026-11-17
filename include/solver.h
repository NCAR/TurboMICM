#ifndef CHEM_SOLVER_H
#define CHEM_SOLVER_H
#include "state.h"
#include <stdint.h>

typedef enum { CHEM_SOLVER_ROSENBROCK, CHEM_SOLVER_BACKWARD_EULER } ChemSolverMethod;
typedef enum { CHEM_BUILD_OK, CHEM_BUILD_MISSING_SYSTEM, CHEM_BUILD_UNUSED_SPECIES, CHEM_BUILD_INVALID, CHEM_BUILD_NO_MEMORY } ChemBuildStatus;
typedef struct {
    const ChemSystem *system;
    const ChemMechanism *mechanism;
    const ChemConstraint *constraints;
    size_t constraint_count;
    const ChemExternalModel *models;
    size_t model_count;
    size_t group_width; /* zero selects unbounded standard storage */
    int ignore_unused_species;
    ChemSolverMethod method;
    ChemRosenbrockParameters rosenbrock;
    ChemBackwardEulerParameters backward_euler;
} ChemSolverBuilder;
typedef struct {
    ChemSolverBuilder options;
    ChemWorkspace *workspace;
} ChemSolver;
ChemSolverBuilder chem_solver_builder(ChemSolverMethod);
ChemBuildStatus chem_solver_build(const ChemSolverBuilder *, ChemSolver *);
void chem_solver_free(ChemSolver *);
size_t chem_solver_maximum_cells(const ChemSolver *);
/* Selects the default parameters when override is null. A non-null override
 * is copied into the solver and persists for subsequent calls. */
int chem_solver_solve(ChemSolver *, double dt, ChemState *, const void *override, ChemStatistics *);
#endif
