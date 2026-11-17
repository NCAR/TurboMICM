#ifndef CHEM_CPU_SOLVER_H
#define CHEM_CPU_SOLVER_H
#include "process.h"
#include "state.h"
#include "backward_euler.h"

/* A solver owns its compiled mechanism and state schema. A zero width allows
   any number of cells; a positive width limits each solve to that many cells.
   Option pointers other than absolute_tolerance are borrowed. */
typedef enum { CHEM_SOLVER_ROSENBROCK, CHEM_SOLVER_BACKWARD_EULER } ChemSolverKind;
typedef struct {
    ChemBuiltMechanism built;
    ChemState schema;
    ChemSolverKind kind;
    size_t group_width;
    ChemRosenbrockMethod method;
    ChemRosenbrockOptions rosenbrock;
    ChemBackwardEulerOptions backward_euler;
} ChemCpuSolver;
/* Initialize solver with {0}; compilation is transactional. Parameter labels
   are indexed as in each user/surface rate descriptor. Callback rates are
   rejected by this wrapper; use the composition API for callbacks. */
ChemBuildStatus chem_cpu_solver_init(ChemCpuSolver *solver,const ChemSystem *system,
    size_t process_count,const ChemProcess *processes,int ignore_unused,
    ChemSolverKind kind,size_t group_width,size_t parameter_count,const char *const *parameter_names);
void chem_cpu_solver_destroy(ChemCpuSolver *solver);
size_t chem_cpu_solver_maximum_cells(const ChemCpuSolver *solver);
ChemStateStatus chem_cpu_solver_state(const ChemCpuSolver *solver,size_t cells,ChemState *state);
/* Rates are updated explicitly, allowing caller-set rates to be used in solve. */
ChemStateStatus chem_cpu_solver_update_parameters(const ChemCpuSolver *solver,ChemState *state);
/* NULL options use saved settings; a non-NULL override persists. The state
   supplies tolerances. Multiple cells share a step and convergence decision. */
ChemSolveStatus chem_cpu_solver_solve(ChemCpuSolver *solver,ChemState *state,double duration,
    const ChemRosenbrockOptions *rosenbrock_override,const ChemBackwardEulerOptions *backward_override,
    ChemSolverStats *stats);
#endif
