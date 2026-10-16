#ifndef CHEM_SOLVER_H
#define CHEM_SOLVER_H
#include "chem_process.h"
#include "chem_sparse.h"

typedef enum { CHEM_SOLVER_ROSENBROCK, CHEM_SOLVER_BACKWARD_EULER } chem_solver_method;
/* Preferred cells per vector group. States may contain multiple groups. */
typedef struct { size_t vector_length; chem_sparse_order sparse_order; } chem_solver_layout;
typedef struct {
    chem_status (*finalize)(const chem_model *,const chem_solver_layout *,void *);
    chem_status (*finalize_sparse)(const chem_model *,chem_sparse_order,size_t,void *);
    void *context;
} chem_setup_callback;
typedef struct {
    chem_solver_method method;
    chem_ros_parameters rosenbrock;
    chem_be_parameters backward_euler;
    chem_lu_algorithm lu;
    chem_solver_layout layout;
    int preserve_species_order, reject_unused;
    size_t n_setups;
    const chem_setup_callback *setups;
} chem_solver_options;
chem_solver_options chem_solver_options_default(chem_solver_method);
typedef struct { chem_process_model chemistry; chem_solver_options options; } chem_solver;
chem_status chem_solver_create(chem_solver *,const chem_system *,size_t,const chem_process *,const chem_solver_options *,chem_error *);
void chem_solver_destroy(chem_solver *);
size_t chem_solver_max_cells(const chem_solver *);
chem_status chem_solver_state(const chem_solver *,chem_state *,size_t);
/* A solver can be shared by threads that each own their state. Explicit
 * parameter overrides below mutate the solver and require exclusive access. */
chem_result chem_solver_solve(const chem_solver *,chem_state *,double);
chem_result chem_solver_solve_ros(chem_solver *,chem_state *,double,const chem_ros_parameters *,chem_step_controller,void *);
chem_result chem_solver_solve_be(chem_solver *,chem_state *,double,const chem_be_parameters *);
#endif
