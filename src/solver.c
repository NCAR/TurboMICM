#include "chem_solver.h"
#include <stdio.h>

chem_solver_options chem_solver_options_default(chem_solver_method method)
{
    return (chem_solver_options){.method=method,.rosenbrock=chem_ros_parameters_default(CHEM_ROS3),
        .backward_euler=chem_be_parameters_default(),.lu=CHEM_LU_DOOLITTLE,.layout={0,CHEM_CSR}};
}
void chem_solver_destroy(chem_solver *s)
{ if (s) { chem_process_model_destroy(&s->chemistry); *s=(chem_solver){0}; } }
chem_status chem_solver_create(chem_solver *s,const chem_system *system,size_t count,const chem_process *processes,
    const chem_solver_options *options,chem_error *error)
{
    if (!s || !options || (options->n_setups && !options->setups) ||
        options->method>CHEM_SOLVER_BACKWARD_EULER || options->layout.sparse_order>CHEM_CSC) return CHEM_INVALID;
    *s=(chem_solver){0}; s->options=*options;
    chem_status status=chem_process_model_create(&s->chemistry,system,count,processes,options->preserve_species_order,
        options->n_external ? 0 : options->reject_unused,options->lu,error);
    if (status!=CHEM_OK) return status;
    if (options->n_external) {
        chem_model_destroy(&s->chemistry.model);
        s->chemistry.mechanism.n_external=options->n_external; s->chemistry.mechanism.external=options->external;
        status=chem_model_create_checked(&s->chemistry.model,&s->chemistry.mechanism,options->lu,error);
        if (status==CHEM_OK && options->reject_unused) {
            const chem_model *model=&s->chemistry.model;
            for (size_t i=0;i<model->mechanism->n_species;++i) if (!chem_model_species_used(model,i)) {
                status=CHEM_INVALID;
                if (error) { error->code=status; snprintf(error->message,sizeof(error->message),"Unused species: %s",model->mechanism->species[i].name); }
                break;
            }
        }
        if (status!=CHEM_OK) { chem_solver_destroy(s); return status; }
        /* Extension descriptors have been consumed; model bindings own copies. */
        s->chemistry.mechanism.n_external=0; s->chemistry.mechanism.external=NULL;
        s->options.n_external=0; s->options.external=NULL;
    }
    for (size_t i=0;i<options->n_setups;++i) {
        const chem_setup_callback *setup=options->setups+i;
        if (setup->finalize) status=setup->finalize(&s->chemistry.model,&options->layout,setup->context);
        else if (setup->finalize_sparse) status=setup->finalize_sparse(&s->chemistry.model,options->layout.sparse_order,options->layout.vector_length,setup->context);
        else status=CHEM_INVALID;
        if (status!=CHEM_OK) {
            if (error) { error->code=status; snprintf(error->message,sizeof(error->message),"Model setup callback failed: %s",chem_status_string(status)); }
            chem_solver_destroy(s); return status;
        }
    }
    /* Setup callbacks run once, so their descriptor array need not outlive creation. */
    s->options.n_setups=0; s->options.setups=NULL;
    return CHEM_OK;
}
size_t chem_solver_max_cells(const chem_solver *s)
{ return s->options.layout.vector_length ? s->options.layout.vector_length : SIZE_MAX; }
chem_status chem_solver_state(const chem_solver *s,chem_state *state,size_t cells)
{
    if (!s || !state || !s->chemistry.model.mechanism) return CHEM_INVALID;
    return chem_state_create(state,&s->chemistry.model,cells);
}
chem_result chem_solver_solve_ros(chem_solver *solver,chem_state *state,double dt,const chem_ros_parameters *parameters,
    chem_step_controller controller,void *context)
{
    if (!solver || !state || !parameters || state->model!=&solver->chemistry.model)
        return (chem_result){.status=CHEM_INVALID};
    solver->options.rosenbrock=*parameters;
    chem_result result=chem_ros_solve(state,dt,&solver->options.rosenbrock,controller,context);
    chem_state_clamp(state); return result;
}
chem_result chem_solver_solve_be(chem_solver *solver,chem_state *state,double dt,const chem_be_parameters *parameters)
{
    if (!solver || !state || !parameters || state->model!=&solver->chemistry.model)
        return (chem_result){.status=CHEM_INVALID};
    solver->options.backward_euler=*parameters;
    chem_result result=chem_be_solve(state,dt,&solver->options.backward_euler);
    chem_state_clamp(state); return result;
}
chem_result chem_solver_solve(const chem_solver *solver,chem_state *state,double dt)
{
    if (!solver || !state || state->model!=&solver->chemistry.model) return (chem_result){.status=CHEM_INVALID};
    chem_result result=solver->options.method==CHEM_SOLVER_BACKWARD_EULER ?
        chem_be_solve(state,dt,&solver->options.backward_euler) :
        chem_ros_solve(state,dt,&solver->options.rosenbrock,NULL,NULL);
    chem_state_clamp(state); return result;
}
