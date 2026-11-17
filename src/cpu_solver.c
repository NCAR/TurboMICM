#include "cpu_solver.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void chem_cpu_solver_destroy(ChemCpuSolver *s)
{
    if (!s) return;
    chem_built_mechanism_destroy(&s->built); chem_state_destroy(&s->schema); *s=(ChemCpuSolver){0};
}
ChemBuildStatus chem_cpu_solver_init(ChemCpuSolver *s,const ChemSystem *system,size_t count,
    const ChemProcess *processes,int ignore_unused,ChemSolverKind kind,size_t width,
    size_t parameter_count,const char *const *parameter_names)
{
    if (!s || (kind!=CHEM_SOLVER_ROSENBROCK && kind!=CHEM_SOLVER_BACKWARD_EULER) ||
        (parameter_count && !parameter_names) || (count && !processes)) return CHEM_BUILD_INVALID;
    for (size_t i=0;i<count;++i) {
        const ChemRate *r=&processes[i].rate;
        size_t need=r->type==CHEM_RATE_SURFACE ? 2 : r->type==CHEM_RATE_USER ? 1 : 0;
        if (r->type==CHEM_RATE_CALLBACK || (need && (r->custom_index>parameter_count || need>parameter_count-r->custom_index)))
            return CHEM_BUILD_INVALID;
    }
    ChemCpuSolver copy={0};
    ChemBuildStatus result=chem_build_mechanism(&copy.built,system,count,processes,ignore_unused);
    if (result!=CHEM_BUILD_OK) return result;
    ChemStateConfig config={0,copy.built.mechanism.species_count,parameter_count,count,
        copy.built.state_names,parameter_names,1e-6,copy.built.absolute_tolerances};
    ChemStateStatus state_result=chem_state_init(&copy.schema,&config);
    if (state_result!=CHEM_STATE_OK) {
        chem_cpu_solver_destroy(&copy);
        return state_result==CHEM_STATE_NO_MEMORY || state_result==CHEM_STATE_SIZE ? CHEM_BUILD_NO_MEMORY : CHEM_BUILD_INVALID;
    }
    copy.kind=kind; copy.group_width=width; copy.method=chem_rosenbrock_three_stage();
    copy.rosenbrock=chem_rosenbrock_defaults(NULL); copy.backward_euler=chem_backward_euler_defaults(NULL);
    chem_cpu_solver_destroy(s); *s=copy; return CHEM_BUILD_OK;
}
size_t chem_cpu_solver_maximum_cells(const ChemCpuSolver *s)
{ return !s ? 0 : s->group_width ? s->group_width : SIZE_MAX; }
ChemStateStatus chem_cpu_solver_state(const ChemCpuSolver *s,size_t cells,ChemState *state)
{
    if (!s || !state || !s->built.system.gas.name || !cells || cells>chem_cpu_solver_maximum_cells(s)) return CHEM_STATE_INVALID;
    ChemStateConfig config={cells,s->schema.species_count,s->schema.parameter_count,s->schema.rate_count,
        (const char *const *)s->schema.species_names,(const char *const *)s->schema.parameter_names,
        s->schema.relative_tolerance,s->schema.absolute_tolerances};
    return chem_state_init(state,&config);
}
static int compatible(const ChemCpuSolver *s,const ChemState *state)
{
    if (!s || !state || !s->built.system.gas.name || !state->cells || state->cells>chem_cpu_solver_maximum_cells(s) ||
        state->species_count!=s->schema.species_count || state->parameter_count!=s->schema.parameter_count ||
        state->rate_count!=s->schema.rate_count || !state->conditions || !state->concentrations ||
        (state->rate_count && !state->rates) || (state->parameter_count && !state->parameters) || !state->absolute_tolerances)
        return 0;
    for (size_t i=0;i<state->species_count;++i)
        if (!state->species_names || !state->species_names[i] || strcmp(state->species_names[i],s->schema.species_names[i])) return 0;
    for (size_t i=0;i<state->parameter_count;++i)
        if (!state->parameter_names || !state->parameter_names[i] || strcmp(state->parameter_names[i],s->schema.parameter_names[i])) return 0;
    return 1;
}
ChemStateStatus chem_cpu_solver_update_parameters(const ChemCpuSolver *s,ChemState *state)
{
    if (!compatible(s,state)) return CHEM_STATE_INVALID;
    for (size_t cell=0;cell<state->cells;++cell)
        chem_rates_parameterized(state->rate_count,s->built.rates,s->built.multipliers,s->built.functions,
            state->conditions[cell],state->parameter_count ? state->parameters+cell*state->parameter_count : NULL,
            state->rate_count ? state->rates+cell*state->rate_count : NULL);
    return CHEM_STATE_OK;
}
typedef struct { const ChemMechanism *mechanism; const double *rates; } SolverCell;
static void rhs(void *v,const double *y,double *f)
{ SolverCell *c=v; chem_rhs(c->mechanism,y,c->rates,f); }
static void jacobian(void *v,const double *y,double *j)
{ SolverCell *c=v; chem_jacobian(c->mechanism,y,c->rates,j); }
ChemSolveStatus chem_cpu_solver_solve(ChemCpuSolver *s,ChemState *state,double duration,
    const ChemRosenbrockOptions *ro,const ChemBackwardEulerOptions *be,ChemSolverStats *stats)
{
    if (stats) memset(stats,0,sizeof(*stats));
    if (!compatible(s,state) || (ro && be) || (ro && s->kind!=CHEM_SOLVER_ROSENBROCK) ||
        (be && s->kind!=CHEM_SOLVER_BACKWARD_EULER)) return CHEM_SOLVE_INVALID;
    if (ro) { s->rosenbrock=*ro; s->rosenbrock.absolute_tolerance=NULL; }
    if (be) { s->backward_euler=*be; s->backward_euler.absolute_tolerance=NULL; }
    size_t cells=state->cells;
    if (cells>SIZE_MAX/sizeof(SolverCell) || cells>SIZE_MAX/sizeof(void *)) return CHEM_SOLVE_ALLOCATION;
    SolverCell *data=malloc(cells*sizeof(*data)); void **contexts=malloc(cells*sizeof(*contexts));
    if (!data || !contexts) { free(data); free(contexts); return CHEM_SOLVE_ALLOCATION; }
    for (size_t i=0;i<cells;++i) {
        data[i]=(SolverCell){&s->built.mechanism,state->rate_count ? state->rates+i*state->rate_count : NULL};
        contexts[i]=data+i;
    }
    ChemSolveStatus result;
    if (s->kind==CHEM_SOLVER_ROSENBROCK) {
        ChemRosenbrockOptions options=s->rosenbrock;
        options.absolute_tolerance=state->absolute_tolerances; options.relative_tolerance=state->relative_tolerance;
        result=cells==1 ? chem_rosenbrock_solve(state->species_count,state->concentrations,duration,rhs,jacobian,data,&s->method,&options,stats) :
            chem_rosenbrock_solve_grouped(cells,state->species_count,state->concentrations,duration,rhs,jacobian,contexts,NULL,&s->method,&options,stats);
    } else {
        ChemBackwardEulerOptions options=s->backward_euler;
        options.absolute_tolerance=state->absolute_tolerances; options.relative_tolerance=state->relative_tolerance;
        result=cells==1 ? chem_backward_euler_solve(state->species_count,state->concentrations,duration,rhs,jacobian,data,&options,stats) :
            chem_backward_euler_solve_grouped(cells,state->species_count,state->concentrations,duration,rhs,jacobian,contexts,NULL,&options,stats);
    }
    free(data); free(contexts); return result;
}
