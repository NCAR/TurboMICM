#include "internal.h"
#include <omp.h>

ChemStatus chem_step_rosenbrock(void *solver, double dt, ChemState *state, ChemStats *stats) {
  return chem_solve(solver,dt,state,stats);
}
ChemStatus chem_step_euler(void *parameters, double dt, ChemState *state, ChemStats *stats) {
  return chem_backward_euler(parameters,dt,state,stats);
}
ChemStatus chem_advance(ChemStepFunction step, void *parameters, double dt, size_t substeps, ChemState *state, ChemStats *stats) {
  if (!step || !parameters || !state || !stats || !substeps || !(dt >= 0.0) || !isfinite(dt)) return CHEM_INVALID;
  memset(stats,0,sizeof *stats);
  double present = 0.0;
  for (size_t i = 0; i < substeps; ++i) {
    double target = dt*((double)(i+1)/(double)substeps);
    while (present < target) {
      ChemStats one = {0};
      ChemStatus status = step(parameters,target-present,state,&one);
      stats->internal_steps += one.internal_steps; stats->accepted += one.accepted; stats->rejected += one.rejected;
      stats->function_calls += one.function_calls; stats->jacobian_updates += one.jacobian_updates;
      stats->decompositions += one.decompositions; stats->solves += one.solves;
      stats->constraint_init_iterations += one.constraint_init_iterations;
      if (!(one.final_time > 0.0) || present+one.final_time == present) return status == CHEM_OK ? CHEM_SMALL_STEP : status;
      present += one.final_time; stats->final_time = present;
      if (status != CHEM_OK && status != CHEM_MAX_STEPS && status != CHEM_UNCONVERGED) return status;
    }
  }
  return CHEM_OK;
}

ChemStatus chem_advance_many(ChemStepFunction step, void *parameters, double dt, size_t substeps,
                            ChemState **states, size_t count, int threads, ChemStats *stats) {
  if (!states || !stats || !count || threads < 0 || count > SIZE_MAX/sizeof(ChemStatus)) return CHEM_INVALID;
  ChemStatus *status = malloc(count*sizeof *status);
  if (!status) return CHEM_MEMORY;
  int workers = threads ? threads : omp_get_max_threads();
  #pragma omp parallel for schedule(static) num_threads(workers)
  for (size_t i = 0; i < count; ++i) status[i] = chem_advance(step,parameters,dt,substeps,states[i],stats+i);
  ChemStatus result = CHEM_OK;
  for (size_t i = 0; i < count; ++i) if (status[i]) { result = status[i]; break; }
  free(status);
  return result;
}
