#include "chem.h"
#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef CHEM_SPECIALIZED
#include "chem_specialized.h"
extern const chem_specialized ts1_specialized;
#define TS1_MECHANISM (*ts1_specialized.mechanism)
#else
extern const chem_mechanism ts1_mechanism;
#define TS1_MECHANISM ts1_mechanism
#endif
#ifndef CHEM_GROUP_SIZE
#define CHEM_GROUP_SIZE 32
#endif

#include "bench_io.h"

int main(int argc, char **argv)
{
    size_t threads=(size_t)omp_get_num_procs(), group_size=CHEM_GROUP_SIZE;
    int first=1, generic=0;
    omp_sched_t schedule=omp_sched_static; size_t chunk=0;
    while (first<argc && argv[first][0]=='-') {
        if (!strcmp(argv[first],"--generic")) { generic=1; ++first; continue; }
        if (!strcmp(argv[first],"--schedule")) {
            if (first+1>=argc) goto usage;
            if (!strcmp(argv[first+1],"static")) schedule=omp_sched_static;
            else if (!strcmp(argv[first+1],"dynamic")) schedule=omp_sched_dynamic;
            else if (!strcmp(argv[first+1],"guided")) schedule=omp_sched_guided;
            else goto usage;
            first+=2; continue;
        }
        size_t value;
        if (first+1>=argc || !integer(argv[first+1],&value,0)) goto usage;
        if (!strcmp(argv[first],"--threads")) threads=value;
        else if (!strcmp(argv[first],"--group-size")) group_size=value;
        else if (!strcmp(argv[first],"--chunk")) chunk=value;
        else goto usage;
        first+=2;
    }
    size_t cells,steps;
    if (argc-first!=4 || !integer(argv[first],&cells,0) || !integer(argv[first+1],&steps,1) || threads>INT_MAX || chunk>INT_MAX || !group_size) goto usage;
    char *end; errno=0; double dt=strtod(argv[first+2],&end);
    if (*end || errno || !isfinite(dt) || dt<=0) goto usage;
    omp_set_num_threads((int)threads);
    omp_set_schedule(schedule,(int)chunk);
    chem_model model;
    chem_status status=chem_model_create(&model,&TS1_MECHANISM,CHEM_LU_MOZART);
    if (status!=CHEM_OK) { fprintf(stderr,"Model: %s\n",chem_status_string(status)); return 1; }
    if (generic) { model.forcing=chem_forcing; model.jacobian=chem_jacobian; }
    size_t groups=cells/group_size+(cells%group_size!=0);
    chem_state *states=calloc(groups,sizeof(chem_state));
    chem_result *results=calloc(groups,sizeof(chem_result));
    int success=0;
    if (!states || !results) goto done;
    for (size_t g=0; g<groups; ++g) {
        size_t count=cells-g*group_size;
        if (count>group_size) count=group_size;
        status=chem_state_create(&states[g],&model,count);
        if (status!=CHEM_OK) { fprintf(stderr,"State: %s\n",chem_status_string(status)); goto done; }
    }
    if (!read_cells(argv[first+3],states,group_size,cells)) goto done;
    chem_ros_parameters parameters=chem_ros_parameters_default(CHEM_ROS3);
    #pragma omp parallel for schedule(runtime)
    for (size_t g=0; g<groups; ++g) {
#ifdef CHEM_SPECIALIZED
        if (!generic) ts1_specialized.rates(&states[g]);
        else
#endif
        chem_calculate_rates(&states[g]);
        results[g]=
#ifdef CHEM_SPECIALIZED
            !generic ? ts1_specialized.solve(&states[g],dt) :
#endif
            chem_ros_solve(&states[g],dt,&parameters,NULL,NULL);
        chem_state_clamp(&states[g]);
    }
    for (size_t g=0; g<groups; ++g) if (results[g].status!=CHEM_OK) {
        fprintf(stderr,"Warm-up group %zu: %s\n",g,chem_status_string(results[g].status)); goto done;
    }
    memset(results,0,groups*sizeof(chem_result));
    double start=omp_get_wtime();
    #pragma omp parallel for schedule(runtime)
    for (size_t g=0; g<groups; ++g) {
        for (size_t step=0; step<steps; ++step) {
            chem_result r=
#ifdef CHEM_SPECIALIZED
                !generic ? ts1_specialized.solve(&states[g],dt) :
#endif
                chem_ros_solve(&states[g],dt,&parameters,NULL,NULL);
            chem_state_clamp(&states[g]);
            chem_stats_add(&results[g].stats,&r.stats);
            if (r.status!=CHEM_OK) { results[g].status=r.status; break; }
        }
    }
    double elapsed=(omp_get_wtime()-start)*1000.0;
    chem_stats stats={0};
    for (size_t g=0; g<groups; ++g) {
        if (results[g].status!=CHEM_OK) { fprintf(stderr,"Group %zu: %s\n",g,chem_status_string(results[g].status)); goto done; }
        chem_stats_add(&stats,&results[g].stats);
    }
    printf("elapsed_ms=%.6f\n",elapsed);
    printf("internal_steps=%" PRIu64 " accepted=%" PRIu64 " rejected=%" PRIu64 " function_calls=%" PRIu64
           " jacobian_updates=%" PRIu64 " decompositions=%" PRIu64 " solves=%" PRIu64 "\n",stats.internal_steps,
           stats.accepted,stats.rejected,stats.function_calls,stats.jacobian_updates,stats.decompositions,stats.solves);
    success=write_final(states,groups);
done:
    if (states) for (size_t g=0; g<groups; ++g) chem_state_destroy(&states[g]);
    free(states); free(results); chem_model_destroy(&model); return success ? 0 : 1;
usage:
    fprintf(stderr,"Usage: %s [--threads N] [--group-size G] [--schedule static|dynamic|guided] [--chunk K] [--generic] CELLS STEPS DT CELLS_CSV\n",argv[0]);
    return 2;
}
