#include "ts1_driver.h"
int main(int argc,char **argv)
{
    int threads=omp_get_max_threads(),arg=1;
    if (argc>2 && !strcmp(argv[1],"--threads")) {
        size_t n=positive_integer(argv[2]);
        if (!n || n>INT_MAX) { fprintf(stderr,"Invalid thread count\n");return 1; }
        threads=(int)n;arg=3;
    }
    if (argc-arg!=4) {
        fprintf(stderr,"Usage: %s [--threads N] CELLS STEPS DT CELLS_CSV\n",argv[0]);return 1;
    }
    size_t cells=positive_integer(argv[arg]),steps=positive_integer(argv[arg+1]);
    char *end; double dt=strtod(argv[arg+2],&end);
    if (!cells || !steps || *end || !(dt>0) || !isfinite(dt)) { fprintf(stderr,"Invalid run dimensions\n");return 1; }
    size_t group_size=CHEM_GROUP_SIZE;
    const char *setting=getenv("CHEM_GROUP_SIZE");
    if (setting) group_size=positive_integer(setting);
    if (!group_size) { fprintf(stderr,"Invalid group size\n");return 1; }
    if (group_size>cells) group_size=cells;
    size_t ng=1+(cells-1)/group_size;
    ChemMechanism *m=NULL;
    Group *groups=calloc(ng,sizeof(Group));
    int success=0;
    if (!groups || create_mechanism(&m)) goto finish;
    for (size_t g=0;g<ng;++g) {
        size_t count=cells-g*group_size;
        if (count>group_size) count=group_size;
        if (chem_state_create(m,count,&groups[g].state) || ts1_initialize(groups[g].state) ||
            create_solver(m,count,&groups[g].solver)) goto finish;
    }
    if (!load_cells(argv[arg+3],groups,cells,group_size,m)) goto finish;
    omp_set_num_threads(threads);
    if (!getenv("OMP_SCHEDULE")) omp_set_schedule(omp_sched_static,0);
    ChemRosenbrockParameters parameters=chem_rosenbrock_parameters(CHEM_ROS3);
    #pragma omp parallel for schedule(runtime)
    for (size_t g=0;g<ng;++g) {
        chem_update_rates(groups[g].state);
        groups[g].status=chem_solve(groups[g].solver,groups[g].state,dt,&parameters,&groups[g].stats);
    }
    for (size_t g=0;g<ng;++g) if (groups[g].status) {
        fprintf(stderr,"Warm-up: %s\n",chem_status_string(groups[g].status));goto finish;
    }
    double start=omp_get_wtime();
    #pragma omp parallel for schedule(runtime)
    for (size_t g=0;g<ng;++g) {
        memset(&groups[g].stats,0,sizeof(ChemSolverStats));
        for (size_t step=0;step<steps;++step) {
            ChemSolverStats stats;
            chem_update_rates(groups[g].state);
            groups[g].status=chem_solve(groups[g].solver,groups[g].state,dt,&parameters,&stats);
            add_stats(&groups[g].stats,&stats);
            if (groups[g].status) break;
        }
    }
    double elapsed=1000*(omp_get_wtime()-start);
    ChemSolverStats stats={0};
    for (size_t g=0;g<ng;++g) {
        if (groups[g].status) { fprintf(stderr,"Solve: %s\n",chem_status_string(groups[g].status));goto finish; }
        add_stats(&stats,&groups[g].stats);
    }
    printf("elapsed_ms=%.6f\n",elapsed);
    printf("internal_steps=%zu accepted=%zu rejected=%zu function_calls=%zu jacobian_updates=%zu decompositions=%zu solves=%zu\n",
        stats.internal_steps,stats.accepted,stats.rejected,stats.function_calls,stats.jacobian_updates,stats.decompositions,stats.solves);
    FILE *output=fopen("ts1_final.csv","w");
    if (!output) { perror("ts1_final.csv");goto finish; }
    fprintf(output,"temperature,pressure");
    for (size_t v=0;v<chem_species_count(m);++v) fprintf(output,",%s",chem_species_name(m,v));
    fputc('\n',output);
    for (size_t cell=0;cell<cells;++cell) {
        ChemState *s=groups[cell/group_size].state; size_t c=cell%group_size;
        fprintf(output,"%.17g,%.17g",s->conditions[c].temperature,s->conditions[c].pressure);
        for (size_t v=0;v<chem_species_count(m);++v) fprintf(output,",%.17g",s->values[v*s->cells+c]);
        fputc('\n',output);
    }
    int io_error=ferror(output); if (fclose(output)) io_error=1;
    success=!io_error;
finish:
    if (groups) for (size_t g=0;g<ng;++g) { chem_solver_free(groups[g].solver);chem_state_free(groups[g].state); }
    free(groups);chem_mechanism_free(m);
    return success?0:1;
}
