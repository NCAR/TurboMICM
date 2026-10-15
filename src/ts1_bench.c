/* SPDX-License-Identifier: Apache-2.0 */
#include "chem_solver.h"
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

ChemStatus ts1_create(ChemMechanism **);
ChemStatus ts1_initialize(ChemState *);
#ifndef CHEM_GROUP_SIZE
#define CHEM_GROUP_SIZE 32
#endif

typedef struct { ChemState *state; ChemSolver *solver; ChemSolverStats stats; ChemStatus status; } Group;
static size_t positive_integer(const char *text)
{
    char *end;
    errno=0;
    unsigned long long v=strtoull(text,&end,10);
    if (errno || *end || end==text || text[0]=='-' || !v || v>SIZE_MAX) return 0;
    return (size_t)v;
}
static char *read_line(FILE *f)
{
    size_t n=0,capacity=256;
    char *line=malloc(capacity);
    if (!line) return NULL;
    int c;
    while ((c=fgetc(f))!=EOF && c!='\n') {
        if (n+1>=capacity) {
            if (capacity>SIZE_MAX/2) { free(line); return NULL; }
            capacity*=2;
            char *next=realloc(line,capacity);
            if (!next) { free(line); return NULL; }
            line=next;
        }
        if (c!='\r') line[n++]=(char)c;
    }
    if (c==EOF && !n) { free(line); return NULL; }
    line[n]='\0'; return line;
}
/* The supplied input format has unquoted, comma-separated names and numbers. */
static char *field(char **cursor)
{
    if (!*cursor) return NULL;
    char *value=*cursor,*comma=strchr(value,',');
    if (comma) { *comma='\0'; *cursor=comma+1; } else *cursor=NULL;
    return value;
}
static int load_cells(const char *filename,Group *groups,size_t cells,size_t group_size,const ChemMechanism *m)
{
    FILE *f=fopen(filename,"r");
    if (!f) { perror(filename); return 0; }
    char *header=read_line(f);
    if (!header) { fclose(f); return 0; }
    size_t columns=1;
    for (char *p=header;*p;++p) if (*p==',') ++columns;
    size_t *indices=calloc(columns,sizeof(size_t));
    int *types=calloc(columns,sizeof(int));
    if (!indices || !types) { free(header);free(indices);free(types);fclose(f);return 0; }
    char *cursor=header;
    int ok=1;
    for (size_t j=0;j<columns;++j) {
        char *name=field(&cursor);
        if (!strcmp(name,"temperature")) types[j]=0;
        else if (!strcmp(name,"pressure")) types[j]=1;
        else if (chem_species_index(m,name,&indices[j])==CHEM_OK) types[j]=2;
        else if (chem_custom_index(m,name,&indices[j])==CHEM_OK) types[j]=3;
        else { fprintf(stderr,"Unknown CSV column: %s\n",name); ok=0; break; }
    }
    for (size_t cell=0;ok && cell<cells;++cell) {
        char *line=read_line(f);
        if (!line) { fprintf(stderr,"Missing CSV row %zu\n",cell+1); ok=0; break; }
        cursor=line;
        ChemState *s=groups[cell/group_size].state;
        size_t c=cell%group_size;
        for (size_t j=0;j<columns;++j) {
            char *value=field(&cursor),*end;
            if (!value) { ok=0; break; }
            errno=0; double x=strtod(value,&end);
            if (*end || end==value || errno || !isfinite(x)) { ok=0; break; }
            switch(types[j]) {
            case 0: s->conditions[c].temperature=x; break;
            case 1: s->conditions[c].pressure=x; break;
            case 2: s->values[indices[j]*s->cells+c]=x; break;
            case 3: s->custom[indices[j]*s->cells+c]=x; break;
            }
        }
        if (cursor) ok=0;
        if (!(s->conditions[c].temperature>0) || !(s->conditions[c].pressure>0)) ok=0;
        chem_ideal_air_density(&s->conditions[c]);
        free(line);
    }
    free(header);free(indices);free(types);fclose(f);return ok;
}
static void add_stats(ChemSolverStats *a,const ChemSolverStats *b)
{
    a->internal_steps+=b->internal_steps; a->accepted+=b->accepted; a->rejected+=b->rejected;
    a->function_calls+=b->function_calls; a->jacobian_updates+=b->jacobian_updates;
    a->decompositions+=b->decompositions; a->solves+=b->solves;
}
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
    if (!groups || ts1_create(&m) || chem_mechanism_reorder(m)) goto finish;
    for (size_t g=0;g<ng;++g) {
        size_t count=cells-g*group_size;
        if (count>group_size) count=group_size;
        if (chem_state_create(m,count,&groups[g].state) || ts1_initialize(groups[g].state) ||
            chem_solver_create(m,count,&groups[g].solver)) goto finish;
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
