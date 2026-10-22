/* SPDX-License-Identifier: Apache-2.0 */
#include "ts1_driver.h"
#include <ctype.h>

typedef struct {
    char *cells, *photolysis;
    size_t steps, snapshot_every;
    double dt;
} Scenario;
typedef struct { double time, factor; } Light;

static char *trim(char *s)
{
    while (isspace((unsigned char)*s)) ++s;
    char *end=s+strlen(s);
    while (end>s && isspace((unsigned char)end[-1])) --end;
    *end=0; return s;
}
static double positive_real(const char *s)
{
    char *end; errno=0; double x=strtod(s,&end);
    return errno || end==s || *end || !isfinite(x) || x<=0 ? 0 : x;
}
static char *relative_path(const char *scenario,const char *name)
{
    const char *slash=strrchr(scenario,'/');
    size_t prefix=name[0]=='/' || !slash ? 0 : (size_t)(slash-scenario)+1;
    if (prefix>SIZE_MAX-strlen(name)-1) return NULL;
    char *path=malloc(prefix+strlen(name)+1);
    if (path) { memcpy(path,scenario,prefix); strcpy(path+prefix,name); }
    return path;
}
static int read_scenario(const char *path,Scenario *s)
{
    FILE *f=fopen(path,"r"); if (!f) { perror(path); return 0; }
    int ok=1; char *line;
    while ((line=read_line(f))) {
        char *comment=strchr(line,'#'); if (comment) *comment=0;
        char *key=trim(line);
        if (!*key) { free(line); continue; }
        char *eq=strchr(key,'=');
        if (!eq) { free(line); ok=0; break; }
        *eq=0; key=trim(key); char *value=trim(eq+1);
        if (!strcmp(key,"cells") && !s->cells) s->cells=relative_path(path,value);
        else if (!strcmp(key,"photolysis") && !s->photolysis) s->photolysis=relative_path(path,value);
        else if (!strcmp(key,"steps") && !s->steps) { s->steps=positive_integer(value); ok=s->steps!=0; }
        else if (!strcmp(key,"dt") && !s->dt) { s->dt=positive_real(value); ok=s->dt!=0; }
        else if (!strcmp(key,"snapshot_every") && !s->snapshot_every) { s->snapshot_every=positive_integer(value); ok=s->snapshot_every!=0; }
        else ok=0;
        free(line); if (!ok) break;
    }
    if (ferror(f)) ok=0;
    fclose(f);
    return ok && s->cells && s->photolysis && s->steps && s->dt && isfinite(s->steps*s->dt);
}
static size_t count_cells(const char *path)
{
    FILE *f=fopen(path,"r"); if (!f) { perror(path); return 0; }
    char *line=read_line(f); free(line);
    size_t count=0;
    while ((line=read_line(f))) { ++count; free(line); }
    if (ferror(f)) count=0;
    fclose(f); return count;
}
static int read_light(const char *path,Light **out,size_t *count)
{
    FILE *f=fopen(path,"r"); if (!f) { perror(path); return 0; }
    char *line=read_line(f); int ok=line && !strcmp(line,"time_s,factor"); free(line);
    size_t capacity=0;
    while (ok && (line=read_line(f))) {
        char *end; errno=0;
        double time=strtod(line,&end), factor=0;
        ok=end!=line && *end==',' && isfinite(time) && !errno;
        if (ok) {
            char *value=end+1; factor=strtod(value,&end);
            ok=end!=value && !*end && isfinite(factor) && factor>=0 && !errno;
        }
        free(line);
        if (!ok || (*count && time<=(*out)[*count-1].time)) { ok=0; break; }
        if (*count==capacity) {
            if (capacity>SIZE_MAX/2/sizeof(Light)) { ok=0; break; }
            capacity=capacity?capacity*2:256;
            Light *next=realloc(*out,capacity*sizeof(Light));
            if (!next) { ok=0; break; } *out=next;
        }
        (*out)[(*count)++]=(Light){time,factor};
    }
    if (ferror(f)) ok=0;
    fclose(f); return ok && *count;
}
static double light_factor(const Light *light,size_t count,double t)
{
    if (t<=light[0].time) return light[0].factor;
    if (t>=light[count-1].time) return light[count-1].factor;
    size_t lo=0,hi=count-1;
    while (hi-lo>1) { size_t mid=lo+(hi-lo)/2; if (light[mid].time<=t) lo=mid; else hi=mid; }
    double fraction=(t-light[lo].time)/(light[hi].time-light[lo].time);
    return (1.0-fraction)*light[lo].factor+fraction*light[hi].factor;
}
static void write_header(FILE *f,const ChemMechanism *m,int snapshots)
{
    if (snapshots) fprintf(f,"step,cell,");
    fprintf(f,"temperature,pressure");
    for (size_t v=0;v<chem_species_count(m);++v) fprintf(f,",%s",chem_species_name(m,v));
    fputc('\n',f);
}
static int write_state(FILE *f,const Group *groups,size_t cells,size_t group_size,const ChemMechanism *m,size_t step)
{
    for (size_t cell=0;cell<cells;++cell) {
        const ChemState *s=groups[cell/group_size].state; size_t c=cell%group_size;
        if (step) fprintf(f,"%zu,%zu,",step,cell);
        fprintf(f,"%.17g,%.17g",s->conditions[c].temperature,s->conditions[c].pressure);
        for (size_t v=0;v<chem_species_count(m);++v) {
            double value=s->values[v*s->cells+c];
            if (!isfinite(value) || value < -1e-20) return 0;
            fprintf(f,",%.17g",value);
        }
        fputc('\n',f);
    }
    return !ferror(f);
}
int main(int argc,char **argv)
{
    int arg=1,threads=omp_get_max_threads();
    size_t group_size=CHEM_GROUP_SIZE;
    double rtol=0,atol=0;
    const char *setting=getenv("CHEM_GROUP_SIZE");
    if (setting) group_size=positive_integer(setting);
    while (arg+1<argc && argv[arg][0]=='-') {
        const char *option=argv[arg++],*value=argv[arg++];
        if (!strcmp(option,"--threads")) {
            size_t n=positive_integer(value); if (!n || n>INT_MAX) goto usage; threads=(int)n;
        } else if (!strcmp(option,"--group-size")) { group_size=positive_integer(value); }
        else if (!strcmp(option,"--rtol")) { rtol=positive_real(value); if (!rtol) goto usage; }
        else if (!strcmp(option,"--atol")) { atol=positive_real(value); if (!atol) goto usage; }
        else goto usage;
    }
    if (arg+1!=argc || !group_size) goto usage;
    Scenario scenario={0}; Light *light=NULL; size_t light_count=0,ng=0;
    ChemMechanism *m=NULL; Group *groups=NULL; double **base=NULL;
    FILE *snapshots=NULL,*final=NULL; int success=0;
    if (!read_scenario(argv[arg],&scenario) || !read_light(scenario.photolysis,&light,&light_count)) goto finish;
    size_t cells=count_cells(scenario.cells); if (!cells) goto finish;
    if (group_size>cells) group_size=cells;
    ng=1+(cells-1)/group_size;
    groups=calloc(ng,sizeof(Group)); base=calloc(ng,sizeof(double *));
    if (!groups || !base || create_mechanism(&m)) goto finish;
    for (size_t g=0;g<ng;++g) {
        size_t count=cells-g*group_size; if (count>group_size) count=group_size;
        if (chem_state_create(m,count,&groups[g].state) || ts1_initialize(groups[g].state) ||
            create_solver(m,count,&groups[g].solver)) goto finish;
        ChemState *s=groups[g].state;
        if (rtol) s->relative_tolerance=rtol;
        if (atol) for (size_t v=0;v<chem_species_count(m);++v) s->absolute_tolerance[v]=atol;
        base[g]=malloc(chem_custom_count(m)*count*sizeof(double));
        if (!base[g]) goto finish;
    }
    if (!load_cells(scenario.cells,groups,cells,group_size,m)) goto finish;
    for (size_t g=0;g<ng;++g) memcpy(base[g],groups[g].state->custom,chem_custom_count(m)*groups[g].state->cells*sizeof(double));
    snapshots=fopen("ts1_snapshots.csv","w"); if (!snapshots) goto finish;
    write_header(snapshots,m,1);
    omp_set_num_threads(threads);
    if (!getenv("OMP_SCHEDULE")) omp_set_schedule(omp_sched_static,0);
    ChemRosenbrockParameters parameters=chem_rosenbrock_parameters(CHEM_ROS3);
    double elapsed=0;
    for (size_t step=0;step<scenario.steps;++step) {
        double start=omp_get_wtime();
        double factor=light_factor(light,light_count,step*scenario.dt);
        #pragma omp parallel for schedule(runtime)
        for (size_t g=0;g<ng;++g) {
            ChemState *s=groups[g].state;
            for (size_t j=0;j<chem_custom_count(m);++j) if (chem_custom_name(m,j)[0]=='j')
                for (size_t c=0;c<s->cells;++c) s->custom[j*s->cells+c]=base[g][j*s->cells+c]*factor;
            chem_update_rates(s);
            ChemSolverStats stats;
            groups[g].status=chem_solve(groups[g].solver,s,scenario.dt,&parameters,&stats);
            add_stats(&groups[g].stats,&stats);
        }
        elapsed+=1000*(omp_get_wtime()-start);
        for (size_t g=0;g<ng;++g) if (groups[g].status) {
            fprintf(stderr,"Step %zu, group %zu: %s\n",step+1,g,chem_status_string(groups[g].status)); goto finish;
        }
        if ((scenario.snapshot_every && (step+1)%scenario.snapshot_every==0) || step+1==scenario.steps)
            if (!write_state(snapshots,groups,cells,group_size,m,step+1)) goto finish;
    }
    ChemSolverStats stats={0};
    for (size_t g=0;g<ng;++g) add_stats(&stats,&groups[g].stats);
    printf("elapsed_ms=%.6f\n",elapsed);
    printf("internal_steps=%zu accepted=%zu rejected=%zu function_calls=%zu jacobian_updates=%zu decompositions=%zu solves=%zu\n",
           stats.internal_steps,stats.accepted,stats.rejected,stats.function_calls,stats.jacobian_updates,stats.decompositions,stats.solves);
    final=fopen("ts1_final.csv","w"); if (!final) goto finish;
    write_header(final,m,0);
    success=write_state(final,groups,cells,group_size,m,0);
finish:
    if (final && fclose(final)) success=0;
    if (snapshots && fclose(snapshots)) success=0;
    if (groups) for (size_t g=0;g<ng;++g) { chem_state_free(groups[g].state); chem_solver_free(groups[g].solver); }
    if (base) for (size_t g=0;g<ng;++g) free(base[g]);
    free(groups); free(base); free(light); free(scenario.cells); free(scenario.photolysis); chem_mechanism_free(m);
    if (!success) fprintf(stderr,"Scenario failed\n");
    return success?0:1;
usage:
    fprintf(stderr,"Usage: %s [--threads N] [--group-size N] [--rtol R] [--atol A] SCENARIO\n",argv[0]);
    return 1;
}
