#include "chem_specialized.h"
#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bench_io.h"

extern const chem_specialized ts1_specialized;
#ifndef CHEM_GROUP_SIZE
#define CHEM_GROUP_SIZE 32
#endif

typedef struct {
    char *cells_path, *light_path;
    size_t steps, snapshot_every;
    double dt;
} scenario;
typedef struct { double time, factor; } light_point;

static int real_number(const char *text, double *value)
{
    char *end;
    errno=0; *value=strtod(text,&end);
    return end!=text && !*end && !errno && isfinite(*value);
}
static char *trim(char *text)
{
    while (isspace((unsigned char)*text)) ++text;
    size_t n=strlen(text);
    while (n && isspace((unsigned char)text[n-1])) text[--n]='\0';
    return text;
}
static char *relative_path(const char *scenario_path, const char *name)
{
    const char *slash=strrchr(scenario_path,'/');
    size_t prefix=name[0]=='/' || !slash ? 0 : (size_t)(slash-scenario_path)+1;
    size_t length=strlen(name);
    if (length>SIZE_MAX-prefix-1) return NULL;
    char *path=malloc(prefix+length+1);
    if (path) { memcpy(path,scenario_path,prefix); memcpy(path+prefix,name,length+1); }
    return path;
}
static int read_scenario(const char *path, scenario *sc)
{
    FILE *file=fopen(path,"r");
    if (!file) return 0;
    unsigned seen=0;
    int ok=1;
    char *line;
    while ((line=read_line(file))) {
        char *comment=strchr(line,'#');
        if (comment) *comment='\0';
        char *key=trim(line), *equal=strchr(key,'=');
        if (!*key) { free(line); continue; }
        if (!equal) { free(line); ok=0; break; }
        *equal='\0'; key=trim(key);
        char *value=trim(equal+1);
        unsigned bit=0;
        if (!strcmp(key,"cells")) {
            bit=1; if (!(seen&bit)) sc->cells_path=relative_path(path,value);
            ok=*value && sc->cells_path;
        } else if (!strcmp(key,"photolysis")) {
            bit=2; if (!(seen&bit)) sc->light_path=relative_path(path,value);
            ok=*value && sc->light_path;
        } else if (!strcmp(key,"dt")) { bit=4; ok=real_number(value,&sc->dt) && sc->dt>0; }
        else if (!strcmp(key,"steps")) { bit=8; ok=integer(value,&sc->steps,0); }
        else if (!strcmp(key,"snapshot_every")) { bit=16; ok=integer(value,&sc->snapshot_every,0); }
        else ok=0;
        if (seen&bit) ok=0;
        seen|=bit; free(line);
        if (!ok) break;
    }
    ok=ok && !ferror(file) && (seen&15)==15;
    fclose(file);
    return ok && isfinite(sc->dt*(double)sc->steps);
}
static size_t count_cells(const char *path)
{
    FILE *file=fopen(path,"r");
    if (!file) return 0;
    char *line=read_line(file); free(line);
    size_t count=0;
    while ((line=read_line(file))) { if (*line) ++count; free(line); }
    if (ferror(file)) count=0;
    fclose(file); return count;
}
static int read_light(const char *path, light_point **points, size_t *count)
{
    FILE *file=fopen(path,"r");
    if (!file) return 0;
    char *line=read_line(file);
    int ok=line && !strcmp(line,"time_s,factor");
    free(line);
    size_t capacity=0;
    while (ok && (line=read_line(file))) {
        if (!*line) { free(line); continue; }
        char *fields[2]; double t,f;
        ok=split_csv(line,fields,2)==2 && real_number(fields[0],&t) && real_number(fields[1],&f) && f>=0;
        if (ok && *count && t<=(*points)[*count-1].time) ok=0;
        if (ok && *count==capacity) {
            size_t next=capacity ? capacity*2 : 128;
            if (next<capacity || next>SIZE_MAX/sizeof(light_point)) ok=0;
            else {
                light_point *grown=realloc(*points,next*sizeof(light_point));
                if (!grown) ok=0;
                else { *points=grown; capacity=next; }
            }
        }
        if (ok) (*points)[(*count)++]=(light_point){t,f};
        free(line);
    }
    ok=ok && *count && !ferror(file);
    fclose(file); return ok;
}
static double light_at(const light_point *p, size_t count, double t)
{
    if (t<=p[0].time) return p[0].factor;
    if (t>=p[count-1].time) return p[count-1].factor;
    size_t lo=0,hi=count-1;
    while (hi-lo>1) { size_t mid=lo+(hi-lo)/2; if (p[mid].time>t) hi=mid; else lo=mid; }
    double weight=(t-p[lo].time)/(p[hi].time-p[lo].time);
    return (1.0-weight)*p[lo].factor+weight*p[hi].factor;
}
static int snapshot(FILE *file, chem_state *states, size_t groups, size_t step)
{
    size_t row=0;
    for (size_t g=0; g<groups; ++g) {
        chem_state *s=&states[g];
        for (size_t c=0; c<s->cells; ++c,++row) {
            fprintf(file,"%zu,%zu,%.17g,%.17g",step,row,s->conditions[c].temperature,s->conditions[c].pressure);
            for (size_t i=0; i<s->model->n_species; ++i) fprintf(file,",%.17g",s->y[i*s->cells+c]);
            fputc('\n',file);
        }
    }
    return !ferror(file);
}

int main(int argc, char **argv)
{
    size_t threads=(size_t)omp_get_num_procs(), group_size=CHEM_GROUP_SIZE;
    double rtol=1e-6, atol=0;
    int first=1, generic=0;
    while (first<argc && argv[first][0]=='-') {
        const char *key=argv[first++];
        if (!strcmp(key,"--generic")) { generic=1; continue; }
        if (first==argc) goto usage;
        const char *value=argv[first++];
        if (!strcmp(key,"--threads")) { if (!integer(value,&threads,0)) goto usage; }
        else if (!strcmp(key,"--group-size")) { if (!integer(value,&group_size,0)) goto usage; }
        else if (!strcmp(key,"--rtol")) { if (!real_number(value,&rtol) || rtol<=0) goto usage; }
        else if (!strcmp(key,"--atol")) { if (!real_number(value,&atol) || atol<=0) goto usage; }
        else goto usage;
    }
    if (argc-first!=1 || threads>INT_MAX || !group_size) goto usage;
    scenario sc={0}; chem_model model={0}; chem_state *states=NULL;
    chem_result *results=NULL; double **base=NULL;
    light_point *light=NULL; size_t light_count=0, groups=0;
    FILE *snapshots=NULL; int success=0;
    if (!read_scenario(argv[first],&sc)) { fprintf(stderr,"Invalid scenario\n"); goto done; }
    size_t cells=count_cells(sc.cells_path);
    if (!cells || !read_light(sc.light_path,&light,&light_count)) { fprintf(stderr,"Invalid scenario inputs\n"); goto done; }
    chem_status status=chem_model_create(&model,ts1_specialized.mechanism,CHEM_LU_MOZART);
    if (status!=CHEM_OK) { fprintf(stderr,"Model: %s\n",chem_status_string(status)); goto done; }
    groups=cells/group_size+(cells%group_size!=0);
    states=calloc(groups,sizeof(*states)); results=calloc(groups,sizeof(*results)); base=calloc(groups,sizeof(*base));
    if (!states || !results || !base) goto done;
    for (size_t g=0; g<groups; ++g) {
        size_t count=cells-g*group_size; if (count>group_size) count=group_size;
        status=chem_state_create(&states[g],&model,count);
        if (status!=CHEM_OK) goto done;
        states[g].rtol=rtol;
        if (atol) for (size_t i=0; i<model.n_species; ++i) states[g].atol[i]=atol;
        base[g]=malloc((model.n_custom+1)*count*sizeof(double));
        if (!base[g]) goto done;
    }
    if (!read_cells(sc.cells_path,states,group_size,cells)) goto done;
    for (size_t g=0; g<groups; ++g) memcpy(base[g],states[g].custom,model.n_custom*states[g].cells*sizeof(double));
    snapshots=fopen("ts1_snapshots.csv","w");
    if (!snapshots) goto done;
    fprintf(snapshots,"step,cell,temperature,pressure");
    for (size_t i=0; i<model.n_species; ++i) fprintf(snapshots,",%s",model.mechanism->species[model.var_to_spec[i]].name);
    fputc('\n',snapshots);
    omp_set_num_threads((int)threads);
    double elapsed=0;
    chem_ros_parameters parameters=chem_ros_parameters_default(CHEM_ROS3);
    size_t interval=sc.snapshot_every ? sc.snapshot_every : sc.steps;
    for (size_t begin=0; begin<sc.steps;) {
        size_t end=interval>sc.steps-begin ? sc.steps : begin+interval;
        double start=omp_get_wtime();
        #pragma omp parallel for schedule(static)
        for (size_t g=0; g<groups; ++g) {
            chem_state *s=&states[g];
            for (size_t step=begin; step<end; ++step) {
                double factor=light_at(light,light_count,(double)step*sc.dt);
                for (size_t j=0; j<model.n_custom; ++j) if (model.custom_names[j][0]=='j')
                    for (size_t c=0; c<s->cells; ++c) s->custom[j*s->cells+c]=base[g][j*s->cells+c]*factor;
                if (generic) chem_calculate_rates(s); else ts1_specialized.rates(s);
                chem_result r=generic ? chem_ros_solve(s,sc.dt,&parameters,NULL,NULL) : ts1_specialized.solve(s,sc.dt);
                chem_state_clamp(s);
                chem_stats_add(&results[g].stats,&r.stats);
                if (r.status!=CHEM_OK) { results[g].status=r.status; break; }
            }
        }
        elapsed+=(omp_get_wtime()-start)*1000;
        for (size_t g=0; g<groups; ++g) if (results[g].status!=CHEM_OK) {
            fprintf(stderr,"Group %zu: %s\n",g,chem_status_string(results[g].status)); goto done;
        }
        if (!snapshot(snapshots,states,groups,end)) goto done;
        begin=end;
    }
    chem_stats stats={0};
    for (size_t g=0; g<groups; ++g) chem_stats_add(&stats,&results[g].stats);
    printf("elapsed_ms=%.6f\n",elapsed);
    printf("internal_steps=%" PRIu64 " accepted=%" PRIu64 " rejected=%" PRIu64 " function_calls=%" PRIu64
           " jacobian_updates=%" PRIu64 " decompositions=%" PRIu64 " solves=%" PRIu64 "\n",stats.internal_steps,
           stats.accepted,stats.rejected,stats.function_calls,stats.jacobian_updates,stats.decompositions,stats.solves);
    success=write_final(states,groups);
done:
    if (snapshots && fclose(snapshots)) success=0;
    for (size_t g=0; g<groups; ++g) { if (states) chem_state_destroy(&states[g]); if (base) free(base[g]); }
    free(states); free(results); free(base); free(light);
    free(sc.cells_path); free(sc.light_path); chem_model_destroy(&model);
    return success ? 0 : 1;
usage:
    fprintf(stderr,"Usage: %s [--threads N] [--group-size G] [--rtol R] [--atol A] [--generic] SCENARIO\n",argv[0]);
    return 2;
}
