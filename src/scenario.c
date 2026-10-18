#include "chem_adaptive.h"
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

extern const chem_mechanism ts1_mechanism;
typedef struct {
    char *cells_path,*photo_path;
    double dt,*times,*factors;
    size_t steps,snapshot_every,cells,n_photo;
} scenario;

static char *trim(char *s)
{
    while (isspace((unsigned char)*s)) ++s;
    char *end=s+strlen(s);
    while (end>s && isspace((unsigned char)end[-1])) --end;
    *end=0; return s;
}
static int real_value(const char *s,double *v)
{
    char *end; errno=0; *v=strtod(s,&end);
    return end!=s && !*end && errno!=ERANGE && isfinite(*v);
}
static char *relative_path(const char *file,const char *value)
{
    const char *slash=strrchr(file,'/');
    size_t prefix=value[0]=='/' || !slash ? 0 : (size_t)(slash-file)+1;
    if (strlen(value)>SIZE_MAX-prefix-1) return NULL;
    char *path=malloc(prefix+strlen(value)+1);
    if (path) { memcpy(path,file,prefix); strcpy(path+prefix,value); }
    return path;
}
static void destroy_scenario(scenario *a)
{ free(a->cells_path); free(a->photo_path); free(a->times); free(a->factors); }

static int load_scenario(const char *path,scenario *a)
{
    FILE *f=fopen(path,"r"); if (!f) return 0;
    char *line; unsigned seen=0; int ok=1;
    while ((line=read_line(f))) {
        char *comment=strchr(line,'#'); if (comment) *comment=0;
        char *key=trim(line),*value=strchr(key,'=');
        unsigned bit=0;
        if (!*key) { free(line); continue; }
        if (!value) { free(line); ok=0; break; }
        *value++=0; key=trim(key); value=trim(value);
        if (!strcmp(key,"cells")) { bit=1; if (!(seen&bit)) a->cells_path=relative_path(path,value); }
        else if (!strcmp(key,"photolysis")) { bit=2; if (!(seen&bit)) a->photo_path=relative_path(path,value); }
        else if (!strcmp(key,"dt")) { bit=4; ok=real_value(value,&a->dt) && a->dt>0; }
        else if (!strcmp(key,"steps")) { bit=8; ok=integer(value,&a->steps,0); }
        else if (!strcmp(key,"snapshot_every")) { bit=16; ok=integer(value,&a->snapshot_every,1); }
        else ok=0;
        if (seen&bit) ok=0;
        seen|=bit; free(line); if (!ok) break;
    }
    fclose(f);
    if (!ok || (seen&15)!=15 || !a->cells_path || !a->photo_path || !isfinite(a->dt*a->steps)) return 0;
    f=fopen(a->cells_path,"r"); if (!f) return 0;
    line=read_line(f); if (!line) { fclose(f); return 0; } free(line);
    while ((line=read_line(f))) { if (*trim(line)) ++a->cells; free(line); }
    fclose(f); if (!a->cells) return 0;
    f=fopen(a->photo_path,"r"); if (!f) return 0;
    line=read_line(f);
    if (!line || strcmp(trim(line),"time_s,factor")) { free(line); fclose(f); return 0; }
    free(line);
    size_t cap=0;
    while ((line=read_line(f))) {
        char *fields[2]; double t,v;
        if (!*trim(line)) { free(line); continue; }
        if (split_csv(line,fields,2)!=2 || !real_value(fields[0],&t) || !real_value(fields[1],&v) ||
            v<0 || (a->n_photo && t<=a->times[a->n_photo-1])) { free(line); ok=0; break; }
        if (a->n_photo==cap) {
            cap=cap ? cap*2 : 64;
            double *new_times=realloc(a->times,cap*sizeof(double));
            if (!new_times) { free(line); ok=0; break; } a->times=new_times;
            double *new_factors=realloc(a->factors,cap*sizeof(double));
            if (!new_factors) { free(line); ok=0; break; } a->factors=new_factors;
        }
        a->times[a->n_photo]=t; a->factors[a->n_photo++]=v; free(line);
    }
    fclose(f); return ok && a->n_photo;
}
static double photo_factor(const scenario *a,double t)
{
    if (t<=a->times[0]) return a->factors[0];
    for (size_t i=1;i<a->n_photo;++i) if (t<a->times[i]) {
        double x=(t-a->times[i-1])/(a->times[i]-a->times[i-1]);
        return a->factors[i-1]+x*(a->factors[i]-a->factors[i-1]);
    }
    return a->factors[a->n_photo-1];
}
static int snapshot(FILE *f,chem_state *states,size_t groups,size_t step)
{
    const chem_model *m=states[0].model;
    if (!step) {
        fprintf(f,"step,cell,temperature,pressure");
        for (size_t i=0;i<m->n_species;++i) fprintf(f,",%s",m->mechanism->species[m->var_to_spec[i]].name);
        fputc('\n',f);
    } else {
        size_t row=0;
        for (size_t g=0;g<groups;++g) for (size_t c=0;c<states[g].cells;++c) {
            chem_state *s=states+g;
            fprintf(f,"%zu,%zu,%.17g,%.17g",step,row++,s->conditions[c].temperature,s->conditions[c].pressure);
            for (size_t i=0;i<m->n_species;++i) fprintf(f,",%.17g",s->y[i*s->cells+c]);
            fputc('\n',f);
        }
    }
    return !ferror(f);
}

static int option(const char *arg,chem_adaptive_parameters *p,chem_ros_parameters *ros,
    double *rtol,double *atol,size_t *group,int *generic)
{
    const char *eq=strchr(arg,'='); if (!eq || eq==arg || (size_t)(eq-arg)>=64) return 0;
    char key[64]; memcpy(key,arg,(size_t)(eq-arg)); key[eq-arg]=0;
    double v; size_t n;
    if (!strcmp(key,"rtol")) return real_value(eq+1,rtol) && *rtol>0;
    if (!strcmp(key,"atol")) return real_value(eq+1,atol) && *atol>0;
    if (!strcmp(key,"h-start")) { if (!real_value(eq+1,&v) || v<0) return 0; p->h_start=v; ros->h_start=v; return 1; }
    if (!strcmp(key,"h-max")) { if (!real_value(eq+1,&v) || v<0) return 0; p->h_max=v; ros->h_max=v; return 1; }
    if (!strcmp(key,"safety")) { if (!real_value(eq+1,&v) || v<=0 || v>=1) return 0; p->safety=v; return 1; }
    if (!strcmp(key,"growth")) { if (!real_value(eq+1,&v) || v<1) return 0; p->growth=v; return 1; }
    if (!strcmp(key,"reuse-window")) { if (!real_value(eq+1,&v) || v<0 || v>=1) return 0; p->reuse_window=v; return 1; }
    if (!strcmp(key,"newton-tol")) { if (!real_value(eq+1,&v) || v<=0 || v>=1) return 0; p->newton_tolerance=v; return 1; }
    if (!integer(eq+1,&n,1)) return 0;
    if (!strcmp(key,"max-steps") && n) p->max_steps=ros->max_steps=n;
    else if (!strcmp(key,"jacobian-age") && n) p->jacobian_age=n;
    else if (!strcmp(key,"newton-iterations") && n) p->newton_iterations=n;
    else if (!strcmp(key,"order") && n>=2 && n<=4) p->extrap_order=n;
    else if (!strcmp(key,"group-size") && n) *group=n;
    else if (!strcmp(key,"norm") && n<=2) p->error_norm=(int)n;
    else if (!strcmp(key,"controller") && n<=1) p->controller=(int)n;
    else if (!strcmp(key,"clip-stages") && n<=1) p->clip_stages=(int)n;
    else if (!strcmp(key,"warm-start") && n<=1) p->warm_start=(int)n;
    else if (!strcmp(key,"adaptive-order") && n<=1) p->adaptive_order=(int)n;
    else if (!strcmp(key,"generic") && n<=1) *generic=(int)n;
    else return 0;
    return 1;
}

int main(int argc,char **argv)
{
    size_t threads=(size_t)omp_get_num_procs(),group=32;
    const char *method="micm",*path=NULL; int adaptive=0,generic=0;
    chem_adaptive_parameters p=chem_adaptive_defaults();
    chem_ros_parameters ros=chem_ros_parameters_default(CHEM_ROS3);
    double rtol=0,atol=0;
    /* Method is read first so option order cannot change its defaults. */
    for (int i=1;i<argc;++i) if (!strcmp(argv[i],"--method")) {
        if (++i==argc) goto usage;
        method=argv[i];
    }
    const char *names[]={"ros2","ros3","ros4","rodas3","rodas4"};
    const char *old_names[]={"micm2","micm","micm4","micm4dae","micm6"};
    int found=0;
    for (size_t i=0;i<5;++i) {
        if (!strcmp(method,names[i])) { adaptive=1; p.ros=(chem_ros_method)i; found=1; }
        if (!strcmp(method,old_names[i])) { ros=chem_ros_parameters_default((chem_ros_method)i); found=1; }
    }
    if (!strcmp(method,"extrap")) { adaptive=1; p.method=CHEM_ADAPT_EXTRAP; found=1; }
    if (!strcmp(method,"sdirk")) { adaptive=1; p.method=CHEM_ADAPT_SDIRK; found=1; }
    if (!found) goto usage;
    for (int i=1;i<argc;++i) {
        if (!strcmp(argv[i],"--method")) { ++i; continue; }
        if (!strcmp(argv[i],"--threads")) {
            if (++i==argc || !integer(argv[i],&threads,0) || threads>INT_MAX) goto usage;
        } else if (!strcmp(argv[i],"--option")) {
            if (++i==argc || !option(argv[i],&p,&ros,&rtol,&atol,&group,&generic)) goto usage;
        } else if (argv[i][0]=='-' || path) goto usage;
        else path=argv[i];
    }
    if (!path) goto usage;
    scenario a={0}; chem_model model={0}; chem_state *states=NULL;
    chem_result *results=NULL; chem_adaptive_workspace *work=NULL; double **base=NULL;
    FILE *snap=NULL; size_t groups=0; int ok=0;
    if (!load_scenario(path,&a)) { fprintf(stderr,"Invalid scenario: %s\n",path); goto done; }
    chem_status rc=chem_model_create(&model,&ts1_mechanism,CHEM_LU_MOZART);
    if (rc!=CHEM_OK) { fprintf(stderr,"Model: %s\n",chem_status_string(rc)); goto done; }
    if (generic) { model.forcing=chem_forcing; model.jacobian=chem_jacobian; }
    groups=a.cells/group+(a.cells%group!=0);
    states=calloc(groups,sizeof(*states)); results=calloc(groups,sizeof(*results));
    work=calloc(groups,sizeof(*work)); base=calloc(groups,sizeof(*base));
    if (!states || !results || !work || !base) goto done;
    for (size_t g=0;g<groups;++g) {
        size_t n=a.cells-g*group; if (n>group) n=group;
        if (chem_state_create(states+g,&model,n)!=CHEM_OK) goto done;
        if (adaptive && chem_adaptive_create(work+g,states+g)!=CHEM_OK) goto done;
        if (rtol || adaptive) states[g].rtol=rtol ? rtol : 1e-4;
        if (atol || adaptive) for (size_t i=0;i<model.n_species;++i) states[g].atol[i]=atol ? atol : 1e-18;
        base[g]=malloc(model.n_custom*n*sizeof(double)); if (!base[g]) goto done;
    }
    if (!read_cells(a.cells_path,states,group,a.cells)) goto done;
    for (size_t g=0;g<groups;++g)
        memcpy(base[g],states[g].custom,model.n_custom*states[g].cells*sizeof(double));
    if (a.snapshot_every) {
        snap=fopen("ts1_snapshots.csv","w");
        if (!snap || !snapshot(snap,states,groups,0)) goto done;
    }
    omp_set_num_threads((int)threads);
    double elapsed=0; size_t completed=0;
    for (size_t step=0;step<a.steps;++step) {
        double start=omp_get_wtime();
        double f=photo_factor(&a,step*a.dt);
        #pragma omp parallel for schedule(static)
        for (size_t g=0;g<groups;++g) {
            chem_state *s=states+g;
            for (size_t j=0;j<model.n_custom;++j) if (model.custom_names[j][0]=='j')
                for (size_t c=0;c<s->cells;++c) s->custom[j*s->cells+c]=base[g][j*s->cells+c]*f;
            chem_calculate_rates(s);
            chem_result r=adaptive ? chem_adaptive_solve(s,a.dt,&p,work+g) : chem_ros_solve(s,a.dt,&ros,NULL,NULL);
            if (!adaptive) chem_state_clamp(s);
            chem_stats_add(&results[g].stats,&r.stats); results[g].status=r.status;
            for (size_t i=0;i<model.n_species*s->cells;++i) if (!isfinite(s->y[i]) || s->y[i]<-1e-20)
                results[g].status=CHEM_NAN;
        }
        elapsed+=(omp_get_wtime()-start)*1000;
        int failed=0;
        for (size_t g=0;g<groups;++g) if (results[g].status!=CHEM_OK) {
            if (!failed) fprintf(stderr,"Step %zu group %zu: %s\n",step+1,g,chem_status_string(results[g].status));
            failed=1;
        }
        if (failed) break;
        completed=step+1;
        if (snap && (completed%a.snapshot_every==0 || completed==a.steps))
            if (!snapshot(snap,states,groups,completed)) goto done;
    }
    chem_stats stats={0};
    uint64_t jr=0,fr=0,pr=0,nr=0,cv=0,ni=0,nj=0,eo[3]={0}; double mass=0;
    for (size_t g=0;g<groups;++g) {
        chem_stats_add(&stats,&results[g].stats);
        jr+=work[g].jacobian_reuses; fr+=work[g].factor_reuses; pr+=work[g].positivity_rejections;
        nr+=work[g].nonfinite_rejections; cv+=work[g].clipped_values; mass+=work[g].clipped_mass;
        ni+=work[g].nonlinear_iterations; nj+=work[g].nonlinear_rejections;
        for (size_t j=0;j<3;++j) eo[j]+=work[g].extrap_orders[j];
    }
    printf("elapsed_ms=%.6f\n",elapsed);
    printf("internal_steps=%" PRIu64 " accepted=%" PRIu64 " rejected=%" PRIu64 " function_calls=%" PRIu64
        " jacobian_updates=%" PRIu64 " decompositions=%" PRIu64 " solves=%" PRIu64 "\n",stats.internal_steps,
        stats.accepted,stats.rejected,stats.function_calls,stats.jacobian_updates,stats.decompositions,stats.solves);
    printf("jacobian_reuses=%" PRIu64 " factor_reuses=%" PRIu64 " positivity_rejections=%" PRIu64
        " nonfinite_rejections=%" PRIu64 " clipped_values=%" PRIu64 " clipped_mass=%.17g completed_steps=%zu\n",
        jr,fr,pr,nr,cv,mass,completed);
    printf("nonlinear_iterations=%" PRIu64 " nonlinear_rejections=%" PRIu64 "\n",ni,nj);
    printf("extrap_order2=%" PRIu64 " extrap_order3=%" PRIu64 " extrap_order4=%" PRIu64 "\n",eo[0],eo[1],eo[2]);
    if (completed==a.steps) ok=write_final(states,groups);
done:
    if (snap && fclose(snap)) ok=0;
    for (size_t g=0;g<groups;++g) {
        if (states) chem_state_destroy(states+g);
        if (work) chem_adaptive_destroy(work+g);
        if (base) free(base[g]);
    }
    free(states); free(results); free(work); free(base); chem_model_destroy(&model); destroy_scenario(&a);
    return ok ? 0 : 1;
usage:
    fprintf(stderr,"Usage: %s [--threads N] [--method micm|micm2|micm4|micm4dae|micm6|ros2|ros3|ros4|rodas3|rodas4|extrap|sdirk] [--option KEY=VALUE ...] SCENARIO\n",argv[0]);
    return 2;
}
