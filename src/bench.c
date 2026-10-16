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

extern const chem_mechanism ts1_mechanism;
#ifndef CHEM_GROUP_SIZE
#define CHEM_GROUP_SIZE 32
#endif

static char *read_line(FILE *file)
{
    size_t length=0, capacity=4096;
    char *s=malloc(capacity);
    if (!s) return NULL;
    int c;
    while ((c=fgetc(file))!=EOF && c!='\n') {
        if (length+1==capacity) {
            if (capacity>SIZE_MAX/2) { free(s); return NULL; }
            capacity*=2;
            char *next=realloc(s,capacity);
            if (!next) { free(s); return NULL; }
            s=next;
        }
        s[length++]=(char)c;
    }
    if (!length && c==EOF) { free(s); return NULL; }
    if (length && s[length-1]=='\r') --length;
    s[length]='\0'; return s;
}

/* CSV fields are compacted in place; double quotes and escaped quotes are accepted. */
static size_t split_csv(char *s, char **fields, size_t limit)
{
    size_t count=0;
    char *read=s,*write=s;
    for (;;) {
        if (count==limit) return CHEM_NONE;
        fields[count++]=write;
        int quoted=*read=='"';
        if (quoted) ++read;
        while (*read) {
            if (quoted && *read=='"') {
                if (read[1]=='"') { *write++='"'; read+=2; continue; }
                ++read; quoted=0; break;
            }
            if (!quoted && *read==',') break;
            *write++=*read++;
        }
        if (quoted || (*read && *read!=',')) return CHEM_NONE;
        int more=*read==',';
        if (more) ++read;
        *write++='\0';
        if (!more) break;
    }
    return count;
}

static int read_cells(const char *path, chem_state *states, size_t group_size, size_t cells)
{
    FILE *file=fopen(path,"r");
    if (!file) { fprintf(stderr,"Cannot open input: %s\n",path); return 0; }
    const chem_model *m=states[0].model;
    size_t limit=m->n_species+m->n_custom+2;
    char **fields=calloc(limit,sizeof(char *)); size_t *indices=calloc(limit,sizeof(size_t));
    int *kind=calloc(limit,sizeof(int)); char *line=read_line(file);
    int success=0;
    if (!fields || !indices || !kind || !line) goto done;
    size_t columns=split_csv(line,fields,limit);
    if (columns==CHEM_NONE || !columns) goto done;
    for (size_t i=0; i<columns; ++i) {
        for (size_t j=0; j<i; ++j) if (!strcmp(fields[i],fields[j])) goto done;
        if (!strcmp(fields[i],"temperature")) kind[i]=0;
        else if (!strcmp(fields[i],"pressure")) kind[i]=1;
        else if ((indices[i]=chem_species_index(m,fields[i]))!=CHEM_NONE) kind[i]=2;
        else if ((indices[i]=chem_custom_index(m,fields[i]))!=CHEM_NONE) kind[i]=3;
        else { fprintf(stderr,"Unknown input column: %s\n",fields[i]); goto done; }
    }
    free(line); line=NULL;
    for (size_t row=0; row<cells;) {
        line=read_line(file);
        if (!line) { fprintf(stderr,"Input has fewer than %zu cells\n",cells); goto done; }
        if (!*line) { free(line); line=NULL; continue; }
        if (split_csv(line,fields,limit)!=columns) goto done;
        chem_state *s=&states[row/group_size]; size_t c=row%group_size;
        for (size_t i=0; i<columns; ++i) {
            char *end; errno=0; double value=strtod(fields[i],&end);
            while (isspace((unsigned char)*end)) ++end;
            if (fields[i]==end || *end || errno==ERANGE || !isfinite(value)) goto done;
            if (kind[i]==0) s->conditions[c].temperature=value;
            else if (kind[i]==1) s->conditions[c].pressure=value;
            else if (kind[i]==2) s->y[indices[i]*s->cells+c]=value;
            else {
                for (size_t j=indices[i]; j<m->n_custom; ++j)
                    if (!strcmp(m->custom_names[j],m->custom_names[indices[i]])) s->custom[j*s->cells+c]=value;
            }
        }
        if (s->conditions[c].temperature<=0 || s->conditions[c].pressure<=0) goto done;
        chem_ideal_air_density(&s->conditions[c]);
        ++row; free(line); line=NULL;
    }
    success=1;
done:
    free(line); free(fields); free(indices); free(kind); fclose(file);
    if (!success) fprintf(stderr,"Invalid CSV input\n");
    return success;
}

static int write_final(chem_state *states, size_t groups)
{
    FILE *file=fopen("ts1_final.csv","w");
    if (!file) return 0;
    const chem_model *m=states[0].model;
    fprintf(file,"temperature,pressure");
    for (size_t i=0; i<m->n_species; ++i) fprintf(file,",%s",m->mechanism->species[m->var_to_spec[i]].name);
    fputc('\n',file);
    for (size_t g=0; g<groups; ++g) {
        chem_state *s=&states[g];
        for (size_t c=0; c<s->cells; ++c) {
            fprintf(file,"%.17g,%.17g",s->conditions[c].temperature,s->conditions[c].pressure);
            for (size_t i=0; i<m->n_species; ++i) fprintf(file,",%.17g",s->y[i*s->cells+c]);
            fputc('\n',file);
        }
    }
    int success=!ferror(file);
    if (fclose(file)) success=0;
    return success;
}

static int integer(const char *s, size_t *value, int zero_ok)
{
    if (*s=='-' || !*s) return 0;
    char *end; errno=0; uintmax_t v=strtoumax(s,&end,10);
    if (*end || errno || v>SIZE_MAX || (!zero_ok && !v)) return 0;
    *value=(size_t)v; return 1;
}

int main(int argc, char **argv)
{
    size_t threads=(size_t)omp_get_num_procs(), group_size=CHEM_GROUP_SIZE;
    int first=1, generic=0;
    while (first<argc && argv[first][0]=='-') {
        if (!strcmp(argv[first],"--generic")) { generic=1; ++first; continue; }
        size_t value;
        if (first+1>=argc || !integer(argv[first+1],&value,0)) goto usage;
        if (!strcmp(argv[first],"--threads")) threads=value;
        else if (!strcmp(argv[first],"--group-size")) group_size=value;
        else goto usage;
        first+=2;
    }
    size_t cells,steps;
    if (argc-first!=4 || !integer(argv[first],&cells,0) || !integer(argv[first+1],&steps,1) || threads>INT_MAX || !group_size) goto usage;
    char *end; errno=0; double dt=strtod(argv[first+2],&end);
    if (*end || errno || !isfinite(dt) || dt<=0) goto usage;
    omp_set_num_threads((int)threads);
    chem_model model;
    chem_status status=chem_model_create(&model,&ts1_mechanism,CHEM_LU_MOZART);
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
    #pragma omp parallel for schedule(static)
    for (size_t g=0; g<groups; ++g) {
        chem_calculate_rates(&states[g]);
        results[g]=chem_ros_solve(&states[g],dt,&parameters,NULL,NULL);
    }
    for (size_t g=0; g<groups; ++g) if (results[g].status!=CHEM_OK) {
        fprintf(stderr,"Warm-up group %zu: %s\n",g,chem_status_string(results[g].status)); goto done;
    }
    memset(results,0,groups*sizeof(chem_result));
    double start=omp_get_wtime();
    #pragma omp parallel for schedule(static)
    for (size_t g=0; g<groups; ++g) {
        for (size_t step=0; step<steps; ++step) {
            chem_result r=chem_ros_solve(&states[g],dt,&parameters,NULL,NULL);
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
    fprintf(stderr,"Usage: %s [--threads N] [--group-size G] [--generic] CELLS STEPS DT CELLS_CSV\n",argv[0]);
    return 2;
}
