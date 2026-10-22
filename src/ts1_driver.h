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
#ifdef CHEM_SPECIALIZED
ChemStatus ts1_specialized_create(ChemMechanism **);
void ts1_specialized_attach(ChemSolver *);
void ts1_update_rates(ChemState *);
ChemStatus ts1_solve(ChemSolver *,ChemState *,double,const ChemRosenbrockParameters *,ChemSolverStats *);
#define chem_update_rates ts1_update_rates
#define chem_solve ts1_solve
#endif
static ChemStatus create_mechanism(ChemMechanism **out)
{
#ifdef CHEM_SPECIALIZED
    return ts1_specialized_create(out);
#else
    ChemStatus status=ts1_create(out);
    return status?status:chem_mechanism_reorder(*out);
#endif
}
static ChemStatus create_solver(const ChemMechanism *m,size_t cells,ChemSolver **out)
{
    ChemStatus status=chem_solver_create(m,cells,out);
#ifdef CHEM_SPECIALIZED
    if (!status) ts1_specialized_attach(*out);
#endif
    return status;
}
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
