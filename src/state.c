#include "state.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void chem_state_destroy(ChemState *s)
{
    if (!s) return;
    if (s->species_names) for (size_t i=0; i<s->species_count; ++i) free(s->species_names[i]);
    if (s->parameter_names) for (size_t i=0; i<s->parameter_count; ++i) free(s->parameter_names[i]);
    free(s->species_names); free(s->parameter_names); free(s->conditions);
    free(s->concentrations); free(s->parameters); free(s->rates); free(s->absolute_tolerances);
    *s = (ChemState){0};
}
static ChemStateStatus copy_names(char ***out, const char *const *names, size_t n)
{
    if (!n) return CHEM_STATE_OK;
    if (!names) return CHEM_STATE_INVALID;
    if (n > SIZE_MAX / sizeof(**out)) return CHEM_STATE_SIZE;
    *out = calloc(n, sizeof(**out));
    if (!*out) return CHEM_STATE_NO_MEMORY;
    for (size_t i=0; i<n; ++i) {
        if (!names[i]) return CHEM_STATE_INVALID;
        for (size_t j=0; j<i; ++j)
            if (!strcmp(names[i],names[j])) return CHEM_STATE_DUPLICATE_NAME;
        size_t len = strlen(names[i]);
        (*out)[i] = malloc(len + 1);
        if (!(*out)[i]) return CHEM_STATE_NO_MEMORY;
        memcpy((*out)[i],names[i],len + 1);
    }
    return CHEM_STATE_OK;
}
static int array_size(size_t a, size_t b, size_t element_size)
{
    return !a || (b <= SIZE_MAX / a && a*b <= SIZE_MAX / element_size);
}
static void *zero_array(size_t n, size_t size)
{
    return n ? calloc(n,size) : NULL;
}
ChemStateStatus chem_state_init(ChemState *state, const ChemStateConfig *c)
{
    if (!state || !c) return CHEM_STATE_INVALID;
    if (!array_size(c->cells,c->species_count,sizeof(double)) ||
        !array_size(c->cells,c->parameter_count,sizeof(double)) ||
        !array_size(c->cells,c->rate_count,sizeof(double)) ||
        !array_size(c->cells,1,sizeof(ChemConditions)) ||
        !array_size(c->species_count,1,sizeof(double))) return CHEM_STATE_SIZE;
    ChemState s = { .cells=c->cells, .species_count=c->species_count,
        .parameter_count=c->parameter_count, .rate_count=c->rate_count,
        .relative_tolerance=c->relative_tolerance };
    ChemStateStatus result = copy_names(&s.species_names,c->species_names,c->species_count);
    if (result != CHEM_STATE_OK) goto fail;
    result = copy_names(&s.parameter_names,c->parameter_names,c->parameter_count);
    if (result != CHEM_STATE_OK) goto fail;
    s.conditions=zero_array(c->cells,sizeof(*s.conditions));
    s.concentrations=zero_array(c->cells*c->species_count,sizeof(double));
    s.parameters=zero_array(c->cells*c->parameter_count,sizeof(double));
    s.rates=zero_array(c->cells*c->rate_count,sizeof(double));
    s.absolute_tolerances=zero_array(c->species_count,sizeof(double));
    if ((c->cells && !s.conditions) || (c->cells && c->species_count && !s.concentrations) ||
        (c->cells && c->parameter_count && !s.parameters) || (c->cells && c->rate_count && !s.rates) ||
        (c->species_count && !s.absolute_tolerances)) { result=CHEM_STATE_NO_MEMORY; goto fail; }
    if (c->absolute_tolerances && c->species_count)
        memcpy(s.absolute_tolerances,c->absolute_tolerances,c->species_count*sizeof(double));
    chem_state_destroy(state);
    *state=s;
    return CHEM_STATE_OK;
fail:
    chem_state_destroy(&s);
    return result;
}
ChemStateStatus chem_state_clone(ChemState *destination, const ChemState *source)
{
    if (!destination || !source) return CHEM_STATE_INVALID;
    if (destination == source) return CHEM_STATE_OK;
    ChemState copy = {0};
    ChemStateConfig c = {source->cells,source->species_count,source->parameter_count,source->rate_count,
        (const char *const *)source->species_names,(const char *const *)source->parameter_names,
        source->relative_tolerance,source->absolute_tolerances};
    ChemStateStatus status=chem_state_init(&copy,&c);
    if (status != CHEM_STATE_OK) return status;
    if (source->cells) memcpy(copy.conditions,source->conditions,source->cells*sizeof(*copy.conditions));
    if (source->cells && source->species_count) memcpy(copy.concentrations,source->concentrations,source->cells*source->species_count*sizeof(double));
    if (source->cells && source->parameter_count) memcpy(copy.parameters,source->parameters,source->cells*source->parameter_count*sizeof(double));
    if (source->cells && source->rate_count) memcpy(copy.rates,source->rates,source->cells*source->rate_count*sizeof(double));
    chem_state_destroy(destination); *destination=copy;
    return CHEM_STATE_OK;
}
void chem_state_move(ChemState *destination, ChemState *source)
{
    if (!destination || !source || destination == source) return;
    chem_state_destroy(destination); *destination=*source; *source=(ChemState){0};
}
static size_t find_name(char *const *names, size_t n, const char *name)
{
    if (name) for (size_t i=0;i<n;++i) if (!strcmp(names[i],name)) return i;
    return n;
}
static ChemStateStatus set_column(size_t cells, size_t columns, char *const *names, double *matrix,
                                  const char *name, const double *values, size_t count)
{
    size_t col=find_name(names,columns,name);
    if (col==columns) return CHEM_STATE_UNKNOWN_NAME;
    if (count!=cells) return CHEM_STATE_SIZE;
    if (count && !values) return CHEM_STATE_INVALID;
    for (size_t cell=0;cell<cells;++cell) matrix[cell*columns+col]=values[cell];
    return CHEM_STATE_OK;
}
ChemStateStatus chem_state_set_concentration(ChemState *s, const char *name, const double *v, size_t count)
{
    return s ? set_column(s->cells,s->species_count,s->species_names,s->concentrations,name,v,count) : CHEM_STATE_INVALID;
}
ChemStateStatus chem_state_set_parameter(ChemState *s, const char *name, const double *v, size_t count)
{
    return s ? set_column(s->cells,s->parameter_count,s->parameter_names,s->parameters,name,v,count) : CHEM_STATE_INVALID;
}
static ChemStateStatus set_columns(size_t cells, size_t n, char *const *names, double *matrix,
                                   const ChemStateColumn *columns, size_t count)
{
    if (count && !columns) return CHEM_STATE_INVALID;
    for (size_t i=0;i<count;++i) {
        if (find_name(names,n,columns[i].name)==n) return CHEM_STATE_UNKNOWN_NAME;
        if (columns[i].count!=cells) return CHEM_STATE_SIZE;
        if (cells && !columns[i].values) return CHEM_STATE_INVALID;
    }
    for (size_t i=0;i<count;++i)
        set_column(cells,n,names,matrix,columns[i].name,columns[i].values,columns[i].count);
    return CHEM_STATE_OK;
}
ChemStateStatus chem_state_set_concentrations(ChemState *s, const ChemStateColumn *columns, size_t count)
{
    return s ? set_columns(s->cells,s->species_count,s->species_names,s->concentrations,columns,count) : CHEM_STATE_INVALID;
}
ChemStateStatus chem_state_set_parameters(ChemState *s, const ChemStateColumn *columns, size_t count)
{
    return s ? set_columns(s->cells,s->parameter_count,s->parameter_names,s->parameters,columns,count) : CHEM_STATE_INVALID;
}
double *chem_state_concentration(ChemState *s, size_t cell, const char *name)
{
    if (!s || cell>=s->cells) return NULL;
    size_t col=find_name(s->species_names,s->species_count,name);
    return col==s->species_count ? NULL : &s->concentrations[cell*s->species_count+col];
}
ChemStateStatus chem_state_set_parameter_table(ChemState *s, const double *values, size_t cells, size_t parameters)
{
    if (!s) return CHEM_STATE_INVALID;
    if (cells!=s->cells || parameters!=s->parameter_count) return CHEM_STATE_SIZE;
    if (cells && parameters && !values) return CHEM_STATE_INVALID;
    if (cells && parameters) memmove(s->parameters,values,cells*parameters*sizeof(double));
    return CHEM_STATE_OK;
}
