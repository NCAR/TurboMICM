#include "state.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static ptrdiff_t name_index(const char *const *names,size_t count,const char *name)
{
    if(!name || (count && !names)) return -1;
    for(size_t i=0;i<count;++i) if(names[i] && !strcmp(names[i],name)) return (ptrdiff_t)i;
    return -1;
}
int chem_state_init(ChemState *s,const ChemStateParameters *p,size_t cells,size_t width)
{
    if(!s || !p || !cells || !width || (p->variable_count && !p->variable_names) || (p->parameter_count && !p->parameter_names) ||
       cells>SIZE_MAX/sizeof(ChemConditions) || p->variable_count>SIZE_MAX/sizeof(double)) return 1;
    ChemState t={0}; t.cells=cells; t.parameters=*p;
    t.conditions=calloc(cells,sizeof(*t.conditions));
    t.absolute_tolerances=p->variable_count?calloc(p->variable_count,sizeof(double)):NULL;
    if(!t.conditions || (p->variable_count && !t.absolute_tolerances) ||
       chem_matrix_init(&t.variables,cells,p->variable_count,width,width!=1,0) ||
       chem_matrix_init(&t.custom_parameters,cells,p->parameter_count,width,width!=1,0) ||
       chem_matrix_init(&t.rate_constants,cells,p->rate_count,width,width!=1,0)) { chem_state_free(&t); return 1; }
    if(p->variable_count && p->absolute_tolerances) memcpy(t.absolute_tolerances,p->absolute_tolerances,p->variable_count*sizeof(double));
    else for(size_t i=0;i<p->variable_count;++i) t.absolute_tolerances[i]=1e-3;
    t.parameters.absolute_tolerances=t.absolute_tolerances;
    *s=t; return 0;
}
void chem_state_free(ChemState *s)
{
    if(!s) return;
    chem_matrix_free(&s->variables); chem_matrix_free(&s->custom_parameters); chem_matrix_free(&s->rate_constants);
    free(s->conditions); free(s->absolute_tolerances); *s=(ChemState){0};
}
int chem_state_copy(ChemState *out,const ChemState *in)
{
    if(!out || !in) return 1;
    if(out==in) return 0;
    if(!in->cells) { chem_state_free(out); return 0; }
    ChemState copy={0}; if(chem_state_init(&copy,&in->parameters,in->cells,in->variables.width)) return 1;
    memcpy(copy.conditions,in->conditions,in->cells*sizeof(*in->conditions));
    chem_matrix_copy(&copy.variables,&in->variables); chem_matrix_copy(&copy.custom_parameters,&in->custom_parameters);
    chem_matrix_copy(&copy.rate_constants,&in->rate_constants);
    chem_state_free(out); *out=copy; return 0;
}
void chem_state_move(ChemState *out,ChemState *in)
{
    if(!out || !in || out==in) return;
    chem_state_free(out); *out=*in; *in=(ChemState){0};
}
ptrdiff_t chem_state_variable_index(const ChemState *s,const char *name)
{
    return s?name_index(s->parameters.variable_names,s->parameters.variable_count,name):-1;
}
ptrdiff_t chem_state_parameter_index(const ChemState *s,const char *name)
{
    return s?name_index(s->parameters.parameter_names,s->parameters.parameter_count,name):-1;
}
double *chem_state_variable(ChemState *s,size_t cell,const char *name)
{
    ptrdiff_t i=chem_state_variable_index(s,name); return i<0?NULL:chem_matrix_at(&s->variables,cell,(size_t)i);
}
const double *chem_state_const_variable(const ChemState *s,size_t cell,const char *name)
{
    ptrdiff_t i=chem_state_variable_index(s,name); return i<0?NULL:chem_matrix_const_at(&s->variables,cell,(size_t)i);
}
static int set_named(ChemState *s,const char *name,const double *values,size_t count,int parameter)
{
    if(!s || count!=s->cells || (count && !values)) return 1;
    ptrdiff_t col=parameter?chem_state_parameter_index(s,name):chem_state_variable_index(s,name);
    if(col<0) return 1;
    return chem_matrix_vector_to_column(parameter?&s->custom_parameters:&s->variables,(size_t)col,values,count);
}
int chem_state_set_concentration(ChemState *s,const char *n,const double *v,size_t c) { return set_named(s,n,v,c,0); }
int chem_state_set_parameter(ChemState *s,const char *n,const double *v,size_t c) { return set_named(s,n,v,c,1); }
static int set_many(ChemState *s,const ChemNamedValues *items,size_t count,int parameter)
{
    if(!s || (count && !items)) return 1;
    for(size_t i=0;i<count;++i) if(items[i].count!=s->cells || !items[i].values ||
        (parameter?chem_state_parameter_index(s,items[i].name):chem_state_variable_index(s,items[i].name))<0) return 1;
    for(size_t i=0;i<count;++i) if(set_named(s,items[i].name,items[i].values,items[i].count,parameter)) return 1;
    return 0;
}
int chem_state_set_concentrations(ChemState *s,const ChemNamedValues *v,size_t n) { return set_many(s,v,n,0); }
int chem_state_set_parameters(ChemState *s,const ChemNamedValues *v,size_t n) { return set_many(s,v,n,1); }
int chem_state_set_parameter_rows(ChemState *s,const double *const *rows,const size_t *sizes,size_t count)
{
    if(!s || count!=s->cells || (count && (!rows || !sizes))) return 1;
    for(size_t i=0;i<count;++i) if(sizes[i]!=s->parameters.parameter_count || (sizes[i] && !rows[i])) return 1;
    for(size_t i=0;i<count;++i) if(chem_matrix_set_row(&s->custom_parameters,i,rows[i],sizes[i])) return 1;
    return 0;
}
int chem_state_update_rate_constants(ChemState *s,const ChemMechanism *m,const ChemCallback *callbacks)
{
    if(!s || !m || s->parameters.rate_count!=m->reaction_count || (m->parameter_count && !m->parameter_names) ||
       m->parameter_count>SIZE_MAX-m->reaction_count || m->parameter_count+m->reaction_count>SIZE_MAX/sizeof(double)) return 1;
    size_t n=m->parameter_count+m->reaction_count;
    double *data=malloc((n?n:1)*sizeof(double)); if(!data) return 1;
    double *params=data, *rates=data+m->parameter_count; int status=0;
    for(size_t cell=0;cell<s->cells;++cell) {
        for(size_t i=0;i<m->parameter_count;++i) {
            ptrdiff_t col=chem_state_parameter_index(s,m->parameter_names[i]);
            if(col<0) { status=1; break; }
            params[i]=*chem_matrix_at(&s->custom_parameters,cell,(size_t)col);
        }
        if(status || chem_rate_constants_with_callbacks(m,s->conditions[cell],params,callbacks,rates) ||
           chem_matrix_set_row(&s->rate_constants,cell,rates,m->reaction_count)) { status=1; break; }
    }
    free(data); return status;
}
int chem_state_apply_phase_tolerances(ChemState *s,ChemPhase phase)
{
    if(!s || (phase.count && !phase.species)) return 1;
    for(size_t i=0;i<phase.count;++i) {
        const ChemSpecies *species=phase.species[i].species;
        if(!species) return 1;
        if(species->is_parameterized) continue;
        ptrdiff_t index=chem_state_variable_index(s,species->name);
        if(index<0) return 1;
        ChemPropertyValue tolerance;
        if(chem_species_get(species,"absolute tolerance",CHEM_PROPERTY_REAL,&tolerance)==CHEM_PROPERTY_OK) {
            if(!(tolerance.real>0)) return 1;
            s->absolute_tolerances[index]=tolerance.real;
        }
    }
    return 0;
}
int chem_state_update_constraint_parameters(ChemState *s,ChemConstraintSet set)
{
    if(!s || chem_constraint_set_validate(set,s->parameters.variable_count)) return 1;
    for(size_t i=0;i<set.count;++i) if(set.items[i].type==CHEM_EQUILIBRIUM_CONSTRAINT) {
        ptrdiff_t column=chem_state_parameter_index(s,set.items[i].name);
        if(column<0) return 1;
        for(size_t cell=0;cell<s->cells;++cell) {
            double value;
            if(chem_constraint_equilibrium_constant(&set.items[i],s->conditions[cell].temperature,&value)) return 1;
            *chem_matrix_at(&s->custom_parameters,cell,(size_t)column)=value;
        }
    }
    return 0;
}
