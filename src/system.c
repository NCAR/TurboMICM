#include "chem_system.h"
#include <stdlib.h>
#include <string.h>

static char *copy_string(const char *s)
{
    size_t n=strlen(s)+1; char *copy=malloc(n);
    if (copy) memcpy(copy,s,n);
    return copy;
}
static int valid_type(chem_property_type type)
{ return type>=CHEM_PROPERTY_DOUBLE && type<=CHEM_PROPERTY_BOOL; }
void chem_species_info_destroy(chem_species_info *s)
{
    if (!s) return;
    free((void *)s->species.name);
    for (size_t i=0;i<s->n_properties;++i) {
        free((void *)s->properties[i].key);
        if (s->properties[i].type==CHEM_PROPERTY_STRING) free((void *)s->properties[i].value.string);
    }
    free(s->properties); *s=(chem_species_info){0};
}
chem_status chem_species_property_set(chem_species_info *s, const char *key, chem_property_type type, chem_property_value value)
{
    if (!s || !key || !valid_type(type) || (type==CHEM_PROPERTY_STRING && !value.string)) return CHEM_INVALID;
    if (type==CHEM_PROPERTY_STRING) {
        value.string=copy_string(value.string);
        if (!value.string) return CHEM_MEMORY;
    }
    size_t i=0;
    for (;i<s->n_properties;++i) if (s->properties[i].type==type && !strcmp(s->properties[i].key,key)) break;
    if (i==s->n_properties) {
        char *name=copy_string(key);
        chem_property *properties=name ? realloc(s->properties,(i+1)*sizeof(*properties)) : NULL;
        if (!properties) {
            free(name); if (type==CHEM_PROPERTY_STRING) free((void *)value.string);
            return CHEM_MEMORY;
        }
        s->properties=properties; s->properties[i]=(chem_property){name,type,{0}}; ++s->n_properties;
    } else if (type==CHEM_PROPERTY_STRING) free((void *)s->properties[i].value.string);
    s->properties[i].value=value;
    if (type==CHEM_PROPERTY_DOUBLE) {
        if (!strcmp(key,"absolute tolerance")) s->species.absolute_tolerance=value.real;
        else if (!strcmp(key,"molecular weight [kg mol-1]")) s->species.molecular_weight=value.real;
        else if (!strcmp(key,"diffusion coefficient [m2 s-1]")) s->species.diffusion=value.real;
    }
    return CHEM_OK;
}
chem_status chem_species_property_get(const chem_species_info *s, const char *key, chem_property_type type, chem_property_value *value)
{
    if (!s || !key || !value || !valid_type(type)) return CHEM_INVALID;
    for (size_t i=0;i<s->n_properties;++i) if (s->properties[i].type==type && !strcmp(s->properties[i].key,key)) {
        *value=s->properties[i].value; return CHEM_OK;
    }
    return CHEM_NOT_FOUND;
}
size_t chem_species_property_count(const chem_species_info *s, chem_property_type type)
{
    size_t count=0; for (size_t i=0;i<s->n_properties;++i) count+=s->properties[i].type==type;
    return count;
}
int chem_species_has_property(const chem_species_info *s, const char *key)
{
    if (!s || !key) return 0;
    for (size_t i=0;i<s->n_properties;++i) if (!strcmp(s->properties[i].key,key)) return 1;
    return 0;
}
chem_status chem_species_info_create(chem_species_info *s, const char *name, size_t count, const chem_property *properties)
{
    if (!s || !name || (count && !properties)) return CHEM_INVALID;
    *s=(chem_species_info){0}; s->species.name=copy_string(name);
    if (!s->species.name) return CHEM_MEMORY;
    for (size_t i=0;i<count;++i) {
        chem_status status=chem_species_property_set(s,properties[i].key,properties[i].type,properties[i].value);
        if (status!=CHEM_OK) { chem_species_info_destroy(s); return status; }
    }
    return CHEM_OK;
}
chem_status chem_species_info_clone(chem_species_info *dest, const chem_species_info *src)
{
    if (!dest || !src || dest==src) return CHEM_INVALID;
    chem_status status=chem_species_info_create(dest,src->species.name,src->n_properties,src->properties);
    if (status!=CHEM_OK) return status;
    const char *name=dest->species.name; dest->species=src->species; dest->species.name=name;
    return CHEM_OK;
}
chem_status chem_species_info_replace(chem_species_info *dest, const chem_species_info *src)
{
    if (!dest || !src) return CHEM_INVALID;
    if (dest==src) return CHEM_OK;
    chem_species_info copy={0}; chem_status status=chem_species_info_clone(&copy,src);
    if (status==CHEM_OK) { chem_species_info_destroy(dest); *dest=copy; }
    return status;
}
void chem_species_parameterize(chem_species_info *s, chem_parameterization p)
{
    s->species.type=CHEM_SPECIES_PARAMETER; s->species.parameter=NULL; s->species.context=NULL;
    s->species.affine=p;
}
void chem_species_third_body(chem_species_info *s)
{
    s->species.type=CHEM_SPECIES_THIRD_BODY; s->species.parameter=NULL; s->species.context=NULL;
    s->species.affine=(chem_parameterization){0,0,0,1};
}
double chem_species_evaluate(const chem_species *s, const chem_conditions *c)
{
    if (s->type==CHEM_SPECIES_THIRD_BODY) return c->air_density;
    if (s->parameter) return s->parameter(c,NULL,s->context);
    return s->affine.constant+s->affine.temperature*c->temperature+s->affine.pressure*c->pressure+s->affine.air_density*c->air_density;
}
chem_status chem_phase_species_create(chem_phase_species *dest, const chem_species_info *src, const double *diffusion, const double *density)
{
    if (!dest || !src) return CHEM_INVALID;
    *dest=(chem_phase_species){0}; chem_status status=chem_species_info_clone(&dest->info,src);
    if (status!=CHEM_OK) return status;
    if (diffusion) chem_phase_species_set_diffusion(dest,*diffusion);
    if (density) chem_phase_species_set_density(dest,*density);
    return CHEM_OK;
}
void chem_phase_species_destroy(chem_phase_species *s)
{ if (s) { chem_species_info_destroy(&s->info); *s=(chem_phase_species){0}; } }
void chem_phase_species_set_diffusion(chem_phase_species *s,double value) { s->has_diffusion=1; s->diffusion=value; }
void chem_phase_species_set_density(chem_phase_species *s,double value) { s->has_density=1; s->density=value; }
void chem_phase_destroy(chem_phase *p)
{
    if (!p) return;
    for (size_t i=0;i<p->n_species;++i) chem_phase_species_destroy(&p->species[i]);
    free(p->species); free(p->name); *p=(chem_phase){0};
}
chem_status chem_phase_create(chem_phase *p,const char *name,size_t n,const chem_phase_species *species)
{
    if (!p || !name || (n && !species)) return CHEM_INVALID;
    *p=(chem_phase){0}; p->name=copy_string(name); p->species=calloc(n ? n : 1,sizeof(*p->species));
    if (!p->name || !p->species) { chem_phase_destroy(p); return CHEM_MEMORY; }
    for (size_t i=0;i<n;++i) {
        chem_status status=chem_phase_species_create(p->species+i,&species[i].info,
            species[i].has_diffusion ? &species[i].diffusion : NULL,species[i].has_density ? &species[i].density : NULL);
        if (status!=CHEM_OK) { chem_phase_destroy(p); return status; }
        ++p->n_species;
    }
    return CHEM_OK;
}
size_t chem_phase_state_size(const chem_phase *p)
{ size_t n=0; for (size_t i=0;i<p->n_species;++i) n+=p->species[i].info.species.type==CHEM_SPECIES_VARIABLE; return n; }
void chem_names_destroy(chem_names *names)
{
    if (!names) return;
    for (size_t i=0;i<names->count;++i) free(names->names[i]);
    free(names->names); *names=(chem_names){0};
}
chem_status chem_phase_names(const chem_phase *p,int prefix,chem_names *names)
{
    if (!p || !names) return CHEM_INVALID;
    *names=(chem_names){0}; size_t count=chem_phase_state_size(p);
    names->names=calloc(count ? count : 1,sizeof(*names->names));
    if (!names->names) return CHEM_MEMORY;
    for (size_t i=0;i<p->n_species;++i) if (p->species[i].info.species.type==CHEM_SPECIES_VARIABLE) {
        const char *s=p->species[i].info.species.name; size_t n=prefix ? strlen(p->name)+1 : 0;
        char *name=malloc(n+strlen(s)+1);
        if (!name) { chem_names_destroy(names); return CHEM_MEMORY; }
        if (prefix) { memcpy(name,p->name,n-1); name[n-1]='.'; }
        strcpy(name+n,s); names->names[names->count++]=name;
    }
    return CHEM_OK;
}
chem_status chem_system_create(chem_system *s,const chem_phase *p)
{ return (!s || !p) ? CHEM_INVALID : chem_phase_create(&s->gas,p->name,p->n_species,p->species); }
void chem_system_destroy(chem_system *s) { if (s) chem_phase_destroy(&s->gas); }
size_t chem_system_state_size(const chem_system *s) { return chem_phase_state_size(&s->gas); }
chem_status chem_system_names(const chem_system *s,const size_t *order,chem_names *names)
{
    if (!s || !names) return CHEM_INVALID;
    chem_status status=chem_phase_names(&s->gas,0,names);
    if (status!=CHEM_OK || !order) return status;
    size_t n=names->count; char **sorted=calloc(n ? n : 1,sizeof(*sorted));
    if (!sorted) { chem_names_destroy(names); return CHEM_MEMORY; }
    for (size_t i=0;i<n;++i) {
        if (order[i]>=n) { free(sorted); chem_names_destroy(names); return CHEM_INVALID; }
        for (size_t j=0;j<i;++j) if (order[i]==order[j]) {
            free(sorted); chem_names_destroy(names); return CHEM_INVALID;
        }
        sorted[i]=names->names[order[i]];
    }
    free(names->names); names->names=sorted; return CHEM_OK;
}
