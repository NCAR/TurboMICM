#include "mechanism_metadata.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void free_names(char **names,size_t count)
{
    if (names) for (size_t i=0;i<count;++i) free(names[i]);
    free(names);
}
void chem_mechanism_metadata_destroy(ChemMechanismMetadata *m)
{
    if (!m) return;
    free_names(m->reaction_names,m->reaction_count);
    free_names(m->photolysis_names,m->photolysis_count);
    free_names(m->species_names,m->species_count);
    *m=(ChemMechanismMetadata){0};
}
static char *copy_name(const char *s)
{
    size_t n=strlen(s);
    if (n==SIZE_MAX) return NULL;
    char *copy=malloc(n+1);
    if (copy) memcpy(copy,s,n+1);
    return copy;
}
ChemMetadataStatus chem_mechanism_metadata_init(ChemMechanismMetadata *m,const ChemMechanism *mechanism,
    const char *const *reactions,const char *const *species,const size_t *photolysis,size_t count)
{
    if (!m || !mechanism || (mechanism->reaction_count && !reactions) ||
        (mechanism->species_count && !species) || (count && !photolysis)) return CHEM_METADATA_INVALID;
    size_t nr=mechanism->reaction_count,ns=mechanism->species_count;
    if (nr>SIZE_MAX/sizeof(char *) || ns>SIZE_MAX/sizeof(char *) || count>SIZE_MAX/sizeof(char *))
        return CHEM_METADATA_NO_MEMORY;
    for (size_t i=0;i<nr;++i) if (!reactions[i]) return CHEM_METADATA_INVALID;
    for (size_t i=0;i<ns;++i) if (!species[i]) return CHEM_METADATA_INVALID;
    for (size_t i=0;i<count;++i) if (photolysis[i]>=nr) return CHEM_METADATA_INVALID;
    ChemMechanismMetadata copy={.reaction_count=nr,.species_count=ns,.photolysis_count=count};
    copy.reaction_names=nr ? calloc(nr,sizeof(char *)) : NULL;
    copy.species_names=ns ? calloc(ns,sizeof(char *)) : NULL;
    copy.photolysis_names=count ? calloc(count,sizeof(char *)) : NULL;
    if ((nr && !copy.reaction_names) || (ns && !copy.species_names) || (count && !copy.photolysis_names)) goto fail;
    for (size_t i=0;i<nr;++i) if (!(copy.reaction_names[i]=copy_name(reactions[i]))) goto fail;
    for (size_t i=0;i<ns;++i) if (!(copy.species_names[i]=copy_name(species[i]))) goto fail;
    for (size_t i=0;i<count;++i) if (!(copy.photolysis_names[i]=copy_name(reactions[photolysis[i]]))) goto fail;
    chem_mechanism_metadata_destroy(m); *m=copy; return CHEM_METADATA_OK;
fail:
    chem_mechanism_metadata_destroy(&copy); return CHEM_METADATA_NO_MEMORY;
}
