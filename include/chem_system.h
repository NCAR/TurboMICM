#ifndef CHEM_SYSTEM_H
#define CHEM_SYSTEM_H
#include "chem.h"

/* Owned metadata. All create/clone destinations must initially be zeroed;
 * replace is transactional and accepts an already initialized destination. */
typedef enum { CHEM_PROPERTY_DOUBLE, CHEM_PROPERTY_STRING, CHEM_PROPERTY_INT, CHEM_PROPERTY_BOOL } chem_property_type;
typedef union { double real; const char *string; int integer; int boolean; } chem_property_value;
typedef struct { const char *key; chem_property_type type; chem_property_value value; } chem_property;
typedef struct { chem_species species; size_t n_properties; chem_property *properties; } chem_species_info;
chem_status chem_species_info_create(chem_species_info *, const char *, size_t, const chem_property *);
chem_status chem_species_info_clone(chem_species_info *, const chem_species_info *);
chem_status chem_species_info_replace(chem_species_info *, const chem_species_info *);
void chem_species_info_destroy(chem_species_info *);
chem_status chem_species_property_set(chem_species_info *, const char *, chem_property_type, chem_property_value);
chem_status chem_species_property_get(const chem_species_info *, const char *, chem_property_type, chem_property_value *);
size_t chem_species_property_count(const chem_species_info *, chem_property_type);
int chem_species_has_property(const chem_species_info *, const char *);
void chem_species_parameterize(chem_species_info *, chem_parameterization);
void chem_species_third_body(chem_species_info *);

typedef struct {
    chem_species_info info;
    int has_diffusion, has_density;
    double diffusion, density;
} chem_phase_species;
chem_status chem_phase_species_create(chem_phase_species *, const chem_species_info *, const double *, const double *);
void chem_phase_species_destroy(chem_phase_species *);
void chem_phase_species_set_diffusion(chem_phase_species *, double);
void chem_phase_species_set_density(chem_phase_species *, double);
typedef struct { char *name; size_t n_species; chem_phase_species *species; } chem_phase;
chem_status chem_phase_create(chem_phase *, const char *, size_t, const chem_phase_species *);
void chem_phase_destroy(chem_phase *);
size_t chem_phase_state_size(const chem_phase *);
typedef struct { size_t count; char **names; } chem_names;
void chem_names_destroy(chem_names *);
chem_status chem_phase_names(const chem_phase *, int prefix, chem_names *);
typedef struct { chem_phase gas; } chem_system;
chem_status chem_system_create(chem_system *, const chem_phase *);
void chem_system_destroy(chem_system *);
size_t chem_system_state_size(const chem_system *);
/* Optional permutation maps each output position to an unpermuted position. */
chem_status chem_system_names(const chem_system *, const size_t *, chem_names *);
#endif
