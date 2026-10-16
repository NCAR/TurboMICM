#ifndef CHEM_PROCESS_H
#define CHEM_PROCESS_H
#include "chem_system.h"

typedef struct { const char *name; double yield; } chem_named_product;
/* A process owns its phase, rate label, and reaction lists. Callback contexts
 * remain borrowed. Species indices refer to this process's owned phase. */
typedef struct { chem_phase phase; chem_reaction reaction; } chem_process;
chem_status chem_process_create(chem_process *,const chem_phase *,size_t,const char *const *,size_t,const chem_named_product *,chem_rate);
chem_status chem_process_clone(chem_process *,const chem_process *);
chem_status chem_process_replace(chem_process *,const chem_process *);
void chem_process_destroy(chem_process *);
chem_status chem_surface_rate(chem_rate *,const char *,const chem_phase_species *,double);
int chem_rate_type_order(chem_rate_type);
void chem_reactions_sort(chem_reaction *,size_t);
size_t chem_model_rate_count(const chem_model *,chem_rate_type);
/* Cumulative offsets have the reference meaning for type-sorted reactions. */
size_t chem_model_rate_offset(const chem_model *,chem_rate_type);
size_t chem_model_parameterized_count(const chem_model *);
size_t chem_model_parameterized_index(const chem_model *,size_t);
double chem_model_parameterized_multiplier(const chem_model *,size_t,const chem_conditions *);
int chem_model_species_used(const chem_model *,size_t species);
/* Kinetic sparsity, before algebraic masking, fill, and added solver diagonals. */
chem_status chem_model_kinetic_pattern(const chem_model *,unsigned char *,size_t size);

/* Owns the entire system and reaction description used by its compiled model.
 * Process-local species indices are resolved by name against the supplied gas
 * system; reaction types are stably grouped as in the reference builder. */
typedef struct {
    chem_system system;
    chem_species *species;
    chem_reaction *reactions;
    chem_mechanism mechanism;
    chem_model model;
} chem_process_model;
chem_status chem_process_model_create(chem_process_model *,const chem_system *,size_t,const chem_process *,
    int preserve_species_order,int reject_unused,chem_lu_algorithm,chem_error *);
void chem_process_model_destroy(chem_process_model *);
#endif
