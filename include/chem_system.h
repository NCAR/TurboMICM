/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_SYSTEM_H
#define CHEM_SYSTEM_H
#include "chem.h"
#include "chem_external.h"

typedef enum { CHEM_PROPERTY_REAL, CHEM_PROPERTY_STRING, CHEM_PROPERTY_INT, CHEM_PROPERTY_BOOL } ChemPropertyType;
typedef union { double real;const char *string;int integer;int boolean; } ChemPropertyValue;
typedef struct { const char *name;ChemPropertyType type;ChemPropertyValue value; } ChemProperty;
typedef struct { double constant,temperature,pressure,air_density;int enabled; } ChemParameterization;
typedef struct ChemSpeciesInfo ChemSpeciesInfo;
/* Owning constructors copy all names/properties. Assignment replaces an owned
 * destination transactionally; property getters return borrowed string values. */
ChemStatus chem_species_info_create(const char *name,size_t property_count,const ChemProperty *,ChemSpeciesInfo **);
ChemStatus chem_species_info_copy(const ChemSpeciesInfo *,ChemSpeciesInfo **);
ChemStatus chem_species_info_assign(ChemSpeciesInfo **,const ChemSpeciesInfo *);
void chem_species_info_free(ChemSpeciesInfo *);
const char *chem_species_info_name(const ChemSpeciesInfo *);
size_t chem_species_property_count(const ChemSpeciesInfo *,ChemPropertyType);
int chem_species_has_property(const ChemSpeciesInfo *,const char *);
ChemStatus chem_species_set_property(ChemSpeciesInfo *,ChemProperty);
ChemStatus chem_species_get_property(const ChemSpeciesInfo *,const char *,ChemPropertyType,ChemPropertyValue *);
void chem_species_parameterize(ChemSpeciesInfo *,ChemParameterization);
void chem_species_set_third_body(ChemSpeciesInfo *);
int chem_species_is_parameterized(const ChemSpeciesInfo *);
double chem_species_parameter_value(const ChemSpeciesInfo *,const ChemConditions *);

typedef struct {
    ChemSpeciesInfo *species;
    double diffusion,density;
    int has_diffusion,has_density;
} ChemPhaseSpecies;
ChemStatus chem_phase_species_create(const ChemSpeciesInfo *,ChemPhaseSpecies *);
void chem_phase_species_free(ChemPhaseSpecies *);
void chem_phase_species_diffusion(ChemPhaseSpecies *,double);
void chem_phase_species_density(ChemPhaseSpecies *,double);
typedef struct ChemPhase ChemPhase;
ChemStatus chem_phase_create(const char *name,size_t count,const ChemPhaseSpecies *,ChemPhase **);
ChemStatus chem_phase_copy(const ChemPhase *,ChemPhase **);
void chem_phase_free(ChemPhase *);
const char *chem_phase_name(const ChemPhase *);
size_t chem_phase_count(const ChemPhase *);
size_t chem_phase_state_size(const ChemPhase *);
const ChemPhaseSpecies *chem_phase_species(const ChemPhase *,size_t);
typedef struct { size_t count;char **values; } ChemNames;
void chem_names_free(ChemNames *);
ChemStatus chem_phase_names(const ChemPhase *,int include_phase,ChemNames *);
typedef struct ChemSystem ChemSystem;
ChemStatus chem_system_create(const ChemPhase *gas_phase,ChemSystem **);
void chem_system_free(ChemSystem *);
const ChemPhase *chem_system_gas_phase(const ChemSystem *);
size_t chem_system_state_size(const ChemSystem *);
typedef const char *(*ChemNameOrder)(const char *const *names,size_t count,size_t i,void *context);
ChemStatus chem_system_names(const ChemSystem *,ChemNameOrder,void *context,ChemNames *);

typedef struct { const char *name;double yield; } ChemNamedProduct;
typedef struct ChemReactionInfo ChemReactionInfo;
ChemStatus chem_reaction_info_create(const ChemPhase *,size_t reactant_count,const char *const *reactants,
    size_t product_count,const ChemNamedProduct *products,ChemRate,ChemReactionInfo **);
ChemStatus chem_reaction_info_copy(const ChemReactionInfo *,ChemReactionInfo **);
ChemStatus chem_reaction_info_assign(ChemReactionInfo **,const ChemReactionInfo *);
void chem_reaction_info_free(ChemReactionInfo *);
const ChemPhase *chem_reaction_phase(const ChemReactionInfo *);
const ChemReaction *chem_reaction_definition(const ChemReactionInfo *);
const char *chem_reaction_reactant(const ChemReactionInfo *,size_t);
const char *chem_reaction_product(const ChemReactionInfo *,size_t);
/* The resulting mechanism owns a system snapshot, including parameterization
 * coefficients. Input metadata and reaction definitions may be freed afterward. */
ChemStatus chem_system_mechanism(const ChemSystem *,size_t reaction_count,
    const ChemReactionInfo *const *reactions,ChemMechanism **);
ChemStatus chem_system_mechanism_with_models(const ChemSystem *,size_t,const ChemReactionInfo *const *,
    size_t,const ChemExternalModel *,ChemMechanism **);
const ChemSystem *chem_mechanism_system(const ChemMechanism *);
#endif
