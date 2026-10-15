#ifndef CHEM_BUILDER_H
#define CHEM_BUILDER_H
#include "chem.h"

/* Declarative external state registration precedes model configuration. The
   configure callback resolves names against the assembled descriptor and writes
   one process model plus exactly constraint_count constraint descriptors.
   Arrays referenced by those outputs and all callback contexts remain borrowed.
   finalize runs once after the sparse plan is available. */
typedef ChemStatus (*ChemModelConfigure)(const ChemMechanismData *,
    ChemExternalModel *process, ChemConstraint *constraints, void *context);
typedef struct {
    size_t variable_count, parameter_count, constraint_count;
    const char *const *variable_names, *const *parameter_names;
    const double *initial_variables, *initial_parameters, *absolute_tolerances;
    size_t used_species_count;
    const char *const *used_species_names;
    ChemModelConfigure configure;
    void *context;
} ChemModelDefinition;
typedef struct {
    ChemMechanismData base;
    size_t model_count;
    const ChemModelDefinition *models;
    int ignore_unused_species;
} ChemBuilder;
/* Reusable builder: every successful call creates independent merged descriptor
   tables owned by the resulting mechanism. Model equations and reactions remain
   caller-owned. The base species are the required gas system. */
ChemStatus chem_builder_build(const ChemBuilder *, ChemMechanism **);
size_t chem_data_species_index(const ChemMechanismData *, const char *name);
size_t chem_data_parameter_index(const ChemMechanismData *, const char *name);
#endif
