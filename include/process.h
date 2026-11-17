#ifndef CHEM_PROCESS_H
#define CHEM_PROCESS_H
#include "chem.h"
#include "species.h"

typedef struct { ChemSpecies species; double yield; } ChemStoichSpecies;
typedef struct {
    size_t reactant_count, product_count;
    ChemSpecies *reactants;
    ChemStoichSpecies *products;
    ChemPhase phase;
    ChemRate rate;
} ChemProcess;
typedef struct {
    size_t reactant_count, product_count;
    const ChemSpecies *const *reactants, *const *products;
    const double *yields;
    const ChemPhase *phase;
    ChemRate rate;
} ChemProcessConfig;
/* Initialize to zero; the process owns independent copies of its inputs. */
ChemSpeciesStatus chem_process_init(ChemProcess *process, const ChemProcessConfig *config);
ChemSpeciesStatus chem_process_clone(ChemProcess *destination, const ChemProcess *source);
void chem_process_destroy(ChemProcess *process);

typedef enum { CHEM_BUILD_OK, CHEM_BUILD_INVALID, CHEM_BUILD_NO_MEMORY,
               CHEM_BUILD_MISSING_SYSTEM, CHEM_BUILD_UNKNOWN_SPECIES,
               CHEM_BUILD_UNUSED_SPECIES, CHEM_BUILD_BAD_TOLERANCE } ChemBuildStatus;
typedef struct {
    ChemSystem system;
    ChemMechanism mechanism;
    ChemReaction *reactions;
    size_t *reactants;
    ChemProduct *products;
    ChemRate *rates;
    ChemRateMultiplier *multipliers;
    ChemParameterizedSpecies *functions;
    const char **state_names; /* borrowed from system; owned pointer array */
    double *absolute_tolerances;
} ChemBuiltMechanism;
/* Compile reactions into indexed mass-action tables. Parameterized reactants
   become condition-dependent rate multipliers; parameterized products are
   excluded from state. Returns an error without modifying destination. */
ChemBuildStatus chem_build_mechanism(ChemBuiltMechanism *destination, const ChemSystem *system,
                                     size_t count, const ChemProcess *processes, int ignore_unused);
void chem_built_mechanism_destroy(ChemBuiltMechanism *mechanism);
#endif
