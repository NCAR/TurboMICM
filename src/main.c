#include "chemistry.h"

#include <stdio.h>
#include <string.h>

extern const chem_mechanism chem_generated_mechanism;

int main(int argc, char **argv) {
  if (argc == 2 && strcmp(argv[1], "--inspect") == 0) {
    printf("species=%zu reactions=%zu\n", chem_generated_mechanism.species_count,
           chem_generated_mechanism.reaction_count);
    return 0;
  }
  fprintf(stderr, "TS1 integrator is not implemented yet. Usage: %s --inspect\n", argv[0]);
  return 2;
}
