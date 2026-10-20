#include "chemistry.h"

#define CHEM_STRINGIFY_INNER(value) #value
#define CHEM_STRINGIFY(value) CHEM_STRINGIFY_INNER(value)

const char *chem_version(void) {
  return CHEM_STRINGIFY(CHEM_VERSION_MAJOR) "." \
         CHEM_STRINGIFY(CHEM_VERSION_MINOR) "." \
         CHEM_STRINGIFY(CHEM_VERSION_PATCH) "+" \
         CHEM_STRINGIFY(CHEM_VERSION_TWEAK);
}

int chem_version_major(void) { return CHEM_VERSION_MAJOR; }
int chem_version_minor(void) { return CHEM_VERSION_MINOR; }
int chem_version_patch(void) { return CHEM_VERSION_PATCH; }
int chem_version_tweak(void) { return CHEM_VERSION_TWEAK; }
