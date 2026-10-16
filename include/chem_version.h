#ifndef CHEM_VERSION_H
#define CHEM_VERSION_H
/* Compatibility version reported by the supplied source's version.hpp. */
#define CHEM_VERSION_MAJOR 3u
#define CHEM_VERSION_MINOR 14u
#define CHEM_VERSION_PATCH 0u
#define CHEM_VERSION_TWEAK 0u
const char *chem_version(void);
unsigned chem_version_major(void);
unsigned chem_version_minor(void);
unsigned chem_version_patch(void);
unsigned chem_version_tweak(void);
#endif
