#ifndef CHEM_VERSION_H
#define CHEM_VERSION_H
/* Version of the reference API and mechanism semantics targeted by this port. */
const char *chem_compatibility_version(void);
unsigned chem_compatibility_version_major(void);
unsigned chem_compatibility_version_minor(void);
unsigned chem_compatibility_version_patch(void);
unsigned chem_compatibility_version_tweak(void);
#endif
