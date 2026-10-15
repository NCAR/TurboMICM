/* SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_VERSION_H
#define CHEM_VERSION_H
/* Compatibility source version from the supplied micm/version.hpp. */
const char *chem_micm_version(void);
unsigned chem_micm_version_major(void);
unsigned chem_micm_version_minor(void);
unsigned chem_micm_version_patch(void);
unsigned chem_micm_version_tweak(void);
#endif
