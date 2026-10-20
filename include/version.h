#ifndef CHEM_VERSION_H
#define CHEM_VERSION_H

#define CHEM_VERSION_MAJOR 3u
#define CHEM_VERSION_MINOR 14u
#define CHEM_VERSION_PATCH 0u
#define CHEM_VERSION_TWEAK 0u
#define CHEM_VERSION_STRING "3.14.0"

static inline const char *chem_get_version(void)
{
  return CHEM_VERSION_STRING;
}

static inline unsigned chem_get_version_major(void)
{
  return CHEM_VERSION_MAJOR;
}

static inline unsigned chem_get_version_minor(void)
{
  return CHEM_VERSION_MINOR;
}

static inline unsigned chem_get_version_patch(void)
{
  return CHEM_VERSION_PATCH;
}

static inline unsigned chem_get_version_tweak(void)
{
  return CHEM_VERSION_TWEAK;
}

#endif
