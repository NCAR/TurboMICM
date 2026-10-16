#include "test.h"
#include "chem_version.h"
TEST(version,Version,FullVersion) { CHECK(!strcmp(chem_version(),"3.14.0")); }
TEST(version,Version,VersionMajor) { CHECK(chem_version_major()==3); }
TEST(version,Version,VersionMinor) { CHECK(chem_version_minor()==14); }
TEST(version,Version,VersionPatch) { CHECK(chem_version_patch()==0); }
TEST(version,Version,VersionTweak) { CHECK(chem_version_tweak()==0); }
