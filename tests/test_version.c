#undef NDEBUG
#include <assert.h>
#include <string.h>
#include "version.h"
int main(void)
{
    assert(!strcmp(chem_compatibility_version(),"3.14.0"));
    assert(chem_compatibility_version_major()==3 && chem_compatibility_version_minor()==14);
    assert(chem_compatibility_version_patch()==0 && chem_compatibility_version_tweak()==0);
    return 0;
}
