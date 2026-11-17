#include "chem.h"
#include <stdio.h>
#include <string.h>

int main(void)
{
    int checks[]={strcmp(chem_version(),"3.14.0")==0,chem_version_major()==3,chem_version_minor()==14,
                  chem_version_patch()==0,chem_version_tweak()==0};
    const char *names[]={"FullVersion","VersionMajor","VersionMinor","VersionPatch","VersionTweak"};
    int failed=0;
    for (size_t i=0;i<5;++i) {
        printf("version:Version.%s %s %s\n",names[i],checks[i]?"PASS":"FAIL",checks[i]?"reference version matched":"reference version differed");
        failed+=!checks[i];
    }
    return failed ? 1 : 0;
}
