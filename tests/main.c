#include "test.h"
#include "test_registry.inc"
int test_failures=0;
int main(int argc,char **argv)
{
    const char *filter=argc>1 ? argv[1] : "";
    FILE *file=fopen("tests/cases.txt","r");
    if (!file) return 2;
    char line[1024], id[960], kind[32];
    int failures=0, matched=0;
    while (fgets(line,sizeof(line),file)) {
        if (line[0]=='#' || sscanf(line,"%959s %31s",id,kind)!=2 || !strstr(id,filter)) continue;
        ++matched;
        size_t i;
        for (i=0; i<sizeof(cases)/sizeof(*cases); ++i) if (!strcmp(cases[i].id,id)) break;
        if (i==sizeof(cases)/sizeof(*cases)) { printf("%s FAIL not yet ported\n",id); ++failures; continue; }
        if (cases[i].na_reason) { printf("%s N/A %s\n",id,cases[i].na_reason); continue; }
        test_failures=0; cases[i].function();
        printf("%s %s %s\n",id,test_failures ? "FAIL" : "PASS",test_failures ? "assertion failure" : "reference assertions");
        if (test_failures) ++failures;
    }
    fclose(file);
    fprintf(stderr,"%d cases selected; %d failed\n",matched,failures);
    return failures || !matched ? 1 : 0;
}
