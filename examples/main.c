#include "tutorials.h"
int main(int argc,char **argv)
{
    if (argc!=2) { fprintf(stderr,"Usage: %s TUTORIAL_NAME\n",argv[0]); return 2; }
    return example_run(argv[1],stdout);
}
