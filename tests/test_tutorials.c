#include "test.h"
#include "tutorials.h"
#include "config_equivalence.h"

static void compare_program(const char *name)
{
    char path[256]; snprintf(path,sizeof(path),"tests/reference/%s.txt",name);
    FILE *expected=fopen(path,"r"),*actual=tmpfile(); REQUIRE(expected && actual);
    CHECK(example_run(name,actual)==0); rewind(actual);
    char a[4096],b[4096]; size_t line=0; unsigned announced=0;
    while (fgets(b,sizeof(b),expected)) {
        ++line; if (!fgets(a,sizeof(a),actual)) { CHECK(0); break; }
        if (!strncmp(b,"Total solve time: ",18)) {
            double duration=-1; int used=0;
            CHECK(sscanf(a,"Total solve time: %lf nanoseconds%n",&duration,&used)==1);
            CHECK(duration>=0 && isfinite(duration)); CHECK(used>0 && !strcmp(a+used,"\n"));
        } else if (!strncmp(b,"Running solver on thread ",25)) {
            unsigned thread=99; int used=0;
            CHECK(sscanf(a,"Running solver on thread %u%n",&thread,&used)==1);
            CHECK(thread<3); if (thread<3) { CHECK(!(announced & (1u<<thread))); announced|=1u<<thread; }
            CHECK(used>0 && !strcmp(a+used,"\n"));
        } else if (strcmp(a,b)) {
            fprintf(stderr,"%s:%zu expected: %sactual: %s",name,line,b,a); CHECK(0);
        }
    }
    CHECK(!fgets(a,sizeof(a),actual)); CHECK(!ferror(expected) && !ferror(actual));
    if (!strcmp(name,"openmp_tutorial")) CHECK(announced==7);
    fclose(expected); fclose(actual);
}
PROGRAM(README_example) { compare_program("README_example"); }
PROGRAM(solve_results) { compare_program("solve_results"); configuration_equivalence(EXAMPLE_ROBERTSON,&robertson_mechanism); }
PROGRAM(multiple_grid_cells) { compare_program("multiple_grid_cells"); }
PROGRAM(rate_constants_no_user_defined_example_by_hand) { compare_program("rate_constants_no_user_defined_example_by_hand"); configuration_equivalence(EXAMPLE_RATES,&rates_mechanism); }
PROGRAM(rate_constants_user_defined_example_by_hand) { compare_program("rate_constants_user_defined_example_by_hand"); configuration_equivalence(EXAMPLE_USER_RATES,&user_rates_mechanism); }
PROGRAM(solver_configuration) { compare_program("solver_configuration"); }
PROGRAM(vectorized_matrix_solver) { compare_program("vectorized_matrix_solver"); }
PROGRAM(openmp_tutorial) { compare_program("openmp_tutorial"); }
