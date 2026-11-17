#include "chem.h"
#include "rosenbrock.h"
#include <stdio.h>
#include <time.h>

typedef struct { const ChemMechanism *mechanism; const double *rates; } System;
static void rhs(void *v,const double *y,double *out)
{ System *s=v; chem_rhs(s->mechanism,y,s->rates,out); }
static void jacobian(void *v,const double *y,double *out)
{ System *s=v; chem_jacobian(s->mechanism,y,s->rates,out); }
static void print_state(double time,const double *y)
{
    printf("%6g,",time);
    for (size_t i=0;i<3;++i) printf(i ? ",%11.2e" : "%11.2e",y[i]);
    putchar('\n');
}
static double nanoseconds(struct timespec start,struct timespec end)
{ return (double)(end.tv_sec-start.tv_sec)*1e9+(double)(end.tv_nsec-start.tv_nsec); }
int main(void)
{
    const size_t reactants[]={0,1,1,1,2};
    const ChemProduct products[]={{1,1},{1,1},{2,1},{0,1},{2,1}};
    const ChemReaction reactions[]={{0,1,0,1},{1,2,1,2},{3,2,3,2}};
    const ChemMechanism mechanism={3,3,reactions,reactants,products};
    const double rates[]={.04,3e7,1e4},atol[]={.001,.001,.001};
    System system={&mechanism,rates};
    ChemRosenbrockMethod methods[]={chem_rosenbrock_two_stage(),chem_rosenbrock_three_stage(),
        chem_rosenbrock_four_stage(),chem_rosenbrock_four_stage_da(),chem_rosenbrock_six_stage_da()};
    const char *labels[]={"Two stages", "Three stages", "Four stages", "Four stages differential algebraic",
                          "Six stages differential algebraic"};
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
    for (size_t method=0;method<5;++method) {
        double y[]={1,0,0},elapsed_ns=0;
        ChemSolverStats total={0};
        printf(method ? "\n%s: \n" : "%s: \n",labels[method]);
        printf("  time,          A,          B,          C\n");
        print_state(0,y);
        for (size_t step=0;step<10;++step) {
            double elapsed=0;
            while (elapsed<200) {
                struct timespec start,end;
                timespec_get(&start,TIME_UTC);
                ChemSolverStats stats;
                ChemSolveStatus status=chem_rosenbrock_solve(3,y,200-elapsed,rhs,jacobian,
                                                             &system,methods+method,&options,&stats);
                timespec_get(&end,TIME_UTC);
                elapsed_ns+=nanoseconds(start,end);
                total.function_calls+=stats.function_calls; total.jacobian_updates+=stats.jacobian_updates;
                total.steps+=stats.steps; total.accepted+=stats.accepted; total.rejected+=stats.rejected;
                total.decompositions+=stats.decompositions; total.solves+=stats.solves;
                if ((status!=CHEM_SOLVE_OK && status!=CHEM_SOLVE_MAX_STEPS) || stats.final_time<=0) {
                    fprintf(stderr,"Solver made no progress: %s\n",chem_solve_status_string(status)); return 1;
                }
                elapsed+=stats.final_time;
            }
            print_state(200*(step+1),y);
        }
        printf("Total solve time: %g nanoseconds\naccepted: %zu\nfunction_calls: %zu\n"
               "jacobian_updates: %zu\nnumber_of_steps: %zu\naccepted: %zu\nrejected: %zu\n"
               "decompositions: %zu\nsolves: %zu\n",elapsed_ns,total.accepted,total.function_calls,
               total.jacobian_updates,total.steps,total.accepted,total.rejected,total.decompositions,total.solves);
    }
    return 0;
}
