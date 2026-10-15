/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * C ports of the eight standalone CPU tutorials. */
#include "chem_solver.h"
#include "chem_matrix.h"
#include "chem_sparse_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "tutorial_mechanisms.h"

static void require(ChemStatus status)
{ if (status) { fprintf(stderr,"Tutorial: %s\n",chem_status_string(status));exit(EXIT_FAILURE); } }
typedef struct { ChemMechanism *m;ChemSolver *solver;ChemState *s; } Tutorial;
static Tutorial create(int mechanism,size_t cells)
{
    Tutorial t={0};
    if (mechanism==0) {
        ChemSpecies sp[]={{.name="Foo"},{.name="Bar"},{.name="Baz"}};
        size_t lhs[][2]={{0},{0,1}};ChemProduct rhs[][2]={{{1,.8},{2,.2}},{{2,1}}};
        ChemRate r1=chem_rate_default(CHEM_ARRHENIUS),r2=r1;r1.data.arrhenius.A=1e-3;
        r2.data.arrhenius.A=1e-5;r2.data.arrhenius.C=110;
        ChemReaction reactions[]={{1,2,lhs[0],rhs[0],r1},{2,1,lhs[1],rhs[1],r2}};
        require(chem_mechanism_create(3,sp,2,reactions,&t.m));
    } else if (mechanism==1) require(hand_robertson(&t.m));
    else require(hand_rates(mechanism==3,&t.m));
    require(chem_mechanism_reorder(t.m));require(chem_state_create(t.m,cells,&t.s));require(chem_solver_create(t.m,cells,&t.solver));
    return t;
}
static void release(Tutorial *t)
{ chem_solver_free(t->solver);chem_state_free(t->s);chem_mechanism_free(t->m); }
static double concentration(const Tutorial *t,const char *name,size_t cell)
{ size_t index;require(chem_species_index(t->m,name,&index));return t->s->values[index*t->s->cells+cell]; }
static void header(const Tutorial *t,size_t n,const char *const *names)
{
    printf("%6s","time");if (t->s->cells>1) printf(",%6s","grid");
    for (size_t i=0;i<n;++i) printf(",%11s",names[i]);
    putchar('\n');
}
static void print_state(const Tutorial *t,size_t n,const char *const *names,double time)
{
    for (size_t c=0;c<t->s->cells;++c) {
        printf("%6g,",time);if (t->s->cells>1) printf("%6zu,",c);
        for (size_t i=0;i<n;++i) printf("%s%11.2e",i?",":"",concentration(t,names[i],c));
        putchar('\n');
    }
}
static ChemStatus solve(Tutorial *t,ChemRosenbrockMethod method,double dt,ChemSolverStats *stats)
{ ChemRosenbrockParameters p=chem_rosenbrock_parameters(method);return chem_solve(t->solver,t->s,dt,&p,stats); }
static void add_stats(ChemSolverStats *a,const ChemSolverStats *b)
{
    a->internal_steps+=b->internal_steps;a->accepted+=b->accepted;a->rejected+=b->rejected;
    a->function_calls+=b->function_calls;a->jacobian_updates+=b->jacobian_updates;
    a->decompositions+=b->decompositions;a->solves+=b->solves;
}
static void interval(Tutorial *t,ChemRosenbrockMethod method,double dt,ChemSolverStats *total)
{
    double elapsed=0;
    while (elapsed<dt) {
        ChemSolverStats stats;chem_update_rates(t->s);ChemStatus status=solve(t,method,dt-elapsed,&stats);
        if (status!=CHEM_OK && status!=CHEM_MAX_STEPS) require(status);
        if (!(stats.final_time>0)) require(CHEM_STEP_TOO_SMALL);
        elapsed+=stats.final_time;if (total) add_stats(total,&stats);
    }
}
static void print_stats(const ChemSolverStats *s)
{
    printf("accepted: %zu\nfunction_calls: %zu\njacobian_updates: %zu\nnumber_of_steps: %zu\naccepted: %zu\nrejected: %zu\ndecompositions: %zu\nsolves: %zu\n",
        s->accepted,s->function_calls,s->jacobian_updates,s->internal_steps,s->accepted,s->rejected,s->decompositions,s->solves);
}
static void init_robertson(Tutorial *t,int multiple)
{
    double initial[]={1,2,.5};
    for (size_t c=0;c<t->s->cells;++c) {
        t->s->conditions[c]=(ChemConditions){272.5,101253.3,1e6};
        size_t index;require(chem_species_index(t->m,"A",&index));t->s->values[index*t->s->cells+c]=multiple?initial[c]:1;
    }
    require(chem_state_fill_custom(t->s,"r1",.04));require(chem_state_fill_custom(t->s,"r2",3e7));
    require(chem_state_fill_custom(t->s,"r3",1e4));
}
static void readme(void)
{
    Tutorial t=create(0,1);const char *names[]={"Foo","Bar","Baz"};
    t.s->conditions[0]=(ChemConditions){287.45,101319.9,0};chem_ideal_air_density(t.s->conditions);
    require(chem_state_fill_concentration(t.s,"Foo",20));header(&t,3,names);
    for (size_t i=0;i<10;++i) {
        ChemSolverStats stats;chem_update_rates(t.s);require(solve(&t,CHEM_ROS3,500,&stats));print_state(&t,3,names,i*500);
    }
    release(&t);
}
static void results(void)
{
    Tutorial t=create(1,3);init_robertson(&t,1);chem_update_rates(t.s);ChemSolverStats stats;
    ChemStatus status=solve(&t,CHEM_ROS3,200,&stats);
    printf("Solver state: %s\n",status==CHEM_OK?"Converged":chem_status_string(status));print_stats(&stats);
    printf("final simulation time: %g\n",stats.final_time);require(status);release(&t);
}
static void multiple(void)
{
    Tutorial t=create(1,3);init_robertson(&t,1);const char *names[]={"A","B","C"};
    header(&t,3,names);print_state(&t,3,names,0);
    for (size_t i=0;i<10;++i) { interval(&t,CHEM_ROS3,200,NULL);print_state(&t,3,names,200*(i+1)); }
    release(&t);
}
static void rates(int user)
{
    Tutorial t=create(user?3:2,1);const char *names[]={"A","B","C","D","E","F","G"};
    t.s->conditions[0]=(ChemConditions){287.45,101319.9,0};chem_ideal_air_density(t.s->conditions);
    require(chem_state_fill_concentration(t.s,"A",1));
    require(chem_state_fill_custom(t.s,"C surface.effective radius [m]",1e-7));
    require(chem_state_fill_custom(t.s,"C surface.particle number concentration [# m-3]",2.5e6));
    header(&t,7,names);print_state(&t,7,names,0);
    double photo=1e-10;
    if (user) { require(chem_state_fill_custom(t.s,"my emission rate",1e-20));require(chem_state_fill_custom(t.s,"my loss rate",1e-20*1e-3)); }
    for (size_t i=0;i<10;++i) {
        if (user) require(chem_state_fill_custom(t.s,"my photolysis rate",photo));
        interval(&t,CHEM_ROS3,500,NULL);print_state(&t,7,names,500*(i+1));photo*=1.5;
    }
    release(&t);
}
static double seconds(void)
{ struct timespec t;if (!timespec_get(&t,TIME_UTC)) return 0;return (double)t.tv_sec+1e-9*(double)t.tv_nsec; }
static void configuration(void)
{
    const char *titles[]={"Two stages","Three stages","Four stages","Four stages differential algebraic","Six stages differential algebraic"};
    const char *names[]={"A","B","C"};
    for (int method=0;method<5;++method) {
        Tutorial t=create(1,1);init_robertson(&t,0);ChemSolverStats stats={0};double elapsed=0;
        printf("%s%s: \n",method?"\n":"",titles[method]);header(&t,3,names);print_state(&t,3,names,0);
        for (size_t i=0;i<10;++i) {
            double start=seconds();interval(&t,(ChemRosenbrockMethod)method,200,&stats);elapsed+=seconds()-start;
            print_state(&t,3,names,200*(i+1));
        }
        printf("Total solve time: %g nanoseconds\n",elapsed*1e9);print_stats(&stats);release(&t);
    }
}
static void vectorized(void)
{
    const char *names[]={"A","B","C"};Tutorial t[2]={create(1,3),create(1,3)};ChemMatrix layouts[2]={{0}};
    for (size_t p=0;p<2;++p) {
        require(chem_matrix_create(3,3,p?3:0,CHEM_MATRIX_REAL,0,layouts+p));
        for (size_t v=0;v<3;++v) for (size_t c=0;c<3;++c) {
            double value=(double)c+1+.1*(v+1);size_t index;require(chem_species_index(t[p].m,names[v],&index));
            t[p].s->values[index*3+c]=value;require(chem_matrix_set(layouts+p,c,v,value));
        }
        for (size_t i=0;i<layouts[p].size;++i) printf("%g\n",((double *)layouts[p].data)[i]);
        putchar('\n');
        /* The original Solve loop initializes Converged and tests != Converged
         * before its first solve, so it only updates rates ten times. */
        for (size_t c=0;c<3;++c) t[p].s->conditions[c]=(ChemConditions){272.5,101253.3,1e6};
        require(chem_state_fill_custom(t[p].s,"r1",.04));require(chem_state_fill_custom(t[p].s,"r2",3e7));
        require(chem_state_fill_custom(t[p].s,"r3",1e4));
        for (size_t i=0;i<10;++i) chem_update_rates(t[p].s);
    }
    for (size_t c=0;c<3;++c) {
        printf("Cell %zu\n%10s%20s%20s\n",c,"Species","Regular","Vectorized");
        for (size_t v=0;v<3;++v) printf("%10s%20g%20g\n",names[v],concentration(t,names[v],c),concentration(t+1,names[v],c));
        putchar('\n');
    }
    ChemSparseBuilder *builder=NULL;require(chem_sparse_builder_create(4,&builder));
    size_t rows[]={0,2,2,3},columns[]={1,1,3,2};
    for (size_t i=0;i<4;++i) require(chem_sparse_builder_add(builder,rows[i],columns[i]));
    require(chem_sparse_builder_blocks(builder,3));
    for (size_t p=0;p<2;++p) {
        ChemSparseMatrix matrix={0};require(chem_sparse_matrix_build(builder,p?3:0,CHEM_SPARSE_ROW,CHEM_MATRIX_INT,&matrix));
        /* Integer IDs retain the tutorial's string labels without C++ strings. */
        for (size_t c=0;c<3;++c) for (size_t i=0;i<4;++i) require(chem_sparse_matrix_set(&matrix,c,rows[i],columns[i],100*c+10*rows[i]+columns[i]));
        printf("\nSparse matrix %s ordering elements\n",p?"vector":"standard");
        for (size_t i=0;i<matrix.values.size;++i) { int label=((int *)matrix.values.data)[i];printf("%d.%d.%d\n",label/100,label/10%10,label%10); }
        chem_sparse_matrix_free(&matrix);chem_matrix_free(layouts+p);release(t+p);
    }
    chem_sparse_builder_free(builder);
}
static void openmp_tutorial(void)
{
    Tutorial base=create(1,1);double results[3][3]={{0}};
    #pragma omp parallel num_threads(3)
    {
        int thread=omp_get_thread_num();Tutorial t={.m=base.m};
        require(chem_state_create(t.m,1,&t.s));require(chem_solver_create(t.m,1,&t.solver));init_robertson(&t,0);
        #pragma omp critical
        printf("Running solver on thread %d\n",thread);
        for (size_t i=0;i<10;++i) interval(&t,CHEM_ROS3,200,NULL);
        memcpy(results[thread],t.s->values,3*sizeof(double));chem_solver_free(t.solver);chem_state_free(t.s);
    }
    for (size_t thread=0;thread<3;++thread) {
        printf("Thread %zu\n%10s,%10s,%10s\n",thread+1,"A","B","C");
        printf("%10.2e,%10.2e,%10.2e\n",results[thread][0],results[thread][1],results[thread][2]);
    }
    release(&base);
}
int main(int argc,char **argv)
{
    if (argc!=2) return 2;
    if (!strcmp(argv[1],"README_example")) readme();
    else if (!strcmp(argv[1],"solve_results")) results();
    else if (!strcmp(argv[1],"multiple_grid_cells")) multiple();
    else if (!strcmp(argv[1],"rate_constants_no_user_defined_example_by_hand")) rates(0);
    else if (!strcmp(argv[1],"rate_constants_user_defined_example_by_hand")) rates(1);
    else if (!strcmp(argv[1],"solver_configuration")) configuration();
    else if (!strcmp(argv[1],"vectorized_matrix_solver")) vectorized();
    else if (!strcmp(argv[1],"openmp_tutorial")) openmp_tutorial();
    else return 2;
    return 0;
}
