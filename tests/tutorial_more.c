#include "chem.h"
#include "chem_matrix.h"
#include "chem_sparse.h"

#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

static const ChemSpecies robertson_species[]={{"A",0,1e-3,0,0},{"B",0,1e-3,0,0},{"C",0,1e-3,0,0}};
static const int robertson_reactants[]={0,1,1,1,2};
static const ChemProduct robertson_products[]={{1,1},{1,1},{2,1},{0,1},{2,1}};
static const ChemReaction robertson_reactions[]={
    {0,1,0,1,0,{.kind=CHEM_USER,.parameter=0,.a=1}},
    {1,2,1,2,0,{.kind=CHEM_USER,.parameter=1,.a=1}},
    {3,2,3,2,0,{.kind=CHEM_USER,.parameter=2,.a=1}}};
static const char *const robertson_names[]={"r1","r2","r3"};
static const double robertson_defaults[]={.04,3e7,1e4};
static const ChemMechanism robertson={3,3,3,robertson_species,robertson_reactions,robertson_reactants,robertson_products,
                                     robertson_names,robertson_defaults,272.5,101253.3};
static ChemState *robertson_state(size_t cells)
{
    ChemState *s=chem_state_create(&robertson,cells);
    if (!s) return NULL;
    for (size_t c=0;c<cells;++c) if (chem_set_air_density(s,c,1e6)!=CHEM_OK) { chem_state_free(s); return NULL; }
    return s;
}
static void print_header(const ChemMechanism *m,size_t cells)
{
    printf("%6s","time"); if (cells>1) printf(",%6s","grid");
    for (size_t i=0;i<m->species_count;++i) printf(",%11s",m->species[i].name);
    putchar('\n');
}
static void print_state(ChemState *s,const ChemMechanism *m,int time)
{
    size_t cells=chem_cell_count(s);
    for (size_t c=0;c<cells;++c) {
        printf("%6d,",time); if (cells>1) printf("%6zu,",c);
        for (size_t i=0;i<m->species_count;++i) printf("%s%11.2e",i ? "," : "",chem_species_values(s,i)[c]);
        putchar('\n');
    }
}
static int advance(ChemState *s,const ChemRosenbrock *p,double dt,ChemStats *total)
{
    double elapsed=0;
    while (elapsed<dt) {
        ChemStats stats; int status=chem_update_rates(s);
        if (!status) status=chem_solve(s,p,dt-elapsed,&stats);
        if ((status!=CHEM_OK && status!=CHEM_TOO_MANY_STEPS) || stats.final_time<=0) return 0;
        elapsed+=stats.final_time;
        if (total) {
            total->function_calls+=stats.function_calls; total->jacobian_updates+=stats.jacobian_updates;
            total->internal_steps+=stats.internal_steps; total->accepted+=stats.accepted; total->rejected+=stats.rejected;
            total->decompositions+=stats.decompositions; total->solves+=stats.solves; total->final_time+=stats.final_time;
        }
    }
    return 1;
}
static int readme(void)
{
    const ChemSpecies species[]={{"Foo",20,1e-3,0,0},{"Bar",0,1e-3,0,0},{"Baz",0,1e-3,0,0}};
    const int reactants[]={0,0,1}; const ChemProduct products[]={{1,.8},{2,.2},{2,1}};
    const ChemReaction reactions[]={
        {0,1,0,2,0,{.kind=CHEM_ARRHENIUS,.a=1e-3,.d=300}},
        {1,2,2,1,0,{.kind=CHEM_ARRHENIUS,.a=1e-5,.c=110,.d=300}}};
    const ChemMechanism m={3,2,0,species,reactions,reactants,products,NULL,NULL,287.45,101319.9};
    ChemState *s=chem_state_create(&m,1); if (!s) return 1;
    ChemRosenbrock p=chem_rosenbrock_three_stage(); print_header(&m,1);
    for (int i=0;i<10;++i) {
        if (!advance(s,&p,500,NULL)) { chem_state_free(s); return 1; }
        print_state(s,&m,i*500);
    }
    chem_state_free(s); return 0;
}
static void print_stats(ChemStats s,int final_time)
{
    printf("accepted: %lu\nfunction_calls: %lu\njacobian_updates: %lu\nnumber_of_steps: %lu\n"
           "accepted: %lu\nrejected: %lu\ndecompositions: %lu\nsolves: %lu\n",
           s.accepted,s.function_calls,s.jacobian_updates,s.internal_steps,s.accepted,s.rejected,s.decompositions,s.solves);
    if (final_time) printf("final simulation time: %g\n",s.final_time);
}
static int solver_results(void)
{
    ChemState *s=robertson_state(3); if (!s) return 1;
    double initial[]={1,2,.5}; for (size_t c=0;c<3;++c) chem_species_values(s,0)[c]=initial[c];
    ChemRosenbrock p=chem_rosenbrock_three_stage(); ChemStats stats={0};
    int ok=advance(s,&p,200,&stats);
    if (ok) { puts("Solver state: Converged"); print_stats(stats,1); }
    chem_state_free(s); return ok ? 0 : 1;
}
static int multiple_cells(void)
{
    ChemState *s=robertson_state(3); if (!s) return 1;
    double initial[]={1,2,.5}; for (size_t c=0;c<3;++c) chem_species_values(s,0)[c]=initial[c];
    ChemRosenbrock p=chem_rosenbrock_three_stage(); print_header(&robertson,3); print_state(s,&robertson,0);
    for (int i=1;i<=10;++i) {
        if (!advance(s,&p,200,NULL)) { chem_state_free(s); return 1; }
        print_state(s,&robertson,200*i);
    }
    chem_state_free(s); return 0;
}
static int solver_configuration(void)
{
    const char *headings[]={"Two stages: ","Three stages: ","Four stages: ","Four stages differential algebraic: ",
                            "Six stages differential algebraic: "};
    for (int i=0;i<5;++i) {
        ChemState *s=robertson_state(1); if (!s) return 1;
        chem_species_values(s,0)[0]=1;
        ChemRosenbrock p=i==0 ? chem_rosenbrock_two_stage() : i==1 ? chem_rosenbrock_three_stage() :
                         i==2 ? chem_rosenbrock_four_stage() : i==3 ? chem_rosenbrock_four_stage_dae() :
                         chem_rosenbrock_six_stage_dae();
        if (i) putchar('\n');
        puts(headings[i]); print_header(&robertson,1); print_state(s,&robertson,0);
        ChemStats total={0}; clock_t elapsed=0;
        for (int step=1;step<=10;++step) {
            clock_t start=clock(); int ok=advance(s,&p,200,&total); elapsed+=clock()-start;
            if (!ok) { chem_state_free(s); return 1; }
            print_state(s,&robertson,step*200);
        }
        printf("Total solve time: %.0f nanoseconds\n",(double)elapsed*1e9/CLOCKS_PER_SEC);
        print_stats(total,0); chem_state_free(s);
    }
    return 0;
}
static int vectorized(void)
{
    ChemState *regular=robertson_state(3),*grouped=robertson_state(3);
    if (!regular || !grouped) { chem_state_free(regular); chem_state_free(grouped); return 1; }
    ChemMatrix a={0},b={0}; double zero=0;
    int ok=chem_matrix_init(&a,3,3,0,sizeof(double),&zero)==CHEM_OK &&
           chem_matrix_init(&b,3,3,3,sizeof(double),&zero)==CHEM_OK;
    for (size_t c=0;c<3 && ok;++c) for (size_t i=0;i<3;++i) {
        double x=1.1+c+.1*i; chem_species_values(regular,i)[c]=chem_species_values(grouped,i)[c]=x;
        *(double *)chem_matrix_at(&a,c,i)=x; *(double *)chem_matrix_at(&b,c,i)=x;
    }
    if (ok) {
        for (size_t i=0;i<a.slots;++i) printf("%g\n",((double *)a.data)[i]);
        putchar('\n');
        for (size_t i=0;i<b.slots;++i) printf("%g\n",((double *)b.data)[i]);
        putchar('\n');
        /* The reference tutorial's solve loop starts in the converged state,
         * so it updates parameters ten times without entering the loop. */
        for (int step=0;step<10;++step) ok &= chem_update_rates(regular)==CHEM_OK && chem_update_rates(grouped)==CHEM_OK;
        for (size_t c=0;c<3;++c) {
            printf("Cell %zu\n%10s%20s%20s\n",c,"Species","Regular","Vectorized");
            for (size_t i=0;i<3;++i)
                printf("%10s%20g%20g\n",robertson.species[i].name,chem_species_values(regular,i)[c],chem_species_values(grouped,i)[c]);
            putchar('\n');
        }
    }
    ChemSparseCoordinate coords[]={{0,1},{2,1},{2,3},{3,2}}; char empty[16]={0}; ChemSparseMatrix standard={0},vector={0};
    if (ok) ok=chem_sparse_init(&standard,4,3,coords,4,CHEM_SPARSE_ROW,0,sizeof(empty),empty)==CHEM_OK &&
               chem_sparse_init(&vector,4,3,coords,4,CHEM_SPARSE_ROW,3,sizeof(empty),empty)==CHEM_OK;
    if (ok) {
        for (size_t c=0;c<3;++c) for (size_t i=0;i<4;++i) {
            char str[16]; snprintf(str,sizeof(str),"%zu.%zu.%zu",c,coords[i].row,coords[i].column);
            strcpy(chem_sparse_at(&standard,c,coords[i].row,coords[i].column),str);
            strcpy(chem_sparse_at(&vector,c,coords[i].row,coords[i].column),str);
        }
        puts("\nSparse matrix standard ordering elements");
        for (size_t i=0;i<standard.storage.slots;++i) puts((char *)standard.storage.data+i*sizeof(empty));
        puts("\nSparse matrix vector ordering elements");
        for (size_t i=0;i<vector.storage.slots;++i) puts((char *)vector.storage.data+i*sizeof(empty));
    }
    chem_sparse_destroy(&standard); chem_sparse_destroy(&vector); chem_matrix_destroy(&a); chem_matrix_destroy(&b);
    chem_state_free(regular); chem_state_free(grouped); return ok ? 0 : 1;
}
static int openmp_tutorial(void)
{
    ChemState *base=robertson_state(1); if (!base) return 1;
    double result[3][3]={{0}}; int good[3]={0};
    omp_set_dynamic(0);
    #pragma omp parallel num_threads(3)
    {
        int thread=omp_get_thread_num(); ChemState *s=chem_state_create_like(base,1); int ok=s!=NULL;
        #pragma omp critical
        printf("Running solver on thread %d\n",thread);
        if (s) {
            chem_species_values(s,0)[0]=1; ok &= chem_set_air_density(s,0,1e6)==CHEM_OK;
            ChemRosenbrock p=chem_rosenbrock_three_stage();
            for (int i=0;i<10 && ok;++i) ok=advance(s,&p,200,NULL);
            for (size_t i=0;i<3;++i) result[thread][i]=chem_species_values(s,i)[0];
            chem_state_free(s);
        }
        good[thread]=ok;
    }
    int ok=1;
    for (int thread=0;thread<3;++thread) {
        ok &= good[thread];
        printf("Thread %d\n%10s,%10s,%10s\n%10.2e,%10.2e,%10.2e\n",thread+1,"A","B","C",
               result[thread][0],result[thread][1],result[thread][2]);
    }
    chem_state_free(base); return ok ? 0 : 1;
}
int main(int argc,char **argv)
{
    if (argc!=2) return 2;
    if (!strcmp(argv[1],"README_example")) return readme();
    if (!strcmp(argv[1],"solve_results")) return solver_results();
    if (!strcmp(argv[1],"multiple_grid_cells")) return multiple_cells();
    if (!strcmp(argv[1],"solver_configuration")) return solver_configuration();
    if (!strcmp(argv[1],"vectorized_matrix_solver")) return vectorized();
    if (!strcmp(argv[1],"openmp_tutorial")) return openmp_tutorial();
    return 2;
}
