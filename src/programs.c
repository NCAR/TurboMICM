/* CPU tutorial ports. Mechanism equations are data; all integration uses the
   same public C library as the benchmark. */
#include "chem.h"
#include "chem_array.h"
#include "chem_sparse.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void header(const char *const *names, size_t n, size_t cells)
{
    printf("%6s", "time");
    if (cells > 1) printf(",%6s", "grid");
    for (size_t i = 0; i < n; ++i) printf(",%11s", names[i]);
    putchar('\n');
}
static void print_state(ChemState *s, size_t n, size_t cells, double time)
{
    const double *y = chem_concentrations(s);
    for (size_t c = 0; c < cells; ++c) {
        printf("%6g,", time);
        if (cells > 1) printf("%6zu,", c);
        for (size_t i = 0; i < n; ++i) printf("%s%11.2e", i ? "," : "", y[i * cells + c]);
        putchar('\n');
    }
}
static ChemMechanism *robertson(void)
{
    static const char *const names[] = {"A", "B", "C"}, *const params[] = {"r1", "r2", "r3"};
    static const int rr[] = {0, 1, 1, 1, 2};
    static const ChemProduct pp[] = {{1, 1}, {1, 1}, {2, 1}, {0, 1}, {2, 1}};
    static ChemReaction rx[3];
    static ChemMechanismData data = {.species_count=3, .reaction_count=3, .parameter_count=3,
        .species_names=names, .parameter_names=params, .temperature=272.5, .pressure=101253.3, .reactions=rx};
    for (size_t i = 0; i < 3; ++i) {
        rx[i] = (ChemReaction){.reactant_count=i ? 2 : 1, .product_count=i ? 2 : 1,
            .reactants=rr+(i ? 2*i-1 : 0), .products=pp+(i ? 2*i-1 : 0), .rate=chem_rate_default(CHEM_USER)};
        rx[i].rate.parameter=(int)i;
    }
    ChemMechanism *m = NULL;
    if (chem_mechanism_create(&data,&m) != CHEM_OK) return NULL;
    return m;
}
static ChemState *initial_state(const ChemMechanism *m, size_t cells)
{
    ChemState *s = NULL;
    if (chem_state_create(m,cells,&s) != CHEM_OK) return NULL;
    for (size_t c = 0; c < cells; ++c) {
        chem_concentrations(s)[c] = cells == 1 ? 1. : (double[]){1.,2.,.5}[c];
        chem_air_densities(s)[c] = 1e6;
        chem_parameters(s)[c]=.04; chem_parameters(s)[cells+c]=3e7; chem_parameters(s)[2*cells+c]=1e4;
    }
    return s;
}
static void add_stats(ChemStats *sum, const ChemStats *part)
{
    sum->internal_steps+=part->internal_steps; sum->accepted+=part->accepted; sum->rejected+=part->rejected;
    sum->function_calls+=part->function_calls; sum->jacobian_updates+=part->jacobian_updates;
    sum->decompositions+=part->decompositions; sum->solves+=part->solves;
}
static int step(ChemState *s, const ChemRosenbrock *p, double dt, ChemStats *total)
{
    double elapsed = 0;
    while (elapsed < dt) {
        ChemStats stats;
        if (chem_update_rates(s) != CHEM_OK) return 1;
        ChemStatus status=chem_solve(s,p,dt-elapsed,&stats);
        if (status != CHEM_OK && status != CHEM_MAX_STEPS) return 1;
        if (stats.final_time <= 0) return 1;
        elapsed+=stats.final_time;
        if (total) add_stats(total,&stats);
    }
    return 0;
}
static void print_stats(const ChemStats *s)
{
    printf("accepted: %zu\nfunction_calls: %zu\njacobian_updates: %zu\nnumber_of_steps: %zu\naccepted: %zu\nrejected: %zu\ndecompositions: %zu\nsolves: %zu\n",
        s->accepted,s->function_calls,s->jacobian_updates,s->internal_steps,s->accepted,s->rejected,s->decompositions,s->solves);
}
static int readme(void)
{
    const char *names[]={"Foo","Bar","Baz"};
    const int rr[]={0,0,1}; const ChemProduct pp[]={{1,.8},{2,.2},{2,1}};
    ChemReaction rx[]={{.reactant_count=1,.product_count=2,.reactants=rr,.products=pp,.rate=chem_rate_default(CHEM_ARRHENIUS)},
        {.reactant_count=2,.product_count=1,.reactants=rr+1,.products=pp+2,.rate=chem_rate_default(CHEM_ARRHENIUS)}};
    rx[0].rate.p[0]=1e-3;rx[1].rate.p[0]=1e-5;rx[1].rate.p[2]=110;
    double initial[]={20,0,0};
    ChemMechanismData data={.species_count=3,.reaction_count=2,.species_names=names,.initial_concentrations=initial,
        .temperature=287.45,.pressure=101319.9,.reactions=rx};
    ChemMechanism *m=NULL;ChemState *s=NULL;
    if (chem_mechanism_create(&data,&m) || chem_state_create(m,1,&s)) return 1;
    if (chem_calculate_air_density(s)) return 1;
    ChemRosenbrock p=chem_rosenbrock_parameters(3);header(names,3,1);
    for (size_t i=0;i<10;++i) {if(step(s,&p,500,NULL))return 1;print_state(s,3,1,500.*i);}
    chem_state_free(s);chem_mechanism_free(m);return 0;
}
static int integration_program(const char *mode)
{
    ChemMechanism *m=robertson(); if(!m)return 1;
    size_t cells=!strcmp(mode,"solver_configuration") ? 1 : 3;
    int methods[]={2,3,4,-4,6};
    const char *labels[]={"Two stages: ","Three stages: ","Four stages: ","Four stages differential algebraic: ","Six stages differential algebraic: "};
    int config=!strcmp(mode,"solver_configuration"), results=!strcmp(mode,"solve_results");
    for(size_t variant=0;variant<(config?5u:1u);++variant) {
        ChemState *s=initial_state(m,cells);if(!s)return 1;
        ChemRosenbrock p=chem_rosenbrock_parameters(config?methods[variant]:3);
        if(results) {
            ChemStats stats;if(chem_update_rates(s))return 1;
            ChemStatus status=chem_solve(s,&p,200,&stats);
            printf("Solver state: %s\n",status==CHEM_OK?"Converged":chem_status_string(status));
            print_stats(&stats);printf("final simulation time: %g\n",stats.final_time);
        } else {
            if(config)printf("%s%s\n",variant?"\n":"",labels[variant]);
            header(chem_mechanism_data(m)->species_names,3,cells);print_state(s,3,cells,0);
            ChemStats total={0};double start=omp_get_wtime();
            for(size_t i=0;i<10;++i) {if(step(s,&p,200,&total))return 1;print_state(s,3,cells,200.*(i+1));}
            if(config) {printf("Total solve time: %.0f nanoseconds\n",(omp_get_wtime()-start)*1e9);print_stats(&total);}
        }
        chem_state_free(s);
    }
    chem_mechanism_free(m);return 0;
}
static int parallel_program(int threads, int check)
{
    ChemMechanism *m=robertson();if(!m)return 1;
    double values[8][3]={{0}};int failed=0;
    #pragma omp parallel num_threads(threads) reduction(|:failed)
    {
        int thread=omp_get_thread_num();
        ChemState *s=initial_state(m,1);ChemRosenbrock p=chem_rosenbrock_parameters(3);
        if(!s)failed=1;
        else {
            for(int i=0;i<10;++i)failed|=step(s,&p,200,NULL);
            memcpy(values[thread],chem_concentrations(s),3*sizeof(double));chem_state_free(s);
        }
    }
    if(check) {
        for(int t=1;t<threads;++t)for(int i=0;i<3;++i)if(values[t][i]!=values[0][i])failed=1;
    } else {
        /* Thread log lines are serialized; their order is immaterial in the
           reference, whose unsynchronized streams may even interleave bytes. */
        for(int t=0;t<threads;++t)printf("Running solver on thread %d\n",t);
        for(int t=0;t<threads;++t)printf("Thread %d\n%10s,%10s,%10s\n%10.2e,%10.2e,%10.2e\n",t+1,"A","B","C",values[t][0],values[t][1],values[t][2]);
    }
    chem_mechanism_free(m);return failed;
}
static int vector_program(void)
{
    double zero=0; ChemArray *arrays[2]={NULL,NULL};
    for(size_t policy=0;policy<2;++policy) {
        if(chem_array_create(3,3,policy?3:0,CHEM_REAL,&zero,&arrays[policy]))return 1;
        for(size_t row=0;row<3;++row)for(size_t col=0;col<3;++col) {
            double value=(double)row+1.1+.1*col;if(chem_array_set(arrays[policy],row,col,&value))return 1;
        }
        const double *data=chem_array_const_data(arrays[policy]);
        for(size_t i=0;i<9;++i)printf("%g\n",data[i]);
        putchar('\n');
    }
    /* The source's Solve loop has a false initial condition, so it updates
       parameters ten times and prints the original concentrations. */
    ChemMechanism *m=robertson();ChemState *s=initial_state(m,3);if(!s)return 1;
    for(int i=0;i<10;++i)if(chem_update_rates(s))return 1;
    for(size_t row=0;row<3;++row) {
        printf("Cell %zu\n%10s%20s%20s\n",row,"Species","Regular","Vectorized");
        for(size_t col=0;col<3;++col) {
            double a,b;chem_array_get(arrays[0],row,col,&a);chem_array_get(arrays[1],row,col,&b);
            printf("%10s%20g%20g\n",chem_mechanism_data(m)->species_names[col],a,b);
        }putchar('\n');
    }
    chem_state_free(s);chem_mechanism_free(m);
    ChemMatrixIndex coords[]={{0,1},{2,1},{2,3},{3,2}};int empty=0;
    for(size_t policy=0;policy<2;++policy) {
        ChemSparse *sparse=NULL;
        if(chem_sparse_create(4,3,policy?3:0,CHEM_CSR,CHEM_INTEGER,4,coords,&empty,&sparse))return 1;
        for(size_t cell=0;cell<3;++cell)for(size_t k=0;k<4;++k) {
            int label=(int)(cell*100+coords[k].row*10+coords[k].column);
            if(chem_sparse_set(sparse,cell,coords[k].row,coords[k].column,&label))return 1;
        }
        printf("\nSparse matrix %s ordering elements\n",policy?"vector":"standard");
        const int *data=chem_array_const_data(chem_sparse_values(sparse));
        for(size_t k=0;k<12;++k)printf("%d.%d.%d\n",data[k]/100,data[k]/10%10,data[k]%10);
        chem_sparse_free(sparse);chem_array_free(arrays[policy]);
    }
    return 0;
}
int main(int argc,char **argv)
{
    if(argc!=2)return 2;
    if(!strcmp(argv[1],"README_example"))return readme();
    if(!strcmp(argv[1],"openmp_tutorial"))return parallel_program(3,0);
    if(!strcmp(argv[1],"openmp_check"))return parallel_program(8,1);
    if(!strcmp(argv[1],"vectorized_matrix_solver"))return vector_program();
    if(!strcmp(argv[1],"multiple_grid_cells")||!strcmp(argv[1],"solve_results")||!strcmp(argv[1],"solver_configuration"))return integration_program(argv[1]);
    return 2;
}
