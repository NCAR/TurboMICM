#include "chem.h"
#include "matrix.h"
#include "sparse.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static const char *names[]={"A","B","C"},*parameter_names[]={"r1","r2","r3"};
static const ChemTerm reactants[]={{0,1},{1,1},{1,1},{1,1},{2,1}},products[]={{1,1},{1,1},{2,1},{0,1},{2,1}};
static const ChemReaction reactions[]={
    {CHEM_USER,0,1,0,1,0,{1}},{CHEM_USER,1,2,1,2,1,{1}},{CHEM_USER,3,2,3,2,2,{1}}
};
static const ChemMechanism mechanism={3,3,3,names,parameter_names,NULL,NULL,reactions,reactants,products,0,NULL};
static const ChemConditions conditions={272.5,101253.3,1e6,1};
static const double rates[]={.04,3e7,1e4};
static void print_header(int cells)
{
    if(cells>1)puts("  time,  grid,          A,          B,          C");
    else puts("  time,          A,          B,          C");
}
static void print_state(int time,int cells,const double *y)
{
    for(int cell=0;cell<cells;++cell) {
        if(cells>1)printf("%6d, %5d, %10.2e, %10.2e, %10.2e\n",time,cell,y[3*cell],y[3*cell+1],y[3*cell+2]);
        else printf("%6d, %10.2e, %10.2e, %10.2e\n",time,y[0],y[1],y[2]);
    }
}
static int reference_readme(void)
{
    const char *sp[]={"Foo","Bar","Baz"};const ChemTerm rs[]={{0,1},{0,1},{1,1}},ps[]={{1,.8},{2,.2},{2,1}};
    const ChemReaction rx[]={{CHEM_ARRHENIUS,0,1,0,2,-1,{1e-3,0,0,300}},{CHEM_ARRHENIUS,1,2,2,1,-1,{1e-5,0,110,300}}};
    const ChemMechanism m={3,2,0,sp,NULL,NULL,NULL,rx,rs,ps,0,NULL};ChemConditions c={287.45,101319.9,0,0};
    ChemWorkspace *w=chem_workspace_create(&m);double y[]={20,0,0},k[2];ChemRosenbrockParameters p=chem_rosenbrock_three_stage();
    if(!w || chem_rate_constants(&m,c,NULL,k))return 1;
    puts("  time,        Foo,        Bar,        Baz");
    for(int i=0;i<10;++i) {ChemStatistics s;if(chem_rosenbrock_solve(&m,w,&p,c,k,500,y,&s))return 1;
        printf("%6d, %10.2e, %10.2e, %10.2e\n",i*500,y[0],y[1],y[2]);}
    chem_workspace_free(w);return 0;
}
static int multi(int results_only)
{
    ChemWorkspace *w[3]={0};ChemConditions c[]={conditions,conditions,conditions};double k[]={.04,3e7,1e4,.04,3e7,1e4,.04,3e7,1e4};
    double y[]={1,0,0,2,0,0,.5,0,0};for(int i=0;i<3;++i)if(!(w[i]=chem_workspace_create(&mechanism)))return 1;
    ChemRosenbrockParameters p=chem_rosenbrock_three_stage();ChemStatistics stats;
    if(!results_only){print_header(3);print_state(0,3,y);}
    for(int step=0;step<(results_only?1:10);++step) {
        if(chem_rosenbrock_solve_group(&mechanism,w,3,&p,c,k,200,y,&stats))return 1;
        if(!results_only)print_state(200*(step+1),3,y);
    }
    if(results_only) {
        puts("Solver state: Converged");printf("accepted: %lu\nfunction_calls: %lu\njacobian_updates: %lu\nnumber_of_steps: %lu\naccepted: %lu\nrejected: %lu\ndecompositions: %lu\nsolves: %lu\nfinal simulation time: %g\n",stats.accepted,stats.function_calls,stats.jacobian_updates,stats.internal_steps,stats.accepted,stats.rejected,stats.decompositions,stats.solves,stats.final_time);
    }
    for(int i=0;i<3;++i)chem_workspace_free(w[i]);
    return 0;
}
static void add_stats(ChemStatistics *a,ChemStatistics b)
{
    a->accepted+=b.accepted;a->function_calls+=b.function_calls;a->jacobian_updates+=b.jacobian_updates;
    a->internal_steps+=b.internal_steps;a->rejected+=b.rejected;a->decompositions+=b.decompositions;a->solves+=b.solves;
}
static int configuration(void)
{
    const char *titles[]={"Two stages: ","Three stages: ","Four stages: ","Four stages differential algebraic: ","Six stages differential algebraic: "};
    for(int method=0;method<5;++method) {
        if(method)putchar('\n');
        puts(titles[method]);
        ChemRosenbrockParameters p=method==0?chem_rosenbrock_two_stage():method==1?chem_rosenbrock_three_stage():method==2?chem_rosenbrock_four_stage():method==3?chem_rosenbrock_four_stage_dae():chem_rosenbrock_six_stage_dae();
        ChemWorkspace *w=chem_workspace_create(&mechanism);if(!w)return 1;double y[]={1,0,0};ChemStatistics total={0};clock_t elapsed=0;
        print_header(1);print_state(0,1,y);
        for(int step=0;step<10;++step){ChemStatistics stats;clock_t start=clock();
            if(chem_rosenbrock_solve(&mechanism,w,&p,conditions,rates,200,y,&stats))return 1;
            elapsed+=clock()-start;add_stats(&total,stats);print_state(200*(step+1),1,y);}
        double ns=1e9*(double)elapsed/CLOCKS_PER_SEC;
        printf("Total solve time: %.0f nanoseconds\naccepted: %lu\nfunction_calls: %lu\njacobian_updates: %lu\nnumber_of_steps: %lu\naccepted: %lu\nrejected: %lu\ndecompositions: %lu\nsolves: %lu\n",ns,total.accepted,total.function_calls,total.jacobian_updates,total.internal_steps,total.accepted,total.rejected,total.decompositions,total.solves);
        chem_workspace_free(w);
    }
    return 0;
}
static int vectorized(void)
{
    ChemMatrix a={0},b={0};if(chem_matrix_init(&a,3,3,1,0,0)||chem_matrix_init(&b,3,3,3,1,0))return 1;
    for(size_t r=0;r<3;++r)for(size_t c=0;c<3;++c)*chem_matrix_at(&a,r,c)=*chem_matrix_at(&b,r,c)=1.1+r+.1*c;
    for(size_t i=0;i<a.storage_size;++i)printf("%g\n",a.data[i]);
    putchar('\n');
    for(size_t i=0;i<b.storage_size;++i)printf("%g\n",b.data[i]);
    putchar('\n');
    for(size_t r=0;r<3;++r){printf("Cell %zu\n",r);printf("%10s%20s%20s\n","Species","Regular","Vectorized");
        for(size_t c=0;c<3;++c)printf("%10s%20g%20g\n",names[c],*chem_matrix_at(&a,r,c),*chem_matrix_at(&b,r,c));
        putchar('\n');}
    ChemSparseBuilder builder={0};ChemSparseMatrix s={0},v={0};if(chem_sparse_builder_init(&builder,4))return 1;
    size_t points[][2]={{0,1},{2,1},{2,3},{3,2}};for(size_t i=0;i<4;++i)chem_sparse_builder_add(&builder,points[i][0],points[i][1]);
    if(chem_sparse_init(&s,&builder,3,1,0,0)||chem_sparse_init(&v,&builder,3,3,0,0))return 1;
    chem_sparse_builder_free(&builder);
    for(size_t cell=0;cell<3;++cell)for(size_t i=0;i<4;++i)*chem_sparse_at(&s,cell,points[i][0],points[i][1])=*chem_sparse_at(&v,cell,points[i][0],points[i][1])=(double)(100*cell+10*points[i][0]+points[i][1]);
    puts("\nSparse matrix standard ordering elements");
    for(size_t i=0;i<s.storage_size;++i){int x=(int)s.data[i];printf("%d.%d.%d\n",x/100,(x/10)%10,x%10);}
    puts("\nSparse matrix vector ordering elements");
    for(size_t i=0;i<v.storage_size;++i){int x=(int)v.data[i];printf("%d.%d.%d\n",x/100,(x/10)%10,x%10);}
    chem_sparse_free(&s);chem_sparse_free(&v);chem_matrix_free(&a);chem_matrix_free(&b);return 0;
}
static int openmp_example(void)
{
    double result[3][3]={{0}};int status[3]={0};
    omp_set_dynamic(0);
#pragma omp parallel num_threads(3)
    {
        int i=omp_get_thread_num();
#pragma omp critical
        printf("Running solver on thread %d\n",i);
        ChemWorkspace *w=chem_workspace_create(&mechanism);ChemRosenbrockParameters p=chem_rosenbrock_three_stage();
        double y[]={1,0,0};if(!w)status[i]=1;
        for(int step=0;step<10 && !status[i];++step){ChemStatistics stats;if(chem_rosenbrock_solve(&mechanism,w,&p,conditions,rates,200,y,&stats))status[i]=1;}
        memcpy(result[i],y,sizeof y);chem_workspace_free(w);
    }
    for(int i=0;i<3;++i){if(status[i])return 1;printf("Thread %d\n",i+1);puts("         A,         B,         C");printf("%10.2e,%10.2e,%10.2e\n",result[i][0],result[i][1],result[i][2]);}
    return 0;
}
static int openmp_unit(void)
{
    double result[8][3]={{0}};int status[8]={0};int threads=0;
    const ChemRosenbrockParameters p=chem_rosenbrock_three_stage();
    omp_set_dynamic(0);
#pragma omp parallel num_threads(8)
    {
        int i=omp_get_thread_num();
#pragma omp single
        threads=omp_get_num_threads();
        ChemWorkspace *w=chem_workspace_create(&mechanism);double y[]={1,0,0};
        if(!w)status[i]=1;
        for(int step=0;step<10 && !status[i];++step){ChemStatistics stats;if(chem_rosenbrock_solve(&mechanism,w,&p,conditions,rates,200,y,&stats))status[i]=1;}
        memcpy(result[i],y,sizeof y);chem_workspace_free(w);
    }
    if(threads!=8)return 1;
    for(int i=0;i<8;++i){if(status[i])return 1;for(int j=0;j<3;++j)if(result[i][j]!=result[0][j])return 1;}
    return 0;
}
int main(int argc,char **argv)
{
    if(argc!=2)return 1;
    if(!strcmp(argv[1],"README_example"))return reference_readme();
    if(!strcmp(argv[1],"solve_results"))return multi(1);
    if(!strcmp(argv[1],"multiple_grid_cells"))return multi(0);
    if(!strcmp(argv[1],"solver_configuration"))return configuration();
    if(!strcmp(argv[1],"vectorized_matrix_solver"))return vectorized();
    if(!strcmp(argv[1],"openmp_tutorial"))return openmp_example();
    if(!strcmp(argv[1],"openmp_unit"))return openmp_unit();
    return 1;
}
