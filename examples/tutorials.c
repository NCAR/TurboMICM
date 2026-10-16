#include "tutorials.h"
#include "chem_process.h"
#include "chem_sparse.h"
#include <inttypes.h>
#include <omp.h>
#include <string.h>

static void reaction(example_chemistry *e,size_t i,size_t nr,size_t r1,size_t r2,
                     size_t np,size_t p1,double y1,size_t p2,double y2,chem_rate rate)
{
    e->reactants[i][0]=r1; e->reactants[i][1]=r2;
    e->products[i][0]=(chem_product){p1,y1}; e->products[i][1]=(chem_product){p2,y2};
    e->reactions[i]=(chem_reaction){nr,np,e->reactants[i],e->products[i],rate};
}
chem_status example_create(example_chemistry *e,example_kind kind)
{
    memset(e,0,sizeof(*e));
    static const char *names[]={"A","B","C","D","E","F","G"};
    static const char *readme[]={"Foo","Bar","Baz"};
    size_t n=kind>=EXAMPLE_RATES ? 7 : 3;
    for (size_t i=0;i<n;++i) e->species[i]=(chem_species){.name=kind==EXAMPLE_README ? readme[i] : names[i]};
    chem_rate r=chem_rate_default(CHEM_ARRHENIUS);
    size_t count=0;
    if (kind==EXAMPLE_README) {
        r.p.arrhenius.A=1e-3; reaction(e,0,1,0,0,2,1,.8,2,.2,r);
        r.p.arrhenius.A=1e-5; r.p.arrhenius.C=110; reaction(e,1,2,0,1,1,2,1,0,0,r); count=2;
    } else if (kind==EXAMPLE_ROBERTSON) {
        r=chem_rate_default(CHEM_USER); r.label="r1"; reaction(e,0,1,0,0,1,1,1,0,0,r);
        r.label="r2"; reaction(e,1,2,1,1,2,1,1,2,1,r);
        r.label="r3"; reaction(e,2,2,1,2,2,0,1,2,1,r); count=3;
    } else {
        const double factor=1e-6*6.02214076e23;
        e->species[2].molecular_weight=.025; e->species[2].diffusion=230;
        r.p.arrhenius.A=2.15e-4; r.p.arrhenius.C=110; reaction(e,0,1,0,0,1,1,1,0,0,r);
        r=chem_rate_default(CHEM_BRANCHED); r.p.branched=(chem_branched){.X=1.2,.Y=204.3,.a0=.001,.n=2};
        reaction(e,1,1,1,0,1,2,1,0,0,r); r.p.branched.nitrate=1; reaction(e,2,1,1,0,1,3,1,0,0,r);
        r=chem_rate_default(CHEM_SURFACE); r.label="C"; r.p.surface=(chem_surface){.diffusion=230,.molecular_weight=.025,.probability=.9};
        reaction(e,3,1,2,0,1,4,1,0,0,r);
        r=chem_rate_default(CHEM_TERNARY); r.p.falloff=(chem_falloff){1.2,2.3,302.3,2.6/factor,-3.1,402.1,.9,1.2};
        reaction(e,4,1,3,0,1,5,2,0,0,r);
        r=chem_rate_default(CHEM_TROE); r.p.falloff=(chem_falloff){1.8e-30*factor*factor,-3,0,2.8e-11*factor,0,0,.6,1};
        reaction(e,5,2,4,4,1,6,1,0,0,r);
        r=chem_rate_default(CHEM_TUNNELING); r.p.tunneling=(chem_tunneling){1.2,2.3,302.3};
        reaction(e,6,1,5,0,1,6,1,0,0,r); count=7;
        if (kind==EXAMPLE_USER_RATES) {
            r=chem_rate_default(CHEM_USER); r.label="my photolysis rate"; reaction(e,7,1,2,0,1,6,1,0,0,r);
            r.label="my emission rate"; reaction(e,8,0,0,0,1,0,1,0,0,r);
            r.label="my loss rate"; reaction(e,9,1,1,0,0,0,0,0,0,r); count=10;
        }
    }
    chem_reactions_sort(e->reactions,count);
    e->mechanism=(chem_mechanism){.n_species=n,.species=e->species,.n_reactions=count,.reactions=e->reactions};
    return chem_model_create(&e->model,&e->mechanism,CHEM_LU_DOOLITTLE);
}
chem_status example_initialize(chem_state *s,example_kind kind,int config)
{
    chem_status status=CHEM_OK;
    for (size_t c=0;c<s->cells;++c) {
        s->conditions[c]=kind==EXAMPLE_ROBERTSON ? (chem_conditions){272.5,101253.3,1e6} : (chem_conditions){287.45,101319.9,0};
        if (kind!=EXAMPLE_ROBERTSON) chem_ideal_air_density(s->conditions+c);
    }
    double initial=kind==EXAMPLE_README ? 20 : 1;
    status=chem_state_set_species(s,kind==EXAMPLE_README ? "Foo" : "A",&initial,1);
    if (status!=CHEM_OK) return status;
    const char *labels[5]; double values[5]; size_t count=0;
    if (kind==EXAMPLE_ROBERTSON) {
        labels[0]="r1"; labels[1]="r2"; labels[2]="r3";
        values[0]=.04; values[1]=3e7; values[2]=1e4; count=3;
    } else if (kind>=EXAMPLE_RATES) {
        labels[0]=config ? "C surface.effective radius [m]" : "C.effective radius [m]";
        labels[1]=config ? "C surface.particle number concentration [# m-3]" : "C.particle number concentration [# m-3]";
        values[0]=1e-7; values[1]=2.5e6; count=2;
        if (kind==EXAMPLE_USER_RATES) {
            labels[2]="my photolysis rate"; labels[3]="my emission rate"; labels[4]="my loss rate";
            values[2]=1e-10; values[3]=1e-20; values[4]=1e-20*1e-3; count=5;
        }
    }
    for (size_t i=0;i<count;++i) if ((status=chem_state_set_custom(s,labels[i],values+i,1))!=CHEM_OK) return status;
    return CHEM_OK;
}
static int print_width(const chem_state *s)
{
    size_t width=0;
    for (size_t i=0;i<s->model->mechanism->n_species;++i) {
        size_t n=strlen(s->model->mechanism->species[i].name); if (n>width) width=n;
    }
    return width<10 ? 11 : (int)width+2;
}
static void print_header(FILE *out,const chem_state *s)
{
    fprintf(out,"%6s","time"); if (s->cells>1) fprintf(out,",%6s","grid");
    for (size_t i=0;i<s->model->mechanism->n_species;++i) fprintf(out,",%*s",print_width(s),s->model->mechanism->species[i].name);
    fputc('\n',out);
}
static void print_state(FILE *out,const chem_state *s,double time)
{
    for (size_t c=0;c<s->cells;++c) {
        fprintf(out,"%6g",time); if (s->cells>1) fprintf(out,",%6zu",c);
        for (size_t i=0;i<s->model->mechanism->n_species;++i) {
            size_t v=s->model->spec_to_var[i]; fprintf(out,",%*.2e",print_width(s),s->y[v*s->cells+c]);
        }
        fputc('\n',out);
    }
}
static void print_stats(FILE *out,const chem_stats *s)
{
    fprintf(out,"accepted: %"PRIu64"\nfunction_calls: %"PRIu64"\njacobian_updates: %"PRIu64"\nnumber_of_steps: %"PRIu64
                "\naccepted: %"PRIu64"\nrejected: %"PRIu64"\ndecompositions: %"PRIu64"\nsolves: %"PRIu64"\n",
            s->accepted,s->function_calls,s->jacobian_updates,s->internal_steps,s->accepted,s->rejected,s->decompositions,s->solves);
}
static int advance(chem_state *s,double dt,const chem_ros_parameters *p,chem_stats *stats)
{
    double advanced=0;
    while (advanced<dt) {
        chem_result result=chem_ros_solve(s,dt-advanced,p,NULL,NULL); chem_state_clamp(s);
        if ((result.status!=CHEM_OK && result.status!=CHEM_MAX_STEPS) || advanced+result.stats.final_time<=advanced) return 1;
        advanced+=result.stats.final_time; if (stats) chem_stats_add(stats,&result.stats);
    }
    return 0;
}
static int trajectory(FILE *out,example_kind kind,size_t cells,int result_only,chem_ros_method method,int configuration)
{
    example_chemistry e; chem_state s={0}; int error=1;
    if (example_create(&e,kind)!=CHEM_OK) return 1;
    if (chem_state_create(&s,&e.model,cells)!=CHEM_OK || example_initialize(&s,kind,0)!=CHEM_OK) goto done;
    if (cells==3 && chem_state_assign_species(&s,"A",(double[]){1,2,.5},3)!=CHEM_OK) goto done;
    chem_ros_parameters p=chem_ros_parameters_default(method); chem_stats stats={0}; double elapsed=0;
    if (result_only) {
        chem_calculate_rates(&s); chem_result r=chem_ros_solve(&s,200,&p,NULL,NULL); chem_state_clamp(&s);
        fprintf(out,"Solver state: %s\n",chem_status_string(r.status)); print_stats(out,&r.stats);
        fprintf(out,"final simulation time: %g\n",r.stats.final_time); error=r.status!=CHEM_OK; goto done;
    }
    print_header(out,&s); if (kind!=EXAMPLE_README) print_state(out,&s,0);
    double dt=kind==EXAMPLE_ROBERTSON ? 200 : 500, photo=1e-10;
    for (size_t i=0;i<10;++i) {
        if (kind==EXAMPLE_USER_RATES && chem_state_set_custom(&s,"my photolysis rate",&photo,1)!=CHEM_OK) goto done;
        double start=omp_get_wtime(); chem_calculate_rates(&s);
        if (advance(&s,dt,&p,&stats)) goto done;
        elapsed+=omp_get_wtime()-start;
        print_state(out,&s,dt*(kind==EXAMPLE_README ? i : i+1)); photo*=1.5;
    }
    if (configuration) { fprintf(out,"Total solve time: %g nanoseconds\n",elapsed*1e9); print_stats(out,&stats); }
    error=ferror(out)!=0;
done:
    chem_state_destroy(&s); chem_model_destroy(&e.model); return error;
}
static int vectorized(FILE *out)
{
    chem_array regular={0},vector={0}; chem_sparse_builder builder={0}; chem_sparse sparse[2]={{0}}; int error=1;
    if (chem_array_create(&regular,3,3,0,CHEM_DOUBLE,0)!=CHEM_OK || chem_array_create(&vector,3,3,3,CHEM_DOUBLE,0)!=CHEM_OK) goto done;
    for (size_t c=0;c<3;++c) for (size_t j=0;j<3;++j) {
        double value=(double)c+1+(double)(j+1)/10;
        if (chem_array_set(&regular,c,j,value)!=CHEM_OK || chem_array_set(&vector,c,j,value)!=CHEM_OK) goto done;
    }
    for (size_t layout=0;layout<2;++layout) {
        const double *data=layout ? vector.data : regular.data;
        for (size_t i=0;i<9;++i) fprintf(out,"%g\n",data[i]);
        fputc('\n',out);
    }
    /* The supplied tutorial's inner solve condition is false from the start:
     * solver_state is Converged and the condition requires != Converged. */
    for (size_t c=0;c<3;++c) {
        fprintf(out,"Cell %zu\n%10s%20s%20s\n",c,"Species","Regular","Vectorized");
        const char *names[]={"A","B","C"};
        for (size_t j=0;j<3;++j) fprintf(out,"%10s%20g%20g\n",names[j],chem_array_get(&regular,c,j),chem_array_get(&vector,c,j));
        fputc('\n',out);
    }
    if (chem_sparse_builder_create(&builder,4)!=CHEM_OK) goto done;
    builder.blocks=3;
    size_t rows[]={0,2,2,3},cols[]={1,1,3,2}; char strings[12][32];
    for (size_t j=0;j<4;++j) if (chem_sparse_builder_add(&builder,rows[j],cols[j])!=CHEM_OK) goto done;
    for (size_t layout=0;layout<2;++layout) {
        if (chem_sparse_create(sparse+layout,&builder,CHEM_CSR,layout ? 3 : 0,CHEM_INDEX)!=CHEM_OK) goto done;
        size_t *handles=sparse[layout].values.data;
        for (size_t c=0;c<3;++c) for (size_t j=0;j<4;++j) {
            size_t offset; if (chem_sparse_index(sparse+layout,c,rows[j],cols[j],&offset)!=CHEM_OK) goto done;
            size_t id=c*4+j; snprintf(strings[id],sizeof(strings[id]),"%zu.%zu.%zu",c,rows[j],cols[j]); handles[offset]=id;
        }
        fprintf(out,"\nSparse matrix %s ordering elements\n",layout ? "vector" : "standard");
        for (size_t j=0;j<sparse[layout].values.size;++j) fprintf(out,"%s\n",strings[handles[j]]);
    }
    error=ferror(out)!=0;
done:
    chem_array_destroy(&regular); chem_array_destroy(&vector); chem_sparse_builder_destroy(&builder);
    for (size_t j=0;j<2;++j) chem_sparse_destroy(sparse+j);
    return error;
}
static int openmp_example(FILE *out)
{
    example_chemistry e; if (example_create(&e,EXAMPLE_ROBERTSON)!=CHEM_OK) return 1;
    double results[3][3]={{0}}; int errors[3]={0},seen[3]={0};
    #pragma omp parallel num_threads(3)
    {
        int thread=omp_get_thread_num(); seen[thread]=1;
        #pragma omp critical
        fprintf(out,"Running solver on thread %d\n",thread);
        chem_state s={0}; chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3);
        if (chem_state_create(&s,&e.model,1)!=CHEM_OK || example_initialize(&s,EXAMPLE_ROBERTSON,0)!=CHEM_OK) errors[thread]=1;
        for (size_t step=0;step<10 && !errors[thread];++step) { chem_calculate_rates(&s); errors[thread]=advance(&s,200,&p,NULL); }
        if (!errors[thread]) memcpy(results[thread],s.y,sizeof(results[thread]));
        chem_state_destroy(&s);
    }
    int error=0;
    for (size_t thread=0;thread<3;++thread) {
        error|=errors[thread] || !seen[thread];
        fprintf(out,"Thread %zu\n%10s,%10s,%10s\n%10.2e,%10.2e,%10.2e\n",thread+1,"A","B","C",results[thread][0],results[thread][1],results[thread][2]);
    }
    chem_model_destroy(&e.model); return error || ferror(out);
}
int example_run(const char *name,FILE *out)
{
    if (!strcmp(name,"README_example")) return trajectory(out,EXAMPLE_README,1,0,CHEM_ROS3,0);
    if (!strcmp(name,"solve_results")) return trajectory(out,EXAMPLE_ROBERTSON,3,1,CHEM_ROS3,0);
    if (!strcmp(name,"multiple_grid_cells")) return trajectory(out,EXAMPLE_ROBERTSON,3,0,CHEM_ROS3,0);
    if (!strcmp(name,"rate_constants_no_user_defined_example_by_hand")) return trajectory(out,EXAMPLE_RATES,1,0,CHEM_ROS3,0);
    if (!strcmp(name,"rate_constants_user_defined_example_by_hand")) return trajectory(out,EXAMPLE_USER_RATES,1,0,CHEM_ROS3,0);
    if (!strcmp(name,"vectorized_matrix_solver")) return vectorized(out);
    if (!strcmp(name,"openmp_tutorial")) return openmp_example(out);
    if (!strcmp(name,"solver_configuration")) {
        const char *headings[]={"Two stages","Three stages","Four stages","Four stages differential algebraic","Six stages differential algebraic"};
        for (size_t i=0;i<5;++i) {
            fprintf(out,"%s%s: \n",i ? "\n" : "",headings[i]);
            if (trajectory(out,EXAMPLE_ROBERTSON,1,0,(chem_ros_method)i,1)) return 1;
        }
        return ferror(out)!=0;
    }
    fprintf(stderr,"Unknown tutorial: %s\n",name); return 2;
}
