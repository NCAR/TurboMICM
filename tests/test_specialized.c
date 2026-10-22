/* SPDX-License-Identifier: Apache-2.0 */
#include "solver_internal.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define DECLARE(name) \
ChemStatus name##_create(ChemMechanism **); \
ChemStatus name##_initialize(ChemState *); \
ChemStatus name##_specialized_create(ChemMechanism **); \
void name##_specialized_attach(ChemSolver *); \
void name##_update_rates(ChemState *); \
ChemStatus name##_solve(ChemSolver *,ChemState *,double,const ChemRosenbrockParameters *,ChemSolverStats *);
DECLARE(ts1)
DECLARE(no_user)
DECLARE(with_user)
DECLARE(robertson)
DECLARE(extra)
DECLARE(mutated)
typedef struct {
    const char *name;
    ChemStatus (*create)(ChemMechanism **),(*initialize)(ChemState *),(*specialized_create)(ChemMechanism **);
    void (*attach)(ChemSolver *),(*rates)(ChemState *);
    ChemStatus (*solve)(ChemSolver *,ChemState *,double,const ChemRosenbrockParameters *,ChemSolverStats *);
} Case;
#define CASE(n) {#n,n##_create,n##_initialize,n##_specialized_create,n##_specialized_attach,n##_update_rates,n##_solve}
static int failures;
static void compare(const char *what,const double *a,const double *b,size_t n)
{
    for (size_t i=0;i<n;++i) if (!isfinite(a[i]) || !isfinite(b[i]) || fabs(a[i]-b[i])>1e-10*fabs(b[i])) {
        fprintf(stderr,"%s[%zu]: %.17g != %.17g\n",what,i,a[i],b[i]); ++failures; return;
    }
}
#define REQUIRE(expr) do { if (!(expr)) { fprintf(stderr,"line %d: %s\n",__LINE__,#expr); ++failures; goto finish; } } while (0)
static void run(const Case *test,size_t cells,int tight)
{
    ChemMechanism *m=NULL,*generic=NULL; ChemState *a=NULL,*b=NULL; ChemSolver *sa=NULL,*sb=NULL;
    REQUIRE(test->specialized_create(&m)==CHEM_OK);
    REQUIRE(test->create(&generic)==CHEM_OK);
    REQUIRE(chem_mechanism_reorder(generic)==CHEM_OK);
    REQUIRE(chem_state_create(m,cells,&a)==CHEM_OK);
    REQUIRE(test->initialize(a)==CHEM_OK);
    REQUIRE(chem_solver_create(m,cells,&sa)==CHEM_OK);
    REQUIRE(chem_solver_create(generic,cells,&sb)==CHEM_OK);
    test->attach(sa);
    size_t nv=chem_species_count(m),nr=chem_reaction_count(m),nc=chem_custom_count(m);
    for (size_t v=0;v<nv;++v) REQUIRE(!strcmp(chem_species_name(m,v),chem_species_name(generic,v)));
    for (size_t c=0;c<cells;++c) {
        a->conditions[c]=(ChemConditions){270.0+1.3*c,80000+731.0*c,0};
        chem_ideal_air_density(a->conditions+c);
        /* Air density is independent runtime data after initialization. */
        a->conditions[c].air_density *= .95;
        for (size_t v=0;v<nv;++v) {
            if (strcmp(test->name,"ts1")) a->values[v*cells+c]=1e-8*(v+1)*(1+.01*c);
            else a->values[v*cells+c]*=1+.01*c;
        }
        for (size_t j=0;j<nc;++j) {
            const char *name=chem_custom_name(m,j);
            double value=a->custom[j*cells+c];
            if (!strcmp(test->name,"ts1") && (!strcmp(name,"jo2_a") || !strcmp(name,"jh2o_b") || !strcmp(name,"jh2o_c") || !strcmp(name,"jch4_a") || !strcmp(name,"jch4_b") || !strcmp(name,"jhf"))) value=0;
            if (strstr(name,"effective radius")) value=1e-7;
            else if (strstr(name,"particle number")) value=2.5e6;
            else if (strcmp(test->name,"ts1")) value=1e-5*(j+1);
            a->custom[j*cells+c]=value*(1+.013*c);
        }
    }
    if (!strcmp(test->name,"robertson")) {
        REQUIRE(chem_state_fill_custom(a,"r1",.04)==CHEM_OK);
        REQUIRE(chem_state_fill_custom(a,"r2",3e7)==CHEM_OK);
        REQUIRE(chem_state_fill_custom(a,"r3",1e4)==CHEM_OK);
        for (size_t c=0;c<cells;++c) a->values[c]=1;
    }
    if (tight) for (size_t v=0;v<nv;++v) a->absolute_tolerance[v]=1e-15;
    REQUIRE(chem_state_copy(a,&b)==CHEM_OK); b->mechanism=generic;
    test->rates(a); chem_update_rates(b);
    compare("rates",a->rates,b->rates,nr*cells);
    /* Exercise cache hits and each independent invalidation input. */
    test->rates(a); compare("cached rates",a->rates,b->rates,nr*cells);
    for (int input=0;input<4;++input) {
        if (input==0) a->conditions[0].temperature+=.25;
        else if (input==1) a->conditions[0].pressure*=.99;
        else if (input==2) a->conditions[0].air_density*=.98;
        else for (size_t j=0;j<nc;++j) a->custom[j*cells]*=1.01;
        memcpy(b->conditions,a->conditions,cells*sizeof(ChemConditions));
        memcpy(b->custom,a->custom,nc*cells*sizeof(double));
        test->rates(a); chem_update_rates(b);
        compare("invalidated rates",a->rates,b->rates,nr*cells);
    }
    chem_update_rates(a); test->rates(a);
    compare("generic cache invalidation",a->rates,b->rates,nr*cells);
    /* Test the whole symbolic factor/solve plan independently of chemistry. */
    size_t nnz=chem_sparse_elements(sa->sparse);
    REQUIRE(nnz==chem_sparse_elements(sb->sparse));
    for (size_t i=0;i<nnz*cells;++i) sa->matrix[i]=.00001*(1+i%17);
    for (size_t v=0;v<nv;++v) for (size_t c=0;c<cells;++c) {
        sa->matrix[chem_sparse_index(sa->sparse,v,v)*cells+c]+=10;
        sa->initial[v*cells+c]=.1*(v+1)+c;
    }
    memcpy(sb->matrix,sa->matrix,nnz*cells*sizeof(double)); memcpy(sb->initial,sa->initial,nv*cells*sizeof(double));
    REQUIRE(sa->linear.factor(sa->sparse,cells,sa->matrix,NULL)==CHEM_OK);
    REQUIRE(chem_lu_factor(sb->sparse,cells,sb->matrix,CHEM_LU_MOZART)==CHEM_OK);
    compare("LU",sa->matrix,sb->matrix,nnz*cells);
    sa->linear.solve(sa->sparse,cells,sa->matrix,sa->initial,NULL);
    chem_lu_solve(sb->sparse,cells,sb->matrix,sb->initial);
    compare("triangular",sa->initial,sb->initial,nv*cells);
    if (nv) {
        size_t diagonal=chem_sparse_index(sa->sparse,0,0)*cells;
        sa->matrix[diagonal]=0;
        REQUIRE(sa->linear.factor(sa->sparse,cells,sa->matrix,NULL)==CHEM_SINGULAR);
        sa->matrix[diagonal]=NAN;
        REQUIRE(sa->linear.factor(sa->sparse,cells,sa->matrix,NULL)==CHEM_SINGULAR);
    }
    ChemRosenbrockParameters parameters=chem_rosenbrock_parameters(CHEM_ROS3);
    for (size_t step=0;step<12;++step) {
        a->relative_tolerance=b->relative_tolerance=1e-6*(1+.01*step);
        for (size_t c=0;c<cells;++c) { a->conditions[c].temperature+=.1; b->conditions[c]=a->conditions[c]; }
        for (size_t j=0;j<nc;++j) if (!strstr(chem_custom_name(m,j),"radius") && !strstr(chem_custom_name(m,j),"concentration"))
            for (size_t c=0;c<cells;++c) { a->custom[j*cells+c]*=.97; b->custom[j*cells+c]=a->custom[j*cells+c]; }
        test->rates(a); chem_update_rates(b); compare("updated rates",a->rates,b->rates,nr*cells);
        ChemSolverStats as,bs;
        double dt=!strcmp(test->name,"ts1")?30:.1;
        ChemStatus ast=test->solve(sa,a,dt,&parameters,&as), bst=chem_solve(sb,b,dt,&parameters,&bs);
        if (ast || bst) fprintf(stderr,"%s cells=%zu tight=%d step=%zu: specialized=%s generic=%s\n",test->name,cells,tight,step,chem_status_string(ast),chem_status_string(bst));
        REQUIRE(ast==CHEM_OK && bst==CHEM_OK);
        compare("integration",a->values,b->values,nv*cells);
        REQUIRE(as.internal_steps==bs.internal_steps && as.rejected==bs.rejected && as.solves==bs.solves);
    }
finish:
    chem_solver_free(sa); chem_solver_free(sb); chem_state_free(a); chem_state_free(b);
    chem_mechanism_free(m); chem_mechanism_free(generic);
}
int main(void)
{
    const Case cases[]={CASE(ts1),CASE(no_user),CASE(with_user),CASE(robertson),CASE(extra),CASE(mutated)};
    const size_t sizes[]={1,7,16,32,33};
    for (size_t i=0;i<sizeof(cases)/sizeof(*cases);++i) {
        int before=failures;
        for (size_t j=0;j<sizeof(sizes)/sizeof(*sizes);++j) for (int tight=0;tight<2;++tight) run(cases+i,sizes[j],tight);
        printf("%s specialized equivalence: %s\n",cases[i].name,failures==before?"PASS":"FAIL");
    }
    return failures?1:0;
}
