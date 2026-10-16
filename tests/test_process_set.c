#include "test.h"
#include "chem_process.h"
#include "chem_verify.h"
#include "test_random.h"

static void real_equal(double actual,double expected)
{
    double lower=expected,upper=expected;
    for (int i=0;i<4;++i) { lower=nextafter(lower,-INFINITY); upper=nextafter(upper,INFINITY); }
    CHECK(actual>=lower && actual<=upper);
}
static void process_set_case(size_t width)
{
    chem_species sp[]={{.name="foo"},{.name="bar"},{.name="baz"},{.name="quz"},{.name="quuz"},{.name="corge"},
        {.name="qux",.type=CHEM_SPECIES_PARAMETER,.affine={.air_density=.72}}};
    const size_t r0[]={0,2},r1[]={1,6},r2[]={3},r3[]={2,6};
    const chem_product p0[]={{1,1},{4,2.4}},p1[]={{0,1},{3,1.4}},p3[]={{1,1},{3,2.5}};
    chem_rate rate=chem_rate_default(CHEM_ARRHENIUS); rate.p.arrhenius.A=12.2; rate.p.arrhenius.C=300;
    chem_reaction rx[]={{2,2,r0,p0,rate},{2,2,r1,p1,rate},{1,0,r2,NULL,rate},{2,2,r3,p3,rate}};
    chem_mechanism mechanism={.n_species=7,.species=sp,.n_reactions=4,.reactions=rx,.preserve_order=1}; chem_model m={0};
    REQUIRE(chem_model_create(&m,&mechanism,CHEM_LU_MOZART)==CHEM_OK); CHECK(m.n_species==6);
    for (size_t i=0;i<7;++i) CHECK(chem_model_species_used(&m,i)==(i!=5));
    unsigned char pattern[36]; REQUIRE(chem_model_kinetic_pattern(&m,pattern,36)==CHEM_OK);
    const size_t pairs[][2]={{0,0},{0,1},{0,2},{1,0},{1,1},{1,2},{2,0},{2,2},{3,1},{3,2},{3,3},{4,0},{4,2}};
    size_t next=0;
    for (size_t i=0;i<6;++i) for (size_t j=0;j<6;++j) if (pattern[i*6+j]) {
        REQUIRE(next<13); CHECK(i==pairs[next][0]); CHECK(j==pairs[next][1]); ++next;
    }
    CHECK(next==13);
    const double values[][6]={{.1,.2,.3,.4,.5,0},{1.1,1.2,1.3,1.4,1.5,0}};
    if (!width) width=2;
    for (size_t start=0;start<2;start+=width) {
        size_t cells=2-start<width ? 2-start : width; chem_state s={0}; REQUIRE(chem_state_create(&s,&m,cells)==CHEM_OK);
        double f[12]; for (size_t i=0;i<6*cells;++i) f[i]=1000;
        for (size_t i=0;i<m.lu.nnz*cells;++i) s.jacobian[i]=100;
        for (size_t c=0;c<cells;++c) {
            size_t gc=start+c; double rho=gc ? 80 : 70;
            for (size_t i=0;i<6;++i) s.y[i*cells+c]=values[gc][i];
            s.rates[c]=gc ? 110 : 10; s.rates[cells+c]=(gc ? 120 : 20)*rho*.72;
            s.rates[2*cells+c]=gc ? 130 : 30; s.rates[3*cells+c]=(gc ? 140 : 40)*rho*.72;
        }
        chem_forcing(&m,cells,s.rates,s.y,f); chem_jacobian(&m,cells,s.rates,s.y,s.jacobian);
        for (size_t c=0;c<cells;++c) {
            size_t gc=start+c; double rho=gc ? 80 : 70,k1=gc ? 110 : 10,k2=gc ? 120 : 20,k3=gc ? 130 : 30,k4=gc ? 140 : 40;
            double a=values[gc][0],b=values[gc][1],z=values[gc][2],q=values[gc][3];
            double expected[]={1000-k1*a*z+k2*rho*.72*b,
                1000+k1*a*z-k2*b*rho*.72+k4*rho*.72*z,1000-k1*a*z-k4*rho*.72*z,
                1000+k2*rho*.72*b*1.4-k3*q+k4*rho*.72*2.5*z,1000+k1*a*z*2.4,1000};
            for (size_t i=0;i<6;++i) {
                if (gc==1 && i==1) NEAR(f[i*cells+c],expected[i],expected[i]*1e-12);
                else real_equal(f[i*cells+c],expected[i]);
            }
            double jac[]={100+k1*z,100-k2*rho*.72,100+k1*a,100-k1*z,100+k2*rho*.72,
                100-k1*a-k4*rho*.72,100+k1*z,100+k1*a+k4*rho*.72,100-1.4*k2*rho*.72,
                100-2.5*k4*rho*.72,100+k3,100-2.4*k1*z,100-2.4*k1*a};
            for (size_t i=0;i<13;++i) real_equal(s.jacobian[m.lu.map[pairs[i][0]*6+pairs[i][1]]*cells+c],jac[i]);
        }
        chem_state_destroy(&s);
    }
    chem_model_destroy(&m);
}
TEST(process_set,ProcessSet,Matrix) { process_set_case(0); }
TEST(process_set,ProcessSet,VectorMatrix) { for (size_t width=1;width<=4;++width) process_set_case(width); }
static void masking_case(size_t width)
{
    chem_species sp[]={{.name="A"},{.name="B"},{.name="C"},{.name="D"}}; size_t reactants[]={0,1}; chem_product products[]={{2,1},{3,1}};
    chem_reaction rx={2,2,reactants,products,chem_rate_default(CHEM_ARRHENIUS)};
    chem_mechanism mechanism={.n_species=4,.species=sp,.n_reactions=1,.reactions=&rx,.preserve_order=1}; chem_model m={0};
    REQUIRE(chem_model_create(&m,&mechanism,CHEM_LU_MOZART)==CHEM_OK); m.algebraic[1]=m.algebraic[3]=1;
    for (size_t start=0;start<2;start+=width) {
        size_t cells=2-start<width ? 2-start : width; chem_state s={0}; REQUIRE(chem_state_create(&s,&m,cells)==CHEM_OK);
        double f[8]; for (size_t i=0;i<4*cells;++i) f[i]=1000;
        for (size_t i=0;i<m.lu.nnz*cells;++i) s.jacobian[i]=500;
        for (size_t c=0;c<cells;++c) {
            for (size_t i=0;i<4;++i) s.y[i*cells+c]=1+i+4*(start+c);
            s.rates[c]=10*(start+c+1);
        }
        chem_forcing(&m,cells,s.rates,s.y,f); chem_jacobian(&m,cells,s.rates,s.y,s.jacobian);
        for (size_t c=0;c<cells;++c) {
            double rate=start+c ? 600 : 20,da=start+c ? 120 : 20,db=start+c ? 100 : 10;
            real_equal(f[c],1000-rate); real_equal(f[cells+c],1000); real_equal(f[2*cells+c],1000+rate); real_equal(f[3*cells+c],1000);
            for (size_t row=0;row<4;++row) for (size_t col=0;col<2;++col) {
                double expected=500+(row==0 ? 1 : row==2 ? -1 : 0)*(col==0 ? da : db);
                real_equal(s.jacobian[m.lu.map[row*4+col]*cells+c],expected);
            }
        }
        chem_state_destroy(&s);
    }
    chem_model_destroy(&m);
}
TEST(process_set,ProcessSetAlgebraicVariables,CudaMatrix) { for (size_t width=1;width<=4;++width) masking_case(width); }
static void fd_forcing(void *context,const double *y,double *f)
{ chem_state *s=context; chem_forcing(s->model,s->cells,s->rates,y,f); }
TEST(process_set,ProcessSetFiniteDifferenceJacobian,Matrix)
{
    chem_species sp[]={{.name="A"},{.name="B"},{.name="C"}};
    size_t r0[]={0,1},r1[]={2}; chem_product p0[]={{2,2}},p1[]={{0,1}}; chem_rate rate=chem_rate_default(CHEM_ARRHENIUS);
    chem_reaction rx[]={{2,1,r0,p0,rate},{1,1,r1,p1,rate}};
    chem_mechanism mechanism={.n_species=3,.species=sp,.n_reactions=2,.reactions=rx,.preserve_order=1}; chem_model m={0}; chem_state s={0};
    REQUIRE(chem_model_create(&m,&mechanism,CHEM_LU_MOZART)==CHEM_OK); REQUIRE(chem_state_create(&s,&m,2)==CHEM_OK);
    memcpy(s.y,(double[]){1,.5,2,4,3,1.5},6*sizeof(double)); memcpy(s.rates,(double[]){5,8,10,12},4*sizeof(double));
    chem_jacobian(&m,2,s.rates,s.y,s.jacobian); double fd[18];
    REQUIRE(chem_finite_difference(3,2,s.y,fd_forcing,&s,1e-8,fd)==CHEM_OK);
    CHECK(chem_compare_jacobian(3,2,m.lu.map,s.jacobian,fd,1e-7,1e-7).passed);
    CHECK(chem_check_sparsity(3,2,m.lu.map,fd,1e-10).passed);
    chem_state_destroy(&s); chem_model_destroy(&m);
}
static uint32_t integer(random_state *r,uint32_t upper)
{
    uint32_t range=upper+1,mask=1; while (mask<range) mask<<=1;
    uint32_t limit=2147483646u/mask*mask,value;
    do { do { value=draw(r)-1; } while (value>=limit); value&=mask-1; } while (value>=range);
    return value;
}
static void random_process_set(size_t cells,size_t reactions,size_t species)
{
    chem_species *sp=calloc(species,sizeof(*sp)); chem_reaction *rx=calloc(reactions,sizeof(*rx));
    REQUIRE(sp && rx); char (*names)[24]=calloc(species,sizeof(*names)); REQUIRE(names);
    for (size_t i=0;i<species;++i) { snprintf(names[i],24,"%zu",i); sp[i].name=names[i]; }
    random_state nr,np,ids,values; random_init(&nr,1,0); random_init(&np,1,0); random_init(&ids,1,0); random_init(&values,1,0);
    for (size_t i=0;i<reactions;++i) {
        rx[i].n_reactants=integer(&nr,3); size_t *rr=calloc(rx[i].n_reactants+1,sizeof(*rr)); REQUIRE(rr); rx[i].reactants=rr;
        for (size_t j=0;j<rx[i].n_reactants;++j) rr[j]=integer(&ids,(uint32_t)species-1);
        rx[i].n_products=integer(&np,10); chem_product *pp=calloc(rx[i].n_products+1,sizeof(*pp)); REQUIRE(pp); rx[i].products=pp;
        for (size_t j=0;j<rx[i].n_products;++j) pp[j]=(chem_product){integer(&ids,(uint32_t)species-1),1.2};
        rx[i].rate=chem_rate_default(CHEM_ARRHENIUS); rx[i].rate.p.arrhenius.A=12.2; rx[i].rate.p.arrhenius.C=300;
    }
    chem_mechanism mechanism={.n_species=species,.species=sp,.n_reactions=reactions,.reactions=rx,.preserve_order=1}; chem_model m={0}; chem_state s={0};
    REQUIRE(chem_model_create(&m,&mechanism,CHEM_LU_MOZART)==CHEM_OK); REQUIRE(chem_state_create(&s,&m,cells)==CHEM_OK);
    for (size_t c=0;c<cells;++c) for (size_t i=0;i<species;++i) s.y[i*cells+c]=lognormal(&values,4);
    for (size_t c=0;c<cells;++c) for (size_t i=0;i<reactions;++i) s.rates[i*cells+c]=lognormal(&values,4);
    double *f=malloc(species*cells*sizeof(*f)),*expected=malloc(species*sizeof(*expected)); REQUIRE(f && expected);
    for (size_t i=0;i<species*cells;++i) f[i]=1000;
    chem_forcing(&m,cells,s.rates,s.y,f);
    /* The reference checks successful execution; also compare a scalar cell-major
     * evaluator with the compiled reaction traversal, retaining accumulation order. */
    for (size_t c=0;c<cells;++c) {
        for (size_t i=0;i<species;++i) expected[i]=1000;
        for (size_t r=0;r<reactions;++r) {
            double rate=s.rates[r*cells+c];
            for (size_t j=0;j<rx[r].n_reactants;++j) rate*=s.y[rx[r].reactants[j]*cells+c];
            for (size_t j=0;j<rx[r].n_reactants;++j) expected[rx[r].reactants[j]]-=rate;
            for (size_t j=0;j<rx[r].n_products;++j) expected[rx[r].products[j].species]+=1.2*rate;
        }
        for (size_t i=0;i<species;++i) real_equal(f[i*cells+c],expected[i]);
    }
    free(f); free(expected); chem_state_destroy(&s); chem_model_destroy(&m);
    for (size_t i=0;i<reactions;++i) { free((void *)rx[i].reactants); free((void *)rx[i].products); }
    free(rx); free(sp); free(names);
}
TEST(process_set,RandomProcessSet,Matrix)
{ random_process_set(200,50,40); random_process_set(300,30,20); random_process_set(400,100,80); }
