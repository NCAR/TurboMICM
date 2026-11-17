#include "chem.h"
#include "chem_system.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int close_value(double a,double b,double tolerance)
{ return fabs(a-b)<=tolerance*fmax(1.0,fabs(b)); }

static int process_rates(size_t cells)
{
    const ChemSpecies species[]={{"foo",0,1e-3,.025,230},{"bar",0,1e-3,0,0}};
    const ChemParameterizedFunction param={.air_density=.82};
    const int reactants[]={0,0};
    const ChemReaction reactions[]={
        {0,1,0,0,0,{.kind=CHEM_ARRHENIUS,.a=12.2,.c=300,.d=300,.parameterization_count=1,.parameterizations=&param}},
        {1,0,0,0,0,{.kind=CHEM_USER,.parameter=0,.a=1,.parameterization_count=1,.parameterizations=&param}},
        {1,1,0,0,0,{.kind=CHEM_SURFACE,.parameter=1,.diffusion=230,.molecular_weight=.025,.probability=1}}};
    const char *names[]={"bar_user","foo_surf.effective radius [m]","foo_surf.particle number concentration [# m-3]"};
    const double defaults[]={0,0,0};
    const ChemMechanism m={2,3,3,species,reactions,reactants,NULL,names,defaults,285,101100};
    ChemState *s=chem_state_create(&m,cells); ChemRateIndex index={0};
    if (!s) return 0;
    int ok=chem_rate_index_init(&index,&m)==CHEM_OK && index.offsets[CHEM_USER]==1 && index.offsets[CHEM_SURFACE]==2;
    for (size_t c=0;c<cells;++c) {
        chem_temperatures(s)[c]=285+1.1*c; chem_pressures(s)[c]=101100+30*c;
        ok &= chem_set_air_density(s,c,10+.1*c)==CHEM_OK;
        chem_parameter_values(s,0)[c]=.01+.001*c; chem_parameter_values(s,1)[c]=1e-8*(1+.01*c);
        chem_parameter_values(s,2)[c]=1e5*(1+.02*c);
    }
    ok &= chem_update_rates(s)==CHEM_OK;
    const double gas=1.380649e-23*6.02214076e23,pi=3.14159265358979323846;
    for (size_t c=0;c<cells;++c) {
        double t=chem_temperatures(s)[c],density=10+.1*c,radius=chem_parameter_values(s,1)[c];
        double speed=sqrt(8*gas/(pi*.025)*t);
        double expected_arr=12.2*exp(300/t)*density*.82;
        double expected_user=chem_parameter_values(s,0)[c]*density*.82;
        double expected_surface=4*chem_parameter_values(s,2)[c]*pi*radius*radius/(radius/230+4/speed);
        ok &= fabs(chem_rate_constants(s,0)[c]-expected_arr)<=1e-12*expected_arr &&
              fabs(chem_rate_constants(s,1)[c]-expected_user)<=1e-12*expected_user &&
              fabs(chem_rate_constants(s,2)[c]-expected_surface)<=1e-10*expected_surface;
    }
    chem_rate_index_destroy(&index); chem_state_free(s); return ok;
}
static int reaction_build(int copy)
{
    const ChemSpecies species[]={{"O3",0,1e-3,0,0},{"NO",0,1e-3,0,0},{"NO2",0,1e-3,0,0},{"O2",0,1e-3,0,0}};
    const int reactants[]={0,1}; const ChemProduct products[]={{2,1},{3,1}};
    const ChemPhase phase={"gas",NULL,0};
    ChemProcessDefinition original={{0,2,0,2,0,{.kind=CHEM_ARRHENIUS,.a=12.2,.c=300,.d=300}},&phase};
    ChemProcessDefinition other=original; if (copy) original.reaction.rate.a=42;
    const ChemMechanism m={4,1,0,species,&other.reaction,reactants,products,NULL,NULL,298,101325};
    return other.reaction.reactant_count==2 && other.reaction.product_count==2 &&
           !strcmp(m.species[m.products[0].species].name,"NO2") && !strcmp(other.phase->name,"gas") &&
           !strcmp(m.species[m.reactants[0]].name,"O3") && other.reaction.rate.kind==CHEM_ARRHENIUS && other.reaction.rate.a==12.2;
}
static int explicit_set(void)
{
    const ChemSpecies species[]={{"foo",0,1e-3,0,0},{"bar",0,1e-3,0,0},{"baz",0,1e-3,0,0},
                                 {"quz",0,1e-3,0,0},{"quuz",0,1e-3,0,0},{"corge",0,1e-3,0,0}};
    const int reactants[]={0,2,1,3,2};
    const ChemProduct products[]={{1,1},{4,2.4},{0,1},{3,1.4},{1,1},{3,2.5}};
    const ChemParameterizedFunction param={.air_density=.72};
    const ChemReaction reactions[]={
        {0,2,0,2,0,{.kind=CHEM_USER,.parameter=0,.a=1}},
        {2,1,2,2,0,{.kind=CHEM_USER,.parameter=1,.a=1,.parameterization_count=1,.parameterizations=&param}},
        {3,1,4,0,0,{.kind=CHEM_USER,.parameter=2,.a=1}},
        {4,1,4,2,0,{.kind=CHEM_USER,.parameter=3,.a=1,.parameterization_count=1,.parameterizations=&param}}};
    const char *names[]={"r1","r2","r3","r4"}; const double defaults[]={0,0,0,0};
    const ChemMechanism m={6,4,4,species,reactions,reactants,products,names,defaults,300,101325};
    ChemState *s=chem_state_create(&m,2); if (!s) return 0;
    for (size_t c=0;c<2;++c) {
        for (size_t i=0;i<5;++i) chem_species_values(s,i)[c]=.1*(i+1)+c;
        for (size_t r=0;r<4;++r) chem_parameter_values(s,r)[c]=10*(r+1)+100*c;
        chem_set_air_density(s,c,70+10*c);
    }
    int ok=chem_update_rates(s)==CHEM_OK; unsigned char pattern[36];
    const unsigned char expected_pattern[36]={
        1,1,1,0,0,0, 1,1,1,0,0,0, 1,0,1,0,0,0, 0,1,1,1,0,0, 1,0,1,0,0,0, 0,0,0,0,0,0};
    ok &= chem_reaction_pattern(&m,NULL,pattern)==CHEM_OK && !memcmp(pattern,expected_pattern,sizeof pattern);
    for (size_t c=0;c<2;++c) {
        double x[6],f[6],j[36];
        for (size_t i=0;i<6;++i) { x[i]=chem_species_values(s,i)[c]; f[i]=1000; }
        for (size_t i=0;i<36;++i) j[i]=100;
        ok &= chem_add_reaction_terms(s,c,x,NULL,f,j)==CHEM_OK;
        double k1=10+100*c,k2=(20+100*c)*(70+10*c)*.72,k3=30+100*c,k4=(40+100*c)*(70+10*c)*.72;
        double expected_f[]={1000-k1*x[0]*x[2]+k2*x[1],1000+k1*x[0]*x[2]-k2*x[1]+k4*x[2],
            1000-k1*x[0]*x[2]-k4*x[2],1000+1.4*k2*x[1]-k3*x[3]+2.5*k4*x[2],1000+2.4*k1*x[0]*x[2],1000};
        double expected_j[36]; for (size_t i=0;i<36;++i) expected_j[i]=100;
        expected_j[0*6+0]+=k1*x[2]; expected_j[0*6+1]-=k2; expected_j[0*6+2]+=k1*x[0];
        expected_j[1*6+0]-=k1*x[2]; expected_j[1*6+1]+=k2; expected_j[1*6+2]-=k1*x[0]+k4;
        expected_j[2*6+0]+=k1*x[2]; expected_j[2*6+2]+=k1*x[0]+k4;
        expected_j[3*6+1]-=1.4*k2; expected_j[3*6+2]-=2.5*k4; expected_j[3*6+3]+=k3;
        expected_j[4*6+0]-=2.4*k1*x[2]; expected_j[4*6+2]-=2.4*k1*x[0];
        for (size_t i=0;i<6;++i) ok &= close_value(f[i],expected_f[i],i==1 && c==1 ? 1e-12 : 8*DBL_EPSILON);
        for (size_t i=0;i<36;++i) ok &= close_value(j[i],expected_j[i],8*DBL_EPSILON);
    }
    chem_state_free(s); return ok;
}
static int masked(void)
{
    const ChemSpecies species[]={{"A",0,1e-3,0,0},{"B",0,1e-3,0,0},{"C",0,1e-3,0,0},{"D",0,1e-3,0,0}};
    const int reactants[]={0,1}; const ChemProduct products[]={{2,1},{3,1}};
    const ChemReaction reaction={0,2,0,2,0,{.kind=CHEM_USER,.parameter=0,.a=1}};
    const char *names[]={"rate"}; const double defaults[]={0};
    const ChemMechanism m={4,1,1,species,&reaction,reactants,products,names,defaults,298,101325};
    ChemState *s=chem_state_create(&m,2); if (!s) return 0;
    chem_parameter_values(s,0)[0]=10; chem_parameter_values(s,0)[1]=20;
    int ok=chem_update_rates(s)==CHEM_OK; const unsigned char algebraic[]={0,1,0,1};
    unsigned char pattern[16]; ok &= chem_reaction_pattern(&m,algebraic,pattern)==CHEM_OK;
    for (size_t c=0;c<2;++c) {
        const double y0[]={1,2,3,4},y1[]={5,6,7,8}; const double *y=c ? y1 : y0;
        double f[]={1000,1000,1000,1000},j[16]; for (size_t i=0;i<16;++i) j[i]=500;
        ok &= chem_add_reaction_terms(s,c,y,algebraic,f,j)==CHEM_OK;
        double k=c ? 20 : 10, rate=k*y[0]*y[1];
        ok &= f[0]==1000-rate && f[1]==1000 && f[2]==1000+rate && f[3]==1000;
        for (size_t row=0;row<4;++row) for (size_t col=0;col<4;++col) {
            double expected=500;
            if ((row==0 || row==2) && col<2) expected+=(row==0 ? 1 : -1)*k*y[1-col];
            ok &= j[row*4+col]==expected && pattern[row*4+col]==((row==0 || row==2) && col<2);
        }
    }
    chem_state_free(s); return ok;
}
static int finite_difference(void)
{
    const ChemSpecies species[]={{"A",0,1e-3,0,0},{"B",0,1e-3,0,0},{"C",0,1e-3,0,0}};
    const int reactants[]={0,1,2}; const ChemProduct products[]={{2,2},{0,1}};
    const ChemReaction reactions[]={
        {0,2,0,1,0,{.kind=CHEM_USER,.parameter=0,.a=1}},{2,1,1,1,0,{.kind=CHEM_USER,.parameter=1,.a=1}}};
    const char *names[]={"r1","r2"}; const double defaults[]={0,0};
    const ChemMechanism m={3,2,2,species,reactions,reactants,products,names,defaults,298,101325};
    ChemState *s=chem_state_create(&m,2); if (!s) return 0;
    chem_parameter_values(s,0)[0]=5; chem_parameter_values(s,0)[1]=8;
    chem_parameter_values(s,1)[0]=10; chem_parameter_values(s,1)[1]=12;
    unsigned char pattern[9]; int ok=chem_update_rates(s)==CHEM_OK && chem_reaction_pattern(&m,NULL,pattern)==CHEM_OK;
    for (size_t c=0;c<2;++c) {
        double values[3]={c ? .5 : 1,c ? 4 : 2,c ? 1.5 : 3}, f[3]={0},j[9]={0};
        ok &= chem_add_reaction_terms(s,c,values,NULL,f,j)==CHEM_OK;
        for (size_t col=0;col<3;++col) {
            double plus[3]={0},minus[3]={0},unused[9]={0},h=1e-6;
            values[col]+=h; ok &= chem_add_reaction_terms(s,c,values,NULL,plus,unused)==CHEM_OK;
            values[col]-=2*h; memset(unused,0,sizeof unused); ok &= chem_add_reaction_terms(s,c,values,NULL,minus,unused)==CHEM_OK;
            values[col]+=h;
            for (size_t row=0;row<3;++row) {
                double fd=(plus[row]-minus[row])/(2*h);
                ok &= fabs(fd+j[row*3+col])<1e-7 && (fabs(fd)<1e-7 || pattern[row*3+col]);
            }
        }
    }
    chem_state_free(s); return ok;
}
static unsigned seed=987654321;
static unsigned next_random(void) { seed=1664525u*seed+1013904223u; return seed; }
static int random_system(size_t cells,size_t nr,size_t ns)
{
    ChemSpecies *species=calloc(ns,sizeof(*species)); ChemReaction *reactions=calloc(nr,sizeof(*reactions));
    int *reactants=calloc(3*nr,sizeof(*reactants)); ChemProduct *products=calloc(10*nr,sizeof(*products));
    char (*species_names)[24]=calloc(ns,sizeof(*species_names)),(*parameter_names)[24]=calloc(nr,sizeof(*parameter_names));
    const char **names=calloc(nr,sizeof(*names)); double *defaults=calloc(nr,sizeof(*defaults));
    double *y=calloc(ns,sizeof(*y)),*force=calloc(ns,sizeof(*force)),*expected=calloc(ns,sizeof(*expected)),*jac=calloc(ns*ns,sizeof(*jac));
    ChemState *s=NULL; int ok=species && reactions && reactants && products && species_names && parameter_names && names && defaults && y && force && expected && jac;
    if (!ok) goto done;
    for (size_t i=0;i<ns;++i) { snprintf(species_names[i],24,"s%zu",i); species[i]=(ChemSpecies){species_names[i],0,1e-3,0,0}; }
    size_t ri=0,pi=0;
    for (size_t r=0;r<nr;++r) {
        size_t nreact=next_random()%4,nprod=next_random()%11;
        reactions[r]=(ChemReaction){ri,nreact,pi,nprod,0,{.kind=CHEM_USER,.parameter=(int)r,.a=1}};
        for (size_t j=0;j<nreact;++j) reactants[ri++]=(int)(next_random()%ns);
        for (size_t j=0;j<nprod;++j) products[pi++]=(ChemProduct){(int)(next_random()%ns),1.2};
        snprintf(parameter_names[r],24,"k%zu",r); names[r]=parameter_names[r];
    }
    ChemMechanism m={ns,nr,nr,species,reactions,reactants,products,names,defaults,298,101325};
    s=chem_state_create(&m,cells); if (!s) { ok=0; goto done; }
    for (size_t c=0;c<cells;++c) {
        for (size_t i=0;i<ns;++i) chem_species_values(s,i)[c]=.01+(next_random()%1000)/1000.0;
        for (size_t r=0;r<nr;++r) chem_parameter_values(s,r)[c]=.01+(next_random()%1000)/1000.0;
    }
    ok=chem_update_rates(s)==CHEM_OK;
    for (size_t c=0;c<cells && ok;++c) {
        for (size_t i=0;i<ns;++i) { y[i]=chem_species_values(s,i)[c]; force[i]=expected[i]=1000; }
        memset(jac,0,ns*ns*sizeof(*jac));
        ok &= chem_add_reaction_terms(s,c,y,NULL,force,jac)==CHEM_OK;
        for (size_t r=0;r<nr;++r) {
            const ChemReaction *rx=&reactions[r]; double flux=chem_parameter_values(s,r)[c];
            for (size_t j=0;j<rx->reactant_count;++j) flux*=y[reactants[rx->reactant_start+j]];
            for (size_t j=0;j<rx->reactant_count;++j) expected[reactants[rx->reactant_start+j]]-=flux;
            for (size_t j=0;j<rx->product_count;++j) { ChemProduct p=products[rx->product_start+j]; expected[p.species]+=p.yield*flux; }
        }
        for (size_t i=0;i<ns;++i) ok &= close_value(force[i],expected[i],8*DBL_EPSILON);
    }
done:
    chem_state_free(s); free(species); free(reactions); free(reactants); free(products); free(species_names);
    free(parameter_names); free(names); free(defaults); free(y); free(force); free(expected); free(jac); return ok;
}
static int failures;
static void report(const char *id,int ok)
{ printf("%s %s %s\n",id,ok ? "PASS" : "FAIL",ok ? "reference reaction behavior matched" : "reaction behavior differed"); failures+=!ok; }
int main(void)
{
    report("process:Process.Matrix",process_rates(5));
    report("process:Process.VectorMatrix",process_rates(1) && process_rates(2) && process_rates(3) && process_rates(4) && process_rates(5));
    report("process:Process.BuildsChemicalReaction",reaction_build(0));
    report("process:Process.ChemicalReactionCopyAssignmentSucceeds",reaction_build(1));
    report("process_set:ProcessSet.Matrix",explicit_set());
    report("process_set:ProcessSet.VectorMatrix",explicit_set());
    report("process_set:RandomProcessSet.Matrix",random_system(200,50,40) && random_system(300,30,20) && random_system(400,100,80));
    report("process_set:ProcessSetAlgebraicVariables.CudaMatrix",masked());
    report("process_set:ProcessSetFiniteDifferenceJacobian.Matrix",finite_difference());
    return failures ? 1 : 0;
}
