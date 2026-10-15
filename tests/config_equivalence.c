/* SPDX-License-Identifier: Apache-2.0 */
#include "test.h"
#include "chem_solver.h"
#include <stdlib.h>
int test_failures;
int total_failures;
ChemStatus no_user_create(ChemMechanism **);
ChemStatus with_user_create(ChemMechanism **);
ChemStatus robertson_create(ChemMechanism **);
ChemStatus ts1_create(ChemMechanism **);
ChemStatus ts1_initialize(ChemState *);

static ChemStatus hand_rates(int user,ChemMechanism **out)
{
    ChemSpecies sp[7]={{.name="A"},{.name="B"},{.name="C",.molecular_weight=.025,.diffusion=2.3e2},
                       {.name="D"},{.name="E"},{.name="F"},{.name="G"}};
    size_t lhs[10][2]={{0},{1},{1},{2},{3},{4,4},{5},{2},{0},{1}};
    ChemProduct rhs[10][1]={{{1,1}},{{2,1}},{{3,1}},{{4,1}},{{5,2}},{{6,1}},{{6,1}},{{6,1}},{{0,1}},{{0,0}}};
    ChemRateType types[]={CHEM_ARRHENIUS,CHEM_BRANCHED,CHEM_BRANCHED,CHEM_SURFACE,CHEM_TERNARY,CHEM_TROE,CHEM_TUNNELING,CHEM_USER,CHEM_USER,CHEM_USER};
    ChemReaction r[10]={0};
    for (size_t i=0;i<10;++i) {
        r[i].reactant_count=i==8?0:i==5?2:1;r[i].reactants=lhs[i];
        r[i].product_count=i==9?0:1;r[i].products=rhs[i];r[i].rate=chem_rate_default(types[i]);
    }
    r[0].rate.data.arrhenius=(ChemArrhenius){2.15e-4,0,110,300,0};
    r[1].rate.data.branched=(ChemBranched){.X=1.2,.Y=204.3,.a0=1e-3,.n=2,.branch=0};
    r[2].rate.data.branched=r[1].rate.data.branched;r[2].rate.data.branched.branch=1;
    r[3].rate.label="C surface";r[3].rate.data.surface=(ChemSurface){.diffusion=2.3e2,.molecular_weight=.025,.probability=.9};
    double conversion=1e-6*CHEM_AVOGADRO;
    r[4].rate.data.falloff=(ChemFalloff){1.2,2.3,302.3,2.6/conversion,-3.1,402.1,.9,1.2};
    r[5].rate.data.falloff=(ChemFalloff){1.8e-30*conversion*conversion,-3,0,2.8e-11*conversion,0,0,.6,1};
    r[6].rate.data.tunneling=(ChemTunneling){1.2,2.3,302.3};
    r[7].rate.label="my photolysis rate";r[8].rate.label="my emission rate";r[9].rate.label="my loss rate";
    return chem_mechanism_create(7,sp,user?10:7,r,out);
}
static ChemStatus hand_robertson(ChemMechanism **out)
{
    ChemSpecies sp[]={{.name="A"},{.name="B"},{.name="C"}};
    const size_t lhs[][2]={{0},{1,1},{1,2}};
    const ChemProduct rhs[][2]={{{1,1}},{{1,1},{2,1}},{{0,1},{2,1}}};
    const char *labels[]={"r1","r2","r3"};
    ChemReaction r[3]={0};
    for (size_t i=0;i<3;++i) {
        r[i].reactants=lhs[i];r[i].reactant_count=i?2:1;
        r[i].products=rhs[i];r[i].product_count=i?2:1;
        r[i].rate=chem_rate_default(CHEM_USER);r[i].rate.label=labels[i];
    }
    return chem_mechanism_create(3,sp,3,r,out);
}
static void initialize(ChemState *s,int which)
{
    for (size_t c=0;c<s->cells;++c) {
        s->conditions[c]=(ChemConditions){.temperature=287.45+c,.pressure=101319.9};
        chem_ideal_air_density(&s->conditions[c]);
        for (size_t v=0;v<chem_species_count(s->mechanism);++v)
            s->values[v*s->cells+c]=v?0:1+.001*c;
    }
    if (which==2) {
        const char *labels[]={"r1","r2","r3"};double rates[]={.04,3e7,1e4};
        for (size_t i=0;i<3;++i) CHECK(chem_state_fill_custom(s,labels[i],rates[i])==CHEM_OK);
    } else {
        double radius=1e-7,number=2.5e6;
        CHECK(chem_state_fill_custom(s,"C surface.effective radius [m]",radius)==CHEM_OK);
        CHECK(chem_state_fill_custom(s,"C surface.particle number concentration [# m-3]",number)==CHEM_OK);
        if (which) {
            double photo=1e-10,emission=1e-20,loss=1e-23;
            CHECK(chem_state_fill_custom(s,"my photolysis rate",photo)==CHEM_OK);
            CHECK(chem_state_fill_custom(s,"my emission rate",emission)==CHEM_OK);
            CHECK(chem_state_fill_custom(s,"my loss rate",loss)==CHEM_OK);
        }
    }
    chem_update_rates(s);
}
static void compare(int which,size_t cells)
{
    ChemMechanism *a=NULL,*b=NULL;
    CHECK((which==0?no_user_create(&a):which==1?with_user_create(&a):robertson_create(&a))==CHEM_OK);
    CHECK((which==2?hand_robertson(&b):hand_rates(which,&b))==CHEM_OK);
    if (!a || !b) { chem_mechanism_free(a);chem_mechanism_free(b);return; }
    CHECK(chem_species_count(a)==chem_species_count(b));CHECK(chem_reaction_count(a)==chem_reaction_count(b));
    CHECK(chem_custom_count(a)==chem_custom_count(b));
    ChemState *sa=NULL,*sb=NULL;ChemSolver *as=NULL,*bs=NULL;
    CHECK(chem_state_create(a,cells,&sa)==CHEM_OK);CHECK(chem_state_create(b,cells,&sb)==CHEM_OK);
    CHECK(chem_solver_create(a,cells,&as)==CHEM_OK);CHECK(chem_solver_create(b,cells,&bs)==CHEM_OK);
    if (!sa || !sb || !as || !bs) goto finish;
    initialize(sa,which);initialize(sb,which);
    for (size_t r=0;r<chem_reaction_count(a);++r) for (size_t c=0;c<cells;++c) {
        size_t br=which==1 && (r==7 || r==8)?15-r:r;
        REL(sa->rates[r*cells+c],sb->rates[br*cells+c],1e-13);
    }
    ChemRosenbrockParameters p=chem_rosenbrock_parameters(CHEM_ROS3);
    for (size_t step=0;step<10;++step) {
        ChemSolverStats astats,bstats;
        CHECK(chem_solve(as,sa,which==2?.1:500,&p,&astats)==CHEM_OK);
        CHECK(chem_solve(bs,sb,which==2?.1:500,&p,&bstats)==CHEM_OK);
        for (size_t i=0;i<cells*chem_species_count(a);++i)
            NEAR(sa->values[i],sb->values[i],1e-30+1e-10*fabs(sb->values[i]));
    }
finish:
    chem_solver_free(as);chem_solver_free(bs);chem_state_free(sa);chem_state_free(sb);
    chem_mechanism_free(a);chem_mechanism_free(b);
}
int main(void)
{
    for (int which=0;which<3;++which) { compare(which,1);compare(which,7); }
    ChemMechanism *m=NULL;ChemState *s=NULL;
    CHECK(ts1_create(&m)==CHEM_OK);
    if (m) {
        CHECK(chem_species_count(m)==209);CHECK(chem_reaction_count(m)==547);
        CHECK(chem_state_create(m,1,&s)==CHEM_OK);
        if (s) CHECK(ts1_initialize(s)==CHEM_OK);
    }
    chem_state_free(s);chem_mechanism_free(m);
    printf("Configuration equivalence: %s\n",test_failures?"FAIL":"PASS");
    return test_failures?1:0;
}
