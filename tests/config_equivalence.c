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

#include "tutorial_mechanisms.h"
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
static int dump_defaults(void)
{
    ChemMechanism *m=NULL;ChemState *s=NULL;
    CHECK(ts1_create(&m)==CHEM_OK);CHECK(chem_state_create(m,1,&s)==CHEM_OK);CHECK(ts1_initialize(s)==CHEM_OK);
    printf("temperature,pressure");
    for (size_t v=0;v<chem_species_count(m);++v) printf(",%s",chem_species_name(m,v));
    for (size_t v=0;v<chem_custom_count(m);++v) printf(",%s",chem_custom_name(m,v));
    printf("\n%.17g,%.17g",s->conditions[0].temperature,s->conditions[0].pressure);
    for (size_t v=0;v<chem_species_count(m);++v) printf(",%.17g",s->values[v]);
    for (size_t v=0;v<chem_custom_count(m);++v) printf(",%.17g",s->custom[v]);
    putchar('\n');chem_state_free(s);chem_mechanism_free(m);return test_failures?1:0;
}
int main(int argc,char **argv)
{
    if (argc==2 && !strcmp(argv[1],"--defaults")) return dump_defaults();
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
