#include "chem.h"

#include <math.h>
#include <stdio.h>

extern const ChemMechanism tutorial_no_user;
extern const ChemMechanism tutorial_user;
extern const ChemMechanism tutorial_robertson;
extern const ChemMechanism tutorial_ordered;
extern const ChemMechanism tutorial_unordered;
extern const ChemMechanism tutorial_parameterized;

#define CHECK(condition) do { if (!(condition)) { \
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); return 1; \
} } while (0)
#define FACTOR (1e-6 * 6.02214076e23)

static const ChemSpecies species[] = {
    {"A", 0, 1e-3, 0, 0}, {"B", 0, 1e-3, 0, 0}, {"C", 0, 1e-3, .025, 230},
    {"D", 0, 1e-3, 0, 0}, {"E", 0, 1e-3, 0, 0}, {"F", 0, 1e-3, 0, 0}, {"G", 0, 1e-3, 0, 0}
};
static const int reactants[] = {0, 1, 1, 2, 3, 4, 4, 5, 2, 1};
static const ChemProduct products[] = {{1,1}, {2,1}, {3,1}, {4,1}, {5,2}, {6,1}, {6,1}, {6,1}, {0,1}};
static const char *const names[] = {"C.effective radius [m]", "C.particle number concentration [# m-3]",
                                    "my photolysis rate", "my loss rate", "my emission rate"};
static const double defaults[] = {0,0,0,0,0};
static const ChemReaction reactions[] = {
    {0,1,0,1,0,{.kind=CHEM_ARRHENIUS, .a=2.15e-4, .c=110, .d=300}},
    {1,1,1,1,0,{.kind=CHEM_BRANCHED, .branch_x=1.2, .branch_y=204.3, .branch_a0=1e-3, .branch_n=2}},
    {2,1,2,1,0,{.kind=CHEM_BRANCHED, .branch_x=1.2, .branch_y=204.3, .branch_a0=1e-3, .branch_n=2, .nitrate_branch=1}},
    {3,1,3,1,0,{.kind=CHEM_SURFACE, .parameter=0, .diffusion=230, .molecular_weight=.025, .probability=.9}},
    {4,1,4,1,0,{.kind=CHEM_TERNARY, .k0_a=1.2, .k0_b=2.3, .k0_c=302.3,
                 .kinf_a=2.6/FACTOR, .kinf_b=-3.1, .kinf_c=402.1, .fc=.9, .n=1.2}},
    {5,2,5,1,0,{.kind=CHEM_TROE, .k0_a=1.8e-30*FACTOR*FACTOR, .k0_b=-3,
                 .kinf_a=2.8e-11*FACTOR, .fc=.6, .n=1}},
    {7,1,6,1,0,{.kind=CHEM_TUNNELING, .a=1.2, .b=2.3, .c=302.3}},
    {8,1,7,1,0,{.kind=CHEM_USER, .parameter=2, .a=1}},
    {9,1,8,0,0,{.kind=CHEM_USER, .parameter=3, .a=1}},
    {10,0,8,1,0,{.kind=CHEM_USER, .parameter=4, .a=1}}
};
static const ChemMechanism hand_no_user = {7,7,2,species,reactions,reactants,products,names,defaults,0,0};
static const ChemMechanism hand_user = {7,10,5,species,reactions,reactants,products,names,defaults,0,0};

static int close_enough(double a, double b)
{ return fabs(a-b) <= 1e-28 + 5e-12*fmax(fabs(a), fabs(b)); }

static int set_parameter(ChemState *s, const ChemMechanism *m, const char *name, double value)
{
    int i = chem_find_parameter(m, name);
    if (i < 0) return 0;
    chem_parameter_values(s, (size_t)i)[0] = value;
    return 1;
}

static int check_pair(const ChemMechanism *generated, const ChemMechanism *hand, int user)
{
    ChemState *a = chem_state_create(generated, 1), *b = chem_state_create(hand, 1);
    CHECK(a && b);
    chem_temperatures(a)[0] = chem_temperatures(b)[0] = 287.45;
    chem_pressures(a)[0] = chem_pressures(b)[0] = 101319.9;
    CHECK(set_parameter(a, generated, "C surface.effective radius [m]", 1e-7));
    CHECK(set_parameter(a, generated, "C surface.particle number concentration [# m-3]", 2.5e6));
    CHECK(set_parameter(b, hand, "C.effective radius [m]", 1e-7));
    CHECK(set_parameter(b, hand, "C.particle number concentration [# m-3]", 2.5e6));
    if (user) {
        CHECK(set_parameter(a, generated, "my loss rate", 1e-23));
        CHECK(set_parameter(b, hand, "my loss rate", 1e-23));
        CHECK(set_parameter(a, generated, "my emission rate", 1e-20));
        CHECK(set_parameter(b, hand, "my emission rate", 1e-20));
    }
    for (size_t i = 0; i < 7; ++i) chem_species_values(a,i)[0] = chem_species_values(b,i)[0] = (i+1)*0.01;
    CHECK(chem_update_rates(a) == CHEM_OK && chem_update_rates(b) == CHEM_OK);
    double fa[7], fb[7], ja[49], jb[49];
    CHECK(chem_evaluate(a, 0, fa, ja) == CHEM_OK && chem_evaluate(b, 0, fb, jb) == CHEM_OK);
    for (size_t i = 0; i < 7; ++i) CHECK(close_enough(fa[i], fb[i]));
    for (size_t i = 0; i < 49; ++i) CHECK(close_enough(ja[i], jb[i]));
    for (size_t i = 0; i < 7; ++i) chem_species_values(a,i)[0] = chem_species_values(b,i)[0] = i == 0 ? 1 : 0;
    ChemRosenbrock p = chem_rosenbrock_three_stage();
    for (size_t step = 0; step < 10; ++step) {
        if (user) {
            double photo = 1e-10*pow(1.5, (double)step);
            CHECK(set_parameter(a, generated, "my photolysis rate", photo));
            CHECK(set_parameter(b, hand, "my photolysis rate", photo));
        }
        CHECK(chem_update_rates(a) == CHEM_OK && chem_update_rates(b) == CHEM_OK);
        ChemStats sa, sb;
        CHECK(chem_solve(a, &p, 500, &sa) == CHEM_OK && chem_solve(b, &p, 500, &sb) == CHEM_OK);
        for (size_t i = 0; i < 7; ++i) CHECK(close_enough(chem_species_values(a,i)[0], chem_species_values(b,i)[0]));
    }
    chem_state_free(a); chem_state_free(b);
    return 0;
}

static int check_robertson(void)
{
    const ChemSpecies sp[] = {{"A",1,1e-3,0,0},{"B",0,1e-3,0,0},{"C",0,1e-3,0,0}};
    const int re[] = {0,1,1,1,2};
    const ChemProduct pr[] = {{1,1},{1,1},{2,1},{0,1},{2,1}};
    const char *const pn[] = {"r1", "r2", "r3"};
    const double pd[] = {.04,3e7,1e4};
    const ChemReaction rx[] = {
        {0,1,0,1,0,{.kind=CHEM_USER,.parameter=0,.a=1}},
        {1,2,1,2,0,{.kind=CHEM_USER,.parameter=1,.a=1}},
        {3,2,3,2,0,{.kind=CHEM_USER,.parameter=2,.a=1}}
    };
    const ChemMechanism hand = {3,3,3,sp,rx,re,pr,pn,pd,300,101325};
    ChemState *a = chem_state_create(&tutorial_robertson, 1), *b = chem_state_create(&hand, 1);
    CHECK(a && b);
    chem_temperatures(a)[0] = 300; chem_pressures(a)[0] = 101325;
    chem_species_values(a,0)[0] = 1;
    for (size_t i = 0; i < 3; ++i) CHECK(set_parameter(a, &tutorial_robertson, pn[i], pd[i]));
    CHECK(chem_update_rates(a) == CHEM_OK && chem_update_rates(b) == CHEM_OK);
    ChemRosenbrock p = chem_rosenbrock_three_stage();
    for (int i = 0; i < 10; ++i) {
        ChemStats sa, sb;
        CHECK(chem_solve(a,&p,.1,&sa) == CHEM_OK && chem_solve(b,&p,.1,&sb) == CHEM_OK);
        for (size_t j = 0; j < 3; ++j) CHECK(close_enough(chem_species_values(a,j)[0], chem_species_values(b,j)[0]));
    }
    chem_state_free(a); chem_state_free(b);
    return 0;
}

static int check_generated_ordering(void)
{
    const ChemMechanism *a=&tutorial_unordered, *b=&tutorial_ordered;
    ChemState *sa=chem_state_create(a,3), *sb=chem_state_create(b,3);
    CHECK(sa && sb && a->species_count==b->species_count);
    int mapping[3], different=0;
    for (size_t i=0;i<a->species_count;++i) {
        mapping[i]=chem_find_species(b,a->species[i].name);
        CHECK(mapping[i]>=0);
        different |= mapping[i]!=(int)i;
        for (size_t cell=0;cell<3;++cell)
            chem_species_values(sa,i)[cell]=chem_species_values(sb,(size_t)mapping[i])[cell]=.01*(i+1)*(cell+1);
    }
    CHECK(different);
    for (size_t cell=0;cell<3;++cell) {
        chem_temperatures(sa)[cell]=chem_temperatures(sb)[cell]=287.45+cell;
        chem_pressures(sa)[cell]=chem_pressures(sb)[cell]=101319.9;
    }
    CHECK(chem_update_rates(sa)==CHEM_OK && chem_update_rates(sb)==CHEM_OK);
    for (size_t cell=0;cell<3;++cell) {
        double fa[3],fb[3],ja[9],jb[9];
        CHECK(chem_evaluate(sa,cell,fa,ja)==CHEM_OK && chem_evaluate(sb,cell,fb,jb)==CHEM_OK);
        for (size_t i=0;i<3;++i) {
            CHECK(close_enough(fa[i],fb[mapping[i]]));
            for (size_t j=0;j<3;++j) CHECK(close_enough(ja[i*3+j],jb[(size_t)mapping[i]*3+(size_t)mapping[j]]));
        }
    }
    ChemRosenbrock p=chem_rosenbrock_three_stage();
    for (size_t step=0;step<5;++step) {
        ChemStats sta,stb;
        CHECK(chem_solve(sa,&p,1,&sta)==CHEM_OK && chem_solve(sb,&p,1,&stb)==CHEM_OK);
        for (size_t i=0;i<3;++i) for (size_t cell=0;cell<3;++cell)
            CHECK(close_enough(chem_species_values(sa,i)[cell],chem_species_values(sb,(size_t)mapping[i])[cell]));
    }
    chem_state_free(sa); chem_state_free(sb);
    return 0;
}

static int check_parameterized(void)
{
    const ChemSpecies species[]={{"A",0,1e-3,0,0},{"B",0,1e-3,0,0}};
    const int reactants[]={0,0}; const ChemProduct products[]={{1,1},{1,1}};
    const ChemParameterizedFunction factors[]={{1,.01,.00001,2},{0,.1,0,0},{0,.1,0,0}};
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=2,.d=300,.parameterization_count=1,.parameterizations=factors}},
        {1,1,1,1,0,{.kind=CHEM_ARRHENIUS,.a=3,.d=300,.parameterization_count=2,.parameterizations=factors+1}}};
    const ChemMechanism hand={2,2,0,species,reactions,reactants,products,NULL,NULL,300,101325};
    ChemState *a=chem_state_create(&tutorial_parameterized,2), *b=chem_state_create(&hand,2);
    CHECK(a && b && tutorial_parameterized.species_count==2);
    for (size_t i=0;i<2;++i) {
        chem_temperatures(a)[i]=chem_temperatures(b)[i]=300+10*i;
        chem_pressures(a)[i]=chem_pressures(b)[i]=101325+100*i;
        CHECK(chem_set_air_density(a,i,5+5*i)==CHEM_OK && chem_set_air_density(b,i,5+5*i)==CHEM_OK);
        chem_species_values(a,0)[i]=chem_species_values(b,0)[i]=.01*(i+1);
    }
    CHECK(chem_update_rates(a)==CHEM_OK && chem_update_rates(b)==CHEM_OK);
    for (size_t i=0;i<2;++i) {
        double t=300+10*i, pressure=101325+100*i, density=5+5*i;
        CHECK(close_enough(chem_rate_constants(a,0)[i],2*(1+.01*t+.00001*pressure+2*density)));
        CHECK(close_enough(chem_rate_constants(a,1)[i],3*(.1*t)*(.1*t)));
        double fa[2],fb[2],ja[4],jb[4];
        CHECK(chem_evaluate(a,i,fa,ja)==CHEM_OK && chem_evaluate(b,i,fb,jb)==CHEM_OK);
        for (size_t j=0;j<2;++j) CHECK(close_enough(fa[j],fb[j]));
        for (size_t j=0;j<4;++j) CHECK(close_enough(ja[j],jb[j]));
    }
    ChemRosenbrock p=chem_rosenbrock_three_stage(); ChemStats sa,sb;
    CHECK(chem_solve(a,&p,.001,&sa)==CHEM_OK && chem_solve(b,&p,.001,&sb)==CHEM_OK);
    for (size_t i=0;i<2;++i) for (size_t j=0;j<2;++j)
        CHECK(close_enough(chem_species_values(a,i)[j],chem_species_values(b,i)[j]));
    chem_state_free(a); chem_state_free(b); return 0;
}

int main(void)
{
    CHECK(check_pair(&tutorial_no_user, &hand_no_user, 0) == 0);
    CHECK(check_pair(&tutorial_user, &hand_user, 1) == 0);
    CHECK(check_robertson() == 0);
    CHECK(check_generated_ordering() == 0);
    CHECK(check_parameterized() == 0);
    puts("YAML, JSON, and Robertson configurations match hand-built mechanisms: PASS");
    return 0;
}
