#include "chem.h"

#include <float.h>
#include <math.h>
#include <stdio.h>

extern const ChemMechanism constrained_original, constrained_ordered;
extern const ChemConstraint constrained_original_constraints[], constrained_ordered_constraints[];
extern const size_t constrained_original_constraint_count, constrained_ordered_constraint_count;

static int failures;
static void report(const char *name,int ok)
{
    printf("equilibrium_constraint:EquilibriumIntegration.%s %s %s\n",name,ok ? "PASS" : "FAIL",
           ok ? "reference quantities matched" : "reference quantity differed");
    failures+=!ok;
}
static const ChemSpecies species[]={{"A",0,1e-3,0,0},{"B",0,1e-3,0,0},{"C",0,1e-3,0,0},
                                    {"D",0,1e-3,0,0},{"E",0,1e-3,0,0},{"F",0,1e-3,0,0}};
static const ChemTerm a2[]={{0,2}}, b[]={{1,1}}, c[]={{2,1}}, d[]={{3,1}}, e[]={{4,1}}, f[]={{5,1}};
static double expected_constant(double k, double enthalpy, double temperature)
{
    double gas=1.380649e-23*6.02214076e23;
    return k*exp((enthalpy/gas)*(1/temperature-1/298.15));
}
static ChemMechanism mechanism(size_t n, const ChemReaction *reactions, size_t nreactions,
                               const int *reactants, const ChemProduct *products)
{ return (ChemMechanism){n,nreactions,0,species,reactions,reactants,products,NULL,NULL,298.15,101325}; }

static int single_api(void)
{
    const int reactants[]={0}; const ChemProduct products[]={{1,1}};
    const ChemReaction rx[]={{0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=.1,.d=300}}};
    ChemMechanism m=mechanism(3,rx,1,reactants,products);
    ChemConstraint eq={"B_C_eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,b,c,1,1,0,.034,-2400,0,0};
    ChemState *s=chem_state_create_constrained(&m,&eq,1,1); if (!s) return 0;
    int ok=chem_find_species(&m,"A")==0 && chem_find_species(&m,"B")==1 && chem_find_species(&m,"C")==2 &&
           chem_find_constraint(s,"B_C_eq")==0 && !chem_is_algebraic(s,0) && !chem_is_algebraic(s,1) && chem_is_algebraic(s,2);
    chem_species_values(s,0)[0]=1; chem_species_values(s,1)[0]=.1; chem_species_values(s,2)[0]=.0034;
    chem_temperatures(s)[0]=300;
    double expected=expected_constant(.034,-2400,300);
    ok &= chem_update_rates(s)==CHEM_OK && fabs(chem_constraint_constants(s,0)[0]-expected)<=fabs(expected)*100*DBL_EPSILON;
    chem_state_free(s); return ok;
}

static int multiple_api(void)
{
    const int reactants[]={0,3}; const ChemProduct products[]={{1,1},{4,1}};
    const ChemReaction rx[]={{0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=.5,.d=300}},
                             {1,1,1,1,0,{.kind=CHEM_ARRHENIUS,.a=.2,.d=300}}};
    ChemMechanism m=mechanism(6,rx,2,reactants,products);
    ChemConstraint eq[]={{"B_C_eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,b,c,1,1,0,.034,-2400,0,0},
                         {"E_F_eq",CHEM_EQUILIBRIUM_CONSTRAINT,5,e,f,1,1,0,.012,-2000,0,0}};
    ChemState *s=chem_state_create_constrained(&m,eq,2,1); if (!s) return 0;
    int ok=chem_find_constraint(s,"B_C_eq")==0 && chem_find_constraint(s,"E_F_eq")==1;
    for (size_t i=0;i<6;++i) ok &= chem_is_algebraic(s,i)==(i==2 || i==5);
    const double y[]={1,.1,.0034,.5,.05,.0006};
    for (size_t i=0;i<6;++i) chem_species_values(s,i)[0]=y[i];
    chem_temperatures(s)[0]=310;
    double k1=expected_constant(.034,-2400,310), k2=expected_constant(.012,-2000,310);
    ok &= chem_update_rates(s)==CHEM_OK && fabs(chem_constraint_constants(s,0)[0]-k1)<=1e-10 &&
          fabs(chem_constraint_constants(s,1)[0]-k2)<=fabs(k2)*100*DBL_EPSILON;
    chem_state_free(s); return ok;
}

static int solve_linear_equilibria(int count)
{
    const int reactants[]={0}; const ChemProduct products[]={{1,1}};
    const ChemReaction rx[]={{0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=1,.d=300}}};
    ChemMechanism m=mechanism(count==1 ? 3 : 4,rx,1,reactants,products);
    ChemConstraint eq[]={{"B_C_eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,b,c,1,1,0,count==1 ? 2 : 3,-2400,0,0},
                         {"B_D_eq",CHEM_EQUILIBRIUM_CONSTRAINT,3,b,d,1,1,0,5,-2400,0,0}};
    ChemState *s=chem_state_create_constrained(&m,eq,(size_t)count,1); if (!s) return 0;
    chem_species_values(s,0)[0]=1; chem_temperatures(s)[0]=count==1 ? 270 : 300;
    int ok=chem_update_rates(s)==CHEM_OK;
    for (int i=0;i<count;++i) {
        double expected=expected_constant(eq[i].equilibrium_constant,-2400,chem_temperatures(s)[0]);
        ok &= fabs(chem_constraint_constants(s,(size_t)i)[0]-expected)<=1e-10;
    }
    ChemRosenbrock p=chem_rosenbrock_three_stage();
    for (int step=0;step<100 && ok;++step) {
        ChemStats stats;
        ok &= chem_update_rates(s)==CHEM_OK && chem_solve(s,&p,.001,&stats)==CHEM_OK;
        for (int i=0;i<count && ok;++i)
            ok &= fabs(chem_constraint_constants(s,(size_t)i)[0]*chem_species_values(s,1)[0]-
                       chem_species_values(s,(size_t)(2+i))[0])<=1e-6;
    }
    ok &= fabs(chem_species_values(s,0)[0]+chem_species_values(s,1)[0]-1)<=.01;
    chem_state_free(s); return ok;
}

static int nonunit(void)
{
    const int reactants[]={2}; const ChemProduct products[]={{0,1}};
    const ChemReaction rx[]={{0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=.5,.d=300}}};
    ChemMechanism m=mechanism(3,rx,1,reactants,products);
    ChemConstraint eq={"A2_B_eq",CHEM_EQUILIBRIUM_CONSTRAINT,1,a2,b,1,1,0,10,-2400,0,0};
    ChemState *s=chem_state_create_constrained(&m,&eq,1,1); if (!s) return 0;
    chem_species_values(s,0)[0]=.1; chem_species_values(s,1)[0]=.1; chem_species_values(s,2)[0]=1;
    chem_temperatures(s)[0]=298;
    int ok=chem_update_rates(s)==CHEM_OK && fabs(chem_constraint_constants(s,0)[0]-expected_constant(10,-2400,298))<=1e-10;
    ChemRosenbrock p=chem_rosenbrock_three_stage();
    for (int step=0;step<50 && ok;++step) {
        ChemStats stats;
        ok &= chem_update_rates(s)==CHEM_OK && chem_solve(s,&p,.001,&stats)==CHEM_OK;
        double av=chem_species_values(s,0)[0], bv=chem_species_values(s,1)[0];
        ok &= fabs(chem_constraint_constants(s,0)[0]*av*av-bv)<=1e-5;
    }
    chem_state_free(s); return ok;
}

static int reordered(void)
{
    const ChemMechanism *mechanisms[]={&constrained_original,&constrained_ordered};
    const ChemConstraint *constraints[]={constrained_original_constraints,constrained_ordered_constraints};
    const size_t counts[]={constrained_original_constraint_count,constrained_ordered_constraint_count};
    ChemState *states[2]={NULL,NULL}; int ai[2],bi[2],ci[2],ok=1;
    for (size_t i=0;i<2;++i) {
        ai[i]=chem_find_species(mechanisms[i],"A"); bi[i]=chem_find_species(mechanisms[i],"B"); ci[i]=chem_find_species(mechanisms[i],"C");
        states[i]=chem_state_create_constrained(mechanisms[i],constraints[i],counts[i],1);
        if (!states[i] || ai[i]<0 || bi[i]<0 || ci[i]<0 || counts[i]!=1) { ok=0; break; }
        chem_species_values(states[i],(size_t)ai[i])[0]=1;
        ok &= chem_update_rates(states[i])==CHEM_OK &&
              fabs(chem_constraint_constants(states[i],0)[0]-expected_constant(2,-2400,400))<=1e-10;
    }
    if (ok) ok &= ai[0]!=ai[1] || bi[0]!=bi[1] || ci[0]!=ci[1];
    ChemRosenbrock p=chem_rosenbrock_three_stage();
    for (int step=0;step<100 && ok;++step) {
        for (size_t i=0;i<2;++i) {
            ChemStats stats;
            ok &= chem_update_rates(states[i])==CHEM_OK && chem_solve(states[i],&p,.001,&stats)==CHEM_OK;
            double bv=chem_species_values(states[i],(size_t)bi[i])[0], cv=chem_species_values(states[i],(size_t)ci[i])[0];
            ok &= fabs(chem_constraint_constants(states[i],0)[0]*bv-cv)<=1e-6;
        }
    }
    if (ok) for (size_t i=0;i<2;++i)
        ok &= fabs(chem_species_values(states[i],(size_t)ai[i])[0]+chem_species_values(states[i],(size_t)bi[i])[0]-1)<=.01;
    if (ok) for (int species_index=0;species_index<3;++species_index) {
        int indices[]={species_index==0 ? ai[0] : species_index==1 ? bi[0] : ci[0],
                       species_index==0 ? ai[1] : species_index==1 ? bi[1] : ci[1]};
        ok &= fabs(chem_species_values(states[0],(size_t)indices[0])[0]-chem_species_values(states[1],(size_t)indices[1])[0])<1e-11;
    }
    chem_state_free(states[0]); chem_state_free(states[1]);
    return ok;
}

int main(void)
{
    report("SetConstraintsAPIWorks",single_api());
    report("SetConstraintsAPIMultipleConstraints",multiple_api());
    report("DAESolveWithConstraint",solve_linear_equilibria(1));
    report("DAESolveWithConstraintAndReorderState",reordered());
    report("DAESolveWithTwoCoupledConstraints",solve_linear_equilibria(2));
    report("DAESolveWithNonUnitStoichiometry",nonunit());
    return failures ? 1 : 0;
}
