#include "chem.h"

#include <math.h>
#include <stdio.h>

static int failures;
static void report(const char *prefix, const char *name, int ok)
{
    printf("%s:%s %s %s\n", prefix, name, ok ? "PASS" : "FAIL",
           ok ? "reference quantities matched" : "reference quantity differed");
    failures += !ok;
}

static int conservation(ChemRosenbrock parameters)
{
    const ChemSpecies species[]={{"A",.9e-6,1e-12,0,0},{"B",0,1e-12,0,0},{"C",.1e-6,1e-12,0,0}};
    const int reactants[]={0}; const ChemProduct products[]={{1,1}};
    const ChemReaction reaction={0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=1e4,.d=300}};
    const ChemTerm terms[]={{0,1},{1,1},{2,1}};
    const ChemConstraint constraint={"mass",CHEM_LINEAR_CONSTRAINT,2,terms,NULL,3,0,1e-6,0,0,0,0};
    const ChemMechanism mechanism={3,1,0,species,&reaction,reactants,products,NULL,NULL,298,101325};
    ChemState *s=chem_state_create_constrained(&mechanism,&constraint,1,1);
    if (!s) return 0;
    parameters.relative_tolerance=1e-6;
    int ok=chem_update_rates(s)==CHEM_OK;
    double advanced=0;
    while (ok && advanced<30) {
        ChemStats stats;
        ok=chem_solve(s,&parameters,30-advanced,&stats)==CHEM_OK && stats.final_time>0;
        advanced+=stats.final_time;
        double a=chem_species_values(s,0)[0], b=chem_species_values(s,1)[0], c=chem_species_values(s,2)[0];
        ok &= fabs(a+b+c-1e-6)<=1e-12 && c>=-1e-18;
    }
    ok &= chem_species_values(s,0)[0]<1e-12 && fabs(chem_species_values(s,1)[0]-.9e-6)<=1e-10 &&
          fabs(chem_species_values(s,2)[0]-.1e-6)<=1e-10;
    chem_state_free(s); return ok;
}

static int equilibrium_conservation(void)
{
    const ChemSpecies species[]={{"A_gas",1e-6/11,1e-8,0,0},{"A_aq",1e-5/11,1e-8,0,0},{"P",0,1e-12,0,0}};
    const int reactants[]={1}; const ChemProduct products[]={{2,1}};
    const ChemReaction reaction={0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=1e3,.d=300}};
    const ChemTerm gas[]={{0,1}}, aq[]={{1,1}}, mass[]={{1,1},{2,1},{0,1}};
    const ChemConstraint constraints[]={
        {"gas_aq",CHEM_EQUILIBRIUM_CONSTRAINT,1,gas,aq,1,1,0,10,0,0,0},
        {"mass",CHEM_LINEAR_CONSTRAINT,0,mass,NULL,3,0,1e-6,0,0,0,0}};
    const ChemMechanism mechanism={3,1,0,species,&reaction,reactants,products,NULL,NULL,298,101325};
    ChemState *s=chem_state_create_constrained(&mechanism,constraints,2,1);
    if (!s) return 0;
    ChemRosenbrock p=chem_rosenbrock_four_stage_dae();
    int ok=chem_update_rates(s)==CHEM_OK;
    double advanced=0;
    while (ok && advanced<30) {
        ChemStats stats;
        ok=chem_solve(s,&p,30-advanced,&stats)==CHEM_OK && stats.final_time>0;
        advanced+=stats.final_time;
        double g=chem_species_values(s,0)[0], a=chem_species_values(s,1)[0], product=chem_species_values(s,2)[0];
        ok &= fabs(g+a+product-1e-6)<=1e-12 && g>=-1e-18 && a>=-1e-18;
    }
    ok &= fabs(chem_species_values(s,2)[0]-1e-6)<=1e-8 && chem_species_values(s,0)[0]>=-1e-18;
    chem_state_free(s); return ok;
}

typedef struct { double gas, product, minimum_gas; unsigned long accepted; int ok; } CascadeResult;
static CascadeResult cascade(double balance_atol, double rate, double k1, double k2)
{
    double initial=1e-6/(1+k1+k1*k2);
    const ChemSpecies species[]={{"A_gas",initial,balance_atol,0,0},{"A_aq",k1*initial,1e-8,0,0},
                                {"B_aq",k1*k2*initial,1e-8,0,0},{"P",0,1e-12,0,0}};
    const int reactants[]={2}; const ChemProduct products[]={{3,1}};
    const ChemReaction reaction={0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=rate,.d=300}};
    const ChemTerm gas[]={{0,1}}, aq[]={{1,1}}, bq[]={{2,1}}, mass[]={{1,1},{2,1},{3,1},{0,1}};
    const ChemConstraint constraints[]={
        {"eq1",CHEM_EQUILIBRIUM_CONSTRAINT,1,gas,aq,1,1,0,k1,0,0,0},
        {"eq2",CHEM_EQUILIBRIUM_CONSTRAINT,2,aq,bq,1,1,0,k2,0,0,0},
        {"mass",CHEM_LINEAR_CONSTRAINT,0,mass,NULL,4,0,1e-6,0,0,0,0}};
    const ChemMechanism mechanism={4,1,0,species,&reaction,reactants,products,NULL,NULL,298,101325};
    CascadeResult result={0,0,1e99,0,0};
    ChemState *s=chem_state_create_constrained(&mechanism,constraints,3,1);
    if (!s) return result;
    ChemRosenbrock p=chem_rosenbrock_four_stage_dae();
    result.ok=chem_update_rates(s)==CHEM_OK;
    double advanced=0;
    while (result.ok && advanced<30) {
        ChemStats stats;
        result.ok=chem_solve(s,&p,30-advanced,&stats)==CHEM_OK && stats.final_time>0;
        advanced+=stats.final_time; result.accepted+=stats.accepted;
        result.minimum_gas=fmin(result.minimum_gas,chem_species_values(s,0)[0]);
    }
    result.gas=chem_species_values(s,0)[0]; result.product=chem_species_values(s,3)[0];
    chem_state_free(s); return result;
}

static int terminator_robertson(void)
{
    const ChemSpecies species[]={{"Cl2",1.8e-10,1e-20,0,0},{"Cl",1.2e-6,1e-20,0,0},
                                 {"A",1,1e-3,0,0},{"B",0,1e-3,0,0},{"C",0,1e-3,0,0}};
    const int reactants[]={0,1,1,2,3,3,3,4};
    const ChemProduct products[]={{1,2},{0,1},{3,1},{3,1},{4,1},{2,1},{4,1}};
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_USER,.parameter=0,.a=1}},
        {1,2,1,1,0,{.kind=CHEM_ARRHENIUS,.a=1,.d=300}},
        {3,1,2,1,0,{.kind=CHEM_USER,.parameter=1,.a=1}},
        {4,2,3,2,0,{.kind=CHEM_USER,.parameter=2,.a=1}},
        {6,2,5,2,0,{.kind=CHEM_USER,.parameter=3,.a=1}}};
    const char *names[]={"terminator_k1","robertson_r1","robertson_r2","robertson_r3"};
    const double rates[]={0,.04,3e7,1e4};
    const ChemTerm mass[]={{2,1},{3,1},{4,1}};
    const ChemConstraint constraint={"mass",CHEM_LINEAR_CONSTRAINT,4,mass,NULL,3,0,1,0,0,0,0};
    const ChemMechanism mechanism={5,5,4,species,reactions,reactants,products,names,rates,298,101300};
    ChemState *s=chem_state_create_constrained(&mechanism,&constraint,1,1);
    if (!s) return 0;
    ChemRosenbrock p=chem_rosenbrock_four_stage_dae(); p.relative_tolerance=1e-8;
    int ok=chem_set_air_density(s,0,42)==CHEM_OK && chem_update_rates(s)==CHEM_OK;
    const double chlorine=2*1.8e-10+1.2e-6;
    double dt=1;
    for (int step=0;step<12 && ok;++step,dt*=10) {
        double advanced=0;
        while (advanced<dt && ok) {
            ChemStats stats;
            int status=chem_solve(s,&p,dt-advanced,&stats);
            ok=(status==CHEM_OK || status==CHEM_TOO_MANY_STEPS) && stats.final_time>0 && advanced+stats.final_time>advanced;
            advanced+=stats.final_time;
        }
        double cl2=chem_species_values(s,0)[0], cl=chem_species_values(s,1)[0];
        double a=chem_species_values(s,2)[0], b=chem_species_values(s,3)[0], c=chem_species_values(s,4)[0];
        ok &= fabs(a+b+c-1)<=1e-10 && fabs(2*cl2+cl-chlorine)<=1e-6*chlorine;
        for (size_t i=0;i<5;++i) { double v=chem_species_values(s,i)[0]; ok &= isfinite(v) && v>=0; }
        if (step==0) ok &= a<.99 && b>0 && c>1e-3;
    }
    ok &= chem_species_values(s,4)[0]>.95 && chem_species_values(s,2)[0]<.05 &&
          chem_species_values(s,0)[0]>.4*chlorine && chem_species_values(s,1)[0]<.1*chlorine;
    chem_state_free(s); return ok;
}

int main(void)
{
    report("dae_constraint_overshoot","DAEConstraintOvershoot.AlgebraicVariableStaysNonNegative",conservation(chem_rosenbrock_four_stage_dae()));
    report("dae_constraint_overshoot","DAEConstraintOvershoot.EquilibriumPlusConservation",equilibrium_conservation());
    int all=conservation(chem_rosenbrock_three_stage()) && conservation(chem_rosenbrock_four_stage()) &&
            conservation(chem_rosenbrock_four_stage_dae()) && conservation(chem_rosenbrock_six_stage_dae());
    report("dae_constraint_overshoot","DAEConstraintOvershoot.AllRosenbrockOrdersConstrained",all);
    CascadeResult loose=cascade(1e-3,100,2,2), tight=cascade(1e-8,100,2,2);
    report("dae_algebraic_error_insensitivity","DAEAlgebraicError.ErrorSensitiveToBalanceAtol",
           loose.ok && tight.ok && tight.accepted>loose.accepted);
    CascadeResult stiff=cascade(1e-8,1e4,100,50);
    report("dae_algebraic_error_insensitivity","DAEAlgebraicError.AlgebraicVariableDoesNotOvershootDeeply",
           stiff.ok && stiff.minimum_gas>=-1e-8 && fabs(stiff.product-1e-6)<=1e-8 && fabs(stiff.gas)<=1e-10);
    report("linear_constraint","DAESolveWithConstraint.TerminatorAndRobertson",terminator_robertson());
    return failures ? 1 : 0;
}
