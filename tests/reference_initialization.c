#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures;
static void report(const char *name, int ok)
{
    printf("constraint_initialization:ConstraintInitialization.%s %s %s\n",name,ok ? "PASS" : "FAIL",
           ok ? "reference quantities matched" : "reference quantity differed");
    failures += !ok;
}
static const ChemSpecies species[]={{"A",1,1e-3,0,0},{"B",.5,1e-3,0,0},{"C",1,1e-3,0,0}};
static const int reactants[]={0}; static const ChemProduct products[]={{1,1}};
static const ChemReaction reactions[]={{0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=.5,.d=300}}};
static const ChemMechanism mechanism={3,1,0,species,reactions,reactants,products,NULL,NULL,298.15,101325};
static const ChemTerm from_b[]={{1,1}}, to_c[]={{2,1}}, to_c_squared[]={{2,2}};
static const ChemConstraint eq={"B_C_eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,from_b,to_c,1,1,0,2,0,0,0};
static const ChemConstraint squared={"squared",CHEM_EQUILIBRIUM_CONSTRAINT,2,from_b,to_c_squared,1,1,0,2,0,0,0};

static ChemState *state(const ChemConstraint *constraint, size_t cells)
{
    ChemState *s=chem_state_create_constrained(&mechanism,constraint,1,cells);
    if (s && chem_update_rates(s)!=CHEM_OK) { chem_state_free(s); return NULL; }
    return s;
}

static int integrate_guess(double guess, int require_one)
{
    ChemState *s=state(&eq,1); if (!s) return 0;
    chem_species_values(s,2)[0]=guess;
    ChemRosenbrock p=chem_rosenbrock_three_stage(); ChemStats stats;
    int ok=chem_solve(s,&p,.001,&stats)==CHEM_OK &&
           (require_one ? stats.constraint_init_iterations<=1 : stats.constraint_init_iterations>0) &&
           fabs(2*chem_species_values(s,1)[0]-chem_species_values(s,2)[0])<=1e-6;
    chem_state_free(s); return ok;
}

static int scaled_solve(void)
{
    for (int k=0;k<3;++k) {
        double scale=k==0 ? 1 : k==1 ? 1e-12 : 1e12;
        ChemTerm terms[]={{0,scale},{1,scale},{2,scale}};
        ChemConstraint c={"mass",CHEM_LINEAR_CONSTRAINT,2,terms,NULL,3,0,1.3*scale,0,0,0,0};
        ChemState *s=state(&c,1); if (!s) return 0;
        chem_species_values(s,0)[0]=.7; chem_species_values(s,1)[0]=.15; chem_species_values(s,2)[0]=5;
        ChemRosenbrock p=chem_rosenbrock_four_stage_dae(); ChemStats stats;
        int ok=chem_solve(s,&p,.001,&stats)==CHEM_OK && fabs(chem_species_values(s,2)[0]-.45)<=4.5e-6;
        chem_state_free(s); if (!ok) return 0;
    }
    return 1;
}

static int pure_ode(void)
{
    ChemState *s=chem_state_create(&mechanism,1); if (!s) return 0;
    ChemRosenbrock p=chem_rosenbrock_three_stage(); ChemStats stats;
    int ok=chem_update_rates(s)==CHEM_OK && chem_solve(s,&p,.01,&stats)==CHEM_OK && stats.constraint_init_iterations==0;
    chem_state_free(s); return ok;
}

static int multiple_cells(void)
{
    ChemState *s=state(&eq,3); if (!s) return 0;
    const double a[]={1,2,.5}, b[]={.5,.3,.1}, c[]={1,5,100};
    for (int i=0;i<3;++i) {
        chem_species_values(s,0)[i]=a[i]; chem_species_values(s,1)[i]=b[i]; chem_species_values(s,2)[i]=c[i];
    }
    ChemRosenbrock p=chem_rosenbrock_three_stage(); ChemStats stats;
    int ok=chem_solve(s,&p,.001,&stats)==CHEM_OK;
    for (int i=0;i<3;++i) ok &= fabs(2*chem_species_values(s,1)[i]-chem_species_values(s,2)[i])<=1e-5;
    chem_state_free(s); return ok;
}

static int subsequent(void)
{
    ChemState *s=state(&eq,1); if (!s) return 0;
    ChemRosenbrock p=chem_rosenbrock_three_stage(); ChemStats stats;
    int ok=chem_solve(s,&p,.001,&stats)==CHEM_OK;
    chem_species_values(s,2)[0]=999;
    ok &= chem_update_rates(s)==CHEM_OK && chem_solve(s,&p,.001,&stats)==CHEM_OK &&
          fabs(2*chem_species_values(s,1)[0]-chem_species_values(s,2)[0])<=1e-6;
    chem_state_free(s); return ok;
}

static int failed_rollback(void)
{
    ChemState *s=state(&squared,1); if (!s) return 0;
    chem_species_values(s,2)[0]=100;
    ChemRosenbrock p=chem_rosenbrock_three_stage(); p.constraint_max_iterations=1; ChemStats stats;
    int ok=chem_initialize_constraints(s,&p,&stats)==CHEM_CONSTRAINT_INITIALIZATION_FAILED &&
           chem_species_values(s,1)[0]==.5 && chem_species_values(s,2)[0]==100;
    chem_state_free(s); return ok;
}

static int successful_rollback(void)
{
    ChemState *s=state(&squared,1); if (!s) return 0;
    ChemRosenbrock p=chem_rosenbrock_three_stage(); ChemStats stats;
    chem_species_values(s,2)[0]=1e-6;
    int ok=chem_initialize_constraints(s,&p,&stats)==CHEM_OK && fabs(chem_species_values(s,2)[0]-1)<=1e-3;
    for (int i=0;i<2;++i) {
        p.constraint_max_backtracks=i ? 2 : 0; chem_species_values(s,2)[0]=1e-6;
        ok &= chem_initialize_constraints(s,&p,&stats)==CHEM_CONSTRAINT_INITIALIZATION_FAILED &&
              chem_species_values(s,1)[0]==.5 && chem_species_values(s,2)[0]==1e-6;
    }
    chem_state_free(s); return ok;
}

static int scaled_correction(void)
{
    for (int k=0;k<3;++k) {
        double scale=k==0 ? 1e-12 : k==1 ? 1 : 1e12;
        ChemTerm terms[]={{0,scale},{2,-scale}};
        ChemConstraint c={"scaled",CHEM_LINEAR_CONSTRAINT,2,terms,NULL,2,0,0,0,0,0,0};
        ChemState *s=state(&c,1); if (!s) return 0;
        chem_species_values(s,2)[0]=4;
        for (int i=0;i<3;++i) chem_absolute_tolerances(s)[i]=1e-12;
        ChemRosenbrock p=chem_rosenbrock_three_stage(); p.relative_tolerance=1e-8; ChemStats stats;
        int ok=chem_initialize_constraints(s,&p,&stats)==CHEM_OK && fabs(chem_species_values(s,2)[0]-1)<=1e-12;
        chem_state_free(s); if (!ok) return 0;
    }
    return 1;
}

static int weighted_tolerances(void)
{
    ChemTerm terms[]={{0,1},{2,-1}};
    ChemConstraint c={"scaled",CHEM_LINEAR_CONSTRAINT,2,terms,NULL,2,0,0,0,0,0,0};
    ChemRosenbrock p=chem_rosenbrock_three_stage();
    int ok=p.constraint_tolerance==.1;
    for (int tight=0;tight<2;++tight) {
        ChemState *s=state(&c,1); if (!s) return 0;
        double initial=1+5e-8; chem_species_values(s,2)[0]=initial;
        for (int i=0;i<3;++i) chem_absolute_tolerances(s)[i]=tight ? 1e-12 : 1e-8;
        p.relative_tolerance=tight ? 1e-8 : 1e-6; ChemStats stats;
        ok &= chem_initialize_constraints(s,&p,&stats)==CHEM_OK;
        if (tight) ok &= fabs(chem_species_values(s,2)[0]-1)<=1e-12 && stats.constraint_init_iterations==2;
        else ok &= chem_species_values(s,2)[0]==initial && stats.constraint_init_iterations==1;
        chem_state_free(s);
    }
    return ok;
}

int main(void)
{
    report("ConsistentICsUnchanged",integrate_guess(1,1));
    report("MildlyInconsistentICsCorrected",integrate_guess(5,0));
    report("SeverelyInconsistentICsConverge",integrate_guess(1000,0));
    report("ScaledConstraintRowsSolveIdenticallyThroughSolve",scaled_solve());
    report("PureODESystemUnaffected",pure_ode());
    report("MultiCellSystems",multiple_cells());
    report("SubsequentSolveCallsReinitialize",subsequent());
    report("SolverStateToStringNewValue",!strcmp(chem_error_string(CHEM_CONSTRAINT_INITIALIZATION_FAILED),"Constraint Initialization Failed"));
    report("FailedInitializationRestoresCallerState",failed_rollback());
    report("SuccessfulInitializationKeepsTheCorrection",successful_rollback());
    report("WeightedCorrectionIsInvariantToConstraintRowScaling",scaled_correction());
    report("WeightedCorrectionUsesStateTolerances",weighted_tolerances());
    return failures ? 1 : 0;
}
