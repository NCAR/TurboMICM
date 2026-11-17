#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures;
static void report(const char *prefix, const char *name, int ok)
{
    printf("%s:%s %s %s\n",prefix,name,ok ? "PASS" : "FAIL",ok ? "reference quantities matched" : "reference quantity differed");
    failures += !ok;
}
static int near(double a, double b, double tol) { return isfinite(a) && fabs(a-b)<=tol; }
static const ChemSpecies sp[]={{"A",0,1e-3,0,0},{"B",0,1e-3,0,0},{"C",0,1e-3,0,0}};

static int evaluate(const ChemConstraint *c, size_t n, const double *y, double *f, double *jac)
{
    ChemMechanism m={n,0,0,sp,NULL,NULL,NULL,NULL,NULL,298.15,101325};
    ChemState *s=chem_state_create_constrained(&m,c,1,1);
    if (!s) return 0;
    for (size_t i=0;i<n;++i) chem_species_values(s,i)[0]=y[i];
    int ok=chem_update_rates(s)==CHEM_OK && chem_evaluate(s,0,f,jac)==CHEM_OK;
    chem_state_free(s);
    return ok;
}

static int finite_difference(const ChemConstraint *c, size_t n, const double *original)
{
    double y[3], f[3], jac[9], plus[3], minus[3], scratch[9];
    memcpy(y,original,n*sizeof(*y));
    if (!evaluate(c,n,y,f,jac)) return 0;
    for (size_t col=0;col<n;++col) {
        double h=1e-8*fmax(1,fabs(y[col])), saved=y[col];
        y[col]=saved+h; if (!evaluate(c,n,y,plus,scratch)) return 0;
        y[col]=saved-h; if (!evaluate(c,n,y,minus,scratch)) return 0;
        y[col]=saved;
        for (size_t row=0;row<n;++row) {
            double fd=(plus[row]-minus[row])/(2*h), analytic=-jac[row*n+col];
            if (!near(fd,analytic,1e-7+1e-7*fmax(fabs(fd),fabs(analytic)))) return 0;
        }
    }
    return 1;
}

static void linear_cases(void)
{
    ChemTerm ab[]={{0,1},{1,1}}, abc[]={{0,1},{1,1},{2,1}}, weighted[]={{0,2},{1,3},{2,-1}};
    ChemTerm balance[]={{0,1},{1,-1}}, fractional[]={{0,.5},{1,1.5}}, pair_weighted[]={{0,2},{1,3}};
    ChemConstraint c={"A_B_conservation",CHEM_LINEAR_CONSTRAINT,1,ab,NULL,2,0,1,0,0,0,0};
    ChemMechanism m={3,0,0,sp,NULL,NULL,NULL,NULL,NULL,298.15,101325};
    report("linear","LinearConstraint.Construction",!strcmp(c.name,"A_B_conservation") && c.constant==1 &&
           c.reactant_count==2 && c.reactants[0].species==0 && c.reactants[1].species==1 &&
           chem_validate_constraints(&m,&c,1)==CHEM_OK);
    c.name="A_B_C_conservation"; c.algebraic_species=2; c.reactants=abc; c.reactant_count=3; c.constant=10;
    report("linear","LinearConstraint.AlgebraicSpecies",c.algebraic_species==2 && chem_validate_constraints(&m,&c,1)==CHEM_OK);
    c.name="weighted_sum"; c.reactants=weighted; c.constant=5;
    report("linear","LinearConstraint.WeightedTerms",!strcmp(c.name,"weighted_sum") && c.reactant_count==3 && c.constant==5 &&
           c.reactants[0].coefficient==2 && c.reactants[1].coefficient==3 && c.reactants[2].coefficient==-1);
    c.name="A_equals_B"; c.algebraic_species=1; c.reactants=balance; c.reactant_count=2; c.constant=0;
    report("linear","LinearConstraint.ZeroConstant",!strcmp(c.name,"A_equals_B") && c.constant==0 && c.reactant_count==2);
    c.name="fractional_conservation"; c.reactants=fractional; c.constant=2;
    report("linear","LinearConstraint.FractionalCoefficients",!strcmp(c.name,"fractional_conservation") && c.constant==2 &&
           c.reactants[0].coefficient==.5 && c.reactants[1].coefficient==1.5);
    double f[3], j[9], y[]={.3,.7,0};
    c.name="A_B_conservation"; c.reactants=ab; c.constant=1;
    int ok=evaluate(&c,2,y,f,j) && near(f[1],0,1e-10);
    y[0]=.5; y[1]=.6; ok &= evaluate(&c,2,y,f,j) && near(f[1],.1,1e-10);
    report("linear","LinearConstraint.ResidualComputationThroughConstraintSet",ok);
    y[0]=.3; y[1]=.7;
    report("linear","LinearConstraint.JacobianComputationThroughConstraintSet",evaluate(&c,2,y,f,j) &&
           near(j[2],-1,1e-10) && near(j[3],-1,1e-10));
    c.algebraic_species=2; c.reactants=weighted; c.reactant_count=3; c.constant=5;
    y[0]=1; y[1]=2; y[2]=3;
    report("linear","LinearConstraint.WeightedSumResidualAndJacobian",evaluate(&c,3,y,f,j) && near(f[2],0,1e-10) &&
           near(j[6],-2,1e-10) && near(j[7],-3,1e-10) && near(j[8],1,1e-10));
    c.reactants=abc; c.constant=10; y[0]=2; y[1]=3; y[2]=5;
    ok=evaluate(&c,3,y,f,j) && near(f[2],0,1e-10);
    y[2]=6; ok &= evaluate(&c,3,y,f,j) && near(f[2],1,1e-10);
    report("linear","LinearConstraint.ThreeSpeciesConservationResidual",ok);
    c.algebraic_species=1; c.reactants=balance; c.reactant_count=2; c.constant=0; y[0]=.5; y[1]=.5;
    ok=evaluate(&c,2,y,f,j) && near(f[1],0,1e-10);
    y[0]=.6; ok &= evaluate(&c,2,y,f,j) && near(f[1],.1,1e-10);
    report("linear","LinearConstraint.ZeroConstantResidual",ok);
    c.reactants=fractional; c.constant=2; y[0]=y[1]=1;
    report("linear","LinearConstraint.FractionalCoefficientsResidualAndJacobian",evaluate(&c,2,y,f,j) && near(f[1],0,1e-10) &&
           near(j[2],-.5,1e-10) && near(j[3],-1.5,1e-10));
    c.reactants=pair_weighted; c.constant=5; y[0]=.1; y[1]=.9;
    ok=evaluate(&c,2,y,f,j) && near(j[2],-2,1e-10) && near(j[3],-3,1e-10);
    y[0]=10; y[1]=20; ok &= evaluate(&c,2,y,f,j) && near(j[2],-2,1e-10) && near(j[3],-3,1e-10);
    report("linear","LinearConstraint.JacobianIndependentOfConcentrations",ok);
    c.reactants=ab; c.constant=1; double y1[]={.3,.7}, y2[]={.8,.2};
    report("linear","LinearConstraint.FiniteDifferenceJacobianSimpleConservation",finite_difference(&c,2,y1) && finite_difference(&c,2,y2));
    c.algebraic_species=2; c.reactants=weighted; c.reactant_count=3; c.constant=5; double y3[]={1.5,.8,2};
    report("linear","LinearConstraint.FiniteDifferenceJacobianWeightedSum",finite_difference(&c,3,y3));
}

static void equilibrium_cases(void)
{
    ChemTerm ab[]={{0,1},{1,1}}, product[]={{2,1}}, single[]={{0,1}}, b[]={{1,1}}, twice[]={{0,2}}, bc[]={{1,1},{2,1}};
    ChemConstraint c={"A_B_equilibrium",CHEM_EQUILIBRIUM_CONSTRAINT,2,ab,product,2,1,0,1000,-2400,0,0};
    ChemMechanism m={3,0,0,sp,NULL,NULL,NULL,NULL,NULL,298.15,101325};
    report("equilibrium","EquilibriumConstraint.Construction",!strcmp(c.name,"A_B_equilibrium") && c.reactant_count==2 &&
           c.product_count==1 && c.reactants[0].species==0 && c.reactants[1].species==1 && c.products[0].species==2 &&
           chem_validate_constraints(&m,&c,1)==CHEM_OK);
    report("equilibrium","EquilibriumConstraint.AlgebraicSpecies",c.algebraic_species==2);
    c.name="A_B_simple"; c.algebraic_species=1; c.reactants=single; c.reactant_count=1; c.products=b; c.equilibrium_constant=10;
    report("equilibrium","EquilibriumConstraint.SingleReactantSingleProduct",!strcmp(c.name,"A_B_simple") && c.algebraic_species==1 &&
           c.reactant_count==1 && c.product_count==1 && c.reactants[0].species==0 && c.products[0].species==1);
    c.name="dissociation"; c.reactants=twice; c.products=bc; c.product_count=2; c.equilibrium_constant=100;
    report("equilibrium","EquilibriumConstraint.MultipleReactantsAndProducts",!strcmp(c.name,"dissociation") && c.algebraic_species==1 &&
           c.reactant_count==1 && c.product_count==2 && c.reactants[0].coefficient==2 && c.products[1].species==2);
    c.reactants=single; c.products=b; c.product_count=1;
    c.equilibrium_constant=-1; int ok=chem_validate_constraints(&m,&c,1)==CHEM_INVALID;
    c.equilibrium_constant=0; ok &= chem_validate_constraints(&m,&c,1)==CHEM_INVALID;
    report("equilibrium","EquilibriumConstraint.InvalidEquilibriumConstant",ok);
    c.equilibrium_constant=1; c.reactant_count=0;
    report("equilibrium","EquilibriumConstraint.EmptyReactantsThrows",chem_validate_constraints(&m,&c,1)==CHEM_INVALID);
    c.reactant_count=1; c.product_count=0;
    report("equilibrium","EquilibriumConstraint.EmptyProductsThrows",chem_validate_constraints(&m,&c,1)==CHEM_INVALID);
    c.product_count=1;
    single[0].coefficient=0; ok=chem_validate_constraints(&m,&c,1)==CHEM_INVALID;
    single[0].coefficient=-1; ok &= chem_validate_constraints(&m,&c,1)==CHEM_INVALID;
    single[0].coefficient=1; b[0].coefficient=0; ok &= chem_validate_constraints(&m,&c,1)==CHEM_INVALID;
    b[0].coefficient=-2; ok &= chem_validate_constraints(&m,&c,1)==CHEM_INVALID;
    report("equilibrium","EquilibriumConstraint.InvalidStoichiometryThrows",ok);
    c.algebraic_species=2; c.reactants=ab; c.reactant_count=2; c.products=product; c.product_count=1; c.equilibrium_constant=1000; c.enthalpy=0;
    double f[3], j[9], y[]={.001,.001,.001};
    ok=evaluate(&c,3,y,f,j) && near(f[2],0,1e-10);
    y[0]=.02; y[1]=.03; y[2]=.05; ok &= evaluate(&c,3,y,f,j) && near(f[2],.55,1e-10);
    report("equilibrium","EquilibriumConstraint.ResidualComputationThroughConstraintSet",ok);
    y[0]=.01; y[1]=.02; y[2]=.05;
    report("equilibrium","EquilibriumConstraint.JacobianComputationThroughConstraintSet",evaluate(&c,3,y,f,j) &&
           near(j[6],-20,1e-8) && near(j[7],-10,1e-8) && near(j[8],1,1e-8));
    c.algebraic_species=1; c.reactants=twice; c.reactant_count=1; c.products=bc; c.product_count=2; c.equilibrium_constant=100;
    y[0]=.1; y[1]=.5; y[2]=2;
    report("equilibrium","EquilibriumConstraint.ComplexStoichiometryResidual",evaluate(&c,3,y,f,j) && near(f[1],0,1e-10));
    c.algebraic_species=2; c.reactants=ab; c.reactant_count=2; c.products=product; c.product_count=1; c.equilibrium_constant=1000;
    const double v1[]={.02,.03,.05}, v2[]={.1,.2,5};
    report("equilibrium","EquilibriumConstraint.FiniteDifferenceJacobianSimple",finite_difference(&c,3,v1) && finite_difference(&c,3,v2));
    c.algebraic_species=1; c.reactants=twice; c.reactant_count=1; c.products=bc; c.product_count=2; c.equilibrium_constant=100;
    const double v3[]={.15,.8,1.5};
    report("equilibrium","EquilibriumConstraint.FiniteDifferenceJacobianComplexStoichiometry",finite_difference(&c,3,v3));
}

int main(void)
{
    linear_cases(); equilibrium_cases();
    return failures ? 1 : 0;
}
