#include "chem.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static double rate_at(ChemRateType type, const double *p, size_t n, double temperature)
{
    ChemReaction r = {0};
    r.type = type;
    for (size_t i = 0; i < n; ++i) r.p[i] = p[i];
    ChemMechanism m = {0};
    m.reaction_count = 1; m.reactions = &r;
    double value;
    if (chem_rate_constants(&m, (ChemConditions){temperature, 0, 0, 0}, NULL, &value)) return NAN;
    return value;
}

static double rate(ChemRateType type, const double *p, size_t n)
{
    return rate_at(type, p, n, 301.24);
}

static int reversible(double a, double c, double reverse, double temperature)
{
    double p[] = {a, c, reverse};
    double expected = a*exp(c/temperature)*reverse;
    return fabs(rate_at(CHEM_REVERSIBLE, p, 3, temperature)-expected) <= 1e-13*(expected == 0 ? 1 : fabs(expected));
}

static int user_rate(double value, double scale)
{
    ChemReaction reaction = {0}; reaction.type = CHEM_USER; reaction.parameter = 0; reaction.p[0] = scale;
    ChemMechanism m = {0}; m.reaction_count = 1; m.parameter_count = 1; m.reactions = &reaction;
    double result;
    return !chem_rate_constants(&m, (ChemConditions){301.24, 0, 0, 0}, &value, &result) && result == value*scale;
}

static int taylor(void)
{
    double p[22] = {1, 0, 0, 300, 0, 1, 1};
    if (fabs(rate(CHEM_TAYLOR, p, 22)-1) > 1e-13) return 0;
    p[0] = 1;
    if (fabs(rate(CHEM_TAYLOR, p, 22)-1) > 1e-13) return 0;
    p[0] = 2.2e-10;
    if (fabs(rate(CHEM_TAYLOR, p, 22)-p[0]) > 1e-13*p[0]) return 0;
    p[0] = 3e-11; p[2] = -200; p[5] = 3; p[6] = 12.5; p[7] = 1.3e-2; p[8] = 5.2e-4;
    double expected = p[0]*exp(p[2]/301.24)*(p[6]+p[7]*301.24+p[8]*pow(301.24,2));
    if (fabs(rate(CHEM_TAYLOR, p, 22)-expected) > 1e-13*expected) return 0;
    p[0] = 2.6e-12; p[2] = -350; p[6] = 1; p[7] = 4.3e-1; p[8] = 7.3e-3;
    expected = p[0]*exp(p[2]/301.24)*(p[6]+p[7]*301.24+p[8]*pow(301.24,2));
    return fabs(rate(CHEM_TAYLOR, p, 22)-expected) <= 1e-13*expected;
}

static double rate_density(ChemRateType type, const double *p, size_t n, double t, double air)
{
    ChemReaction reaction = {0}; reaction.type = type;
    for (size_t i = 0; i < n; ++i) reaction.p[i] = p[i];
    ChemMechanism mechanism = {0}; mechanism.reaction_count = 1; mechanism.reactions = &reaction;
    double result;
    if (chem_rate_constants(&mechanism, (ChemConditions){t, 0, air, 1}, NULL, &result)) return NAN;
    return result;
}

static int falloff(ChemRateType type, int mode)
{
    double p[] = {1,0,0,1,0,0,.6,1};
    if (mode == 1) { double q[] = {1.2,2.3,302.3,2.6,-3.1,402.1,.9,1.2}; memcpy(p,q,sizeof p); }
    if (mode == 2) p[0] = 4e-18;
    if (mode == 3) { double q[] = {1.2e-12,167,3,136,5,24,.9,.8}; memcpy(p,q,sizeof p); }
    double t = 301.24, air = 42.2;
    double k0 = p[0]*exp(p[2]/t)*pow(t/300,p[1]);
    double ki = p[3]*exp(p[5]/t)*pow(t/300,p[4]);
    double expected = k0*(type == CHEM_TROE ? air : 1.0)/(1+k0*air/ki)
        * pow(p[6],1.0/(1.0+(1.0/p[7])*pow(log10(k0*air/ki),2)));
    double actual = rate_density(type,p,8,t,air);
    if (mode == 3) return fabs(actual-expected)/fmax(fabs(actual),fabs(expected)) <= 1e-14;
    return actual == expected;
}

static int jpl(int zero)
{
    double t = zero ? 250 : 301.24, air = zero ? 30 : 42.2;
    double p[] = {1.2,2.3,302.3,298,2.6,-3.1,402.1,310,.7,.5,120,290,.9,1.2};
    if (zero) { double q[] = {.003,-2,0,298,1.5,-1,0,298,0,0,0,298,.6,1}; memcpy(p,q,sizeof p); }
    double k0 = p[0]*exp(p[2]/t)*pow(t/p[3],p[1]);
    double ki = p[4]*exp(p[6]/t)*pow(t/p[7],p[5]);
    double kint = p[8]*exp(p[10]/t)*pow(t/p[11],p[9]);
    double ratio = k0*air/ki;
    double expected = k0*air/(1+ratio)*pow(p[12],p[13]/(p[13]+pow(log10(ratio),2)));
    expected += kint*(1-expected/ki);
    return fabs(rate_density(CHEM_TERNARY_JPL19,p,14,t,air)-expected) <= (zero ? 1e-30 : 0);
}

static int branched(int nitrate)
{
    double p[] = {1.2,204.3,1e-3,2}, t = 301.24, air = 42.2;
    double a = 2e-22*exp(2)*2.45e19, b = .43*pow(293.0/298,-8);
    double z = a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)))*(1-1e-3)/1e-3;
    a = 2e-22*exp(2)*(air*6.02214076e23*1e-6); b = .43*pow(t/298,-8);
    double af = a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)));
    double expected = 1.2*exp(-204.3/t)*(nitrate ? af/(z+af) : z/(z+af));
    double actual = rate_density(nitrate ? CHEM_BRANCH_NITRATE : CHEM_BRANCH_ALKOXY,p,4,t,air);
    return fabs(actual-expected) <= 1e-13*expected;
}

static int surface(double probability)
{
    ChemReaction reaction = {0}; reaction.type = CHEM_SURFACE; reaction.parameter = 0;
    reaction.p[0] = 230; reaction.p[1] = .025; reaction.p[2] = probability;
    ChemMechanism m = {0}; m.reaction_count = 1; m.parameter_count = 2; m.reactions = &reaction;
    double custom[] = {1e-7, 2.5e6}, actual;
    double gas = 1.380649e-23 * 6.02214076e23, pi = 3.14159265358979323846;
    double expected = 4*2.5e6*pi*pow(1e-7,2) /
        (1e-7/230 + 4/(sqrt(8*gas*273.65/(pi*.025))*probability));
    return !chem_rate_constants(&m, (ChemConditions){273.65, 0, 0, 0}, custom, &actual)
        && fabs(actual-expected) <= 1e-13*expected;
}

static int missing_surface(int weight)
{
    const char *names[] = {"foo"};
    ChemTerm reactant[] = {{0, 1}};
    ChemReaction reaction = {0}; reaction.type = CHEM_SURFACE; reaction.reactant_count = 1;
    reaction.p[0] = weight ? 230 : 0; reaction.p[1] = weight ? 0 : .025;
    ChemMechanism m = {0}; m.species_count = 1; m.species_names = names;
    m.reaction_count = 1; m.reactions = &reaction; m.reactants = reactant;
    char message[128];
    ChemValidation result = chem_validate_mechanism(&m, message, sizeof message);
    return result == (weight ? CHEM_MISSING_MOLECULAR_WEIGHT : CHEM_MISSING_DIFFUSION)
        && strstr(message, "foo") && strstr(message, weight ? "molecular weight" : "Diffusion coefficient");
}

static int backward_convergence(void)
{
    double residual[12] = {0}, y[12] = {0}, atol[] = {1e-6,1e-6,1e-6};
    double small = 1e-6, rtol = 1e-3;
    if (!chem_backward_euler_converged(4,3,residual,y,small,rtol,atol)) return 0;
    residual[1] = 1e-5;
    if (chem_backward_euler_converged(4,3,residual,y,small,rtol,atol)) return 0;
    small = 1e-4;
    if (!chem_backward_euler_converged(4,3,residual,y,small,rtol,atol)) return 0;
    residual[11] = 1e-3;
    if (chem_backward_euler_converged(4,3,residual,y,small,rtol,atol)) return 0;
    y[11] = 10;
    if (!chem_backward_euler_converged(4,3,residual,y,small,rtol,atol)) return 0;
    residual[11] = .1;
    if (chem_backward_euler_converged(4,3,residual,y,small,rtol,atol)) return 0;
    atol[2] = 1;
    return chem_backward_euler_converged(4,3,residual,y,small,rtol,atol);
}

static int backward_solve(void)
{
    const char *names[] = {"A","B","C"};
    const ChemTerm react[] = {{0,1},{1,1}}, product[] = {{1,1},{2,1}};
    const ChemReaction reactions[] = {
        {CHEM_ARRHENIUS,0,1,0,1,-1,{2.15e-11,0,110,300,0}},
        {CHEM_ARRHENIUS,1,1,1,1,-1,{3.3e-11,0,55,300,0}}
    };
    const ChemMechanism m = {3,2,0,names,NULL,NULL,NULL,reactions,react,product};
    ChemConditions c = {272.5,101253.3,1e6,1};
    ChemBackwardEulerParameters p = chem_backward_euler_defaults();
    p.absolute_tolerance = 1e-6;
    ChemWorkspace *work = chem_workspace_create(&m);
    double rates[2], y[] = {1,0,0};
    ChemStatistics stats;
    int ok = work && !chem_rate_constants(&m,c,NULL,rates) && !chem_backward_euler_solve(&m,work,&p,c,rates,1,y,&stats);
    chem_workspace_free(work);
    return ok;
}

static int linear_case(int which)
{
    ChemTerm terms[] = {{0,1},{1,1},{2,1}};
    ChemConstraint c = {"conservation", CHEM_LINEAR_CONSTRAINT, 1, terms, NULL, 2, 0, 1, 0, 0, 0};
    double y[] = {.3,.7,0}, forcing[] = {7,7,7}, jac[9] = {0};
    size_t n = 2;
    if (which == 0) return !strcmp(c.name,"conservation") && c.constant == 1 && c.reactant_count == 2 && terms[0].species == 0 && terms[1].species == 1 && chem_constraint_validate(&c,n) == CHEM_CONSTRAINT_OK;
    if (which == 1) { c.row = 2; c.reactant_count = 3; c.constant = 10; return c.row == 2 && chem_constraint_validate(&c,3) == CHEM_CONSTRAINT_OK; }
    if (which == 2 || which == 7 || which == 13) {
        n=3; c.row=2; c.reactant_count=3; c.constant=5;
        terms[0].yield=2; terms[1].yield=3; terms[2].yield=-1;
        y[0]=1; y[1]=2; y[2]=3;
        if (which == 2) return c.constant == 5 && c.reactant_count == 3 && terms[0].yield == 2 && terms[1].yield == 3 && terms[2].yield == -1;
    }
    if (which == 3 || which == 9) { c.constant=0; terms[1].yield=-1; y[0]=y[1]=.5; if (which == 3) return c.constant == 0 && c.reactant_count == 2; }
    if (which == 4 || which == 10) { c.constant=2; terms[0].yield=.5; terms[1].yield=1.5; y[0]=y[1]=1; if (which == 4) return terms[0].yield == .5 && terms[1].yield == 1.5 && c.constant == 2; }
    if (which == 8) { n=3; c.row=2; c.reactant_count=3; c.constant=10; y[0]=2; y[1]=3; y[2]=5; }
    if (which == 11) { terms[0].yield=2; terms[1].yield=3; y[0]=.1; y[1]=.2; }
    if (which == 13) { y[0]=1.5; y[1]=.8; y[2]=2; }
    if (chem_constraint_apply(&c,n,300,y,forcing,jac) != CHEM_CONSTRAINT_OK) return 0;
    if (which == 5) {
        if (fabs(forcing[1]) > 1e-10) return 0;
        y[0]=.5; y[1]=.6;
        return !chem_constraint_apply(&c,n,300,y,forcing,jac) && fabs(forcing[1]-.1) < 1e-10;
    }
    if (which == 6) return fabs(jac[2]+1)<1e-10 && fabs(jac[3]+1)<1e-10;
    if (which == 7) return fabs(forcing[2]) < 1e-10 && fabs(jac[6]+2)<1e-10 && fabs(jac[7]+3)<1e-10 && fabs(jac[8]-1)<1e-10;
    if (which == 8) {
        if (fabs(forcing[2])>1e-10) return 0;
        y[2]=6; return !chem_constraint_apply(&c,n,300,y,forcing,jac) && fabs(forcing[2]-1)<1e-10;
    }
    if (which == 9) {
        if (fabs(forcing[1])>1e-10) return 0;
        y[0]=.6; return !chem_constraint_apply(&c,n,300,y,forcing,jac) && fabs(forcing[1]-.1)<1e-10;
    }
    if (which == 10) return fabs(forcing[1])<1e-10 && fabs(jac[2]+.5)<1e-10 && fabs(jac[3]+1.5)<1e-10;
    if (which == 11) {
        if (fabs(jac[2]+2)>1e-10 || fabs(jac[3]+3)>1e-10) return 0;
        y[0]=10; y[1]=20;
        return !chem_constraint_apply(&c,n,300,y,forcing,jac) && fabs(jac[2]+2)<1e-10 && fabs(jac[3]+3)<1e-10;
    }
    if (which == 12 || which == 13) {
        int cells = which == 12 ? 2 : 1;
        for (int cell=0; cell<cells; ++cell) {
            if (which == 12 && cell == 1) { y[0]=.8; y[1]=.2; }
            if (chem_constraint_apply(&c,n,300,y,forcing,jac)) return 0;
            double base = forcing[c.row];
            for (size_t col=0; col<n; ++col) {
                double old=y[col], f[3]={0}, j[9]={0}, h=1e-6;
                y[col] += h;
                if (chem_constraint_apply(&c,n,300,y,f,j)) return 0;
                y[col] = old;
                if (fabs((f[c.row]-base)/h + jac[c.row*n+col]) > 1e-8) return 0;
            }
        }
        return 1;
    }
    return 0;
}

static int equilibrium_case(int which)
{
    ChemTerm reactants[] = {{0,1},{1,1}}, products[] = {{2,1},{2,1}};
    ChemConstraint c = {"A_B_equilibrium", CHEM_EQUILIBRIUM_CONSTRAINT, 2, reactants, products, 2, 1, 1000, -2400, 8.31446261815324, 298.15};
    if (which == 0) return !strcmp(c.name,"A_B_equilibrium") && c.reactant_count == 2 && c.product_count == 1 && reactants[0].species == 0 && reactants[1].species == 1 && products[0].species == 2 && chem_constraint_validate(&c,3) == CHEM_CONSTRAINT_OK;
    if (which == 1) return c.row == 2 && chem_constraint_validate(&c,3) == CHEM_CONSTRAINT_OK;
    if (which == 2) { c.name="A_B_simple"; c.row=1; c.reactant_count=1; products[0].species=1; c.constant=10;
        return !strcmp(c.name,"A_B_simple") && c.row==1 && c.reactant_count==1 && c.product_count==1 && reactants[0].species==0 && products[0].species==1 && chem_constraint_validate(&c,2)==CHEM_CONSTRAINT_OK; }
    if (which == 3 || which == 10 || which == 12) {
        c.name="dissociation"; c.row=1; c.reactant_count=1; c.product_count=2;
        reactants[0].yield=2; products[0].species=1; products[1].species=2; c.constant=100;
        if (which == 3) return !strcmp(c.name,"dissociation") && c.row==1 && c.reactant_count==1 && c.product_count==2 && reactants[0].yield==2 && reactants[0].species==0 && products[0].species==1 && products[1].species==2 && chem_constraint_validate(&c,3)==CHEM_CONSTRAINT_OK;
    }
    if (which == 4) { c.constant=-1; if (chem_constraint_validate(&c,3)!=CHEM_CONSTRAINT_BAD_EQUILIBRIUM) return 0; c.constant=0; return chem_constraint_validate(&c,3)==CHEM_CONSTRAINT_BAD_EQUILIBRIUM; }
    if (which == 5) { c.reactant_count=0; return chem_constraint_validate(&c,3)==CHEM_CONSTRAINT_EMPTY_REACTANTS; }
    if (which == 6) { c.product_count=0; return chem_constraint_validate(&c,3)==CHEM_CONSTRAINT_EMPTY_PRODUCTS; }
    if (which == 7) {
        reactants[0].yield=0; if (chem_constraint_validate(&c,3)!=CHEM_CONSTRAINT_BAD_STOICHIOMETRY) return 0;
        reactants[0].yield=-1; if (chem_constraint_validate(&c,3)!=CHEM_CONSTRAINT_BAD_STOICHIOMETRY) return 0;
        reactants[0].yield=1; products[0].yield=0; if (chem_constraint_validate(&c,3)!=CHEM_CONSTRAINT_BAD_STOICHIOMETRY) return 0;
        products[0].yield=-2; return chem_constraint_validate(&c,3)==CHEM_CONSTRAINT_BAD_STOICHIOMETRY;
    }
    double y[3]={.001,.001,.001}, f[3]={0}, j[9]={0};
    if (which == 9) { y[0]=.01; y[1]=.02; y[2]=.05; }
    if (which == 10) { y[0]=.1; y[1]=.5; y[2]=2; }
    if (which == 11) { y[0]=.02; y[1]=.03; y[2]=.05; }
    if (which == 12) { y[0]=.15; y[1]=.8; y[2]=1.5; }
    if (chem_constraint_apply(&c,3,298.15,y,f,j)) return 0;
    if (which == 8) {
        if (fabs(f[2])>1e-10) return 0;
        y[0]=.02; y[1]=.03; y[2]=.05;
        return !chem_constraint_apply(&c,3,298.15,y,f,j) && fabs(f[2]-.55)<1e-10;
    }
    if (which == 9) return fabs(j[6]+20)<1e-8 && fabs(j[7]+10)<1e-8 && fabs(j[8]-1)<1e-8;
    if (which == 10) return fabs(f[1])<1e-10;
    if (which == 11 || which == 12) {
        int cells=which==11 ? 2 : 1;
        for (int cell=0; cell<cells; ++cell) {
            if (cell) { y[0]=.1; y[1]=.2; y[2]=5; }
            if (chem_constraint_apply(&c,3,298.15,y,f,j)) return 0;
            for (int col=0; col<3; ++col) {
                double old=y[col], fp[3]={0}, fm[3]={0}, unused[9]={0}, h=1e-6;
                y[col]=old+h; if (chem_constraint_apply(&c,3,298.15,y,fp,unused)) return 0;
                y[col]=old-h; if (chem_constraint_apply(&c,3,298.15,y,fm,unused)) return 0;
                y[col]=old;
                if (fabs((fp[c.row]-fm[c.row])/(2*h)+j[c.row*3+col]) > 1e-7) return 0;
            }
        }
        return 1;
    }
    return 0;
}

static int arrhenius(void)
{
    double p[] = {1, 0, 0, 300, 0};
    if (fabs(rate(CHEM_ARRHENIUS, p, 5) - 1) > 1e-13) return 0;
    p[0] = 1;
    if (fabs(rate(CHEM_ARRHENIUS, p, 5) - 1) > 1e-13) return 0;
    p[0] = 2.2e-10;
    if (fabs(rate(CHEM_ARRHENIUS, p, 5) - p[0]) > 1e-13*p[0]) return 0;
    p[0] = 3e-11; p[2] = -200;
    double expected = p[0]*exp(-200/301.24);
    if (fabs(rate(CHEM_ARRHENIUS, p, 5)-expected) > 1e-13*expected) return 0;
    p[0] = 2.6e-12; p[2] = -350;
    expected = p[0]*exp(-350/301.24);
    return fabs(rate(CHEM_ARRHENIUS, p, 5)-expected) <= 1e-13*expected;
}

static int tunneling(int full)
{
    double p[] = {1, 0, 0};
    if (full) { p[0] = 1.2; p[1] = 2.3; p[2] = 302.3; }
    double expected = p[0]*exp(-p[1]/301.24)*exp(p[2]/pow(301.24, 3));
    return fabs(rate(CHEM_TUNNELING, p, 3)-expected) <= 1e-13*expected;
}

int main(void)
{
    printf("arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithSystem %d\n", arrhenius());
    printf("arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithPrescribedArguments %d\n", arrhenius());
    printf("tunneling_rate_constant:TunnelingRateConstant.CalculateWithMinimalArguments %d\n", tunneling(0));
    printf("tunneling_rate_constant:TunnelingRateConstant.CalculateWithAllArguments %d\n", tunneling(1));
    printf("user_defined_rate_constant:UserDefinedRateConstant.CalculateWithSystem %d\n", user_rate(.5,1));
    printf("user_defined_rate_constant:UserDefinedRateConstant.ConstructorWithRate %d\n", user_rate(1.1,1));
    printf("user_defined_rate_constant:UserDefinedRateConstant.ConstructorWithRateAndName %d\n", user_rate(1.1,1));
    printf("user_defined_rate_constant:UserDefinedRateConstant.CustomScalingFactor %d\n", user_rate(1.2,2));
    printf("reversible_rate_constant:ReversibleRateConstant.DefaultConstructor %d\n", reversible(1,0,0,298.15));
    printf("reversible_rate_constant:ReversibleRateConstant.CalculateWithSystem %d\n", reversible(1.14e-2,2300,.32,301.24));
    printf("reversible_rate_constant:ReversibleRateConstant.CalculateWithPrescribedArguments %d\n", reversible(1.14e-2,2300,.32,298.15));
    printf("reversible_rate_constant:ReversibleRateConstant.TemperatureDependence %d\n", reversible(1.5,2000,.5,273.15) && reversible(1.5,2000,.5,298.15) && reversible(1.5,2000,.5,323.15));
    printf("reversible_rate_constant:ReversibleRateConstant.ZeroReverseRate %d\n", reversible(5,1000,0,298.15));
    printf("reversible_rate_constant:ReversibleRateConstant.EquilibriumConstantOnly %d\n", reversible(3.5,1800,1,298.15));
    printf("taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithSystem %d\n", taylor());
    printf("taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithPrescribedArguments %d\n", taylor());
    printf("troe_rate_constant:TroeRateConstant.CalculateWithMinimalArguments %d\n", falloff(CHEM_TROE,0));
    printf("troe_rate_constant:TroeRateConstant.CalculateWithAllArguments %d\n", falloff(CHEM_TROE,1));
    printf("troe_rate_constant:TroeRateConstant.AnalyticalTroeExampleAB %d\n", falloff(CHEM_TROE,2));
    printf("troe_rate_constant:TroeRateConstant.AnalyticalTroeExampleBC %d\n", falloff(CHEM_TROE,3));
    printf("ternary_chemical_activation_rate_constant:TernaryChemicalActivationRateConstant.CalculateWithMinimalArguments %d\n", falloff(CHEM_TERNARY,0));
    printf("ternary_chemical_activation_rate_constant:TernaryChemicalActivationRateConstant.CalculateWithAllArguments %d\n", falloff(CHEM_TERNARY,1));
    printf("ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.CalculateWithAllArguments %d\n", jpl(0));
    printf("ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.ZeroKintReducesToTroe %d\n", jpl(1));
    printf("branched_rate_constant:BranchedRateConstant.CalculateAlkoxyBranchWithAllArguments %d\n", branched(0));
    printf("branched_rate_constant:BranchedRateConstant.CalculateNitrateBranchWithAllArguments %d\n", branched(1));
    printf("surface_rate_constant:SurfaceRateConstant.CalculateDefaultProbability %d\n", surface(1));
    printf("surface_rate_constant:SurfaceRateConstant.CalculateSpecifiedProbability %d\n", surface(.74));
    printf("surface_rate_constant:SurfaceRateConstant.DiffusionCoefficientIsMissing %d\n", missing_surface(0));
    printf("surface_rate_constant:SurfaceRateConstant.MolecularWeightIsMissing %d\n", missing_surface(1));
    printf("reaction_rate_store:ReactionRateConstantStore.SurfaceMissingDiffusionCoefficientThrows %d\n", missing_surface(0));
    printf("reaction_rate_store:ReactionRateConstantStore.SurfaceMissingMolecularWeightThrows %d\n", missing_surface(1));
    printf("backward_euler:BackwardEuler.CanCallSolve %d\n", backward_solve());
    printf("backward_euler:BackwardEuler.IsConverged %d\n", backward_convergence());
    printf("linear:LinearConstraint.Construction %d\n", linear_case(0));
    printf("linear:LinearConstraint.AlgebraicSpecies %d\n", linear_case(1));
    printf("linear:LinearConstraint.WeightedTerms %d\n", linear_case(2));
    printf("linear:LinearConstraint.ZeroConstant %d\n", linear_case(3));
    printf("linear:LinearConstraint.FractionalCoefficients %d\n", linear_case(4));
    printf("linear:LinearConstraint.ResidualComputationThroughConstraintSet %d\n", linear_case(5));
    printf("linear:LinearConstraint.JacobianComputationThroughConstraintSet %d\n", linear_case(6));
    printf("linear:LinearConstraint.WeightedSumResidualAndJacobian %d\n", linear_case(7));
    printf("linear:LinearConstraint.ThreeSpeciesConservationResidual %d\n", linear_case(8));
    printf("linear:LinearConstraint.ZeroConstantResidual %d\n", linear_case(9));
    printf("linear:LinearConstraint.FractionalCoefficientsResidualAndJacobian %d\n", linear_case(10));
    printf("linear:LinearConstraint.JacobianIndependentOfConcentrations %d\n", linear_case(11));
    printf("linear:LinearConstraint.FiniteDifferenceJacobianSimpleConservation %d\n", linear_case(12));
    printf("linear:LinearConstraint.FiniteDifferenceJacobianWeightedSum %d\n", linear_case(13));
    printf("equilibrium:EquilibriumConstraint.Construction %d\n", equilibrium_case(0));
    printf("equilibrium:EquilibriumConstraint.AlgebraicSpecies %d\n", equilibrium_case(1));
    printf("equilibrium:EquilibriumConstraint.SingleReactantSingleProduct %d\n", equilibrium_case(2));
    printf("equilibrium:EquilibriumConstraint.MultipleReactantsAndProducts %d\n", equilibrium_case(3));
    printf("equilibrium:EquilibriumConstraint.InvalidEquilibriumConstant %d\n", equilibrium_case(4));
    printf("equilibrium:EquilibriumConstraint.EmptyReactantsThrows %d\n", equilibrium_case(5));
    printf("equilibrium:EquilibriumConstraint.EmptyProductsThrows %d\n", equilibrium_case(6));
    printf("equilibrium:EquilibriumConstraint.InvalidStoichiometryThrows %d\n", equilibrium_case(7));
    printf("equilibrium:EquilibriumConstraint.ResidualComputationThroughConstraintSet %d\n", equilibrium_case(8));
    printf("equilibrium:EquilibriumConstraint.JacobianComputationThroughConstraintSet %d\n", equilibrium_case(9));
    printf("equilibrium:EquilibriumConstraint.ComplexStoichiometryResidual %d\n", equilibrium_case(10));
    printf("equilibrium:EquilibriumConstraint.FiniteDifferenceJacobianSimple %d\n", equilibrium_case(11));
    printf("equilibrium:EquilibriumConstraint.FiniteDifferenceJacobianComplexStoichiometry %d\n", equilibrium_case(12));
    return 0;
}
