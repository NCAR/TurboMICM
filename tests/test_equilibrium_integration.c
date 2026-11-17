#include "state.h"
#include <float.h>
#include <math.h>
#include <stdio.h>

static double eq(double constant,double enthalpy,double temperature)
{
    return constant*exp((enthalpy/(1.380649e-23*6.02214076e23))*(1/temperature-1/298.15));
}
static int mixed_parameter_order(void)
{
    const char *names[]={"A"},*state_labels[]={"constraint","photolysis"},*rate_labels[]={"photolysis"};
    ChemReaction reaction={.type=CHEM_USER,.parameter=0,.p={2}};
    ChemMechanism mechanism={.species_count=1,.reaction_count=1,.parameter_count=1,.species_names=names,.parameter_names=rate_labels,.reactions=&reaction};
    ChemStateParameters params={1,1,2,names,state_labels,1e-6,NULL}; ChemState state={0};
    if(chem_state_init(&state,&params,1,1)) return 0;
    state.conditions[0]=(ChemConditions){300,101325,0,0};
    *chem_matrix_at(&state.custom_parameters,0,0)=100;
    *chem_matrix_at(&state.custom_parameters,0,1)=.3;
    int ok=!chem_state_update_rate_constants(&state,&mechanism,NULL) && *chem_matrix_at(&state.rate_constants,0,0)==.6;
    chem_state_free(&state);return ok;
}
static int run(int kind)
{
    const char *names3[]={"A","B","C"},*names6[]={"A","B","C","D","E","F"},*names4[]={"A","B","C","D"},*reordered[]={"B","C","A"};
    const char *ordered_names[3];
    const char **names=kind==1?names6:kind==4?names4:kind==3?ordered_names:names3;
    size_t n=kind==1?6:kind==4?4:3,nconstraints=kind==1 || kind==4?2:1;
    size_t A=0,B=1,C=2;
    if(kind==3) {
        const unsigned char pattern[]={1,0,1, 1,1,0, 0,0,1}; size_t order[3],inverse[3];
        if(chem_diagonal_markowitz_order(3,pattern,order)) return 0;
        for(size_t i=0;i<3;++i) {ordered_names[i]=reordered[order[i]];inverse[order[i]]=i;}
        A=inverse[2]; B=inverse[0]; C=inverse[1];
    }
    const char *labels[]={"B_C_eq","E_F_eq"}; if(kind==4) labels[1]="B_D_eq";if(kind==5) labels[0]="A2_B_eq";
    double gas=1.380649e-23*6.02214076e23, temperature=kind==1?310:kind==2?270:kind==3?400:kind==5?298:300;
    double k1=kind==0||kind==1?.034:kind==4?3:kind==5?10:2,k2=kind==1?.012:5;
    ChemTerm reactants[2]={{(int)(kind==5?C:A),1},{3,1}},products[2]={{(int)(kind==5?A:B),1},{4,1}};
    ChemReaction reactions[2]={{.type=CHEM_ARRHENIUS,.reactant_count=1,.product_count=1,.parameter=-1,.p={0,0,0,300}},
        {.type=CHEM_ARRHENIUS,.reactant_begin=1,.reactant_count=1,.product_begin=1,.product_count=1,.parameter=-1,.p={.2,0,0,300}}};
    reactions[0].p[0]=kind==0?.1:kind==1?.5:kind==5?.5:1;
    ChemMechanism mechanism={.species_count=n,.reaction_count=kind==1?2u:1u,.species_names=names,.reactions=reactions,.reactants=reactants,.products=products};
    ChemTerm eq_react[2]={{(int)(kind==5?A:B),kind==5?2:1},{kind==1?4:(int)B,1}};
    ChemTerm eq_product[2]={{(int)(kind==5?B:C),1},{kind==1?5:3,1}};
    ChemConstraint constraints[2]={{labels[0],CHEM_EQUILIBRIUM_CONSTRAINT,(size_t)eq_product[0].species,&eq_react[0],&eq_product[0],1,1,k1,-2400,gas,298.15},
        {labels[1],CHEM_EQUILIBRIUM_CONSTRAINT,(size_t)eq_product[1].species,&eq_react[1],&eq_product[1],1,1,k2,kind==1?-2000:-2400,gas,298.15}};
    ChemConstraintSet set={constraints,nconstraints}; unsigned char mass[6];
    ChemStateParameters state_parameters={mechanism.reaction_count,n,nconstraints,names,labels,1e-6,NULL};ChemState state={0};
    if(chem_state_init(&state,&state_parameters,1,4) || chem_constraint_mass_diagonal(set,n,mass)) return 0;
    ChemWorkspace *w=chem_workspace_create_constrained(&mechanism,constraints,nconstraints);if(!w) {chem_state_free(&state);return 0;}
    state.conditions[0]=(ChemConditions){temperature,101325,0,0};
    int ok=state.parameters.variable_count==n;
    for(size_t i=0;i<n;++i) {int algebraic=0;for(size_t j=0;j<nconstraints;++j) if(constraints[j].row==i) algebraic=1;if(mass[i]!=(unsigned char)!algebraic) ok=0;}
    if(kind<=1) {
        *chem_state_variable(&state,0,"A")=1;*chem_state_variable(&state,0,"B")=.1;*chem_state_variable(&state,0,"C")=.1*k1;
        if(kind==1) {*chem_state_variable(&state,0,"D")=.5;*chem_state_variable(&state,0,"E")=.05;*chem_state_variable(&state,0,"F")=.05*k2;}
    } else if(kind==5) { *chem_state_variable(&state,0,"A")=.1;*chem_state_variable(&state,0,"B")=.1;*chem_state_variable(&state,0,"C")=1; }
    else *chem_state_variable(&state,0,"A")=1;
    if(chem_state_update_constraint_parameters(&state,set) || chem_state_update_rate_constants(&state,&mechanism,NULL)) ok=0;
    for(size_t i=0;i<nconstraints;++i) {
        double expected=eq(constraints[i].constant,constraints[i].enthalpy,temperature),actual=*chem_matrix_at(&state.custom_parameters,0,i);
        double tolerance=kind==0 || (kind==1 && i==1)?fabs(expected)*100*DBL_EPSILON:1e-10;
        if(fabs(actual-expected)>tolerance) ok=0;
    }
    if(kind>=2) {
        ChemRosenbrockParameters params=chem_rosenbrock_three_stage();double rates[2]={0},y[6]={0};
        chem_matrix_get_row(&state.rate_constants,0,rates,mechanism.reaction_count);chem_matrix_get_row(&state.variables,0,y,n);
        for(int step=0;step<(kind==5?50:100) && ok;++step) {
            ChemStatistics stats;
            if(chem_state_update_constraint_parameters(&state,set) || chem_rosenbrock_solve(&mechanism,w,&params,state.conditions[0],rates,.001,y,&stats)) {ok=0;break;}
            for(size_t i=0;i<nconstraints;++i) {
                double kval=*chem_matrix_at(&state.custom_parameters,0,i);
                double left=kval*pow(y[eq_react[i].species],eq_react[i].yield),right=y[eq_product[i].species];
                if(fabs(left-right)>(kind==5?1e-5:1e-6)) ok=0;
            }
        }
        if(kind!=5 && fabs(y[A]+y[B]-1)>.01) ok=0;
    }
    chem_workspace_free(w);chem_state_free(&state);return ok;
}
int main(void)
{
    const char *names[]={"SetConstraintsAPIWorks","SetConstraintsAPIMultipleConstraints","DAESolveWithConstraint","DAESolveWithConstraintAndReorderState","DAESolveWithTwoCoupledConstraints","DAESolveWithNonUnitStoichiometry"};
    for(int i=0;i<6;++i) printf("equilibrium_constraint:EquilibriumIntegration.%s %d\n",names[i],run(i) && (i!=0 || mixed_parameter_order()));
    return 0;
}
