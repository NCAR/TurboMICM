#include "chem.h"

#include <math.h>
#include <stdio.h>

#define CHECK(condition) do { if (!(condition)) { \
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); return 1; \
} } while (0)

static const ChemSpecies species[]={{"A",1,1e-3,0,0},{"B",.5,1e-3,0,0},{"C",5,1e-3,0,0}};
static const int reactants[]={0};
static const ChemProduct products[]={{1,1}};
static const ChemReaction reactions[]={{0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=.5,.d=300}}};
static const ChemMechanism mechanism={3,1,0,species,reactions,reactants,products,NULL,NULL,298.15,101325};
static const ChemTerm eq_reactants[]={{1,1}}, eq_products[]={{2,1}};
static const ChemConstraint equilibrium={"B_C_eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,
                                          eq_reactants,eq_products,1,1,0,2,0,0,0};

static int integration(ChemRosenbrock p, size_t cells)
{
    ChemState *s=chem_state_create_constrained(&mechanism,&equilibrium,1,cells);
    CHECK(s);
    for (size_t c=0;c<cells;++c) {
        chem_species_values(s,1)[c]=.5+.1*c;
        chem_species_values(s,2)[c]=5+10*c;
    }
    CHECK(chem_update_rates(s)==CHEM_OK);
    double f[3], j[9];
    CHECK(chem_evaluate(s,0,f,j)==CHEM_OK);
    CHECK(fabs(f[0]+.5)<1e-14 && fabs(f[1]-.5)<1e-14 && fabs(f[2]+4)<1e-14);
    CHECK(j[0]==.5 && j[3]==-.5 && j[7]==-2 && j[8]==1);
    ChemStats stats;
    int status=chem_solve(s,&p,.001,&stats);
    if (status) fprintf(stderr,"constrained solve: %s\n",chem_error_string(status));
    CHECK(status==CHEM_OK && stats.constraint_init_iterations>0);
    for (size_t c=0;c<cells;++c)
        CHECK(fabs(chem_species_values(s,2)[c]-2*chem_species_values(s,1)[c])<1e-6);
    chem_state_free(s);
    return 0;
}

static int scaled_projection(double scale)
{
    ChemTerm terms[]={{0,scale},{1,scale},{2,scale}};
    ChemConstraint constraint={"mass",CHEM_LINEAR_CONSTRAINT,2,terms,NULL,3,0,scale*1.3,0,0,0,0};
    ChemState *s=chem_state_create_constrained(&mechanism,&constraint,1,1);
    CHECK(s); chem_species_values(s,0)[0]=.7; chem_species_values(s,1)[0]=.15;
    CHECK(chem_update_rates(s)==CHEM_OK);
    ChemRosenbrock p=chem_rosenbrock_four_stage_dae(); ChemStats stats;
    CHECK(chem_solve(s,&p,.001,&stats)==CHEM_OK);
    CHECK(fabs(chem_species_values(s,2)[0]-.45)<4.5e-6);
    chem_state_free(s);
    return 0;
}

int main(void)
{
    double y[]={.02,.03,.05}, gradient[3], residual;
    ChemTerm left[]={{0,1},{1,1}}, right[]={{2,1}};
    ChemConstraint c={"eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,left,right,2,1,0,1000,0,0,0};
    CHECK(chem_constraint_evaluate(&c,y,3,298.15,&residual,gradient)==CHEM_OK);
    CHECK(fabs(residual-.55)<1e-12 && fabs(gradient[0]-30)<1e-12 &&
          fabs(gradient[1]-20)<1e-12 && gradient[2]==-1);
    for (size_t i=0;i<3;++i) {
        double h=1e-7, value=y[i], plus, minus;
        y[i]=value+h; CHECK(chem_constraint_evaluate(&c,y,3,298.15,&plus,NULL)==CHEM_OK);
        y[i]=value-h; CHECK(chem_constraint_evaluate(&c,y,3,298.15,&minus,NULL)==CHEM_OK);
        y[i]=value;
        CHECK(fabs((plus-minus)/(2*h)-gradient[i])<1e-7);
    }
    CHECK(integration(chem_rosenbrock_three_stage(),1)==0);
    CHECK(integration(chem_rosenbrock_four_stage_dae(),4)==0);
    CHECK(integration(chem_rosenbrock_six_stage_dae(),2)==0);
    CHECK(scaled_projection(1)==0 && scaled_projection(1e-12)==0 && scaled_projection(1e12)==0);
    ChemState *s=chem_state_create_constrained(&mechanism,&equilibrium,1,1); CHECK(s);
    CHECK(chem_update_rates(s)==CHEM_OK);
    ChemRosenbrock p=chem_rosenbrock_three_stage(); p.constraint_max_iterations=0;
    ChemStats stats;
    CHECK(chem_solve(s,&p,.001,&stats)==CHEM_CONSTRAINT_INITIALIZATION_FAILED);
    CHECK(chem_species_values(s,0)[0]==1 && chem_species_values(s,1)[0]==.5 && chem_species_values(s,2)[0]==5);
    chem_state_free(s);
    puts("constraint residuals, Jacobians, projection, and DAE integration: PASS");
    return 0;
}
