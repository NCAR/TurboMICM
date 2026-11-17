#include "rosenbrock.h"
#include "reorder.h"
#include "state.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct { size_t n, source, product; double rate; } Reaction;
static void rhs(void *context,const double *y,double *f)
{
    Reaction *r=context; memset(f,0,r->n*sizeof(double));
    f[r->source]-=r->rate*y[r->source]; f[r->product]+=r->rate*y[r->source];
}
static void jac(void *context,const double *y,double *j)
{
    Reaction *r=context; (void)y; memset(j,0,r->n*r->n*sizeof(double));
    j[r->source*r->n+r->source]-=r->rate; j[r->product*r->n+r->source]+=r->rate;
}
static void run(ChemConstraintSet *set, Reaction *reaction, double *initial, size_t steps, double dt, double temperature)
{
    size_t n=set->variables;
    unsigned char pattern[16]={0};
    pattern[reaction->source*n+reaction->source]=1;
    pattern[reaction->product*n+reaction->source]=1;
    chem_constraint_sparsity(set,pattern);
    ChemSparsePlan *plan=chem_sparse_plan_from_pattern(n,pattern);
    assert(plan);
    const double atol[]={1e-3,1e-3,1e-3,1e-3};
    ChemRosenbrockMethod method=chem_rosenbrock_three_stage();
    for (size_t sparse=0;sparse<2;++sparse)
      for (size_t algorithm=0;algorithm<2;++algorithm) {
        ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
        options.constraints=set; options.constraint_temperature=temperature;
        options.sparse_plan=sparse ? plan : NULL; options.lu_algorithm=(ChemLuAlgorithm)algorithm;
        double y[4]; memcpy(y,initial,n*sizeof(double));
        for (size_t step=0;step<steps;++step) {
            assert(chem_rosenbrock_solve(n,y,dt,rhs,jac,reaction,&method,&options,NULL)==CHEM_SOLVE_OK);
            for (size_t i=0;i<set->count;++i) {
                const ChemConstraint *c=set->constraints+i;
                double k=c->constant*exp((c->enthalpy/c->gas_constant)*(1/temperature-1/c->reference_temperature));
                double left=k,right=1;
                for (size_t j=0;j<c->reactant_count;++j) left*=pow(y[c->reactants[j].species],c->reactants[j].coefficient);
                for (size_t j=0;j<c->product_count;++j) right*=pow(y[c->products[j].species],c->products[j].coefficient);
                assert(fabs(left-right)<1e-6);
            }
        }
        assert(fabs(y[reaction->source]+y[reaction->product]-initial[reaction->source]-initial[reaction->product])<.01);
      }
    chem_sparse_plan_destroy(plan);
}
typedef struct { const ChemMechanism *mechanism; double rate; } ReorderedReaction;
static void reordered_rhs(void *context,const double *y,double *out)
{ ReorderedReaction *reaction=context; chem_rhs(reaction->mechanism,y,&reaction->rate,out); }
static void reordered_jac(void *context,const double *y,double *out)
{ ReorderedReaction *reaction=context; chem_jacobian(reaction->mechanism,y,&reaction->rate,out); }
static void automatically_reordered(void)
{
    /* Start with B,A,C so the computed ordering changes the state indices. */
    const char *original_names[]={"B","A","C"};
    const ChemConstraintTerm b[]={{0,1}}, c[]={{2,1}};
    const ChemConstraint row={"B_C_eq",CHEM_CONSTRAINT_EQUILIBRIUM,2,b,c,1,1,2,-2400,8.31446261815324,298.15,NULL};
    const ChemConstraintSet original={3,1,&row};
    const size_t reactants[]={1}; const ChemProduct products[]={{0,1}};
    const ChemReaction reactions[]={{0,1,0,1}};
    const ChemMechanism mechanism={3,1,reactions,reactants,products};
    ChemPermutation permutation={0};
    assert(chem_permutation_for_system(&permutation,&mechanism,&original)==CHEM_REORDER_OK);
    assert(permutation.old_to_new[1]!=1);
    for (size_t i=0;i<3;++i) assert(permutation.new_to_old[permutation.old_to_new[i]]==i);
    ChemReorderedMechanism reordered_mechanism={0};
    assert(chem_reordered_mechanism_init(&reordered_mechanism,&mechanism,&permutation)==CHEM_REORDER_OK);
    ChemReorderedConstraints reordered={0};
    assert(chem_reordered_constraints_init(&reordered,&original,&permutation)==CHEM_REORDER_OK);
    const char *names[3];
    for (size_t i=0;i<3;++i) names[i]=original_names[permutation.new_to_old[i]];
    const double atol[]={.001,.001,.001};
    const ChemStateConfig config={1,3,0,0,names,NULL,1e-6,atol};
    ChemState state={0};
    assert(chem_state_init_constrained(&state,&config,0,NULL,&reordered.set)==CHEM_STATE_OK);
    size_t a=chem_state_species_index(&state,"A"),bb=chem_state_species_index(&state,"B"),cc=chem_state_species_index(&state,"C");
    size_t k=chem_state_parameter_index(&state,"B_C_eq");
    assert(a==permutation.old_to_new[1] && bb==permutation.old_to_new[0] && cc==permutation.old_to_new[2] && k==0);
    state.concentrations[a]=1; state.conditions[0]=(ChemConditions){400,101325,0};
    assert(chem_state_update_constraint_parameters(&state,&reordered.set)==CHEM_STATE_OK);
    const double expected=2*exp(-2400/8.31446261815324*(1.0/400-1.0/298.15));
    assert(fabs(state.parameters[k]-expected)<1e-10);
    ReorderedReaction reaction={&reordered_mechanism.mechanism,1};
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(state.absolute_tolerances);
    options.constraints=&reordered.set; options.constraint_parameters=state.parameters; options.constraint_temperature=400;
    ChemRosenbrockMethod method=chem_rosenbrock_three_stage();
    for (size_t step=0;step<100;++step) {
        assert(chem_state_update_constraint_parameters(&state,&reordered.set)==CHEM_STATE_OK);
        assert(chem_rosenbrock_solve(3,state.concentrations,.001,reordered_rhs,reordered_jac,&reaction,&method,&options,NULL)==CHEM_SOLVE_OK);
        assert(fabs(state.parameters[k]*state.concentrations[bb]-state.concentrations[cc])<1e-6);
    }
    assert(fabs(state.concentrations[a]+state.concentrations[bb]-1)<.01);
    chem_state_destroy(&state); chem_reordered_constraints_destroy(&reordered);
    chem_reordered_mechanism_destroy(&reordered_mechanism); chem_permutation_destroy(&permutation);
}
int main(void)
{
    automatically_reordered();
    const ChemConstraintTerm b[]={{1,1}}, c[]={{2,1}}, d[]={{3,1}};
    ChemConstraint constraints[]={
        {"B_C",CHEM_CONSTRAINT_EQUILIBRIUM,2,b,c,1,1,2,-2400,8.31446261815324,298.15,NULL},
        {"B_D",CHEM_CONSTRAINT_EQUILIBRIUM,3,b,d,1,1,5,-2400,8.31446261815324,298.15,NULL}
    };
    ChemConstraintSet set={3,1,constraints}; Reaction reaction={3,0,1,1};
    double y[]={1,0,0}; run(&set,&reaction,y,100,.001,270);
    const ChemConstraintTerm reordered_b[]={{2,1}}, reordered_c[]={{0,1}};
    constraints[0]=(ChemConstraint){"B_C",CHEM_CONSTRAINT_EQUILIBRIUM,0,reordered_b,reordered_c,1,1,2,-2400,8.31446261815324,298.15,NULL};
    reaction=(Reaction){3,1,2,1}; double reordered[]={0,1,0}; run(&set,&reaction,reordered,100,.001,400);
    constraints[0]=(ChemConstraint){"B_C",CHEM_CONSTRAINT_EQUILIBRIUM,2,b,c,1,1,3,-2400,8.31446261815324,298.15,NULL};
    reaction=(Reaction){3,0,1,1}; set.variables=4; set.count=2; reaction.n=4;
    double y4[]={1,0,0,0}; run(&set,&reaction,y4,100,.001,300);
    const ChemConstraintTerm a2[]={{0,2}}, bproduct[]={{1,1}};
    constraints[0]=(ChemConstraint){"A2_B",CHEM_CONSTRAINT_EQUILIBRIUM,1,a2,bproduct,1,1,10,-2400,8.31446261815324,298.15,NULL};
    set=(ChemConstraintSet){3,1,constraints}; reaction=(Reaction){3,2,0,.5};
    double yn[]={.1,.1,1}; run(&set,&reaction,yn,50,.001,298);
    puts("equilibrium DAE integration PASS"); return 0;
}
