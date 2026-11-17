#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "model.h"
#include "rosenbrock.h"

typedef struct { size_t n,source,target,solvent; double k; int use_solvent; size_t finalized; } Transfer;
static const ChemState *finalized_state;
static void forcing(void *v,const double *y,double *f)
{
    Transfer *t=v; double rate=t->k*y[t->source]*(t->use_solvent ? y[t->solvent] : 1);
    f[t->source]-=rate; f[t->target]+=rate;
}
static void jacobian(void *v,const double *y,double *j)
{
    Transfer *t=v; double d=t->k*(t->use_solvent ? y[t->solvent] : 1);
    j[t->source*t->n+t->source]-=d; j[t->target*t->n+t->source]+=d;
    if (t->use_solvent) {
        d=t->k*y[t->source]; j[t->source*t->n+t->solvent]-=d; j[t->target*t->n+t->solvent]+=d;
    }
}
static void finalize(void *v,const ChemState *state,const ChemSparseMatrix *structure)
{
    Transfer *t=v; ++t->finalized;
    finalized_state=state;
    assert(state->species_count==t->n && structure->dimension==t->n);
    assert(chem_sparse_matrix_entry(structure,t->source,t->source)!=SIZE_MAX);
    assert(chem_sparse_matrix_entry(structure,t->target,t->source)!=SIZE_MAX);
    if (t->use_solvent) assert(chem_sparse_matrix_entry(structure,t->target,t->solvent)!=SIZE_MAX);
}
static void finalize_sparse_only(void *v,const ChemState *state,const ChemSparseMatrix *structure)
{
    Transfer *t=v; (void)state; ++t->finalized;
    assert(structure->dimension==t->n && structure->values.width==4);
    assert(chem_sparse_matrix_entry(structure,t->target,t->source)!=SIZE_MAX);
}
static ChemOdeComponent component(Transfer *t)
{ return (ChemOdeComponent){t,NULL,forcing,jacobian}; }
static void solve(ChemModelSystem *m,double *y,size_t steps,double dt)
{
    ChemRosenbrockMethod method=chem_rosenbrock_four_stage_da();
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(m->state.absolute_tolerances);
    options.constraints=m->constraints.count ? &m->constraints : NULL;
    options.constraint_temperature=298; options.sparse_plan=m->plan;
    for (size_t i=0;i<steps;++i)
        assert(chem_rosenbrock_solve(m->state.species_count,y,dt,chem_composite_rhs,
                    chem_composite_jacobian,&m->ode,&method,&options,NULL)==CHEM_SOLVE_OK);
}
typedef struct { size_t reactant,product; double k; } Equilibrium;
static double eq_residual(void *v,const double *y,double t,const double *p)
{ Equilibrium *e=v; (void)t; (void)p; return e->k*y[e->reactant]-y[e->product]; }
static void eq_jac(void *v,const double *y,double t,const double *p,double *row)
{ Equilibrium *e=v; (void)y; (void)t; (void)p; row[e->reactant]-=e->k; row[e->product]+=1; }

static void metadata_and_patterns(void)
{
    const char *const gas[]={"A_GAS"},*const aerosol[]={"AEROSOL.A_AQ"};
    const double base_atol[]={.001},external_atol[]={.001};
    ChemStateConfig base={1,1,0,0,gas,NULL,1e-6,base_atol};
    const ChemConstraintTerm mass[]={{0,1},{1,1}};
    ChemConstraint conservation={.name="mass",.type=CHEM_CONSTRAINT_LINEAR,.row=1,
        .reactants=mass,.reactant_count=2,.constant=1};
    const ChemSparseEntry entries[]={{0,0},{1,0}};
    Transfer transfer={2,0,1,0,.01,0,0};
    ChemExternalModel model={{1,0,aerosol,NULL,external_atol},component(&transfer),2,entries,1,&conservation,finalize};
    for (size_t width=0;width<=4;width+=4) for (int order=0;order<2;++order) {
        ChemModelSystem built={0}; transfer.finalized=0;
        assert(chem_model_system_init(&built,&base,NULL,NULL,NULL,1,&model,width,(ChemSparseOrdering)order)==CHEM_MODEL_OK);
        assert(transfer.finalized==1 && built.state.species_count==2 && built.state.constraint_count==1);
        assert(finalized_state==&built.state);
        assert(built.structure.values.width==width && built.structure.ordering==(ChemSparseOrdering)order);
        assert(chem_state_species_index(&built.state,"A_GAS")==0 && chem_state_species_index(&built.state,"AEROSOL.A_AQ")==1);
        assert(!built.state.algebraic[0] && built.state.algebraic[1]);
        chem_model_system_destroy(&built);
        model.constraint_count=0; model.constraints=NULL;
        assert(chem_model_system_init(&built,&base,NULL,NULL,NULL,1,&model,width,(ChemSparseOrdering)order)==CHEM_MODEL_OK);
        assert(built.state.constraint_count==0 && !built.state.algebraic);
        chem_model_system_destroy(&built);
        model.constraint_count=1; model.constraints=&conservation;
    }
    {
        ChemModelSystem built={0}; transfer.finalized=0; model.finalize=finalize_sparse_only;
        assert(chem_model_system_init(&built,&base,NULL,NULL,NULL,1,&model,4,CHEM_SPARSE_CSR)==CHEM_MODEL_OK);
        assert(transfer.finalized==1); chem_model_system_destroy(&built); model.finalize=finalize;
    }
    {
        ChemModelSystem built={0}; base.absolute_tolerances=NULL; model.state.absolute_tolerances=NULL;
        assert(chem_model_system_init(&built,&base,NULL,NULL,NULL,1,&model,0,CHEM_SPARSE_CSR)==CHEM_MODEL_OK);
        assert(built.state.absolute_tolerances[0]==.001 && built.state.absolute_tolerances[1]==.001);
        const ChemSparseEntry invalid[]={{5,0}};
        size_t calls=transfer.finalized; double *old=built.state.concentrations;
        model.jacobian_entries=invalid; model.jacobian_count=1;
        assert(chem_model_system_init(&built,&base,NULL,NULL,NULL,1,&model,0,CHEM_SPARSE_CSR)==CHEM_MODEL_INVALID);
        assert(built.state.concentrations==old && transfer.finalized==calls);
        model.jacobian_entries=entries; model.jacobian_count=2;
        chem_model_system_destroy(&built); base.absolute_tolerances=base_atol; model.state.absolute_tolerances=external_atol;
    }
    const char *const names4[]={"A_GAS","B","C"};
    const double atol4[]={.001,.001,.001};
    base=(ChemStateConfig){1,3,0,0,names4,NULL,1e-6,atol4};
    const ChemConstraintTerm eq_r[]={{1,1}},eq_p[]={{2,1}},mass4[]={{0,1},{3,1}};
    ChemConstraint eq={.name="B_C_eq",.type=CHEM_CONSTRAINT_EQUILIBRIUM,.row=2,.reactants=eq_r,.products=eq_p,
        .reactant_count=1,.product_count=1,.constant=5,.enthalpy=0,.gas_constant=8.31446261815324,.reference_temperature=298.15};
    ChemConstraintSet builtin={4,1,&eq};
    conservation=(ChemConstraint){.name="mass",.type=CHEM_CONSTRAINT_LINEAR,.row=3,.reactants=mass4,.reactant_count=2,.constant=1};
    const ChemSparseEntry four_entries[]={{0,0},{3,0}};
    transfer=(Transfer){4,0,3,0,.01,0,0};
    model=(ChemExternalModel){{1,0,aerosol,NULL,external_atol},component(&transfer),2,four_entries,1,&conservation,finalize};
    Transfer kinetic={4,0,1,0,.05,0,0}; ChemOdeComponent base_process=component(&kinetic);
    unsigned char kinetic_pattern[16]={1,0,0,0,1,0,0,0};
    ChemModelSystem combined={0};
    assert(chem_model_system_init(&combined,&base,&base_process,kinetic_pattern,&builtin,1,&model,0,CHEM_SPARSE_CSR)==CHEM_MODEL_OK);
    assert(combined.state.constraint_count==2 && combined.state.algebraic[2] && combined.state.algebraic[3]);
    chem_model_system_destroy(&combined);
}

static void algebraic_pattern(void)
{
    const char *const gas[]={"A_GAS","S"},*const aerosol[]={"AEROSOL.A_AQ"};
    const double atol[]={.001,.001},ext_atol[]={.001};
    ChemStateConfig base={1,2,0,0,gas,NULL,1e-6,atol};
    const ChemConstraintTerm mass[]={{0,1},{2,1}};
    ChemConstraint conservation={.name="mass",.type=CHEM_CONSTRAINT_LINEAR,.row=2,.reactants=mass,.reactant_count=2,.constant=1};
    const ChemSparseEntry entries[]={{0,0},{0,1},{2,0},{2,1}};
    Transfer transfer={3,0,2,1,.1,1,0};
    ChemExternalModel model={{1,0,aerosol,NULL,ext_atol},component(&transfer),4,entries,1,&conservation,finalize};
    ChemModelSystem built={0};
    assert(chem_model_system_init(&built,&base,NULL,NULL,NULL,1,&model,0,CHEM_SPARSE_CSR)==CHEM_MODEL_OK);
    assert(built.state.algebraic[2] && chem_sparse_matrix_entry(&built.structure,2,1)!=SIZE_MAX);
    double y[]={.9,1,.1};
    for (size_t step=0;step<10;++step) { solve(&built,y,1,10); assert(fabs(y[0]+y[2]-1)<1e-4); }
    chem_model_system_destroy(&built);
    const char *const names[]={"A_GAS","B"}; base.species_names=names;
    const ChemConstraintTerm r[]={{1,1}},p[]={{2,1}};
    ChemConstraint eq={.name="B_AQ_eq",.type=CHEM_CONSTRAINT_EQUILIBRIUM,.row=2,.reactants=r,.products=p,
        .reactant_count=1,.product_count=1,.constant=5,.gas_constant=8.31446261815324,.reference_temperature=298.15};
    ChemConstraintSet builtin={3,1,&eq};
    const ChemSparseEntry simple[]={{0,0},{2,0}};
    transfer=(Transfer){3,0,2,0,.1,0,0};
    model=(ChemExternalModel){{1,0,aerosol,NULL,ext_atol},component(&transfer),2,simple,0,NULL,finalize};
    assert(chem_model_system_init(&built,&base,NULL,NULL,&builtin,1,&model,0,CHEM_SPARSE_CSR)==CHEM_MODEL_OK);
    assert(built.state.constraint_count==1 && built.state.algebraic[2]);
    assert(chem_sparse_matrix_entry(&built.structure,2,0)!=SIZE_MAX);
    chem_model_system_destroy(&built);
}
static void model_order_and_reuse(void)
{
    const char *const names[]={"A","B","C","D"}; const double atol[]={.001,.001,.001,.001};
    ChemStateConfig base={1,3,0,0,names,NULL,1e-6,atol};
    const size_t dependencies[]={1,2}; Equilibrium eq={1,2,5};
    ChemCustomConstraint custom={&eq,dependencies,2,eq_residual,eq_jac};
    ChemConstraint external={.name="B_C",.type=CHEM_CONSTRAINT_CUSTOM,.row=2,.custom=&custom};
    ChemExternalModel models[2]={0}; models[0].constraint_count=1; models[0].constraints=&external;
    Transfer drive={3,0,1,0,.1,0,0}; ChemOdeComponent process=component(&drive);
    unsigned char base_pattern[9]={1,0,0,1,0,0,0,0,0};
    double results[2][3];
    for (size_t reverse=0;reverse<2;++reverse) {
        if (reverse) { ChemExternalModel swap=models[0]; models[0]=models[1]; models[1]=swap; }
        ChemModelSystem built={0};
        assert(chem_model_system_init(&built,&base,&process,base_pattern,NULL,2,models,0,CHEM_SPARSE_CSR)==CHEM_MODEL_OK);
        assert(built.state.constraint_count==1 && built.state.algebraic[2]);
        double y[]={1,0,0}; solve(&built,y,20,1);
        assert(y[1]>0 && fabs(5*y[1]-y[2])<1e-6);
        memcpy(results[reverse],y,sizeof(y)); chem_model_system_destroy(&built);
    }
    assert(results[0][1]==results[1][1] && results[0][2]==results[1][2]);

    base.species_count=4; drive.n=4;
    const ChemConstraintTerm ar[]={{0,1}},dp[]={{3,1}};
    ChemConstraint built_in={.name="A_D_eq",.type=CHEM_CONSTRAINT_EQUILIBRIUM,.row=3,.reactants=ar,.products=dp,
        .reactant_count=1,.product_count=1,.constant=2,.gas_constant=8.31446261815324,.reference_temperature=298.15};
    ChemConstraintSet set={4,1,&built_in};
    unsigned char pattern4[16]={1,0,0,0,1,0,0,0};
    ChemExternalModel single={0}; single.constraint_count=1; single.constraints=&external;
    ChemModelSystem first={0},second={0};
    assert(chem_model_system_init(&first,&base,&process,pattern4,&set,1,&single,0,CHEM_SPARSE_CSR)==CHEM_MODEL_OK);
    assert(chem_model_system_init(&second,&base,&process,pattern4,&set,1,&single,0,CHEM_SPARSE_CSR)==CHEM_MODEL_OK);
    assert(first.state.constraint_count==2 && second.state.constraint_count==2);
    double a[]={1,0,0,2},b[]={1,0,0,2}; solve(&first,a,20,1); solve(&second,b,20,1);
    assert(fabs(2*a[0]-a[3])<1e-6 && fabs(5*a[1]-a[2])<1e-6);
    for (size_t i=0;i<4;++i) assert(a[i]==b[i]);
    chem_model_system_destroy(&first); chem_model_system_destroy(&second);
}
int main(void) { metadata_and_patterns(); algebraic_pattern(); model_order_and_reuse(); return 0; }
