#include "chem.h"
#include "jacobian_verify.h"
#include "matrix.h"
#include "sparse.h"
#include "state.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

enum { MAX_N=4 };
typedef struct {
    size_t n,source,target,row,solvent;int kinetic,has_solvent,mode,temperature;
    double rate,total,ratio,updated_ratio;
    unsigned char kinetic_pattern[MAX_N*MAX_N],rows[MAX_N],constraint_pattern[MAX_N*MAX_N];
} External;
static double ratio(const External *c,ChemConditions conditions)
{ (void)conditions; return c->temperature ? c->updated_ratio : c->ratio; }
static int process_f(void *opaque,ChemConditions cond,const double *y,double *f)
{
    (void)cond;External *c=opaque;double v=c->rate*y[c->source]*(c->has_solvent?y[c->solvent]:1);
    f[c->source]-=v;f[c->target]+=v;return 0;
}
static int process_j(void *opaque,ChemConditions cond,const double *y,double *j)
{
    (void)cond;External *c=opaque;double v=c->rate*(c->has_solvent?y[c->solvent]:1);
    j[c->source*c->n+c->source]-=v;j[c->target*c->n+c->source]+=v;
    if(c->has_solvent) {v=c->rate*y[c->source];j[c->source*c->n+c->solvent]-=v;j[c->target*c->n+c->solvent]+=v;}return 0;
}
static int constraint_f(void *opaque,ChemConditions cond,const double *y,double *f)
{
    External *c=opaque;
    if(c->mode==1) f[c->row]=y[c->source]+y[c->target]-c->total;
    else if(c->mode==2) f[c->row]=ratio(c,cond)*y[c->source]-y[c->target];
    else if(c->mode==3) {f[c->row]=-c->total;for(size_t i=0;i<c->n;++i) f[c->row]+=y[i];}
    return 0;
}
static int constraint_j(void *opaque,ChemConditions cond,const double *y,double *j)
{
    (void)y;External *c=opaque;
    if(c->mode==1) {j[c->row*c->n+c->source]+=1;j[c->row*c->n+c->target]+=1;}
    else if(c->mode==2) {j[c->row*c->n+c->source]+=ratio(c,cond);j[c->row*c->n+c->target]-=1;}
    else if(c->mode==3) for(size_t i=0;i<c->n;++i) j[c->row*c->n+i]+=1;
    return 0;
}
static int update_parameter(void *opaque,ChemConditions cond,double *p,size_t count)
{ External *c=opaque;if(!count || !(cond.temperature>0))return 1;p[0]=c->ratio*exp(3000*(1/298.15-1/cond.temperature));c->updated_ratio=p[0];return 0; }
static ChemExternalModel make_model(External *c)
{
    ChemExternalModel model={0};
    if(c->kinetic) {
        c->kinetic_pattern[c->source*c->n+c->source]=c->kinetic_pattern[c->target*c->n+c->source]=1;
        if(c->has_solvent) c->kinetic_pattern[c->source*c->n+c->solvent]=c->kinetic_pattern[c->target*c->n+c->solvent]=1;
        model.jacobian_pattern=c->kinetic_pattern;model.add_forcing=process_f;model.add_jacobian=process_j;
    }
    if(c->mode) {
        c->rows[c->row]=1;
        if(c->mode==3) for(size_t j=0;j<c->n;++j)c->constraint_pattern[c->row*c->n+j]=1;
        else c->constraint_pattern[c->row*c->n+c->source]=c->constraint_pattern[c->row*c->n+c->target]=1;
        model.constraint_rows=c->rows;model.constraint_pattern=c->constraint_pattern;model.constraint_residual=constraint_f;model.constraint_jacobian=constraint_j;
    }
    model.context=c;if(c->temperature)model.update_parameters=update_parameter;return model;
}
static ChemMechanism empty_mechanism(size_t n,const char *const *names)
{ ChemMechanism m={.species_count=n,.species_names=names};return m; }
static int basic(int kind)
{
    const char *names[]={"A_GAS","AEROSOL.A_AQ"};ChemMechanism m=empty_mechanism(2,names);
    int mode=kind==1 || kind==5?0:1;External c={.n=2,.source=0,.target=1,.row=1,.kinetic=1,.mode=mode,.rate=kind==3?10:.1,.total=1};
    ChemExternalModel model=make_model(&c);ChemWorkspace *w=chem_workspace_create_external(&m,NULL,0,&model,1);if(!w)return 0;
    unsigned char mass[2];int ok=!chem_workspace_mass_diagonal(w,mass,2) && mass[0]==1 && mass[1]==(mode?0:1);
    ChemConditions conditions={298,101325,0,0};double dummy[]={0},y[]={.9,.1};ChemRosenbrockParameters p=chem_rosenbrock_four_stage_dae();
    if(kind==3) {
        double f[2],j[4];y[0]=.7;y[1]=.2;
        if(chem_workspace_evaluate(w,&m,conditions,dummy,y,f,j)) ok=0;
        if(fabs(f[0]+7)>1e-12 || fabs(f[1]+.1)>1e-12 || fabs(j[0]-10)>1e-12 || fabs(j[2]+1)>1e-12 || fabs(j[3]+1)>1e-12) ok=0;
    } else if(kind==2 || kind==6) {
        if(kind==6){y[0]=.8;y[1]=.2;}
        for(int step=0;step<(kind==2?20:1);++step) {ChemStatistics s;if(chem_rosenbrock_solve(&m,w,&p,conditions,dummy,kind==2?10:50,y,&s) || fabs(y[0]+y[1]-1)>1e-4)ok=0;}
    } else if(kind==7) {
        ChemWorkspace *work[]={w,w,w}; /* allocate separate workspaces below */
        ChemWorkspace *w2=chem_workspace_create_external(&m,NULL,0,&model,1),*w3=chem_workspace_create_external(&m,NULL,0,&model,1);
        work[1]=w2;work[2]=w3;ChemConditions cs[]={conditions,conditions,conditions};double states[]={.9,.1,.7,.3,.5,.5};ChemStatistics s;
        if(!w2 || !w3 || chem_rosenbrock_solve_group(&m,work,3,&p,cs,dummy,50,states,&s))ok=0;
        for(int i=0;i<3;++i) if(fabs(states[2*i]+states[2*i+1]-1)>1e-4)ok=0;
        chem_workspace_free(w2);chem_workspace_free(w3);
    }
    chem_workspace_free(w);return ok;
}
static int combined(void)
{
    const char *names[]={"A_GAS","B","C","AEROSOL.A_AQ"};ChemTerm rs[]={{0,1}},ps[]={{1,1}},er[]={{1,1}},ep[]={{2,1}};
    ChemReaction rx={.type=CHEM_ARRHENIUS,.reactant_count=1,.product_count=1,.p={.05,0,0,300}};
    ChemMechanism m={.species_count=4,.reaction_count=1,.species_names=names,.reactions=&rx,.reactants=rs,.products=ps};
    ChemConstraint constraint={"B_C_eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,er,ep,1,1,5,0,8.31446261815324,298.15};
    External c={.n=4,.source=0,.target=3,.row=3,.kinetic=1,.mode=1,.rate=.01,.total=1};ChemExternalModel model=make_model(&c);
    ChemWorkspace *w=chem_workspace_create_external(&m,&constraint,1,&model,1);if(!w)return 0;
    unsigned char mass[4];int ok=!chem_workspace_mass_diagonal(w,mass,4) && mass[0] && mass[1] && !mass[2] && !mass[3];
    chem_workspace_free(w);return ok;
}
static int surviving_pattern(int builtin)
{
    const char *names[]={"A_GAS","S","AEROSOL.A_AQ"};ChemMechanism m=empty_mechanism(3,names);
    External c={.n=3,.source=0,.target=2,.row=2,.solvent=1,.kinetic=1,.has_solvent=!builtin,.mode=builtin?0:1,.rate=.1,.total=1};ChemExternalModel model=make_model(&c);
    ChemTerm er[]={{1,1}},ep[]={{2,1}};ChemConstraint constraint={"eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,er,ep,1,1,5,0,8.31446261815324,298.15};
    ChemWorkspace *w=chem_workspace_create_external(&m,builtin?&constraint:NULL,builtin?1:0,&model,1);if(!w)return 0;
    unsigned char mass[3],pattern[9];int ok=!chem_workspace_mass_diagonal(w,mass,3) && !mass[2] && !chem_workspace_jacobian_pattern(w,pattern,9) && pattern[2*3+0];
    if(!builtin && !pattern[2*3+1]) ok=0;
    if(!builtin) {double y[]={.9,1,.1},dummy[]={0};ChemConditions cond={298,101325,0,0};ChemRosenbrockParameters p=chem_rosenbrock_four_stage_dae();
        for(int i=0;i<10;++i) {ChemStatistics stats;if(chem_rosenbrock_solve(&m,w,&p,cond,dummy,10,y,&stats) || fabs(y[0]+y[2]-1)>1e-4)ok=0;}}
    chem_workspace_free(w);return ok;
}
typedef struct {External *model;int constraint;} FDContext;
static void fd_force(size_t cell,const double *y,double *f,void *opaque)
{ (void)cell;FDContext *c=opaque;ChemConditions cond={300,101325,0,0};if(c->constraint)constraint_f(c->model,cond,y,f);else process_f(c->model,cond,y,f); }
static int finite_difference(int kind)
{
    External c={.n=2,.source=0,.target=1,.row=1,.kinetic=kind==0,.mode=kind==0?0:kind==1?1:2,.rate=3.5,.total=1,.ratio=2.5};ChemExternalModel model=make_model(&c);
    size_t cells=kind==2?1:2;ChemMatrix y={0},fd={0};ChemSparseBuilder b={0};ChemSparseMatrix jac={0};int ok=1;
    if(chem_matrix_init(&y,cells,2,1,0,0)||chem_matrix_init(&fd,cells,4,1,0,0)||chem_sparse_builder_init(&b,2))return 0;
    const unsigned char *pattern=kind==0?model.jacobian_pattern:model.constraint_pattern;
    for(size_t i=0;i<2;++i) for(size_t j=0;j<2;++j)if(pattern[i*2+j])chem_sparse_builder_add(&b,i,j);
    if(chem_sparse_init(&jac,&b,cells,1,0,0))ok=0;
    chem_sparse_builder_free(&b);
    double input[2][2]={{kind==2?3:.6,kind==2?5:.4},{.2,.8}};
    if(kind==0){input[0][0]=.8;input[0][1]=.2;input[1][0]=.3;input[1][1]=.7;}
    for(size_t cell=0;cell<cells;++cell) {
        double j[4]={0};chem_matrix_set_row(&y,cell,input[cell],2);ChemConditions cond={300,101325,0,0};
        if(kind==0)process_j(&c,cond,input[cell],j);else constraint_j(&c,cond,input[cell],j);
        for(size_t i=0;i<2;++i)for(size_t k=0;k<2;++k)if(pattern[i*2+k])*chem_sparse_at(&jac,cell,i,k)=-j[i*2+k];
    }
    FDContext ctx={&c,kind!=0};ChemJacobianComparison comparison,sparsity;
    if(chem_finite_difference_jacobian(&y,1e-8,fd_force,&ctx,&fd)||chem_compare_jacobian(&jac,&fd,1e-7,1e-7,&comparison)||!comparison.passed)ok=0;
    if((kind==0 || kind==2) && (chem_check_jacobian_sparsity(&jac,&fd,1e-6,&sparsity)||!sparsity.passed))ok=0;
    chem_sparse_free(&jac);chem_matrix_free(&y);chem_matrix_free(&fd);return ok;
}
typedef struct {External *first,*second;} Compound;
static int compound_f(void *opaque,ChemConditions cond,const double *y,double *f)
{Compound *c=opaque;return constraint_f(c->first,cond,y,f) || constraint_f(c->second,cond,y,f);}
static int compound_j(void *opaque,ChemConditions cond,const double *y,double *j)
{Compound *c=opaque;return constraint_j(c->first,cond,y,j) || constraint_j(c->second,cond,y,j);}
static int convergence_run(int mode,int four,double *result)
{
    const char *names[]={"A","B","C","D"};size_t n=four?4:3,nr=mode==0?(four?5:3):1;
    ChemTerm rs[]={{0,1},{1,1},{2,1},{1,1},{3,1}},ps[]={{1,1},{2,1},{1,1},{3,1},{1,1}};
    double rates[]={.1,10,2,8,8.0/3};ChemReaction reactions[5]={0};
    for(size_t i=0;i<nr;++i) reactions[i]=(ChemReaction){.type=CHEM_ARRHENIUS,.reactant_begin=i,.reactant_count=1,.product_begin=i,.product_count=1,.p={1,0,0,300}};
    ChemMechanism m={.species_count=n,.reaction_count=nr,.species_names=names,.reactions=reactions,.reactants=rs,.products=ps};
    External c[3]={{.n=n,.source=1,.target=2,.row=2,.mode=2,.ratio=5},
        {.n=n,.source=1,.target=3,.row=3,.mode=2,.ratio=3}, {.n=n,.row=1,.mode=3,.total=1}};
    ChemExternalModel models[3]={{0}};size_t count=0;unsigned char rows[MAX_N]={0},pattern[MAX_N*MAX_N]={0};Compound compound={&c[0],&c[2]};
    if(mode==1 && !four) {
        ChemExternalModel eq_model=make_model(&c[0]),mass_model=make_model(&c[2]);
        for(size_t i=0;i<n;++i) rows[i]=eq_model.constraint_rows[i] || mass_model.constraint_rows[i];
        for(size_t i=0;i<n*n;++i) pattern[i]=eq_model.constraint_pattern[i] || mass_model.constraint_pattern[i];
        models[count++]=(ChemExternalModel){.context=&compound,.constraint_rows=rows,.constraint_pattern=pattern,.constraint_residual=compound_f,.constraint_jacobian=compound_j};
    } else if(mode) {models[count++]=make_model(&c[0]);if(four) models[count++]=make_model(&c[1]);if(mode==1)models[count++]=make_model(&c[2]);}
    ChemWorkspace *w=chem_workspace_create_external(&m,NULL,0,models,count);if(!w)return 0;
    ChemRosenbrockParameters p=mode?chem_rosenbrock_four_stage_dae():chem_rosenbrock_three_stage();
    ChemConditions conditions={298,101325,0,0};double y[4]={1,0,0,0};int ok=1;
    for(int step=0;step<200;++step) {ChemStatistics stats;if(chem_rosenbrock_solve(&m,w,&p,conditions,rates,1,y,&stats)) {ok=0;break;}}
    memcpy(result,y,n*sizeof(double));chem_workspace_free(w);return ok;
}
static int convergence(int mode)
{
    int four=mode==2;double kinetic[4]={0},constrained[4]={0};
    int ok=convergence_run(0,four,kinetic) && convergence_run(mode==1?2:1,four,constrained);
    if(mode==1) {
        if(fabs(kinetic[2]/kinetic[1]-5)>1e-3 || fabs(constrained[2]/constrained[1]-5)>1e-3 ||
           fabs(kinetic[0]+kinetic[1]+kinetic[2]-1)>1e-6 || fabs(constrained[0]+constrained[1]+constrained[2]-6)>.1)ok=0;
    } else {
        double expected[]={0,1.0/(four?9:6),5.0/(four?9:6),3.0/9};
        for(size_t i=0;i<(four?4u:3u);++i) {
            double tolerance=i==0?1e-6:1e-3;
            if(fabs(kinetic[i]-expected[i])>tolerance || fabs(constrained[i]-expected[i])>tolerance || fabs(kinetic[i]-constrained[i])>tolerance)ok=0;
        }
        double mass1=0,mass2=0;for(size_t i=0;i<(four?4u:3u);++i){mass1+=kinetic[i];mass2+=constrained[i];}
        if(fabs(mass1-1)>1e-6 || fabs(mass2-1)>1e-3 || fabs(5*constrained[1]-constrained[2])>1e-6 || (four && fabs(3*constrained[1]-constrained[3])>1e-6))ok=0;
    }
    return ok;
}
static int builtin_external_comparison(void)
{
    const char *names[]={"A","B","C"};ChemTerm rs[]={{0,1}},ps[]={{1,1}},er[]={{1,1}},ep[]={{2,1}};
    ChemReaction rx={.type=CHEM_ARRHENIUS,.reactant_count=1,.product_count=1,.p={.1,0,0,300}};
    ChemMechanism m={.species_count=3,.reaction_count=1,.species_names=names,.reactions=&rx,.reactants=rs,.products=ps};
    ChemConstraint c={"eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,er,ep,1,1,5,0,8.31446261815324,298.15};
    External ctx={.n=3,.source=1,.target=2,.row=2,.mode=2,.ratio=5};ChemExternalModel model=make_model(&ctx);
    ChemWorkspace *a=chem_workspace_create_constrained(&m,&c,1),*b=chem_workspace_create_external(&m,NULL,0,&model,1);
    if(!a || !b) {chem_workspace_free(a);chem_workspace_free(b);return 0;}
    ChemRosenbrockParameters p=chem_rosenbrock_four_stage_dae();ChemConditions cond={298,101325,0,0};double rates[]={.1},x[]={1,0,0},y[]={1,0,0};int ok=1;
    for(int step=0;step<100;++step) {ChemStatistics sa,sb;
        if(chem_rosenbrock_solve(&m,a,&p,cond,rates,.5,x,&sa)||chem_rosenbrock_solve(&m,b,&p,cond,rates,.5,y,&sb)) {ok=0;break;}
        for(size_t i=0;i<3;++i)if(fabs(x[i]-y[i])>1e-8)ok=0;
        if(fabs(5*x[1]-x[2])>1e-6 || fabs(5*y[1]-y[2])>1e-6)ok=0;
    }
    chem_workspace_free(a);chem_workspace_free(b);return ok;
}
static int no_op(void *context,ChemConditions conditions,const double *y,double *out)
{(void)context;(void)conditions;(void)y;(void)out;return 0;}
static int model_order(void)
{
    const char *names[]={"A","B","C"};ChemTerm rs[]={{0,1}},ps[]={{1,1}};ChemReaction rx={.type=CHEM_ARRHENIUS,.reactant_count=1,.product_count=1};
    ChemMechanism m={.species_count=3,.reaction_count=1,.species_names=names,.reactions=&rx,.reactants=rs,.products=ps};
    External ctx={.n=3,.source=1,.target=2,.row=2,.mode=2,.ratio=5};ChemExternalModel eq_model=make_model(&ctx);
    unsigned char zero[9]={0};ChemExternalModel noop={.jacobian_pattern=zero,.add_forcing=no_op,.add_jacobian=no_op};
    ChemExternalModel first[]={eq_model,noop},second[]={noop,eq_model};ChemWorkspace *a=chem_workspace_create_external(&m,NULL,0,first,2),*b=chem_workspace_create_external(&m,NULL,0,second,2);
    if(!a || !b){chem_workspace_free(a);chem_workspace_free(b);return 0;}
    ChemRosenbrockParameters p=chem_rosenbrock_four_stage_dae();ChemConditions cond={298,101325,0,0};double rates[]={.1},x[]={1,0,0},y[]={1,0,0};int ok=1;
    for(int i=0;i<20;++i){ChemStatistics sa,sb;if(chem_rosenbrock_solve(&m,a,&p,cond,rates,1,x,&sa)||chem_rosenbrock_solve(&m,b,&p,cond,rates,1,y,&sb)){ok=0;break;}}
    if(!(x[1]>0) || fabs(5*x[1]-x[2])>1e-6 || fabs(5*y[1]-y[2])>1e-6 || x[1]!=y[1] || x[2]!=y[2])ok=0;
    chem_workspace_free(a);chem_workspace_free(b);return ok;
}
static int multiple_workspaces(void)
{
    const char *names[]={"A","B","C","D"};ChemTerm rs[]={{0,1}},ps[]={{1,1}},er[]={{0,1}},ep[]={{3,1}};
    ChemReaction rx={.type=CHEM_ARRHENIUS,.reactant_count=1,.product_count=1};ChemMechanism m={.species_count=4,.reaction_count=1,.species_names=names,.reactions=&rx,.reactants=rs,.products=ps};
    ChemConstraint c={"eq",CHEM_EQUILIBRIUM_CONSTRAINT,3,er,ep,1,1,2,0,8.31446261815324,298.15};
    External ctx={.n=4,.source=1,.target=2,.row=2,.mode=2,.ratio=5};ChemExternalModel model=make_model(&ctx);
    ChemWorkspace *a=chem_workspace_create_external(&m,&c,1,&model,1),*b=chem_workspace_create_external(&m,&c,1,&model,1);
    if(!a || !b){chem_workspace_free(a);chem_workspace_free(b);return 0;}
    ChemRosenbrockParameters p=chem_rosenbrock_four_stage_dae();ChemConditions cond={298,101325,0,0};double rates[]={.1},x[]={1,0,0,2},y[]={1,0,0,2};int ok=1;
    for(int i=0;i<20;++i){ChemStatistics sa,sb;if(chem_rosenbrock_solve(&m,a,&p,cond,rates,1,x,&sa)||chem_rosenbrock_solve(&m,b,&p,cond,rates,1,y,&sb)){ok=0;break;}}
    if(fabs(2*x[0]-x[3])>1e-6 || fabs(5*x[1]-x[2])>1e-6 || memcmp(x,y,sizeof x))ok=0;
    chem_workspace_free(a);chem_workspace_free(b);return ok;
}
static int temperature_constraint(void)
{
    const char *names[]={"A","B","C"},*labels[]={"temperature_equilibrium"};ChemTerm rs[]={{0,1}},ps[]={{1,1}};ChemReaction rx={.type=CHEM_ARRHENIUS,.reactant_count=1,.product_count=1};
    ChemMechanism m={.species_count=3,.reaction_count=1,.species_names=names,.reactions=&rx,.reactants=rs,.products=ps};
    External ctx={.n=3,.source=1,.target=2,.row=2,.mode=2,.ratio=2,.temperature=1};ChemExternalModel model=make_model(&ctx);
    ChemStateParameters sp={1,3,1,names,labels,1e-6,NULL};ChemState state={0};if(chem_state_init(&state,&sp,1,1))return 0;
    ChemWorkspace *w=chem_workspace_create_external(&m,NULL,0,&model,1);if(!w){chem_state_free(&state);return 0;}
    ChemRosenbrockParameters p=chem_rosenbrock_four_stage_dae();double rates[]={.1};int ok=1;
    for(int t=0;t<2;++t) {
        double temperature=t?350:298.15;state.conditions[0]=(ChemConditions){temperature,101325,0,0};double y[]={1,0,0};
        for(int step=0;step<200;++step) {ChemStatistics stats;
            if(chem_state_update_external_parameters(&state,&model,1) || chem_rosenbrock_solve(&m,w,&p,state.conditions[0],rates,1,y,&stats)){ok=0;break;}}
        double expected=2*exp(3000*(1/298.15-1/temperature));
        if(!(y[1]>0) || (t && !(expected>2)) || fabs(y[2]/y[1]-expected)>1e-4 || *chem_matrix_at(&state.custom_parameters,0,0)!=expected)ok=0;
    }
    chem_workspace_free(w);chem_state_free(&state);return ok;
}

int main(void)
{
    const char *names[]={"AddExternalModelWithConstraints","AddExternalModelProcessOnly","DAESolveEnforcesConservation","AlgebraicRowHoldsOnlyConstraintTerms","CombinedBuiltInAndExternalConstraints","AddExternalModelOnlyStandardRosenbrock","AddExternalModelConstraintsOnly","MultiGridCell"};
    for(int i=0;i<8;++i)printf("external_model_constraints:ExternalModelConstraints.%s %d\n",names[i],i==4?combined():basic(i));
    printf("external_model_constraints:ExternalModelConstraints.ProcessJacobianElementInAlgebraicRowSurvivesFiltering %d\n",surviving_pattern(0));
    printf("external_model_constraints:ExternalModelConstraints.ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering %d\n",surviving_pattern(1));
    const char *fdnames[]={"ProcessForcingJacobian","ConstraintResidualJacobian","EquilibriumConstraintModelJacobian"};
    for(int i=0;i<3;++i)printf("external_model_constraints:ExternalModelFiniteDifferenceJacobian.%s %d\n",fdnames[i],finite_difference(i));
    printf("external_model_constraints:ExternalModelConstraints.KineticVsConservativeConstraintConvergence %d\n",convergence(0));
    printf("external_model_constraints:ExternalModelConstraints.SimpleConstraintPreservesRatioNotMass %d\n",convergence(1));
    printf("external_model_constraints:ExternalModelConstraints.MultiEquilibriumKineticVsComposedConstraints %d\n",convergence(2));
    printf("external_model_constraints:ExternalModelConstraints.BuiltInVsExternalModelConstraintStepByStep %d\n",builtin_external_comparison());
    printf("external_model_constraints:ExternalModelConstraints.TemperatureDependentConstraintParameter %d\n",temperature_constraint());
    printf("external_model_constraints:ExternalModelConstraints.ConstraintActiveWhenProcessOnlyModelAddedFirst %d\n",model_order());
    printf("external_model_constraints:ExternalModelConstraints.BuilderCanBuildMultipleSolvers %d\n",multiple_workspaces());
    return 0;
}
