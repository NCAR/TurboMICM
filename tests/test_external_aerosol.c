#include "chem.h"
#include "state.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

enum {N=22};
static const char *names[]={"FO2","BAR","STUB1.MODE1.QUUX.BAZ","STUB1.MODE1.QUUX.QUX","STUB1.MODE2.QUUX.BAZ","STUB1.MODE2.QUUX.QUX","STUB1.MODE2.CORGE.FO2","STUB1.MODE2.CORGE.BAZ","STUB1.MODE2.CORGE.QUX","STUB2.MODE1.NUMBER","STUB2.MODE1.QUUX.BAZ","STUB2.MODE1.QUUX.QUX","STUB2.MODE2.NUMBER","STUB2.MODE2.CORGE.FO2","STUB2.MODE2.CORGE.BAZ","STUB2.MODE2.CORGE.QUX","STUB2.MODE3.NUMBER","STUB2.MODE3.QUUX.BAZ","STUB2.MODE3.QUUX.QUX","STUB2.MODE3.CORGE.FO2","STUB2.MODE3.CORGE.BAZ","STUB2.MODE3.CORGE.QUX"};
static const char *params[]={"STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT","STUB2.PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT"};
typedef struct {int model;double rate;} Context;
static int forcing(void *opaque,ChemConditions conditions,const double *y,double *f)
{
    Context *c=opaque;
    if(c->model==1) {double a=1e-3*y[0],b=2e-3*y[2];f[0]-=a;f[6]+=a;f[2]-=b;f[4]+=b;}
    else {double a=c->rate*y[13],b=.005*conditions.temperature*y[17];f[13]-=a;f[14]+=a;f[17]-=b;f[18]+=b;}
    return 0;
}
static int jacobian(void *opaque,ChemConditions conditions,const double *y,double *j)
{
    (void)y;Context *c=opaque;
    if(c->model==1) {j[0*N+0]-=1e-3;j[6*N+0]+=1e-3;j[2*N+2]-=2e-3;j[4*N+2]+=2e-3;}
    else {j[13*N+13]-=c->rate;j[14*N+13]+=c->rate;j[17*N+17]-=.005*conditions.temperature;j[18*N+17]+=.005*conditions.temperature;}
    return 0;
}
static int update(void *opaque,ChemConditions conditions,double *parameters,size_t count)
{
    (void)opaque;if(count<2) return 1;parameters[1]=.005*conditions.temperature;return 0;
}
static void initialize(ChemExternalModel *models,Context *contexts,unsigned char patterns[2][N*N])
{
    memset(patterns,0,2*N*N);int row[][4]={{0,6,2,4},{13,14,17,18}}, col[][4]={{0,0,2,2},{13,13,17,17}};
    for(int m=0;m<2;++m) for(int i=0;i<4;++i) patterns[m][row[m][i]*N+col[m][i]]=1;
    contexts[0]=(Context){1,0};contexts[1]=(Context){2,0};
    models[0]=(ChemExternalModel){.jacobian_pattern=patterns[0],.add_forcing=forcing,.add_jacobian=jacobian,.context=&contexts[0],.variable_names=names+2,.variable_count=7};
    models[1]=(ChemExternalModel){.jacobian_pattern=patterns[1],.add_forcing=forcing,.add_jacobian=jacobian,.context=&contexts[1],.variable_names=names+9,.variable_count=13,.parameter_names=params,.parameter_count=2,.update_parameters=update};
}
static int make_state(ChemState *state,size_t cells,size_t width,ChemExternalModel *models)
{
    const char *all[N],*all_params[2];size_t n=0,np=0;
    all[n++]=names[0];all[n++]=names[1];
    for(size_t m=0;m<2;++m) {
        for(size_t j=0;j<models[m].variable_count;++j) all[n++]=models[m].variable_names[j];
        for(size_t j=0;j<models[m].parameter_count;++j) all_params[np++]=models[m].parameter_names[j];
    }
    if(n!=N || np!=2) return 1;
    /* The API borrows names, so use the persistent equivalent arrays. */
    for(size_t i=0;i<n;++i) if(strcmp(all[i],names[i])) return 1;
    for(size_t i=0;i<np;++i) if(strcmp(all_params[i],params[i])) return 1;
    ChemStateParameters p={0,N,2,names,params,1e-6,NULL};return chem_state_init(state,&p,cells,width);
}
static int metadata(size_t width)
{
    ChemExternalModel models[2];Context context[2];unsigned char patterns[2][N*N];initialize(models,context,patterns);
    ChemState s={0};if(make_state(&s,1,width,models)) return 0;
    int ok=s.parameters.variable_count==2+models[0].variable_count+models[1].variable_count && models[0].parameter_count==0 && models[1].parameter_count==2;
    for(size_t i=0;i<N;++i) if(chem_state_variable_index(&s,names[i])<0) ok=0;
    if(chem_state_variable_index(&s,"STUB1.MODE1.QUUX.QUX")<0 || chem_state_variable_index(&s,"STUB1.MODE2.CORGE.FO2")<0 ||
       chem_state_variable_index(&s,"STUB2.MODE1.NUMBER")<0 || chem_state_variable_index(&s,"STUB2.MODE2.CORGE.BAZ")<0) ok=0;
    chem_state_free(&s);return ok;
}
static int state_updates(size_t cells,size_t width)
{
    ChemExternalModel models[2];Context context[2];unsigned char patterns[2][N*N];initialize(models,context,patterns);
    ChemState s={0};if(make_state(&s,cells,width,models)) return 0;
    const size_t columns[]={0,1,2,6,21,20,9,12}; const double start[]={1.34,2.53,.75,1.23,.42,.33,1000,500};
    const double increment[]={.01,.01,.01,.01,.01,.01,1,1}; int ok=1;
    for(size_t j=0;j<8;++j) {double values[3];for(size_t i=0;i<cells;++i) values[i]=start[j]+i*increment[j];if(chem_state_set_concentration(&s,names[columns[j]],values,cells)) ok=0;}
    for(size_t i=0;i<cells;++i) for(size_t j=0;j<8;++j) {
        double expected=start[j]+i*increment[j];
        if(*chem_state_const_variable(&s,i,names[columns[j]])!=expected || *chem_matrix_const_at(&s.variables,i,columns[j])!=expected) ok=0;
    }
    chem_state_free(&s);return ok;
}
static int derivatives(int mode)
{
    ChemExternalModel models[2];Context context[2];unsigned char patterns[2][N*N];initialize(models,context,patterns);
    double y[N]={0},out[N*N]={0};y[0]=1;y[2]=.5;ChemConditions c={275,101325,0,0};int ok;
    if(mode==0) {ok=!models[0].add_forcing(models[0].context,c,y,out);ok &= out[0]==-1e-3 && out[6]==1e-3 && out[2]==-1e-3 && out[4]==1e-3;}
    else {ok=!models[0].add_jacobian(models[0].context,c,y,out);ok &= out[0*N+0]==-1e-3 && out[6*N+0]==1e-3 && out[2*N+2]==-2e-3 && out[4*N+2]==2e-3;}
    return ok;
}
static int solve_case(int backward,int method,int algorithm,size_t width,size_t cells,int both)
{
    ChemExternalModel models[2];Context context[2];unsigned char patterns[2][N*N];initialize(models,context,patterns);
    ChemState s={0};if(make_state(&s,cells,width,models)) return 0;
    ChemMechanism m={.species_count=N,.species_names=names};
    ChemWorkspace *w=chem_workspace_create_external(&m,NULL,0,models,2);if(!w){chem_state_free(&s);return 0;}
    for(size_t cell=0;cell<cells;++cell) {
        *chem_matrix_at(&s.variables,cell,0)=1+.5*cell;*chem_matrix_at(&s.variables,cell,1)=2+.5*cell;
        *chem_matrix_at(&s.variables,cell,2)=.5+.2*cell;*chem_matrix_at(&s.variables,cell,6)=.8+.2*cell;
        *chem_matrix_at(&s.variables,cell,9)=1000+100*cell;*chem_matrix_at(&s.variables,cell,12)=500+50*cell;
        *chem_matrix_at(&s.variables,cell,20)=.2+.05*cell;*chem_matrix_at(&s.variables,cell,21)=.3+.1*cell;
        if(both) { *chem_matrix_at(&s.variables,cell,13)=.4+.2*cell;*chem_matrix_at(&s.variables,cell,14)=.1+.05*cell;
            *chem_matrix_at(&s.variables,cell,17)=.2+.1*cell;*chem_matrix_at(&s.variables,cell,18)=.3+.1*cell;
            *chem_matrix_at(&s.custom_parameters,cell,0)=.01+.005*cell; }
        s.conditions[cell]=(ChemConditions){both?275+10*cell:0,101325,0,0};
    }
    int ok=!chem_state_update_external_parameters(&s,models,2);
    double grouped[3][N];
    if(!backward && cells>1 && ok) {
        ChemExternalModel cell_models[3][2];Context cell_context[3][2];unsigned char cell_patterns[3][2][N*N];
        ChemWorkspace *group_work[3]={0};
        for(size_t cell=0;cell<cells;++cell) {
            initialize(cell_models[cell],cell_context[cell],cell_patterns[cell]);
            cell_context[cell][1].rate=*chem_matrix_at(&s.custom_parameters,cell,0);
            group_work[cell]=chem_workspace_create_external(&m,NULL,0,cell_models[cell],2);
            if(!group_work[cell] || chem_matrix_get_row(&s.variables,cell,grouped[cell],N)) ok=0;
        }
        ChemRosenbrockParameters gp=method==0?chem_rosenbrock_two_stage():method==1?chem_rosenbrock_three_stage():method==2?chem_rosenbrock_four_stage():method==3?chem_rosenbrock_four_stage_dae():chem_rosenbrock_six_stage_dae();
        gp.linear_solver=(ChemLUAlgorithm)algorithm; double dummy[]={0};ChemStatistics stats;
        if(ok && chem_rosenbrock_solve_group(&m,group_work,cells,&gp,s.conditions,dummy,10,&grouped[0][0],&stats)) ok=0;
        for(size_t cell=0;cell<cells;++cell) chem_workspace_free(group_work[cell]);
    }
    for(size_t cell=0;cell<cells && ok;++cell) {
        context[1].rate=*chem_matrix_at(&s.custom_parameters,cell,0);
        if(*chem_matrix_at(&s.custom_parameters,cell,1)!=.005*s.conditions[cell].temperature) ok=0;
        double y[N],initial[N],dummy[1]={0};ChemStatistics stats;chem_matrix_get_row(&s.variables,cell,y,N);memcpy(initial,y,sizeof y);
        ChemRosenbrockParameters rp=method==0?chem_rosenbrock_two_stage():method==1?chem_rosenbrock_three_stage():method==2?chem_rosenbrock_four_stage():method==3?chem_rosenbrock_four_stage_dae():chem_rosenbrock_six_stage_dae();
        ChemBackwardEulerParameters bp=chem_backward_euler_defaults();rp.linear_solver=bp.linear_solver=(ChemLUAlgorithm)algorithm;
        int status=0;
        if(!backward && cells>1) memcpy(y,grouped[cell],sizeof y);
        else status=backward?chem_backward_euler_solve(&m,w,&bp,s.conditions[cell],dummy,10,y,&stats):chem_rosenbrock_solve(&m,w,&rp,s.conditions[cell],dummy,10,y,&stats);
        if(status) {ok=0;break;}
        double tolerance=backward?(both?.1:2e-4):both && method==0?5e-4:5e-5;
        double delta1=initial[0]*(1-exp(-.001*10)),delta2=initial[2]*(1-exp(-.002*10));
        if(fabs(y[0]-(initial[0]-delta1))>tolerance*initial[0] || y[1]!=initial[1] ||
           fabs(y[2]-(initial[2]-delta2))>tolerance*initial[2] || fabs(y[4]-delta2)>tolerance*initial[2] ||
           fabs(y[6]-(initial[6]+delta1))>tolerance*initial[0] || y[9]!=initial[9] || y[12]!=initial[12]) ok=0;
        if(both) {
            double delta3=initial[13]*(1-exp(-context[1].rate*10)),delta4=initial[17]*(1-exp(-.005*s.conditions[cell].temperature*10));
            if(fabs(y[13]-(initial[13]-delta3))>tolerance*initial[13] || fabs(y[14]-(initial[14]+delta3))>tolerance*initial[13] ||
               fabs(y[17]-(initial[17]-delta4))>tolerance*initial[17] || fabs(y[18]-(initial[18]+delta4))>tolerance*initial[17]) {
                fprintf(stderr,"external solve backward=%d method=%d alg=%d width=%zu cell=%zu errors=%.3g %.3g %.3g %.3g tol=%.3g\n",backward,method,algorithm,width,cell,
                    fabs(y[13]-(initial[13]-delta3)),fabs(y[14]-(initial[14]+delta3)),fabs(y[17]-(initial[17]-delta4)),fabs(y[18]-(initial[18]+delta4)),tolerance);
                ok=0;
            }
        } else if(y[20]!=initial[20] || y[21]!=initial[21]) ok=0;
    }
    chem_workspace_free(w);chem_state_free(&s);return ok;
}
int main(void)
{
    const char *names_out[]={"StateIncludesStubAerosolModel","CanUpdateStateWithStubAerosolModel","CanUpdateMultiCellStateWithStubAerosolModel","CanCalculateSingleGridCellForcingWithStubAerosolModel","CanCalculateSingleGridCellJacobianWithStubAerosolModel","CanSolveSingleGridCellWithStubAerosolModel1","CanSolveSingleGridCellWithTwoStubAerosolModels","CanSolveMultiGridCellWithStubAerosolModel1","CanSolveMultiGridCellWithTwoStubAerosolModels"};
    for(int backward=0;backward<2;++backward) for(int kind=0;kind<9;++kind) {
        int ok=1;
        for(size_t width=1;width<=4 && ok;++width) {
            if(kind==0) ok &= metadata(width);
            else if(kind==1 || kind==2) ok &= state_updates(kind==1?1:3,width);
            else if(kind==3 || kind==4) ok &= derivatives(kind-3);
            else for(int method=0;method<(backward?1:5) && ok;++method) for(int alg=0;alg<2 && ok;++alg) ok &= solve_case(backward,method,alg,width,kind>=7?3:1,kind%2==0);
        }
        printf("%s:AerosolModelIntegration.%s %d\n",backward?"aerosol_model_backward_euler":"aerosol_model_rosenbrock",names_out[kind],ok);
    }
    return 0;
}
