#include "chem.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>

static ChemRosenbrockParameters method_parameters(int method,int algorithm)
{
    ChemRosenbrockParameters p=method==0?chem_rosenbrock_two_stage():method==1?chem_rosenbrock_three_stage():method==2?chem_rosenbrock_four_stage():method==3?chem_rosenbrock_four_stage_dae():chem_rosenbrock_six_stage_dae();
    p.linear_solver=(ChemLUAlgorithm)algorithm;return p;
}
static int advance(const ChemMechanism *m,ChemWorkspace *w,const ChemRosenbrockParameters *rp,const ChemBackwardEulerParameters *bp,
                   const double *rates,double dt,double *y,int backward,ChemMatrix *view)
{
    ChemConditions c={298,101325,0,0};double elapsed=0;
    while(elapsed<dt){ChemStatistics stats;int status;
        if(view) {
            ChemWorkspace *group[]={w};
            if(chem_matrix_set_row(view,0,y,m->species_count) || chem_matrix_get_row(view,0,y,m->species_count))return 0;
            status=backward?chem_backward_euler_solve_group(m,group,1,bp,&c,rates,dt-elapsed,y,&stats):
                            chem_rosenbrock_solve_group(m,group,1,rp,&c,rates,dt-elapsed,y,&stats);
            if(chem_matrix_set_row(view,0,y,m->species_count) || chem_matrix_get_row(view,0,y,m->species_count))return 0;
        } else status=backward?chem_backward_euler_solve(m,w,bp,c,rates,dt-elapsed,y,&stats):chem_rosenbrock_solve(m,w,rp,c,rates,dt-elapsed,y,&stats);
        if((status && status!=2) || !(stats.final_time>0) || elapsed+stats.final_time<=elapsed)return 0;
        elapsed+=stats.final_time;}
    return 1;
}
static int hires(int backward,int method,int algorithm,size_t width)
{
    const char *names[]={"Y0","Y1","Y2","Y3","Y4","Y5","Y6","Y7"},*params[]={"r1","r2","r3","r4","r5","r6","r7","r8","r9"};
    const ChemTerm r[]={{0,1},{1,1},{2,1},{3,1},{4,1},{5,1},{6,1},{5,1},{7,1}};
    const ChemTerm q[]={{1,1},{0,.43/8.75},{3,8.32/8.75},{0,8.32/10.03},{3,1.71/10.03},{0,1},{2,.43/1.12},{5,.69/1.12},{2,.035/1.745},{5,1.71/1.745},{4,1},{4,.43/1.81},{5,.69/1.81},{7,1},{6,1}};
    const ChemReaction rx[]={
        {CHEM_USER,0,1,0,1,0,{1}},{CHEM_USER,1,1,1,2,1,{1}},{CHEM_USER,2,1,3,2,2,{1}},
        {CHEM_USER,3,0,5,1,3,{1}},{CHEM_USER,3,1,6,2,4,{1}},{CHEM_USER,4,1,8,2,5,{1}},
        {CHEM_USER,5,1,10,1,6,{1}},{CHEM_USER,6,1,11,3,7,{1}},{CHEM_USER,7,2,14,1,8,{1}}
    };
    const ChemMechanism m={8,9,9,names,params,NULL,NULL,rx,r,q,0,NULL};
    const double expected[2][8]={
        {.000737131257332567,.000144248572631618,.000058887297409676,.001175651343283149,.002386356198831330,.006238968252742796,.002849998395185769,.002850001604814231},
        {.000670305503581864,.000130996846986347,.000046862231597733,.001044668020551705,.000594883830951485,.001399628833942774,.001014492757718480,.004685507242281520}
    };
    double rates[]={1.71,8.75,10.03,.0007,1.12,1.745,.43,1.81,280},y[]={1,0,0,0,0,0,0,.0057};
    ChemWorkspace *w=chem_workspace_create(&m);if(!w)return 0;
    ChemMatrix matrix={0};if(width && chem_matrix_init(&matrix,1,8,width,1,0)){chem_workspace_free(w);return 0;}
    ChemRosenbrockParameters p=method_parameters(method,algorithm);p.relative_tolerance=1e-6;p.absolute_tolerance=1e-8;
    ChemBackwardEulerParameters b=chem_backward_euler_defaults();b.relative_tolerance=1e-6;b.absolute_tolerance=1e-8;b.linear_solver=(ChemLUAlgorithm)algorithm;
    double tolerance=backward?.1:(method==1 || method==2?1e-7:1e-6),dt=321.8122;int ok=1;
    for(int i=0;i<2 && ok;++i){if(!advance(&m,w,&p,&b,rates,dt,y,backward,width?&matrix:NULL))ok=0;
        for(int j=0;j<8;++j)if(!isfinite(y[j]) || fabs(y[j]-expected[i][j])>tolerance){fprintf(stderr,"HIRES method %d backward %d step %d species %d got %.9g expected %.9g\n",method,backward,i,j,y[j],expected[i][j]);ok=0;}
        dt+=100;}
    chem_matrix_free(&matrix);chem_workspace_free(w);return ok;
}
static int e5(int method,int algorithm,size_t width)
{
    const char *names[]={"A1","A2","A3","A4","A5","A6"},*params[]={"r1","r2","r3","r4"};
    const ChemTerm r[]={{0,1},{1,1},{2,1},{0,1},{2,1},{3,1}},q[]={{1,1},{2,1},{4,1},{3,1},{2,1},{5,1}};
    const ChemReaction rx[]={{CHEM_USER,0,1,0,2,0,{1}},{CHEM_USER,1,2,2,1,1,{1}},{CHEM_USER,3,2,3,1,2,{1}},{CHEM_USER,5,1,4,2,3,{1}}};
    const ChemMechanism m={6,4,4,names,params,NULL,NULL,rx,r,q,0,NULL};
    const double expected[7][4]={
        {1.7599259497677897e-3,1.3846281519376516e-11,7.637003853007391e-13,1.3082581134075777e-11},
        {1.6180769999072943e-3,1.3822370304983735e-10,8.251573500683834e-12,1.2997212954915352e-10},
        {7.481320822429222e-6,2.3734781561205975e-12,2.2123586689581664e-12,1.6111948716243114e-13},
        {4.715033363040163e-10,1.8188895860807022e-14,1.8188812376786725e-14,8.34840202963217e-20},
        {3.1317148329357e-14,1.4840957952870064e-16,1.4840957948345691e-16,4.5243728279782625e-26},
        {3.813903518978709e-49,1.0192582567660293e-20,1.0192582567660293e-20,3.784493550748622e-65},
        {0,8.861233497626378e-23,8.861233497626378e-23,0}
    };
    double rates[]={7.89e-10,1.13e9,1.1e7,1.13e3},y[]={1.76e-3,0,0,0,0,0},atol[]={1e-7,1e-17,1e-17,1e-17,1e-7,1e-7};
    ChemWorkspace *w=chem_workspace_create(&m);if(!w)return 0;
    ChemMatrix matrix={0};if(width && chem_matrix_init(&matrix,1,6,width,1,0)){chem_workspace_free(w);return 0;}
    ChemRosenbrockParameters p=method_parameters(method,algorithm);p.relative_tolerance=1e-13;p.absolute_tolerances=atol;
    ChemBackwardEulerParameters b=chem_backward_euler_defaults();double target=10,now=0;int ok=1;
    for(int i=0;i<7 && ok;++i){if(!advance(&m,w,&p,&b,rates,target-now,y,0,width?&matrix:NULL))ok=0;
        for(int j=0;j<4;++j){double error=fabs(y[j]-expected[i][j]),absolute=j?1e-13:1e-6;if(!isfinite(y[j]) || (error>=absolute && error>=1e-10*fabs(expected[i][j]))){fprintf(stderr,"E5 method %d step %d species %d got %.9g expected %.9g\n",method,i,j,y[j],expected[i][j]);ok=0;}}
        now=target;target*=100;}
    chem_matrix_free(&matrix);chem_workspace_free(w);return ok;
}
static int oregonator(int backward,int method,int algorithm,size_t width)
{
    const char *names[]={"X","Y","Z","P","Q"},*params[]={"r1","r2","r3","r4","r5"};
    const ChemTerm r[]={{1,1},{0,1},{1,1},{0,1},{0,1},{0,1},{2,1}},q[]={{0,1},{3,1},{2,1},{0,2},{4,1},{1,1}};
    const ChemReaction rx[]={{CHEM_USER,0,1,0,1,0,{1}},{CHEM_USER,1,2,1,1,1,{1}},{CHEM_USER,3,1,2,2,2,{1}},{CHEM_USER,4,2,4,1,3,{1}},{CHEM_USER,6,1,5,1,4,{1}}};
    const ChemMechanism m={5,5,5,names,params,NULL,NULL,rx,r,q,0,NULL};
    const double expected[12][3]={
        {1.000661467180497,1512.778937348249,10358.54312767229},{1.000874625199626,1144.336972384497,83.72149966624639},
        {1.001890368438751,529.9926232295553,1.662279579042420},{1.004118022612645,243.8326079910346,1.008822224048647},
        {1.008995416634061,112.1664388662539,1.007783229065319},{1.019763472537298,51.59761322947535,1.016985778956374},
        {1.043985088527474,23.73442027531524,1.037691843544522},{1.100849071667922,10.91533805469020,1.085831969810860},
        {1.249102130020572,5.013945178605446,1.208326626237875},{1.779724751937019,2.281852385542403,1.613754023671725},
        {1.000889326903503,1125.438585746596,16410.49483777168},{1.000814870318523,1228.178521549889,132.0554942846513}
    };
    const double s=77.27,qconst=8.375e-6,wconst=.161,k2=1.6e9;
    const double alpha=1/(s*wconst*k2),eta=s/(wconst*k2),rho=1/(wconst*wconst*k2);
    double rates[]={1/(s*wconst),k2,s/wconst,s*s*qconst*k2/2,1},y[]={alpha,2*eta,3*rho,0,0};
    double atol[]={alpha*1e-11,eta*1e-11,rho*1e-11,eta,eta};
    ChemWorkspace *work=chem_workspace_create(&m);if(!work)return 0;
    ChemMatrix matrix={0};if(width && chem_matrix_init(&matrix,1,5,width,1,0)){chem_workspace_free(work);return 0;}
    ChemRosenbrockParameters p=method_parameters(method,algorithm);p.relative_tolerance=1e-9;p.absolute_tolerances=atol;
    ChemBackwardEulerParameters b=chem_backward_euler_defaults();b.relative_tolerance=1e-9;b.absolute_tolerances=atol;b.linear_solver=(ChemLUAlgorithm)algorithm;
    const double dt=30*wconst,scales[]={alpha,eta,rho},tolerance=backward?5e-3:1e-6;int ok=1,substeps=backward?18000:1;
    for(int i=0;i<12 && ok;++i){double elapsed=0;
        while(elapsed<dt){double step=fmin(dt/substeps,dt-elapsed);if(!advance(&m,work,&p,&b,rates,step,y,backward,width?&matrix:NULL)){ok=0;break;}if(elapsed+step<=elapsed){ok=0;break;}elapsed+=step;}
        for(int j=0;j<3;++j){double ref=expected[i][j]*scales[j];if(!isfinite(y[j]) || fabs(y[j]-ref)>=tolerance*fabs(ref)){fprintf(stderr,"Oregonator method %d backward %d step %d species %d got %.9g expected %.9g\n",method,backward,i,j,y[j],ref);ok=0;}}
    }
    chem_matrix_free(&matrix);chem_workspace_free(work);return ok;
}
int main(void)
{
    int hires_r=1,hires_b=1,e5_r=1,oregon_r=1,oregon_b=1;
    for(int algorithm=0;algorithm<2;++algorithm)for(size_t width=0;width<=4;++width){
        hires_b &= hires(1,0,algorithm,width);oregon_b &= oregonator(1,0,algorithm,width);
        for(int method=0;method<5;++method){hires_r &= hires(0,method,algorithm,width);e5_r &= e5(method,algorithm,width);oregon_r &= oregonator(0,method,algorithm,width);}
    }
    printf("analytical_rosenbrock_integration:AnalyticalExamples.HIRES %d\n",hires_r);
    printf("analytical_backward_euler:AnalyticalExamples.HIRES %d\n",hires_b);
    printf("analytical_rosenbrock_integration:AnalyticalExamples.E5 %d\n",e5_r);
    printf("analytical_rosenbrock_integration:AnalyticalExamples.Oregonator %d\n",oregon_r);
    printf("analytical_backward_euler:AnalyticalExamples.Oregonator %d\n",oregon_b);
    return 0;
}
