#include "chem.h"
#include "rates.h"
#include "backward_euler.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
extern const ChemMechanism generated_mechanism;
extern const ChemRate generated_rates[];
extern const ChemRateMultiplier generated_multipliers[];
extern const ChemParameterizedSpecies generated_multiplier_functions[];
extern const char *const generated_callback_labels[];
extern const size_t generated_callback_count;
typedef struct { const double *rates; } Context;
static double callback(ChemConditions c,void *context) { return *(double *)context*c.temperature; }
static double rebound(ChemConditions c,void *context) { (void)context; return 2e-3*c.temperature; }
static void rhs(void *context,const double *y,double *out)
{ chem_rhs(&generated_mechanism,y,((Context *)context)->rates,out); }
static void jac(void *context,const double *y,double *out)
{ chem_jacobian(&generated_mechanism,y,((Context *)context)->rates,out); }
int main(void)
{
    assert(generated_mechanism.reaction_count==2 && generated_callback_count==1);
    assert(!strcmp(generated_callback_labels[0],"lambda_rc"));
    assert(generated_rates[0].type==CHEM_RATE_ARRHENIUS && generated_rates[1].type==CHEM_RATE_CALLBACK);
    assert(generated_rates[1].custom_index==0);
    ChemConditions c={300,101253.3,1e6}; double scale=2, rates[2];
    ChemRateBinding bindings[]={{callback,&scale}};
    chem_rates_parameterized_with_callbacks(2,generated_rates,generated_multipliers,generated_multiplier_functions,bindings,1,c,NULL,rates);
    assert(fabs(rates[1]-600)<1e-13*600);
    scale=1e-3;
    chem_rates_parameterized_with_callbacks(2,generated_rates,generated_multipliers,generated_multiplier_functions,bindings,1,c,NULL,rates);
    assert(fabs(rates[1]-.3)<1e-14);
    bindings[0].function=rebound;
    chem_rates_parameterized_with_callbacks(2,generated_rates,generated_multipliers,generated_multiplier_functions,bindings,1,c,NULL,rates);
    assert(fabs(rates[1]-.6)<1e-12);
    rates[0]=0; Context context={rates};
    const double atol[]={1e-6,1e-6}; double y[]={1,0}; ChemSolverStats stats;
    ChemBackwardEulerOptions options=chem_backward_euler_defaults(atol);
    assert(chem_backward_euler_solve(2,y,1,rhs,jac,&context,&options,&stats)==CHEM_SOLVE_OK);
    assert(stats.final_time==1 && fabs(y[0]-1/1.6)<1e-14 && fabs(y[1]-.6/1.6)<1e-14);
    chem_rates_parameterized(2,generated_rates,generated_multipliers,generated_multiplier_functions,c,NULL,rates);
    assert(isnan(rates[1]));
    puts("callback rates and rebinding PASS"); return 0;
}
