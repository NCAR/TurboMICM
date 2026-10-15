/* SPDX-License-Identifier: Apache-2.0 */
#include "test.h"
#include "chem_solver.h"
int test_failures,total_failures;
ChemStatus extra_create(ChemMechanism **);
int main(void)
{
    ChemMechanism *m=NULL;ChemState *s=NULL;CHECK(extra_create(&m)==CHEM_OK);
    CHECK(chem_state_create(m,3,&s)==CHEM_OK);
    for (size_t i=0;i<3;++i) {
        double t=270+20*i,p=80000+9000*i;s->conditions[i]=(ChemConditions){t,p,0};chem_ideal_air_density(s->conditions+i);
    }
    chem_update_rates(s);
    for (size_t i=0;i<3;++i) {
        double t=s->conditions[i].temperature,p=s->conditions[i].pressure,density=s->conditions[i].air_density;
        double k0=.013*exp(2/t)*pow(t/290,1.2),ki=.7*exp(8/t)*pow(t/280,.3),km=.002*exp(9/t)*pow(t/270,.1);
        double ratio=k0*density/ki,kf=k0*density/(1+ratio)*pow(.8,1.1/(1.1+pow(log10(ratio),2)));
        REL(s->rates[i],kf+km*(1-kf/ki),1e-14);
        REL(s->rates[3+i],(1+.001*t+.000001*t*t)*.004*exp(14/t)*pow(t/280,.7)*(1+.00001*p),1e-14);
        REL(s->rates[6+i],2*exp(-30/t)*.8,1e-14);
    }
    chem_state_free(s);chem_mechanism_free(m);
    printf("Additional generated rate types: %s\n",test_failures?"FAIL":"PASS");return test_failures?1:0;
}
