#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <string.h>
#include "rate_store.h"

int main(void)
{
    ChemRateDescriptor d[8]={0};
    d[0].rate.type=CHEM_RATE_SURFACE; d[0].rate.p[0]=230; d[0].rate.p[1]=.025; d[0].rate.p[2]=.74;
    d[0].has_diffusion_coefficient=d[0].has_molecular_weight=1;
    d[1].rate.type=CHEM_RATE_ARRHENIUS;
    d[1].rate.p[0]=2.15e-4; d[1].rate.p[1]=1.2; d[1].rate.p[2]=110; d[1].rate.p[3]=300;
    d[2].rate.type=CHEM_RATE_USER; d[2].rate.p[0]=2;
    d[3].rate.type=CHEM_RATE_TROE;
    d[4].rate.type=CHEM_RATE_ARRHENIUS; d[4].rate.p[0]=42;
    d[5].rate.type=CHEM_RATE_TUNNELING;
    d[6].rate.type=CHEM_RATE_USER; d[6].rate.p[0]=.5;
    d[7].rate.type=CHEM_RATE_CALLBACK; d[7].rate.custom_index=9;
    ChemRateStore s={0};
    assert(chem_rate_store_init(&s,8,d)==CHEM_RATE_STORE_OK);
    const size_t offsets[]={0,2,3,3,3,3,4,4,4,6,7,8};
    assert(s.count==8 && s.custom_count==4);
    for (size_t i=0;i<12;++i) assert(s.offset[i]==offsets[i]);
    const size_t order[]={1,4,3,5,2,6,0,7};
    for (size_t i=0;i<8;++i) assert(s.entries[i].source_index==order[i]);
    assert(s.entries[0].rate.p[0]==2.15e-4 && s.entries[0].rate.p[1]==1.2 &&
           s.entries[0].rate.p[2]==110 && s.entries[0].rate.p[3]==300);
    assert(s.entries[4].rate.custom_index==0 && s.entries[4].rate.p[0]==2);
    assert(s.entries[5].rate.custom_index==1 && s.entries[5].rate.p[0]==.5);
    assert(s.entries[6].rate.custom_index==2 && s.entries[6].rate.p[0]==230 &&
           s.entries[6].rate.p[2]==.74 && s.entries[7].rate.custom_index==9);
    const double avogadro=6.02214076e23, gas=1.380649e-23*avogadro, pi=3.14159265358979323846;
    assert(fabs(s.entries[6].mean_free_speed_factor-8*gas/(pi*.025))<1e-14);
    double custom[]={10,20,1e-6,2e8};
    assert(chem_rate(&s.entries[4].rate,(ChemConditions){300,101325,1},custom)==20);
    assert(chem_rate(&s.entries[5].rate,(ChemConditions){300,101325,1},custom)==10);
    assert(isfinite(chem_rate(&s.entries[6].rate,(ChemConditions){300,101325,1},custom)));
    chem_rate_store_destroy(&s);

    /* Reference offset fixture: two Arrhenius, then Troe, tunneling,
       user-defined, and surface. */
    ChemRateDescriptor six[]={d[1],d[4],d[3],d[5],d[2],d[0]};
    const size_t expected[]={0,2,3,3,3,3,4,4,4,5,6,6};
    assert(chem_rate_store_init(&s,6,six)==CHEM_RATE_STORE_OK);
    for (size_t i=0;i<12;++i) assert(s.offset[i]==expected[i]);
    assert(s.entries[5].rate.custom_index==1 && s.custom_count==3);
    chem_rate_store_destroy(&s);
    assert(chem_rate_store_init(&s,1,d)==CHEM_RATE_STORE_OK);
    assert(s.entries[0].rate.custom_index==0 && s.entries[0].rate.p[0]==230 &&
           s.entries[0].rate.p[2]==.74 && fabs(s.entries[0].mean_free_speed_factor-8*gas/(pi*.025))<1e-14);
    chem_rate_store_destroy(&s);

    ChemRateDescriptor branched={0};
    branched.rate.type=CHEM_RATE_ALKOXY; branched.rate.p[0]=1; branched.rate.p[2]=.5; branched.rate.p[3]=3;
    assert(chem_rate_store_init(&s,1,&branched)==CHEM_RATE_STORE_OK);
    double k0=2e-22*avogadro*1e-6*exp(3), a=k0*(2.45e19/avogadro*1e6);
    double b=.43*pow(293.0/298,-8), z=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)));
    assert(fabs(s.entries[0].branch_k0-k0)<=1e-10*k0);
    assert(fabs(s.entries[0].branch_z-z)<=1e-10*fabs(z));
    chem_rate_store_destroy(&s);

    d[0].has_diffusion_coefficient=0;
    size_t bad=123; char message[160];
    ChemRateStoreStatus status=chem_rate_store_init_diagnostic(&s,8,d,&bad);
    assert(status==CHEM_RATE_STORE_MISSING_DIFFUSION && bad==0);
    assert(chem_rate_store_error_message(status,"foo",message,sizeof(message))>0);
    assert(strstr(message,"Diffusion coefficient for species 'foo' is not defined"));
    assert(!s.entries && !s.count);
    d[0].has_diffusion_coefficient=1; d[0].has_molecular_weight=0;
    status=chem_rate_store_init_diagnostic(&s,8,d,&bad);
    assert(status==CHEM_RATE_STORE_MISSING_MOLECULAR_WEIGHT && bad==0);
    assert(chem_rate_store_error_message(status,"foo",message,sizeof(message))>0);
    assert(strstr(message,"foo") && strstr(message,"molecular weight"));
    assert(!s.entries && !s.count);
    assert(chem_rate_store_init(&s,0,NULL)==CHEM_RATE_STORE_OK);
    chem_rate_store_destroy(&s);
    return 0;
}
