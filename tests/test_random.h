#ifndef TEST_RANDOM_H
#define TEST_RANDOM_H
/* Test-only deterministic reproduction of the reference's libc++ generators:
 * minstd_rand, MT19937, two-draw canonical doubles and the polar normal method.
 * No C++ headers or runtime are part of the build. Fixed sequences make these
 * ill-conditioned reference inputs identical across C compilers/platforms. */
typedef struct {
    uint32_t state, mt[624]; size_t next;
    int use_mt, saved; double normal;
} random_state;
static inline void random_init(random_state *r,uint32_t seed,int mt)
{
    memset(r,0,sizeof(*r));r->state=seed;r->use_mt=mt;r->next=624;
    r->mt[0]=seed;
    for(size_t i=1;i<624;++i)r->mt[i]=1812433253u*(r->mt[i-1]^(r->mt[i-1]>>30))+(uint32_t)i;
}
static inline uint32_t draw(random_state *r)
{
    if(!r->use_mt){r->state=(uint32_t)(((uint64_t)r->state*48271)%2147483647);return r->state;}
    if(r->next==624){
        for(size_t i=0;i<624;++i){
            uint32_t y=(r->mt[i]&0x80000000u)|(r->mt[(i+1)%624]&0x7fffffffu);
            r->mt[i]=r->mt[(i+397)%624]^(y>>1)^((y&1u) ? 0x9908b0dfu : 0);
        }
        r->next=0;
    }
    uint32_t x=r->mt[r->next++];x^=x>>11;x^=(x<<7)&0x9d2c5680u;x^=(x<<15)&0xefc60000u;x^=x>>18;return x;
}
static inline double uniform(random_state *r)
{
    double range=r->use_mt ? 4294967296.0 : 2147483646.0;
    uint32_t minimum=r->use_mt ? 0 : 1;
    double first=draw(r)-minimum;
    first+=(draw(r)-minimum)*range;
    return first/(range*range);
}
static inline int coin(random_state *r)
{return r->use_mt ? uniform(r)<.5 : ((draw(r)-1)&1u)!=0;}
static inline double lognormal_at(random_state *r,double mean,double sigma)
{
    double z;
    if(r->saved) {r->saved=0;z=r->normal;}
    else {
        double x,y,q;do{x=2*uniform(r)-1;y=2*uniform(r)-1;q=x*x+y*y;}while(q>1 || q==0);
        double factor=sqrt(-2*log(q)/q);z=x*factor;r->normal=y*factor;r->saved=1;
    }
    return exp(mean+sigma*z);
}
static inline double lognormal(random_state *r,double sigma) { return lognormal_at(r,-2,sigma); }
#endif
