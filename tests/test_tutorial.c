#include "chem.h"
#include "rates.h"
#include "rosenbrock.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

extern const ChemMechanism generated_mechanism;
extern const ChemRate generated_rates[];
extern const size_t generated_custom_count;
typedef struct { const ChemMechanism *m; double k[10]; } System;
static void rhs(void *c, const double *y, double *f) { System *s=c; chem_rhs(s->m,y,s->k,f); }
static void jac(void *c, const double *y, double *j) { System *s=c; chem_jacobian(s->m,y,s->k,j); }
static void close_value(double a, double b) { assert(fabs(a-b) <= 1e-12 * fmax(1.0, fabs(b))); }

int main(void)
{
    const double conv = 6.02214076e17;
    const size_t reactants[] = {0,1,1,2,3,4,4,5,2,1};
    const ChemProduct products[] = {{1,1},{2,1},{3,1},{4,1},{5,2},{6,1},{6,1},{6,1},{0,1}};
    const ChemReaction reactions[] = {
        {0,1,0,1},{1,1,1,1},{2,1,2,1},{3,1,3,1},{4,1,4,1},
        {5,2,5,1},{7,1,6,1},{8,1,7,1},{9,1,8,0},{10,0,8,1}};
    const ChemRate rates[] = {
        {CHEM_RATE_ARRHENIUS,{2.15e-4,0,110,300,0},0},
        {CHEM_RATE_ALKOXY,{1.2,204.3,1e-3,2},0},
        {CHEM_RATE_NITRATE,{1.2,204.3,1e-3,2},0},
        {CHEM_RATE_SURFACE,{2.3e2,.025,.9},0},
        {CHEM_RATE_TERNARY,{1.2,2.3,302.3,2.6/conv,-3.1,402.1,.9,1.2},0},
        {CHEM_RATE_TROE,{1.8e-30*conv*conv,-3,0,2.8e-11*conv,0,0,.6,1},0},
        {CHEM_RATE_TUNNELING,{1.2,2.3,302.3},0},
        {CHEM_RATE_USER,{1},2},{CHEM_RATE_USER,{1},3},{CHEM_RATE_USER,{1},4}};
#ifdef WITH_CUSTOM
    const size_t count = 10;
#else
    const size_t count = 7;
#endif
    ChemMechanism hand = {7,count,reactions,reactants,products};
    ChemConditions conditions = {287.45, 101319.9, 101319.9/(1.380649e-23*6.02214076e23*287.45)};
    const double custom[] = {1e-7,2.5e6,.01,.02,.03};
    assert(generated_mechanism.reaction_count == count);
    assert(generated_custom_count == (count == 7 ? 2 : 5));
    System a = {&generated_mechanism,{0}}, b = {&hand,{0}};
    chem_rates(count, generated_rates, conditions, custom, a.k);
    chem_rates(count, rates, conditions, custom, b.k);
    for (size_t i=0; i<count; ++i) close_value(a.k[i], b.k[i]);
    double y[] = {1,.2,.3,.4,.5,.6,.7}, f[7], g[7], j[49], k[49];
    rhs(&a,y,f); rhs(&b,y,g); jac(&a,y,j); jac(&b,y,k);
    for (size_t i=0; i<7; ++i) close_value(f[i],g[i]);
    for (size_t i=0; i<49; ++i) close_value(j[i],k[i]);
    const double atol[] = {1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3};
    double ya[] = {1,0,0,0,0,0,0}, yb[] = {1,0,0,0,0,0,0};
    ChemRosenbrockMethod method = chem_rosenbrock_three_stage();
    ChemRosenbrockOptions options = chem_rosenbrock_defaults(atol);
#ifndef WITH_CUSTOM
    /* Reference tutorial output at its printed, three-significant-digit precision. */
    const double reference[10][7] = {
        {.854,.000457,.144,.000155,0,1.23e-22,.000644},
        {.730,.000390,.265,.000289,8.18e-7,2.28e-22,.00244},
        {.623,.000333,.366,.000402,9.66e-7,3.18e-22,.00520},
        {.532,.000285,.449,.000500,1.07e-6,3.95e-22,.00877},
        {.455,.000243,.518,.000583,1.15e-6,4.61e-22,.0130},
        {.388,.000208,.575,.000654,1.21e-6,5.17e-22,.0178},
        {.332,.000177,.621,.000714,1.26e-6,5.65e-22,.0230},
        {.283,.000152,.659,.000766,1.30e-6,6.06e-22,.0286},
        {.242,.000129,.688,.000810,1.33e-6,6.41e-22,.0345},
        {.207,.000111,.711,.000848,1.35e-6,6.71e-22,.0406}
    };
#endif
    for (int step=0; step<10; ++step) {
        assert(chem_rosenbrock_solve(7,ya,500,rhs,jac,&a,&method,&options,NULL)==CHEM_SOLVE_OK);
        assert(chem_rosenbrock_solve(7,yb,500,rhs,jac,&b,&method,&options,NULL)==CHEM_SOLVE_OK);
        for (size_t i=0; i<7; ++i) {
            close_value(ya[i], yb[i]);
#ifndef WITH_CUSTOM
            if (reference[step][i] == 0) assert(ya[i] == 0);
            else assert(fabs(ya[i]-reference[step][i]) <= .0051*fabs(reference[step][i]));
#endif
        }
#ifdef PRINT_STATE
        printf("%d", 500*(step+1));
        for (size_t i=0; i<7; ++i) printf(", %.2e", ya[i]);
        putchar('\n');
#endif
    }
    puts("generated tutorial matches hand-built mechanism PASS");
    return 0;
}
