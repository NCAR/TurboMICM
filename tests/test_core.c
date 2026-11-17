#include "chem.h"

#include <math.h>
#include <stdio.h>

#define CHECK(condition) do { if (!(condition)) { \
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); return 1; \
} } while (0)

static double fixed_step(double error, double step, const ChemRosenbrock *p, void *context)
{
    (void)error; (void)step; (void)p;
    ++*(unsigned *)context;
    return 1;
}

int main(void)
{
    const ChemSpecies species[] = {{"A", 1, 1e-12, 0, 0}, {"B", 0, 1e-12, 0, 0}};
    const int reactants[] = {0};
    const ChemProduct products[] = {{1, 1}};
    const char *const names[] = {"rate"};
    const double defaults[] = {0.1};
    const ChemReaction reactions[] = {
        {0, 1, 0, 1, 0, {.kind = CHEM_USER, .parameter = 0, .a = 1}}
    };
    const ChemMechanism mechanism = {2, 1, 1, species, reactions, reactants, products, names, defaults, 300, 101325};
    ChemState *s = chem_state_create(&mechanism, 3);
    CHECK(s);
    CHECK(chem_find_species(&mechanism, "B") == 1);
    CHECK(chem_find_parameter(&mechanism, "rate") == 0);
    CHECK(chem_find_species(&mechanism, "missing") == -1);
    chem_species_values(s, 0)[1] = 2;
    chem_species_values(s, 0)[2] = 3;
    CHECK(chem_update_rates(s) == CHEM_OK);
    double f[2], j[4];
    CHECK(chem_evaluate(s, 1, f, j) == CHEM_OK);
    CHECK(fabs(f[0]+0.2) < 1e-15 && fabs(f[1]-0.2) < 1e-15);
    CHECK(fabs(j[0]-0.1) < 1e-15 && j[1] == 0 && fabs(j[2]+0.1) < 1e-15 && j[3] == 0);
    ChemRosenbrock p = chem_rosenbrock_three_stage();
    ChemStats stats;
    CHECK(chem_solve(s, &p, 1, &stats) == CHEM_OK);
    CHECK(stats.final_time == 1 && stats.accepted > 0 && stats.solves == 3*stats.internal_steps);
    for (size_t c = 0; c < 3; ++c) {
        double a = chem_species_values(s, 0)[c], b = chem_species_values(s, 1)[c];
        CHECK(fabs(a - (c+1)*exp(-0.1)) < 1e-7);
        CHECK(fabs(a+b-(c+1)) < 1e-12);
    }
    const ChemRosenbrock variants[] = {chem_rosenbrock_two_stage(), chem_rosenbrock_four_stage(),
                                      chem_rosenbrock_four_stage_dae(), chem_rosenbrock_six_stage_dae()};
    for (size_t v = 0; v < sizeof variants / sizeof variants[0]; ++v) {
        for (size_t c = 0; c < 3; ++c) {
            chem_species_values(s, 0)[c] = (double)c+1;
            chem_species_values(s, 1)[c] = 0;
        }
        CHECK(chem_solve(s, &variants[v], 1, &stats) == CHEM_OK);
        for (size_t c = 0; c < 3; ++c) {
            double a = chem_species_values(s, 0)[c], b = chem_species_values(s, 1)[c];
            CHECK(fabs(a - (c+1)*exp(-0.1)) < 1e-7);
            CHECK(fabs(a+b-(c+1)) < 1e-12);
        }
    }
    unsigned calls=0;
    p=chem_rosenbrock_three_stage(); p.h_start=.1; p.step_controller=fixed_step; p.step_controller_context=&calls;
    for (size_t c=0;c<3;++c) { chem_species_values(s,0)[c]=1; chem_species_values(s,1)[c]=0; }
    CHECK(chem_solve(s,&p,1,&stats)==CHEM_OK && calls==stats.internal_steps && calls>=9);
    ChemLinearSolver doolittle=chem_doolittle_linear_solver();
    CHECK(chem_set_linear_solver(s,&doolittle,NULL)==CHEM_OK);
    p=chem_rosenbrock_three_stage();
    for (size_t c=0;c<3;++c) { chem_species_values(s,0)[c]=1; chem_species_values(s,1)[c]=0; }
    CHECK(chem_solve(s,&p,1,&stats)==CHEM_OK);
    for (size_t c=0;c<3;++c) CHECK(fabs(chem_species_values(s,0)[c]-exp(-.1))<1e-7);
    CHECK(chem_set_linear_solver(s,NULL,NULL)==CHEM_OK);
    chem_state_free(s);
    puts("core analytic decay and Jacobian: PASS");
    return 0;
}
