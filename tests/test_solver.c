#include "chem.h"
#include <math.h>
#include <stdio.h>

int main(void)
{
    /* A -> B, with k = 2/s; total concentration is conserved. */
    const char *names[] = {"A", "B"};
    const double init[] = {1, 0};
    const ChemTerm reactants[] = {{0, 1}}, products[] = {{1, 1}};
    const ChemReaction reactions[] = {{CHEM_ARRHENIUS, 0, 1, 0, 1, -1, {2, 0, 0, 300, 0}}};
    const ChemMechanism mechanism = {2, 1, 0, names, NULL, init, NULL, reactions, reactants, products, 0, NULL};
    ChemConditions conditions = {300, 100000, 0, 0};
    ChemWorkspace *w = chem_workspace_create(&mechanism);
    double rates[1];
    if (!w || chem_rate_constants(&mechanism, conditions, NULL, rates)) return 1;
    for (int algorithm = 0; algorithm < 2; ++algorithm) for (int method = 0; method < 3; ++method) {
        ChemRosenbrockParameters p = method == 0 ? chem_rosenbrock_two_stage() :
                                     method == 1 ? chem_rosenbrock_three_stage() : chem_rosenbrock_four_stage();
        p.absolute_tolerance = 1e-9; p.relative_tolerance = 1e-6;
        p.linear_solver = (ChemLUAlgorithm)algorithm;
        double y[] = {1, 0};
        ChemStatistics stats;
        if (chem_rosenbrock_solve(&mechanism, w, &p, conditions, rates, .1, y, &stats)) return 1;
        if (fabs(y[0] - exp(-.2)) > 1e-6 || fabs(y[0]+y[1]-1) > 1e-12 ||
            stats.accepted == 0 || stats.solves != (unsigned)p.stages*stats.internal_steps) {
            fprintf(stderr, "unexpected method %d solution %.17g %.17g\n", method, y[0], y[1]); return 2;
        }
        if (!chem_rosenbrock_solve(&mechanism, w, &p, conditions, rates, 0, y, &stats)) return 3;
    }
    ChemBackwardEulerParameters be = chem_backward_euler_defaults();
    for (int algorithm = 0; algorithm < 2; ++algorithm) {
        be.linear_solver = (ChemLUAlgorithm)algorithm;
        double backward_y[] = {1, 0};
        ChemStatistics backward_stats;
        if (chem_backward_euler_solve(&mechanism, w, &be, conditions, rates, .1, backward_y, &backward_stats) ||
            fabs(backward_y[0] - 1.0/1.2) > 1e-14 || fabs(backward_y[0]+backward_y[1]-1) > 1e-14) return 4;
    }
    chem_workspace_free(w);
    puts("Rosenbrock and backward Euler decay and conservation: PASS");
    return 0;
}
