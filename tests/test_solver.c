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
    const ChemMechanism mechanism = {2, 1, 0, names, NULL, init, NULL, reactions, reactants, products};
    ChemConditions conditions = {300, 100000};
    ChemRosenbrockParameters p = chem_rosenbrock_three_stage();
    p.absolute_tolerance = 1e-12; p.relative_tolerance = 1e-9;
    ChemWorkspace *w = chem_workspace_create(&mechanism);
    double rates[1], y[] = {1, 0};
    ChemStatistics stats;
    if (!w || chem_rate_constants(&mechanism, conditions, NULL, rates) ||
        chem_rosenbrock_solve(&mechanism, w, &p, conditions, rates, 1, y, &stats)) return 1;
    if (fabs(y[0] - exp(-2)) > 1e-8 || fabs(y[0]+y[1]-1) > 1e-12 ||
        stats.accepted == 0 || stats.solves != 3*stats.internal_steps) {
        fprintf(stderr, "unexpected solution %.17g %.17g\n", y[0], y[1]); return 2;
    }
    if (!chem_rosenbrock_solve(&mechanism, w, &p, conditions, rates, 0, y, &stats)) return 3;
    chem_workspace_free(w);
    puts("Rosenbrock decay and conservation: PASS");
    return 0;
}
