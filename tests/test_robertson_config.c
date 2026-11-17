#include "chem.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern const ChemMechanism generated_mechanism;

int main(void)
{
    const char *names[] = {"A", "B", "C"}, *parameter_names[] = {"r1", "r2", "r3"};
    const double defaults[] = {0, 0, 0};
    const ChemTerm reactants[] = {{0,1}, {1,1}, {1,1}, {1,1}, {2,1}};
    const ChemTerm products[] = {{1,1}, {1,1}, {2,1}, {0,1}, {2,1}};
    const ChemReaction reactions[] = {
        {CHEM_USER, 0,1,0,1,0,{1}}, {CHEM_USER,1,2,1,2,1,{1}}, {CHEM_USER,3,2,3,2,2,{1}}
    };
    const ChemMechanism hand = {3,3,3,names,parameter_names,defaults,defaults,reactions,reactants,products,0,NULL};
    const ChemMechanism *mechanisms[] = {&hand, &generated_mechanism};
    double solutions[2][3];
    for (size_t m = 0; m < 2; ++m) {
        const ChemMechanism *mechanism = mechanisms[m];
        ChemWorkspace *work = chem_workspace_create(mechanism);
        double y[3] = {0}, rates[3], params[3];
        int index[3] = {-1,-1,-1};
        if (!work) return 1;
        for (size_t i = 0; i < 3; ++i) {
            int species = mechanism->species_names[i][0] - 'A';
            if (species < 0 || species >= 3) return 1;
            index[species] = (int)i;
            if (!strcmp(mechanism->parameter_names[i], "r1")) params[i] = .01;
            else if (!strcmp(mechanism->parameter_names[i], "r2")) params[i] = .02;
            else if (!strcmp(mechanism->parameter_names[i], "r3")) params[i] = .03;
            else return 1;
        }
        y[index[0]] = 1;
        ChemConditions conditions = {300, 100000, 0, 0};
        ChemRosenbrockParameters p = chem_rosenbrock_three_stage();
        p.absolute_tolerance = 1e-12; p.relative_tolerance = 1e-9;
        if (chem_rate_constants(mechanism, conditions, params, rates)) return 1;
        for (int step = 0; step < 10; ++step) {
            ChemStatistics stats;
            if (chem_rosenbrock_solve(mechanism, work, &p, conditions, rates, 1, y, &stats)) return 1;
        }
        for (size_t i = 0; i < 3; ++i) solutions[m][i] = y[index[i]];
        chem_workspace_free(work);
    }
    for (size_t i = 0; i < 3; ++i) {
        if (fabs(solutions[0][i]-solutions[1][i]) > 1e-14 + 1e-12*fabs(solutions[0][i])) return 2;
    }
    puts("Robertson configuration and C-built mechanism: PASS");
    return 0;
}
