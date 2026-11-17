#include "chem.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

extern const ChemMechanism generated_mechanism;

/* The same mechanism constructed in C, using SI rate constants. */
static ChemMechanism hand_mechanism(int user)
{
    static const char *names[] = {"A","B","C","D","E","F","G"};
    static const char *params[] = {"C surface.effective radius [m]", "C surface.particle number concentration [# m-3]", "my photolysis rate", "my loss rate", "my emission rate"};
    static const double initial[] = {0,0,0,0,0,0,0}, defaults[] = {0,0,0,0,0};
    static const ChemTerm react[] = {{0,1},{1,1},{1,1},{2,1},{3,1},{4,1},{4,1},{5,1},{2,1},{1,1}};
    static const ChemTerm prod[] = {{1,1},{2,1},{3,1},{4,1},{5,2},{6,1},{6,1},{6,1},{0,1}};
    static ChemReaction r[10];
    const double factor = 1e-6 * 6.02214076e23;
    r[0] = (ChemReaction){CHEM_ARRHENIUS,0,1,0,1,-1,{2.15e-4,0,110,300,0}};
    r[1] = (ChemReaction){CHEM_BRANCH_ALKOXY,1,1,1,1,-1,{1.2,204.3,1e-3,2}};
    r[2] = (ChemReaction){CHEM_BRANCH_NITRATE,2,1,2,1,-1,{1.2,204.3,1e-3,2}};
    r[3] = (ChemReaction){CHEM_SURFACE,3,1,3,1,0,{230,.025,.9}};
    r[4] = (ChemReaction){CHEM_TERNARY,4,1,4,1,-1,{1.2,2.3,302.3,2.6/factor,-3.1,402.1,.9,1.2}};
    r[5] = (ChemReaction){CHEM_TROE,5,2,5,1,-1,{1.8e-30*factor*factor,-3,0,2.8e-11*factor,0,0,.6,1}};
    r[6] = (ChemReaction){CHEM_TUNNELING,7,1,6,1,-1,{1.2,2.3,302.3}};
    r[7] = (ChemReaction){CHEM_USER,8,1,7,1,2,{1}};
    r[8] = (ChemReaction){CHEM_USER,9,1,8,0,3,{1}};
    r[9] = (ChemReaction){CHEM_USER,10,0,8,1,4,{1}};
    return (ChemMechanism){7, user ? 10 : 7, user ? 5 : 2, names, params, initial, defaults, r, react, prod};
}

static int run(const ChemMechanism *m, int user, double output[11][7])
{
    ChemConditions cond = {287.45, 101319.9, 0, 0};
    ChemRosenbrockParameters p = chem_rosenbrock_three_stage();
    ChemWorkspace *work = chem_workspace_create(m);
    double *rates = calloc(m->reaction_count, sizeof(double));
    double *params = calloc(m->parameter_count, sizeof(double));
    double *y = calloc(m->species_count, sizeof(double));
    if (!work || !rates || !params || !y) return 1;
    size_t map[7];
    for (int j = 0; j < 7; ++j) {
        map[j] = m->species_count;
        for (size_t i = 0; i < m->species_count; ++i) if (m->species_names[i][0] == 'A'+j && !m->species_names[i][1]) map[j] = i;
        if (map[j] == m->species_count) return 1;
    }
    y[map[0]] = 1;
    for (size_t i = 0; i < m->parameter_count; ++i) {
        if (strstr(m->parameter_names[i], "effective radius")) params[i] = 1e-7;
        else if (strstr(m->parameter_names[i], "particle number")) params[i] = 2.5e6;
        else if (!strcmp(m->parameter_names[i], "my emission rate")) params[i] = 1e-20;
        else if (!strcmp(m->parameter_names[i], "my loss rate")) params[i] = 1e-23;
    }
    double photo = 1e-10;
    for (int step = 0; step <= 10; ++step) {
        for (int j = 0; j < 7; ++j) output[step][j] = y[map[j]];
        if (step == 10) break;
        if (user) for (size_t i = 0; i < m->parameter_count; ++i)
            if (!strcmp(m->parameter_names[i], "my photolysis rate")) params[i] = photo;
        ChemStatistics stats;
        if (chem_rate_constants(m, cond, params, rates) || chem_rosenbrock_solve(m, work, &p, cond, rates, 500, y, &stats)) return 1;
        photo *= 1.5;
    }
    chem_workspace_free(work); free(rates); free(params); free(y);
    return 0;
}

int main(int argc, char **argv)
{
    int user = argc > 1 && !strcmp(argv[1], "user");
    ChemMechanism hand = hand_mechanism(user);
    double a[11][7], b[11][7];
    if (run(&generated_mechanism, user, a) || run(&hand, user, b)) return 1;
    /* Check the printed values as well as configuration vs hand-built state. */
    puts("  time,          A,          B,          C,          D,          E,          F,          G");
    for (int i = 0; i <= 10; ++i) {
        printf("%6d", i*500);
        for (int j = 0; j < 7; ++j) {
            if (fabs(a[i][j] - b[i][j]) > 1e-25 + 1e-12*fabs(b[i][j])) {
                fprintf(stderr, "configuration and hand-built states differ at %d,%d: %.17g %.17g\n", i, j, a[i][j], b[i][j]); return 1;
            }
            char left[40], right[40];
            snprintf(left, sizeof left, "%11.2e", a[i][j]);
            snprintf(right, sizeof right, "%11.2e", b[i][j]);
            if (strcmp(left, right)) { fprintf(stderr, "configuration and hand-built differ at %d,%d\n", i, j); return 1; }
            printf(",%s", left);
        }
        putchar('\n');
    }
    return 0;
}
