#include "chem.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern const ChemMechanism generated_mechanism;
extern const ChemConditions generated_mechanism_initial_conditions;

static int near(double a, double b, double f) { return fabs(a - b) <= 2e-5 * (1 + fabs(b)) + 1e-8 * fabs(f); }

int main(void)
{
    const ChemMechanism *m = &generated_mechanism;
    double *rates = calloc(m->reaction_count, sizeof *rates);
    double *y = calloc(m->species_count, sizeof *y);
    double *f = calloc(m->species_count, sizeof *f);
    double *fp = calloc(m->species_count, sizeof *fp);
    double *fm = calloc(m->species_count, sizeof *fm);
    double *j = calloc(m->species_count * m->species_count, sizeof *j);
    if (!rates || !y || !f || !fp || !fm || !j) return 1;
    ChemConditions c = generated_mechanism_initial_conditions;
    if (c.temperature!=287.45 || c.pressure!=101319.9 || c.has_air_density) return 8;
    if (chem_rate_constants(m, c, m->initial_parameters, rates)) return 2;
    for (size_t i = 0; i < m->species_count; ++i) y[i] = 0.1 + i * 0.001;
    chem_forcing(m, c, rates, y, f);
    chem_jacobian(m, c, rates, y, j);
    for (size_t col = 0; col < m->species_count; ++col) {
        double original = y[col], h = 1e-4;
        y[col] = original + h;
        chem_forcing(m, c, rates, y, fp);
        y[col] = original - h;
        chem_forcing(m, c, rates, y, fm);
        y[col] = original;
        for (size_t row = 0; row < m->species_count; ++row) {
            double fd = (fp[row] - fm[row]) / (2 * h);
            if (!near(fd, j[row * m->species_count + col], f[row])) {
                fprintf(stderr, "Jacobian mismatch row=%zu col=%zu: %.17g vs %.17g\n", row, col, fd, j[row*m->species_count+col]);
                return 3;
            }
        }
    }
    /* Exercise both LU algorithms on the generated mechanism, including fill. */
    ChemWorkspace *work = chem_workspace_create(m);
    if (!work) return 4;
    for (size_t i = 0; i < m->species_count; ++i) y[i] = fp[i] = m->initial_species[i];
    ChemRosenbrockParameters parameters = chem_rosenbrock_three_stage();
    ChemStatistics stats;
    if (chem_rosenbrock_solve(m, work, &parameters, c, rates, 0.001, y, &stats)) return 5;
    parameters.linear_solver = CHEM_LU_DOOLITTLE;
    if (chem_rosenbrock_solve(m, work, &parameters, c, rates, 0.001, fp, &stats)) return 6;
    for (size_t i = 0; i < m->species_count; ++i)
        if (fabs(y[i]-fp[i]) > 1e-30 + 1e-10*fabs(y[i])) {
            fprintf(stderr, "LU algorithms disagree for %s: %.17g %.17g\n", m->species_names[i], y[i], fp[i]);
            return 7;
        }
    chem_workspace_free(work);
    free(rates); free(y); free(f); free(fp); free(fm); free(j);
    puts("forcing, finite-difference Jacobian, and LU algorithms: PASS");
    return 0;
}
