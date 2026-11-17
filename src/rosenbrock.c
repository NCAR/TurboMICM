#include "chem.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct ChemWorkspace {
    size_t n;
    unsigned char *pattern;
    size_t *lower_offset, *upper_offset, *lower_rows, *upper_cols;
    double *matrix, *jacobian, *initial, *trial, *error, *stage;
};

ChemRosenbrockParameters chem_rosenbrock_three_stage(void)
{
    ChemRosenbrockParameters p = {0};
    p.stages = 3; p.new_function[0] = p.new_function[1] = 1; p.max_steps = 1000;
    p.a[0] = 1; p.a[1] = 1;
    p.c[0] = -0.10156171083877702091975600115545e+01;
    p.c[1] = 0.40759956452537699824805835358067e+01;
    p.c[2] = 0.92076794298330791242156818474003e+01;
    p.m[0] = 1; p.m[1] = 0.61697947043828245592553615689730e+01;
    p.m[2] = -0.42772256543218573326238373806514;
    p.e[0] = 0.5; p.e[1] = -0.29079558716805469821718236208017e+01;
    p.e[2] = 0.22354069897811569627360909276199;
    p.gamma = 0.43586652150845899941601945119356; p.order = 3;
    p.factor_min = 0.2; p.factor_max = 6; p.rejection_factor = 0.1; p.safety = 0.9;
    p.relative_tolerance = 1e-6; p.absolute_tolerance = 1e-3;
    return p;
}

ChemRosenbrockParameters chem_rosenbrock_two_stage(void)
{
    ChemRosenbrockParameters p = chem_rosenbrock_three_stage();
    double g = 1.0 + 1.0/sqrt(2.0);
    memset(p.a, 0, sizeof p.a); memset(p.c, 0, sizeof p.c);
    memset(p.m, 0, sizeof p.m); memset(p.e, 0, sizeof p.e);
    p.stages = 2; p.a[0] = 1.0/g; p.c[0] = -2.0/g;
    p.m[0] = 3.0/(2.0*g); p.m[1] = 1.0/(2.0*g);
    p.e[0] = p.e[1] = 1.0/(2.0*g); p.gamma = g; p.order = 2;
    p.new_function[0] = p.new_function[1] = 1;
    return p;
}

ChemRosenbrockParameters chem_rosenbrock_four_stage(void)
{
    ChemRosenbrockParameters p = chem_rosenbrock_three_stage();
    memset(p.a, 0, sizeof p.a); memset(p.c, 0, sizeof p.c);
    memset(p.m, 0, sizeof p.m); memset(p.e, 0, sizeof p.e);
    p.stages = 4; p.new_function[0] = p.new_function[1] = p.new_function[2] = 1; p.new_function[3] = 0;
    p.a[0] = 2.0; p.a[1] = p.a[3] = 1.867943637803922; p.a[2] = p.a[4] = .2344449711399156;
    p.c[0] = -7.137615036412310; p.c[1] = 2.580708087951457; p.c[2] = .6515950076447975;
    p.c[3] = -2.137148994382534; p.c[4] = -.3214669691237626; p.c[5] = -.6949742501781779;
    p.m[0] = 2.255570073418735; p.m[1] = .2870493262186792; p.m[2] = .4353179431840180; p.m[3] = 1.093502252409163;
    p.e[0] = -.2815431932141155; p.e[1] = -.07276199124938920; p.e[2] = -.1082196201495311; p.e[3] = -1.093502252409163;
    p.gamma = .57282; p.order = 4;
    return p;
}

ChemWorkspace *chem_workspace_create(const ChemMechanism *m)
{
    if (!m || chem_validate_mechanism(m, NULL, 0) != CHEM_VALID || m->species_count > (size_t)-1 / m->species_count) return NULL;
    size_t n = m->species_count, nn = n*n;
    ChemWorkspace *w = calloc(1, sizeof *w);
    if (!w) return NULL;
    w->n = n;
    w->pattern = calloc(nn, 1);
    w->lower_offset = calloc(n+1, sizeof(size_t)); w->upper_offset = calloc(n+1, sizeof(size_t));
    w->lower_rows = malloc(nn*sizeof(size_t)); w->upper_cols = malloc(nn*sizeof(size_t));
    w->matrix = calloc(nn, sizeof(double)); w->jacobian = calloc(nn, sizeof(double));
    w->initial = calloc(n, sizeof(double)); w->trial = calloc(n, sizeof(double));
    w->error = calloc(n, sizeof(double)); w->stage = calloc(6*n, sizeof(double));
    if (!w->pattern || !w->lower_offset || !w->upper_offset || !w->lower_rows || !w->upper_cols || !w->matrix || !w->jacobian || !w->initial || !w->trial || !w->error || !w->stage) {
        chem_workspace_free(w); return NULL;
    }
    for (size_t i = 0; i < n; ++i) w->pattern[i*n+i] = 1;
    for (size_t r = 0; r < m->reaction_count; ++r) {
        const ChemReaction *rx = &m->reactions[r];
        for (size_t q = 0; q < rx->reactant_count; ++q) {
            int col = m->reactants[rx->reactant_begin+q].species;
            if (col < 0) continue;
            for (size_t k = 0; k < rx->reactant_count; ++k) {
                int row = m->reactants[rx->reactant_begin+k].species;
                if (row >= 0) w->pattern[(size_t)row*n+(size_t)col] = 1;
            }
            for (size_t k = 0; k < rx->product_count; ++k) {
                int row = m->products[rx->product_begin+k].species;
                w->pattern[(size_t)row*n+(size_t)col] = 1;
            }
        }
    }
    /* Symbolic fill for the in-place Mozart elimination. */
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i+1; j < n; ++j) if (w->pattern[j*n+i])
            for (size_t k = i+1; k < n; ++k) if (w->pattern[i*n+k]) w->pattern[j*n+k] = 1;
    size_t nl = 0, nu = 0;
    for (size_t i = 0; i < n; ++i) {
        w->lower_offset[i] = nl; w->upper_offset[i] = nu;
        for (size_t j = i+1; j < n; ++j) if (w->pattern[j*n+i]) w->lower_rows[nl++] = j;
        for (size_t k = i+1; k < n; ++k) if (w->pattern[i*n+k]) w->upper_cols[nu++] = k;
    }
    w->lower_offset[n] = nl; w->upper_offset[n] = nu;
    return w;
}

void chem_workspace_free(ChemWorkspace *w)
{
    if (!w) return;
    free(w->pattern); free(w->lower_offset); free(w->upper_offset); free(w->lower_rows); free(w->upper_cols); free(w->matrix); free(w->jacobian); free(w->initial);
    free(w->trial); free(w->error); free(w->stage); free(w);
}

static int factor(ChemWorkspace *w, double alpha, ChemLUAlgorithm algorithm)
{
    size_t n = w->n;
    double *a = w->matrix;
    for (size_t i = 0; i < n*n; ++i) a[i] = -w->jacobian[i];
    for (size_t i = 0; i < n; ++i) a[i*n+i] += alpha;
    if (algorithm == CHEM_LU_DOOLITTLE) {
        const unsigned char *pattern = w->pattern;
        for (size_t i = 0; i < n; ++i) {
            for (size_t k = i; k < n; ++k) if (pattern[i*n+k])
                for (size_t j = 0; j < i; ++j) if (pattern[i*n+j] && pattern[j*n+k])
                    a[i*n+k] -= a[i*n+j]*a[j*n+k];
            if (a[i*n+i] == 0 || !isfinite(a[i*n+i])) return 1;
            for (size_t k = i+1; k < n; ++k) if (pattern[k*n+i]) {
                for (size_t j = 0; j < i; ++j) if (pattern[k*n+j] && pattern[j*n+i])
                    a[k*n+i] -= a[k*n+j]*a[j*n+i];
                a[k*n+i] /= a[i*n+i];
            }
        }
        return 0;
    }
    if (algorithm != CHEM_LU_MOZART) return 1;
    for (size_t i = 0; i < n; ++i) {
        if (a[i*n+i] == 0 || !isfinite(a[i*n+i])) return 1;
        double inverse = 1.0 / a[i*n+i];
        size_t lb = w->lower_offset[i], le = w->lower_offset[i+1];
        for (size_t j = lb; j < le; ++j) a[w->lower_rows[j]*n+i] *= inverse;
        for (size_t u = w->upper_offset[i]; u < w->upper_offset[i+1]; ++u) {
            size_t k = w->upper_cols[u];
            for (size_t j = lb; j < le; ++j) {
                size_t row = w->lower_rows[j]*n;
                a[row+k] -= a[row+i] * a[i*n+k];
            }
        }
    }
    return 0;
}

static void substitute(const ChemWorkspace *w, double *x)
{
    size_t n = w->n;
    const double *a = w->matrix;
    const unsigned char *s = w->pattern;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < i; ++j) if (s[i*n+j]) x[i] -= a[i*n+j]*x[j];
    for (size_t i = n; i-- > 0;) {
        for (size_t j = i+1; j < n; ++j) if (s[i*n+j]) x[i] -= a[i*n+j]*x[j];
        x[i] /= a[i*n+i];
    }
}

int chem_rosenbrock_solve(const ChemMechanism *m, ChemWorkspace *w,
                          const ChemRosenbrockParameters *p, ChemConditions cond,
                          const double *rates, double dt, double *y, ChemStatistics *stats)
{
    if (!m || !w || !p || !rates || !y || !stats || w->n != m->species_count ||
        !(dt > 0) || p->stages < 1 || p->stages > 6 || p->gamma <= 0 || p->order <= 0 ||
        p->absolute_tolerance <= 0 || p->relative_tolerance < 0) return 1;
    memset(stats, 0, sizeof *stats);
    size_t n = w->n;
    double hmin = p->h_min == 0 ? 1e-15*dt : p->h_min;
    double hmax = p->h_max == 0 ? dt : fmin(dt, p->h_max);
    double h = p->h_start == 0 ? 1e-6*dt : fmin(hmax, p->h_start);
    h = fmin(fmax(hmin, fabs(h)), fabs(hmax));
    double time = 0;
    int reject_last = 0, reject_more = 0;
    while (time - dt + DBL_EPSILON <= 0) {
        if (stats->internal_steps > (unsigned)p->max_steps || time + 0.1*h == time || h <= fmin(DBL_EPSILON, hmin)) return 2;
        h = fmin(h, fabs(dt-time));
        chem_forcing(m, cond, rates, y, w->initial); ++stats->function_calls;
        chem_jacobian(m, cond, rates, y, w->jacobian); ++stats->jacobian_updates;
        int accepted = 0;
        while (!accepted) {
            if (factor(w, 1.0/(h*p->gamma), p->linear_solver)) return 3;
            ++stats->decompositions;
            for (int stage = 0; stage < p->stages; ++stage) {
                double *ks = w->stage + (size_t)stage*n;
                int offset = stage*(stage-1)/2;
                if (stage == 0) memcpy(ks, w->initial, n*sizeof(double));
                else if (p->new_function[stage]) {
                    memcpy(w->trial, y, n*sizeof(double));
                    for (int j = 0; j < stage; ++j)
                        for (size_t i = 0; i < n; ++i) w->trial[i] += p->a[offset+j]*w->stage[(size_t)j*n+i];
                    chem_forcing(m, cond, rates, w->trial, ks); ++stats->function_calls;
                }
                if (stage+1 < p->stages && !p->new_function[stage+1]) memcpy(ks+n, ks, n*sizeof(double));
                for (int j = 0; j < stage; ++j)
                    for (size_t i = 0; i < n; ++i) ks[i] += (p->c[offset+j]/h)*w->stage[(size_t)j*n+i];
                substitute(w, ks); ++stats->solves;
            }
            memcpy(w->trial, y, n*sizeof(double)); memset(w->error, 0, n*sizeof(double));
            for (int stage = 0; stage < p->stages; ++stage)
                for (size_t i = 0; i < n; ++i) {
                    w->trial[i] += p->m[stage]*w->stage[(size_t)stage*n+i];
                    w->error[i] += p->e[stage]*w->stage[(size_t)stage*n+i];
                }
            double err = 0;
            for (size_t i = 0; i < n; ++i) {
                double atol = p->absolute_tolerances ? p->absolute_tolerances[i] : p->absolute_tolerance;
                if (!(atol > 0)) return 1;
                double scale = atol + p->relative_tolerance*fmax(fabs(y[i]), fabs(w->trial[i]));
                double e = w->error[i]/scale; err += e*e;
            }
            err = sqrt(err/n);
            if (!isfinite(err)) return 4;
            err = fmax(err, 1e-10);
            double fac = fmin(p->factor_max, fmax(p->factor_min, p->safety/pow(err, 1/p->order)));
            double next = h*fac; ++stats->internal_steps;
            if (err < 1 || h < hmin) {
                ++stats->accepted; time += h; stats->final_time = time; memcpy(y, w->trial, n*sizeof(double));
                next = fmax(hmin, fmin(next, hmax));
                if (reject_last) next = fmin(next, h);
                reject_last = reject_more = 0; h = next; accepted = 1;
            } else {
                if (reject_more) next = h*p->rejection_factor;
                reject_more = reject_last; reject_last = 1; h = next;
                if (stats->accepted) ++stats->rejected;
                ++stats->jacobian_updates;
                if (h <= fmin(DBL_EPSILON, hmin)) return 2;
            }
        }
    }
    for (size_t i = 0; i < n; ++i) y[i] = y[i] > 0.0 ? y[i] : 0.0;
    return 0;
}

ChemBackwardEulerParameters chem_backward_euler_defaults(void)
{
    ChemBackwardEulerParameters p = {0};
    p.small = 1e-40; p.max_iterations = 11;
    p.reductions[0] = p.reductions[1] = p.reductions[2] = p.reductions[3] = 0.5;
    p.reductions[4] = 0.1;
    p.relative_tolerance = 1e-6; p.absolute_tolerance = 1e-3;
    return p;
}

int chem_backward_euler_converged(size_t cells, size_t species, const double *residual,
                                  const double *y, double small, double rtol, const double *atol)
{
    for (size_t c = 0; c < cells; ++c) for (size_t i = 0; i < species; ++i) {
        double r = fabs(residual[c*species+i]);
        if (!isfinite(r) || (r > small && r > atol[i] && r > rtol*fabs(y[c*species+i]))) return 0;
    }
    return 1;
}

int chem_backward_euler_solve(const ChemMechanism *m, ChemWorkspace *w,
                              const ChemBackwardEulerParameters *p, ChemConditions cond,
                              const double *rates, double dt, double *y, ChemStatistics *stats)
{
    if (!m || !w || !p || !rates || !y || !stats || w->n != m->species_count ||
        !(dt > 0) || p->max_iterations < 1 || p->small < 0 || p->relative_tolerance < 0 ||
        (!p->absolute_tolerances && p->absolute_tolerance < 0)) return 1;
    memset(stats, 0, sizeof *stats);
    size_t n = w->n;
    double h = p->h_start == 0 ? dt : p->h_start, time = 0;
    if (!(h > 0)) return 1;
    int successes = 0, failures = 0;
    memcpy(w->initial, y, n*sizeof(double));
    while (time < dt) {
        int converged = 0;
        for (int iteration = 0; iteration < p->max_iterations; ++iteration) {
            ++stats->internal_steps;
            chem_forcing(m, cond, rates, y, w->trial); ++stats->function_calls;
            chem_jacobian(m, cond, rates, y, w->jacobian); ++stats->jacobian_updates;
            if (factor(w, 1/h, p->linear_solver)) return 3;
            ++stats->decompositions;
            for (size_t i = 0; i < n; ++i) w->trial[i] -= (y[i]-w->initial[i])/h;
            substitute(w, w->trial); ++stats->solves;
            for (size_t i = 0; i < n; ++i) {
                double next = y[i]+w->trial[i];
                y[i] = 0.0 > next ? 0.0 : next;
            }
            if (iteration == 0) continue;
            converged = 1;
            for (size_t i = 0; i < n; ++i) {
                double a = p->absolute_tolerances ? p->absolute_tolerances[i] : p->absolute_tolerance;
                double residual = fabs(w->trial[i]);
                if (!isfinite(residual) || (residual > p->small && residual > a && residual > p->relative_tolerance*fabs(y[i]))) {
                    converged = 0; break;
                }
            }
            if (converged) break;
        }
        if (!converged) {
            ++stats->rejected; successes = 0;
            if (failures >= 5) { stats->final_time = time+h; return 2; }
            memcpy(y, w->initial, n*sizeof(double));
            h *= p->reductions[failures++];
            if (!(h > 0) || time+h == time) return 2;
        } else {
            ++stats->accepted; time += h; stats->final_time = time;
            memcpy(w->initial, y, n*sizeof(double));
            if (++successes >= 2) { successes = 0; h *= 2; }
        }
        h = fmin(h, dt-time);
    }
    return 0;
}
