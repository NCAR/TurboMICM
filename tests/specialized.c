/* Independent generic-versus-generated kernel and trajectory checks.
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <stdio.h>
#include <omp.h>
#include SPEC_HEADER
#define JOIN_(a,b) a##b
#define JOIN(a,b) JOIN_(a,b)
#define SPEC(a) JOIN(a,SPEC_NAME)

static int compare(const char *kind, const double *actual, const double *expected, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    if (!isfinite(actual[i]) || !isfinite(expected[i]) ||
        fabs(actual[i]-expected[i]) > 1e-30+1e-10*fabs(expected[i])) {
      fprintf(stderr,"%s[%zu]: %.17g != %.17g\n",kind,i,actual[i],expected[i]);
      return 1;
    }
  }
  return 0;
}
static void initialize(ch_state *s) {
  SPEC(ch_init_)(s);
  size_t nc = s->cells;
  for (size_t c = 0; c < nc; ++c) {
    s->conditions[c] = (ch_conditions){260.0+1.3*c, 80000.0+211.0*c, 0};
    ch_ideal_air_density(&s->conditions[c]);
    for (size_t i = 0; i < s->model->nv; ++i) {
      double *y = &s->y[i*nc+c];
      *y = (*y == 0 ? 1e-10 : *y)*(1.0+.01*c);
    }
    for (size_t i = 0; i < s->model->nparam; ++i) {
      const char *name = ch_model_parameter_name(s->model,i);
      double *v = &s->custom[i*nc+c];
      if (strstr(name,"radius")) *v = 1e-7;
      else if (strstr(name,"concentration")) *v = 1e6;
      else if (s->model->nv > 100) { if (name[0] == 'j' && *v >= 1) *v = 0; }
      else if (*v == 0 || *v >= 1) *v = .001;
      *v *= 1+.02*c;
    }
  }
}
static int check(size_t nc, int tight) {
  ch_model *m = SPEC(ch_build_)();
  ch_state *a = ch_state_create(m,nc), *b = ch_state_create(m,nc);
  if (!m || !a || !b) { fprintf(stderr,"allocation failed\n"); return 1; }
  initialize(a);
  if (tight) {
    double atol = 1e-15;
    if (ch_state_tolerances(a,1e-6,&atol,1)) return 1;
  }
  ch_state_copy(b,a);
  ch_rosenbrock_parameters p = ch_rosenbrock_default(3);
  /* Large initial trials exercise regeneration after rejected steps. */
  if (tight) p.h_start = 30.0;
  size_t nv = m->nv*nc, nr = m->nr*nc, nnz = m->sparse->nnz*nc;
  int fail = 0;
  for (int step = 0; step < 5 && !fail; ++step) {
    for (size_t c = 0; c < nc; ++c) {
      a->conditions[c].temperature += 0.7;
      a->conditions[c].pressure += 19;
      ch_ideal_air_density(&a->conditions[c]);
      /* Density is an independent runtime datum too. */
      a->conditions[c].air_density *= 1+.001*step;
      b->conditions[c] = a->conditions[c];
    }
    for (size_t i = 0; i < m->nparam*nc; ++i) b->custom[i] = (a->custom[i] *= .93);
    fail |= ch_update_rates(a) != CH_OK || SPEC(ch_rates_)(b) != CH_OK;
    fail |= compare("rates",b->rates,a->rates,nr);
    ch_force(a,a->y,a->forcing); SPEC(ch_force_)(b,a->y,b->forcing);
    fail |= compare("forcing",b->forcing,a->forcing,nv);
    ch_jacobian(a,a->y,a->jacobian); SPEC(ch_jacobian_)(b,a->y,b->jacobian);
    fail |= compare("Jacobian",b->jacobian,a->jacobian,nnz);
    ch_alpha_minus_jacobian(a,1000.0,a->jacobian); ch_factor_state(a);
    SPEC(ch_factor_)(b,1000.0);
    fail |= compare("LU",b->jacobian,a->jacobian,nnz);
    for (size_t i = 0; i < nv; ++i) a->forcing[i] = b->forcing[i] = 1e-6*(1+i%17);
    ch_linear_state(a,a->forcing); SPEC(ch_linear_)(b,b->forcing);
    fail |= compare("linear solve",b->forcing,a->forcing,nv);
    ch_stats sa, sb;
    int ra = ch_solve(a,30,&p,&sa), rb = SPEC(ch_solve_)(b,30,&p,&sb);
    if (ra || rb || sa.internal_steps != sb.internal_steps || sa.rejected != sb.rejected) {
      fprintf(stderr,"solve: %s / %s; steps %zu/%zu rejects %zu/%zu\n",
              ch_status_string(ra),ch_status_string(rb),sa.internal_steps,sb.internal_steps,sa.rejected,sb.rejected);
      fail = 1;
    }
    fail |= compare("trajectory",b->y,a->y,nv);
  }
  ch_state_free(a); ch_state_free(b); ch_model_free(m);
  return fail;
}
static int profile(size_t nc) {
  ch_model *m = SPEC(ch_build_)();
  ch_state *s = ch_state_create(m,nc);
  if (!s) return 1;
  initialize(s); ch_update_rates(s);
  const int repeats = 2000;
  for (int specialized = 0; specialized < 2; ++specialized) {
    double start = omp_get_wtime();
    for (int i = 0; i < repeats; ++i)
      if (specialized) SPEC(ch_force_)(s,s->y,s->forcing); else ch_force(s,s->y,s->forcing);
    double force = omp_get_wtime()-start;
    start = omp_get_wtime();
    for (int i = 0; i < repeats; ++i)
      if (specialized) SPEC(ch_jacobian_)(s,s->y,s->jacobian); else ch_jacobian(s,s->y,s->jacobian);
    double jac = omp_get_wtime()-start;
    double *saved = malloc(m->sparse->nnz*nc*sizeof(double));
    if (!saved) return 1;
    memcpy(saved,s->jacobian,m->sparse->nnz*nc*sizeof(double));
    start = omp_get_wtime();
    for (int i = 0; i < repeats; ++i) {
      memcpy(s->jacobian,saved,m->sparse->nnz*nc*sizeof(double));
      if (specialized) SPEC(ch_factor_)(s,1000); else { ch_alpha_minus_jacobian(s,1000,s->jacobian); ch_factor_state(s); }
    }
    double factor = omp_get_wtime()-start;
    start = omp_get_wtime();
    for (int i = 0; i < repeats; ++i) {
      memcpy(s->forcing,s->y,m->nv*nc*sizeof(double));
      if (specialized) SPEC(ch_linear_)(s,s->forcing); else ch_linear_state(s,s->forcing);
    }
    printf("%s cells=%zu force_ms=%.3f jacobian_ms=%.3f factor_ms=%.3f linear_ms=%.3f\n",
           specialized?"specialized":"generic",nc,1000*force,1000*jac,1000*factor,1000*(omp_get_wtime()-start));
    free(saved);
  }
  ch_state_free(s); ch_model_free(m); return 0;
}
int main(int argc, char **argv) {
  if (argc == 2) return profile((size_t)strtoul(argv[1],NULL,10));
  const size_t sizes[] = {1,3,7,32,33};
  for (size_t i = 0; i < sizeof(sizes)/sizeof(*sizes); ++i) for (int tight = 0; tight < 2; ++tight)
    if (check(sizes[i],tight)) { fprintf(stderr,"FAIL %s cells=%zu tight=%d\n",SPEC_HEADER,sizes[i],tight); return 1; }
  printf("PASS %s: rates, forcing, Jacobian, LU, solve and adaptive trajectories (1e-10)\n",SPEC_HEADER);
  return 0;
}
