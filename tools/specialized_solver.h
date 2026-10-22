/* Three-stage Rosenbrock integration specialized to generated chemistry.
 * Coefficients, error estimate, and controller follow MICM and src/rosenbrock.c.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Included once in each generated translation unit, with CH_* names defined.
 */
int CH_SOLVE(ch_state *s, double dt, const ch_rosenbrock_parameters *p, ch_stats *stats) {
  if (!s || !p || !stats || !isfinite(dt) || dt < 0 || p->stages != 3 ||
      p->gamma[0] == 0 || p->estimator_order <= 0 || !p->controller ||
      !p->new_function[1] || p->new_function[2]) return CH_INVALID;
  *stats = (ch_stats){0};
  const size_t n = CH_SPECIES*s->cells;
  const double h_min = p->h_min == 0 ? 1e-15*dt : p->h_min;
  const double h_max = p->h_max == 0 ? dt : fmin(dt,p->h_max);
  const double h_start = p->h_start == 0 ? 1e-6*dt : fmin(h_max,p->h_start);
  double h = fmin(fmax(h_min,fabs(h_start)),fabs(h_max)), time = 0.0;
  int reject_last = 0, reject_more = 0, status = CH_OK;
  double *restrict k0 = s->stages, *restrict k1 = s->stages+n, *restrict k2 = s->stages+2*n;
  while (time-dt+p->round_off <= 0) {
    if (stats->internal_steps > (size_t)p->max_steps) { status = CH_MAX_STEPS; break; }
    if (time+0.1*h == time || h <= fmin(p->round_off,h_min)) { status = CH_SMALL_STEP; break; }
    h = fmin(h,fabs(dt-time));
    CH_FORCE(s,s->y,s->forcing); ++stats->function_calls;
    CH_JACOBIAN(s,s->y,s->jacobian); ++stats->jacobian_updates;
    int accepted = 0;
    while (!accepted) {
      CH_FACTOR(s,1.0/(h*p->gamma[0])); ++stats->decompositions;
      memcpy(k0,s->forcing,n*sizeof(double));
      CH_LINEAR(s,k0);
      double *restrict yn = s->new_y;
      const double *restrict y = s->y;
      #pragma omp simd
      for (size_t i = 0; i < n; ++i) yn[i] = y[i]+p->a[0]*k0[i];
      CH_FORCE(s,yn,k1); ++stats->function_calls;
      const double c0 = p->c[0]/h, c1 = p->c[1]/h, c2 = p->c[2]/h;
      #pragma omp simd
      for (size_t i = 0; i < n; ++i) {
        k2[i] = k1[i];
        k1[i] += c0*k0[i];
      }
      CH_LINEAR(s,k1);
      #pragma omp simd
      for (size_t i = 0; i < n; ++i) {
        k2[i] += c1*k0[i];
        k2[i] += c2*k1[i];
      }
      CH_LINEAR(s,k2); stats->solves += 3;
      double *restrict errors = s->error;
      #pragma omp simd
      for (size_t i = 0; i < n; ++i) {
        double value = y[i];
        value += p->m[0]*k0[i]; value += p->m[1]*k1[i]; value += p->m[2]*k2[i];
        yn[i] = value;
        double error = 0.0;
        error += p->e[0]*k0[i]; error += p->e[1]*k1[i]; error += p->e[2]*k2[i];
        errors[i] = error;
      }
      /* Keep the reference's species-major, ordered reduction. */
      double error = ch_normalized_error(s,y,yn,errors);
      double hnew = h*p->controller(error,p,p->controller_context);
      ++stats->internal_steps;
      if (isnan(error)) { status = CH_NAN; goto done; }
      if (isinf(error)) { status = CH_INF; goto done; }
      if (error < 1 || h < h_min) {
        ++stats->accepted; time += h;
        double *swap = s->y; s->y = s->new_y; s->new_y = swap;
        hnew = fmax(h_min,fmin(hnew,h_max));
        if (reject_last) hnew = fmin(hnew,h);
        reject_last = reject_more = 0; h = hnew; accepted = 1;
      } else {
        if (reject_more) hnew = h*p->rejection_factor;
        reject_more = reject_last; reject_last = 1; h = hnew;
        if (stats->accepted >= 1) ++stats->rejected;
        CH_JACOBIAN(s,s->y,s->jacobian); ++stats->jacobian_updates;
        if (h == 0 || !isfinite(h)) { status = CH_SMALL_STEP; goto done; }
      }
    }
  }
done:
  stats->final_time = time;
  ch_post_solve(s);
  return status;
}
