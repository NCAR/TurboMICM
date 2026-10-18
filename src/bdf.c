#include "adaptive_internal.h"

/* Arbitrary-step BDF coefficients, obtained by differentiating the Lagrange
 * interpolant at the new endpoint. Nodes are scaled by the current step. */
static void coefficients(int q, const double *times, double h, double *d) {
  double x[6] = {0}; d[0] = 0.0;
  for (int j = 1; j <= q; ++j) { x[j] = (times[j-1]-h)/h; d[0] -= 1.0/x[j]; }
  for (int j = 1; j <= q; ++j) {
    double value = 1.0/x[j];
    for (int k = 1; k <= q; ++k) if (k != j) value *= -x[k]/(x[j]-x[k]);
    d[j] = value;
  }
}
/* Endpoint interpolation and its first omitted term share the same nodes.
 * Multiplying the (q+1)st divided difference by prod(h-t_i)/d0 estimates
 * the BDF local solution error. This is an asymptotic estimator, tested by
 * tolerance refinement and independent scenario references. */
static void extrapolate(const ChemState *s, const ChemWorkspace *w, int points, double h, double *out) {
  size_t n = s->model->nvariables*s->stride;
  memset(out,0,n*sizeof(double));
  for (int j = 0; j < points; ++j) {
    double a = 1.0;
    for (int k = 0; k < points; ++k) if (j != k) a *= (h-w->times[k])/(w->times[j]-w->times[k]);
    chem_axpy(n,out,a,w->history+(size_t)j*n);
  }
}
static double error_estimate(const ChemAdaptive *p, const ChemState *s, ChemWorkspace *w,
                             int q, double h, const double *candidate, double *error) {
  size_t n = s->model->nvariables*s->stride;
  if (w->history_count < (size_t)q+1) return INFINITY;
  extrapolate(s,w,q+1,h,error);
  double d[6]; coefficients(q,w->times,h,d);
  double factor = 1.0/(d[0]*(1.0-w->times[q]/h));
  for (size_t a = 0; a < n; ++a) error[a] = factor*(candidate[a]-error[a]);
  return chem_adaptive_norm(p,s,s->variables,candidate,error);
}
static void factor(ChemState *s, ChemWorkspace *w, double alpha, ChemAdaptiveStats *stats) {
  chem_jacobian(s,s->variables,w->jacobian); ++stats->base.jacobian_updates;
  memcpy(s->jacobian,w->jacobian,s->model->nnz*s->stride*sizeof(double));
  chem_alpha_minus_jacobian(s,alpha);
  chem_lu_mozart(s->model,s->cells,s->stride,s->jacobian); ++stats->base.decompositions;
  w->factor_h = alpha; w->factor_valid = 1; w->jacobian_age = 0;
}
ChemStatus chem_bdf_solve(const ChemAdaptive *p, double dt, ChemState *s, ChemWorkspace *w, ChemAdaptiveStats *stats) {
  size_t n = s->model->nvariables*s->stride;
  double *constant = w->scratch, *predictor = constant+n, *delta = predictor+n;
  double *error_vec = delta+n;
  if (!w->history_count) {
    memcpy(w->history,s->variables,n*sizeof(double)); w->times[0] = 0.0;
    w->history_count = 1; w->order = 1; w->factor_valid = 0;
  }
  double t = 0.0, h = w->next_h > 0.0 ? w->next_h : p->h_start;
  double hmax = p->h_max > 0.0 ? fmin(p->h_max,dt) : dt;
  int successes = 0, rejected = 0;
  while (t < dt) {
    if (stats->base.internal_steps >= p->max_steps) return CHEM_MAX_STEPS;
    h = fmin(h,fmin(hmax,dt-t));
    if (!(h > 0.0) || t+h == t) return CHEM_SMALL_STEP;
    int q = w->order;
    if (q > p->max_order) q = p->max_order;
    if ((size_t)q > w->history_count) q = (int)w->history_count;
    double d[6]; coefficients(q,w->times,h,d); double alpha = d[0]/h;
    memset(constant,0,n*sizeof(double));
    for (int j = 1; j <= q; ++j) chem_axpy(n,constant,-d[j]/h,w->history+(size_t)(j-1)*n);
    if (w->history_count == 1) memcpy(predictor,s->variables,n*sizeof(double));
    else extrapolate(s,w,q < (int)w->history_count ? q+1 : q,h,predictor);
    /* A nonnegative predictor avoids evaluating polynomial chemistry at a large
     * extrapolated negative concentration. Newton still solves the full system. */
    for (size_t a = 0; a < n; ++a) s->new_y[a] = fmax(0.0,predictor[a]);
    int updated = 0;
    if (!p->reuse || !w->factor_valid || w->jacobian_age >= p->jacobian_age ||
        fabs(alpha/w->factor_h-1.0) > p->reuse_ratio) {
      factor(s,w,alpha,stats); updated = 1;
    } else ++stats->matrix_reuses;
    int converged = 0;
    double last = INFINITY;
    /* A failed chord iteration first retries with a current matrix, then cuts h. */
    for (int retry = 0; retry < 2 && !converged; ++retry) {
      for (int iteration = 0; iteration < 8; ++iteration) {
        chem_forcing(s,s->new_y,delta); ++stats->base.function_calls;
        for (size_t a = 0; a < n; ++a) delta[a] += constant[a]-alpha*s->new_y[a];
        chem_linear_solve(s->model,s->cells,s->stride,s->jacobian,delta);
        ++stats->base.solves; ++stats->nonlinear_iterations;
        double norm = chem_adaptive_norm(p,s,s->variables,s->new_y,delta);
        if (!isfinite(norm)) break;
        chem_axpy(n,s->new_y,1.0,delta);
        if (norm < 0.03) { converged = 1; break; }
        if (iteration && norm > last*0.9) break;
        last = norm;
      }
      if (converged || updated) break;
      factor(s,w,alpha,stats); updated = 1; last = INFINITY;
      for (size_t a = 0; a < n; ++a) s->new_y[a] = fmax(0.0,predictor[a]);
    }
    ++stats->base.internal_steps;
    double error = INFINITY;
    if (converged) {
      if (w->history_count == 1) {
        /* Backward Euler vs its explicit tangent predictor for startup. */
        chem_forcing(s,s->variables,error_vec); ++stats->base.function_calls;
        for (size_t a = 0; a < n; ++a) error_vec[a] = 0.5*(s->new_y[a]-s->variables[a]-h*error_vec[a]);
        error = chem_adaptive_norm(p,s,s->variables,s->new_y,error_vec);
      } else error = error_estimate(p,s,w,q,h,s->new_y,error_vec);
    } else { ++stats->failed_iterations; w->factor_valid = 0; }
    int positive = chem_adaptive_positive(p,s,s->variables,s->new_y);
    if (converged && positive && error <= 1.0) {
      ++stats->base.accepted; ++w->jacobian_age; ++successes;
      size_t keep = w->history_count < 7 ? w->history_count : 6;
      for (size_t j = keep; j > 0; --j) {
        memcpy(w->history+j*n,w->history+(j-1)*n,n*sizeof(double)); w->times[j] = w->times[j-1]-h;
      }
      for (size_t a = 0; a < n; ++a) s->variables[a] = fmax(0.0,s->new_y[a]);
      memcpy(w->history,s->variables,n*sizeof(double)); w->times[0] = 0.0;
      if (w->history_count < 7) ++w->history_count;
      t = h == dt-t ? dt : t+h; stats->base.final_time = t;
      double ratio = p->safety*pow(fmax(error,1e-12),-1.0/(q+1)); int next_order = q;
      if (successes >= q+1 && !rejected) {
        if (q > 1) {
          extrapolate(s,w,q+1,h,predictor);
          double low = error_estimate(p,s,w,q-1,h,predictor,error_vec);
          double other = p->safety*pow(fmax(low,1e-12),-1.0/q);
          if (other > ratio*1.2) { ratio = other; next_order = q-1; }
        }
        /* Raising the order needs an extra historical point for its estimator.
         * The candidate is assessed at a future endpoint using history alone. */
        if (q < p->max_order && w->history_count >= (size_t)q+3) {
          extrapolate(s,w,q+3,h,predictor);
          double high = error_estimate(p,s,w,q+1,h,predictor,error_vec);
          double other = p->safety*pow(fmax(high,1e-12),-1.0/(q+2));
          if (other > ratio*1.1) { ratio = other; next_order = q+1; }
        }
        if (next_order != q) successes = 0;
      }
      ratio = fmin(p->growth,fmax(p->shrink,ratio));
      if (rejected) ratio = fmin(1.0,ratio);
      /* A deadband keeps matrices reusable and avoids excessive step-ratio noise. */
      if (p->reuse && ratio >= 0.8 && ratio <= 1.3) ratio = 1.0;
      w->order = next_order; w->next_h = h*ratio; h = w->next_h; rejected = 0;
    } else {
      ++stats->base.rejected;
      if (!positive) ++stats->positivity_rejections;
      double ratio = isfinite(error) ? p->safety*pow(fmax(error,1e-12),-1.0/(q+1)) : 0.25;
      h *= fmax(p->shrink,fmin(0.5,ratio)); rejected = 1; successes = 0;
      if (q > 1) w->order = q-1;
    }
  }
  return CHEM_OK;
}
