#include "adaptive_internal.h"

ChemAdaptive chem_adaptive_parameters(int method) {
  ChemAdaptive p = {0};
  p.method = method; p.norm = 1; p.max_order = 5;
  p.jacobian_age = 1; p.max_steps = 200000;
  p.h_start = 1.0; p.safety = 0.85; p.growth = 3.0; p.shrink = 0.1;
  p.positivity = 1.0; p.reuse = 1; p.reuse_ratio = 0.2;
  p.step_hold = 1.0;
  p.controller_b = 1.0; p.controller_k = 2.0;
  return p;
}
ChemWorkspace *chem_workspace_create(const ChemState *s) {
  if (!s || s->model->nconstraints || s->model->nvariables > SIZE_MAX/s->stride/(16*sizeof(double))) return NULL;
  ChemWorkspace *w = calloc(1, sizeof *w);
  if (!w) return NULL;
  size_t n = s->model->nvariables*s->stride;
  w->owner = s;
  w->jacobian = malloc(s->model->nnz*s->stride*sizeof(double));
  w->history = calloc(7*n, sizeof(double)); w->scratch = calloc(9*n, sizeof(double));
  if (!w->jacobian || !w->history || !w->scratch) { chem_workspace_free(w); return NULL; }
  w->previous_error = 1.0; w->order = 1;
  return w;
}
void chem_workspace_free(ChemWorkspace *w) {
  if (w) { free(w->jacobian); free(w->history); free(w->scratch); free(w); }
}
void chem_workspace_reset(ChemWorkspace *w) {
  if (w) {
    w->history_count = 0; w->factor_valid = 0; w->jacobian_age = 0; w->order = 1;
    w->previous_error = 1.0; w->previous_h = 0.0;
  }
}
double chem_adaptive_norm(const ChemAdaptive *p, const ChemState *s, const double *y,
                          const double *yn, const double *e) {
  double total = 0.0, largest = 0.0;
  for (size_t c = 0; c < s->cells; ++c) {
    double sum = 0.0;
    for (size_t i = 0; i < s->model->nvariables; ++i) {
      size_t a = i*s->stride+c;
      double scale = s->absolute_tolerance[i]+s->relative_tolerance*fmax(fabs(y[a]),fabs(yn[a]));
      double z = scale > 0.0 ? fabs(e[a])/scale : (e[a] == 0.0 ? 0.0 : INFINITY);
      if (!isfinite(z)) return INFINITY;
      if (p->norm == 2) { if (z > largest) largest = z; }
      else sum += z*z;
    }
    if (p->norm != 2) {
      total += sum;
      if (sum > largest) largest = sum;
    }
  }
  if (p->norm == 2) return largest;
  return sqrt((p->norm ? largest : total/s->cells)/s->model->nvariables);
}
int chem_adaptive_positive(const ChemAdaptive *p, const ChemState *s, const double *y, const double *yn) {
  for (size_t i = 0; i < s->model->nvariables; ++i) for (size_t c = 0; c < s->cells; ++c) {
    size_t a = i*s->stride+c;
    double scale = s->absolute_tolerance[i]+s->relative_tolerance*fabs(y[a]);
    if (!isfinite(yn[a]) || yn[a] < -p->positivity*scale) return 0;
  }
  return 1;
}
static ChemStatus rosenbrock(const ChemAdaptive *p, double dt, ChemState *s, ChemWorkspace *w, ChemAdaptiveStats *stats) {
  const ChemModel *m = s->model; ChemRosenbrock coefficients = chem_rosenbrock_parameters(p->method);
  const ChemRosenbrock *r = &coefficients;
  size_t n = m->nvariables*s->stride, nj = m->nnz*s->stride;
  double t = 0.0, h = w->next_h > 0.0 ? w->next_h : p->h_start;
  double hmax = p->h_max > 0.0 ? fmin(p->h_max,dt) : dt;
  int fresh = 1, rejected = 0;
  w->factor_valid = 0;
  while (t < dt) {
    if (stats->base.internal_steps >= p->max_steps) return CHEM_MAX_STEPS;
    h = fmin(h,fmin(hmax,dt-t));
    if (!(h > 0.0) || t+h == t) return CHEM_SMALL_STEP;
    if (fresh) {
      chem_forcing(s,s->variables,s->forcing); ++stats->base.function_calls;
      if (!w->jacobian_age || w->jacobian_age >= p->jacobian_age) {
        chem_jacobian(s,s->variables,w->jacobian); ++stats->base.jacobian_updates;
        w->jacobian_age = 0; w->factor_valid = 0;
      }
      fresh = 0;
    }
    if (!p->reuse || !w->factor_valid || w->factor_h != h) {
      memcpy(s->jacobian,w->jacobian,nj*sizeof(double));
      chem_alpha_minus_jacobian(s,1.0/(h*r->gamma[0]));
      chem_lu_mozart(m,s->cells,s->stride,s->jacobian); ++stats->base.decompositions;
      w->factor_h = h; w->factor_valid = 1;
    } else ++stats->matrix_reuses;
    for (int stage = 0; stage < r->stages; ++stage) {
      size_t off = (size_t)(stage*(stage-1)/2);
      if (!stage) memcpy(s->stages[stage],s->forcing,n*sizeof(double));
      else if (r->new_function[stage]) {
        memcpy(s->new_y,s->variables,n*sizeof(double));
        for (int j = 0; j < stage; ++j) chem_axpy(n,s->new_y,r->a[off+j],s->stages[j]);
        chem_forcing(s,s->new_y,s->stages[stage]); ++stats->base.function_calls;
      }
      if (stage+1 < r->stages && !r->new_function[stage+1]) memcpy(s->stages[stage+1],s->stages[stage],n*sizeof(double));
      for (int j = 0; j < stage; ++j) chem_axpy(n,s->stages[stage],r->c[off+j]/h,s->stages[j]);
      chem_linear_solve(m,s->cells,s->stride,s->jacobian,s->stages[stage]); ++stats->base.solves;
    }
    memcpy(s->new_y,s->variables,n*sizeof(double)); memset(s->error,0,n*sizeof(double));
    for (int j = 0; j < r->stages; ++j) {
      chem_axpy(n,s->new_y,r->m[j],s->stages[j]); chem_axpy(n,s->error,r->e[j],s->stages[j]);
    }
    ++stats->base.internal_steps;
    double error = chem_adaptive_norm(p,s,s->variables,s->new_y,s->error);
    int positive = chem_adaptive_positive(p,s,s->variables,s->new_y);
    double factor = p->shrink;
    if (isfinite(error)) {
      if (p->controller == 2 && !rejected && w->previous_h > 0.0) {
        double exponent = 1.0/(p->controller_b*p->controller_k);
        factor = p->safety*pow(fmax(error,1e-12),-exponent)*pow(w->previous_error,-exponent)*
          pow(h/w->previous_h,-1.0/p->controller_b);
      } else if (p->controller == 1 && !rejected)
        factor = p->safety*pow(fmax(error,1e-12),-0.7/r->order)*pow(w->previous_error,0.4/r->order);
      else factor = p->safety*pow(fmax(error,1e-12),-1.0/r->order);
      factor = fmin(p->growth,fmax(p->shrink,factor));
    }
    if (positive && error <= 1.0) {
      for (size_t a = 0; a < n; ++a) {
        if (s->new_y[a] < 0.0) {
          ++stats->clipped_values;
          stats->largest_clip = fmax(stats->largest_clip,-s->new_y[a]);
        }
        s->variables[a] = fmax(0.0,s->new_y[a]);
      }
      t = h == dt-t ? dt : t+h; stats->base.final_time = t;
      ++stats->base.accepted; ++w->jacobian_age; fresh = 1;
      w->previous_error = fmax(error,1e-12);
      w->previous_h = h;
      if (rejected) factor = fmin(1.0,factor);
      if (factor >= 1.0/p->step_hold && factor <= p->step_hold) factor = 1.0;
      w->next_h = h*factor; h = w->next_h; rejected = 0;
    } else {
      ++stats->base.rejected;
      if (!positive) { ++stats->positivity_rejections; factor = fmin(0.5,factor); }
      h *= fmin(0.8,factor); rejected = 1;
      /* Refresh an approximate Jacobian before persisting with a failed trial. */
      if (w->jacobian_age) { w->jacobian_age = 0; fresh = 1; }
      w->factor_valid = 0;
    }
  }
  return CHEM_OK;
}
ChemStatus chem_adaptive_solve(const ChemAdaptive *p, double dt, ChemState *s, ChemWorkspace *w, ChemAdaptiveStats *stats) {
  if (!p || !s || !w || w->owner != s || !stats || s->model->nconstraints || !isfinite(dt) || dt < 0.0 ||
      !((p->method >= 2 && p->method <= 6) || p->method == 10) || p->controller < 0 || p->controller > 2 ||
      p->norm < 0 || p->norm > 2 || p->max_order < 1 || p->max_order > 5 || !p->max_steps || !p->jacobian_age ||
      !(p->h_start > 0.0) || !isfinite(p->h_start) || !(p->h_max >= 0.0) || !isfinite(p->h_max) ||
      !(p->safety > 0.0) || !isfinite(p->safety) || !(p->growth >= 1.0) || !isfinite(p->growth) ||
      !(p->shrink > 0.0 && p->shrink < 1.0) || !(p->positivity >= 0.0) || !isfinite(p->positivity) ||
      !(p->reuse_ratio >= 0.0 && p->reuse_ratio < 1.0) || !(p->step_hold >= 1.0) || !isfinite(p->step_hold) ||
      !(p->controller_b >= 1.0) || !isfinite(p->controller_b) || !(p->controller_k > 0.0) || !isfinite(p->controller_k)) return CHEM_INVALID;
  memset(stats,0,sizeof *stats);
  if (dt == 0.0) return CHEM_OK;
  if (p->method == 10) return chem_bdf_solve(p,dt,s,w,stats);
  w->jacobian_age = 0;
  return rosenbrock(p,dt,s,w,stats);
}
