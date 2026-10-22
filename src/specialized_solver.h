/* Three-stage Rosenbrock from MICM (Apache-2.0; see NOTICE).
   Included by generated C with SP_PREFIX, SP_N, SP_NNZ and sp_diagonal.
   Only stage arithmetic is fused; the error reduction and controller retain
   the reference order, initial step, rejection rules and stopping conditions. */
#include <float.h>
#define SP_CAT_(a, b) a##b
#define SP_CAT(a, b) SP_CAT_(a, b)
#define SP_FN(suffix) SP_CAT(SP_PREFIX, suffix)

static ck_status SP_FN(_solve)(ck_state *s, double dt, ck_stats *stats)
{
  if (!s || !stats || !(dt >= 0) || !isfinite(dt))
    return CK_INVALID;
  memset(stats, 0, sizeof *stats);
  if (dt == 0)
  {
    ck_post_solve_clamp(s);
    return CK_OK;
  }
  ck_workspace *w = s->workspace;
  const size_t cells = s->cells, size = SP_N * cells;
  const double h_min = 1e-15 * dt;
  double h = fmin(fmax(h_min, fabs(1e-6 * dt)), fabs(dt)), time = 0;
  int reject_last = 0, reject_more = 0;
  ck_status status = CK_OK;
  while ((time - dt + DBL_EPSILON) <= 0)
  {
    if (stats->internal_steps > 1000)
    {
      status = CK_MAX_STEPS;
      break;
    }
    if (time + 0.1 * h == time || h <= fmin(DBL_EPSILON, h_min))
    {
      status = CK_STEP_SMALL;
      break;
    }
    h = fmin(h, fabs(dt - time));
#if SP_FUSE
    SP_FN(_both)(s, s->y, w->initial, w->jacobian);
#else
    SP_FN(_forcing)(s, s->y, w->initial);
    SP_FN(_jacobian)(s, s->y, w->jacobian);
#endif
    ++stats->function_calls;
    ++stats->jacobian_updates;
    int accepted = 0;
    while (!accepted)
    {
      double alpha = 1.0 / (h * 0.43586652150845899941601945119356);
      for (size_t j = 0; j < SP_N; ++j)
      {
        double *restrict diagonal = w->jacobian + sp_diagonal[j] * cells;
#pragma omp simd
        for (size_t c = 0; c < cells; ++c)
          diagonal[c] += alpha;
      }
      SP_FN(_factor)(cells, w->jacobian);
      ++stats->decompositions;
      double *restrict k0 = w->stage[0], *restrict k1 = w->stage[1], *restrict k2 = w->stage[2];
      double *restrict next = w->next;
      const double *restrict y = s->y;
      memcpy(k0, w->initial, size * sizeof(double));
      SP_FN(_linear)(cells, w->jacobian, k0);
#pragma omp simd
      for (size_t i = 0; i < size; ++i)
        next[i] = y[i] + 1.0 * k0[i];
      SP_FN(_forcing)(s, next, k1);
      ++stats->function_calls;
      double c10 = -0.10156171083877702091975600115545e+01 / h;
#pragma omp simd
      for (size_t i = 0; i < size; ++i)
      {
        k2[i] = k1[i];
        k1[i] += c10 * k0[i];
      }
      SP_FN(_linear)(cells, w->jacobian, k1);
      double c20 = 0.40759956452537699824805835358067e+01 / h;
      double c21 = 0.92076794298330791242156818474003e+01 / h;
#pragma omp simd
      for (size_t i = 0; i < size; ++i)
      {
        k2[i] += c20 * k0[i];
        k2[i] += c21 * k1[i];
      }
      SP_FN(_linear)(cells, w->jacobian, k2);
      stats->solves += 3;
      double *restrict errors = w->error;
#pragma omp simd
      for (size_t i = 0; i < size; ++i)
      {
        double value = y[i];
        value += 0.1e+01 * k0[i];
        value += 0.61697947043828245592553615689730e+01 * k1[i];
        value += -0.42772256543218573326238373806514 * k2[i];
        next[i] = value;
        double error = 0;
        error += 0.5 * k0[i];
        error += -0.29079558716805469821718236208017e+01 * k1[i];
        error += 0.22354069897811569627360909276199 * k2[i];
        errors[i] = error;
      }
      double error = ck_normalized_error(s, y, next, errors);
      double factor = fmin(6.0, fmax(0.2, 0.9 / pow(error, 1.0 / 3.0)));
      double next_h = h * factor;
      ++stats->internal_steps;
      if (isnan(error))
      {
        status = CK_NAN;
        goto done;
      }
      if (isinf(error))
      {
        status = CK_INF;
        goto done;
      }
      if (error < 1 || h < h_min)
      {
        ++stats->accepted;
        time += h;
        double *swap = s->y;
        s->y = w->next;
        w->next = swap;
        next_h = fmax(h_min, fmin(next_h, dt));
        if (reject_last)
          next_h = fmin(next_h, h);
        reject_last = reject_more = 0;
        h = next_h;
        accepted = 1;
      }
      else
      {
        if (reject_more)
          next_h = h * 0.1;
        reject_more = reject_last;
        reject_last = 1;
        h = next_h;
        if (stats->accepted >= 1)
          ++stats->rejected;
        SP_FN(_jacobian)(s, s->y, w->jacobian);
        ++stats->jacobian_updates;
      }
    }
  }
done:
  stats->final_time = time;
  double *restrict y = s->y;
#pragma omp simd
  for (size_t i = 0; i < size; ++i)
    y[i] = y[i] > 0 ? y[i] : 0;
  return status;
}

#undef SP_FN
#undef SP_CAT
#undef SP_CAT_
#undef SP_PREFIX
#undef SP_N
#undef SP_NNZ
#undef SP_FUSE
