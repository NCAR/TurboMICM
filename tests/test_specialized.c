#include "internal.h"
#include "chem_specialized.h"
#include <stdio.h>

extern const ck_specialized ts1_specialized, sp_tutorial_no, sp_tutorial_user, sp_tutorial_robertson,
    sp_extended, sp_emission, sp_alternate;

static int compare(const char *what, const double *reference, const double *actual, size_t size)
{
  for (size_t i = 0; i < size; ++i)
    if (!isfinite(actual[i]) || !isfinite(reference[i]) ||
        fabs(actual[i] - reference[i]) > 1e-30 + 1e-10 * fabs(reference[i]))
    {
      fprintf(stderr, "%s[%zu]: reference %.17g, generated %.17g\n", what, i, reference[i], actual[i]);
      return 0;
    }
  return 1;
}
#define REQUIRE(condition)                                                                                   \
  do                                                                                                         \
  {                                                                                                          \
    if (!(condition))                                                                                        \
    {                                                                                                        \
      fprintf(stderr, "line %d: %s\n", __LINE__, #condition);                                                \
      return 0;                                                                                              \
    }                                                                                                        \
  } while (0)

static int check(const ck_specialized *sp, size_t cells, int tight)
{
  ck_model *m = NULL;
  ck_state *a = NULL, *b = NULL;
  REQUIRE(ck_model_create(sp->mechanism, &m) == CK_OK);
  REQUIRE(ck_state_create(m, cells, &a) == CK_OK);
  int ts1_case = sp == &ts1_specialized;
  for (size_t c = 0; c < cells; ++c)
  {
    a->conditions[c] = (ck_conditions){270 + 0.7 * (c % 17), 100000 + 100 * c, 0};
    ck_ideal_air_density(&a->conditions[c]);
    a->conditions[c].air_density *= 1.001;
    for (size_t v = 0; v < m->n; ++v)
      a->y[v * cells + c] =
          ts1_case ? a->y[v * cells + c] * (0.8 + 0.01 * c) : 1e-12 * (v + 1) * (1 + 0.03 * c);
    for (size_t j = 0; j < m->nc; ++j)
    {
      const char *name = ck_custom_name(m, j);
      double value = a->custom[j * cells + c];
      if (ts1_case && (!strcmp(name, "jo2_a") || !strcmp(name, "jh2o_b") || !strcmp(name, "jh2o_c") ||
                       !strcmp(name, "jch4_a") || !strcmp(name, "jch4_b") || !strcmp(name, "jhf")))
        value = 0;
      else if (strstr(name, "effective radius"))
        value = 1e-7;
      else if (strstr(name, "particle number"))
        value = 2.5e6;
      else if (!ts1_case)
        value = strstr(name, "emission") ? 1e-20 : 1e-12;
      a->custom[j * cells + c] = value * (0.9 + 0.01 * c);
    }
    if (tight)
      for (size_t v = 0; v < m->n; ++v)
        a->absolute_tolerance[v] = 1e-15;
  }
  REQUIRE(ck_state_clone(a, &b) == CK_OK);
  REQUIRE(ck_update_rates(a) == CK_OK && sp->update_rates(b) == CK_OK);
  REQUIRE(compare("rates", a->rates, b->rates, m->nr * cells));
  size_t ny = m->n * cells, na = m->nnz * cells;
  double *fa = calloc(ny, sizeof(double)), *fb = calloc(ny, sizeof(double));
  double *ja = calloc(na, sizeof(double)), *jb = calloc(na, sizeof(double));
  REQUIRE(fa && fb && ja && jb);
  ck_forcing(a, a->y, fa);
  sp->forcing(b, b->y, fb);
  REQUIRE(compare("forcing", fa, fb, ny));
  ck_jacobian(a, a->y, ja);
  sp->jacobian(b, b->y, jb);
  REQUIRE(compare("jacobian", ja, jb, na));
  ck_alpha_minus_jacobian(m, cells, ja, 1e4);
  ck_alpha_minus_jacobian(m, cells, jb, 1e4);
  ck_factor(m, cells, ja, CK_LU_MOZART);
  sp->factor(cells, jb);
  REQUIRE(compare("factorization", ja, jb, na));
  ck_linear_solve(m, cells, ja, fa);
  sp->linear_solve(cells, jb, fb);
  REQUIRE(compare("linear solve", fa, fb, ny));
  ck_rosenbrock_parameters p;
  ck_rosenbrock_parameters_init(3, &p);
  for (size_t step = 0; step < 6; ++step)
  {
    /* Change every runtime input between steps, including darkness. */
    for (size_t j = 0; j < m->nc * cells; ++j)
      a->custom[j] = b->custom[j] = a->custom[j] * (step == 3 ? 0 : 1.03);
    for (size_t c = 0; c < cells; ++c)
    {
      a->conditions[c].temperature += 0.1;
      a->conditions[c].pressure += 10;
      a->conditions[c].air_density += 0.01;
      b->conditions[c] = a->conditions[c];
    }
    REQUIRE(ck_update_rates(a) == CK_OK && sp->update_rates(b) == CK_OK);
    REQUIRE(compare("updated rates", a->rates, b->rates, m->nr * cells));
    ck_stats sa, sb;
    double dt = ts1_case ? (tight ? 0.01 : 30) : 0.1;
    ck_status status = ck_solve(a, dt, &p, &sa);
    if (status)
      fprintf(stderr, "generic step %zu: %s\n", step, ck_status_string(status));
    REQUIRE(status == CK_OK);
    REQUIRE(sp->solve(b, dt, &sb) == CK_OK);
    REQUIRE(sa.internal_steps == sb.internal_steps && sa.accepted == sb.accepted &&
            sa.rejected == sb.rejected && sa.jacobian_updates == sb.jacobian_updates &&
            sa.function_calls == sb.function_calls && sa.solves == sb.solves &&
            sa.decompositions == sb.decompositions && sa.final_time == sb.final_time);
    REQUIRE(compare("trajectory", a->y, b->y, ny));
  }
  ck_stats zero;
  REQUIRE(sp->solve(b, 0, &zero) == CK_OK && zero.internal_steps == 0);
  REQUIRE(sp->solve(b, -1, &zero) == CK_INVALID);
  free(fa);
  free(fb);
  free(ja);
  free(jb);
  ck_state_free(a);
  ck_state_free(b);
  ck_model_free(m);
  return 1;
}

int main(void)
{
  const ck_specialized *cases[] = {&ts1_specialized,       &sp_tutorial_no, &sp_tutorial_user,
                                   &sp_tutorial_robertson, &sp_extended,    &sp_emission,
                                   &sp_alternate};
  const size_t widths[] = {1, 3, 32, 37, 64};
  for (size_t i = 0; i < sizeof cases / sizeof *cases; ++i)
    for (size_t w = 0; w < sizeof widths / sizeof *widths; ++w)
      for (int tight = 0; tight < 2; ++tight)
      {
        if (!check(cases[i], widths[w], tight))
        {
          fprintf(stderr, "FAIL %s cells=%zu tight=%d\n", cases[i]->mechanism->name, widths[w], tight);
          return 1;
        }
        printf("PASS %s cells=%zu tight=%d\n", cases[i]->mechanism->name, widths[w], tight);
      }
  return 0;
}
