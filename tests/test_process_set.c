#include "test.h"
#include "internal.h"
#include "chem_verify.h"
#include <float.h>

static void independent(const ck_mechanism *d, const ck_state *s, double *f, double *j)
{
  size_t n = ck_model_nspecies(s->model), cells = s->cells;
  for (size_t c = 0; c < cells; ++c)
    for (size_t r = 0; r < d->nreactions; ++r)
    {
      const ck_reaction *rx = d->reactions + r;
      double rate = s->rates[r * cells + c];
      for (size_t a = 0; a < rx->nreactants; ++a)
        if (!d->species[rx->reactants[a]].parameterized)
          rate *= s->y[rx->reactants[a] * cells + c];
      for (size_t a = 0; a < rx->nreactants; ++a)
        if (!d->species[rx->reactants[a]].parameterized && !ck_is_algebraic(s->model, rx->reactants[a]))
          f[rx->reactants[a] * cells + c] -= rate;
      for (size_t a = 0; a < rx->nproducts; ++a)
        if (!ck_is_algebraic(s->model, rx->products[a].species))
          f[rx->products[a].species * cells + c] += rx->products[a].yield * rate;
      for (size_t a = 0; a < rx->nreactants; ++a)
      {
        size_t col = rx->reactants[a];
        if (d->species[col].parameterized)
          continue;
        double partial = s->rates[r * cells + c];
        for (size_t b = 0; b < rx->nreactants; ++b)
          if (a != b && !d->species[rx->reactants[b]].parameterized)
            partial *= s->y[rx->reactants[b] * cells + c];
        for (size_t b = 0; b < rx->nreactants; ++b)
        {
          size_t row = rx->reactants[b];
          if (d->species[row].parameterized || ck_is_algebraic(s->model, row))
            continue;
          j[(row * n + col) * cells + c] += partial;
        }
        for (size_t b = 0; b < rx->nproducts; ++b)
        {
          size_t row = rx->products[b].species;
          if (ck_is_algebraic(s->model, row))
            continue;
          j[(row * n + col) * cells + c] -= rx->products[b].yield * partial;
        }
      }
    }
}
static int fixed_system(int mask, int difference)
{
  ck_species sp[] = {{.name = "foo"},
                     {.name = "bar"},
                     {.name = "baz"},
                     {.name = "quz"},
                     {.name = "quuz"},
                     {.name = "corge"},
                     {.name = "qux", .parameterized = 1}};
  size_t re[][2] = {{0, 2}, {1, 6}, {3, 0}, {2, 6}};
  ck_product pr[][2] = {{{1, 1}, {4, 2.4}}, {{0, 1}, {3, 1.4}}, {{0, 0}, {0, 0}}, {{1, 1}, {3, 2.5}}};
  ck_reaction rx[4] = {0};
  size_t nr = 4, n = 6;
  for (size_t r = 0; r < 4; ++r)
    rx[r] = (ck_reaction){.nreactants = r == 2 ? 1 : 2,
                          .reactants = re[r],
                          .nproducts = r == 2 ? 0 : 2,
                          .products = pr[r],
                          .rate = ck_rate_defaults(CK_ARRHENIUS)};
  if (mask)
  {
    n = 4;
    nr = 1;
    sp[0].name = "A";
    sp[1].name = "B";
    sp[2].name = "C";
    sp[3].name = "D";
    re[0][1] = 1;
    pr[0][0] = (ck_product){2, 1};
    pr[0][1] = (ck_product){3, 1};
  }
  if (difference)
  {
    n = 3;
    nr = 2;
    sp[0].name = "A";
    sp[1].name = "B";
    sp[2].name = "C";
    re[0][1] = 1;
    pr[0][0] = (ck_product){2, 2};
    rx[0].nproducts = 1;
    re[1][0] = 2;
    rx[1].nreactants = 1;
    pr[1][0] = (ck_product){0, 1};
    rx[1].nproducts = 1;
  }
  ck_mechanism d = {.name = "process set",
                    .nspecies = mask || difference ? n : 7,
                    .species = sp,
                    .nreactions = nr,
                    .reactions = rx};
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  CHECK(ck_state_create(m, 2, &s) == CK_OK);
  CHECK(m->n == n);
  if (!mask && !difference)
  {
    CHECK(ck_model_validate_species_used(&d) == CK_INVALID);
    for (size_t v = 0; v < 6; ++v)
    {
      size_t index;
      CHECK(ck_species_index(m, sp[v].name, &index) == CK_OK);
      CHECK(index == v);
    }
  }
  if (mask)
  {
    m->algebraic[1] = m->algebraic[3] = 1;
  }
  double values[] = {.1, .2, .3, .4, .5, 0, 1.1, 1.2, 1.3, 1.4, 1.5, 0};
  if (mask)
  {
    double v[] = {1, 2, 3, 4, 5, 6, 7, 8};
    memcpy(values, v, sizeof v);
  }
  if (difference)
  {
    double v[] = {1, 2, 3, .5, 4, 1.5};
    memcpy(values, v, sizeof v);
  }
  for (size_t c = 0; c < 2; ++c)
    for (size_t v = 0; v < n; ++v)
      s->y[v * 2 + c] = values[c * n + v];
  for (size_t c = 0; c < 2; ++c)
    for (size_t r = 0; r < nr; ++r)
    {
      s->rates[r * 2 + c] = mask         ? 10 + 10 * c
                            : difference ? (r ? 10 + 2 * c : 5 + 3 * c)
                                         : (r + 1) * 10 + 100 * c;
      if (!mask && !difference && (r == 1 || r == 3))
        s->rates[r * 2 + c] *= (70 + 10 * c) * .72;
    }
  size_t nnz = ck_jacobian_size(m);
  double *j = malloc(nnz * 2 * sizeof(double)), *f = malloc(n * 2 * sizeof(double));
  double *expected_j = malloc(n * n * 2 * sizeof(double)), *expected_f = malloc(n * 2 * sizeof(double));
  CHECK(j && f && expected_j && expected_f);
  double initial_j = mask ? 500 : 100;
  for (size_t i = 0; i < nnz * 2; ++i)
    j[i] = initial_j;
  for (size_t i = 0; i < n * n * 2; ++i)
    expected_j[i] = initial_j;
  for (size_t i = 0; i < n * 2; ++i)
    expected_f[i] = f[i] = 1000;
  ck_add_forcing(s, s->y, f);
  ck_subtract_jacobian(s, s->y, j);
  independent(&d, s, expected_f, expected_j);
  for (size_t i = 0; i < n * 2; ++i)
    NEAR(f[i], expected_f[i], 4 * (nextafter(fabs(expected_f[i]), INFINITY) - fabs(expected_f[i])));
  for (size_t row = 0; row < n; ++row)
    for (size_t col = 0; col < n; ++col)
    {
      size_t off = ck_jacobian_index(m, row, col);
      if (off != SIZE_MAX)
        for (size_t c = 0; c < 2; ++c)
          NEAR(j[off * 2 + c], expected_j[(row * n + col) * 2 + c],
               4 * (nextafter(fabs(expected_j[(row * n + col) * 2 + c]), INFINITY) -
                    fabs(expected_j[(row * n + col) * 2 + c])));
      else
        for (size_t c = 0; c < 2; ++c)
          CHECK(expected_j[(row * n + col) * 2 + c] == initial_j);
    }
  if (difference)
  {
    ck_jacobian(s, s->y, j);
    double plus[6], minus[6];
    for (size_t col = 0; col < n; ++col)
      for (size_t c = 0; c < 2; ++c)
      {
        size_t index = col * 2 + c;
        double base = s->y[index], h = 1e-8 * fmax(1, fabs(base));
        s->y[index] = base + h;
        ck_forcing(s, s->y, plus);
        s->y[index] = base - h;
        ck_forcing(s, s->y, minus);
        s->y[index] = base;
        for (size_t row = 0; row < n; ++row)
        {
          double fd = (plus[row * 2 + c] - minus[row * 2 + c]) / (2 * h);
          size_t off = ck_jacobian_index(m, row, col);
          if (off == SIZE_MAX)
            CHECK(fabs(fd) <= 1e-6);
          else
            NEAR(-j[off * 2 + c], fd, 1e-7 + 1e-7 * fmax(fabs(fd), fabs(j[off * 2 + c])));
        }
      }
  }
  free(j);
  free(f);
  free(expected_j);
  free(expected_f);
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
static uint32_t random_state = 78513;
static uint32_t random_value(void)
{
  random_state = 1664525 * random_state + 1013904223;
  return random_state;
}
static int random_system(size_t cells, size_t nr, size_t n)
{
  ck_species *sp = calloc(n, sizeof *sp);
  char (*names)[24] = calloc(n, sizeof *names);
  ck_reaction *rx = calloc(nr, sizeof *rx);
  size_t (*reactants)[3] = calloc(nr, sizeof *reactants);
  ck_product(*products)[10] = calloc(nr, sizeof *products);
  CHECK(sp && names && rx && reactants && products);
  for (size_t v = 0; v < n; ++v)
  {
    snprintf(names[v], 24, "%zu", v);
    sp[v].name = names[v];
  }
  for (size_t r = 0; r < nr; ++r)
  {
    rx[r] = (ck_reaction){.nreactants = random_value() % 4,
                          .reactants = reactants[r],
                          .nproducts = random_value() % 11,
                          .products = products[r],
                          .rate = ck_rate_defaults(CK_ARRHENIUS)};
    for (size_t a = 0; a < rx[r].nreactants; ++a)
      reactants[r][a] = random_value() % n;
    for (size_t a = 0; a < rx[r].nproducts; ++a)
      products[r][a] = (ck_product){random_value() % n, 1.2};
  }
  ck_mechanism d = {.name = "random", .nspecies = n, .species = sp, .nreactions = nr, .reactions = rx};
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  CHECK(ck_state_create(m, cells, &s) == CK_OK);
  for (size_t i = 0; i < n * cells; ++i)
    s->y[i] = exp(-2 + 4 * ((double)(random_value() % 10000) / 10000 - .5));
  for (size_t i = 0; i < nr * cells; ++i)
    s->rates[i] = exp(-2 + 4 * ((double)(random_value() % 10000) / 10000 - .5));
  double *f = malloc(n * cells * sizeof(double)), *ef = malloc(n * cells * sizeof(double)),
         *ej = calloc(n * n * cells, sizeof(double));
  CHECK(f && ef && ej);
  for (size_t i = 0; i < n * cells; ++i)
    f[i] = ef[i] = 1000;
  ck_add_forcing(s, s->y, f);
  independent(&d, s, ef, ej);
  for (size_t i = 0; i < n * cells; ++i)
    NEAR(f[i], ef[i], 4 * (nextafter(fabs(ef[i]), INFINITY) - fabs(ef[i])));
  free(f);
  free(ef);
  free(ej);
  ck_state_free(s);
  ck_model_free(m);
  free(sp);
  free(names);
  free(rx);
  free(reactants);
  free(products);
  return 0;
}
int test_process_set(const char *suite, const char *test)
{
  (void)test;
  if (!strcmp(suite, "ProcessSet"))
    return fixed_system(0, 0);
  if (!strcmp(suite, "ProcessSetAlgebraicVariables"))
    return fixed_system(1, 0);
  if (!strcmp(suite, "ProcessSetFiniteDifferenceJacobian"))
    return fixed_system(0, 1);
  if (!strcmp(suite, "RandomProcessSet"))
  {
    CHECK(random_system(200, 50, 40) == 0);
    CHECK(random_system(300, 30, 20) == 0);
    return random_system(400, 100, 80);
  }
  return 2;
}
