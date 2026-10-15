/* External process and constraint integration cases, including composed models. */
#include "test.h"

typedef struct
{
  size_t n, nrows;
  const char *names[4], *row_names[3];
  double coefficients[3][4], constant[3], rate;
  int process, solvent, temperature;
  size_t indices[4], rows[3], jacobian[3][4], parameter;
} affine_model;
static void *affine_clone(const void *source)
{
  affine_model *copy = malloc(sizeof *copy);
  if (copy)
    *copy = *(const affine_model *)source;
  return copy;
}
static ck_status affine_finalize(const ck_model *m, void *context)
{
  affine_model *a = context;
  for (size_t j = 0; j < a->n; ++j)
    if (ck_species_index(m, a->names[j], a->indices + j))
      return CK_UNKNOWN;
  for (size_t r = 0; r < a->nrows; ++r)
  {
    if (ck_species_index(m, a->row_names[r], a->rows + r))
      return CK_UNKNOWN;
    for (size_t j = 0; j < a->n; ++j)
    {
      a->jacobian[r][j] = ck_jacobian_index(m, a->rows[r], a->indices[j]);
      if (a->coefficients[r][j] && a->jacobian[r][j] == SIZE_MAX)
        return CK_INVALID;
    }
  }
  if (a->temperature && ck_custom_index(m, "equilibrium_K", &a->parameter))
    return CK_UNKNOWN;
  return CK_OK;
}
static ck_status affine_update(ck_state *s, void *context)
{
  affine_model *a = context;
  if (a->temperature)
    for (size_t c = 0; c < s->cells; ++c)
      s->custom[a->parameter * s->cells + c] =
          2 * exp(3000 * (1 / 298.15 - 1 / s->conditions[c].temperature));
  return CK_OK;
}
static void affine_residual(const ck_state *s, const double *y, double *f, void *context)
{
  const affine_model *a = context;
  for (size_t r = 0; r < a->nrows; ++r)
    for (size_t c = 0; c < s->cells; ++c)
    {
      double value = 0;
      for (size_t j = 0; j < a->n; ++j)
        value += (a->temperature && j == 0 ? s->custom[a->parameter * s->cells + c] : a->coefficients[r][j]) *
                 y[a->indices[j] * s->cells + c];
      f[a->rows[r] * s->cells + c] = value + a->constant[r];
    }
}
static void affine_constraint_jacobian(const ck_state *s, const double *y, double *jac, void *context)
{
  (void)y;
  const affine_model *a = context;
  for (size_t r = 0; r < a->nrows; ++r)
    for (size_t j = 0; j < a->n; ++j)
      if (a->coefficients[r][j])
        for (size_t c = 0; c < s->cells; ++c)
          jac[a->jacobian[r][j] * s->cells + c] -=
              a->temperature && j == 0 ? s->custom[a->parameter * s->cells + c] : a->coefficients[r][j];
}
static void affine_process(const ck_state *s, const double *y, double *f, void *context)
{
  const affine_model *a = context;
  if (!a->process)
    return;
  for (size_t c = 0; c < s->cells; ++c)
  {
    double rate = a->rate * y[a->indices[0] * s->cells + c];
    if (a->solvent)
      rate *= y[a->indices[2] * s->cells + c];
    f[a->indices[0] * s->cells + c] -= rate;
    f[a->indices[1] * s->cells + c] += rate;
  }
}
static void affine_process_jacobian(const ck_state *s, const double *y, double *jac, void *context)
{
  const affine_model *a = context;
  if (!a->process)
    return;
  for (size_t c = 0; c < s->cells; ++c)
  {
    double derivative = a->rate * (a->solvent ? y[a->indices[2] * s->cells + c] : 1);
    jac[ck_jacobian_index(s->model, a->indices[0], a->indices[0]) * s->cells + c] += derivative;
    jac[ck_jacobian_index(s->model, a->indices[1], a->indices[0]) * s->cells + c] -= derivative;
    if (a->solvent)
    {
      derivative = a->rate * y[a->indices[0] * s->cells + c];
      jac[ck_jacobian_index(s->model, a->indices[0], a->indices[2]) * s->cells + c] += derivative;
      jac[ck_jacobian_index(s->model, a->indices[1], a->indices[2]) * s->cells + c] -= derivative;
    }
  }
}
static ck_external_model external_descriptor(affine_model *a, ck_named_jacobian_entry *pattern)
{
  size_t count = 0;
  for (size_t r = 0; r < a->nrows; ++r)
    for (size_t j = 0; j < a->n; ++j)
      if (a->coefficients[r][j])
        pattern[count++] = (ck_named_jacobian_entry){a->row_names[r], a->names[j], 1};
  if (a->process)
    for (size_t r = 0; r < 2; ++r)
    {
      pattern[count++] = (ck_named_jacobian_entry){a->names[r], a->names[0], 0};
      if (a->solvent)
        pattern[count++] = (ck_named_jacobian_entry){a->names[r], a->names[2], 0};
    }
  static const char *const custom[] = {"equilibrium_K"};
  return (ck_external_model){.name = "affine external model",
                             .ncustom = a->temperature ? 1 : 0,
                             .custom_names = custom,
                             .nalgebraic = a->nrows,
                             .algebraic_names = a->row_names,
                             .njacobian = count,
                             .jacobian = pattern,
                             .context = a,
                             .clone_context = affine_clone,
                             .destroy_context = free,
                             .finalize = affine_finalize,
                             .update_constraints = affine_update,
                             .forcing = affine_process,
                             .jacobian_terms = affine_process_jacobian,
                             .constraint_residual = affine_residual,
                             .constraint_jacobian = affine_constraint_jacobian};
}
static int external_mass(const char *test)
{
  int processonly = IS("AddExternalModelProcessOnly") || IS("AddExternalModelOnlyStandardRosenbrock") ||
                    IS("ProcessForcingJacobian") ||
                    IS("ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering");
  int solvent = IS("ProcessJacobianElementInAlgebraicRowSurvivesFiltering"),
      combined = IS("CombinedBuiltInAndExternalConstraints"),
      builtinrow = IS("ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering");
  double k = IS("AlgebraicRowHoldsOnlyConstraintTerms")                                           ? 10
             : IS("ProcessForcingJacobian")                                                       ? 3.5
             : IS("DAESolveEnforcesConservation") || IS("MultiGridCell") || solvent || builtinrow ? .1
                                                                                                  : .01;
  affine_model a = {.n = solvent ? 3 : 2,
                    .names = {"A_GAS", "AEROSOL.A_AQ", "S"},
                    .nrows = processonly ? 0 : 1,
                    .row_names = {"AEROSOL.A_AQ"},
                    .coefficients = {{1, 1, 0}},
                    .constant = {-1},
                    .rate = k,
                    .process = 1,
                    .solvent = solvent};
  ck_named_jacobian_entry pattern[16];
  ck_external_model e = external_descriptor(&a, pattern);
  ck_species extra = {.name = "AEROSOL.A_AQ"};
  e.nspecies = 1;
  e.species = &extra;
  ck_species gas[] = {{.name = "A_GAS"}, {.name = solvent ? "S" : "B"}, {.name = "C"}};
  ck_mechanism d = {.name = "external mass",
                    .nspecies = combined                ? 3
                                : solvent || builtinrow ? 2
                                                        : 1,
                    .species = gas,
                    .nexternal = 1,
                    .external = &e,
                    .temperature = 298,
                    .pressure = 101325};
  ck_product reactant = {1, 1}, product = {builtinrow ? 2 : 2, 1};
  ck_constraint constraint = {.type = CK_EQUILIBRIUM_CONSTRAINT,
                              .name = "B_C_eq",
                              .algebraic_species = 2,
                              .nreactants = 1,
                              .reactants = &reactant,
                              .nproducts = 1,
                              .products = &product,
                              .equilibrium_ref = 5};
  size_t react = 0;
  ck_product prod = {1, 1};
  ck_reaction rx = {.nreactants = 1,
                    .reactants = &react,
                    .nproducts = 1,
                    .products = &prod,
                    .rate = ck_rate_defaults(CK_ARRHENIUS)};
  rx.rate.A = .05;
  if (combined || builtinrow)
  {
    d.nconstraints = 1;
    d.constraints = &constraint;
  }
  if (combined)
  {
    d.nreactions = 1;
    d.reactions = &rx;
  }
  ck_model *m;
  ck_state *s;
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  size_t cells =
      IS("MultiGridCell") || IS("ProcessForcingJacobian") || IS("ConstraintResidualJacobian") ? 3 : 1;
  if (IS("ProcessForcingJacobian") || IS("ConstraintResidualJacobian"))
    cells = 2;
  CHECK(ck_state_create(m, cells, &s) == CK_OK);
  size_t g, q;
  CHECK(ck_species_index(m, "A_GAS", &g) == CK_OK);
  CHECK(ck_species_index(m, "AEROSOL.A_AQ", &q) == CK_OK);
  CHECK(ck_model_nconstraints(m) == (size_t)(combined ? 2 : processonly && !builtinrow ? 0 : 1));
  CHECK(!ck_is_algebraic(m, g));
  CHECK(ck_is_algebraic(m, q) == (!processonly || builtinrow));
  if (builtinrow)
    CHECK(ck_jacobian_index(m, q, g) != SIZE_MAX);
  if (combined)
  {
    size_t c;
    CHECK(ck_species_index(m, "C", &c) == CK_OK);
    CHECK(ck_is_algebraic(m, c));
  }
  for (size_t c = 0; c < cells; ++c)
  {
    s->y[g * cells + c] = IS("AddExternalModelConstraintsOnly") ? .8 : .9 - .2 * c;
    s->y[q * cells + c] = 1 - s->y[g * cells + c];
    if (solvent)
    {
      size_t index;
      CHECK(ck_species_index(m, "S", &index) == CK_OK);
      s->y[index] = 1;
    }
  }
  ck_rosenbrock_parameters p;
  CHECK(ck_rosenbrock_parameters_init(43, &p) == CK_OK);
  ck_stats stats;
  if (IS("AlgebraicRowHoldsOnlyConstraintTerms"))
  {
    s->y[g] = .7;
    s->y[q] = .2;
    double f[2];
    double *j = calloc(ck_jacobian_size(m), sizeof(double));
    CHECK(j);
    ck_forcing(s, s->y, f);
    ck_constraint_forcing(s, s->y, f);
    ck_jacobian(s, s->y, j);
    ck_constraint_jacobian(s, s->y, j);
    NEAR(f[g], -7, 1e-12);
    NEAR(f[q], -.1, 1e-12);
    NEAR(j[ck_jacobian_index(m, g, g)], 10, 1e-12);
    NEAR(j[ck_jacobian_index(m, q, g)], -1, 1e-12);
    NEAR(j[ck_jacobian_index(m, q, q)], -1, 1e-12);
    free(j);
  }
  else if (IS("DAESolveEnforcesConservation") || IS("MultiGridCell") ||
           IS("AddExternalModelConstraintsOnly") || solvent)
  {
    size_t steps = IS("DAESolveEnforcesConservation") ? 20 : solvent ? 10 : 1;
    double dt = steps == 1 ? 50 : 10;
    for (size_t t = 0; t < steps; ++t)
    {
      CHECK(ck_solve(s, dt, &p, &stats) == CK_OK);
      for (size_t c = 0; c < cells; ++c)
        NEAR(s->y[g * cells + c] + s->y[q * cells + c], 1, 1e-4);
    }
  }
  else if (IS("ProcessForcingJacobian") || IS("ConstraintResidualJacobian"))
  {
    int process = IS("ProcessForcingJacobian");
    double values[] = {process ? .8 : .6, process ? .2 : .4, process ? .3 : .2, process ? .7 : .8};
    for (size_t c = 0; c < 2; ++c)
    {
      s->y[g * 2 + c] = values[c * 2];
      s->y[q * 2 + c] = values[c * 2 + 1];
    }
    size_t n = 2, nnz = ck_jacobian_size(m);
    double *j = calloc(nnz * 2, sizeof(double)), plus[4] = {0}, minus[4] = {0};
    CHECK(j);
    if (process)
      ck_jacobian(s, s->y, j);
    else
      ck_constraint_jacobian(s, s->y, j);
    for (size_t col = 0; col < n; ++col)
      for (size_t c = 0; c < 2; ++c)
      {
        size_t index = col * 2 + c;
        double base = s->y[index], h = 1e-8 * fmax(1, fabs(base));
        s->y[index] = base + h;
        memset(plus, 0, sizeof plus);
        memset(minus, 0, sizeof minus);
        if (process)
          ck_forcing(s, s->y, plus);
        else
          ck_constraint_forcing(s, s->y, plus);
        s->y[index] = base - h;
        if (process)
          ck_forcing(s, s->y, minus);
        else
          ck_constraint_forcing(s, s->y, minus);
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
    free(j);
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
static int coupled_external(const char *test)
{
  int multi = IS("MultiEquilibriumKineticVsComposedConstraints"),
      built = IS("BuiltInVsExternalModelConstraintStepByStep"),
      rebuild = IS("BuilderCanBuildMultipleSolvers"),
      ordered = IS("ConstraintActiveWhenProcessOnlyModelAddedFirst"),
      temp = IS("TemperatureDependentConstraintParameter"), fd = IS("EquilibriumConstraintModelJacobian");
  int conservative = IS("KineticVsConservativeConstraintConvergence") || multi;
  size_t n = multi || rebuild ? 4 : fd ? 2 : 3;
  ck_species sp[] = {{.name = fd ? "A" : "A"}, {.name = fd ? "B" : "B"}, {.name = "C"}, {.name = "D"}};
  size_t reactants[] = {0, 1, 2, 1, 3};
  ck_product products[] = {{1, 1}, {2, 1}, {1, 1}, {3, 1}, {1, 1}};
  ck_reaction rx[5] = {0};
  double rate[] = {.1, 10, 2, 8, 8. / 3};
  for (size_t r = 0; r < 5; ++r)
  {
    rx[r] = (ck_reaction){.nreactants = 1,
                          .reactants = reactants + r,
                          .nproducts = 1,
                          .products = products + r,
                          .rate = ck_rate_defaults(CK_ARRHENIUS)};
    rx[r].rate.A = rate[r];
  }
  affine_model a = {.n = 2,
                    .names = {fd ? "A" : "B", fd ? "B" : "C"},
                    .nrows = 1,
                    .row_names = {fd ? "B" : "C"},
                    .coefficients = {{fd     ? 2.5
                                      : temp ? 2
                                             : 5,
                                      -1}},
                    .temperature = temp};
  affine_model b = {.n = 2, .names = {"B", "D"}, .nrows = 1, .row_names = {"D"}, .coefficients = {{3, -1}}};
  affine_model mass = {.n = n,
                       .names = {"A", "B", "C", "D"},
                       .nrows = 1,
                       .row_names = {"B"},
                       .coefficients = {{1, 1, 1, multi ? 1 : 0}},
                       .constant = {-1}};
  ck_named_jacobian_entry patterns[3][16];
  ck_external_model ext[3] = {external_descriptor(&a, patterns[0]), external_descriptor(&b, patterns[1]),
                              external_descriptor(&mass, patterns[2])};
  size_t ne = conservative ? (multi ? 3 : 2) : 1;
  if (conservative && !multi)
    ext[1] = ext[2];
  ck_mechanism d = {.name = "coupled",
                    .nspecies = n,
                    .species = sp,
                    .nreactions = fd ? 0 : 1,
                    .reactions = rx,
                    .nexternal = ne,
                    .external = ext,
                    .temperature = 298,
                    .pressure = 101325};
  ck_product cr[] = {{built ? 1 : 0, 1}}, cp[] = {{rebuild ? 3 : 2, 1}};
  ck_constraint constraint = {.type = CK_EQUILIBRIUM_CONSTRAINT,
                              .name = "builtin",
                              .algebraic_species = rebuild ? 3 : 2,
                              .nreactants = 1,
                              .reactants = cr,
                              .nproducts = 1,
                              .products = cp,
                              .equilibrium_ref = rebuild ? 2 : 5};
  ck_model *models[2] = {NULL, NULL};
  ck_state *states[2] = {NULL, NULL};
  size_t count =
      conservative || IS("SimpleConstraintPreservesRatioNotMass") || built || ordered || rebuild || temp ? 2
                                                                                                         : 1;
  for (size_t k = 0; k < count; ++k)
  {
    ck_mechanism desc = d;
    ck_external_model ordered_ext[2];
    if ((conservative || IS("SimpleConstraintPreservesRatioNotMass")) && k == 0)
    {
      desc.nexternal = 0;
      desc.nreactions = multi ? 5 : 3;
    }
    if (built && k == 0)
    {
      desc.nexternal = 0;
      desc.nconstraints = 1;
      desc.constraints = &constraint;
    }
    if (rebuild)
    {
      desc.nconstraints = 1;
      desc.constraints = &constraint;
    }
    if (ordered)
    {
      ordered_ext[0] = (ck_external_model){0};
      ordered_ext[1] = ext[0];
      if (k == 0)
      {
        ordered_ext[0] = ext[0];
        ordered_ext[1] = (ck_external_model){0};
      }
      desc.nexternal = 2;
      desc.external = ordered_ext;
    }
    CHECK(ck_model_create_with_order(&desc, 0, models + k) == CK_OK);
    CHECK(ck_state_create(models[k], 1, states + k) == CK_OK);
    states[k]->y[0] = 1;
    if (rebuild)
    {
      states[k]->y[3] = 2;
      CHECK(ck_model_nconstraints(models[k]) == 2);
    }
    if (temp)
      states[k]->conditions[0].temperature = k ? 350 : 298.15;
    if (ordered)
    {
      CHECK(ck_model_nconstraints(models[k]) == 1);
      CHECK(ck_is_algebraic(models[k], 2));
    }
  }
  if (fd)
  {
    ck_state *s = states[0];
    s->y[0] = 3;
    s->y[1] = 5;
    double f[2] = {0}, j[4] = {0};
    ck_constraint_forcing(s, s->y, f);
    ck_constraint_jacobian(s, s->y, j);
    NEAR(f[1], 2.5, 1e-12);
    for (size_t col = 0; col < 2; ++col)
    {
      double base = s->y[col], h = 1e-8 * fmax(1, base), plus[2] = {0}, minus[2] = {0};
      s->y[col] = base + h;
      ck_constraint_forcing(s, s->y, plus);
      s->y[col] = base - h;
      ck_constraint_forcing(s, s->y, minus);
      s->y[col] = base;
      for (size_t row = 0; row < 2; ++row)
      {
        double fdval = (plus[row] - minus[row]) / (2 * h);
        size_t off = ck_jacobian_index(models[0], row, col);
        if (off == SIZE_MAX)
          CHECK(fabs(fdval) <= 1e-6);
        else
          NEAR(-j[off], fdval, 1e-7 + 1e-7 * fmax(fabs(fdval), fabs(j[off])));
      }
    }
  }
  else
  {
    size_t steps = built ? 100 : ordered || rebuild ? 20 : 200;
    double dt = built ? .5 : 1;
    for (size_t t = 0; t < steps; ++t)
      for (size_t k = 0; k < count; ++k)
      {
        CHECK(ck_update_rates(states[k]) == CK_OK);
        ck_rosenbrock_parameters p;
        ck_stats stats;
        int method = (conservative || IS("SimpleConstraintPreservesRatioNotMass")) && k == 0 ? 3 : 43;
        CHECK(ck_rosenbrock_parameters_init(method, &p) == CK_OK);
        CHECK(ck_solve(states[k], dt, &p, &stats) == CK_OK);
        if (built && k == 1)
        {
          for (size_t j = 0; j < 3; ++j)
            NEAR(states[0]->y[j], states[1]->y[j], 1e-8);
          for (size_t l = 0; l < 2; ++l)
            NEAR(5 * states[l]->y[1] - states[l]->y[2], 0, 1e-6);
        }
      }
    if (conservative)
      for (size_t k = 0; k < 2; ++k)
      {
        double *y = states[k]->y;
        double B = 1. / (multi ? 9 : 6);
        NEAR(y[0], 0, 1e-6);
        NEAR(y[1], B, 1e-3);
        NEAR(y[2], 5 * B, 1e-3);
        if (multi)
          NEAR(y[3], 3 * B, 1e-3);
        double sum = 0;
        for (size_t j = 0; j < n; ++j)
          sum += y[j];
        NEAR(sum, 1, k ? 1e-3 : 1e-6);
        for (size_t j = 0; j < n; ++j)
          NEAR(states[0]->y[j], states[1]->y[j], j == 0 ? 1e-6 : 1e-3);
        if (multi && k == 1)
        {
          NEAR(5 * y[1] - y[2], 0, 1e-6);
          NEAR(3 * y[1] - y[3], 0, 1e-6);
        }
      }
    if (IS("SimpleConstraintPreservesRatioNotMass"))
    {
      for (size_t k = 0; k < 2; ++k)
        NEAR(states[k]->y[2] / states[k]->y[1], 5, 1e-3);
      NEAR(states[0]->y[0] + states[0]->y[1] + states[0]->y[2], 1, 1e-6);
      NEAR(states[1]->y[0] + states[1]->y[1] + states[1]->y[2], 6, .1);
    }
    if (ordered || rebuild)
    {
      for (size_t k = 0; k < 2; ++k)
      {
        CHECK(states[k]->y[1] > 0);
        NEAR(5 * states[k]->y[1] - states[k]->y[2], 0, 1e-6);
        if (rebuild)
          NEAR(2 * states[k]->y[0] - states[k]->y[3], 0, 1e-6);
      }
      for (size_t j = 0; j < n; ++j)
        CHECK(states[0]->y[j] == states[1]->y[j]);
    }
    if (temp)
      for (size_t k = 0; k < 2; ++k)
      {
        double K = 2 * exp(3000 * (1 / 298.15 - 1 / states[k]->conditions[0].temperature));
        CHECK(states[k]->y[1] > 0);
        if (k)
          CHECK(K > 2);
        NEAR(states[k]->y[2] / states[k]->y[1], K, 1e-4);
      }
  }
  for (size_t k = 0; k < count; ++k)
  {
    ck_state_free(states[k]);
    ck_model_free(models[k]);
  }
  return 0;
}
int test_external(const char *id, const char *test)
{
  if (strncmp(id, "external_model_constraints:", 27))
    return 2;
  if (IS("KineticVsConservativeConstraintConvergence") || IS("SimpleConstraintPreservesRatioNotMass") ||
      IS("BuiltInVsExternalModelConstraintStepByStep") ||
      IS("MultiEquilibriumKineticVsComposedConstraints") || IS("TemperatureDependentConstraintParameter") ||
      IS("ConstraintActiveWhenProcessOnlyModelAddedFirst") || IS("BuilderCanBuildMultipleSolvers") ||
      IS("EquilibriumConstraintModelJacobian"))
    return coupled_external(test);
  return external_mass(test);
}
