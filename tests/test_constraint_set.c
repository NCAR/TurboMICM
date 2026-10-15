/* ConstraintSet cases use the bound constraint collection owned by ck_model. */
#include "test.h"

static int collection_case(const char *test)
{
  ck_species species[] = {{.name = "A"}, {.name = "B"}, {.name = "AB"}, {.name = "C"}, {.name = "D"}};
  ck_product r[] = {{0, 1}, {1, 1}, {3, 1}}, p[] = {{2, 1}, {4, 1}};
  ck_constraint constraints[] = {{.type = CK_EQUILIBRIUM_CONSTRAINT,
                                  .name = "A_B_eq",
                                  .algebraic_species = 2,
                                  .nreactants = 2,
                                  .reactants = r,
                                  .nproducts = 1,
                                  .products = p,
                                  .equilibrium_ref = .033,
                                  .enthalpy = -24000},
                                 {.type = CK_EQUILIBRIUM_CONSTRAINT,
                                  .name = "C_D_eq",
                                  .algebraic_species = 4,
                                  .nreactants = 1,
                                  .reactants = r + 2,
                                  .nproducts = 1,
                                  .products = p + 1,
                                  .equilibrium_ref = .033,
                                  .enthalpy = -24000}};
  size_t n = 3, nc = 1, cells = 1;
  double values[15] = {.01, .02, .05}, expected[15] = {0};
  if (IS("MultipleConstraints"))
  {
    n = 5;
    nc = 2;
  }
  if (IS("EmptyConstraintSet"))
    nc = 0;
  if (IS("ReplaceStateRowsMapsToAlgebraicSpecies"))
  {
    r[0].species = 1;
    constraints[0].nreactants = 1;
  }
  if (IS("UnknownSpeciesThrows"))
    r[0].species = 99;
  if (IS("AddForcingTerms"))
  {
    cells = 2;
    double v[] = {.2, .4, .6, .1, .3, .7};
    memcpy(values, v, sizeof v);
    expected[2] = -.59736;
    expected[5] = -.69901;
  }
  if (IS("ThreeDStateOneConstraint"))
  {
    r[0].species = 0;
    p[0].species = 1;
    constraints[0].nreactants = 1;
    constraints[0].algebraic_species = 1;
    cells = 2;
    double v[] = {10, .2, .5, 10, .33, .3};
    memcpy(values, v, sizeof v);
    expected[1] = .13;
  }
  if (IS("FourDStateTwoConstraints"))
  {
    n = 4;
    nc = 2;
    cells = 3;
    r[0].species = 0;
    p[0].species = 1;
    constraints[0].nreactants = 1;
    constraints[0].algebraic_species = 1;
    ck_product rr[] = {{2, 1}, {3, 1}};
    memcpy(r + 1, rr, sizeof rr);
    p[1].species = 0;
    constraints[1].reactants = r + 1;
    constraints[1].nreactants = 2;
    constraints[1].algebraic_species = 0;
    double v[] = {.33, .01089, 10, 1, .33, .01089, 5, 1, .33, .02, 5, 1};
    memcpy(values, v, sizeof v);
    expected[4] = expected[8] = -.165;
    expected[9] = -.00911;
  }
  if (IS("CoupledConstraintsSharedSpecies"))
  {
    nc = 2;
    r[0].species = 0;
    p[0].species = 1;
    constraints[0].nreactants = 1;
    constraints[0].algebraic_species = 1;
    r[2].species = 0;
    p[1].species = 2;
    constraints[1].algebraic_species = 2;
    values[0] = .1;
    values[1] = values[2] = .0033;
  }
  if (IS("VectorizedMatricesRespectGridCellIndexing"))
  {
    cells = 3;
    double v[] = {.01, .02, .05, .03, .01, .2, .001, .002, .004};
    memcpy(values, v, sizeof v);
    for (size_t c = 0; c < cells; ++c)
      expected[c * n + 2] = .033 * values[c * n] * values[c * n + 1] - values[c * n + 2];
  }
  ck_mechanism d = {.name = "constraint collection",
                    .nspecies = n,
                    .species = species,
                    .nconstraints = nc,
                    .constraints = constraints,
                    .temperature = 298.15,
                    .pressure = 101325};
  ck_model *m = NULL;
  ck_state *s = NULL;
  if (IS("UnknownSpeciesThrows"))
  {
    CHECK(ck_model_create_with_order(&d, 0, &m) == CK_UNKNOWN);
    CHECK(!m);
    return 0;
  }
  CHECK(ck_model_create_with_order(&d, 0, &m) == CK_OK);
  CHECK(ck_model_nconstraints(m) == nc);
  CHECK(ck_state_create(m, cells, &s) == CK_OK);
  CHECK(ck_update_rates(s) == CK_OK);
  for (size_t i = 0; i < n; ++i)
    for (size_t c = 0; c < cells; ++c)
      s->y[i * cells + c] = values[c * n + i];
  size_t nnz = ck_jacobian_size(m);
  double *j = calloc(nnz * cells, sizeof(double)), f[16] = {0};
  CHECK(j);
  if (IS("EmptyConstraintSet"))
    for (size_t i = 0; i < 4; ++i)
      f[i] = 1;
  ck_constraint_forcing(s, s->y, f);
  ck_constraint_jacobian(s, s->y, j);
  if (IS("EmptyConstraintSet"))
  {
    for (size_t i = 0; i < 4; ++i)
      CHECK(f[i] == 1);
  }
  else
  {
    for (size_t k = 0; k < nc; ++k)
    {
      ck_constraint *q = constraints + k;
      size_t row = q->algebraic_species;
      CHECK(ck_is_algebraic(m, row));
      size_t dependency_count = 0;
      for (size_t col = 0; col < n; ++col)
      {
        int present = col == row;
        for (size_t a = 0; a < q->nreactants; ++a)
          present |= q->reactants[a].species == col;
        for (size_t a = 0; a < q->nproducts; ++a)
          present |= q->products[a].species == col;
        if (present)
        {
          ++dependency_count;
          CHECK(ck_jacobian_index(m, row, col) != SIZE_MAX);
        }
      }
      CHECK(dependency_count == q->nreactants + q->nproducts);
      for (size_t c = 0; c < cells; ++c)
      {
        double residual = .033;
        for (size_t a = 0; a < q->nreactants; ++a)
          residual *= s->y[q->reactants[a].species * cells + c];
        residual -= s->y[q->products[0].species * cells + c];
        NEAR(f[row * cells + c], residual, 1e-12);
        for (size_t a = 0; a < q->nreactants; ++a)
        {
          double derivative = -.033;
          for (size_t b = 0; b < q->nreactants; ++b)
            if (a != b)
              derivative *= s->y[q->reactants[b].species * cells + c];
          NEAR(j[ck_jacobian_index(m, row, q->reactants[a].species) * cells + c], derivative, 1e-12);
        }
        NEAR(j[ck_jacobian_index(m, row, q->products[0].species) * cells + c], 1, 1e-12);
      }
    }
    if (IS("AddForcingTerms") || IS("ThreeDStateOneConstraint") || IS("FourDStateTwoConstraints") ||
        IS("CoupledConstraintsSharedSpecies") || IS("VectorizedMatricesRespectGridCellIndexing"))
      for (size_t i = 0; i < n; ++i)
        for (size_t c = 0; c < cells; ++c)
          NEAR(f[i * cells + c], expected[c * n + i], 1e-12);
  }
  free(j);
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
int test_constraint_set(const char *suite, const char *test)
{
  if (strcmp(suite, "ConstraintSet"))
    return 2;
  if (IS("VectorMatrix1") || IS("VectorMatrix2") || IS("VectorMatrix3"))
  {
    CHECK(collection_case("Construction") == 0);
    CHECK(collection_case("NonZeroJacobianElements") == 0);
    CHECK(collection_case("AddForcingTerms") == 0);
    return collection_case("SubtractJacobianTerms");
  }
  return collection_case(test);
}
