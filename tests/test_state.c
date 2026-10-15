#include "test.h"
static const ck_species species[] = {{.name = "foo"}, {.name = "bar"}, {.name = "baz"}, {.name = "quz"}};
int test_state(const char *suite, const char *test)
{
  if (strcmp(suite, "State"))
    return 2;
  if (IS("DefaultConstructor"))
  {
    ck_state state = {0}, states[10] = {{0}};
    CHECK(state.y == NULL && state.cells == 0);
    for (size_t i = 0; i < 10; ++i)
      CHECK(states[i].model == NULL && states[i].custom == NULL);
    return 0;
  }
  const char *labels[] = {"O1", "O2", "O3", "AAA", "BBB"};
  int copy = IS("CopyConstructor") || IS("CopyAssignmentOperator") || IS("MoveConstructor") ||
             IS("MoveAssignmentOperator"),
      constructor = IS("Constructor");
  if (copy || constructor)
  {
    labels[0] = "quux";
    labels[1] = "corge";
  }
  ck_reaction reactions[10];
  for (size_t j = 0; j < 10; ++j)
  {
    reactions[j] =
        (ck_reaction){.rate = ck_rate_defaults(j < (copy || constructor ? 2u : 5u) ? CK_USER : CK_ARRHENIUS)};
    if (j < (copy || constructor ? 2u : 5u))
      reactions[j].rate.label = labels[j];
  }
  ck_mechanism mechanism = {
      .name = "state", .nspecies = 4, .species = species, .nreactions = 10, .reactions = reactions};
  ck_model *m;
  CHECK(ck_model_create_with_order(&mechanism, 0, &m) == CK_OK);
  int single = IS("SetConcentrationByElementSingleValue") || IS("SetStateWithSquareBracketOperator") ||
               IS("SettingCustomRateParameterWithInvalidLabelThrowsException") ||
               IS("SetCustomRateParameter") || IS("UnsafelySetCustomRateParameterOneCell");
  size_t cells = single ? 1
                 : IS("UnsafelySetCustomRateParameterCatchesTooFewGridCells") ||
                         IS("UnsafelySetCustomRateParameterCatchesTooParameters")
                     ? 2
                     : 3;
  ck_state *s;
  CHECK(ck_state_create(m, cells, &s) == CK_OK);
  if (constructor || copy)
  {
    CHECK(s->cells == 3 && ck_model_nspecies(m) == 4 && ck_model_nreactions(m) == 10 &&
          ck_model_ncustom(m) == 2);
    for (size_t j = 0; j < 4; ++j)
    {
      size_t index;
      CHECK(ck_species_index(m, species[j].name, &index) == CK_OK && index == j);
    }
    for (size_t j = 0; j < 2; ++j)
    {
      size_t index;
      CHECK(ck_custom_index(m, labels[j], &index) == CK_OK && index == j);
    }
    if (copy)
    {
      double atol[] = {1e-10, 1e-10, 1e-10, 1e-10};
      CHECK(ck_state_tolerances(s, 1e-5, atol, 4) == CK_OK);
      ck_state *other = NULL;
      int move = strstr(test, "Move") != NULL;
      if (move)
      {
        ck_state *old = s;
        CHECK(ck_state_move(&s, &other) == CK_OK);
        CHECK(s == NULL);
        s = other;
        other = NULL;
        CHECK(s == old);
        CHECK(ck_state_clone(s, &other) == CK_OK);
      }
      else
        CHECK(ck_state_clone(s, &other) == CK_OK);
      CHECK(other->cells == s->cells && other->model == s->model &&
            other->relative_tolerance == s->relative_tolerance);
      for (size_t j = 0; j < 4; ++j)
        NEAR(other->absolute_tolerance[j], s->absolute_tolerance[j], 0);
      for (size_t i = 0; i < 4 * cells; ++i)
        NEAR(other->y[i], s->y[i], 0);
      for (size_t i = 0; i < 2 * cells; ++i)
        NEAR(other->custom[i], s->custom[i], 0);
      for (size_t i = 0; i < 10 * cells; ++i)
        NEAR(other->rates[i], s->rates[i], 0);
      other->y[0] = 42;
      NEAR(s->y[0], 0, 0);
      ck_state_free(other);
    }
  }
  else if (IS("SettingSingleConcentrationWithInvalidArgumentsThowsException"))
  {
    double values[] = {1, 2, 3, 4};
    CHECK(ck_state_set(s, "foo", values, 1) == CK_INVALID);
    CHECK(ck_state_set(s, "foo", values, 2) == CK_INVALID);
    CHECK(ck_state_set(s, "not foo", values, 1) != CK_OK);
    CHECK(ck_state_set(s, "not foo", values, 4) != CK_OK);
  }
  else if (IS("SetSingleConcentration") || IS("SetConcentrationByElementVector"))
  {
    double values[] = {12, 42, 35.2};
    const char *name = IS("SetSingleConcentration") ? "bar" : "foo";
    size_t index;
    CHECK(ck_species_index(m, name, &index) == CK_OK);
    CHECK(ck_state_set(s, name, values, 3) == CK_OK);
    for (size_t c = 0; c < 3; ++c)
      NEAR(s->y[index * 3 + c], values[c], 0);
    if (IS("SetSingleConcentration"))
    {
      ck_state *one;
      CHECK(ck_state_create(m, 1, &one) == CK_OK);
      double value = 324.2;
      CHECK(ck_state_set(one, "bar", &value, 1) == CK_OK);
      NEAR(one->y[1], value, 0);
      ck_state_free(one);
    }
  }
  else if (IS("SetConcentrationByElementSingleValue"))
  {
    double x = 42, y = 12;
    CHECK(ck_state_set(s, "baz", &x, 1) == CK_OK);
    CHECK(ck_state_set(s, "quz", &y, 1) == CK_OK);
    NEAR(s->y[2], 42, 0);
    NEAR(s->y[3], 12, 0);
  }
  else if (IS("SetConcentrations") || IS("SettingConcentrationsWithInvalidArguementsThrowsException"))
  {
    double values[][3] = {{0.9, 0.99, 0.999}, {0.2, 0.22, 0.222}, {0.3, 0.33, 0.333}, {0.4, 0.44, 0.444}};
    ck_named_values named[4];
    for (size_t j = 0; j < 4; ++j)
      named[j] = (ck_named_values){species[j].name, values[j], 3};
    if (IS("SettingConcentrationsWithInvalidArguementsThrowsException"))
    {
      named[0].name = "FUU";
      CHECK(ck_state_set_many(s, named, 4, 0) != CK_OK);
    }
    else
    {
      CHECK(ck_state_set_many(s, named, 4, 0) == CK_OK);
      for (size_t j = 0; j < 4; ++j)
        for (size_t c = 0; c < 3; ++c)
          NEAR(s->y[j * 3 + c], values[j][c], 0);
    }
  }
  else if (IS("SetStateWithSquareBracketOperator"))
  {
    double values[] = {42, 12, 35.2, 24.2};
    for (size_t j = 0; j < 4; ++j)
    {
      CHECK(ck_state_set(s, species[j].name, &values[j], 1) == CK_OK);
      NEAR(s->y[j], values[j], 0);
    }
    s->y[0] += 8;
    NEAR(s->y[0], 50, 0);
    s->y[0] -= 20;
    NEAR(s->y[0], 30, 0);
    s->y[0] *= 2;
    NEAR(s->y[0], 60, 0);
    s->y[0] /= 3;
    NEAR(s->y[0], 20, 0);
  }
  else if (IS("SetMultiCellStateWithSquareBrackets"))
  {
    double foo[] = {1, 2, 3}, bar[] = {4, 5, 6};
    CHECK(ck_state_set(s, "foo", foo, 3) == CK_OK);
    CHECK(ck_state_set(s, "bar", bar, 3) == CK_OK);
    for (size_t c = 0; c < 3; ++c)
    {
      NEAR(s->y[c], foo[c], 0);
      NEAR(s->y[3 + c], bar[c], 0);
    }
  }
  else if (IS("SettingCustomRateParameterWithInvalidVectorSizeThrowsException"))
  {
    double value = 0.3;
    CHECK(ck_state_set_custom(s, "O3", &value, 1) == CK_INVALID);
  }
  else if (IS("SettingCustomRateParameterWithInvalidLabelThrowsException"))
  {
    double value = 0.7;
    CHECK(ck_state_set_custom(s, "CCC", &value, 1) == CK_UNKNOWN);
  }
  else if (IS("SetCustomRateParameter"))
  {
    double value = 42.3;
    CHECK(ck_state_set_custom(s, "O2", &value, 1) == CK_OK);
    NEAR(s->custom[1], 42.3, 0);
  }
  else if (IS("SetCustomRateParameters"))
  {
    double values[][3] = {
        {0.1, 0.11, 0.111}, {0.5, 0.55, 0.555}, {0.3, 0.33, 0.333}, {0.5, 0.55, 0.555}, {0.7, 0.77, 0.777}};
    ck_named_values named[5];
    for (size_t j = 0; j < 5; ++j)
      named[j] = (ck_named_values){labels[j], values[j], 3};
    CHECK(ck_state_set_many(s, named, 5, 1) == CK_OK);
    for (size_t j = 0; j < 5; ++j)
      for (size_t c = 0; c < 3; ++c)
        NEAR(s->custom[j * 3 + c], values[j][c], 0);
  }
  else if (!strncmp(test, "UnsafelySetCustomRateParameter", 29))
  {
    double rows[15];
    for (size_t c = 0; c < 3; ++c)
      for (size_t j = 0; j < 5; ++j)
        rows[c * 5 + j] = (j + 1) * 0.1;
    if (IS("UnsafelySetCustomRateParameterCatchesTooFewGridCells"))
      CHECK(ck_state_set_custom_rows(s, rows, 1, 5) == CK_INVALID);
    else if (IS("UnsafelySetCustomRateParameterCatchesTooParameters"))
      CHECK(ck_state_set_custom_rows(s, rows, 1, 3) == CK_INVALID);
    else
    {
      CHECK(ck_state_set_custom_rows(s, rows, cells, 5) == CK_OK);
      for (size_t c = 0; c < cells; ++c)
        for (size_t j = 0; j < 5; ++j)
          NEAR(s->custom[j * cells + c], rows[c * 5 + j], 0);
    }
  }
  else
  {
    ck_state_free(s);
    ck_model_free(m);
    return 2;
  }
  ck_state_free(s);
  ck_model_free(m);
  return 0;
}
