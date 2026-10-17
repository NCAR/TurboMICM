/* MICM test/unit/version/test_version.cpp */
#include "harness.h"

TMC_TEST(full_version)
{
  const char *v = tmc_version();
  TMC_EXPECT_TRUE(v != NULL && strlen(v) > 0);
  return TMC_TEST_PASS;
}

TMC_TEST(major)
{
  unsigned v = tmc_version_major();
  TMC_EXPECT_TRUE(v == 3);
  return TMC_TEST_PASS;
}

TMC_TEST(minor)
{
  unsigned v = tmc_version_minor();
  TMC_EXPECT_TRUE(v == 7);
  return TMC_TEST_PASS;
}

TMC_TEST(patch)
{
  unsigned v = tmc_version_patch();
  TMC_EXPECT_TRUE(v < 1000u);
  return TMC_TEST_PASS;
}

TMC_TEST(tweak)
{
  /* MICM's tweak number counts commits past the tag; the port has no build-time tweak. */
  TMC_EXPECT_TRUE(tmc_version_patch() == 0);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("version:Version.FullVersion", full_version),
  TMC_CASE("version:Version.VersionMajor", major),
  TMC_CASE("version:Version.VersionMinor", minor),
  TMC_CASE("version:Version.VersionPatch", patch),
  TMC_CASE("version:Version.VersionTweak", tweak),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_version(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
