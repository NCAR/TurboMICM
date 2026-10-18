/* Runs the ported MICM cases in the order of the reference case list and
 * prints one line per case. A case in the list with no port prints FAIL.
 *
 *   test_cpu [FILTER]      run only the cases whose id contains FILTER
 *   test_cpu --list        print the registered ids
 *   test_cpu --extra       run the registered checks that are not MICM cases
 */
#include "harness.h"
#include "registry.h"

#include <stdarg.h>

static char reason[256];
static int failures;
static const char *current_id;

void tmc_test_set_reason(const char *r)
{
  snprintf(reason, sizeof reason, "%s", r);
}

void tmc_test_fail(const char *file, int line, const char *format, ...)
{
  ++failures;
  if (failures <= 3)
  {
    fprintf(stderr, "  %s: %s:%d: ", current_id, file, line);
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fputc('\n', stderr);
  }
}

bool tmc_test_double_eq(double a, double b)
{
  if (a == b)
    return true;
  if (isnan(a) || isnan(b))
    return false;
  /* 4 ULPs, as GoogleTest's AlmostEquals */
  double lo = a < b ? a : b, hi = a < b ? b : a;
  for (int i = 0; i < 4; ++i)
  {
    lo = nextafter(lo, hi);
    if (lo >= hi)
      return true;
  }
  return false;
}

static char path_buffer[4][1024];
static int path_slot;

static const char *join(const char *base, const char *relative)
{
  char *p = path_buffer[path_slot++ % 4];
  snprintf(p, sizeof path_buffer[0], "%s/%s", base, relative);
  return p;
}

const char *tmc_test_micm_path(const char *relative)
{
  const char *base = getenv("TMC_TEST_MICM");
  return join(base ? base : "../micm", relative);
}

const char *tmc_test_reference_path(const char *relative)
{
  const char *base = getenv("TMC_TEST_REFERENCE");
  return join(base ? base : "../reference", relative);
}

const char *tmc_test_data_path(const char *relative)
{
  return join("tests/data", relative);
}

/* The case ids of the reference list, in order (generated into case_list.inc). */
static const char *const case_list[] = {
#include "case_list.inc"
};

int main(int argc, char **argv)
{
  const char *filter = argc > 1 && argv[1][0] != '-' ? argv[1] : NULL;
  bool list_only = argc > 1 && strcmp(argv[1], "--list") == 0;
  bool extra_only = argc > 1 && strcmp(argv[1], "--extra") == 0;
  int n_registered = 0;
  const tmc_test_case *const *groups = tmc_test_groups(&n_registered);
  if (list_only)
  {
    for (int g = 0; groups[g]; ++g)
      for (const tmc_test_case *c = groups[g]; c->id; ++c)
        printf("%s\n", c->id);
    return 0;
  }
  const size_t n_cases = sizeof case_list / sizeof case_list[0];
  int pass = 0, fail = 0, na = 0;
  if (extra_only)
  {
    /* checks of this port's own (generated kernels, the TS1 gate) */
    for (int g = 0; groups[g]; ++g)
      for (const tmc_test_case *c = groups[g]; c->id; ++c)
      {
        bool listed = false;
        for (size_t i = 0; i < n_cases && !listed; ++i)
          listed = strcmp(case_list[i], c->id) == 0;
        if (listed)
          continue;
        current_id = c->id;
        failures = 0;
        reason[0] = '\0';
        int status = c->run();
        if (status == TMC_TEST_PASS && failures > 0)
          status = TMC_TEST_FAIL;
        if (status == TMC_TEST_PASS)
          ++pass;
        else if (status == TMC_TEST_NA)
          ++na;
        else
          ++fail;
        printf("%s %s%s%s\n", c->id, status == TMC_TEST_PASS ? "PASS" : status == TMC_TEST_NA ? "N/A" : "FAIL", reason[0] ? " " : "", reason);
      }
    fprintf(stderr, "extra checks: %d PASS, %d FAIL, %d N/A\n", pass, fail, na);
    return fail ? 1 : 0;
  }
  for (size_t i = 0; i < n_cases; ++i)
  {
    const char *id = case_list[i];
    if (filter && !strstr(id, filter))
      continue;
    const tmc_test_case *found = NULL;
    for (int g = 0; groups[g] && !found; ++g)
      for (const tmc_test_case *c = groups[g]; c->id; ++c)
        if (strcmp(c->id, id) == 0)
        {
          found = c;
          break;
        }
    current_id = id;
    failures = 0;
    reason[0] = '\0';
    int status;
    if (!found)
    {
      status = TMC_TEST_FAIL;
      tmc_test_set_reason("not ported");
    }
    else
    {
      status = found->run();
      if (status == TMC_TEST_PASS && failures > 0)
        status = TMC_TEST_FAIL;
    }
    if (status == TMC_TEST_PASS)
    {
      ++pass;
      printf("%s PASS\n", id);
    }
    else if (status == TMC_TEST_NA)
    {
      ++na;
      printf("%s N/A %s\n", id, reason);
    }
    else
    {
      ++fail;
      printf("%s FAIL%s%s\n", id, reason[0] ? " " : "", reason);
    }
    fflush(stdout);
  }
  fprintf(stderr, "%d PASS, %d FAIL, %d N/A\n", pass, fail, na);
  return fail ? 1 : 0;
}
