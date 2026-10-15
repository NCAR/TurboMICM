#include "test.h"
#include "tutorials.h"
int main(int argc, char **argv)
{
  if (argc != 2)
    return 2;
  FILE *cases = fopen(argv[1], "r");
  if (!cases)
  {
    perror(argv[1]);
    return 2;
  }
  char line[512];
  size_t passed = 0, failed = 0, na = 0;
  while (fgets(line, sizeof line, cases))
  {
    if (line[0] == '#' || line[0] == '\n')
      continue;
    char id[400], kind[20];
    if (sscanf(line, "%399s %19s", id, kind) != 2)
      return 2;
    int result = 2;
    char *colon = strchr(id, ':');
    if (!colon)
      result = test_tutorial(id);
    if (colon)
    {
      char suite[200], test[200];
      if (sscanf(colon + 1, "%199[^.].%199s", suite, test) != 2)
        return 2;
      result = test_rate(suite, test);
      if (result == 2)
        result = test_linear(suite, test);
      if (result == 2)
        result = test_constraint(suite, test);
      if (result == 2)
        result = test_integrator(suite, test);
      if (result == 2)
        result = test_initialization(suite, test);
      if (result == 2)
        result = test_matrix(suite, test);
      if (result == 2)
        result = test_sparse(suite, test);
      if (result == 2)
        result = test_state(suite, test);
      if (result == 2)
        result = test_system(suite, test);
      if (result == 2)
        result = test_constraint_set(suite, test);
      if (result == 2)
        result = test_builder(suite, test);
      if (result == 2)
        result = test_verification(suite, test);
      if (result == 2)
        result = test_analytical(id, test);
      if (result == 2)
        result = test_process(suite, test);
      if (result == 2)
        result = test_process_set(suite, test);
      if (result == 2)
        result = test_grouped(suite, test);
      if (result == 2)
        result = test_aerosol(id, test);
      if (result == 2)
        result = test_external(id, test);
      if (result == 2)
        result = test_dae(suite, test);
      if (result == 2)
        result = test_chemistry(id, test);
      if (result == 2)
        result = test_chapman(id, suite, test);
      if (result == 2)
        result = test_stiff(id, suite, test);
    }
    if (result == 0)
    {
      printf("%s PASS assertions verified\n", id);
      ++passed;
    }
    else if (result == 3)
    {
      printf("%s N/A C++ compile-time concepts/overload categories; unified C view has no equivalent\n", id);
      ++na;
    }
    else
    {
      printf("%s FAIL %s\n", id, result == 2 ? "port not implemented" : "assertion failed");
      ++failed;
    }
  }
  fclose(cases);
  fprintf(stderr, "cases=%zu PASS=%zu FAIL=%zu N/A=%zu\n", passed + failed + na, passed, failed, na);
  return failed ? 1 : 0;
}
