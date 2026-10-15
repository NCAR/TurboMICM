/* Test values and tolerances ported from MICM's rate-constant unit tests. */
#include "test.h"
static double twice_temperature(const ck_conditions *c, const double *p, void *context)
{
  (void)p;
  (void)context;
  return 2 * c->temperature;
}
int test_rate(const char *suite, const char *test)
{
  ck_conditions c = {.temperature = 301.24, .air_density = 42.2};
  ck_rate r;
  double k, expected, t = c.temperature;
  if (!strcmp(suite, "ArrheniusRateConstant") || !strcmp(suite, "TaylorSeriesRateConstant"))
  {
    int taylor = !strcmp(suite, "TaylorSeriesRateConstant");
    if (!IS("CalculateWithSystem") && !IS("CalculateWithPrescribedArguments"))
      return 2;
    r = ck_rate_defaults(taylor ? CK_TAYLOR : CK_ARRHENIUS);
    NEAR(ck_rate_calculate(&r, &c, NULL), 1, 1e-13);
    r.A = 1;
    NEAR(ck_rate_calculate(&r, &c, NULL), 1, 1e-13);
    r.A = 2.2e-10;
    NEAR(ck_rate_calculate(&r, &c, NULL), 2.2e-10, 1e-13 * 2.2e-10);
    r.A = 3e-11;
    r.C = -200;
    if (taylor)
    {
      r.ncoefficients = 3;
      r.coefficients[0] = 12.5;
      r.coefficients[1] = 1.3e-2;
      r.coefficients[2] = 5.2e-4;
    }
    expected = 3e-11 * exp(-200 / 301.24) * (taylor ? 12.5 + 1.3e-2 * 301.24 + 5.2e-4 * pow(301.24, 2) : 1);
    NEAR(ck_rate_calculate(&r, &c, NULL), expected, 1e-13 * expected);
    r.A = 2.6e-12;
    r.C = -350;
    if (taylor)
    {
      r.coefficients[0] = 1;
      r.coefficients[1] = 4.3e-1;
      r.coefficients[2] = 7.3e-3;
    }
    expected = 2.6e-12 * exp(-350 / 301.24) * (taylor ? 1 + 4.3e-1 * 301.24 + 7.3e-3 * pow(301.24, 2) : 1);
    NEAR(ck_rate_calculate(&r, &c, NULL), expected, 1e-13 * expected);
    return 0;
  }
  if (!strcmp(suite, "TroeRateConstant") || !strcmp(suite, "TernaryChemicalActivationRateConstant"))
  {
    int troe = !strcmp(suite, "TroeRateConstant");
    r = ck_rate_defaults(troe ? CK_TROE : CK_TERNARY);
    double k0 = 1, kinf = 1;
    if (IS("CalculateWithMinimalArguments"))
    {
      k = ck_rate_calculate(&r, &c, NULL);
      expected = (troe ? 42.2 : 1) * k0 / (1 + 42.2 * k0 / kinf) *
                 pow(0.6, 1.0 / (1 + pow(log10(42.2 * k0 / kinf), 2)));
      NEAR(k, expected, 0);
    }
    else if (IS("CalculateWithAllArguments"))
    {
      r.k0_A = 1.2;
      r.k0_B = 2.3;
      r.k0_C = 302.3;
      r.kinf_A = 2.6;
      r.kinf_B = -3.1;
      r.kinf_C = 402.1;
      r.Fc = 0.9;
      r.N = 1.2;
      k0 = 1.2 * exp(302.3 / t) * pow(t / 300, 2.3);
      kinf = 2.6 * exp(402.1 / t) * pow(t / 300, -3.1);
      expected = (troe ? 42.2 : 1) * k0 / (1 + 42.2 * k0 / kinf) *
                 pow(0.9, 1.0 / (1.0 + 1.0 / 1.2 * pow(log10(42.2 * k0 / kinf), 2)));
      NEAR(ck_rate_calculate(&r, &c, NULL), expected, 0);
    }
    else if (troe && IS("AnalyticalTroeExampleAB"))
    {
      r.k0_A = 4e-18;
      k0 = 4e-18;
      expected = k0 * 42.2 / (1 + k0 * 42.2) * pow(0.6, 1.0 / (1 + pow(log10(k0 * 42.2), 2)));
      NEAR(ck_rate_calculate(&r, &c, NULL), expected, 0);
    }
    else if (troe && IS("AnalyticalTroeExampleBC"))
    {
      r.k0_A = 1.2e-12;
      r.k0_B = 167;
      r.k0_C = 3;
      r.kinf_A = 136;
      r.kinf_B = 5;
      r.kinf_C = 24;
      r.Fc = 0.9;
      r.N = 0.8;
      k0 = 1.2e-12 * exp(3.0 / 301.24) * pow(301.24 / 300, 167);
      kinf = 136 * exp(24.0 / 301.24) * pow(301.24 / 300, 5);
      expected = k0 * 42.2 / (1 + k0 * 42.2 / kinf) *
                 pow(0.9, 1.0 / (1 + (1.0 / 0.8) * pow(log10(k0 * 42.2 / kinf), 2)));
      k = ck_rate_calculate(&r, &c, NULL);
      CHECK(fabs(k - expected) / fmax(fabs(k), fabs(expected)) <= 1e-14);
    }
    else
      return 2;
    return 0;
  }
  if (!strcmp(suite, "TernaryChemicalActivationJPL19RateConstant"))
  {
    r = ck_rate_defaults(CK_TERNARY_JPL19);
    double k0, kinf, ratio;
    if (IS("CalculateWithAllArguments"))
    {
      r.k0_A = 1.2;
      r.k0_B = 2.3;
      r.k0_C = 302.3;
      r.k0_D = 298;
      r.kinf_A = 2.6;
      r.kinf_B = -3.1;
      r.kinf_C = 402.1;
      r.kinf_D = 310;
      r.kint_A = 0.7;
      r.kint_B = 0.5;
      r.kint_C = 120;
      r.kint_D = 290;
      r.Fc = 0.9;
      r.N = 1.2;
      k0 = 1.2 * exp(302.3 / t) * pow(t / 298, 2.3);
      kinf = 2.6 * exp(402.1 / t) * pow(t / 310, -3.1);
      ratio = k0 * c.air_density / kinf;
      double kint = 0.7 * exp(120 / t) * pow(t / 290, 0.5);
      double kf = k0 * c.air_density / (1 + ratio) * pow(0.9, 1.2 / (1.2 + pow(log10(ratio), 2)));
      expected = kf + kint * (1 - kf / kinf);
      NEAR(ck_rate_calculate(&r, &c, NULL), expected, 0);
    }
    else if (IS("ZeroKintReducesToTroe"))
    {
      c.temperature = 250;
      c.air_density = 30;
      r.k0_A = 3e-3;
      r.k0_B = -2;
      r.kinf_A = 1.5;
      r.kinf_B = -1;
      r.kint_A = 0;
      k0 = 3e-3 * pow(250.0 / 298, -2);
      kinf = 1.5 * pow(250.0 / 298, -1);
      ratio = k0 * 30 / kinf;
      expected = k0 * 30 / (1 + ratio) * pow(0.6, 1.0 / (1 + pow(log10(ratio), 2)));
      NEAR(ck_rate_calculate(&r, &c, NULL), expected, 1e-30);
    }
    else
      return 2;
    return 0;
  }
  if (!strcmp(suite, "TunnelingRateConstant"))
  {
    r = ck_rate_defaults(CK_TUNNELING);
    if (IS("CalculateWithMinimalArguments"))
      expected = 1;
    else if (IS("CalculateWithAllArguments"))
    {
      r.A = 1.2;
      r.B = 2.3;
      r.C = 302.3;
      expected = 1.2 * exp(-2.3 / t) * exp(302.3 / pow(t, 3));
    }
    else
      return 2;
    NEAR(ck_rate_calculate(&r, &c, NULL), expected, 1e-13 * expected);
    return 0;
  }
  if (!strcmp(suite, "BranchedRateConstant"))
  {
    int alkoxy = IS("CalculateAlkoxyBranchWithAllArguments");
    if (!alkoxy && !IS("CalculateNitrateBranchWithAllArguments"))
      return 2;
    r = ck_rate_defaults(alkoxy ? CK_BRANCH_ALKOXY : CK_BRANCH_NITRATE);
    r.X = 1.2;
    r.Y = 204.3;
    r.a0 = 1e-3;
    r.n = 2;
    CHECK(ck_rate_prepare(&r) == CK_OK);
    double air = 42.2 * CK_AVOGADRO * 1e-6;
    double a = 2e-22 * exp(2) * 2.45e19, b = 0.43 * pow(293.0 / 298, -8);
    double z = a / (1 + a / b) * pow(0.41, 1.0 / (1 + pow(log10(a / b), 2))) * (1 - 1e-3) / 1e-3;
    a = 2e-22 * exp(2) * air;
    b = 0.43 * pow(t / 298, -8);
    double val = a / (1 + a / b) * pow(0.41, 1.0 / (1 + pow(log10(a / b), 2)));
    expected = 1.2 * exp(-204.3 / t) * (alkoxy ? z / (z + val) : val / (z + val));
    NEAR(ck_rate_calculate(&r, &c, NULL), expected, 1e-13 * expected);
    return 0;
  }
  if (!strcmp(suite, "ReversibleRateConstant"))
  {
    r = ck_rate_defaults(CK_REVERSIBLE);
    if (IS("DefaultConstructor"))
    {
      c.temperature = 298.15;
      NEAR(ck_rate_calculate(&r, &c, NULL), 0, 1e-13);
    }
    else if (IS("CalculateWithSystem") || IS("CalculateWithPrescribedArguments"))
    {
      if (IS("CalculateWithPrescribedArguments"))
        c.temperature = 298.15;
      r.A = 1.14e-2;
      r.C = 2300;
      r.reverse = 0.32;
      expected = 1.14e-2 * exp(2300 / c.temperature) * 0.32;
      NEAR(ck_rate_calculate(&r, &c, NULL), expected, 1e-13 * expected);
    }
    else if (IS("TemperatureDependence"))
    {
      double temperatures[] = {273.15, 298.15, 323.15};
      r.A = 1.5;
      r.C = 2000;
      r.reverse = 0.5;
      for (size_t i = 0; i < 3; ++i)
      {
        c.temperature = temperatures[i];
        expected = 1.5 * exp(2000 / c.temperature) * 0.5;
        NEAR(ck_rate_calculate(&r, &c, NULL), expected, 1e-13 * expected);
      }
    }
    else if (IS("ZeroReverseRate"))
    {
      r.A = 5;
      r.C = 1000;
      r.reverse = 0;
      c.temperature = 298.15;
      NEAR(ck_rate_calculate(&r, &c, NULL), 0, 1e-13);
    }
    else if (IS("EquilibriumConstantOnly"))
    {
      r.A = 3.5;
      r.C = 1800;
      r.reverse = 1;
      c.temperature = 298.15;
      expected = 3.5 * exp(1800 / 298.15);
      NEAR(ck_rate_calculate(&r, &c, NULL), expected, 1e-13 * expected);
    }
    else
      return 2;
    return 0;
  }
  if (!strcmp(suite, "UserDefinedRateConstant"))
  {
    r = ck_rate_defaults(CK_USER);
    double custom = 1.1;
    expected = 1.1;
    if (IS("CalculateWithSystem"))
      custom = expected = 0.5;
    else if (IS("CustomScalingFactor"))
    {
      custom = 1.2;
      r.scaling = 2;
      expected = 1.2 * 2;
    }
    else if (!IS("ConstructorWithRate") && !IS("ConstructorWithRateAndName"))
      return 2;
    NEAR(ck_rate_calculate(&r, &c, &custom), expected, 0);
    return 0;
  }
  if (!strcmp(suite, "SurfaceRateConstant"))
  {
    r = ck_rate_defaults(CK_SURFACE);
    r.diffusion = 2.3e2;
    r.molecular_weight = 0.025;
    if (IS("DiffusionCoefficientIsMissing"))
    {
      r.diffusion = 0;
      CHECK(ck_rate_prepare(&r) == CK_INVALID);
      return 0;
    }
    if (IS("MolecularWeightIsMissing"))
    {
      r.molecular_weight = 0;
      CHECK(ck_rate_prepare(&r) == CK_INVALID);
      return 0;
    }
    if (IS("CalculateSpecifiedProbability"))
      r.probability = 0.74;
    else if (!IS("CalculateDefaultProbability"))
      return 2;
    c.temperature = 273.65;
    double custom[] = {1e-7, 2.5e6};
    expected = 4 * 2.5e6 * CK_PI * pow(1e-7, 2) /
               (1e-7 / 2.3e2 + 4 / (sqrt(8 * CK_GAS_CONSTANT * 273.65 / (CK_PI * 0.025)) * r.probability));
    NEAR(ck_rate_calculate(&r, &c, custom), expected, 1e-13 * expected);
    return 0;
  }
  if (!strcmp(suite, "LambdaRateConstant") && IS("CalculateWithSystem"))
  {
    r = ck_rate_defaults(CK_LAMBDA);
    r.function = twice_temperature;
    c.temperature = 300;
    NEAR(ck_rate_calculate(&r, &c, NULL), 600, 1e-13 * 600);
    return 0;
  }
  return 2;
}
