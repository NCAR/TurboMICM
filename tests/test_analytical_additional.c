#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define COUNT(items) (sizeof(items) / sizeof((items)[0]))
#define PI 3.14159265358979323846264338327950288
#define GAS_CONSTANT 8.31446261815324

static int integrate(
    const char *method, const chem_mechanism *mechanism, double *state,
    double interval, chem_rosenbrock_parameters *rosenbrock) {
  const double temperature[] = { 272.5 };
  const double pressure[] = { 101253.3 };
  if (strcmp(method, "backward_euler") == 0)
    return chem_integrate_group_backward_euler(mechanism, 1u, state, NULL, temperature, pressure, interval, NULL);
  return chem_integrate_group_with_parameters(mechanism, 1u, state, NULL, temperature, pressure,
                                               interval, rosenbrock, NULL);
}

static int near(double actual, double expected, double tolerance) {
  return isfinite(actual) && fabs(actual - expected) <= tolerance;
}

static const chem_species surface_species[] = {
  { "foo", 0.04607, 0.95e-5, 0.0, 1.0e-14, 0u },
  { "bar", 0.0, 0.0, 0.0, 1.0e-14, 0u },
  { "baz", 0.0, 0.0, 0.0, 1.0e-14, 0u },
};
static const chem_parameter surface_parameters[] = {
  { "foo.effective radius [m]", 0.0 },
  { "foo.particle number concentration [# m-3]", 0.0 },
};
static const chem_term surface_reactants[] = { { 0u, 1.0 } };
static const chem_term surface_products[] = { { 1u, 1.0 }, { 2u, 0.4 } };
static const chem_reaction surface_reactions[] = {
  { 0u, 1u, 0u, 2u, { .kind = CHEM_RATE_SURFACE, .custom_parameter = 0u,
      .surface_species = 0u, .p = { [0] = 2.0e-2 } } },
};
static chem_mechanism make_surface_mechanism(void) {
  chem_mechanism mechanism = {
    .name = "analytical surface reaction", .species_count = COUNT(surface_species),
    .reaction_count = COUNT(surface_reactions), .reactant_term_count = COUNT(surface_reactants),
    .product_term_count = COUNT(surface_products), .parameter_count = COUNT(surface_parameters),
    .species = surface_species, .reactions = surface_reactions, .reactant_terms = surface_reactants,
    .product_terms = surface_products, .parameters = surface_parameters,
    .temperature = 272.5, .pressure = 101253.3, .relative_tolerance = 1.0e-6,
  };
  return mechanism;
}

static int surface_reaction(const char *method) {
  const double gmd = 1.0e-6, gsd = 0.1, density = 1000.0;
  const double conc_stuff = 2.0e-3, conc_more_stuff = 3.0e-3;
  const double radius = gmd / 2.0 * exp(5.0 * log(gsd) * log(gsd) / 2.0);
  const double number = 6.0 / (PI * pow(gmd, 3.0) * exp(9.0 / 2.0 * log(gsd) * log(gsd))) *
      (conc_stuff / density + conc_more_stuff / density);
  chem_mechanism mechanism = make_surface_mechanism();
  chem_parameter parameters[] = {
    { "foo.effective radius [m]", radius },
    { "foo.particle number concentration [# m-3]", number },
  };
  const double parameter_values[] = { radius, number };
  mechanism.parameters = parameters;
  double rate = 0.0;
  if (chem_evaluate_rate(&mechanism, 0u, parameter_values, 272.5, 101253.3, &rate) != 0 || rate <= 0.0) return 1;
  double state[] = { 1.0, 0.0, 0.0 };
  chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  const double dt = 0.1 / rate;
  const double tolerance = strcmp(method, "backward_euler") == 0 ? 5.0e-2 : 1.0e-2;
  for (size_t step = 1u; step <= 10u; ++step) {
    if (integrate(method, &mechanism, state, dt, &solver) != 0) return 1;
    const double expected_foo = exp(-rate * (double)step * dt);
    const double expected_bar = 1.0 - expected_foo;
    const double expected_baz = 0.4 * (1.0 - expected_foo);
    const double actual[] = { state[0], state[1], state[2] };
    const double expected[] = { expected_foo, expected_bar, expected_baz };
    for (size_t i = 0; i < 3u; ++i) {
      if (!isfinite(actual[i]) || fabs(actual[i] - expected[i]) > tolerance * fabs(expected[i])) {
        fprintf(stderr, "%s surface step %zu species %zu: %.17g vs %.17g\n",
            method, step, i, actual[i], expected[i]);
        return 1;
      }
    }
  }
  return 0;
}

static const chem_species hires_species[] = {
  { "Y0", 0.0, 0.0, 0.0, 1.0e-8, 0u }, { "Y1", 0.0, 0.0, 0.0, 1.0e-8, 0u },
  { "Y2", 0.0, 0.0, 0.0, 1.0e-8, 0u }, { "Y3", 0.0, 0.0, 0.0, 1.0e-8, 0u },
  { "Y4", 0.0, 0.0, 0.0, 1.0e-8, 0u }, { "Y5", 0.0, 0.0, 0.0, 1.0e-8, 0u },
  { "Y6", 0.0, 0.0, 0.0, 1.0e-8, 0u }, { "Y7", 0.0, 0.0, 0.0, 1.0e-8, 0u },
};
static const chem_parameter hires_parameters[] = {
  { "r1", 1.71 }, { "r2", 8.75 }, { "r3", 10.03 }, { "r4", 0.0007 },
  { "r5", 1.12 }, { "r6", 1.745 }, { "r7", 0.43 }, { "r8", 1.81 }, { "r9", 280.0 },
};
static const chem_term hires_reactants[] = {
  { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 }, { 3u, 1.0 }, { 4u, 1.0 },
  { 5u, 1.0 }, { 6u, 1.0 }, { 5u, 1.0 }, { 7u, 1.0 },
};
static const chem_term hires_products[] = {
  { 1u, 1.0 }, { 0u, 0.43 / 8.75 }, { 3u, 8.32 / 8.75 },
  { 0u, 8.32 / 10.03 }, { 3u, 1.71 / 10.03 }, { 0u, 1.0 },
  { 2u, 0.43 / 1.12 }, { 5u, 0.69 / 1.12 },
  { 2u, 0.035 / 1.745 }, { 5u, 1.71 / 1.745 }, { 4u, 1.0 },
  { 4u, 0.43 / 1.81 }, { 5u, 0.69 / 1.81 }, { 7u, 1.0 }, { 6u, 1.0 },
};
static const chem_reaction hires_reactions[] = {
  { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 0u } },
  { 1u, 1u, 1u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 1u } },
  { 2u, 1u, 3u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 2u } },
  { 3u, 0u, 5u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 3u } },
  { 3u, 1u, 6u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 4u } },
  { 4u, 1u, 8u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 5u } },
  { 5u, 1u, 10u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 6u } },
  { 6u, 1u, 11u, 3u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 7u } },
  { 7u, 2u, 14u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 8u } },
};
static const chem_mechanism hires_mechanism = {
  .name = "analytical HIRES", .species_count = COUNT(hires_species), .reaction_count = COUNT(hires_reactions),
  .reactant_term_count = COUNT(hires_reactants), .product_term_count = COUNT(hires_products),
  .parameter_count = COUNT(hires_parameters), .species = hires_species, .reactions = hires_reactions,
  .reactant_terms = hires_reactants, .product_terms = hires_products, .parameters = hires_parameters,
  .temperature = 272.5, .pressure = 101253.3, .relative_tolerance = 1.0e-6,
};
static const double hires_expected[2][8] = {
  { 0.000737131257332567, 0.000144248572631618, 0.000058887297409676, 0.001175651343283149,
    0.002386356198831330, 0.006238968252742796, 0.002849998395185769, 0.002850001604814231 },
  { 0.000670305503581864, 0.000130996846986347, 0.000046862231597733, 0.001044668020551705,
    0.000594883830951485, 0.001399628833942774, 0.001014492757718480, 0.004685507242281520 },
};

static int hires(const char *method) {
  const double tolerance = strcmp(method, "backward_euler") == 0 ? 1.0e-1 : 1.0e-6;
  double state[] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0057 };
  chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  solver.max_internal_steps = 100000u;
  double interval = 321.8122;
  for (size_t output = 0; output < 2u; ++output) {
    if (integrate(method, &hires_mechanism, state, interval, &solver) != 0) return 1;
    for (size_t species = 0; species < 8u; ++species) {
      if (!near(state[species], hires_expected[output][species], tolerance)) {
        fprintf(stderr, "%s HIRES output %zu species %zu: %.17g vs %.17g\n",
            method, output, species, state[species], hires_expected[output][species]);
        return 1;
      }
    }
    interval += 100.0;
  }
  return 0;
}

static const chem_species robertson_species[] = {
  { "A", 0.0, 0.0, 0.0, 1.0e-12, 0u },
  { "B", 0.0, 0.0, 0.0, 1.0e-12, 0u },
  { "C", 0.0, 0.0, 0.0, 1.0e-12, 0u },
};
static const chem_parameter robertson_parameters[] = {
  { "r1", 0.04 }, { "r2", 3.0e7 }, { "r3", 1.0e4 },
};
static const chem_term robertson_reactants[] = {
  { 0u, 1.0 }, { 1u, 1.0 }, { 1u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 },
};
static const chem_term robertson_products[] = {
  { 1u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 }, { 0u, 1.0 }, { 2u, 1.0 },
};
static const chem_reaction robertson_reactions[] = {
  { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 0u } },
  { 1u, 2u, 1u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 1u } },
  { 3u, 2u, 3u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 2u } },
};
static const chem_mechanism robertson_mechanism = {
  .name = "analytical Robertson system", .species_count = COUNT(robertson_species),
  .reaction_count = COUNT(robertson_reactions), .reactant_term_count = COUNT(robertson_reactants),
  .product_term_count = COUNT(robertson_products), .parameter_count = COUNT(robertson_parameters),
  .species = robertson_species, .reactions = robertson_reactions, .reactant_terms = robertson_reactants,
  .product_terms = robertson_products, .parameters = robertson_parameters,
  .temperature = 272.5, .pressure = 101253.3, .relative_tolerance = 1.0e-10,
};
static const double robertson_expected[13][3] = {
  { 1.0, 0.0, 0.0 },
  { 0.9664597373330035, 0.00003074626578578675, 0.03350951640121071 },
  { 0.8413699238414729, 0.00001623390937990473, 0.1586138422491472 },
  { 0.6172348823960878, 0.000006153591274639123, 0.3827589640126376 },
  { 0.3368745306607069, 0.000002013702318261393, 0.6631234556369748 },
  { 0.1073004285378040, 0.0000004800166972571660, 0.8926990914454987 },
  { 0.01786592114209946, 0.00000007274751468436319, 0.9821340061103859 },
  { 0.002031483924973415, 0.000000008142277783356159, 0.9979685079327488 },
  { 0.0002076093439016395, 0.0000000008306077485067610, 0.9997923898254906 },
  { 0.00002082417512179460, 0.00000000008329841429908955, 0.9999791757415798 },
  { 0.000002083229471647004, 0.000000000008332935037760723, 0.9999979167621954 },
  { 0.0000002083328471883087, 0.0000000000008333315602809495, 0.9999997916663195 },
  { 0.00000002083340149701284, 0.00000000000008333360770334744, 0.9999999791665152 },
};

static int robertson(const char *method) {
  double state[] = { 1.0, 0.0, 0.0 };
  chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  solver.max_internal_steps = 100000u;
  double target_time = 1.0, current_time = 0.0;
  for (size_t output = 1u; output < COUNT(robertson_expected); ++output) {
    if (integrate(method, &robertson_mechanism, state, target_time - current_time, &solver) != 0) return 1;
    for (size_t species = 0; species < 3u; ++species) {
      const double difference = fabs(state[species] - robertson_expected[output][species]);
      const double relative = difference / fabs(robertson_expected[output][species]);
      if (!(difference < 0.3e-3 || relative < 1.0e-6)) {
        fprintf(stderr, "Rosenbrock Robertson output %zu species %zu: %.17g vs %.17g\n",
            output, species, state[species], robertson_expected[output][species]);
        return 1;
      }
    }
    current_time = target_time;
    target_time *= 10.0;
  }
  return 0;
}

static const chem_species e5_species[] = {
  { "A1", 0.0, 0.0, 0.0, 1.0e-7, 0u }, { "A2", 0.0, 0.0, 0.0, 1.0e-17, 0u },
  { "A3", 0.0, 0.0, 0.0, 1.0e-17, 0u }, { "A4", 0.0, 0.0, 0.0, 1.0e-17, 0u },
  { "A5", 0.0, 0.0, 0.0, 1.0e-7, 0u }, { "A6", 0.0, 0.0, 0.0, 1.0e-7, 0u },
};
static const chem_parameter e5_parameters[] = {
  { "r1", 7.89e-10 }, { "r2", 1.13e9 }, { "r3", 1.1e7 }, { "r4", 1.13e3 },
};
static const chem_term e5_reactants[] = {
  { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 }, { 0u, 1.0 }, { 2u, 1.0 }, { 3u, 1.0 },
};
static const chem_term e5_products[] = {
  { 1u, 1.0 }, { 2u, 1.0 }, { 4u, 1.0 }, { 3u, 1.0 }, { 2u, 1.0 }, { 5u, 1.0 },
};
static const chem_reaction e5_reactions[] = {
  { 0u, 1u, 0u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 0u } },
  { 1u, 2u, 2u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 1u } },
  { 3u, 2u, 3u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 2u } },
  { 5u, 1u, 4u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 3u } },
};
static const chem_mechanism e5_mechanism = {
  .name = "analytical E5 system", .species_count = COUNT(e5_species), .reaction_count = COUNT(e5_reactions),
  .reactant_term_count = COUNT(e5_reactants), .product_term_count = COUNT(e5_products),
  .parameter_count = COUNT(e5_parameters), .species = e5_species, .reactions = e5_reactions,
  .reactant_terms = e5_reactants, .product_terms = e5_products, .parameters = e5_parameters,
  .temperature = 272.5, .pressure = 101253.3, .relative_tolerance = 1.0e-13,
};
static const double e5_expected[8][4] = {
  { 1.76e-3, 0.0, 0.0, 0.0 },
  { 1.7599259497677897058e-3, 1.3846281519376516449e-11, 7.6370038530073911180e-13, 1.3082581134075777338e-11 },
  { 1.6180769999072942552e-3, 1.3822370304983735443e-10, 8.2515735006838336088e-12, 1.2997212954915352082e-10 },
  { 7.4813208224292220114e-6, 2.3734781561205975019e-12, 2.2123586689581663654e-12, 1.6111948716243113653e-13 },
  { 4.7150333630401632232e-10, 1.8188895860807021729e-14, 1.8188812376786725407e-14, 8.3484020296321693074e-20 },
  { 3.1317148329356996037e-14, 1.4840957952870064294e-16, 1.4840957948345691466e-16, 4.5243728279782625194e-26 },
  { 3.8139035189787091771e-49, 1.0192582567660293322e-20, 1.0192582567660293322e-20, 3.7844935507486221171e-65 },
  { 0.0, 8.8612334976263783420e-23, 8.8612334976263783421e-23, 0.0 },
};

static int e5(void) {
  double state[] = { 1.76e-3, 0.0, 0.0, 0.0, 0.0, 0.0 };
  chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  solver.max_internal_steps = 1000000u;
  double target_time = 10.0, current_time = 0.0;
  for (size_t output = 1u; output < COUNT(e5_expected); ++output) {
    if (integrate("rosenbrock", &e5_mechanism, state, target_time - current_time, &solver) != 0) return 1;
    for (size_t species = 0; species < 4u; ++species) {
      const double difference = fabs(state[species] - e5_expected[output][species]);
      const double absolute_tolerance = species == 0u ? 1.0e-6 : 1.0e-13;
      const double relative = e5_expected[output][species] == 0.0 ? 0.0 :
          difference / fabs(e5_expected[output][species]);
      if (!(difference < absolute_tolerance || relative < 1.0e-10)) {
        fprintf(stderr, "E5 output %zu species %zu: %.17g vs %.17g\n",
            output, species, state[species], e5_expected[output][species]);
        return 1;
      }
    }
    current_time = target_time;
    target_time *= 100.0;
  }
  return 0;
}

static const double oregonator_dimensionless[13][3] = {
  { 1.0, 2.0, 3.0 },
  { 1.000661467180497, 1512.778937348249, 10358.54312767229 },
  { 1.000874625199626, 1144.336972384497, 83.7214996662463 },
  { 1.001890368438751, 529.9926232295553, 1.662279579042420 },
  { 1.004118022612645, 243.8326079910346, 1.008822224048647 },
  { 1.008995416634061, 112.1664388660253, 1.007783229065319 },
  { 1.019763472537298, 51.59761322947535, 1.016985778956374 },
  { 1.043985088527474, 23.73442027531524, 1.037691843544522 },
  { 1.100849071667922, 10.91533805469022, 1.085831969810860 },
  { 1.249102130020572, 5.013945178605446, 1.208326626237875 },
  { 1.779724751930701, 2.281852385542403, 1.613754023671725 },
  { 1.000889326903503, 1125.438585747168, 16410.49483777168 },
  { 1.000814870318523, 1228.178521546513, 132.0554942846513 },
};

static int oregonator(const char *method) {
  const double s = 77.27, q = 8.375e-6, w = 0.161;
  const double k2 = 1.6e9, k5 = 1.0;
  const double k1 = k5 / (s * w), k3 = s * k5 / w, k4 = s * s * q * k2 / 2.0;
  const double alpha = k5 / (s * w * k2);
  const double eta = s * k5 / (w * k2);
  const double rho = k5 * k5 / (w * w * k2);
  const double tau = w / k5;
  chem_species species[] = {
    { "X", 0.0, 0.0, 0.0, alpha * 1.0e-11, 0u },
    { "Y", 0.0, 0.0, 0.0, eta * 1.0e-11, 0u },
    { "Z", 0.0, 0.0, 0.0, rho * 1.0e-11, 0u },
    { "P", 0.0, 0.0, 0.0, eta, 0u },
    { "Q", 0.0, 0.0, 0.0, eta, 0u },
  };
  const chem_parameter parameters[] = {
    { "r1", k1 }, { "r2", k2 }, { "r3", k3 }, { "r4", k4 }, { "r5", k5 },
  };
  const chem_term reactants[] = {
    { 1u, 1.0 }, { 0u, 1.0 }, { 1u, 1.0 }, { 0u, 1.0 }, { 0u, 1.0 }, { 0u, 1.0 }, { 0u, 1.0 }, { 2u, 1.0 },
  };
  const chem_term products[] = {
    { 0u, 1.0 }, { 3u, 1.0 }, { 2u, 1.0 }, { 0u, 2.0 }, { 4u, 1.0 }, { 1u, 1.0 },
  };
  const chem_reaction reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 0u } },
    { 1u, 2u, 1u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 1u } },
    { 3u, 1u, 2u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 2u } },
    { 4u, 2u, 4u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 3u } },
    { 7u, 1u, 5u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 4u } },
  };
  const chem_mechanism mechanism = {
    .name = "analytical Oregonator", .species_count = COUNT(species), .reaction_count = COUNT(reactions),
    .reactant_term_count = COUNT(reactants), .product_term_count = COUNT(products),
    .parameter_count = COUNT(parameters), .species = species, .reactions = reactions,
    .reactant_terms = reactants, .product_terms = products, .parameters = parameters,
    .temperature = 298.15, .pressure = 101325.0, .relative_tolerance = 1.0e-9,
  };
  double state[] = { alpha, 2.0 * eta, 3.0 * rho, 0.0, 0.0 };
  chem_rosenbrock_parameters solver = chem_rosenbrock_four_stage_parameters();
  solver.max_internal_steps = 1000000u;
  solver.lu_method = CHEM_LU_DOOLITTLE;
  const int backward_euler = strcmp(method, "backward_euler") == 0;
  const size_t substeps = backward_euler ? 18000u : 1u;
  const double tolerance = backward_euler ? 5.0e-3 : 1.0e-6;
  const double interval = 30.0 * tau;
  for (size_t output = 1u; output < COUNT(oregonator_dimensionless); ++output) {
    const double dt = interval / (double)substeps;
    for (size_t substep = 0; substep < substeps; ++substep)
      if (integrate(method, &mechanism, state, dt, &solver) != 0) {
        fprintf(stderr, "%s Oregonator failed at output %zu substep %zu\n", method, output, substep);
        return 1;
      }
    const double expected[] = {
      oregonator_dimensionless[output][0] * alpha,
      oregonator_dimensionless[output][1] * eta,
      oregonator_dimensionless[output][2] * rho,
    };
    for (size_t species_index = 0; species_index < 3u; ++species_index) {
      const double relative = fabs(state[species_index] - expected[species_index]) / fabs(expected[species_index]);
      if (!isfinite(relative) || relative >= tolerance) {
        fprintf(stderr, "%s Oregonator output %zu species %zu: %.17g vs %.17g rel %.3g\n",
            method, output, species_index, state[species_index], expected[species_index], relative);
        return 1;
      }
    }
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc != 3) return 2;
  const char *method = argv[1];
  const char *test = strrchr(argv[2], '.');
  test = test ? test + 1 : argv[2];
  if (strcmp(method, "backward_euler") != 0 && strcmp(method, "rosenbrock") != 0) return 2;
  if (strcmp(test, "SurfaceRxn") == 0) return surface_reaction(method);
  if (strcmp(test, "HIRES") == 0) return hires(method);
  if (strcmp(test, "Oregonator") == 0) return oregonator(method);
  if (strcmp(test, "Robertson") == 0 && strcmp(method, "rosenbrock") == 0) return robertson(method);
  if (strcmp(test, "E5") == 0 && strcmp(method, "rosenbrock") == 0) return e5();
  return 2;
}


