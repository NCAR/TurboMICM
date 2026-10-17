/* MICM's hard-coded Chapman solver and its regression tests' Chapman
 * mechanism, see chapman_fixed.h. Force, DforceDy, the factorization, the
 * triangular solves, the error norm, and the Rosenbrock loop are transcribed
 * from test/regression/RosenbrockChapman/chapman_ode_solver.hpp. */
#include "chapman_fixed.h"

#include <float.h>

enum
{
  M,
  Ar,
  CO2,
  H2O,
  N2,
  O1D,
  O,
  O2,
  O3
};

const char *const cf_species_names[CF_N_SPECIES] = { "M", "Ar", "CO2", "H2O", "N2", "O1D", "O", "O2", "O3" };
const char *const cf_photolysis_names[CF_N_PHOTOLYSIS] = { "O2_1", "O3_1", "O3_2" };
const int cf_elements[CF_N_ELEMENTS][2] = { { M, M },     { O, M },     { O2, M },   { O3, M },   { Ar, Ar },  { CO2, CO2 },
                                            { H2O, H2O }, { N2, N2 },   { O1D, N2 }, { O, N2 },   { O1D, O1D }, { O, O1D },
                                            { O, O },     { O2, O },    { O3, O },   { O1D, O2 }, { O, O2 },    { O2, O2 },
                                            { O3, O2 },   { O1D, O3 },  { O, O3 },   { O2, O3 },  { O3, O3 } };

/* ---- the mechanism's algebra ---------------------------------------------- */

void cf_force(const double *rate_constants, const double *number_densities, double *force)
{
  for (int i = 0; i < CF_N_SPECIES; ++i)
    force[i] = 0;
  /* M, Ar, CO2, H2O, N2 are all zero */
  /* O1D */
  force[5] = force[5] + rate_constants[5] * number_densities[8];
  force[5] = force[5] - rate_constants[0] * number_densities[4] * number_densities[5];
  force[5] = force[5] - rate_constants[1] * number_densities[5] * number_densities[7];
  /* O */
  force[6] = force[6] + 2 * rate_constants[4] * number_densities[7];
  force[6] = force[6] + rate_constants[6] * number_densities[8];
  force[6] = force[6] + rate_constants[0] * number_densities[4] * number_densities[5];
  force[6] = force[6] + rate_constants[1] * number_densities[5] * number_densities[7];
  force[6] = force[6] - rate_constants[2] * number_densities[6] * number_densities[8];
  force[6] = force[6] - rate_constants[3] * number_densities[0] * number_densities[6] * number_densities[7];
  /* O2 */
  force[7] = force[7] - rate_constants[4] * number_densities[7];
  force[7] = force[7] + rate_constants[5] * number_densities[8];
  force[7] = force[7] + rate_constants[6] * number_densities[8];
  force[7] = force[7] + 2 * rate_constants[2] * number_densities[6] * number_densities[8];
  force[7] = force[7] - rate_constants[3] * number_densities[0] * number_densities[6] * number_densities[7];
  /* O3 */
  force[8] = force[8] - rate_constants[5] * number_densities[8];
  force[8] = force[8] - rate_constants[6] * number_densities[8];
  force[8] = force[8] - rate_constants[2] * number_densities[6] * number_densities[8];
  force[8] = force[8] + rate_constants[3] * number_densities[0] * number_densities[6] * number_densities[7];
}

void cf_dforce_dy(const double *rate_constants, const double *number_densities, double *jacobian)
{
  for (int e = 0; e < CF_N_ELEMENTS; ++e)
    jacobian[e] = 0;
  jacobian[1] = jacobian[1] - rate_constants[3] * number_densities[6] * number_densities[7];
  jacobian[2] = jacobian[2] - rate_constants[3] * number_densities[6] * number_densities[7];
  jacobian[3] = jacobian[3] + rate_constants[3] * number_densities[6] * number_densities[7];
  jacobian[8] = jacobian[8] - rate_constants[0] * number_densities[5];
  jacobian[9] = jacobian[9] + rate_constants[0] * number_densities[5];
  jacobian[10] = jacobian[10] - rate_constants[0] * number_densities[4];
  jacobian[10] = jacobian[10] - rate_constants[1] * number_densities[7];
  jacobian[11] = jacobian[11] + rate_constants[0] * number_densities[4];
  jacobian[11] = jacobian[11] + rate_constants[1] * number_densities[7];
  jacobian[12] = jacobian[12] - rate_constants[2] * number_densities[8];
  jacobian[12] = jacobian[12] - rate_constants[3] * number_densities[0] * number_densities[7];
  jacobian[13] = jacobian[13] + 2 * rate_constants[2] * number_densities[8];
  jacobian[13] = jacobian[13] - rate_constants[3] * number_densities[0] * number_densities[7];
  jacobian[14] = jacobian[14] - rate_constants[2] * number_densities[8];
  jacobian[14] = jacobian[14] + rate_constants[3] * number_densities[0] * number_densities[7];
  jacobian[15] = jacobian[15] - rate_constants[1] * number_densities[5];
  jacobian[16] = jacobian[16] + 2 * rate_constants[4];
  jacobian[16] = jacobian[16] + rate_constants[1] * number_densities[5];
  jacobian[16] = jacobian[16] - rate_constants[3] * number_densities[0] * number_densities[6];
  jacobian[17] = jacobian[17] - rate_constants[4];
  jacobian[17] = jacobian[17] - rate_constants[3] * number_densities[0] * number_densities[6];
  jacobian[18] = jacobian[18] + rate_constants[3] * number_densities[0] * number_densities[6];
  jacobian[19] = jacobian[19] + rate_constants[5];
  jacobian[20] = jacobian[20] + rate_constants[6];
  jacobian[20] = jacobian[20] - rate_constants[2] * number_densities[6];
  jacobian[21] = jacobian[21] + rate_constants[5];
  jacobian[21] = jacobian[21] + rate_constants[6];
  jacobian[21] = jacobian[21] + 2 * rate_constants[2] * number_densities[6];
  jacobian[22] = jacobian[22] - rate_constants[5];
  jacobian[22] = jacobian[22] - rate_constants[6];
  jacobian[22] = jacobian[22] - rate_constants[2] * number_densities[6];
}

/* ChapmanODESolver::Factor: an in-place LU with the pivot reciprocals on the diagonal */
static void factor(double *jacobian)
{
  jacobian[0] = 1. / jacobian[0];
  jacobian[1] = jacobian[1] * jacobian[0];
  jacobian[2] = jacobian[2] * jacobian[0];
  jacobian[3] = jacobian[3] * jacobian[0];
  jacobian[4] = 1. / jacobian[4];
  jacobian[5] = 1. / jacobian[5];
  jacobian[6] = 1. / jacobian[6];
  jacobian[7] = 1. / jacobian[7];
  jacobian[8] = jacobian[8] * jacobian[7];
  jacobian[9] = jacobian[9] * jacobian[7];
  jacobian[10] = 1. / jacobian[10];
  jacobian[11] = jacobian[11] * jacobian[10];
  jacobian[16] = jacobian[16] - jacobian[11] * jacobian[15];
  jacobian[20] = jacobian[20] - jacobian[11] * jacobian[19];
  jacobian[12] = 1. / jacobian[12];
  jacobian[13] = jacobian[13] * jacobian[12];
  jacobian[14] = jacobian[14] * jacobian[12];
  jacobian[17] = jacobian[17] - jacobian[13] * jacobian[16];
  jacobian[18] = jacobian[18] - jacobian[14] * jacobian[16];
  jacobian[21] = jacobian[21] - jacobian[13] * jacobian[20];
  jacobian[22] = jacobian[22] - jacobian[14] * jacobian[20];
  jacobian[17] = 1. / jacobian[17];
  jacobian[18] = jacobian[18] * jacobian[17];
  jacobian[22] = jacobian[22] - jacobian[18] * jacobian[21];
  jacobian[22] = 1. / jacobian[22];
}

/* ChapmanODESolver::FactoredAlphaMinusJac */
static void factored_alpha_minus_jac(const double *dforce_dy, double alpha, double *jacobian)
{
  for (int e = 0; e < CF_N_ELEMENTS; ++e)
    jacobian[e] = -dforce_dy[e];
  jacobian[0] = -dforce_dy[0] + alpha;
  jacobian[4] = -dforce_dy[4] + alpha;
  jacobian[5] = -dforce_dy[5] + alpha;
  jacobian[6] = -dforce_dy[6] + alpha;
  jacobian[7] = -dforce_dy[7] + alpha;
  jacobian[10] = -dforce_dy[10] + alpha;
  jacobian[12] = -dforce_dy[12] + alpha;
  jacobian[17] = -dforce_dy[17] + alpha;
  jacobian[22] = -dforce_dy[22] + alpha;
  factor(jacobian);
}

/* ChapmanODESolver::BacksolveLYEqB and BacksolveUXEqB */
static void lin_solve(const double *jacobian, const double *b, double *x)
{
  double y[CF_N_SPECIES];
  y[0] = b[0];
  y[1] = b[1];
  y[2] = b[2];
  y[3] = b[3];
  y[4] = b[4];
  y[5] = b[5];
  y[5] = y[5] - jacobian[8] * y[4];
  y[6] = b[6];
  y[6] = y[6] - jacobian[1] * y[0];
  y[6] = y[6] - jacobian[9] * y[4];
  y[6] = y[6] - jacobian[11] * y[5];
  y[7] = b[7];
  y[7] = y[7] - jacobian[2] * y[0];
  y[7] = y[7] - jacobian[13] * y[6];
  y[8] = b[8];
  y[8] = y[8] - jacobian[3] * y[0];
  y[8] = y[8] - jacobian[14] * y[6];
  y[8] = y[8] - jacobian[18] * y[7];

  double temporary;
  temporary = y[8];
  x[8] = jacobian[22] * temporary;
  temporary = y[7];
  temporary = temporary - jacobian[21] * x[8];
  x[7] = jacobian[17] * temporary;
  temporary = y[6];
  temporary = temporary - jacobian[16] * x[7];
  temporary = temporary - jacobian[20] * x[8];
  x[6] = jacobian[12] * temporary;
  temporary = y[5];
  temporary = temporary - jacobian[15] * x[7];
  temporary = temporary - jacobian[19] * x[8];
  x[5] = jacobian[10] * temporary;
  temporary = y[4];
  x[4] = jacobian[7] * temporary;
  temporary = y[3];
  x[3] = jacobian[6] * temporary;
  temporary = y[2];
  x[2] = jacobian[5] * temporary;
  temporary = y[1];
  x[1] = jacobian[4] * temporary;
  temporary = y[0];
  x[0] = jacobian[0] * temporary;
}

void cf_update_state(double temperature, double pressure, const double *photolysis_rates, double *rate_constants)
{
  tmc_conditions conditions = { temperature, pressure, 0.0 };
  tmc_arrhenius_parameters params = tmc_arrhenius_defaults();
  tmc_rate_constant rc;
  /* N2_O1D_1: N2 + O1D -> O + N2 */
  params.A = 2.15e-11;
  params.C = 110;
  rc = tmc_rate_arrhenius(params);
  rate_constants[0] = tmc_rate_constant_calculate(&rc, &conditions, NULL);
  /* O1D_O2_1: O1D + O2 -> O + O2 */
  params.A = 3.3e-11;
  params.C = 55;
  rc = tmc_rate_arrhenius(params);
  rate_constants[1] = tmc_rate_constant_calculate(&rc, &conditions, NULL);
  /* O_O3_1: O + O3 -> 2 O2 */
  params.A = 8e-12;
  params.C = -2060;
  rc = tmc_rate_arrhenius(params);
  rate_constants[2] = tmc_rate_constant_calculate(&rc, &conditions, NULL);
  /* M_O_O2_1: M + O + O2 -> O3 + M */
  params.A = 6e-34;
  params.B = 2.4;
  params.C = 0;
  rc = tmc_rate_arrhenius(params);
  rate_constants[3] = tmc_rate_constant_calculate(&rc, &conditions, NULL);
  /* O2_1, O3_1, O3_2: the photolysis rates */
  rate_constants[4] = photolysis_rates[0];
  rate_constants[5] = photolysis_rates[1];
  rate_constants[6] = photolysis_rates[2];
}

/* ---- ChapmanODESolver::Solve ---------------------------------------------- */

/* ChapmanODESolver::ErrorNorm, with its tolerances */
static double error_norm(const double *Y, const double *Ynew, const double *errors)
{
  const double absolute_tolerance = 1e-3, relative_tolerance = 1e-4;
  double sum = 0;
  for (int i = 0; i < CF_N_SPECIES; ++i)
  {
    const double scale = absolute_tolerance + relative_tolerance * fmax(fabs(Y[i]), fabs(Ynew[i]));
    sum += pow(errors[i] / scale, 2);
  }
  const double error_min = 1.0e-10;
  return fmax(sqrt(sum / CF_N_SPECIES), error_min);
}

cf_result cf_solve(double time_start, double time_end, const double *y, const double *rate_constants)
{
  /* ChapmanParameters with ThreeStageRosenbrock() */
  enum
  {
    STAGES = 3,
    N = CF_N_SPECIES
  };
  static const double a[3] = { 1, 1, 0 };
  static const double c[3] = { -0.10156171083877702091975600115545e+01, 0.40759956452537699824805835358067e+01, 0.92076794298330791242156818474003e+01 };
  static const bool new_function_evaluation[3] = { true, true, false };
  static const double m[3] = { 0.1e+01, 0.61697947043828245592553615689730e+01, -0.42772256543218573326238373806514 };
  static const double e[3] = { 0.5, -0.29079558716805469821718236208017e+01, 0.22354069897811569627360909276199 };
  const double estimator_of_local_order = 3;
  const double gamma_0 = 0.43586652150845899941601945119356;
  const int max_number_of_steps = 100;
  const double round_off = DBL_EPSILON, factor_min = 0.2, factor_max = 6, rejection_factor_decrease = 0.1, safety_factor = 0.9;
  const double h_min = 0, DELTA_MIN = 1.0e-6;
  const double h_max = time_end - time_start;
  const double h_start = fmax(h_min, DELTA_MIN);

  cf_result result;
  memset(&result, 0, sizeof result);
  double Y[N], K[STAGES][N], forcing[N], initial_forcing[N], new_Y[N], Yerror[N], dforce_dy[CF_N_ELEMENTS], ode_jacobian[CF_N_ELEMENTS];
  memcpy(Y, y, sizeof Y);

  double present_time = time_start;
  double H = fmin(fmax(fabs(h_min), fabs(h_start)), fabs(h_max));
  if (fabs(H) <= 10 * round_off)
    H = DELTA_MIN;
  bool reject_last_h = false, reject_more_h = false;

  while ((present_time - time_end + round_off) <= 0)
  {
    if (result.number_of_steps > max_number_of_steps)
      break; /* ConvergenceExceededMaxSteps */
    if (((present_time + 0.1 * H) == present_time) || (H <= round_off))
      break; /* StepSizeTooSmall */
    H = fmin(H, fabs(time_end - present_time));

    cf_force(rate_constants, Y, initial_forcing);
    result.function_calls += 1;

    bool accepted = false;
    while (!accepted)
    {
      if (result.number_of_steps > max_number_of_steps)
        break;
      /* LinFactor: the factored alpha I - J */
      const double alpha = 1 / (H * gamma_0);
      cf_dforce_dy(rate_constants, Y, dforce_dy);
      factored_alpha_minus_jac(dforce_dy, alpha, ode_jacobian);
      result.decompositions += 1;
      result.jacobian_updates += 1;

      for (int stage = 0; stage < STAGES; ++stage)
      {
        if (stage == 0)
        {
          memcpy(K[0], initial_forcing, sizeof K[0]);
          memcpy(forcing, initial_forcing, sizeof forcing);
        }
        else
        {
          const int stage_combinations = ((stage + 1) - 1) * ((stage + 1) - 2) / 2;
          if (new_function_evaluation[stage])
          {
            memcpy(new_Y, Y, sizeof new_Y);
            for (int j = 0; j < stage; ++j)
              for (int idx = 0; idx < N; ++idx)
                new_Y[idx] += a[stage_combinations + j] * K[j][idx];
            cf_force(rate_constants, new_Y, forcing);
            result.function_calls += 1;
          }
          memcpy(K[stage], forcing, sizeof K[stage]);
          for (int j = 0; j < stage; ++j)
          {
            const double HC = c[stage_combinations + j] / H;
            for (int idx = 0; idx < N; ++idx)
              K[stage][idx] += HC * K[j][idx];
          }
        }
        double solved[N];
        lin_solve(ode_jacobian, K[stage], solved);
        memcpy(K[stage], solved, sizeof solved);
        result.solves += 1;
      }

      memcpy(new_Y, Y, sizeof new_Y);
      for (int stage = 0; stage < STAGES; ++stage)
        for (int idx = 0; idx < N; ++idx)
          new_Y[idx] += m[stage] * K[stage][idx];
      for (int idx = 0; idx < N; ++idx)
        Yerror[idx] = 0;
      for (int stage = 0; stage < STAGES; ++stage)
        for (int idx = 0; idx < N; ++idx)
          Yerror[idx] += e[stage] * K[stage][idx];
      const double error = error_norm(Y, new_Y, Yerror);

      double Hnew = H * fmin(factor_max, fmax(factor_min, safety_factor / pow(error, 1 / estimator_of_local_order)));
      result.number_of_steps += 1;
      if ((error < 1) || (H < h_min))
      {
        result.accepted += 1;
        present_time = present_time + H;
        memcpy(Y, new_Y, sizeof Y);
        Hnew = fmax(h_min, fmin(Hnew, h_max));
        if (reject_last_h)
          Hnew = fmin(Hnew, H);
        reject_last_h = false;
        reject_more_h = false;
        H = Hnew;
        accepted = true;
      }
      else
      {
        if (reject_more_h)
          Hnew = H * rejection_factor_decrease;
        reject_more_h = reject_last_h;
        reject_last_h = true;
        H = Hnew;
        if (result.accepted >= 1)
          result.rejected += 1;
      }
    }
  }
  result.final_time = present_time;
  memcpy(result.y, Y, sizeof result.y);
  return result;
}

/* ---- util.hpp: the mechanism for the general solver ------------------------ */

static tmc_rate_constant arrhenius(double A, double B, double C)
{
  tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
  p.A = A;
  p.B = B;
  p.C = C;
  return tmc_rate_arrhenius(p);
}

tmc_solver *cf_chapman_solver(tmc_builder *builder)
{
  tmc_species s[CF_N_SPECIES];
  tmc_phase gas = tmc_phase_create("gas");
  for (int i = 0; i < CF_N_SPECIES; ++i)
  {
    s[i] = tmc_species_create(cf_species_names[i]);
    tmc_phase_add_species(&gas, &s[i]);
  }
  tmc_reaction r[7];
  for (int i = 0; i < 7; ++i)
  {
    r[i] = tmc_reaction_create();
    tmc_reaction_set_phase(&r[i], &gas);
  }
  /* photo_1: O2 -> 2 O */
  tmc_reaction_add_reactant(&r[0], &s[O2]);
  tmc_reaction_add_product(&r[0], &s[O], 2);
  tmc_reaction_set_rate_constant(&r[0], tmc_rate_user_defined("jO2", 1.0));
  /* photo_2: O3 -> O1D + O2 */
  tmc_reaction_add_reactant(&r[1], &s[O3]);
  tmc_reaction_add_product(&r[1], &s[O1D], 1);
  tmc_reaction_add_product(&r[1], &s[O2], 1);
  tmc_reaction_set_rate_constant(&r[1], tmc_rate_user_defined("jO3a", 1.0));
  /* photo_3: O3 -> O + O2 */
  tmc_reaction_add_reactant(&r[2], &s[O3]);
  tmc_reaction_add_product(&r[2], &s[O], 1);
  tmc_reaction_add_product(&r[2], &s[O2], 1);
  tmc_reaction_set_rate_constant(&r[2], tmc_rate_user_defined("jO3b", 1.0));
  /* r1: O1D + N2 -> O + N2 */
  tmc_reaction_add_reactant(&r[3], &s[O1D]);
  tmc_reaction_add_reactant(&r[3], &s[N2]);
  tmc_reaction_add_product(&r[3], &s[O], 1);
  tmc_reaction_add_product(&r[3], &s[N2], 1);
  tmc_reaction_set_rate_constant(&r[3], arrhenius(2.15e-11, 0, 110));
  /* r2: O1D + O2 -> O + O2 */
  tmc_reaction_add_reactant(&r[4], &s[O1D]);
  tmc_reaction_add_reactant(&r[4], &s[O2]);
  tmc_reaction_add_product(&r[4], &s[O], 1);
  tmc_reaction_add_product(&r[4], &s[O2], 1);
  tmc_reaction_set_rate_constant(&r[4], arrhenius(3.3e-11, 0, 55));
  /* r3: O + O3 -> 2 O2 */
  tmc_reaction_add_reactant(&r[5], &s[O]);
  tmc_reaction_add_reactant(&r[5], &s[O3]);
  tmc_reaction_add_product(&r[5], &s[O2], 2);
  tmc_reaction_set_rate_constant(&r[5], arrhenius(8e-12, 0, -2060));
  /* r4: O + O2 + M -> O3 + M */
  tmc_reaction_add_reactant(&r[6], &s[O]);
  tmc_reaction_add_reactant(&r[6], &s[O2]);
  tmc_reaction_add_reactant(&r[6], &s[M]);
  tmc_reaction_add_product(&r[6], &s[O3], 1);
  tmc_reaction_add_product(&r[6], &s[M], 1);
  tmc_reaction_set_rate_constant(&r[6], arrhenius(6.0e-34, 2.4, 0));

  tmc_system system = tmc_system_create(&gas);
  tmc_builder_set_system(builder, &system);
  tmc_builder_set_reactions(builder, r, 7);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_solver *solver = tmc_builder_build(builder, &err);
  if (!solver)
    tmc_test_fail(__FILE__, __LINE__, "Chapman build failed: %s", err.message);
  tmc_builder_free(builder);
  for (int i = 0; i < 7; ++i)
    tmc_reaction_free(&r[i]);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return solver;
}
