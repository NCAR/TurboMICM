/* Rosenbrock parameter sets (MICM's rosenbrock_solver_parameters.hpp),
 * backward Euler defaults, and the default step-size controller. */
#include "internal.h"

#include <float.h>
#include <math.h>

static tmc_rosenbrock_parameters base(void)
{
  tmc_rosenbrock_parameters p;
  memset(&p, 0, sizeof p);
  p.max_number_of_steps = 1000;
  p.round_off = DBL_EPSILON;
  p.factor_min = 0.2;
  p.factor_max = 6.0;
  p.rejection_factor_decrease = 0.1;
  p.safety_factor = 0.9;
  p.constraint_init_max_iterations = 10;
  p.constraint_init_tolerance = 0.1;
  p.constraint_init_max_backtracks = 24;
  p.constraint_init_backtrack_factor = 0.5;
  p.constraint_init_sufficient_decrease = 1e-4;
  return p;
}

tmc_rosenbrock_parameters tmc_rosenbrock_two_stage(void)
{
  /* an L-stable method, 2 stages, order 2 */
  tmc_rosenbrock_parameters p = base();
  double g = 1.0 + 1.0 / sqrt(2.0);
  p.stages = 2;
  p.a[0] = 1.0 / g;
  p.c[0] = -2.0 / g;
  for (int i = 0; i < TMC_MAX_STAGES; ++i)
    p.new_function_evaluation[i] = true;
  p.m[0] = 3.0 / (2.0 * g);
  p.m[1] = 1.0 / (2.0 * g);
  p.e[0] = 1.0 / (2.0 * g);
  p.e[1] = 1.0 / (2.0 * g);
  p.estimator_of_local_order = 2.0;
  p.alpha[0] = 0.0;
  p.alpha[1] = 1.0;
  p.gamma[0] = g;
  p.gamma[1] = -g;
  return p;
}

tmc_rosenbrock_parameters tmc_rosenbrock_three_stage(void)
{
  /* an L-stable method, 3 stages, order 3, 2 function evaluations
   * (Sandu et al., 1997, Atmospheric Environment 31, 3459-3472) */
  tmc_rosenbrock_parameters p = base();
  p.stages = 3;
  p.a[0] = 1;
  p.a[1] = 1;
  p.a[2] = 0;
  p.c[0] = -0.10156171083877702091975600115545e+01;
  p.c[1] = 0.40759956452537699824805835358067e+01;
  p.c[2] = 0.92076794298330791242156818474003e+01;
  p.new_function_evaluation[0] = true;
  p.new_function_evaluation[1] = true;
  p.new_function_evaluation[2] = false;
  p.m[0] = 0.1e+01;
  p.m[1] = 0.61697947043828245592553615689730e+01;
  p.m[2] = -0.42772256543218573326238373806514;
  p.e[0] = 0.5;
  p.e[1] = -0.29079558716805469821718236208017e+01;
  p.e[2] = 0.22354069897811569627360909276199;
  p.estimator_of_local_order = 3;
  p.alpha[0] = 0;
  p.alpha[1] = 0.43586652150845899941601945119356;
  p.alpha[2] = 0.43586652150845899941601945119356;
  p.gamma[0] = 0.43586652150845899941601945119356;
  p.gamma[1] = 0.24291996454816804366592249683314;
  p.gamma[2] = 0.21851380027664058511513169485832e+01;
  return p;
}

tmc_rosenbrock_parameters tmc_rosenbrock_four_stage(void)
{
  /* L-stable Rosenbrock method of order 4 with 4 stages, embedded order 3
   * (Hairer and Wanner, 1990) */
  tmc_rosenbrock_parameters p = base();
  p.stages = 4;
  p.a[0] = 0.2000000000000000E+01;
  p.a[1] = 0.1867943637803922E+01;
  p.a[2] = 0.2344449711399156;
  p.a[3] = p.a[1];
  p.a[4] = p.a[2];
  p.a[5] = 0.0;
  p.c[0] = -0.7137615036412310E+01;
  p.c[1] = 0.2580708087951457E+01;
  p.c[2] = 0.6515950076447975;
  p.c[3] = -0.2137148994382534E+01;
  p.c[4] = -0.3214669691237626;
  p.c[5] = -0.6949742501781779;
  p.new_function_evaluation[0] = true;
  p.new_function_evaluation[1] = true;
  p.new_function_evaluation[2] = true;
  p.new_function_evaluation[3] = false;
  p.m[0] = 0.2255570073418735E+01;
  p.m[1] = 0.2870493262186792;
  p.m[2] = 0.4353179431840180;
  p.m[3] = 0.1093502252409163E+01;
  p.e[0] = -0.2815431932141155;
  p.e[1] = -0.7276199124938920E-01;
  p.e[2] = -0.1082196201495311;
  p.e[3] = -0.1093502252409163E+01;
  p.estimator_of_local_order = 4.0;
  p.alpha[0] = 0.0;
  p.alpha[1] = 0.1145640000000000E+01;
  p.alpha[2] = 0.6552168638155900;
  p.alpha[3] = p.alpha[2];
  p.gamma[0] = 0.5728200000000000;
  p.gamma[1] = -0.1769193891319233E+01;
  p.gamma[2] = 0.7592633437920482;
  p.gamma[3] = -0.1049021087100450;
  return p;
}

tmc_rosenbrock_parameters tmc_rosenbrock_four_stage_dae(void)
{
  /* a stiffly-stable method, 4 stages, order 3 */
  tmc_rosenbrock_parameters p = base();
  p.stages = 4;
  p.a[0] = 0.0;
  p.a[1] = 2.0;
  p.a[2] = 0.0;
  p.a[3] = 2.0;
  p.a[4] = 0.0;
  p.a[5] = 1.0;
  p.c[0] = 4.0;
  p.c[1] = 1.0;
  p.c[2] = -1.0;
  p.c[3] = 1.0;
  p.c[4] = -1.0;
  p.c[5] = -(8.0 / 3.0);
  p.new_function_evaluation[0] = true;
  p.new_function_evaluation[1] = false;
  p.new_function_evaluation[2] = true;
  p.new_function_evaluation[3] = true;
  p.m[0] = 2.0;
  p.m[2] = 1.0;
  p.m[3] = 1.0;
  p.e[3] = 1.0;
  p.estimator_of_local_order = 3.0;
  p.alpha[2] = 1.0;
  p.alpha[3] = 1.0;
  p.gamma[0] = 0.5;
  p.gamma[1] = 1.5;
  return p;
}

tmc_rosenbrock_parameters tmc_rosenbrock_six_stage_dae(void)
{
  /* stiffly-stable Rosenbrock method of order 4 with 6 stages
   * (Hairer and Wanner, 1996) */
  tmc_rosenbrock_parameters p = base();
  p.stages = 6;
  p.alpha[0] = 0.000;
  p.alpha[1] = 0.386;
  p.alpha[2] = 0.210;
  p.alpha[3] = 0.630;
  p.alpha[4] = 1.000;
  p.alpha[5] = 1.000;
  p.gamma[0] = 0.2500000000000000;
  p.gamma[1] = -0.1043000000000000;
  p.gamma[2] = 0.1035000000000000;
  p.gamma[3] = -0.3620000000000023E-01;
  p.gamma[4] = 0.0;
  p.gamma[5] = 0.0;
  p.a[0] = 0.1544000000000000E+01;
  p.a[1] = 0.9466785280815826;
  p.a[2] = 0.2557011698983284;
  p.a[3] = 0.3314825187068521E+01;
  p.a[4] = 0.2896124015972201E+01;
  p.a[5] = 0.9986419139977817;
  p.a[6] = 0.1221224509226641E+01;
  p.a[7] = 0.6019134481288629E+01;
  p.a[8] = 0.1253708332932087E+02;
  p.a[9] = -0.6878860361058950;
  p.a[10] = p.a[6];
  p.a[11] = p.a[7];
  p.a[12] = p.a[8];
  p.a[13] = p.a[9];
  p.a[14] = 1.0;
  p.c[0] = -0.5668800000000000E+01;
  p.c[1] = -0.2430093356833875E+01;
  p.c[2] = -0.2063599157091915;
  p.c[3] = -0.1073529058151375;
  p.c[4] = -0.9594562251023355E+01;
  p.c[5] = -0.2047028614809616E+02;
  p.c[6] = 0.7496443313967647E+01;
  p.c[7] = -0.1024680431464352E+02;
  p.c[8] = -0.3399990352819905E+02;
  p.c[9] = 0.1170890893206160E+02;
  p.c[10] = 0.8083246795921522E+01;
  p.c[11] = -0.7981132988064893E+01;
  p.c[12] = -0.3152159432874371E+02;
  p.c[13] = 0.1631930543123136E+02;
  p.c[14] = -0.6058818238834054E+01;
  p.m[0] = p.a[6];
  p.m[1] = p.a[7];
  p.m[2] = p.a[8];
  p.m[3] = p.a[9];
  p.m[4] = 1.0;
  p.m[5] = 1.0;
  p.e[5] = 1.0;
  for (int i = 0; i < TMC_MAX_STAGES; ++i)
    p.new_function_evaluation[i] = true;
  p.estimator_of_local_order = 4.0;
  return p;
}

/* A Rosenbrock or W tableau in the standard form of Hairer and Wanner (IV.7):
 *   (I - h gamma J) k_i = h f(y0 + sum_j alpha_ij k_j) + h J sum_j gamma_ij k_j,
 *   y1 = y0 + sum_i b_i k_i, the embedded solution with bhat,
 * converted to MICM's form with the stage vectors u = Gamma k (Gamma lower
 * triangular, gamma on its diagonal): a = alpha Gamma^-1, c_ij = -(Gamma^-1)_ij
 * for i > j, m = b Gamma^-1, and e = (b - bhat) Gamma^-1, so that
 *   (I / (h gamma) - J) u_i = f(y0 + sum_j a_ij u_j) + sum_j (c_ij / h) u_j,
 *   y1 = y0 + sum_i m_i u_i, and the error estimate is sum_i e_i u_i.
 * A stage whose alpha row repeats the previous row's reuses that forcing.
 * The rows of alpha and gamma_ij hold the entries j < i; estimator_order is
 * the order of the embedded pair plus one, as MICM's constructors set it. */
/* alpha and gamma_ij are read only; they are not declared const because C before
 * C23 does not convert a pointer to an array of double into a pointer to an
 * array of const double (gcc -Wpedantic warns). */
static tmc_rosenbrock_parameters from_standard_form(int stages, double alpha[TMC_MAX_STAGES][TMC_MAX_STAGES],
                                                     double gamma_ij[TMC_MAX_STAGES][TMC_MAX_STAGES], double gamma,
                                                     const double *b, const double *bhat, double estimator_order)
{
  tmc_rosenbrock_parameters p = base();
  double inverse[TMC_MAX_STAGES][TMC_MAX_STAGES]; /* Gamma^-1, lower triangular */
  p.stages = stages;
  for (int l = 0; l < stages; ++l)
    for (int i = 0; i < stages; ++i)
    {
      double value = i == l ? 1.0 : 0.0;
      for (int j = l; j < i; ++j)
        value -= gamma_ij[i][j] * inverse[j][l];
      inverse[i][l] = i < l ? 0.0 : value / gamma;
    }
  for (int i = 0; i < stages; ++i)
  {
    const int offset = i * (i - 1) / 2;
    for (int j = 0; j < i; ++j)
    {
      double a = 0.0;
      for (int l = j; l < i; ++l)
        a += alpha[i][l] * inverse[l][j];
      p.a[offset + j] = a;
      p.c[offset + j] = -inverse[i][j];
    }
    bool same_argument = i > 0;
    for (int j = 0; j < i && same_argument; ++j)
      same_argument = alpha[i][j] == alpha[i - 1][j];
    if (i > 0 && alpha[i][i - 1] != 0.0)
      same_argument = false;
    p.new_function_evaluation[i] = !same_argument;
    double row_alpha = 0.0, row_gamma = gamma;
    for (int j = 0; j < i; ++j)
    {
      row_alpha += alpha[i][j];
      row_gamma += gamma_ij[i][j];
    }
    p.alpha[i] = row_alpha;
    p.gamma[i] = row_gamma;
  }
  p.gamma[0] = gamma;
  for (int j = 0; j < stages; ++j)
  {
    double m = 0.0, e = 0.0;
    for (int i = j; i < stages; ++i)
    {
      m += b[i] * inverse[i][j];
      e += (b[i] - bhat[i]) * inverse[i][j];
    }
    p.m[j] = m;
    p.e[j] = e;
  }
  p.estimator_of_local_order = estimator_order;
  return p;
}

tmc_rosenbrock_parameters tmc_rosenbrock_ros34pw2(void)
{
  /* ROS34PW2 (Rang and Angermann, 2005): a four-stage, third-order, L-stable
   * Rosenbrock-W method with an embedded second-order method; it keeps its
   * order with an approximate Jacobian, so a factorization can be reused.
   * Coefficients as printed by Blom et al. (2016, Table 7). */
  const double g = 4.3586652150845900e-01;
  double alpha[TMC_MAX_STAGES][TMC_MAX_STAGES] = { { 0 } }, gamma_ij[TMC_MAX_STAGES][TMC_MAX_STAGES] = { { 0 } };
  alpha[1][0] = 8.7173304301691801e-01;
  alpha[2][0] = 8.4457060015369423e-01;
  alpha[2][1] = -1.1299064236484185e-01;
  alpha[3][0] = 0.0;
  alpha[3][1] = 0.0;
  alpha[3][2] = 1.0;
  gamma_ij[1][0] = -8.7173304301691801e-01;
  gamma_ij[2][0] = -9.0338057013044082e-01;
  gamma_ij[2][1] = 5.4180672388095326e-02;
  gamma_ij[3][0] = 2.4212380706095346e-01;
  gamma_ij[3][1] = -1.2232505839045147e+00;
  gamma_ij[3][2] = 5.4526025533510214e-01;
  const double b[TMC_MAX_STAGES] = { 2.4212380706095346e-01, -1.2232505839045147e+00, 1.5452602553351020e+00, 4.3586652150845900e-01, 0, 0 };
  const double bhat[TMC_MAX_STAGES] = { 3.7810903145819369e-01, -9.6042292212423178e-02, 5.0000000000000000e-01, 2.1793326075422950e-01, 0, 0 };
  return from_standard_form(4, alpha, gamma_ij, g, b, bhat, 3.0);
}

tmc_backward_euler_parameters tmc_backward_euler_defaults(void)
{
  tmc_backward_euler_parameters p;
  memset(&p, 0, sizeof p);
  p.max_number_of_steps = 11;
  p.small = 1.0e-40;
  p.h_start = 0.0;
  p.time_step_reductions_count = 5;
  p.time_step_reductions[0] = 0.5;
  p.time_step_reductions[1] = 0.5;
  p.time_step_reductions[2] = 0.5;
  p.time_step_reductions[3] = 0.5;
  p.time_step_reductions[4] = 0.1;
  return p;
}

/* MICM: fac = min(factor_max, max(factor_min, safety / error^(1/order))) */
double tmc_default_step_controller(double error, const tmc_rosenbrock_parameters *p, void *context)
{
  (void)context;
  double fac = p->safety_factor / pow(error, 1.0 / p->estimator_of_local_order);
  if (fac < p->factor_min)
    fac = p->factor_min;
  if (fac > p->factor_max)
    fac = p->factor_max;
  return fac;
}
