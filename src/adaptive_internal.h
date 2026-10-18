#ifndef CHEM_ADAPTIVE_INTERNAL_H
#define CHEM_ADAPTIVE_INTERNAL_H
#include "internal.h"
#include "chem_adaptive.h"
struct ChemWorkspace {
  const ChemState *owner;
  double *jacobian, *history, *scratch;
  double next_h, previous_error, previous_h, factor_h, times[7];
  size_t history_count, jacobian_age;
  int order, factor_valid;
};
double chem_adaptive_norm(const ChemAdaptive *, const ChemState *, const double *, const double *, const double *);
int chem_adaptive_positive(const ChemAdaptive *, const ChemState *, const double *, const double *);
ChemStatus chem_bdf_solve(const ChemAdaptive *, double, ChemState *, ChemWorkspace *, ChemAdaptiveStats *);
#endif
