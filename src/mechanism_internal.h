#ifndef CHEM_MECHANISM_INTERNAL_H
#define CHEM_MECHANISM_INTERNAL_H

#include "chemistry.h"

typedef struct {
  size_t species_count;
  size_t dependency_count;
  size_t *offsets;
  size_t *reactions;
  size_t *reactant_terms;
} chem_jacobian_plan;

int chem_jacobian_plan_init(chem_jacobian_plan *plan,
                            const chem_mechanism *mechanism);
void chem_jacobian_plan_destroy(chem_jacobian_plan *plan);
int chem_negative_jacobian_with_plan(
    const chem_mechanism *mechanism, const double *concentration,
    const double *rate_constants, const chem_jacobian_plan *plan,
    double *negative_jacobian);
int chem_negative_jacobian_with_plan_pattern(
    const chem_mechanism *mechanism, const double *concentration,
    const double *rate_constants, const chem_jacobian_plan *plan,
    const size_t *row_counts, const size_t *row_cols,
    double *negative_jacobian);

#endif
