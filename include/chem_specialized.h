#ifndef CHEM_SPECIALIZED_H
#define CHEM_SPECIALIZED_H

#include "chem.h"

#ifdef __cplusplus
extern "C" {
#endif

/* This function is emitted beside a generated mechanism. Negative tolerance
   values select the configuration's default; a nonnegative absolute tolerance
   replaces every species tolerance for this call. */
const chem_mechanism *chem_generated_mechanism(void);
int chem_generated_integrate_group_steps(
    size_t cell_count,
    double *y_cell_major,
    const double *parameters_cell_major,
    const double *temperature,
    const double *pressure,
    double interval,
    size_t interval_count,
    double relative_tolerance,
    double absolute_tolerance,
    chem_statistics *statistics);

#ifdef __cplusplus
}
#endif

#endif
