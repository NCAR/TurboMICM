/* Rosenbrock integration and error control ported from the reference.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <float.h>
#include <math.h>
#include <string.h>

double chem_normalized_error(const ChemState *s, const double *old, const double *ynew, const double *error)
{
    double sum = 0;
    size_t nc = s->cells, n = s->mechanism->data->species_count;
    /* Strict species, then cell order, matching the reference's error reduction. */
    for (size_t i = 0; i < n; ++i)
        for (size_t c = 0; c < nc; ++c) {
            size_t j = i * nc + c;
            double ymax = fmax(fabs(old[j]), fabs(ynew[j]));
            double e = error[j] / (s->atol[i] + s->rtol * ymax);
            sum += e * e;
        }
    return fmax(sqrt(sum / (n * nc)), 1e-10);
}
#include "rosenbrock_impl.h"
